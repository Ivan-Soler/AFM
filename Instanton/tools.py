import numpy as np
import struct 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import tarfile

def parse_in_file(file_name):
    f = open(file_name, "r")
    param=[0]*13
    for line in f:
        line=line.replace(" ","")
        line=line.replace("\n","")
        if '#' in line:
            continue
        if 'Directory' in line:
            param[0]=line.split("=")[1]
        if 'File_name' in line:
            param[1]=line.split("=")[1]
        if 'Start_config' in line:
            param[2]=int(line.split("=")[1])
        if 'End_config' in line:
            param[3]=int(line.split("=")[1])
        if 'Read_step' in line:
            param[4]=int(line.split("=")[1])
        if 'Rho' in line:
            param[5]=float(line.split("=")[1])
        if 'Norm' in line:
            param[6]=float(line.split("=")[1])
        if 'Eps_rho' in line:
            param[7]=float(line.split("=")[1])
        if 'Eps_norm' in line:
            param[8]=float(line.split("=")[1])
        if 'Neigh' in line:
            param[9]=int(line.split("=")[1])
        if 'Compare_action' in line:
            param[10]=bool(int(line.split("=")[1]))
        if 'Eps_action' in line:
            param[11]=float(line.split("=")[1])    
        if 'Conf_number' in line:
            param[12]=int(line.split("=")[1])  
    f.close()
    return(param)

def read_top(file):
    sizes=[[]]*4
    precision=0

    count=0
    i=0
    j=0
    k=0
    l=0
    norm=0
    with open(file,"rb") as file:
        for i in range (0,4):
            sizes[i]=int(file.readline())
        precision=int(file.readline())
        colors=int(file.readline())
        elements=sizes[0]*sizes[1]*sizes[2]*sizes[3]
        density=np.zeros((sizes[0],sizes[1],sizes[2],sizes[3]))
        for i in range(0, sizes[0]):
            for j in range(0, sizes[1]):
                for k in range(0, sizes[2]):
                    for l in range(0,sizes[3]):
                        count+=1
                        num=struct.unpack('d', file.read(precision))
                        density[i,j,k,l]=num[0]
                        norm+=num[0]
    if (count==elements):
        return (density,np.array(sizes))

    if (count!=elements):
        print("Four dimensional array wrongly readed")
        return (density,np.array(sizes))
    
def read_top_tony(file):
    sizes=[[]]*4
    precision=0

    count=0
    i=0
    j=0
    k=0
    l=0
    norm=0
    with open(file,"rb") as file:
        colors=int(struct.unpack('i', file.read(4))[0])
        for i in range (0,4):
            sizes[i]=int(struct.unpack('i', file.read(4))[0])
        print(sizes)
        elements=sizes[0]*sizes[1]*sizes[2]*sizes[3]
        density=np.zeros((sizes[0],sizes[1],sizes[2],sizes[3]))
        for i in range(0, sizes[3]):
            for j in range(0, sizes[2]):
                for k in range(0, sizes[1]):
                    for l in range(0,sizes[0]):
                        count+=1
                        num=struct.unpack('d', file.read(8))
                        density[i,j,k,l]=num[0]
                        norm+=num[0]
    if (count==elements):
        return (density,np.array(sizes))

    if (count!=elements):
        print("Four dimensional array wrongly readed")
        return (density,np.array(sizes))
    
def projection_2d(density,sizes):
    #First checks which are the sizes of the smaller directions
    small_sizes=[100,100]
    for i in range(0,len(sizes)):
        if sizes[i]<small_sizes[0]:
            small_sizes[1]=small_sizes[0]
            small_sizes[0]=sizes[i]
            continue
        if sizes[i]<small_sizes[1]:
            small_sizes[1]=sizes[i]
    
    #Then the index of the smaller directions are computed
    index_small=[0,0]
    if small_sizes[0]==small_sizes[1]:
        index_small=np.where(sizes==small_sizes[0])[0]
    else:
        index_small[0]=np.where(sizes==small_sizes[0])[0]
        index_small[1]=np.where(sizes==small_sizes[1])[0]

    #The density is integrated over the small directions
    density_2d=density.sum(axis=(index_small[0],index_small[1]))
    
    #Only the big directions are returned
    sizes_big=np.delete(sizes,index_small)
    
    return(density_2d,sizes_big,index_small)

def local_q(density,sizes,i,j,R):
    #computes the local density around the maxima (i,j)
    q=0

    for ii in range(0,R):
        for jj in range(0,R):
            q+=density[(i+ii-int(R/2))%sizes[0],(j+jj-int(R/2))%sizes[1]]
    return(q)


def condition(density,i,j,sizes,rho,norm,eps_rho,eps_norm,neigh):
    #So far fit only works for positive densities
    fractional=False
    
    if density[i,j]<-0.01:
        try:
            p_fit=fit_inst(-density,[i,j],neigh,sizes)
        except:
            return(fractional,[0,0,0,0])
    elif density[i,j]>0.01:
        try:
            p_fit=fit_inst(density,[i,j],neigh,sizes)  
        except:
            return(fractional,[0,0,0,0])
        #print(p_fit[2],rho,eps_rho,p_fit[2]/rho)
        #print(p_fit[3],norm,eps_norm,p_fit[3]/norm)
    else:
        return(fractional,[0,0,0,0])
    
    if p_fit[2]/rho>1-eps_rho and p_fit[2]/rho<1+eps_rho and p_fit[3]/norm>1-eps_norm and p_fit[3]/norm<1+eps_norm :
        #print('True')
        fractional=True
    #print(fractional)
    return(fractional, p_fit)

def remove_duplicate(maxima):
    
    temp=[]
    for i in range(0,len(maxima)):
        duplicate=False
        for j in range(0,i):
            distance=np.sqrt((maxima[i][2][0]-maxima[j][2][0])**2 + (maxima[i][2][1]-maxima[j][2][1])**2)
            if distance < 1:
                duplicate=True
        if not duplicate:
            temp.append(maxima[i])
    return(temp)

def find_inst_2d(top,act,sizes,rho,norm,eps_rho,eps_norm,norm_action,eps_action,neigh):
    inst=[]
    frac=[]
    a_inst=[]
    a_frac=[]
    total=[]
    
    maxima=find_max_2d(top,sizes)

    for element in maxima:
        i=element[0]
        j=element[1]
        
        condition_q, pfit = condition(top,i,j,sizes,rho,norm,eps_rho,eps_norm,neigh)
        condition_ac, pfit = condition(act,i,j,sizes,rho,norm_action,eps_rho,eps_action,neigh)
        if condition_q and condition_ac:
            if top[i,j]>0:
                frac.append([i,j,pfit])
            else:
                a_frac.append([i,j,pfit])
    
    frac=remove_duplicate(frac)
    a_frac=remove_duplicate(a_frac)
    
    t_frac= a_frac + frac
    t_inst= inst + a_inst 
    total=t_frac+t_inst
                 
    return(inst, a_inst, frac, a_frac, t_frac, t_inst, total)
                 
                 
def find_max_2d(density,sizes):
    #Search for the maxima and minima of the density
    maxima=[]
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            if ((abs(density[i,j]) < abs(density[(i+1)%sizes[0],j])) or (abs(density[i,j]) < abs(density[(i-1)%sizes[0],j]))):
                continue
            if ((abs(density[i,j]) < abs(density[i,(j+1)%sizes[1]])) or (abs(density[i,j]) < abs(density[i,(j-1)%sizes[1]]))):
                continue
            maxima.append([i,j])
                 
    return(maxima)


def plot_dens_2d(file,density_2d,sizes,maxima):

    fig = plt.figure(tight_layout=True,figsize=(10,10))
    ax1 = fig.add_subplot(2, 1, 1, projection="3d") 
    X = np.arange(0,sizes[0])
    Y = np.arange(0,sizes[1])
    X, Y = np.meshgrid(X, Y)

    surf = ax1.plot_surface(X, Y,density_2d, rstride=1, cstride=1,
                    cmap="viridis", edgecolor='none')
    ax1.set_title('Topological charge');
    
    ax2 = fig.add_subplot(2, 1, 2) 
    ax2.set_aspect("equal")
    ax2.contourf(X, Y, density_2d, cmap = "viridis")
    
    x=[]
    y=[]
    for element in maxima:
        x.append(element[0])
        y.append(element[1])
    ax2.scatter(y,x,color="red")

    plt.savefig(file.replace(".dat",".png"))
    
    plt.close(fig)

    return()

def inst(position,maxima_x,maxima_y,rho,norm):
    return(norm/(5*2*np.pi*np.pi/(16*rho**3))*(rho/((position[:,0]-maxima_x)**2 +
                    (position[:,1]-maxima_y)**2+rho**2))**4)

def fit_inst(density_2d,maxima,neigh,sizes):
    x=[]
    y=[]
    data=[]
    for i in range(-neigh,neigh+1):
        for j in range(-neigh,neigh+1):
            #if i!=0 and j!=0:
            x=(int(maxima[0])+i)%sizes[0]
            y=(int(maxima[1])+j)%sizes[1]
            data.append([x,y,density_2d[x,y]])
    data=np.array(data)
    
    popt, pcov = curve_fit(inst, data[:,:2], data[:,2], 
                       p0=np.array((maxima[0],maxima[1],3.7,0.7)))
    
    return(popt)

def compare_fit(list_old,list_new,cap):
    temp_list=list_old.copy()
    founds=[]
    created=[]
    for i in range(0,len(list_new)):
        param_new=list_new[i]
        find=False
        r_min=cap
        for j in range(0,len(list_old)):
            param_old=list_old[j]
            r=np.sqrt((param_old[2][0]-param_new[2][0])**2+(param_old[2][1]-param_new[2][1])**2)
            if r<cap:
                find=True
                if r<r_min:
                    r_min=r
                    ind_temp=[i,j]
                    
            if ind_temp not in founds: #it is not repeated
                temp_list[ind_temp[1]]=[param_old[0], param_old[1],param_new[2]]
                founds.append(ind_temp)
                find=True

            elif ind_temp in founds: # we found two instantons close
                created.append(ind_temp)
                find=True
            
        else:
            temp_list.append(param_new) #it was a new one
    
    #We take care of the ones that were created and have a similar distance to other objects
    print(created)
    print(founds)
    for element_c in created:
        for element_f in founds:
            if element_c[1]==element_f[1]:                  
                created_param=list_new[element_c[0]]
                new_param=temp_list[element_f[0]]
                old_param=temp_list[element_f[1]]
                
                dist_created=(old_param[2][0]-created_param[2][0])**2+(old_param[2][1]-created_param[2][1])**2
                dist_new=(old_param[2][0]-new_param[2][0])**2+(old_param[2][1]-new_param[2][1])**2
                
                if dist_created<dist_new: #If the supposed created is actually the old one
                    print("something")
                    temp_list[element[1]]=[old_param[0], old_param[1],created_param[2]]
                    temp_list.append(new_param)
                else:
                    temp_list.append(created_param)

                    
    #Now we update the ones that were destroyed during the whole flow
    founds=np.array(founds)
    for j in range(0,len(list_old)):
        if j not in founds[:,1]:
            params_destr=[list_old[j][0],list_old[j][1],[list_old[j][2][0],list_old[j][2][1],0,0]]
            temp_list[j]=params_destr
    
    return(temp_list)
    
        