import numpy as np
import struct 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def parse_in_file(file_name):
    f = open(file_name, "r")
    param=[0]*8
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
        if 'rho' in line:
            param[5]=float(line.split("=")[1])
        if 'norm' in line:
            param[6]=float(line.split("=")[1])
        if 'eps_rho' in line:
            param[7]=float(line.split("=")[1])
        if 'eps_norm' in line:
            param[8]=float(line.split("=")[1])
        if 'Size_local' in line:
            param[9]=int(line.split("=")[1])
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


def find_max_2d(density,sizes,rho,norm,eps_rho,eps_norm,size_local,neigh):
    #Search for the maxima and minima of the density and categorize them in fractional or |Q|=1 instantons
    inst=[]
    frac=[]
    a_inst=[]
    a_frac=[]
    total=[]

    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            if ((abs(density[i,j]) < abs(density[(i+1)%sizes[0],j])) or (abs(density[i,j]) < abs(density[(i-1)%sizes[0],j]))):
                continue
            if ((abs(density[i,j]) < abs(density[i,(j+1)%sizes[1]])) or (abs(density[i,j]) < abs(density[i,(j-1)%sizes[1]]))):
                continue
            #So far fit only works for positive densities
            if density[i,j]<-0.01:
                p_fit,q,p_fit2,pfit_3=fit_inst(-density,[i,j],neigh,sizes)  
            elif density[i,j]>0.01:
                p_fit,q,p_fit2,pfit_3=fit_inst(density,[i,j],neigh,sizes)  
            else:
                continue
            print(p_fit[2],p_fit[3])
            if p_fit[2]/rho>1-eps_rho and p_fit[2]/rho<1+eps_rho and p_fit[3]/norm>1-eps_norm and p_fit[3]/norm<1+eps_norm :
                
                total.append([i,j,q])
                if q>0:
                    frac.append([i,j,q])
                else:
                    a_frac.append([i,j,q])

    t_frac= a_frac + frac
    t_inst= t_frac + inst + a_inst
    return(inst, a_inst, frac, a_frac, t_frac, t_inst, total)

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
    return(norm/(5*2*2*np.pi*np.pi/(16*rho**3))*(rho/((position[:,0]-maxima_x)**2 +
                    (position[:,1]-maxima_y)**2+rho**2))**4)

#def inst(position,maxima_x,maxima_y,rho,norm):
#    return( rho*((position[:,0]-maxima_x)**2 + (position[:,1]-maxima_y)**2)
#           + norm)

def inst_frac_norm(position,maxima_x,maxima_y,rho):
    return(1/(5*2*2*np.pi*np.pi/(16*rho**3))*(rho/((position[:,0]-maxima_x)**2 +
                    (position[:,1]-maxima_y)**2+rho**2))**4)

def inst_norm(position,maxima_x,maxima_y,rho):
    return(6/(np.pi*np.pi*rho**4)*(rho**2/((position[:,0]-maxima_x)**2 +
                    (position[:,1]-maxima_y)**2+rho**2))**4)

def inst_height(position,maxima_x,maxima_y,rho):
    return((rho**2/((position[:,0]-maxima_x)**2 +
                    (position[:,1]-maxima_y)**2+rho**2))**4)

def q0(rho):
    #maxim=[23,23]
    return(6/(np.pi*np.pi*rho**4))

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
                       p0=np.array((maxima[0],maxima[1],3,60)))
    q=q0(popt[2])
    
    popt2, pcov2 = curve_fit(inst_frac_norm, data[:,:2], data[:,2], 
                       p0=np.array((maxima[0],maxima[1],2)))
    
    popt3, pcov3 = curve_fit(lambda position, rho:inst_height(position,maxima[0],maxima[1],rho), data[:,:2], data[:,2])
        
    q=q0(popt[2])
    
    return(popt,q,popt2,popt3)

