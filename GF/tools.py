import numpy as np
import struct 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import tarfile

def parse_in_file(file_name):
    f = open(file_name, "r")
    param=[0]*20
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
        if 'Inst_norm' in line:
            param[13]=float(line.split("=")[1])    
        if 'Beta' in line:
            param[14]=str(line.split("=")[1])  
        if 'Ls' in line:
            param[15]=int(line.split("=")[1])   
        if 'Lr' in line:
            param[16]=int(line.split("=")[1])  
        if 'Tau' in line:
            param[17]=str(line.split("=")[1])  
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
        return (np.array(density),np.array(sizes))

    if (count!=elements):
        print("Four dimensional array wrongly readed")
        return (np.array(density),np.array(sizes))
    
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
    #print("\n")
    #print(density[0])
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

def remove_duplicate(maxima):
    
    temp=[]
    if maxima:
        for i in range(0,len(maxima)):
            duplicate=False
            for j in range(0,i):
                distance=np.sqrt((maxima[i][2][0]-maxima[j][2][0])**2 + (maxima[i][2][1]-maxima[j][2][1])**2)
                if distance < 1:
                    duplicate=True
            if not duplicate:
                temp.append(maxima[i])
            #else:
                #temp.append(maxima[i])
    return(temp)

def find_inst_2d(top,en,sizes,norm_frac,norm_inst,neigh):
    inst=[]
    frac=[]
    a_inst=[]
    a_frac=[]
    total=[]
    
    maxima=find_max_2d(top,sizes)

    for element in maxima:
        i=element[0]
        j=element[1]
        
        frac_q, total_q, ppol, pgaus = condition(top,en,i,j,sizes,norm_frac,norm_inst,neigh)
        if frac_q:
            if top[i,j]>0:
                #print(pfit.append(top[i,j]))
                frac.append([i,j,ppol,pgaus])
            else:
                #print(pfit)
                a_frac.append([i,j,ppol,pgaus])
        if total_q:
            if top[i,j]>0:
                inst.append([i,j,ppol,pgaus])
            else:
                a_inst.append([i,j,ppol,pgaus])
           
    frac=remove_duplicate(frac)
    a_frac=remove_duplicate(a_frac)
    inst=remove_duplicate(inst)
    a_inst=remove_duplicate(a_inst)    
    
    t_frac= a_frac + frac
    t_inst= inst + a_inst 
    total=t_frac+t_inst
    return(inst, a_inst, frac, a_frac, t_frac, t_inst, total)


def condition(density,density_en,i,j,sizes,norm_frac,norm_inst,neigh):

    fit_function="Polynomial"
    fractional=False
    total=False
    if density[i,j]<0:
      ppol, pgauss =fit_inst(-density,[i,j],neigh,sizes,fit_function)
            #return(fractional,total,[0,0,0,0,0,0], [0,0,0,0,0,0])
    elif density[i,j]>0:
      ppol, pgauss=fit_inst(density,[i,j],neigh,sizes,fit_function)  
            #return(fractional,total,[0,0,0,0,0,0], [0,0,0,0,0,0])

    Q=local_q(density,sizes,i,j,1)
    S=local_q(density_en,sizes,i,j,1)

    selfdual=8*np.pi*np.pi*Q/(1.0*S)
    fractional=True
    ppol.append(density[i,j])
    ppol.append(selfdual)

    pgauss.append(density[i,j])
    pgauss.append(selfdual)
    
    return(fractional, total, ppol, pgauss)


def fit_inst(density_2d,maxima,neigh,sizes,fit_function):
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

    try:
      popt_pol, pcov_pol = curve_fit(inst_pol, data[:,:2], data[:,2],
      p0=[maxima[0],maxima[1],1,1])
      popt_pol=list(popt_pol)
      popt_pol.append(np.sum(np.sqrt(np.diag(pcov_pol))))
    except RuntimeError:
      popt_pol=[0,0,0,0,10000000]
    #for i in range(0,len(data)):
      #data[i][2]=np.log(data[i][2])
    try:
      popt_gaus, pcov_gaus = curve_fit(inst_gauss, data[:,:2], data[:,2],
      p0=[maxima[0],maxima[1],1,1])
      popt_gaus=list(popt_gaus)
      popt_gaus.append(np.sum(np.sqrt(np.diag(pcov_gaus))))
    except RuntimeError:
      popt_gaus=[0,0,0,0,1000000000]
    #p0=[maxima[0],maxima[1],1,1])
    
  
    return(popt_pol,popt_gaus)

def inst_pol(position,maxima_x,maxima_y,rho,norm):
    return(norm/(2*np.pi)*(rho**2/((position[:,0]-maxima_x)**2 +
                    (position[:,1]-maxima_y)**2+rho**2)**2))

def inst_gauss_log(position,maxima_x,maxima_y,rho,norm):
    return(norm-1/rho**2*((position[:,0]-maxima_x)**2 +
                    (position[:,1]-maxima_y)**2))
def inst_gauss(position,maxima_x,maxima_y,rho,norm):
    return(norm*np.exp(-1/rho**2*((position[:,0]-maxima_x)**2 +
                    (position[:,1]-maxima_y)**2)))

def inst_plot_pol(position,maxima_x,maxima_y,rho,norm):
    return(norm/(2*np.pi)*(rho**2/((position[0]-maxima_x)**2 +
                    (position[1]-maxima_y)**2+rho**2)**2))

def inst_plot_gauss(position,maxima_x,maxima_y,rho,norm):
    return(norm*np.exp(-1/rho**2*((position[0]-maxima_x)**2 +
                    (position[1]-maxima_y)**2)))
  
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


def plot_dens_2d(file,density_2d,sizes,frac,inst):

    fig = plt.figure(tight_layout=True,figsize=(10,10))
    ax1 = fig.add_subplot(2, 1, 1, projection="3d") 
    X = np.arange(0,sizes[0])
    Y = np.arange(0,sizes[1])
    X, Y = np.meshgrid(X, Y)
    
    #ax1.set_zlim([-0.08,0.08])

    surf = ax1.plot_surface(X, Y,density_2d, rstride=1, cstride=1,
                    cmap="viridis", edgecolor='none')
    ax1.set_title(file);
    
    ax2 = fig.add_subplot(2, 1, 2) 
    ax2.set_aspect("equal")
    ax2.contourf(X, Y, density_2d, cmap = "viridis")
    
    x=[]
    y=[]
    for element in frac:
        x.append(element[0])
        y.append(element[1])
    ax2.scatter(y,x,color="red")
    
    x=[]
    y=[]
    for element in inst:
        x.append(element[0])
        y.append(element[1])
    ax2.scatter(y,x,color="black")
    
    plt.savefig(file.replace(".dat",".png"))
    
    plt.close(fig)

    return()


def compare_fit(list_old,list_new,cap):
    temp_list=list_old.copy()
    founds=[]
    created=[]
    for i in range(0,len(list_old)):
        param_old=list_old[i]
        find=False
        r_min=cap
        for j in range(0,len(list_new)):
            param_temp=list_new[j]
            r=np.sqrt((param_old[2][0]-param_temp[2][0])**2+(param_old[2][1]-param_temp[2][1])**2)
            if r<cap:
                find=True
                if r<r_min:
                    r_min=r
                    ind_temp=j
                    param_new=param_temp 
        temp_list[i]=[param_old[0], param_old[1],param_new[2]]
        founds.append(ind_temp)
        if not find: #it was destroyed
            params_destr=[list_old[i][0],list_old[i][1],[list_old[i][2][0],list_old[i][2][1],0,0]]
            #print(params_destr)
            temp_list[i]=params_destr
    
    #We take care of the ones that were created
    inst=[]
    for element in temp_list:
        inst.append([element[0],element[1]])
    for element in list_new:
        if [element[0],element[1]] not in inst:
            temp_list.append(element)
    
    return(temp_list)
    
def plot_inst(sizes,popt,directory,file,col,sign,ax="None"):
    if ax=="None":
        f, ax = plt.subplots()
    data_plot=np.zeros((sizes[0],sizes[1]))
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            data_plot[i,j]=inst_plot([i,j],popt[0],popt[1],popt[2],popt[3])
            
    if sign=="+":
        data_plot_1d=data_plot.sum(axis=1)
    if sign=="-":
        data_plot_1d=-data_plot.sum(axis=1)
    x=np.arange(0,sizes[0])
    ax.set_ylim([-0.15,0.15])
    ax.plot(x,data_plot_1d, color=col)
    return(ax)

def plot_fit(density_2d_top,sizes,ppol,pgaus, directory, file, maxy):

  plt.xlabel("x")
  plt.ylabel("q(x)")

  data_plot=np.zeros((sizes[0],sizes[1]))
  for i in range(0,sizes[0]):
    for j in range(0,sizes[1]):
      data_plot[i,j]=inst_plot_pol([i,j],ppol[0],ppol[1],ppol[2],ppol[3])
  data_plot_1d=data_plot[:,maxy]
  plt.plot(data_plot_1d, color="green", 
           label="Instanton")

  data_plot=np.zeros((sizes[0],sizes[1]))
  for i in range(0,sizes[0]):
    for j in range(0,sizes[1]):
      data_plot[i,j]=inst_plot_gauss([i,j],pgaus[0],pgaus[1],pgaus[2],pgaus[3])
  data_plot_1d=data_plot[:,maxy]
  plt.plot(data_plot_1d, color="blue", 
           label="Gaussian")

  density_1d=density_2d_top[:,maxy]
  plt.plot(density_1d, color="purple", label="Data")
  plt.legend(loc="upper right")
  plt.savefig(directory+file, dpi=150)
  plt.show()

  return








