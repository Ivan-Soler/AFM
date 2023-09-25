import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d, Axes3D
import Read
import Maxima_find
#plt.rcParams.update({'font.size': 12})

from scipy.optimize import curve_fit

def plot_all_peaks(folder,density, sizes, label):

    T = np.arange(0, sizes[0], 1)
    X = np.arange(0, sizes[1], 1)
    Y = np.arange(0, sizes[2], 1)
    Z = np.arange(0, sizes[3], 1)

    #Tt, Xx = np.mgrid[0:sizes[0], 0:sizes[1]]
    Tt, Zz =np.mgrid[0:sizes[0], 0:sizes[3]]
    #Tt, Yy = np.mgrid[0:sizes[0], 0:sizes[2]]
    #Yy, Zz= np.mgrid[0:sizes[2], 0:sizes[3]]
    #Xx, Yy = np.mgrid[0:sizes[1], 0:sizes[2]]
    #Xx, Zz = np.mgrid[0:sizes[1], 0:sizes[3]]

    res=Maxima_find.simple(density,sizes)
    print(res)
    print(len(res[0]))
    for i in range (0, len(res[0])):
        fig = plt.figure(figsize=plt.figaspect(2.))
        ax = fig.add_subplot(2, 1, 1, projection='3d')
        #fig,(ax,ax2)=plt.subplots(2,subplot_kw={"projection": "3d"})
        #ax = Axes3D(fig)
        #ax = plt.axes(projection='3d')
        #ax.plot_surface(Tt, Xx, density[:,:,res[2][i],res[3][i]], rstride=1, cstride=1,
        ax.plot_surface(Tt, Zz, density[:,res[1][i],res[2][i],:], rstride=1, cstride=1,
        #ax.plot_surface(Tt, Yy, density[:,res[1][i],:,res[3][i]], rstride=1, cstride=1,
        #ax.plot_surface(Yy, Zz, density[res[0][i],res[1][i],:,:], rstride=1, cstride=1,
        #ax.plot_surface(Xx, Zz, density[res[0][i],:,res[2][i],:], rstride=1, cstride=1,
        #ax.plot_surface(Xx, Yy, density[res[0][i],:,:,res[3][i]], rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
        plt.title("Maximum at: (" +  str(res[0][i]) + ", "  + str(res[1][i])+", "  + str(res[2][i]) + ", "  + str(res[3][i])  + ")")
        file_name="Maximum_" + str(i)
        
        ax = fig.add_subplot(2, 1, 2)
        ax.set_xlabel("t")
        ax.set_ylabel("x")
        #plt.savefig("AFM"+"_"+label, dpi=150, bbox_inches='tight');
        #plt.show()
        #plt.close()
        fig,ax=plt.subplots(1,1)
        #cp = ax.contourf(Tt, Xx, density[:,:,res[2][i],res[3][i]])
        cp = ax.contourf(Tt, Zz, density[:,res[1][i],res[2][i],:])
        #cp= ax.contourf(Tt, Yy, density[:,res[1][i],:,res[3][i]],)
        #ax.contourf(Yy, Zz, density[res[0][i],res[1][i],:,:])
        #cp = ax.contourf(Xx, Yy, density[res[0][i],:,:,res[3][i]])
        #ax.contourf(Xx, Zz, density[res[0][i],:,res[2][i],:])
        #fig.colorbar(cp) 
        #ax.scatter(res[0][i], res[1][i], c="red")
        ax.scatter(res[0][i], res[3][i], c="red")
        #ax.scatter(res[0][i], res[2][i], c="red")
        #ax.scatter(res[0][i], res[3][i], c="red")
        #ax.scatter(res[1][i], res[2][i], c="red")
        #ax.scatter(res[1][i], res[3][i], c="red")
        plt.savefig(folder+label+"_m"+str(i), dpi=150, bbox_inches='tight');
        plt.show()
        plt.close()
    res = Maxima_find.simple(density,sizes)
    mat = np.matrix(res)
    distance=np.zeros([len(res[0]),len(res[0])])

    for i in range(0,len(res[0])):         
        for j in range(0,len(res[0])):
            for k in range(0,4):
                   distance[i][j]+=np.copy((np.abs(res[k][i]-res[k][j])%sizes[k])*(np.abs(res[k][i]-res[k][j])%sizes[k]))
            distance[i][j]=np.copy(np.sqrt(distance[i][j]))
    repeated=0
    for i in range(0,len(res[0])):         
        for j in range(0,len(res[0])):
            if (distance[i][j]<=2 and i!=j and i<j):
                #print(i,j)
                repeated+=1
    print("repeated maxima = ",repeated)

    #big_max=0
    #for i in range(0,len(res[0])):
    #     if (density[res[0][i],res[1][i],res[2][i],res[3][i]])>10000:
    #            big_max+=1
    #print(big_max)


    with open('position_maxima_25.txt','w') as f:
        for line in mat:
            np.savetxt(f, line, fmt='%.2f')
        print("\n number of repeated maxima= "+str(repeated), file=f)

def plot_all_peaks_2d(file,density, sizes, label):

    T = np.arange(0, sizes[0], 1)
    Z = np.arange(0, sizes[3], 1)

   
    Tt, Zz =np.mgrid[0:sizes[0], 0:sizes[3]]
    res=Maxima_find.simple_2d(density,sizes)

    print(len(res[0]))
    for i in range(0, len(res[0])):
        fig = plt.figure(figsize=plt.figaspect(2.))
        ax = fig.add_subplot(2, 1, 1, projection='3d')
        ax.plot_surface(Tt, Zz, density, rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
        ax = fig.add_subplot(2, 1, 2)
        ax.set_xlabel("t")
        ax.set_ylabel("x")
        cp = ax.contourf(Tt, Zz, density)
        ax.scatter(res[0,:], res[1,:], c="red")
        print(file)
        plt.savefig(file+"_"+str(i)+".png", dpi=150, bbox_inches='tight');
        #plt.show()
        plt.close()
    
def plot_all_peaks_1d(file,density, sizes, label):


    Z = np.arange(0, sizes[3], 1)  
    Tt, Zz =np.mgrid[0:sizes[3]]
    res=Maxima_find.simple_1d(density,sizes)

    print(len(res[0]))
    fig = plt.figure(figsize=plt.figaspect(2.))
    ax = fig.add_subplot(2, 1, 1, projection='3d')
    ax.plot_surface(Tt, Zz, density, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
    ax = fig.add_subplot(2, 1, 2)
    ax.set_xlabel("t")
    ax.set_ylabel("x")
    cp = ax.contourf(Tt, Zz, density)
    ax.scatter(res[0,:], res[1,:], c="red")
    print(file)
    plt.savefig(file+".png", dpi=150, bbox_inches='tight');
    #plt.show()
    plt.close()

def plot_density_2d(file,density, sizes):
    
    XX,TT = np.mgrid[0:sizes[0], 0:sizes[3]]

    fig = plt.figure(figsize=plt.figaspect(2.))
    ax = fig.add_subplot(2, 1, 1, projection='3d')
    ax.plot_surface(XX, TT, density, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
    #ax.view_init(elev=20.)
    #ax.tick_params(axis='x', which='major', pad=-2.5)
    #ax.tick_params(axis='y', which='major', pad=-2.5)
    #ax.tick_params(axis='z', which='major', pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    #ax.set_zlim((0,0.001))
    plt.savefig(file+".png", dpi=150, bbox_inches='tight');
    #plt.show()
    plt.close()

    return(fig)

def plot_density(plot_name, density, sizes, position_max):
    
    X = np.arange(0, sizes[0], 1)
    Y = np.arange(0, sizes[1], 1)
    Z = np.arange(0, sizes[2], 1)
    T = np.arange(0, sizes[3], 1)
    
    #YY, Zz =np.mgrid[0:sizes[2], 0:sizes[3]]
    #YY, TT=np.mgrid[0:sizes[1], 0:sizes[3]]
    #XX, YY=np.mgrid[0:sizes[0], 0:sizes[1]]
    XX,TT = np.mgrid[0:sizes[0], 0:sizes[3]]
    
    
    fig = plt.figure()
    ax = Axes3D(fig)
    #ax.set_zlim([-0.003,0.003])
    #ax.set_ylim([-2,8])
    #ax.set_xlim([0,15])
    #ax.plot_surface(YY, Zz, density[round(position_max[0]),round(position_max[1]),:,:], rstride=1, cstride=1,
    #ax.plot_surface(YY, TT, density[round(position_max[0]),:,round(position_max[2]),:], rstride=1, cstride=1,
    #ax.plot_surface(XX, YY, density[:,:,round(position_max[2]),round(position_max[3])], rstride=1, cstride=1,
    ax.plot_surface(XX, TT, density[:,round(position_max[1]),round(position_max[2]),:], rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
    ax.view_init(elev=20.)
    ax.set_ylabel("x")

    return(fig)

def plot_int_3d(plot_name, density, ax_x, ax_y, sizes):
    print('creting slightly bigger grid')
    approximation="linear"
    max_parameter=40
    
    X = np.arange(0, sizes[0], 1)
    Y = np.arange(0, sizes[1], 1)
    Z = np.arange(0, sizes[2], 1)
    T = np.arange(0, sizes[3], 1)


    interp = RegularGridInterpolator((X, Y, Z, T), density, method=approximation)
    x_points=(sizes[0]-1)*5j+1j
    y_points=(sizes[1]-1)*5j+1j
    z_points=(sizes[2]-1)*5j+1j
    t_points=(sizes[3]-1)*5j+1j
 
    
    points_big=Xbig, Ybig, Zbig, Tbig = np.mgrid[0:sizes[0]-1:x_points, 0:sizes[1]-1:y_points, 0:sizes[2]-1:z_points, 0:sizes[3]-1:t_points]
    print('before interpolating')
    density=interp((Xbig,Ybig,Zbig,Tbig))

    Xx, Yy = np.mgrid[0:sizes[0]-1:x_points, 0:sizes[1]-1:y_points]

    print('searching for maxima next')
    data_max=maximum_filter(density, max_parameter, mode = 'constant', cval=0.0) 
    maxima=(density==data_max)
    res = np.where(density == data_max)

    for i in range (0, len(res[1])):
        fig = plt.figure()
        ax = Axes3D(fig)
        #ax.plot_surface(Xx, Yy, density[:,:,res[2][i],res[3][i]], rstride=1, cstride=1,
        ax.plot_surface(Xx, Yy, density[:,res[1][i],res[2][i],:], rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
        plt.title("Maximum at: (" +  str(res[0][i]/5) + ", "  + str(res[1][i]/5)+", "  + str(res[2][i]/5) + ", "  + str(res[3][i]/5)  + ")")
        #plt.title("Maximum at: (" +  str(res[0][i]) + ", "  + str(res[1][i])+", "  + str(res[2][i]) + ", "  + str(res[3][i])  + ")")
        file_name="Maximum_" + str(i)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        plt.savefig(plot_name+"_"+file_name, dpi=150, bbox_inches='tight')
        
        fig,ax=plt.subplots(1,1)
        cp = ax.contourf(Xx, Yy, density[:,res[1][i],res[2][i],:])
        fig.colorbar(cp) 
        ax.set_xlabel(ax_x)
        ax.set_ylabel(ax_y)
        #plt.scatter(res[0][i]/5, res[1][i]/5, c="red")
        plt.scatter(res[0][i]/5, res[3][i]/5, c="red")
        plt.savefig(plot_name+ "_"+file_name+"_contour", dpi=150, bbox_inches='tight')
        print("The first maxima is located at "  + str(res[0][i]/5) + ", "  + str(res[1][i]/5)+ ", " +str(res[2][i]/5) + ", "  + str(res[3][i]/5) )
       # print("The first maxima is located at "  + str(res[0][i]) + ", "  + str(res[1][i])+ ", " +str(res[2][i]) + ", "  + str(res[3][i]) )
        #print(data_max)
    return

def plot_nint_3d(plot_name, density, ax_x, ax_y, sizes):
    max_parameter=40
    
    X = np.arange(0, sizes[0], 1)
    Y = np.arange(0, sizes[1], 1)
    Z = np.arange(0, sizes[2], 1)
    T = np.arange(0, sizes[3], 1)

    Xx, Yy = np.mgrid[0:sizes[0], 0:sizes[1]]
 
    density=density
    print('searching for maxima next')
    data_max=maximum_filter(density, max_parameter, mode = 'constant', cval=0.0) 
    maxima=(density==data_max)
    res = np.where(density == data_max)

    for i in range (0, len(res[1])):
        ax.plot_surface(Xx, Yy, density[:,:,res[2][i],res[3][i]], rstride=1, cstride=1,
       # ax.plot_surface(Xx, Yy, density[:,:,res[2][i],res[3][i]], rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
        plt.title("Maximum at: (" +  str(res[0][i]) + ", "  + str(res[1][i])+", "  + str(res[2][i]) + ", "  + str(res[3][i])  + ")")
        file_name="Maximum_" + str(i)
        ax.set_xlabel(ax_x)
        ax.set_ylabel(ax_y)
        plt.savefig(plot_name+"_"+file_name, dpi=150, bbox_inches='tight')
        
        fig,ax=plt.subplots(1,1)
        cp = ax.contourf(X, Y, density[:,:,res[2][i],res[3][i]])
        fig.colorbar(cp) 
        ax.set_xlabel(ax_x)
        ax.set_ylabel(ax_y)
        plt.scatter(res[0][i], res[3][i], c="red")
        #plt.scatter(res[0][i], res[3][i], c="red")
        plt.savefig(plot_name+ "_"+file_name+"_contour", dpi=150, bbox_inches='tight')
        print("The first maxima is located at "  + str(res[0][i]) + ", "  + str(res[1][i])+ ", " +str(res[2][i]) + ", "  + str(res[3][i]) )
       
    return

def plot_density_1d(data, fig, color='k', legend='', title='',  xlab='', ylab='',
                line_width=1, alpha=1,
                subplots=False, show_grid=True,row=1,column=1,position=1):
    #plt.legend(loc="upper left", prop={'size': 10})
    xpoints=np.arange(0,len(data))
    #if (subplots):
    #    plt.subplot(row,column,position) # here is where you add the subplot to f
    plt.plot(xpoints, data, linewidth=line_width,
        alpha=alpha, color=color, label=legend)
    plt.xlim([min(xpoints), max(xpoints)])
    plt.xlabel("")
    plt.ylabel("")
    plt.grid(show_grid)
    plt.title(title, size=16)
    #plt.savefig(fig,dpi=150, bbox_inches='tight')
    #plt.show()
    #plt.close()
    return

def shifting_coordinates(density,sizes):
    temp=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    #temp[(i+int(sizes[0]/2))%sizes[0],(j+1)%sizes[1],(k+1)%sizes[2],(l-4)%sizes[3]]=density[i,j,k,l]
                    temp[(i)%sizes[0],j,k,(l+4)%sizes[3]]=density[i,j,k,l]
    return temp

def gif_2d(directory,name,pref, t_init, t_final, t_steps, m):
    i=0
    for t in range(t_init, t_final, t_steps):
        time=str(t/10)
        if not (t%10):
            time=str(int(t/10))
        file_name=directory+name+time+pref

        density,sizes=Read.density_2d(file_name)
        
        res=Maxima_find.simple_2d(density,sizes)

        figure=plot_density_2d(density, sizes)
        figure.savefig(directory+"3d/-"+str(int(t/4))+".png", dpi=150)
            
    t_points=np.linspace(t_init,t_final,int((t_final-t_init)/t_steps))
    
    from PIL import Image
    files = []
    for i in range(t_init,t_final,4):
        seq = str(i*2)
        file_names = "-"+str(int(i/4))+'.png'
        files.append(file_names)

    # Create the frames
    frames = []
    files
    for i in files:
        new_frame = Image.open(directory+"3d/"+i)
        frames.append(new_frame)

    # Save into a GIF file that loops forever   
    frames[0].save('3d_vis.gif', format='GIF',
                   append_images=frames[1:],
                   save_all=True,
                   duration=60, loop=0)


def gif(directory,name,pref, t_init, t_final, t_steps, m):
    i=0
    for t in range(t_init, t_final, t_steps):
        time=str(t/10)
        if not (t%10):
            time=str(int(t/10))
        file_name=directory+name+time+pref
        # plotting line
        density,sizes=Read.density(file_name)
        density=-density
        density=shifting_coordinates(density,sizes)

        res=Maxima_find.simple(density,sizes)
        position_max=[res[0][m],res[1][m],res[2][m],res[3][m]]
        count=0
        for element in position_max:
            if element>0:
                count=count+1
        #print(time)
        #print(position_max)
        #print(count)
        #print()
        i=i+1
        if count!=0:
            figure=plot_density("", density, sizes, position_max)
            figure.savefig("3d/-"+str(int(t/4))+".png", dpi=150)
            
    t_points=np.linspace(t_init,t_final,int((t_final-t_init)/t_steps))
    
    from PIL import Image
    files = []
    for i in range(t_init,t_final,4):
        seq = str(i*2)
        file_names = "-"+str(int(i/4))+'.png'
        files.append(file_names)

    # Create the frames
    frames = []
    files
    for i in files:
        new_frame = Image.open("3d/"+i)
        frames.append(new_frame)

    # Save into a GIF file that loops forever   
    frames[0].save('3d_vis.gif', format='GIF',
                   append_images=frames[1:],
                   save_all=True,
                   duration=60, loop=0)
    

def plot_density_back(plot_name, density, sizes, position_max):
    T = np.arange(0, sizes[0], 1)
    X = np.arange(0, sizes[1], 1)
    Y = np.arange(0, sizes[2], 1)
    Z = np.arange(0, sizes[3], 1)
    #YY, Zz =np.mgrid[0:sizes[2], 0:sizes[3]]
    YY, TT=np.mgrid[0:sizes[2], 0:sizes[0]]
    fig = plt.figure()
    ax = Axes3D(fig)
    #ax.plot_surface(YY, Zz, density[round(position_max[0]),round(position_max[1]),:,:], rstride=1, cstride=1,
    ax.plot_surface(YY, TT, density[:,round(position_max[1]),:,round(position_max[3])], rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
    #plt.title("Maximum at: (" +  str(position_max[0]) + ", "  + str(position_max[1])+", "  + str(position_max[2]) + ", "  + str(position_max[3])  + ")")
    #file_name="Maximum_" + str(i)
    ax.set_xlabel("t")
    ax.set_ylabel("x")
    #plt.savefig("AFM"+"_"+str(i), dpi=150, bbox_inches='tight');
    #plt.show()
    #plt.close()
    #fig,ax=plt.subplots(1,1)
    #cp = ax.contourf(Tt, Zz, density[:,density[1],density[2],:])
    #fig.colorbar(cp) 
    #ax.set_xlabel("t")
    #ax.set_ylabel("x")
    #plt.scatter(res[0][i], res[3][i], c="red")
    #plt.savefig("test", dpi=150, bbox_inches='tight');
    #plt.show()
    #plt.close()

    return(fig)

def gap_C(conf, file, measurement):
    
    ev=[]
    for line in open(file).readlines():
        parts = line.split(":")
       
        if len(parts)>4 and parts[1] == conf and parts[2]==measurement :
            ev.append(float(parts[8]))

    x=np.arange(0,len(ev))
    plt.scatter(x,ev)
    return()


def gap_R(conf, file, measurement):
    
    ev=[]
    for line in open(file).readlines():
        parts = line.split(":")
        #print(parts)
       
        if len(parts)>4 and parts[1] == conf and parts[2]==measurement :
            entry=parts[6]
            entry=entry.replace("(", "")
            entry=entry.replace(")", "")
            Re_Im=entry.split(",")
            ev.append(float(Re_Im[0]))
            

    x=np.arange(0,len(ev))
    plt.scatter(x,ev)
    return()

def gap_M(conf, file, measurement):
    
    ev=[]
    for line in open(file).readlines():
        parts = line.split(":")
        #print(parts)
       
        if len(parts)>4 and parts[1] == conf and parts[2]==measurement :
            entry=parts[6]
            entry=entry.replace("(", "")
            entry=entry.replace(")", "")
            Re_Im=entry.split(",")
            ev.append(float(Re_Im[0]))
            

    x=np.arange(0,len(ev))
    plt.scatter(x,ev)
    return()