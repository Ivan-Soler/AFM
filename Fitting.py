import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import Read
import Maxima_find
import Plotting
import copy
plt.rcParams.update({'font.size': 12})

def instanton_parameter(directory,gf_zm,name,pref,guess_number_maxima,t_init,t_final,t_step,X1,rho1,height1):   

    Cap_number_maxima=guess_number_maxima
    
    #Variables to store the parameters
    Xmax_final=[[[0 for k in range(4)] for j in range(Cap_number_maxima)] for i in range(t_init,t_final,t_step)]
    height_final=[[[0 for k in range(4)] for j in range(Cap_number_maxima)] for i in range(t_init,t_final,t_step)]
    rho_final=[[[0 for k in range(4)] for j in range(Cap_number_maxima)] for i in range(t_init,t_final,t_step)]
    Xmax_previous=[[]]
    Xmax_disordered=[[0 for k in range(4)] for j in range(Cap_number_maxima)]
    height_disordered=[[0 for k in range(4)] for j in range(Cap_number_maxima)]
    rho_disordered=[[0 for k in range(4)] for j in range(Cap_number_maxima)]
    # Just tries to count the time as it appears in the file naming.
    t=0
    l=0
    res1=[[]]
    for i in range(t_init,t_final,t_step):
        time=str(i/10)
        if not (i%10):
            time=str(int(i/10))
        if (gf_zm):
            density,sizes=Read.density(directory+name+time+pref)
        else:
            density,sizes=Read.susy_mode(directory+name+time+pref)
        if ((density.sum())<0):
            density=-density
        res=copy.deepcopy(Maxima_find.simple(density,sizes)) #find the maxima
        #print(len(res[0]))
        #print(res)
        if (len(res[0])<=Cap_number_maxima and len(res[0])>0): # Cap to smooth enough configuration by the number of maxima
            rho=[]
            Xmax=[]
            height=[]
            for j in range(0,len(res[0])): # For every maxima we do a fit.
                Xmax_find=np.array([round(res[0][j]),round(res[1][j]),round(res[2][j]),round(res[3][j])])
                rho_temp=np.array([10,0,0,0])
                Xmax_temp=np.array([10,0,0,0])
                height_temp=np.array([10,0,0,0])
                for d in (0,1,2,3): # the fit is done independently in each direction
                    rho_temp[d], Xmax_temp[d], height_temp[d]=fitting_instanton(density,d,sizes,Xmax_find,0)
                #print(j,Xmax_temp)
                #Xmax_disordered[j]=copy.deepcopy(Xmax_temp)
                Xmax_disordered[j]=np.copy(Xmax_find)
                height_disordered[j]=np.copy(height_temp)
                rho_disordered[j]=np.copy(rho_temp)
            #print(i,Xmax_disordered)
            if (l==0 and not X1): #if it is not the first configuration and we do not have a set of maxima, we define it
                Xmax_final[t]=np.copy(Xmax_disordered)
                height_final[t]=np.copy(height_disordered)
                rho_final[t]=np.copy(rho_disordered)
                res1=res
            elif (l==0):
                order_index=order_maxima(len(res[0]),  [Xmax_disordered,height_disordered,rho_disordered], [X1,height1,rho1], len(X1), sizes)
                #print(order_index)
                for p in range(0, len(res[0])):
                    #print(order_index[p])
                    Xmax_final[t][order_index[p]]=Xmax_disordered[p]
                    height_final[t][order_index[p]]=height_disordered[p]
                    rho_final[t][order_index[p]]=rho_disordered[p]
                Xmax_previous=np.copy(Xmax_final[t])
                res1=res
            else:
                order_index=order_maxima(len(res[0]),  [Xmax_disordered,height_disordered,rho_disordered], [Xmax_final[t-1],height_final[t-1],rho_final[t-1]], len(res1[0]), sizes)
                for p in range(0, len(res[0])):
                    #print(order_index[p])
                    Xmax_final[t][order_index[p]]=Xmax_disordered[p]
                    height_final[t][order_index[p]]=height_disordered[p]
                    rho_final[t][order_index[p]]=rho_disordered[p]
                Xmax_previous=np.copy(Xmax_final[t])
                res1=res
            l=l+1
        t=t+1
    #print(t)
    #print(l)
    return(Xmax_final, rho_final, height_final)

def periodic_ind(shift, size):
    ind=shift%size
    return ind

def fitting_instanton(density, d, sizes, Xmax, plot): #needs the density, the direction of the section d, the position of the maxima Xmax, returns the coefficient in front of x^2
    #1-D parabola
    def func(x, h, Xmax_extr, rho):
        return h - 1/(rho**2)*(x-Xmax_extr)**2
    
    # pick up one maxima
    #Xmax_int=[res[0][0]/5,res[1][0]/5,res[2][0]/5,res[3][0]/5]
    #Xmax=[round(Xmax_int[0]),round(Xmax_int[1]),round(Xmax_int[2]),round(Xmax_int[3])]
    
    # prepare the data: we will use 3 points, the Xmax and the nearest neighbours
    x=np.array([Xmax[d]-1,Xmax[d],Xmax[d]+1])
    y=np.array([Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4]])
    z=np.array([Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4]])
    t=np.array([Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4]])
    shift= { d:x, (d+1)%4:y, (d+2)%4:z, (d+3)%4:t} # dictionary from integer to axis
    data=density[periodic_ind(shift[0],sizes[0]),periodic_ind(shift[1],sizes[1]),
                     periodic_ind(shift[2],sizes[2]),periodic_ind(shift[3],sizes[3])]
    
    try:
        popt,pcov= curve_fit(func, x, data)
        
    except RuntimeError:
        print("error")
        return [0, 0, 0]
    rho=popt[2]
    Xmax_extr=popt[1]
    height=popt[0]
    
    #xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
    #plt.plot(xdata, func(xdata, *popt))
    #plt.plot(x, data)
    #plt.show()
    
    if plot:
        plt.plot(x, data)
        x=np.array([Xmax[d]-2,Xmax[d]-1, Xmax[d],Xmax[d]+1, Xmax[d]+2])
        y=np.array([Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4]])
        z=np.array([Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4]])
        t=np.array([Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4]])
        shift= { d:x, (d+1)%4:y, (d+2)%4:z, (d+3)%4:t} # dictionary from integer to axis
        data=density[periodic_ind(shift[0],sizes[0]),periodic_ind(shift[1],sizes[1]),
                     periodic_ind(shift[2],sizes[2]),periodic_ind(shift[3],sizes[3])]
        
        xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
        plt.plot(xdata, func(xdata, *popt))
        #plt.show()
        
    return popt

def fitting_instanton_impr(density, d, sizes, Xmax, plot): #needs the density, the direction of the section d, the position of the maxima Xmax, returns the coefficient in front of x^2
    #1-D parabola
    def func(x, h, Xmax_extr, rho):
        return h - 1/(rho**2)*(x-Xmax_extr)**2
    
    # prepare the data: we will use 5 points, the Xmax and the nearest neighbours

    x=np.array([Xmax[d]-2,Xmax[d]-1, Xmax[d],Xmax[d]+1, Xmax[d]+2])
    y=np.array([Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4]])
    z=np.array([Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4]])
    t=np.array([Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4]])
    shift= { d:x, (d+1)%4:y, (d+2)%4:z, (d+3)%4:t} # dictionary from integer to axis
    data=density[periodic_ind(shift[0],sizes[0]),periodic_ind(shift[1],sizes[1]),
                 periodic_ind(shift[2],sizes[2]),periodic_ind(shift[3],sizes[3])]

    xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
    #plt.plot(xdata, func(xdata, *popt))
    #plt.plot(x, data)
    #plt.show()
    try:
        popt, pcov = curve_fit(func, x, data)
        #popt = curve_fit(func, x, data)
    except RuntimeError:
        print("error")
        return [0, 0, 0], [0, 0, 0]
        #return 0,0,0,[0, 0, 0]
    rho=popt[2]
    Xmax_extr=popt[1]
    height= popt[0]
    
    #xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
    #plt.plot(xdata, func(xdata, *popt))
    #plt.plot(x, data)
    #plt.show()
    
    if plot:
        x=np.array([Xmax[d]-2,Xmax[d]-1, Xmax[d],Xmax[d]+1, Xmax[d]+2])
        y=np.array([Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4]])
        z=np.array([Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4]])
        t=np.array([Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4]])
        shift= { d:x, (d+1)%4:y, (d+2)%4:z, (d+3)%4:t} # dictionary from integer to axis
        data=density[periodic_ind(shift[0],sizes[0]),periodic_ind(shift[1],sizes[1]),
                     periodic_ind(shift[2],sizes[2]),periodic_ind(shift[3],sizes[3])]
        
        xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
        plt.plot(xdata, func(xdata, *popt))
        plt.plot(x, data)
        #plt.show()
    return (popt, pcov)

    #return (height, Xmax_extr, rho, pcov)

def instanton_parameter_np(directory,gf_zm,name,pref,guess_number_maxima,t_init,t_final,t_step,X1,rho1,height1):   

    Cap_number_maxima=guess_number_maxima*10
    #Variables to store the parameters
    Xmax_final=np.zeros((int((t_final-t_init)/t_step),Cap_number_maxima, 4))
    height_final=np.zeros((int((t_final-t_init)/t_step),Cap_number_maxima, 4))
    rho_final=np.zeros((int((t_final-t_init)/t_step),Cap_number_maxima, 4))
    Xmax_previous=np.zeros((Cap_number_maxima, 4))
    Xmax_disordered=np.zeros((Cap_number_maxima, 4))
    height_disordered=np.zeros((Cap_number_maxima, 4))
    rho_disordered=np.zeros((Cap_number_maxima, 4))
    # Just tries to count the time as it appears in the file naming.
    t=0
    res1=[[]]
    for i in range(t_init,t_final,t_step):
        time=str(i/10)
        if not (i%10):
            time=str(int(i/10))
        if (gf_zm):
            density,sizes=Read.density(directory+name+time+pref)
            density*=-density
        else:
            density,sizes=Read.susy_mode(directory+name+time+pre)
            density*=-density
        if ((density.sum())<0):
            density=-density
        res=Maxima_find.simple(density,sizes) #find the maxima
        #print(len(res[0]))
        #print(res)
        if (len(res[0])>0): # Just if there is some maxima otherwise it breaks
            rho=[]
            Xmax=[]
            height=[]
            for j in range(0,len(res[0])): # For every maxima we do a fit.
                Xmax_find=[round(res[0][j]),round(res[1][j]),round(res[2][j]),round(res[3][j])]
                #print(len(res[0]))
                rho_temp=[10,0,0,0]
                Xmax_temp=[10,0,0,0]
                height_temp=[10,0,0,0]
                for d in (0,1,2,3): # the fit is done independently in each direction
                    rho_temp[d], Xmax_temp[d], height_temp[d]=fitting_instanton(density,d,sizes,Xmax_find,0)
                Xmax_disordered[j]=np.copy(Xmax_temp)
                height_disordered[j]=np.copy(height_temp)
                rho_disordered[j]=np.copy(rho_temp)
            if (t==0 and not np.size(X1)):
                Xmax_final[t]=np.copy(Xmax_disordered)
                height_final[t]=np.copy(height_disordered)
                rho_final[t]=np.copy(rho_disordered)
                res1=np.copy(res)
            elif (t==0 and np.size(X1)):
                order_index=order_maxima(len(res[0]),  [Xmax_disordered,height_disordered,rho_disordered], [X1,height1,rho1], len(X1), sizes)
                #print(order_index)
                for p in range(0, len(res[0])):
                    #print(order_index[p])
                    Xmax_final[t][order_index[p]]=np.copy(Xmax_disordered[p])
                    height_final[t][order_index[p]]=np.copy(height_disordered[p])
                    rho_final[t][order_index[p]]=np.copy(rho_disordered[p])
                Xmax_previous=np.copy(Xmax_final[t])
                res1=res
            else:
                order_index=order_maxima(len(res[0]),  [Xmax_disordered,height_disordered,rho_disordered], [Xmax_final[t-1],height_final[t-1],rho_final[t-1]], len(res1[0]), sizes)
                for p in range(0, len(res[0])):
                    #print(order_index[p])
                    Xmax_final[t][order_index[p]]=np.copy(Xmax_disordered[p])
                    height_final[t][order_index[p]]=np.copy(height_disordered[p])
                    rho_final[t][order_index[p]]=np.copy(rho_disordered[p])
                Xmax_previous=np.copy(Xmax_final[t])
                res1=res
            t=t+1
        #print(t)
            #print(res)
            #print()
    return(Xmax_final, rho_final, height_final)
 

def distance_param(param_in, param_out,sizes):
    #print(param_in)
    return((np.abs(param_in[0]-param_out[0]))%(sizes[0]) + (np.abs(param_in[1]-param_out[1]))%(sizes[1]) + (np.abs(param_in[2]-param_out[2]))%(sizes[2]) + (np.abs(param_in[3]-param_out[3]))%(sizes[3]))

def distance_param2(param_in, param_out,sizes):
    return(np.abs(param_in[0]-param_out[0])+ np.abs(param_in[1]-param_out[1]) + np.abs(param_in[2]-param_out[2]) + np.abs(param_in[3]-param_out[3]))

def order_maxima(num_max_1, param_disordered, param_previous, num_max_2, sizes):
    #print(num_max_1, num_max_2)
    tmp_max_1=max(num_max_1,num_max_2)
    Max_range=np.arange((tmp_max_1))
    tmp_max_2=min(num_max_1,num_max_2)
    index=[[]]*tmp_max_2
    blocked_list=[]
    #print(len(index),np.size(Max_range))
    #print(num_max_1, num_max_2)
    for p in range(0,tmp_max_2): # We compare with the previous parameters to identify the maxima and see how the parameters change
        distance=1000000
        for l in range(0,tmp_max_1):
            distance_temp=distance_param(param_disordered[0][p],param_previous[0][l],sizes)#+distance_param2(param_disordered[1][p],param_previous[1][l],sizes)+distance_param2(param_disordered[2][p],param_previous[2][l],sizes)
            if (distance_temp<=distance and l not in blocked_list): #and Xmax_previous[l][0]>0 and Xmax_disordered[p][0]>0): # if the point is closer we substitute it
                index[p]=l #We store the position on the matrix 
                distance=distance_temp
        blocked_list.append(index[p])
    # We correct for the missing ones
    if (tmp_max_1==num_max_1):
        missing=list((set(Max_range).difference(index)))
        joined_index=index+missing
    else:
        joined_index=index
    #print(joined_index,'\n')
    #if (num_max_1>num_max_2):
    #    index=[0]*num_max_1
    #    print("bla")
    #if (num_max_1==0):  
    #    index=[0]
    #    print("ble")
    #print(index)
    return joined_index

def instanton_evolution_new():
    title=" "
    directory="./q0/gf_unbug/"
    name="profile4dt"
    legend4='wilson'
    pref="c100to.dat"
    max_guess=8
    gf_zm=True
    t_init=10
    t_final=14
    t_steps=4
    for t in range(t_init,t_final,t_steps):
        file=directory+name+str(t)+pref
        density,sizes=Read.density(file)
        density=np.abs(density)
        res=Maxima_find.improve(density,sizes)
        print(res)
        for i in range(0,len(res[0]),1):
            for d in range(0,3,1):
                Xmax_find=[round(res[0][i]),round(res[1][i]),round(res[2][i]),round(res[3][i])]
                rho,Xmax,height=fitting_instanton(density,d,sizes,Xmax_find,0)
    print("finished")
    return(rho,Xmax,height)

def instanton_evolution():
    title=" "
    
    directory="./notwist/wilson/"
    #directory="./heated/heated25/epsilon_1/"
    #directory="./new/epsilon_1/"
    name="profile4dt"
    legend4='wilson'
    pref="c100to.dat"
    max_guess=8
    gf_zm=True
    t_init=40
    t_final=360
    t_steps=4
    X4,rho4,height4=instanton_parameter_np(directory,gf_zm,name,pref,max_guess,t_init,t_final,t_steps,[],[],[])
    print("finished")
    
    #directory="./notwist/wilson/"
    directory="./heated/heated25/wilson/"
    #directory="./new/wilson/"
    legend1='wilson'
    name="profile4dt"
    pref="c100to.dat"

    X1,rho1,height1=instanton_parameter_np(directory,gf_zm,name,pref,max_guess,t_init,t_final,t_steps,[],[],[])
    print("finished")
    
    #directory="./notwist/symanzik/"
    directory="./heated/heated25/symanzik/"
    #directory="./new/symanzik/"
    name="profile4dt"
    legend2='symanzik'
    pref="c100to.dat"
    gf_zm=True
    X2,rho2,height2=instanton_parameter_np(directory,gf_zm,name,pref,max_guess,t_init,t_final,t_steps,[],[],[])
    
    #directory="./notwist/overimproved/"
    directory="./heated/heated25/overimproved/"
    #directory="./new/overimproved/"
    name="profile4dt"
    legend3='too_overimproved'
    pref="c100to.dat"
    gf_zm=True
    X3,rho3,height3=instanton_parameter_np(directory,gf_zm,name,pref,max_guess,t_init,t_final,t_steps,[],[],[])
    
    t_points=[]
    for i in range(t_init,t_final,t_steps):
        t_points.append(np.sqrt(8*i/100))
    for maxima_number in range(0,max_guess):
        f=plt.figure(figsize=(20,15))
        for direction in range(1,2):
            X1_plot=[]
            X2_plot=[]
            X3_plot=[]
            X4_plot=[]
            rho1_plot=[]
            rho2_plot=[]
            rho3_plot=[]
            rho4_plot=[]
            height1_plot=[]
            height2_plot=[]
            height3_plot=[]
            height4_plot=[]

            for t in range(0,len(X1)):
                try:
                    X1_plot.append(X1[t][maxima_number][direction])
                    X2_plot.append(X2[t][maxima_number][direction])
                    rho1_plot.append(rho1[t][maxima_number][direction])
                    rho2_plot.append(rho2[t][maxima_number][direction])
                    height1_plot.append(height1[t][maxima_number][direction])
                    height2_plot.append(height2[t][maxima_number][direction])
                    
                    X3_plot.append(X3[t][maxima_number][direction])
                    rho3_plot.append(rho3[t][maxima_number][direction])
                    height3_plot.append(height3[t][maxima_number][direction])
                    
                    X4_plot.append(X4[t][maxima_number][direction])
                    rho4_plot.append(rho4[t][maxima_number][direction])
                    height4_plot.append(height4[t][maxima_number][direction])
                except IndexError:
                    break 
                    
            
            #.plot_list(t_points,X1_plot,f, "", 
            #                   "flow_time","X",legend1,1,1,"blue",True,True,3,4,direction+1)
            #Plotting.plot_list(t_points,X2_plot,f, "", 
            #                   "flow_time","X",legend2,1,1,"red",True,True,3,4,direction+1)
            #Plotting.plot_list(t_points,X3_plot,f, "", 
            #                   "flow_time","X",legend3,1,1,"green",True,True,3,4,direction+1)
            plt.legend(loc="lower left", prop={'size': 20})
            Plotting.plot_list(t_points,rho1_plot,f, "", 
                               "","",legend1,1,1,"blue",True,True,2,1,2)
            Plotting.plot_list(t_points,rho2_plot,f, "", 
                               "","",legend2,1,1,"red",True,True,2,1,2)
            Plotting.plot_list(t_points,rho4_plot,f, "", 
                               "","",legend4,1,1,"green",True,True,2,1,2)

            #Plotting.plot_list(t_points,height1_plot,f, "",
            #                   "flow_time","height",legend1,1,1,"blue",True,True,3,4,direction+1+8)
            #Plotting.plot_list(t_points,height2_plot,f, "",
            #                   "flow_time","height",legend2,1,1,"red",True,True,3,4,direction+1+8)
            #Plotting.plot_list(t_points,height3_plot,f, "",
            #                   "flow_time","height",legend3,1,1,"green",True,True,3,4,direction+1+8)
        #plt.savefig("Maxima"+str(maxima_number))
        
            
            #Plotting.plot_list(t_points,rho1_plot,f,"", 
            #                   "flow_time","rho",legend1,1,1,"blue",True,True,3,1,2)
            #Plotting.plot_list(t_points,rho2_plot,f, "", 
            #Plotting.plot_list(t_points,rho3_plot,f, "", 
            #                   "flow_time","rho",legend3,1,1,"green",True,True,3,1,2)

            plt.legend(loc="lower left", prop={'size': 20})
            Plotting.plot_list(t_points,height1_plot,f, "",
                               "","",legend1,1,1,"blue",True,True,2,1,1)
            Plotting.plot_list(t_points,height2_plot,f, "",
                               "","",legend2,1,1,"red",True,True,2,1,1)
            #Plotting.plot_list(t_points,height3_plot,f, "",
            #                   "flow_time","height",legend3,1,1,"green",True,True,3,1,2)
            Plotting.plot_list(t_points,height4_plot,f, "",
                               "","",legend4,1,1,"purple",True,True,2,1,1)
            plt.legend(loc="lower left", prop={'size': 20})
            
            #Plotting.plot_list(t_points,X1_plot,f, title, 
            #                   "flow_time","X_0",legend1,1,1,"blue",True,True,3,1,1)
            #Plotting.plot_list(t_points,X2_plot,f, title, 
            #                   "flow_time","X_0",legend2,1,1,"red",True,True,3,1,1)
           # Plotting.plot_list(t_points,X3_plot,f, title, 
           #                    "flow_time","X_0",legend3,1,1,"green",True,True,3,1,1)
            #Plotting.plot_list(t_points,X4_plot,f, title, 
             #                  "flow_time","X_0",legend4,1,1,"purple",True,True,3,1,1)
            plt.legend(loc="lower left", prop={'size': 20})
            plt.savefig(title+"_Maxima"+str(maxima_number)+".png",dpi=100)
            plt.show()
            
    count1=max_count(height1)
    count2=max_count(height2)
    #count3=max_count(height3)
    count4=max_count(height4)
    
    Plotting.plot_list(t_points,count1,f, "",
                               "flow_time","#maxima",legend1,1,1,"blue",True,True,1,1,1)
    Plotting.plot_list(t_points,count2,f, "",
                               "flow_time","#maxima",legend2,1,1,"red",True,True,1,1,1)
   # Plotting.plot_list(t_points,count3,f, "",
    #                           "flow_time","#maxima",legend3,1,1,"green",True,True,1,1,1)
    Plotting.plot_list(t_points,count4,f, "",
                               "flow_time","#maxima",legend4,1,1,"purple",True,True,1,1,1)
    plt.rcParams.update({'font.size': 12})
    lt.legend(loc="lower left", prop={'size': 12})
    plt.savefig("Maxima_evolution.pdf")
    
    return(height1, height2, height3, height4)

def max_count(height):
    count=np.zeros((np.shape(height)[0]))
    for i in range(0, np.shape(height)[0]):
            for j in range(0, np.shape(height)[1]):
                if (height[i][j][0]>1e-12):
                    count[i]+=1
    return(count)