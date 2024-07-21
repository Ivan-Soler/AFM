import os
import matplotlib.pyplot as plt
import sys
import re
import numpy as np
import tools
import pandas as pd
from collections import OrderedDict
import pickle
import copy

def nearest_distance(fractionals):
  distance=[]
  for i in range(len(fractionals)):
    distance_min=100000
    for j in range(len(fractionals)):
      if i!=j:
        distance_tmp=np.sqrt((fractionals[i][0]-fractionals[j][0])**2+(fractionals[i][1]-fractionals[j][1])**2)
        if distance_tmp<distance_min:
          distance_min=distance_tmp
    distance.append(distance_min)

  distance=np.array((distance))
  mean_distance=np.mean(distance)
  return mean_distance
      
def mean_distance(table_ensembles):
  
  for key in table_ensembles:
    distance=[]
    for conf in table_ensembles[key]['pos_afrac']:
      if table_ensembles[key]['pos_afrac'][conf] and table_ensembles[key]['pos_frac'][conf]:
        fractionals=np.vstack((table_ensembles[key]['pos_frac'][conf],table_ensembles[key]['pos_afrac'][conf]))
        distance.append(nearest_distance(fractionals))
        
        #distance.append(nearest_distance(table_ensembles[key]['pos_frac'][conf]))
        #distance.append(nearest_distance(table_ensembles[key]['pos_afrac'][conf]))
    distance=np.array((distance))
    mean_distance=np.mean(distance)
    err_distance=np.std(distance)/np.sqrt(len(distance))
    table_ensembles[key]["means"]=np.hstack((table_ensembles[key]["means"],mean_distance))
    table_ensembles[key]["errors"]=np.hstack((table_ensembles[key]["errors"],0))
    #print(key, mean_distance)
    
  return table_ensembles
        
def rho_smooth(a,b,x):
  return(a+b*x)
  
def plot_scaling_ls(nr_list, nt_list, t_list, beta_list, physical, physical_errors, ylabel, plotname, feature):

    for nr in nr_list:
        for nt in nt_list:
          #for t in t_list:
            data_plot=[]
            data_error=[]
            x=[]
            plot=False
            for i in range(0,len(physical)):
                if physical[i][1]==nt and physical[i][2]==nr and float(physical[i][0])<0.8:# and physical[i][3]==t
                    data_plot.append(physical[i][feature])
                    data_error.append(physical_errors[i][feature])
                    x.append(physical[i][0])
                    plot=True
            if plot:
                plt.errorbar(x,data_plot, yerr=np.transpose(data_error),marker="o", label=str(nt)+", "+str(nr))


    plt.legend(loc="lower right",ncol=2)
    #plt.ylim(0.6,1.6)
    plt.xlabel('$l_s$(fm)')
    plt.ylabel(ylabel)
    plt.savefig(plotname,dpi=800, bbox_inches='tight')
    plt.close()

    return 
  
def plot_scaling_N_even_odd(nr_list, nt_list, t_list, beta_list, physical, physical_errors, ylabel, plotname, feature,beta):

    
    data_plot=[]
    data_error=[]
    x=[]

    data_plot_odd=[]
    data_error_odd=[]
    x_odd=[]
  
    for nr in nr_list:
      #for b in beta_list:
        for nt in nt_list:
          #for t in t_list:
            #data_plot=[]
            #data_plot_odd=[]
            #data_error=[]
            #x=[]
            plot=False
            for i in range(0,len(physical)):
                if int(physical[i][1])%2 and physical[i][2]==nr and physical[i][3]==t and physical[i][4]==beta:
                    a=float(physical[i][0])/float(physical[i][1])
                    data_plot.append(physical[i][feature])
                    data_error.append(physical_errors[i][feature])
                    x.append(int(physical[i][1]))#+float(t)*0.05)
                    #x.append(physical[i][0])#+float(t)*0.05)
                    plot=True
                  
                if not int(physical[i][1])%2 and physical[i][2]==nr and physical[i][3]==t and physical[i][4]==beta:
                    a=float(physical[i][0])/float(physical[i][1])
                    data_plot_odd.append(physical[i][feature])
                    data_error_odd.append(physical_errors[i][feature])
                    x_odd.append(int(physical[i][1]))#+float(t)*0.05)
                    #x_odd.append(physical[i][0])#+float(t)*0.05)
                    plot=True
              #if plot:
                #plt.errorbar(x,data_plot, yerr=np.transpose(data_error), marker="o")#, label=str(nt)+", "+str(nr)+", "+str(t)) 
    xi_even=np.argsort(x)
    xi_odd=np.argsort(x_odd)

    x=np.array((x))[xi_even]
    x_odd=np.array((x_odd))[xi_odd]

    data_plot=np.array((data_plot))[xi_even]
    data_plot_odd=np.array((data_plot_odd))[xi_odd]

    data_error_odd=np.array((data_error_odd))[xi_odd]
    data_error=np.array((data_error))[xi_even]
    
    plt.errorbar(x,data_plot, yerr=np.transpose(data_error), marker="o", label="even")
    plt.errorbar(x_odd,data_plot_odd, yerr=np.transpose(data_error_odd), marker="o", label="odd")

    #x = np.linspace(0, 20, 1000)
    #y=rho_smooth(-0.243684,0.44703,x)
    #plt.plot(x,y)#, label="Smooth", color="black")
    plt.legend(loc="lower right",ncol=2) 
    #plt.ylim(0.6,1.6)  
    plt.xlabel("ls")  
    plt.ylabel(ylabel)
    plt.savefig(plotname,dpi=800, bbox_inches='tight')
    plt.close()

    return 

def plot_scaling_N(nr_list, nt_list, t_list, beta_list, physical, physical_errors, ylabel, plotname, feature,beta):



    for nr in nr_list:
        for nt in nt_list:
          #for t in t_list:
            data_plot=[]
            data_error=[]
            x=[]
            plot=False
            for i in range(0,len(physical)):
                if physical[i][1]==nt and physical[i][2]==nr and float(physical[i][0])<0.8:# and physical[i][3]==t
                    a=float(physical[i][0])/float(physical[i][1])
                    data_plot.append(physical[i][feature]/a)
                    data_error.append(physical_errors[i][feature]/a)
                    x.append(int(physical[i][1]))
                    plot=True
            if plot:
                plt.errorbar(x,data_plot, yerr=np.transpose(data_error),marker="o", label=str(nt)+", "+str(nr))

    x = np.linspace(2, 14, 1000)
    y=rho_smooth(-0.243684,0.44703,x)
    plt.plot(x,y, label="Smooth", color="black",linestyle='dashed')
    plt.legend(loc="upper left",ncol=2)
    plt.ylim(0.5,7)
    plt.xlabel('$N_s$')
    plt.ylabel(ylabel)
    plt.savefig(plotname,dpi=800, bbox_inches='tight')
    plt.close()

    return

  
def plot_scaling_flow(nr_list, nt_list, t_list, beta_list, physical, physical_errors, ylabel, plotname, feature):
    for nr in nr_list:
        for nt in nt_list:
          for b in beta_list:
            for t in t_list:
              data_plot=[]
              data_error=[]
              x=[]
              plot=False
              for i in range(0,len(physical)):
                if physical[i][1]==nt and physical[i][2]==nr and physical[i][4]==b:
                    data_plot.append(physical[i][feature])
                    data_error.append(physical_errors[i][feature])
                    x.append(float(physical[i][3]))
                    plot=True
            if plot:
                plt.errorbar(x,data_plot, yerr=np.transpose(data_error), marker="o")#,label=str(nt)+", "+str(nr) +","+str(b)) 
    #plt.legend(loc="upper left",ncol=2) 
    #plt.ylim(0.6,1.6)
    plt.xlabel("flow time")    
    plt.ylabel(ylabel)
    plt.savefig(plotname)
    plt.close()

    return 
  
nt_list=["4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"]
nr_list=["32", "45", "104"]
nr_list=["104", "64","45"]
t_list=["10","15","25","35","50"]
beta_list=["2.4","2.5","2.6", "2.7"]
t_list=["3.52", "7.34", "15", "27.37", "47.84"]
#scaling_n3=pickle.load(open('scaling_'+str(n)+'.pkl', "rb" ))
#scaling_n8=pickle.load(open('scaling_'+str(n)+'.pkl', "rb" ))
scaling_min=pickle.load(open('scaling_s0.7.pkl', "rb" ))
scaling_med=pickle.load(open('scaling_s0.8.pkl', "rb" ))
scaling_max=pickle.load(open('scaling_s0.9.pkl', "rb" ))

for key in scaling_med:
  print(key, len(scaling_med[key]["pos_afrac"]))
  print(key, len(scaling_min[key]["pos_afrac"]))
  print(key, len(scaling_max[key]["pos_afrac"]))
  
scaling_min=mean_distance(scaling_min)
scaling_med=mean_distance(scaling_med)
scaling_max=mean_distance(scaling_max)

table_ensembles=copy.deepcopy(scaling_med)

for key in scaling_min:  
    media=np.array((scaling_med[key]["means"]))
    table_ensembles[key]["means"]=np.array((media)) 
    higher=np.maximum(scaling_min[key]["means"],scaling_max[key]["means"])
    lower=np.minimum(scaling_min[key]["means"],scaling_max[key]["means"])

    errors=np.stack((lower,higher))
    errors=np.sort(errors,axis=0)
    errors=np.transpose(errors)
    errors=np.array((np.abs(media-errors[:,0]),np.abs(errors[:,1]-media)))
    errors=np.transpose(errors)

    table_ensembles[key]["errors"]=np.dstack((table_ensembles[key]["errors"],table_ensembles[key]["errors"]))[0]+np.array((errors))
    
physical=[]
physical_errors=[]

for key in table_ensembles:
    if table_ensembles[key]["configurations"]>1 and float(table_ensembles[key]["beta"])<2.75 :
        ls=table_ensembles[key]["ls"]
        nr=table_ensembles[key]["nr"]
        nt=table_ensembles[key]["nt"]
        t=table_ensembles[key]["t"]
        beta=table_ensembles[key]["beta"]
        physical.append([ls,nt,nr,t,beta,table_ensembles[key]["means"][0]/table_ensembles[key]["vol"], #density
                  abs(table_ensembles[key]["means"][1])/np.pi/(table_ensembles[key]["means"][2]*table_ensembles[key]["a"])**2, #height_fit=norm/(pi*rho**2)
                  table_ensembles[key]["means"][2]*table_ensembles[key]["a"], #rho
                  abs(table_ensembles[key]["means"][3])/(table_ensembles[key]["a"]**2), #height
                       table_ensembles[key]["means"][1]/(table_ensembles[key]["a"]**2),#norm
                        table_ensembles[key]["means"][4]*table_ensembles[key]["a"]]) #distance

        physical_errors.append([ls,nt,nr,t,beta,table_ensembles[key]["errors"][0]/table_ensembles[key]["vol"], 
                  table_ensembles[key]["errors"][1]/np.pi/(table_ensembles[key]["means"][2]*table_ensembles[key]["a"])**2, #height_fit=norm/(pi*rho**2)
                  table_ensembles[key]["errors"][2]*table_ensembles[key]["a"], #rho
                  table_ensembles[key]["errors"][3]/(table_ensembles[key]["a"]**2), #height
                              table_ensembles[key]["errors"][1]/(table_ensembles[key]["a"]**2),#norm
                               table_ensembles[key]["errors"][4]*table_ensembles[key]["a"]]) #distance

#Sorting according to ls
physical = sorted(physical, key=lambda a_entry: a_entry[0]) 
physical_errors = sorted(physical_errors, key=lambda a_entry: a_entry[0]) 

beta=""

plt.rcParams.update({'font.size': 14})

ylabel='$n_{fi}(1/fm^2)$'
plotname="scaling_dens.pdf"
plot_scaling_ls(nr_list, nt_list, t_list, beta_list,physical, physical_errors, ylabel, plotname, 5)  

ylabel='density$(1/fm^2)$'
plotname="scaling_dens_Nb"+beta+".pdf"
plot_scaling_N(nr_list, nt_list, t_list, beta_list,physical, physical_errors, ylabel, plotname, 5,beta)  

ylabel='height$(fm^2)$'
plotname="scaling_height_fit.pdf"
plot_scaling_ls(nr_list, nt_list, t_list, beta_list,physical, physical_errors, ylabel, plotname, 8)   

ylabel='$\\rho_{fi}(1/fm^2)$'
plotname="scaling_rho.pdf"
plot_scaling_ls(nr_list, nt_list,t_list, beta_list, physical, physical_errors, ylabel, plotname, 7)  

ylabel='$\\rho_{fi}$'
plotname="scaling_rho_Nb"+beta+".pdf"
plot_scaling_N(nr_list, nt_list, t_list, beta_list,physical, physical_errors, ylabel, plotname, 7,beta)   

ylabel="width"
plotname="scaling_rho_flow.pdf"
plot_scaling_flow(nr_list, nt_list, t_list, beta_list, physical, physical_errors, ylabel, plotname, 7)  

ylabel="norm(fm^2)"
plotname="scaling_norm.pdf"
plot_scaling_ls(nr_list, nt_list, t_list, beta_list,physical, physical_errors, ylabel, plotname, 9)   

ylabel="distance(fm)"
plotname="scaling_distance.pdf"
plot_scaling_ls(nr_list, nt_list, t_list, beta_list,physical, physical_errors, ylabel, plotname, 10)   




