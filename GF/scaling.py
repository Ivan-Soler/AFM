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

def plot_scaling(nr_list, nt_list, t_list, physical, physical_errors, ylabel, plotname, feature):
    for nr in nr_list:
        for nt in nt_list:
          for t in t_list:
            data_plot=[]
            data_error=[]
            x=[]
            plot=False
            for i in range(0,len(physical)):
                if physical[i][1]==nt and physical[i][2]==nr and physical[i][3]==t:
                    data_plot.append(physical[i][feature])
                    data_error.append(physical_errors[i][feature])
                    x.append(physical[i][0]+float(t)/100)
                    plot=True
            if plot:
                plt.errorbar(x,data_plot, yerr=np.transpose(data_error), label=str(nt)+", "+str(nr) +","+str(t), marker="o") 
    plt.legend(loc="upper left",ncol=2) 
    #plt.ylim(0.6,1.6)
    plt.xlabel("l_s(fm)")    
    plt.ylabel(ylabel)
    plt.savefig(plotname)
    plt.close()

    return 


def plot_scaling_flow(nr_list, nt_list, t_list, beta_list,physical, physical_errors, ylabel, plotname, feature):
    for nr in nr_list:
        for nt in nt_list:
          for b in beta_list:
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
                plt.errorbar(x,data_plot, yerr=np.transpose(data_error), label=str(nt)+", "+str(nr) +","+str(b), marker="o") 
    plt.legend(loc="upper left",ncol=2) 
    #plt.ylim(0.6,1.6)
    plt.xlabel("flow time")    
    plt.ylabel(ylabel)
    plt.savefig(plotname)
    plt.close()

    return 
  
nt_list=["4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"]
nr_list=["45", "64", "104"]
beta_list=["2.4","2.5","2.6","2.7"]
beta_list=["2.6"]
t_list=["10","15","25","50"]


#scaling_n3=pickle.load(open('scaling_'+str(n)+'.pkl', "rb" ))
#scaling_n8=pickle.load(open('scaling_'+str(n)+'.pkl', "rb" ))
scaling_min=pickle.load(open('scaling_s0.8.pkl', "rb" ))
scaling_med=pickle.load(open('scaling_s0.85.pkl', "rb" ))
scaling_max=pickle.load(open('scaling_s0.9.pkl', "rb" ))

#print(n)

table_ensembles=copy.deepcopy(scaling_med)
#print(table_ensembles)

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
                       table_ensembles[key]["means"][1]])#norm

        physical_errors.append([ls,nt,nr,t,beta,table_ensembles[key]["errors"][0]/table_ensembles[key]["vol"], 
                  table_ensembles[key]["errors"][1]/np.pi/(table_ensembles[key]["means"][2]*table_ensembles[key]["a"])**2, #height_fit=norm/(pi*rho**2)
                  table_ensembles[key]["errors"][2]*table_ensembles[key]["a"], #rho
                  table_ensembles[key]["errors"][3]/(table_ensembles[key]["a"]**2), #height
                              table_ensembles[key]["errors"][1]])#norm

for element in physical:
  print(element[3], element[5])
#Sorting according to ls
physical = sorted(physical, key=lambda a_entry: a_entry[0]) 
physical_errors = sorted(physical_errors, key=lambda a_entry: a_entry[0]) 

ylabel="density(1/fm^2)"
plotname="scaling_dens.pdf"
plot_scaling(nr_list, nt_list, t_list, physical, physical_errors, ylabel, plotname, 5)  

ylabel="density(1/fm^2)"
plotname="scaling_flow.pdf"
plot_scaling_flow(nr_list, nt_list, t_list, beta_list, physical, physical_errors, ylabel, plotname, 5)  

ylabel="height(1/fm^2)"
plotname="scaling_height_fit.png"
plot_scaling(nr_list, nt_list, t_list, physical, physical_errors, ylabel, plotname, 6) 

ylabel="density(1/fm^2)"
plotname="scaling_height_flow.pdf"
plot_scaling_flow(nr_list, nt_list, t_list, beta_list,physical, physical_errors, ylabel, plotname, 6)  

ylabel="width(1/fm^2)"
plotname="scaling_rho.png"
plot_scaling(nr_list, nt_list, t_list, physical, physical_errors, ylabel, plotname, 7) 

ylabel="density(1/fm^2)"
plotname="scaling_rho_flow.pdf"
plot_scaling_flow(nr_list, nt_list, t_list, beta_list,physical, physical_errors, ylabel, plotname, 7)  

ylabel="norm"
plotname="scaling_norm.png"
plot_scaling(nr_list, nt_list, t_list, physical, physical_errors, ylabel, plotname, 8)   
