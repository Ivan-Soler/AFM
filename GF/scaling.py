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

def plot_scaling(nr_list, nt_list, physical, physical_errors, ylabel, plotname, feature):
    for nr in nr_list:
        for nt in nt_list:
            data_plot=[]
            data_error=[]
            x=[]
            plot=False
            for i in range(0,len(physical)):
                if physical[i][1]==nt and physical[i][2]==nr:
                    data_plot.append(physical[i][feature])
                    data_error.append(physical_errors[i][feature])
                    x.append(physical[i][0])
                    plot=True
            if plot:
                plt.errorbar(x,data_plot, yerr=np.transpose(data_error), label=str(nt)+", "+str(nr), marker="o") 
    plt.legend(loc="upper left",ncol=2) 
    #plt.ylim(0.6,1.6)
    plt.xlabel("l_s(fm)")    
    plt.ylabel(ylabel)
    plt.savefig(plotname)
    plt.close()

    return 

nt_list=["4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"]
nr_list=["32", "45", "104"]
nr_list=["64","45", "104"]
#scaling_n3=pickle.load(open('scaling_'+str(n)+'.pkl', "rb" ))
#scaling_n8=pickle.load(open('scaling_'+str(n)+'.pkl', "rb" ))
scaling_min=pickle.load(open('scaling_s0.6.pkl', "rb" ))
scaling_med=pickle.load(open('scaling_s0.7.pkl', "rb" ))
scaling_max=pickle.load(open('scaling_s0.8.pkl', "rb" ))

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
print( table_ensembles[key]["errors"])
physical=[]
physical_errors=[]
for key in table_ensembles:
    if table_ensembles[key]["configurations"]>1 and float(table_ensembles[key]["beta"])<2.75 :
        ls=table_ensembles[key]["ls"]
        nr=table_ensembles[key]["nr"]
        nt=table_ensembles[key]["nt"]
        physical.append([ls,nt,nr,table_ensembles[key]["means"][0]/table_ensembles[key]["vol"], #density
                  abs(table_ensembles[key]["means"][1])/np.pi/(table_ensembles[key]["means"][2]*table_ensembles[key]["a"])**2, #height_fit=norm/(pi*rho**2)
                  table_ensembles[key]["means"][2]*table_ensembles[key]["a"], #rho
                  abs(table_ensembles[key]["means"][3])/(table_ensembles[key]["a"]**2), #height
                       table_ensembles[key]["means"][1]])#norm

        physical_errors.append([ls,nt,nr,table_ensembles[key]["errors"][0]/table_ensembles[key]["vol"], 
                  table_ensembles[key]["errors"][1]/np.pi/(table_ensembles[key]["means"][2]*table_ensembles[key]["a"])**2, #height_fit=norm/(pi*rho**2)
                  table_ensembles[key]["errors"][2]*table_ensembles[key]["a"], #rho
                  table_ensembles[key]["errors"][3]/(table_ensembles[key]["a"]**2), #height
                              table_ensembles[key]["errors"][1]])#norm

#Sorting according to ls
physical = sorted(physical, key=lambda a_entry: a_entry[0]) 
physical_errors = sorted(physical_errors, key=lambda a_entry: a_entry[0]) 

ylabel="density(1/fm^2)"
plotname="scaling_dens.pdf"
plot_scaling(nr_list, nt_list, physical, physical_errors, ylabel, plotname, 3)  

ylabel="height(1/fm^2)"
plotname="scaling_height_fit.png"
plot_scaling(nr_list, nt_list, physical, physical_errors, ylabel, plotname, 4)   

ylabel="width(1/fm^2)"
plotname="scaling_rho.png"
plot_scaling(nr_list, nt_list, physical, physical_errors, ylabel, plotname, 5)   

ylabel="norm"
plotname="scaling_norm.png"
plot_scaling(nr_list, nt_list, physical, physical_errors, ylabel, plotname, 7)   
