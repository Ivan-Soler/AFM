mport matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d, Axes3D
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interpn
import Projecting
import Plotting
import Read
import Maxima_find
import Fitting
import subprocess
import math
import re
from subprocess import call

plt.rcParams.update({'font.size': 12})


def fit_impr(folder,file,zm_top,normalization):
    
    if zm_top:
        density,sizes=Read.ascii_mode(file)
        title="susy"
    else:
        density,sizes=Read.topology(file)
        title="top"
    res=Maxima_find.simple(density,sizes)
    density=density*4*np.pi*np.pi*len(res[0])*normalization
    Plotting.plot_all_peaks(folder,density, sizes, title)

    param_all=[[]]
    param_all.pop(0)
    param_temp=[[]]
    param_temp.pop(0)

    err_all=[[]]
    err_all.pop(0)
    err_temp=[[]]
    err_temp.pop(0)
    for i in range(0,len(res[0])):
        param_temp=[[]]
        param_temp.pop(0)
        for d in range(0,4):
            Xmax=[res[0,i],res[1,i],res[2,i],res[3,i]]
            popt_temp,pcov_temp=Fitting.fitting_instanton_impr(density, d, sizes, Xmax, "0") #h,x,rho
            param_temp.append(popt_temp)
            err_temp.append(np.sqrt(np.diag(pcov_temp)))
        param_all.append(param_temp)
        err_all.append(err_temp)

    param=np.array(param_all)
    err=np.array(err_all)
 
    return density, param, err 

def fit_simple(folder,file,zm_top,normalization):
    if zm_top:
        density,sizes=Read.ascii_mode(file)
        res=Maxima_find.simple(density,sizes)
        density=density*4*np.pi*np.pi*len(res[0])
        title="susy"
    else:
        density,sizes=Read.topology(file)
        res=Maxima_find.simple(density,sizes)
        title="top"
    
    Plotting.plot_all_peaks(folder,density, sizes, title)

    param_all=[[]]
    param_all.pop(0)
    param_temp=[[]]
    param_temp.pop(0)

    for i in range(0,len(res[0])):
        param_temp=[[]]
        param_temp.pop(0)
        for d in range(0,4):
            Xmax=[res[0,i],res[1,i],res[2,i],res[3,i]]
            popt_temp=Fitting.fitting_instanton(density, d, sizes, Xmax, "0")
            param_temp.append(popt_temp)
        param_all.append(param_temp)
    param=np.array(param_all)
    return density, param