import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import Read
import Maxima_find
import Plotting
import copy
plt.rcParams.update({'font.size': 12})

def inst(density, size, Xmax, plot):
    
    def func(x, Xmax_extr, rho, height):
        return (height + (x-Xmax_extr)*(x-Xmax_extr)/(rho*rho))
    
    x=np.array([Xmax-1,Xmax,Xmax+1])
    
    data=density[x]**(-2/5)
    #try:
    popt, pcov= curve_fit(func, x, data)
        
    #except RuntimeError:
    #    print('error')
    #    return 100, 100, 100
    rho=popt[1]
    Xmax_extr=popt[0]
    print(popt)
    
    #Normalization
    N=density[Xmax]**(-2/5)/func(Xmax, *popt)
    #xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
    #plt.plot(xdata, func(xdata, *popt))
    #plt.plot(x, data)
    #plt.show()
    
    if plot:
        x=np.array([(Xmax-2),(Xmax-1), Xmax,(Xmax+1), (Xmax+2)])
        data=density[x]     
        xdata = np.linspace(Xmax-2, Xmax+2, 50)
        plt.plot(xdata, (N*func(xdata, *popt))**(-5/2))
        plt.plot(x, data)
        #plt.scatter(x, (N*func(x, *popt))**(-5/2))
        plt.show()
        
    return popt

def par(density, size, Xmax, plot):
    
    def parabola(x, Xmax_extr, a, c):
        return (a*(x-Xmax_extr)*(x-Xmax_extr)+c)
    
    x=np.array([Xmax-1,Xmax,Xmax+1])
    data=density[x]
    #try:
    popt, pcov = curve_fit(parabola, x, data, bounds=([Xmax-1, -10000,-1], [Xmax+1, 0,1]))
        
    heihgt=popt[2]
    #popt[1]=np.sqrt(popt[1])
    Xmax_extr=popt[0]

    #xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
    #plt.plot(xdata, func(xdata, *popt))
    #plt.plot(x, data)
    #plt.show()
    
    if plot:
        x=np.array([(Xmax-2),(Xmax-1), Xmax,(Xmax+1), (Xmax+2)])
        data=density[x]      
        xdata = np.linspace(Xmax-2, Xmax+2, 50)
        plt.plot(xdata, parabola(xdata, *popt))
        plt.scatter(x, data)
        #plt.scatter(x, (N*func(x, *popt))**(-5/2))
        plt.show()
        
    return popt

