import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d, Axes3D
from scipy import interpolate
from scipy.ndimage.filters import maximum_filter
from scipy.interpolate import RegularGridInterpolator

def d4(data,sizes):
    print('creting grid extrapolated')
    X = np.arange(0, sizes[0], 1)
    Y = np.arange(0, sizes[1], 1)
    Z = np.arange(0, sizes[2], 1)
    T = np.arange(0, sizes[3], 1)

    interp = RegularGridInterpolator((X, Y, Z, T), data)

    Xg, Yg, Zg, Tg = np.meshgrid(X, Y, Z, T, indexing='ij')
    points_big=Tbig, Zbig, Ybig, Xbig = np.mgrid[0:sizes[3]-1:71j, 0:sizes[2]-1:71j, 0:sizes[1]-1:71j, 0:sizes[0]-1:71j]
    
    data_int=interp((Xbig,Ybig,Zbig,Tbig))
    print('searching maxima')
    data_max=maximum_filter(data_int, 3, mode = 'constant', cval=0.0) 
    maxima=(data_int==data_max)
    res = np.where(data_int == data_max)

    print("The total maxima are located at ")
    print(res)
    return(res, data_int)

def d4_backup(data, sizes):
  
    X = np.arange(0, sizes[0], 1)
    Y = np.arange(0, sizes[1], 1)
    Z = np.arange(0, sizes[2], 1)
    X, Y, Z = np.meshgrid(X, Y, Z, indexing='ij')

    Zbig, Ybig, Xbig = np.mgrid[0:sizes[2]-1:71j, 0:sizes[1]-1:71j, 0:sizes[0]-1:71j]
    tck = interpolate.bisplrep(Z, Y, X, data, s=0)
    #data_int= interpolate.bisplev(Zbig[:,0], Ybig[:,0], Xbig[0,:], tck)
    print(interpn(points, values, data_int))
    
    data_max=maximum_filter(data_int, 8, mode = 'constant', cval=0.0) 
    maxima=(data_int==data_max)
    res = np.where(data_int == data_max)
    
    print("The total maxima is located at ")
    print(res[1]/10)
    #print(data_max)
