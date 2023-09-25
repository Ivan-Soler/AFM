import struct
import numpy as np
import sys

import Projecting
import Plotting
import Max_finder

import os

def remaining_axes(dir_1,dir_2):
    if (dir_1==3):
        if (dir_2==2):
            return 'x' , 'y'
        if (dir_2==1):
            return 'x' , 'z'
        if (dir_2==0):
            return 'y' , 'z'
        
    if (dir_1==2):
        if (dir_2==1):
            return 'x' , 't'
        if (dir_2==0):
            return 'y' , 't'
        
    if (dir_1==1):
        if (dir_2==0):
            return 'z' , 't'

    print("direction 1 must be bigger than direction 2")
    return 0

file_name_1=str(sys.argv[1])
#plot_title_1="GF_top"
plot_title_1="AFM"
if (len(sys.argv)>2):
	number_maxima=sys.argv[2]
else:
	number_maxima=8 
if (len(sys.argv)>3):
    sizes=[sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6]]
    
if (len(sys.argv)>7):    
    file_name_2=str(sys.argv[7])


f = open(file_name_1, "rb")
for i in range (0,4):
    sizes[i]=int(f.readline())
precision=int(f.readline())
f.close()
print(sizes)
print(precision)



#density=Projecting.D4(file_name_1, sizes, precision)

#Plotting.plot_int_3d(plot_title_1, density, "x", "y", sizes)
#Plotting.plot_nint_3d(plot_title_1, density, "x", "y", sizes)

            
