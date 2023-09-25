import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import re

import Read
import Maxima_find
import os


matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rcParams.update({'font.size': 11})
matplotlib.rcParams["figure.figsize"] = (8/2.54,6/2.54)
matplotlib.rc('text', usetex=True)

axis=list((map(int,input("Enter the axis to be summed over in ascending order: "))))


directory = os.fsencode("./")
Maxima_file=open('maxima.dat','w')

for file in os.listdir(directory):
    file_name = os.fsdecode(file)
    if file_name.endswith("to.dat"): 
        density,sizes=Read.topology(file_name)
        density_sq=np.sqrt(density*density)

        density=density.sum(axis=axis[0]).sum(axis=(axis[1]-1)).sum(axis=(axis[2]-2))
        res=Maxima_find.simple_1d(density,len(density))
        density_neg=np.array(-density)
        res_min=Maxima_find.simple_1d(density_neg,len(density))
        print(len(density))
        print("Maxima density:\t" + str(np.size(res)))
        print("Minima density:\t" + str(np.size(res_min)))

        #Getting time and conf
        match = re.search(r'dt\d*.*\d*c', file_name)
        if match:
             time= match.group().replace("dt","").replace("c","")
        match = re.search(r'c\d+to', file_name)    
        if match:
             conf= match.group().replace("to","").replace("c","")

        print(conf+"\t"+time+"\t"+str(np.size(res))+"\t"+str(np.size(res_min)),file=Maxima_file)
          
        i=np.arange(0,len(density))
        plt.xlabel(r'Time')
        plt.ylabel(r'q(t)')
        #plt.xlim(0,128)
        #plt.xticks(np.arange(0,128,20))
        plt.plot(i,density)
        plt.savefig(file_name.replace('.dat','')+".png",dpi=150, bbox_inches='tight')
        plt.close()


        
        continue
    
    else:
        continue

Maxima_file.close()
