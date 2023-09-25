import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import re

import Read
import Maxima_find
import os

matplotlib.rcParams.update({'font.size': 20})
plt.rcParams["figure.figsize"] = (7,5)


axis=list((map(int,input("Enter the axis to be summed over in ascending order: "))))



directory = os.fsencode("./")
Maxima_file=open('maxima.dat','w')

for file in os.listdir(directory):
     file_name = os.fsdecode(file)
     if file_name.endswith("to.dat"): 
        density,sizes,colors=Read.topology(file_name)
        density_sq=np.sqrt(density*density)

        density=density.sum(axis=(axis[0],axis[1]))
        sizes=np.delete(np.delete(sizes,axis[0]),axis[1]-1)
        res=Maxima_find.simple_2d(density,sizes)
        density_neg=np.array(-density)
        res_min=Maxima_find.simple_2d(density_neg,sizes)
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
        plt.xlabel("time")
        plt.ylabel("Q(t)")
        plt.plot(i,density)
        plt.savefig(file_name.replace('.dat','')+".png",dpi=150, bbox_inches='tight')
        plt.close()


        
        continue
    
     else:
         continue

Maxima_file.close()
