import matplotlib.pyplot as plt
import numpy as np
import tools
import sys
import tarfile
import os
import pandas as pd
import re
#------------Program to find instantons on a T²xR² lattice-----------------------------#
# Command to run: python3 ./main.py config.txt 
# There is a config.txt template example on the folder

directory=sys.argv[1]
tau=sys.argv[2]
neigh=1

beta=re.search(r"b\S{3}",directory).group().replace("b","")
print(beta)
nr=re.search(r"nr(\d+)",directory).group().replace("nr","")
print(nr)
nt=re.search(r"runns(.+?(?=nt))",directory).group().replace("nt","").replace("runns","")
print(nt)
f = open("identification"+str(nt)+"nt"+str(nr)+"nr"+str(tau)+"t"+str(beta)+"b.txt", "w")
delta_n=0
mean_frac=0
diff_frac=0
count=0

for filename in os.listdir(directory):
    if "to.dat" in filename and "dt"+str(tau) in filename:
        print(filename)
        count+=1
        conf=filename.replace("profile4dt"+str(tau)+"c", "")
        conf=conf.replace("to.dat", "")
        top_density,sizes=tools.read_top(directory+filename)
        density_2d_top,sizes_big,index_smal=tools.projection_2d(top_density,sizes)

        inst, a_inst, frac, a_frac, t_frac, t_inst, total= tools.find_inst_2d(density_2d_top,sizes_big,
                                                                              0,0,neigh)

        Q_top=density_2d_top.sum()
        mean_frac+= len(t_frac)
        diff_frac+= len(frac)-len(a_frac)
        Q_instantons= len(inst)-len(a_inst)+1/2*len(frac)-1/2*len(a_frac)

        f.write(str(conf)+" " + str(len(frac)) + " " +str(len(a_frac))+ " "+ 
                str(len(inst)) + " " +str(len(a_inst))+ " " +
        
        for element in t_frac:
            f.write( " "+ str(element[0])+ " "+str(element[1]))
            for fit in element[2]:
                f.write( " " +str(fit))
        f.write(" \n")
