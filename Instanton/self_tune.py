import matplotlib.pyplot as plt
import numpy as np
import tools
import sys


#------------Program to find instantons on a T²xR² lattice-----------------------------#
# Command to run: python3 ./main.py config.txt 
# There is a config.txt template example on the folder

config=sys.argv[1]
param=tools.parse_in_file(config)

directory=param[0]
config_template=param[1]
start=param[2]
end=param[3]
step=param[4]
size_frac=param[5]
size_inst=param[6]
size_local=param[7]

f = open("tuning.txt", "w")
f.write("size frac \t size inst \t delta Q \n")

for s in range(int(100*size_frac),int(100*size_frac)+10,1):
    for ss in range(int(100*size_inst),int(100*size_inst)+10,1):
        delta_q=0
        for i in range(start,end,step):
            file=directory+config_template+str(i)+"to.dat"
            #print(file)

            density,sizes=tools.read_top(file)

            density_2d,sizes_big,index_smal=tools.projection_2d(density,sizes)

            inst, a_inst, frac, a_frac, t_frac, t_inst, total=tools.find_max_2d(density_2d,sizes_big,s/100,ss/100,size_local)

            Q_top=density_2d.sum()
            Q_instantons=len(inst)-len(a_inst)+1/2*len(frac)-1/2*len(a_frac)
            
            delta_q+=abs(round(Q_top)-Q_instantons)

        f.write(str(s)+"\t "+ str(ss) + "\t" +str(delta_q)+"\n")
        print(s,ss)

f.close()


