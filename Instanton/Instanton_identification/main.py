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

f = open("data_inst.txt", "w")
f.write("Config \t Frac \t Anti_frac \t Inst \t Anti_Inst \t Q_top \t Q_inst  \n \n")

for i in range(start,end,step):
    file=directory+config_template+str(i)+"to.dat"
    print(file)
    
    density,sizes=tools.read_top(file)
    
    density_2d,sizes_big,index_smal=tools.projection_2d(density,sizes)
    
    inst, a_inst, frac, a_frac, t_frac, t_inst, total=tools.find_max_2d(density_2d,sizes_big,size_frac,size_inst)
    
    Q_top=density_2d.sum()
    Q_instantons=len(inst)-len(a_inst)+1/2*len(frac)-1/2*len(a_frac)
    
    f.write(str(i)+"\t \t" + str(len(frac)) + "\t" +str(len(a_frac))+ 
            "\t" +str(len(inst)) + "\t" +str(len(a_inst))+ 
            "\t"+str(Q_top)+"\t"+str(Q_instantons)+"\n")
    
    tools.plot_dens_2d(file,density_2d,sizes_big, t_frac)
    
f.close()
