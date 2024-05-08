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
rho=param[5]
norm=param[6]
eps_rho=param[7]
eps_norm=param[8]
neigh=param[9]
action=param[10]
eps_action=param[11]

print(action)
f = open(directory+"data_inst.txt", "w")
#f.write("Config \t Frac \t Anti_frac \t Inst \t Anti_Inst \t Q_top \t Q_inst \t Tot_fract \n \n")
f.write("Config \t Frac \t Anti_frac \t Tot_fract \n \n")
delta_n=0
for i in range(start,end,step):
    file_top=directory+config_template+str(i)+"to.dat"
    print(file_top)      
    top_density,sizes=tools.read_top(file_top)
    density_2d_top,sizes_big,index_smal=tools.projection_2d(top_density,sizes)
    
    if action:
        file_act=directory+config_template+str(i)+"en.dat"
        act_density=tools.read_top(file_act)
        density_2d_act,sizes_big,index_smal=tools.projection_2d(act_density,sizes)
        norm_action=norm*8*np.pi*np.pi
    else:
        density_2d_act=np.array(density_2d_top)
        eps_action=eps_norm
        norm_action=norm

    #print(rho)
    inst, a_inst, frac, a_frac, t_frac, t_inst, total=    tools.find_inst_2d(density_2d_top,density_2d_act,sizes_big,
                       rho,norm,eps_rho,eps_norm,norm_action,eps_action,neigh)
    
    Q_top=density_2d_top.sum()
    Q_instantons=len(inst)-len(a_inst)+1/2*len(frac)-1/2*len(a_frac)
    
    if (len(frac)+len(a_frac)) % 2:
        delta_n+=1
        
    
    f.write(str(i)+"\t \t" + str(len(frac)) + "\t" +str(len(a_frac))+ "\t"+ str(len(frac)+len(a_frac))+"\n")
    
    tools.plot_dens_2d(file_top,density_2d_top,sizes_big, t_frac)

f.write("Configuration with fractional topological charge: " + str(delta_n))
f.close()
print(str(delta_n))



