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
#eps_action=eps_norm

print(action)
f = open(directory+"data_inst_"+directory.replace("../runns","").replace("/","").replace("nr104","")+"t"+
         config_template.replace("profile4dt","").replace("c","")+".txt", "w")
f.write("Config \t Frac \t Anti_frac \t Inst \t Anti_Inst \t delta_q \t Tot_fract \n \n")
#f.write("Config \t Frac \t Anti_frac \t Tot_fract \n \n")
delta_n=0
mean_frac=0
diff_frac=0
count=0
for i in range(start,end,step):
    count+=1
    file_top=directory+config_template+str(i)+"to.dat"    
    top_density,sizes=tools.read_top(file_top)
    density_2d_top,sizes_big,index_smal=tools.projection_2d(top_density,sizes)
    
    if action:
        file_act=directory+config_template+str(i)+"en.dat"
        act_density,sizes=tools.read_top(file_act)
        #act_density=-act_density
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
    mean_frac+=len(t_frac)
    diff_frac+= len(frac)-len(a_frac)
    Q_instantons=len(inst)-len(a_inst)+1/2*len(frac)-1/2*len(a_frac)
    
    if (len(frac)+len(a_frac)) % 2:
        delta_n+=1
        
    #maxima=tools.find_max_2d(density_2d_top,sizes_big)
    f.write(str(i)+"\t \t" + str(len(frac)) + "\t" +str(len(a_frac))+ "\t"+ 
            str(len(inst)) + "\t" +str(len(a_inst))+ "\t" +
            str(Q_top-Q_instantons) + "\t"+str(len(frac)+len(a_frac))+"\n")
    
    tools.plot_dens_2d(file_top,density_2d_top,sizes_big, t_frac, t_inst)
    tools.plot_dens_2d(file_act,-density_2d_act,sizes_big, t_frac, t_inst)

f.write("Configuration with fractional topological charge: " + str(delta_n) + "\n")
f.write("Mean frac_inst + frac_anti_inst = " + str(mean_frac/count)+ "\n")
f.write("Mean frac_inst - frac_anti_inst = " + str(diff_frac/count)+ "\n")
f.close()
print(count)
print(t_frac)
print(str(delta_n))



