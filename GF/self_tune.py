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
norm_action=norm*8*np.pi*np.pi

eps_rho=param[7]
eps_norm=param[8]
neigh=param[9]
action=param[10]
eps_action=param[11]

print(action)

f = open(directory+"tuning.txt", "w")
f.write("eps rho \t eps norm  \t eps action \t delta n \t delta q \t n \n")

for s in range(int(100*eps_rho),int(100*eps_rho)+1,10):
    for ss in range(int(100*eps_norm),int(100*eps_norm)+1,10):
        for sss in range(int(1000*eps_action),int(1000*eps_action)+11,1):

            n=0
            delta_q=0
            delta_n=0
            for i in range(start,end,step):
                file_top=directory+config_template+str(i)+"to.dat"
                #print(file_top)      
                top_density,sizes=tools.read_top(file_top)
                density_2d_top,sizes_big,index_smal=tools.projection_2d(top_density,sizes)

                if action:
                    file_act=directory+config_template+str(i)+"en.dat"
                    act_density,sizes=tools.read_top(file_act)
                    density_2d_act,sizes_big,index_smal=tools.projection_2d(act_density,sizes)
                else:
                    density_2d_act=np.array(density_2d_top)

                #print(rho)
                inst, a_inst, frac, a_frac, t_frac, t_inst, total=    tools.find_inst_2d(density_2d_top,density_2d_act,sizes_big,
                                   rho,norm,s/100,ss/100,norm_action,sss/100,neigh)

                Q_top=density_2d_top.sum()
                Q_instantons=len(inst)-len(a_inst)+1/2*len(frac)-1/2*len(a_frac)
                n+=len(frac)+len(a_frac)
                if (len(frac)+len(a_frac)) % 2:
                    delta_n+=1
                delta_q+=abs(round(Q_top)-Q_instantons)
            print(s/100,ss/100,sss/100)

            f.write(str(s/100)+"\t "+ str(ss/100) + "\t" + str(sss/100) + "\t" +
                    str(delta_n) + "\t" + str(delta_q)+ "\t" + str(n) +"\n")
        

f.close()


