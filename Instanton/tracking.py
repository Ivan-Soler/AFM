import matplotlib.pyplot as plt
import numpy as np
import tools
import sys
import tarfile


#----------Program to fit instantons on a T²xR² lattice under flow------------------------#
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
conf_number=param[12]

new_t_frac=[]
target=[]
cap=10
start=50
end=51
setp=3

for t in range(start,end,step):
    old_t_frac=new_t_frac.copy()
    file_top=directory+config_template+str(t)+"c"+str(conf_number)+"to.dat"
    top_density,sizes=tools.read_top(file_top)
    density_2d_top,sizes_big,index_smal=tools.projection_2d(top_density,sizes)
    
    if action:
        file_act=directory+config_template+str(t)+"c"+str(conf_number)+"en.dat"
        act_density,sizes=tools.read_top(file_act)
        density_2d_act,sizes_big,index_smal=tools.projection_2d(act_density,sizes)
        norm_action=norm*8*np.pi*np.pi
    else:
        density_2d_act=np.array(density_2d_top)
        eps_action=eps_norm
        norm_action=norm
    #We perform the fit
    inst, a_inst, frac, a_frac, t_frac, t_inst, total=    tools.find_inst_2d(density_2d_top,density_2d_act,sizes_big,
                       rho,norm,eps_rho,eps_norm,norm_action,eps_action,neigh)

    #print(t_frac)
    new_t_frac=tools.compare_fit(old_t_frac, t_frac, cap)
    
    #we track the element we want to check
    #t_frac=tools.find_max_2d(-density_2d_act,sizes_big)
    #print(t_inst)
    tools.plot_dens_2d(file_top,-density_2d_top,sizes_big, t_frac,t_inst)
    tools.plot_dens_2d(file_act,-density_2d_act,sizes_big, t_frac,t_inst)
    for element in new_t_frac:
        file="_fit_t"+str(t)
        if element[0]==28 and element[1]==85:
            target.append(element)
            if element[2][2]!=0:
                ax=tools.plot_inst(sizes_big,element[2],directory,file,"blue","+")
                ax.plot()
        if element[0]==42 and element[1]==91:
            target.append(element)
            if element[2][2]!=0:
                ax2=tools.plot_inst(sizes_big,element[2],directory,file,"green","-",ax)
    plt.ylabel("q(x)")
    plt.xlabel("x")
    plt.title("Fitted fractionals t="+str(t))
    plt.savefig("./fit.png",dpi=150)
    plt.close()
    
    #plt.plot(-density_2d_top[82,:])
    #plt.show()

 
    
f = open("./"+"tracking.txt", "w")
for element in target:
    f.write(str(element))
    f.write("\n")
f.close()

f = open(directory+"fractionals.txt", "w")
for element in new_t_frac:
    f.write(str(element))
    f.write("\n")
f.close()


