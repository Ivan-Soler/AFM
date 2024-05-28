
import matplotlib.pyplot as plt
import numpy as np
import tools
import sys
import tarfile
import os
import pandas as pd
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
norm_frac=param[6]
eps_rho=param[7]
eps_norm=param[8]
neigh=param[9]
action=param[10]
eps_action=param[11]
norm_inst=param[13]
beta=param[14]
lt=param[15]
lr=param[16]
tau=param[17]

#print(action)
f = open("identification"+str(lt)+"nt"+str(lr)+"nr"+str(tau)+"t"+str(beta)+"b.txt", "w")
#f.write(str(beta)+ " " + str(lt) + " "+ str(lr) + "\n")
delta_n=0
mean_frac=0
diff_frac=0
count=0

arclist=pd.read_pickle(directory+'dataprofile.pkl')
flowtime=tau
d=arclist[arclist['FlowTime'] == float(flowtime)]
d.sort_values('ConfigNumber',inplace=True)

lt_height={
   "4":0.07,
   "5":0.05,
   "6":0.035,
   "7":0.025,
   "8":0.02,
   "9":0.015,
   "10":0.01
}

lt_frac={
   "4":0.01,
   "5":0.01,
   "6":0.01,
   "7":0.01,
   "8":0.01,
   "9":0.005,
   "10":0.005
}


#norm_inst=lt_height[str(lt)]
#norm_frac=lt_frac[str(lt)]

norm_frac=0

#print(norm_inst)
#print(norm_frac)

error=open("../error.txt","a")
progress=open("../progress.txt","a")
for index,row in d.iterrows():
   fname=row['ArchiveName']
   filen=row['FileName']
   #print(fname+"/"+filen)
   if "to.dat" in filen and "dt"+str(tau) in filen:
       #print(filen)
       try:
           tar= tarfile.open(directory+fname)
       except PermissionError:
           error.write("permission error for " + directory+fname+"\n")
           continue
       count+=1
       conf=filen.replace("profile4dt"+str(tau)+"c", "")
       conf=conf.replace("to.dat", "")
       tar.extract(filen)
       
       top_density,sizes=tools.read_top(filen)
       #print(filen)
       density_2d_top,sizes_big,index_smal=tools.projection_2d(top_density,sizes)

       inst, a_inst, frac, a_frac, t_frac, t_inst, total= tools.find_inst_2d(density_2d_top,sizes_big,
                                                                              norm_frac,norm_inst,neigh)

       Q_top=density_2d_top.sum()
       mean_frac+= len(t_frac)
       diff_frac+= len(frac)-len(a_frac)
       Q_instantons= len(inst)-len(a_inst)+1/2*len(frac)-1/2*len(a_frac)

       if (len(frac)+len(a_frac)) % 2:
            delta_n+=1
            #maxima=tools.find_max_2d(density_2d_top,sizes_big)
            #tools.plot_dens_2d(filen,density_2d_top,sizes_big, t_frac, t_inst)
            #print(filen)
            
       f.write(str(conf)+" " + str(len(frac)) + " " +str(len(a_frac))+ " "+ 
                str(len(inst)) + " " +str(len(a_inst))+ " " +
                str(Q_top))
       for element in t_frac:
           f.write( " "+ str(element[0])+ " "+str(element[1]))
           for fit in element[2]:
              f.write( " " +str(fit))
           #if element[2][4]>0:
              #print("positive")
           #if element[2][4]<0:
           #   print(element[2])
       f.write(" \n")
       os.remove(filen)
       if not count%20:
          tools.plot_dens_2d(filen,density_2d_top,sizes_big, [], [])

       #print(len(t_frac),len(t_inst))

print("Configurations with fractional topological charge: " + str(delta_n) + "\n")
print("Mean frac_inst + frac_anti_inst = " + str(mean_frac/count)+ "\n")
print("Mean frac_inst - frac_anti_inst = " + str(diff_frac/count)+ "\n")
print("Number of configurations = " + str(count) + "\n")
progress.write(directory+fname+"\n")
progress.close()
error.close()
