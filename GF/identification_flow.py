import pickle 
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
norm_frac=0
tarfiles=os.listdir("profiles")
for directory in tarfiles:
      try:
          tar= tarfile.open("/data/datatopoym/iscalero/T2xR2/flow_13ns/runns13nt64b2.6nr64/profiles/"+directory)
      except Exception as e:
          error.write("permission error for " + directory+fname+"\n")
          continue
      list_files=tar.getmembers()
      for profile in list_files:
        filen=profile.name
        print(filen)
        if "to" in filen and "4dt"+str(tau)+"c" in filen:
          conf=filen.replace("profile4dt"+str(tau)+"c", "")
          conf=filen.replace("to.dat", "")
      
          en_file=filen.replace("to","en")
          tar.extract(filen)
          tar.extract(en_file)
          
      
          top_density,sizes=tools.read_top(filen)
          en_density,sizes=tools.read_top(en_file)
      
          density_2d_top,sizes_big,index_smal=tools.projection_2d(top_density,sizes)
          density_2d_en,sizes_big,index_smal=tools.projection_2d(en_density,sizes)
      
          inst, a_inst, frac, a_frac, t_frac, t_inst, total= tools.find_inst_2d(density_2d_top,density_2d_en,sizes_big,
                                                                    norm_frac,norm_inst,neigh)
      
          Q_top=density_2d_top.sum()
          S_en=density_2d_en.sum()
          f.write(str(conf)+" " + str(len(frac)) + " " +str(len(a_frac))+ " "+ 
          str(len(inst)) + " " +str(len(a_inst))+ " " +
          str(Q_top) + " "+ str(S_en))
          for element in t_frac:
              f.write( " "+ str(element[0])+ " "+str(element[1]))
              for fit in element[2]:
                  f.write( " " +str(fit))
              for fit in element[3]:
                  f.write( " " +str(fit))
          f.write(" \n")
          os.remove(filen)
          os.remove(en_file)