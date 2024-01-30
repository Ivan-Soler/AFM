import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys

Sum_over=(0,1,2)
conf = str(sys.argv[1])
colors = 2
spin_length=4
max_modes_1=int(sys.argv[2])
max_modes_2=int(sys.argv[3])

sizes=[8,8,8,8]
t=np.arange(sizes[3])
susy_mode_1=np.zeros(sizes[3])
for i in range(0,max_modes_1):
    file="./sector_1/SusyMode_bin_"+str(i)+"-"+str(conf)
    density,sizes=Read.bin_mode_1d(file,sizes,colors,spin_length)
    susy_mode_1+=density/2

susy_mode_0=np.zeros(sizes[3])
for i in range(0,max_modes_2):
    file="./sector_0/SusyMode_bin_"+str(i)+"-"+str(conf)
    density,sizes=Read.bin_mode_1d(file,sizes,colors,spin_length)
    susy_mode_0+=-1*density/2

plt.title("Configuration:"+str(conf))
plt.plot(t,susy_mode_1+susy_mode_0, label="susy_mode")

Topology="../gf/profile4dt2c"+str(conf)+"to.dat"
#density_top,sizes=Read.topology_1d(Topology)
#plt.plot(t,density_top, label="top charge")
plt.legend(loc="upper right")
plt.title("Configuration:"+str(conf))
plt.savefig("./Susymode_"+str(conf)+".png", dpi=150, bbox_inches='tight')
plt.close()
