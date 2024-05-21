import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys

Sum_over=(0,1,2)
max_conf = int(sys.argv[1])
Bin = int(sys.argv[2])
Chirality= int(sys.argv[3])
colors = 3
spin_length=4
max_modes=12

if Bin:
    sizes=[8,8,8,64]
    t=np.arange(sizes[3])
    for conf in range(120,max_conf,20):
        susy_mode_1=np.zeros(sizes[3])
        for i in range(0,max_modes):
            file="./SusyMode_bin_"+str(i)+"-"+str(conf)
            density,sizes=Read.bin_mode_1d(file,sizes,colors,spin_length)
            susy_mode_1+=Chirality*density/2
        plt.title("Configuration:"+str(conf))
        plt.plot(t,susy_mode_1, label="susy_mode")
        Topology="../../gf/profile4dt2c"+str(conf)+"to.dat"
        density_top,sizes=Read.topology_1d(Topology)
        plt.plot(t,density_top, label="top charge")
        plt.legend(loc="upper right")
        plt.title("Configuration:"+str(conf))
        plt.savefig("./Susymode_"+str(conf)+".png", dpi=150, bbox_inches='tight')
        plt.close()
else:
    density_mode,sizes=Read.ascii_mode(Mode)  
