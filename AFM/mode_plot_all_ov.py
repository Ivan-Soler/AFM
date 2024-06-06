import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys

Sum_over=(0,1,2)
max_conf = int(sys.argv[1])
All_summed = int(sys.argv[2])
Chirality= int(sys.argv[3])
colors = 3
spin_length=4
max_modes=12

if All_summed:
    sizes=[4,4,4,32]
    t=np.arange(sizes[3])
    for conf in range(10,max_conf,10):
        susy_mode_1=np.zeros(sizes[3])
        file="./OverlapMode"+str(conf)+"_"
        i=0
        susy_mode=np.zeros((spin_length, colors, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
        while i<max_modes:
            mode,density,sizes=Read.ascii_mode(file+str(i))
            susy_mode+=mode/np.sqrt(2)
            i+=1
            
            mode,density,sizes=Read.ascii_mode(file+str(i))
            susy_mode+=mode/np.sqrt(2)
            i+=1
            
        density_susy=Read.mode_to_density(susy_mode,colors,spin_length,sizes)
        plt.title("Configuration:"+str(conf))
        plt.plot(t,Chirality*density_susy, label="susy_mode")
        Topology="../../gf/profile4dt2c"+str(conf)+"to.dat"
        density_top,sizes=Read.topology_1d(Topology)
        plt.plot(t,density_top, label="top charge")
        plt.legend(loc="upper right")
        plt.title("Configuration:"+str(conf))
        plt.savefig("./Susymode_"+str(conf)+".png", dpi=150, bbox_inches='tight')
        plt.close()
else:
    sizes=[4,4,4,32]
    t=np.arange(sizes[3])
    for conf in range(10,max_conf,10):
        susy_mode_1=np.zeros(sizes[3])
        file="./OverlapMode"+str(conf)+"_"
        i=0
        while i<max_modes:
            susy_mode=np.zeros((spin_length, colors, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
            mode,density,sizes=Read.ascii_mode(file+str(i))
            susy_mode+=mode/np.sqrt(2)
            i+=1
            
            mode,density,sizes=Read.ascii_mode(file+str(i))
            susy_mode+=mode/np.sqrt(2)
            
            
            density_susy=Read.mode_to_density(susy_mode,colors,spin_length,sizes)
            plt.title("Configuration:"+str(conf))
            plt.plot(t,Chirality*density_susy, label="susy_mode")
            Topology="../../gf/profile4dt2c"+str(conf)+"to.dat"
            density_top,sizes=Read.topology_1d(Topology)
            plt.plot(t,density_top, label="top charge")
            plt.legend(loc="upper right")
            plt.title("Configuration:"+str(conf))
            plt.savefig("./Overlapmode_"+str(conf)+"_"+str(i)+".png", dpi=150, bbox_inches='tight')
            plt.close()
            
            i+=1