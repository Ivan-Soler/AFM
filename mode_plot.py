import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys

Sum_over=(0,1,2)
Mode = str(sys.argv[1])
Bin = int(sys.argv[2])
Chirality= 1
colors = 3
spin_length=4

if Bin:
    sizes=[8,8,8,64]
    density_mode,sizes=Read.bin_mode(Mode,sizes,colors,spin_length)
else:
    density_mode,sizes=Read.ascii_mode(Mode)  

density_mode=Chirality*density_mode.sum(axis=Sum_over)

t=np.linspace(0,sizes[3],sizes[3])
plt.plot(t,density_mode)
Mode.replace("SusyMode_bin_","")
Mode.split("-")
file_name=("SusyMode_"+Mode[1]+"_"+Mode[0]
plt.savefig(Mode+".png",dpi=150, bbox_inches='tight')
