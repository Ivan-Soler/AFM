import matplotlib.pyplot as plt
import numpy as np
import Plotting
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
    density_mode,sizes=Read.bin_mode(Mode,sizes,colors,spin_length)
else:
    density_mode,sizes=Read.ascii_mode(Mode)  

density_mode=Chirality*density_mode.sum(axis=Sum_over)

Plotting.plot_density_1d(Mode,density_mode)

plt.savefig(Mode+".png",dpi=150, bbox_inches='tight')
