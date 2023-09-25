import numpy as np
import Plotting
import Read
import analyzer
import sys

Sum_over=(1,2)
Mode = str(sys.argv[1])
Bin = int(sys.argv[2])
Chirality= int(sys.argv[3])
Topology = str(sys.argv[4])
colors=int(str(sys.argv[5]))
spin_length=4

density_top,sizess=Read.topology(Topology)

if Bin:
    density_mode,sizes=Read.bin_mode(Mode,sizes,colors,spin_length)
else:
    density_mode,sizes=Read.ascii_mode(Mode)
GM=analyzer.Geom_mean(Chirality*density_mode,density_top)    
print(Mode,GM)

density_mode=Chirality*density_mode.sum(axis=Sum_over)
density_top=density_top.sum(axis=Sum_over)

Plotting.plot_density_2d(Mode,density_mode,sizes)
Plotting.plot_density_2d(Topology,density_top,sizes)


