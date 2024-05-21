import matplotlib.pyplot as plt
import numpy as np
import subprocess
import Read
import re
plt.rcParams.update({'font.size': 12})

#Example of reading and plotting topological summing over spatial dimensions
file="./6x32_su2/gf/profile4dt4c680to.dat"
density_top,sizes=Read.topology_1d(file)
t=np.arange(0,sizes[3])
plt.plot(t, density_top)
plt.savefig("./6x32_su2/gf/profile4dt4c680to.png")
plt.close()
          
    
    
#Example of reading and plotting Susy Mode summing over spatial dimensions
file_susy="./6x32_su2/afm/sector_0/SusyMode0-680"
color_range=3
spin_range=4
density_susy_1d,sizes=Read.bin_mode_1d(file_susy,sizes,color_range,spin_range)
plt.plot(t, density_susy_1d)
plt.savefig("./6x32_su2/afm/sector_0/SusyMode0-680.png")
plt.close()



#Example of reading Susy Mode in 4D. You need to pass by the sizes, color and spin range
file_susy="./6x32_su2/afm/sector_0/SusyMode0-680"
color_range=3
spin_range=4
sizes=[6,6,6,32]
density_susy_4d=Read.bin_mode(file_susy,sizes,color_range,spin_range)



#Example of reading Overlap Modes
file_overlap="./6x32_su2/afm/sector_0/OverlapMode6800"
#Reading and plotting summing over spatial dimensions
density_overlap_1d,sizes=Read.ascii_mode_1d(file_overlap)
plt.plot(t, density_susy_1d)
plt.savefig("./6x32_su2/afm/sector_0/OverlapMode6800.png")
plt.close()
#Reading in 4 D
density_overlap_1d,sizes=Read.ascii_mode(file_overlap)