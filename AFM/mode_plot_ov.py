import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys

Sum_over=(0,1,2)
Mode_name=str(sys.argv[1])
Mode_number=str(sys.argv[2])
Chirality= 1
colors = 3
spin_length=4
sizes=[4,4,4,32]

susy_mode=np.zeros((spin_length, colors, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
for i in range(0,int(Mode_number)):
    mode,density,sizes=Read.ascii_mode(Mode_name+str(i))
    susy_mode+=mode/np.sqrt(2)

density_susy=Read.mode_to_density(susy_mode,colors,spin_length,sizes)

t=np.arange(0,sizes[3])
plt.plot(t,density_susy/2)

plt.savefig(Mode_name+"density.png",dpi=150, bbox_inches='tight')
