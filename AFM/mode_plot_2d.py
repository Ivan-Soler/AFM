import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys
import Maxima
Sum_over=(1,2)
Mode = str(sys.argv[1])
Bin = int(sys.argv[2])
Chirality= int(sys.argv[3])
colors = 8
spin_length=4

if Bin:
    #sizes=[4,4,4,32]
    sizes=[40,6,6,40]
    density_mode,sizes=Read.bin_mode(Mode,sizes,colors,spin_length)
    #Mode=Mode.replace("./SusyMode_bin_","")
    #Mode=Mode.split("-")
    #file_name="SusyMode_"+Mode[1]+"_"+Mode[0]
    file_name=Mode
else:
    zero_mode,density_mode,sizes=Read.ascii_mode(Mode)  
    file_name=Mode

density_2d=Chirality*density_mode.sum(axis=(Sum_over))
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# Make data.
X = np.arange(0,sizes[0])
Y = np.arange(0,sizes[3])
X, Y = np.meshgrid(X, Y)

#ax.set_zlim3d(-0.05,0)
ax.set_box_aspect((1,1,1)) 

# Plot the surface.
surf = ax.plot_surface(X, Y, density_2d, cmap='viridis')

plt.savefig(file_name+".png",dpi=150, bbox_inches='tight')
maxima=Maxima.simple(Chirality*density_mode,sizes)
print(maxima)
 
