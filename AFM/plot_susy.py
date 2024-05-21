import matplotlib.pyplot as plt
import numpy as np
import Read
import Maxima
import re
import os
import sys
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import importlib
importlib.reload(Read)
importlib.reload(Maxima)

sizes=[4,4,4,32]
colors=3
spin=4

if sys.argv[2]:
    sizes=[16,8,8,8]
    density,sizes=Read.bin_mode(sys.argv[1],sizes,colors,spin)
    density_2d=density.sum(axis=(1,2))
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    X=np.arange(0,sizes[3])
    Y=np.arange(0,sizes[0])
    
    X, Y = np.meshgrid(X, Y)
    
    ax.plot_surface(X,Y,density_2d, cmap="viridis")
    plt.savefig("./"+sys.argv[1] + ".png",dpi=150, bbox_inches='tight')
    
else:
    density,sizes=Read.bin_mode_1d(sys.argv[1],sizes,colors,spin)
    t=np.linspace(0,sizes[3],32)
    plt.plot(t,density)
    plt.savefig("./"+sys.argv[1] + ".png",dpi=150, bbox_inches='tight')