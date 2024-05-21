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

density,sizes=Read.topology_1d(sys.argv[1])
t=np.linspace(0,sizes[3],sizes[3])
plt.plot(t,density)
plt.savefig("./"+sys.argv[1] + ".png",dpi=150, bbox_inches='tight')
