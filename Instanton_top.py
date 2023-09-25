import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d, Axes3D
import sys

import Plotting
import Read
import Maxima_find
import Fitting
plt.rcParams.update({'font.size': 12})


file_name="../"+str(sys.argv[1])

density,sizes=Read.density(file_name)

