import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d, Axes3D
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interpn
import Projecting
import Plotting
import Read
import Maxima_find
import Fitting
import Fit_1d
import subprocess
from subprocess import call
import importlib
plt.rcParams.update({'font.size': 12})

folder="./Marga/b2p55"
meas="/gf/"

directory=folder+meas
name="profile4dt"

conf="1001"
pref="c"+conf+"to.dat"

time_max=20
time_step=1

for time in range(0,time_max,time_step):
    file_name=directory+name+str(time)+pref
    density_top,sizes=Read.density_1d(file_name)
    file=directory+"top_density_"+str(time)
    figure=Plotting.plot_list(sizes[3],density_top, file)
    plt.show()