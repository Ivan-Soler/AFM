import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import sys
plt.rcParams.update({'font.size': 16})
#plt.rcParams['text.usetex'] = True

#Param and definitions
folder=str(sys.argv[1])
sizes=[4,4,4,32]
max_modes=8
colors=3
spin_length=4

conf_start=10
conf_end=1000
conf_step=10
conf=np.arange(conf_start,conf_end,conf_step)

lambda_min=0.01
lambda_max=0.20
steps=20
lambdas=np.linspace(lambda_min,lambda_max,num=steps)

RPO_threshold=0.15

#Compare.Index_dic(folder,lambdas, conf)
#Compare.GM_RPO_cut(folder,sizes,max_modes,colors,spin_length,conf,lambdas,RPO_threshold)

ov_max, susy_max=np.loadtxt(folder+"end_spectrum.txt")
lambda_opt=susy_max[0]
analyzer.susy_plot(folder,sizes,colors,spin_length,max_modes,lambda_opt,conf)
