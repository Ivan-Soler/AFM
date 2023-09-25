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
conf_end=500
conf_step=10
conf_tot=int((conf_end-conf_start)/conf_step)

lambda_min=0.01
lambda_max=0.20
steps=1

RPO_threshold=0.15

Compare.Index_dic(folder,lambda_min,lambda_max,steps)
#Compare.GM_RPO_cut(folder,sizes,max_modes,colors,spin_length,conf_start,conf_end,conf_step,lambda_min,lambda_max,steps,RPO_threshold)

#Compare.GM_RPO_modes(folder,sizes,max_modes,colors,spin_length,conf_start,conf_end,conf_step,RPO_threshold)