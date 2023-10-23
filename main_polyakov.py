import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import sys
import Plot_generator
plt.rcParams.update({'font.size': 16})
plt.rcParams['text.usetex'] = False

#Param and definitions
folder_in=str(sys.argv[1]) #./gf_afm_1p5t/ 
tau_compare=str(sys.argv[2]) # 1.5
folder_out=str(sys.argv[3]) #./compare_1p5t/gf_afm_1p5t/         

sizes=[4,4,4,32]
max_modes=8
colors=3
spin_length=4

conf_start=10
conf_end=1000
conf_step=10
configurations=np.arange(conf_start,conf_end,conf_step)
folder_gf="/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/gf/"

top_gauge,conf_read=analyzer.Count_index_gf(folder_gf,configurations)
tony_folder="/home/mi37fud/b2p44_new/cases/"
tony_data=np.loadtxt(tony_folder+"tony_instantons.txt",dtype=int)
susy_read_s1={}
susy_read_s0={}

for element in tony_data:
    susy_read_s1[str(element[2])]=element[0]
    susy_read_s0[str(element[2])]=element[1]

Compare.GM_RPO("/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/"+str(sys.argv[1]),folder_out,sizes,max_modes,colors,spin_length,conf_read,tau_compare,susy_read_s0,susy_read_s1)
