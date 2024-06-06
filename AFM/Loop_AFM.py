import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import sys
import Plot_generator
plt.rcParams.update({'font.size': 16})
plt.rcParams['text.usetex'] = False

folder_in=str(sys.argv[1]) #./gf_afm_1p5t/ 
tau_compare=str(sys.argv[2]) # 1.5
folder_out=str(sys.argv[3]) #./compare_1p5t/gf_afm_1p5t/  

sizes=[4,4,4,32]
max_modes=8
colors=3
spin_length=4
max_modes=8

folder_gf="/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/gf/"
top_gauge,conf_read=analyzer.Count_index_gf(folder_gf,configurations)

tony_folder="../cases/"
tony_data=np.loadtxt(tony_folder+"tony_instantons.txt")
susy_read_s1={}
susy_read_s0={}
for element in tony_data:
    susy_read_s1[element[2]]=element[0]
    susy_read_s0[element[2]]=element[1]
Plot_generator.susy_plot("/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/"+str(sys.argv[1]),folder_out,sizes,colors,spin_length,max_modes,conf_read,
                         susy_read_s0,susy_read_s1,Load=False,Plot=True)

