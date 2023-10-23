import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import sys
import Plot_generator
plt.rcParams.update({'font.size': 16})
plt.rcParams['text.usetex'] = False

#Param and definitions
tau_compare=str(sys.argv[1])
folder_in=str(sys.argv[2])
folder_out=folder_in
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

observable="GM"
measures=["gf_afm_4p0t/", "gf_afm_3p0t/", "gf_afm_2p0t/","gf_afm_1p5t/", 
          "gf_afm_1p125t/", "gf_afm_0p75t/", "gf_afm_0p5t/", "gf_afm_0p25t/","gf_afm_0p0t/"]
time_measures=[4,3,2,1.5,1.125,0.75,0.5,0.25,0]

measures=[ "gf_afm_2p0t/","gf_afm_1p5t/",
          "gf_afm_1p125t/", "gf_afm_0p75t/", "gf_afm_0p5t/", "gf_afm_0p25t/","gf_afm_0p0t/"]
time_measures=[2,1.5,1.125,0.75,0.5,0.25,0]

measures=["gf_afm_1p5t/","gf_afm_1p125t/", "gf_afm_0p75t/", "gf_afm_0p5t/", "gf_afm_0p25t/","gf_afm_0p0t/"] 
time_measures=[1.5,1.125,0.75,0.5,0.25,0]  

lambdas=[0]
Plot_generator.MC_history(folder_out,folder_out,measures,lambdas,observable,Plot=True,Polyakov=True)

t_start=0
t_end=float(tau_compare)
t_step=0.25
RPO_threshold=0.15
Plot_generator.GF_vs_AFM(folder_out, folder_gf, folder_out, conf_read, t_start, t_end, t_step,
                         RPO_threshold,tau_compare,measures,time_measures,observable,Polyakov=True)

tony_folder="/home/mi37fud/b2p44_new/cases/"
tony_data=np.loadtxt(tony_folder+"tony_instantons.txt",dtype=int)
susy_read_s1={}
susy_read_s0={}

for element in tony_data:
    susy_read_s1[str(element[2])]=element[0]
    susy_read_s0[str(element[2])]=element[1]

for measure in measures:
    Plot_generator.susy_plot("./"+measure,folder_out+measure,sizes,colors,spin_length,max_modes,conf_read,
                         susy_read_s0,susy_read_s1,Load=False,Plot=True,Polyakov=True)

