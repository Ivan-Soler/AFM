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

lambda_min=0.01
lambda_max=0.15
steps=30
lambdas=np.linspace(lambda_min,lambda_max,num=steps)
RPO_threshold=0.15

conf_start=10
conf_end=1000
conf_step=10
configurations=np.arange(conf_start,conf_end,conf_step)
folder_gf="/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/gf/"

top_gauge,conf_read=analyzer.Count_index_gf(folder_gf,configurations)
#Compare.GM_RPO_cut(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf_read,lambdas,RPO_threshold,tau_compare,save=True)

tau_compare=4
folder_in="/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/compare_"+str(tau_compare)+"p0t/"
folder_out=folder_in
measures=["gf_afm_4p0t/", "gf_afm_3p0t/", "gf_afm_2p0t/","gf_afm_1p5t/", 
          "gf_afm_1p125t/", "gf_afm_0p75t/", "gf_afm_0p5t/", "gf_afm_0p25t/","gf_afm_0p0t/"]
time_measures=[4,3,2,1.5,1.125,0.75,0.5,0.25,0]


observable="GM"
#Plot_generator.MC_history(folder_out,folder_out,measures,lambdas,observable,Plot=True)
#Plot_generator.Cut_dependence(folder_in,folder_out,measures,observable)

t_start=0
t_end=4
t_step=0.25
RPO_threshold=0.15
#Plot_generator.GF_vs_AFM(folder_in, folder_gf, folder_out, conf_read, t_start, t_end, t_step,
#                         RPO_threshold,tau_compare,measures,time_measures,observable)

folder_in=str(sys.argv[1]) #./gf_afm_1p5t/                                                                                                                                                                 
tau_compare=str(sys.argv[2]) # 1.5                                                                                                                                                                         
folder_out=str(sys.argv[3]) #./compare_1p5t/gf_afm_1p5t/ 


#f=open(folder_out+"lambda_opt.txt",'r')
#lamba_string=f.read().split('\n')
#lambda_opt,index_opt=float(lamba_string[0]), int(float(lamba_string[1]))
#f.close()
#susy_read_s0=analyzer.Count_index(folder_in+"sector_0/Measure.seq",":OverlapFilterModeR:",lambda_opt,conf_read)
#susy_read_s1=analyzer.Count_index(folder_in+"sector_1/Measure.seq",":OverlapFilterModeR:",lambda_opt,conf_read)
#Plot_generator.susy_plot("/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/"+str(sys.argv[1]),folder_out,sizes,colors,spin_length,max_modes,conf_read,
#                         susy_read_s0,susy_read_s1,Load=False,Plot=True)

tony_folder="../cases/"
tony_data=np.loadtxt(tony_folder+"tony_instantons.txt")
susy_read_s1={}
susy_read_s0={}
for element in tony_data:
    susy_read_s1[element[2]]=element[0]
    susy_read_s0[element[2]]=element[1]
Plot_generator.susy_plot("/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/"+str(sys.argv[1]),folder_out,sizes,colors,spin_length,max_modes,conf_read,
                         susy_read_s0,susy_read_s1,Load=False,Plot=True)