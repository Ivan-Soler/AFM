import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import sys
plt.rcParams.update({'font.size': 16})
#plt.rcParams['text.usetex'] = True

#Param and definitions
folder_in=str(sys.argv[1])
tau_compare=str(sys.argv[2])
folder_out=str(sys.argv[3])

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
conf=np.arange(conf_start,conf_end,conf_step)
folder_gf="../4x4x4x32/b2p44_new/gf/"

top_gauge,conf_read=analyzer.Count_index_gf(folder_gf,configurations)
#Compare.GM_RPO_cut(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf,lambdas,RPO_threshold,tao_compare,False)

tao_compare=4
folder_in="../4x4x4x32/b2p44_new/compare_"+str(tao_compare)+"p0t/"
folder_out=folder_in
measures=["gf_afm_4p0t/", "gf_afm_3p0t/", "gf_afm_2p0t/","gf_afm_1p5t/", 
          "gf_afm_1p125t/", "gf_afm_0p75t/", "gf_afm_0p5t/", "gf_afm_0p25t/","gf_afm_0p0t/"]
time_measures=[4,3,2,1.5,1.125,0.75,0.5,0.25,0]
obersvable="GM"
Plot_generator.susy_plot("../4x4x4x32/b2p44_new/"+measures[8],folder_out+measures[8],sizes,colors,spin_length,max_modes,lambda_opt,configurations)
Plot_generator.MC_history(folder_in,folder_out,measures,lambdas,observable)
Plot_generator.Cut_dependence(folder_in,folder_out,measures,observable)
folder_gf="../4x4x4x32/b2p44_new/gf/"
Plot_generator.GF_vs_AFM(folder_in, folder_gf, folder_out, conf_read, t_start, t_end, t_step,
                         #RPO_trehsold,tau_compare,measures,time_measures,observable)