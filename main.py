import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import sys
import Plot_generator
plt.rcParams.update({'font.size': 8})
plt.rcParams['text.usetex'] = False

#Param and definitions
folder_in=str(sys.argv[1]) #./gf_afm_1p5t/ 
tao_compare=str(sys.argv[2]) # 1.5
folder_out=str(sys.argv[3]) #./compare_1p5t/gf_afm_1p5t/         

sizes=[4,4,4,32]
max_modes=8
colors=3
spin_length=4

lambda_min=0.01
lambda_max=0.15
steps=20
lambdas=np.linspace(lambda_min,lambda_max,num=steps)
RPO_threshold=0.15

conf_start=10
conf_end=1000
conf_step=10
conf=np.arange(conf_start,conf_end,conf_step)
folder_gf="/home/mi37fud/b2p44_new/gf/"

top_gauge,conf_read=analyzer.Count_index_gf(folder_gf,conf)

Compare.GM_RPO_cut(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf_read,lambdas,RPO_threshold,tao_compare)

f=open(folder_out+"lambda_opt.txt",'r')
lamba_string=f.read().split('\n')
lambda_opt,index_opt=float(lamba_string[0]), int(float(lamba_string[1]))

susy_read_s0, susy_read_s1=analyzer.Count_index_all(folder_in,"",lambda_opt,conf_read)
#Compare.GM_doublers(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf_read,tao_compare,susy_read_s0,susy_read_s1,save=True)
#Compare.GM_doublers_cut(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf_read,tao_compare,susy_read_s0,susy_read_s1,save=True)
