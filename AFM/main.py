import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import sys
import Plot_generator
plt.rcParams.update({'font.size': 8})
plt.rcParams['text.usetex'] = False

#Definitions parameters
colors=3
spin_length=4
RPO_threshold=0.15
folder_gf="./gf/"

#Param read from screen
folder_in=str(sys.argv[1]) #./gf_afm_1p5t/ 
tao_compare=str(sys.argv[2]) # 1.5
folder_out=str(sys.argv[3])+str(sys.argv[1]) #./compare_1p5t/gf_afm_1p5t/   
operator=str(sys.argv[4]) #OverlapFilterModeC
modes=sys.argv[5] # true, sum over modes, false, sum over densities

#Param read from file
f=open(str(sys.argv[3])+"main_parameters.txt", 'r')
sizes_str=f.readline().replace("\n","").split(" ")
sizes=[int(element) for element in sizes_str]
max_modes=int(f.readline())
method=f.readline().replace("\n","")
lambda_min=float(f.readline())
lambda_max=float(f.readline())
steps=int(f.readline())
lambdas=np.linspace(lambda_min,lambda_max,num=steps)

conf_start=int(f.readline())
conf_end=int(f.readline())
conf_step=int(f.readline())
conf=np.arange(conf_start,conf_end,conf_step)
f.close()

top_gauge,conf_read=analyzer.Count_index_gf(folder_gf,conf)
Compare.GM_RPO_cut(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf_read,lambdas,RPO_threshold,tao_compare,operator,method,modes)

f=open(folder_out+"lambda_opt.txt",'r')
lamba_string=f.read().split('\n')
lambda_opt,index_opt=float(lamba_string[0]), int(float(lamba_string[1]))

if method=="cut":
    susy_read_s0,susy_read_s1,start_spectrum = analyzer.Count_index_cut(folder_in,"",lambda_opt,conf_read,max_modes,operator)
if method=="gap":
    susy_read_s0,susy_read_s1,start_spectrum = analyzer.Count_index_gap(folder_in,"",lambda_opt,conf_read,max_modes,operator)
print(str(sys.argv[1]),"End main")
#print(str(sys.argv[1]),"Doubler contribution")
#Compare.GM_doublers(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf_read,tao_compare,susy_read_s0,susy_read_s1,operator,save=True)
#print(str(sys.argv[1]),"Doubler contribution at cut")
#print(susy_read_s0)
#print(susy_read_s1)
#Compare.GM_doublers_cut(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf_read,tao_compare,susy_read_s0,susy_read_s1,operator,save=True)
