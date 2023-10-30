import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import sys
import Plot_generator
plt.rcParams.update({'font.size': 16})
plt.rcParams['text.usetex'] = False

#Param and definitions
folder_modes=str(sys.argv[1]) #./b2p44_new/ 
tau_compare=str(sys.argv[2]) # 1.5
folder_out=str(sys.argv[3]) #./compare_1p5t/         

sizes=[4,4,4,32]
max_modes=8
colors=3
spin_length=4

lambda_min=0.001
lambda_max=0.1
steps=3
lambdas=np.linspace(lambda_min,lambda_max,num=steps)
RPO_threshold=0.15

conf_start=10
conf_end=1000
conf_step=10
configurations=np.arange(conf_start,conf_end,conf_step)
folder_gf="/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/gf/"
top_gauge,conf_read=analyzer.Count_index_gf(folder_gf,configurations)
print(len(conf_read))

#Plots about GM
measures=np.loadtxt(folder_out+"measures.txt", dtype=str)[0]
time_measures=np.array(np.loadtxt(folder_out+"measures.txt", dtype=str)[1],dtype=float)
print(measures)
observable="GM"
Plot_generator.MC_history(folder_out,folder_out,measures,lambdas,observable,Plot=True)
Plot_generator.Cut_dependence(folder_out,folder_out,measures,observable)

#Plot with the GF
t_start=0
t_end=1.5
t_step=0.25
Plot_generator.GF_vs_AFM(folder_out, folder_gf, folder_out, conf_read, t_start, t_end, t_step,
                         RPO_threshold,tau_compare,measures,time_measures,observable)

#for measure in measures:
#    susy_read_s0, susy_read_s1 = analyzer.Count_index_all("./"+measure,"",1000,conf_read,max_modes)
#    Plot_generator.histogram(folder_out+measure,"GM_doublers", conf_read, max_modes,susy_read_s0, susy_read_s1)
#    f=open(folder_out+measure+"lambda_opt.txt",'r') #mesure[0] to take the threshold at 0p5t
#    lamba_string=f.read().split('\n')
#    lambda_opt,index_opt=float(lamba_string[0]), int(float(lamba_string[1]))
#    f.close()
#    print(measure)
#    susy_read_s0, susy_read_s1 = analyzer.Count_index_all("./"+measure,"",lambda_opt,conf_read,max_modes)
#    Plot_generator.histogram(folder_out+measure,"GM_doublers_cut", conf_read, max_modes,susy_read_s0, susy_read_s1)
#    print(folder_out+measure)
    #Plot_generator.susy_plot(folder_modes+measure,folder_out+measure,sizes,colors,
    #                         spin_length,max_modes,conf_read,susy_read_s0,susy_read_s1)


measure=measures[0]
f=open(folder_out+measure+"lambda_opt.txt",'r') #mesure[0] to take the threshold at 0p5t
lamba_string=f.read().split('\n')
lambda_opt,index_opt=float(lamba_string[0]), int(float(lamba_string[1]))
susy_read_s0, susy_read_s1 = analyzer.Count_index_all("./"+measure,"",lambda_opt,conf_read,max_modes)
Plot_generator.susy_plot(folder_modes+measure,folder_out+measure,sizes,colors,
                             spin_length,max_modes,conf_read,susy_read_s0,susy_read_s1)
