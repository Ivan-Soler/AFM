import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import sys
import Plot_generator
plt.rcParams.update({'font.size': 16})
plt.rcParams['text.usetex'] = False

#Definitions parameters
colors=3
spin_length=4
RPO_threshold=0.15
folder_gf="./gf/"


#Param read from screen
folder_modes=str(sys.argv[1]) #./b2p44_new/ 
tau_compare=str(sys.argv[2]) # 1.5
folder_out=str(sys.argv[3]) #./compare_1p5t/   
operator=str(sys.argv[4]) 

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
print(len(conf_read))

#Plots about GM
f=open(folder_out+"measures.txt", 'r')
measures=f.readline().replace(" \n","").split(" ")
time_measures=f.readline().replace(" \n","").split(" ")
print(measures)
print(time_measures)
f.close()

observable="GM"
print("Monte carlo history of Xi")
print(folder_out)
#Plot_generator.MC_history(folder_out,folder_out,measures,lambdas,observable,Plot=True)
print("Xi dependence with cut")
#Plot_generator.Cut_dependence(folder_out,folder_out,measures,observable)
print(folder_out)
#Plot with the GF
t_start=0
t_end=4
t_step=0.25
print("GF vs AFM")
print(folder_out)
#Plot_generator.GF_vs_AFM(folder_out, folder_gf, folder_out, conf_read, t_start, t_end, t_step,
#                         RPO_threshold,tau_compare,measures,time_measures,observable)

for measure in measures:
    print(measure)
    #susy_read_s0, susy_read_s1 = analyzer.Count_index_cut("./"+measure,"",1000,conf_read,max_modes,operator) #Just to have the maximum number of modes
    print("histogram doublers all modes")
    #Plot_generator.histogram(folder_out+measure,"GM_doublers", conf_read, max_modes,susy_read_s0, susy_read_s1)
    f=open(folder_out+measure+"lambda_opt.txt",'r') #mesure[0] to take the threshold at 0p5t
    lamba_string=f.read().split('\n')
    lambda_opt,index_opt=float(lamba_string[0]), int(float(lamba_string[1]))
    f.close()
    if method=="cut":
        susy_read_s0,susy_read_s1,start_spectrum = analyzer.Count_index_cut("./"+measure,"",lambda_opt,conf_read,max_modes,operator)
    elif method=="gap":
        susy_read_s0,susy_read_s1 = analyzer.Count_index_gap("./"+measure,"",lambda_opt,conf_read,max_modes,operator)
    else:
        print("method "+ method +"not in list")
    #print("histogram doublers cut")
    #Plot_generator.histogram(folder_out+measure,"GM_doublers_cut", conf_read, max_modes,susy_read_s0, susy_read_s1)
    print("susy modes at cut")
    Plot_generator.susy_plot(folder_modes+measure,folder_out+measure,sizes,colors,
                             spin_length,max_modes,conf_read,susy_read_s0,susy_read_s1,operator,"optimal")
    #cut_i=0
    #for cut in lambdas:
    #    print("susy modes at ", str(cut))
    #    lambda_opt,index_opt=float(cut), cut_i
    #    if method=="cut":
    #        susy_read_s0,susy_read_s1,start_spectrum = analyzer.Count_index_cut("./"+measure,"",lambda_opt,conf_read,max_modes,operator)
    #    if method=="gap":
    #        susy_read_s0,susy_read_s1 = analyzer.Count_index_gap("./"+measure,"",lambda_opt,conf_read,max_modes,operator)
    #        Plot_generator.susy_plot(folder_modes+measure,folder_out+measure,sizes,colors,
    #                                 spin_length,max_modes,conf_read,susy_read_s0,susy_read_s1,operator,str(cut_i))
    #    cut_i+=1
    
