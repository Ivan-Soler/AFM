import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import Read
import Maxima
import re
import os
import sys
import pickle
plt.rcParams.update({"text.usetex": False, "font.size": 8})


def susy_plot(folder_in,folder_out,sizes,colors,spin_length,max_modes,conf_read,susy_read_s0, susy_read_s1, Load=False,Plot=True,Polyakov=False):
    
    dictionary_s1=analyzer.Real_eigenvalue(folder_in+"./sector_1/Measure.seq")
    dictionary_s0=analyzer.Real_eigenvalue(folder_in+"./sector_0/Measure.seq")
    
    if not Polyakov:
        print(conf_read)
        for conf in conf_read:
            if ( susy_read_s1[str(conf)]==0) and ( susy_read_s0[str(conf)]==0):
                print(conf)
                ev1=float(dictionary_s1[str(conf)][0])
                ev2=float(dictionary_s0[str(conf)][0])
                if ev1<ev2:
                    print("s1")
                    susy_read_s1[str(conf)]=1
                    susy_read_s0[str(conf)]=0
                else:
                    print("s2")
                    susy_read_s0[str(conf)]=1
                    susy_read_s1[str(conf)]=0
        with open(folder_out+"modes_used.txt", 'w') as f:
            for conf in conf_read:
                print(conf, susy_read_s1[str(conf)],susy_read_s0[str(conf)], file=f)
        print(susy_read_s1)
        print(susy_read_s0)
    for conf in conf_read:
        #Read GF
        Topology_1=folder_in+"../gf/profile4dt0.5c"+str(conf)+"to.dat"
        Topology_2=folder_in+"../gf/profile4dt2c"+str(conf)+"to.dat"
        Topology_3=folder_in+"../gf/profile4dt4c"+str(conf)+"to.dat"
        density_top_1,sizes=Read.topology_1d(Topology_1)
        density_top_2,sizes=Read.topology_1d(Topology_2)
        density_top_3,sizes=Read.topology_1d(Topology_3)
        normalization=np.sum(np.abs(density_top_3))

        #Construct susy mode
        if Load:
            density_susy=np.loadtxt(folder+measure+"susy_mode_"+str(conf)+"c.txt")
        else:
            density_susy=Compare.Construct_susy(folder_in,susy_read_s0[str(conf)],susy_read_s1[str(conf)],
                                                conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes)
            np.savetxt(folder_in+"susy_mode_"+str(conf)+"c.txt", density_susy)

        #density_susy=density_susy*(normalization/np.sum(np.abs(density_susy)))

        #Plot the three densities
        plt.plot(density_top_1, label='top. density t=0.5')
        plt.plot(density_top_2, label='top. density t=2')
        plt.plot(density_top_3, label='top. density t=4')
        plt.plot(density_susy, label="AFM")
        plt.legend(loc="lower left", ncol=2)
        plt.savefig(folder_out+"./susy_mode_"+str(conf)+"c.png",dpi=150, bbox_inches='tight')
        plt.close()      
        np.savetxt(folder_out+"./susy_mode_"+str(conf)+"c.txt",density_susy)
        
    return()

def MC_history(folder_in,folder_out,measures,lambdas,observable_name,Plot=False,Polyakov=False):
    
    for measure in (measures):
        index_lambda=0
        print(measure)
        for element in lambdas:
            observable={}
            data=np.loadtxt(folder_in+measure+"./"+observable_name+"_hist_"+str(index_lambda)+".txt",delimiter=" ", dtype=float)
            for element in data:
                observable[str(element[0])] = element[1]
            #Compute mean, error and variance
            x=[]
            y=[]
            observable_mean=0
            variance=0
            for key in observable:
                observable_mean+=float(observable[key])
            observable_mean/=len(observable)
            for key in observable:
                variance+=(observable[key]-observable_mean)**2
                x.append(float(key)/10)
                y.append(observable[key])
            error=np.sqrt(variance)/len(observable)
            
            #Save the errors
            with open(folder_out+measure+"./"+observable_name+"_error_"+str(index_lambda)+".txt", 'w') as f:
                f.write(str(error))
                print(str(index_lambda))
            if Plot:
                plt.scatter(x,y, marker="x")
                plt.hlines(observable_mean, xmin=0, xmax=100, linestyle="--")
               #plt.xlabel(r'Configuration')
               # plt.ylabel(r'$$ \mbox{\huge $ \Xi$}$$')
                plt.xlabel('Configuration')
                plt.ylabel('Xi')

                plt.xticks(np.arange(0, 120,  step=20))
                plt.savefig(folder_out+measure+""+observable_name+"_history_"+str(index_lambda)+".pdf",dpi=150, bbox_inches='tight')
                plt.close()
            #Save the optimal one
            if not Polyakov:
                f=open(folder_in+measure+"lambda_opt.txt",'r')
                lamba_string=f.read().split('\n')
                lambda_opt,index_opt=float(lamba_string[0]), int(float(lamba_string[1]))
                f.close()
                if index_opt==index_lambda:
                    f_hist=open(folder_out+measure+observable_name+"_history_opt.txt", 'w')
                    for element in observable:
                        print(str(element)+"\t"+str(observable[element]),file=f_hist)
                        plt.savefig(folder_out+measure+""+observable_name+"_history_opt.pdf",dpi=150, bbox_inches='tight')
                    f_hist.close()
                    plt.close()
            index_lambda+=1
            
    return()

def Cut_dependence(folder_in,folder_out,measures,observable):
    ax = plt.gca()
    for measure in (measures):
        time=re.sub("gf_afm", "", measure)
        time=re.sub("t","", time)
        time=re.sub("p",".", time)
        time=re.sub("/","", time)
        time=re.sub("_","", time)
        data=np.loadtxt(folder_in+measure+observable+".txt")
        susy_max=np.loadtxt(folder_in+measure+"end_spectrum.txt")

        #Read the error
        error=[]
        for t in range(0,len(data[0])):
            with open(folder_in+measure+observable+"_error_"+str(t)+".txt", 'r') as f:
                error.append(float(f.readline()))
            f.close()

        color = next(ax._get_lines.prop_cycler)['color']
        plt.errorbar(data[1],data[0], yerr=error, label=r'$\tau$'+"="+time, color=color)
        plt.fill_between(data[1], data[0]-error, data[0]+error, alpha=0.1, color=color)
        plt.scatter(susy_max[0],data[0,int(susy_max[1])], marker="v", color=color)

    #plt.xlabel(r'$$ \mbox{\huge $\lambda$}_{cut} $$')
    #plt.ylabel(r'$$ \mbox{\huge $ \Xi$}$$')
    plt.xlabel('lambda')
    plt.ylabel('Xi')
    plt.legend(loc="upper right", ncol=1)
    plt.ylim([0,1.1])
    plt.xlim([0.0,0.15])

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 1, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc="upper right", ncol=1, bbox_to_anchor=(1.40, 1))

    plt.savefig(folder_out+observable+"_cut.pdf",dpi=150, bbox_inches='tight')
    plt.close()
    
    return()


def Xi_max(data,susy_max):
    maximum=0
    for i in range(0,len(data[0])):
            if (data[0,i] > maximum) and (data[1,i]<susy_max[0]):
                maximum=data[0,i]
    return(maximum)

def find_max(folder,measure,observable):
    data=np.loadtxt(folder+measure+observable+".txt")
    susy_max=np.loadtxt(folder+measure+"end_spectrum.txt")
    maximum=Xi_max(data,susy_max)
    return(maximum)

def GF_vs_AFM(folder_in, folder_gf, folder_out, configurations, t_start, t_end, t_step, 
              RPO_trehsold,tau_compare,measures,time_measures,observable,Polyakov=False):

    GM_GF=Compare.GF_vs_GF(folder_gf, folder_out, configurations, t_start, t_end, t_step, RPO_trehsold, tau_compare)

    maximum=np.zeros((len(measures)))
    error=np.zeros((len(measures)))
    flow_time=0
    for measure in (measures):
        data=np.loadtxt(folder_in+measure+observable+".txt")
        #Reading lambda optimal
        if Polyakov:
            lambda_opt,index_opt=[0,0]
            maximum[flow_time]=data[0]
        else:
            f=open(folder_in+measure+"lambda_opt.txt",'r')
            lamba_string=f.read().split('\n')
            lambda_opt,index_opt=float(lamba_string[0]), int(float(lamba_string[1]))
            f.close()
            maximum[flow_time]=data[0,index_opt]
        #Readinf error
        f=open(folder_in+measure+observable+"_error_"+str(index_opt)+".txt",'r')
        error[flow_time]=float(f.read())
        f.close()
        flow_time+=1

    
    plt.xlabel("t")
    plt.ylabel(r'$\Xi$')
    plt.plot(GM_GF[0],GM_GF[1], label="GF")
    plt.errorbar(time_measures,maximum, yerr=error, color="orange", label="GF + AFM")
    plt.fill_between(time_measures, maximum-error, maximum+error,color="orange",alpha=0.2)
    plt.legend(loc="lower right")
    plt.savefig(folder_out+"./GF_AFM.pdf",dpi=150, bbox_inches='tight')
    plt.close()
    return()

def histogram(folder,file_name):
    
    GM_hist=[]
    with open(folder+file_name) as f:
        for line in f:
            GM_hist.append(float(line.split()[3].replace("[","").replace("]","").replace(",","")))
    print(GM_hist)
    plt.hist(GM_hist, bins='auto')
    plt.savefig(folder+"GM_hist.pdf",dpi=150, bbox_inches='tight')
    
    return()
    