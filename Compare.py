import matplotlib.pyplot as plt
import numpy as np
import Plotting
import Read
import importlib
import analyzer
import re
import Maxima_find
import pickle
plt.rcParams.update({'font.size': 16})
#plt.rcParams['text.usetex'] = True

# Reads the eigenvalue of the supersymmetric operator of each mode
def Real_eigenvalue(file):
    Measure_file = open(file, "r")
    pattern=":OverlapFilterModeR:"
    eigenvectors=[[]]
    count=[[]]
    eigenvectors.pop(0)
    count={}
    dictionary={}
    for line in Measure_file:
        if re.search(pattern,line):
            line_split=line.split(":")
            if (float(line_split[8])<1):
                eigenvectors.append([line_split[1],line_split[4],line_split[8]])
                if line_split[1] in dictionary:
                    dictionary[line_split[1]].append(line_split[8])
                    count[line_split[1]]+=1
                else:
                    dictionary.update({line_split[1]:[line_split[8]]})
                    count.update({line_split[1]:1})
    return dictionary

def GM_matrix(modes1, modes2):
    GM=np.zeros((len(modes1),len(modes2)))
    for i in range(0,len(modes1)):
        for j in range(0,len(modes2)):
            GM[i,j]=analyzer.Geom_mean_1d(modes1(i),modes2(j))
    return(GM)
                
def Construct_susy(modes_s0, modes_s1,top):
    GM_00=Compare.GM_matrix(modes_s0,modes_s0)
    GM_01=Compare.GM_matrix(modes_s0,modes_s1)
    GM_11=Compare.GM_matrix(modes_s1,modes_s1)
    
    modes_s0_used=np.full((len(modes_s0),False))
    modes_s1_used=np.full((len(modes_s1),False))
    for i in range(0,len(modes_s0)):        
        for j in range(0,i):
            if GM_00[i,j]>0.9:
                repeated=True
        if not repeated:
            susy-=modes_s0
            modes_s0_used[i]=True
                
    for i in range(0,len(modes_s1)):        
        for j in range(0,i):
            if GM_11[i,j]>0.9:
                repeated=True
        if not repeated:
            susy+=modes_s1
            modes_s1_used[i]=True
                
    for i in range(0,len(modes_s1)):        
        for j in range(0,i):
            if GM_01[i,j]>0.9:
                repeated=True
        if repeated:
            if top[Maxima_find.simple_1d(modes_s1,sizes)]<0:
                susy-=modes_s1
                modes_s1_used[i]=False
            else:
                susy+=modes_s0
                modes_s0_used[i]=False
    return(susy, modes_s0_used, modes_s1_used)     
                
def GM_RPO(folder,sizes,max_modes,colors,spin_length,conf_start,conf_end,conf_step,
           threshold_start,threshold_end,steps,RPO_threshold):
    
    #Param and definitions
    conf_tot=int((conf_end-conf_start)/conf_step)
    
    GM_tot=np.zeros((steps))
    RPO_tot=np.zeros((steps))      
    param=np.zeros((steps))
    
    threshold_steps=(threshold_end-threshold_start)/steps
    dictionary_s1=Real_eigenvalue(folder+"./sector_1/Measure.seq")
    dictionary_s0=Real_eigenvalue(folder+"./sector_0/Measure.seq")
    end_overlap, end_susy=analyzer.End_spectrum(folder)
    
    for k in range(0, steps):
        
        threshold=threshold_start+threshold_steps*k
        param[k]=threshold
        GM={}
        RPO={}
        j=0
        for conf in range(conf_start,conf_end,conf_step):
            #Read GF
            Topology =folder+"../gf/profile4dt4c"+str(conf)+"to.dat"
            density_top,sizes=Read.topology_1d(Topology)
            #max_top=density_top[Maxima_find(density_top,sizes)]
            #min_top=density_top[Maxima_find(density_top,sizes)]
            #Read supersymmetric modes up to a threshold of the eigenvalue
            density_susy=np.zeros(sizes[3])
            read=False #needed to know if there is some topological object in the measure
            for i in range(0,max_modes):
                ev=float(dictionary_s1[str(conf)][i])
                if (abs(ev)<threshold):
                    read=True
                    Mode = folder+"sector_1/SusyMode_bin_"+str(i)+"-"+str(conf)
                    density_s1,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
                    density_susy+=density_s1

            for i in range(0,max_modes):
                ev=float(dictionary_s0[str(conf)][i])
                if (abs(ev)<threshold):
                    read=True
                    Mode = folder+"sector_0/SusyMode_bin_"+str(i)+"-"+str(conf)
                    density_s0,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
                    density_susy+=-density_s0
            if read:
                GM[conf]=analyzer.Geom_mean_1d(density_susy,density_top)
                RPO[conf]=analyzer.RPO(np.absolute(density_susy),np.absolute(density_top),RPO_threshold)
            print(conf)
        GM_mean=0
        RPO_mean=0
        count_meas=0
        for key in GM:
            GM_mean+=GM[key]
            RPO_mean+=RPO[key]
            count_meas+=1
        RPO_mean=RPO_mean/count_meas  
        GM_mean=GM_mean/count_meas

        GM_tot[k]=GM_mean
        RPO_tot[k]=RPO_mean
        print(k)

    #Plotting and saving the GM and RPO
    np.savetxt(folder+"./GM.txt",np.vstack((GM_tot,param)))
    np.savetxt(folder+"./RPO.txt", np.vstack((RPO_tot,param)))

def GM_RPO_cut(folder,sizes,max_modes,colors,spin_length,conf_start,conf_end,conf_step,
           threshold_start,threshold_end,steps,RPO_threshold):
    
    #Param and definitions
    conf_tot=int((conf_end-conf_start)/conf_step)
    
    GM_tot=np.zeros((steps))
    RPO_tot=np.zeros((steps))      
    param=np.zeros((steps))
    
    threshold_steps=(threshold_end-threshold_start)/steps
    dictionary_s1=Real_eigenvalue(folder+"./sector_1/Measure.seq")
    dictionary_s0=Real_eigenvalue(folder+"./sector_0/Measure.seq")
    
    end_overlap, end_susy=analyzer.End_spectrum(folder)

    ov_top_dif,susy_top_dif,conf_tot,susy_tot=analyzer.Topology_dic(folder+"../gf_afm_4p0t/",0.175) #To know the good configurations, susy tot is a dictionary with the information of the configurations being used
    
    for k in range(0, steps):
        
        threshold=threshold_start+threshold_steps*k
        param[k]=threshold
        GM={}
        RPO={}
        j=0
        for conf in range(conf_start,conf_end,conf_step):
            #Read GF
            Topology =folder+"../gf/profile4dt4c"+str(conf)+"to.dat"
            density_top,sizes=Read.topology_1d(Topology)
            #max_top=density_top[Maxima_find(density_top,sizes)]
            #min_top=density_top[Maxima_find(density_top,sizes)]
            #Read supersymmetric modes up to a threshold of the eigenvalue
            density_susy=np.zeros(sizes[3])
            read=False
            if susy_tot[str(conf)]:
                
                for i in range(0,max_modes):
                    ev=float(dictionary_s1[str(conf)][i])
                    if (abs(ev)<threshold):
                        read=True
                        Mode = folder+"sector_1/SusyMode_bin_"+str(i)+"-"+str(conf)
                        density_s1,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
                        density_susy+=density_s1

                for i in range(0,max_modes):
                    ev=float(dictionary_s0[str(conf)][i])
                    if (abs(ev)<threshold):
                        read=True
                        Mode = folder+"sector_0/SusyMode_bin_"+str(i)+"-"+str(conf)
                        density_s0,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
                        density_susy+=-density_s0
                if not read:
                    ev1=float(dictionary_s1[str(conf)][0])
                    ev2=float(dictionary_s0[str(conf)][0])
                    if ev1<ev2:
                        read=True
                        Mode = folder+"sector_1/SusyMode_bin_"+str(i)+"-"+str(conf)
                        density_s1,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
                        density_susy+=density_s1
                    else:
                        read=True
                        Mode = folder+"sector_0/SusyMode_bin_"+str(i)+"-"+str(conf)
                        density_s0,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
                        density_susy+=-density_s0
                if read:
                    GM[conf]=analyzer.Geom_mean_1d(density_susy,density_top)
                    RPO[conf]=analyzer.RPO(np.absolute(density_susy),np.absolute(density_top),RPO_threshold)
                print(conf)
                
        GM_mean=0
        RPO_mean=0
        count_meas=0
        for key in GM:
            GM_mean+=GM[key]
            RPO_mean+=RPO[key]
            count_meas+=1
        RPO_mean=RPO_mean/count_meas  
        GM_mean=GM_mean/count_meas

        GM_tot[k]=GM_mean
        RPO_tot[k]=RPO_mean
        print(k)

    #Plotting and saving the GM and RPO
    np.savetxt(folder+"./GM.txt",np.vstack((GM_tot,param)))
    np.savetxt(folder+"./RPO.txt", np.vstack((RPO_tot,param)))

    
def GM_RPO_impr(folder,sizes,max_modes,colors,spin_length,conf_start,conf_end,conf_step,
           threshold_start,threshold_end,steps,RPO_threshold):
    
    #Param and definitions
    conf_tot=int((conf_end-conf_start)/conf_step)
    
    GM_tot=np.zeros((steps))
    RPO_tot=np.zeros((steps))      
    param=np.zeros((steps))
    
    threshold_steps=(threshold_end-threshold_start)/steps
    dictionary_s1=Real_eigenvalue(folder+"./sector_1/Measure.seq")
    dictionary_s0=Real_eigenvalue(folder+"./sector_0/Measure.seq")
    end_overlap, end_susy=analyzer.End_spectrum(folder)
    
    table_used_modes={}
    
    for k in range(0, steps):
        
        threshold=threshold_start+threshold_steps*k
        param[k]=threshold
        GM={}
        RPO={}
        j=0
        for conf in range(conf_start,conf_end,conf_step):
            #Read GF
            density_s0=np.zeros((max_modes))
            density_s1=np.zeros((max_modes))
            Topology =folder+"../gf/profile4dt4c"+str(conf)+"to.dat"
            density_top,sizes=Read.topology_1d(Topology)
            #max_top=density_top[Maxima_find(density_top,sizes)]
            #min_top=density_top[Maxima_find(density_top,sizes)]
            #Read supersymmetric modes up to a threshold of the eigenvalue
            density_susy=np.zeros(sizes[3])
            read=False #needed to know if there is some topological object in the measure
            for i in range(0,max_modes):
                ev=float(dictionary_s1[str(conf)][i])
                if (abs(ev)<threshold):
                    read=True
                    Mode = folder+"sector_1/SusyMode_bin_"+str(i)+"-"+str(conf)
                    density_s1[i],sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)

            for i in range(0,max_modes):
                ev=float(dictionary_s0[str(conf)][i])
                if (abs(ev)<threshold):
                    read=True
                    Mode = folder+"sector_0/SusyMode_bin_"+str(i)+"-"+str(conf)
                    density_s0[i],sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
            density_susy, used_modes_s0, used_modes_s1=Compare.Construct_susy(density_s0, density_s1,density_top)
            
            ov_top_dif, susy_top_dif, conf_tot = analyzer.Topology_dic_impr(folder,threshold,used_modes_s0,used_modes_s1)
            for key in ov_top_dif:
                ov_dif_mean[i]+=abs(ov_top_dif[key])
                if abs(ov_top_dif[key]>0.1): ov_dif_count[i]+=1
            for key in susy_top_dif:
                susy_dif_mean[i]+=abs(susy_top_dif[key])
                if abs(susy_top_dif[key]>0.1): susy_dif_count[i]+=1
            
            #Adding the step and the lambda where you run out of eigenmodes
            if not ov_max and threshold > end_overlap:
                ov_max=[threshold,i]
            if not susy_max and threshold > end_susy:
                susy_max=[threshold,i]
            
            if read:
                GM[conf]=analyzer.Geom_mean_1d(density_susy,density_top)
                RPO[conf]=analyzer.RPO(np.absolute(density_susy),np.absolute(density_top),RPO_threshold)
            print(conf)
        GM_mean=0
        RPO_mean=0
        for key in GM:
            GM_mean+=GM[key]
            RPO_mean+=RPO[key]
        RPO_mean=RPO_mean/conf_tot  
        GM_mean=GM_mean/conf_tot

        GM_tot[k]=GM_mean
        RPO_tot[k]=RPO_mean
        print(k)

    #Plotting and saving the GM and RPO
    np.savetxt(folder+"./GM.txt",np.vstack((GM_tot,param)))
    np.savetxt(folder+"./RPO.txt", np.vstack((RPO_tot,param)))

def Index_dic(folder,lambda_min,lambda_max,steps):

    ov_dif_count=np.zeros((steps))
    susy_dif_count=np.zeros((steps)) 
    ov_dif_mean=np.zeros((steps))
    susy_dif_mean=np.zeros((steps)) 
    param=np.zeros((steps))
    ov_max=0
    susy_max=0
    
    threshold_step=(lambda_max-lambda_min)/steps
    #Calculating where do you run out of eigenvalues
    end_overlap, end_susy=analyzer.End_spectrum(folder)
    
    #index configurations
    for i in range(0, steps):
        threshold=lambda_min+threshold_step*i
        param[i]=threshold
        ov_top_dif, susy_top_dif, conf_tot = analyzer.Topology_dic_impr(folder,threshold)
        for key in ov_top_dif:
            ov_dif_mean[i]+=abs(ov_top_dif[key])
            if abs(ov_top_dif[key]>0.1): ov_dif_count[i]+=1
        for key in susy_top_dif:
            susy_dif_mean[i]+=abs(susy_top_dif[key])
            if abs(susy_top_dif[key]>0.1): susy_dif_count[i]+=1
            
        #Adding the step and the lambda where you run out of eigenmodes
        if not ov_max and threshold > end_overlap:
            ov_max=[threshold,i]
        if not susy_max and threshold > end_susy:
            susy_max=[threshold,i]
    
    #Saving the results
    np.savetxt(folder+"./end_spectrum.txt", np.vstack((ov_max,susy_max)))
    np.savetxt(folder+"./Index_ov.txt", np.vstack((ov_dif_count/conf_tot,param)))
    np.savetxt(folder+"./Top_ov.txt", np.vstack((ov_dif_mean,param)))
    np.savetxt(folder+"./Index_susy.txt",  np.vstack((susy_dif_count/conf_tot,param)))
    np.savetxt(folder+"./Top_susy.txt", np.vstack((susy_dif_mean,param)))
    np.savetxt(folder+"./Top_susy.txt", np.vstack((susy_dif_mean,param)))
    
    with open('used_conf.pkl', 'wb') as f:
        pickle.dump(susy, f)
    return()

def GF(folder,conf_start, conf_end, conf_step,t_start,t_end,t_step,RPO_threshold):
    
    steps=int((t_end-t_start)/t_step)
    tot_conf=int((conf_end-conf_start)/conf_step)
    GM_tot=np.zeros((steps+1))
    RPO_tot=np.zeros((steps+1))
    
    k=0
    ov_top_dif,susy_top_dif,conf_tot,susy_tot=analyzer.Topology_dic(folder+"../gf_afm_4p0t/",0.175)
    for i in range(0, steps+1):
        t=t_start+t_step*i
        if int(t)==t: t=int(t)
        GM={}
        RPO={}
        tot_conf=0
        for conf in range(conf_start, conf_end, conf_step):
            if susy_tot[str(conf)]: 
                Topology_smooth=folder+"profile4dt4c"+str(conf)+"to.dat"
                density_smooth,sizes=Read.topology_1d(Topology_smooth)

                Topology_t=folder+"profile4dt"+str(t)+"c"+str(conf)+"to.dat"
                density_t,sizes=Read.topology_1d(Topology_t)

                GM[conf]=analyzer.Geom_mean_1d(density_t,density_smooth)
                RPO[conf]=analyzer.RPO(np.absolute(density_t),np.absolute(density_smooth),RPO_threshold)           
                tot_conf+=1
        GM_mean=0
        RPO_mean=0
        for key in GM:
            GM_mean+=GM[key]
            RPO_mean+=RPO[key]
        GM_tot[k]=GM_mean/tot_conf
        RPO_tot[k]=RPO_mean/tot_conf
        k+=1
        print(tot_conf)
    param=np.arange(t_start, t_end+t_step, t_step)
    np.savetxt(folder+"./GM.txt",np.vstack((GM_tot,param)))
    np.savetxt(folder+"./RPO.txt", np.vstack((RPO_tot,param)))  
    
    return(GM)
            
            
def GM_RPO_modes(folder,sizes,max_modes,colors,spin_length,conf_start,conf_end,conf_step,RPO_threshold):
    
    #Param and definitions
    conf_tot=int((conf_end-conf_start)/conf_step)
    
    GM_tot=np.zeros((max_modes+1))
    RPO_tot=np.zeros((max_modes+1))      
    param=np.zeros((max_modes+1))
    
    #dictionary_s1=Real_eigenvalue(folder+"./sector_1/Measure.seq")
    #dictionary_s0=Real_eigenvalue(folder+"./sector_0/Measure.seq")
    end_overlap, end_susy=analyzer.End_spectrum(folder)
    
    for k in range(1, max_modes+1):
        GM={}
        RPO={}
        j=0
        for conf in range(conf_start,conf_end,conf_step):
            #Read GF
            Topology =folder+"../gf/profile4dt4c"+str(conf)+"to.dat"
            density_top,sizes=Read.topology_1d(Topology)
            #max_top=density_top[Maxima_find(density_top,sizes)]
            #min_top=density_top[Maxima_find(density_top,sizes)]
            #Read supersymmetric modes up to a threshold of the eigenvalue
            density_susy=np.zeros(sizes[3])
            read=False #needed to know if there is some topological object in the measure
            for i in range(0,k):
                Mode = folder+"sector_1/SusyMode_bin_"+str(i)+"-"+str(conf)
                density_s1,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
                density_susy+=density_s1

            for i in range(0,k):
                #ev=float(dictionary_s0[str(conf)][i])
                #if (abs(ev)<threshold):
                Mode = folder+"sector_0/SusyMode_bin_"+str(i)+"-"+str(conf)
                density_s0,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
                density_susy+=-density_s0
            GM[conf]=analyzer.Geom_mean_1d(density_susy,density_top)
            RPO[conf]=analyzer.RPO(np.absolute(density_susy),np.absolute(density_top),RPO_threshold)
            print(conf)
        GM_mean=0
        RPO_mean=0
        count_meas=0
        for key in GM:
            GM_mean+=GM[key]
            RPO_mean+=RPO[key]
            count_meas+=1
        RPO_mean=RPO_mean/count_meas  
        GM_mean=GM_mean/count_meas

        GM_tot[k]=GM_mean
        RPO_tot[k]=RPO_mean
        print(k)

    #Plotting and saving the GM and RPO
    np.savetxt(folder+"./GM.txt",np.vstack((GM_tot,param)))
    np.savetxt(folder+"./RPO.txt", np.vstack((RPO_tot,param)))        
        