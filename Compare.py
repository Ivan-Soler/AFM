import matplotlib.pyplot as plt
import numpy as np
import sys

import Read
import importlib
import analyzer
import re
import pickle
plt.rcParams.update({'font.size': 16})

#plt.rcParams['text.usetex'] = True
                
    
def Construct_susy(folder,susy_read_s0,susy_read_s1,conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes): #It assumes eigenvalues are ordered
    
    #Read supersymmetric modes up to a threshold of the eigenvalue
    density_susy=np.zeros(sizes[3])
    read=False
    for mode in range(0,susy_read_s1):
        if mode<max_modes:
            read=True
            Mode = folder+"sector_1/SusyMode_bin_"+str(mode)+"-"+str(conf)
            density_s1,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
            density_susy+=density_s1
        
    for mode in range(0,susy_read_s0):
        if mode<max_modes:
            read=True
            Mode = folder+"sector_0/SusyMode_bin_"+str(mode)+"-"+str(conf)
            density_s0,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
            density_susy-=density_s0
        
    #If there is no eigenvalue below the threshold we get the lowest one    
    if not read:
        ev1=float(dictionary_s1[str(conf)][0])
        ev2=float(dictionary_s0[str(conf)][0])
        if ev1<ev2:
            read=True
            Mode = folder+"sector_1/SusyMode_bin_"+str(0)+"-"+str(conf)
            density_s1,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
            density_susy+=density_s1
        else:
            read=True
            Mode = folder+"sector_0/SusyMode_bin_"+str(0)+"-"+str(conf)
            density_s0,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
            density_susy-=density_s0
    return(density_susy)   

def Construct_susy_improve(modes_s0, modes_s1,top):
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

def GM_RPO_cut(folder,sizes,max_modes,colors,spin_length,configurations,lambdas,RPO_threshold,save):
    
    #Param and definitions
    conf_tot=len(configurations)
    steps=len(lambdas)
    GM_tot=np.zeros((steps))
    RPO_tot=np.zeros((steps))      
    param=np.zeros((steps))
    
    #Parsing the eigenvalues
    dictionary_s1=analyzer.Real_eigenvalue(folder+"./sector_1/Measure.seq")
    dictionary_s0=analyzer.Real_eigenvalue(folder+"./sector_0/Measure.seq")

    #Checking wich configurations have non-trivial topological content
    top_gauge,conf_read=analyzer.Count_index_gf(folder,configurations)
    
    print(conf_read)
    ov_max=0
    susy_max=0 #To know in which step of lambda one runs out of values
    k=0
    for threshold in lambdas:
        #Check how many modes for each configuration we need to read 
        susy_read_s0=analyzer.Count_index(folder+"sector_0/Measure.seq",":OverlapFilterModeR:",threshold,conf_read)
        susy_read_s1=analyzer.Count_index(folder+"sector_1/Measure.seq",":OverlapFilterModeR:",threshold,conf_read)
        GM={}
        RPO={}
        j=0
        
        for conf in conf_read:
            #Read GF
            Topology =folder+"../gf/profile4dt2c"+str(conf)+"to.dat"
            density_top,sizes=Read.topology_1d(Topology)
            
            #Construct susy mode
            density_susy=Construct_susy(folder,susy_read_s0[conf],susy_read_s1[conf],conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes)
            
            #Compute the distance between the susy and the topological density
            GM[conf]=Geom_mean_1d(density_susy,density_top)
            RPO[conf]=Relative_point(np.absolute(density_susy),np.absolute(density_top),RPO_threshold)
            
        if save:   
            with open(folder+"./GM_hist.txt", 'wb') as f:
                pickle.dump(GM, f)
                sys.exit()

        #Store the means        
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
        print(k) #to see the progress
        if not count_meas==len(conf_read):
            print("GM measures left out")
        print(count_meas)
        print(len(conf_read))
        k+=1
        
        #To check when one runs out of modes
        end_overlap, end_susy=analyzer.End_spectrum(folder,conf_read) #Computes the last value of the spectrum
        i=0
        for threshold in lambdas:
            if not ov_max and threshold > end_overlap: #Computes on which step
                ov_max=[threshold,i]
            if not susy_max and threshold > end_susy:
                susy_max=[threshold,i]
            i+=1
        
    #Plotting and saving the GM and RPO
    np.savetxt(folder+"./end_spectrum.txt", np.vstack((ov_max,susy_max)))
    np.savetxt(folder+"./GM.txt",np.vstack((GM_tot,lambdas)))
    np.savetxt(folder+"./RPO.txt", np.vstack((RPO_tot,lambdas)))


def Index_dic(folder,lambdas,configurations):

    steps=len(lambdas)
    ov_dif_count=np.zeros((steps))
    susy_dif_count=np.zeros((steps)) 
    ov_dif_mean=np.zeros((steps))
    susy_dif_mean=np.zeros((steps)) 
    ov_max=0
    susy_max=0

    #Calculating where do you run out of eigenvalues
    end_overlap, end_susy=analyzer.End_spectrum(folder)
    
    #index configurations
    for threshold in lambdas:
        ov_top_dif, susy_top_dif, conf_tot = analyzer.Topology_dic_impr(folder,threshold,configurations)
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

def GF_vs_GF(folder,configurations, t_start,t_end,t_step,RPO_threshold):
    
    steps=int((t_end-t_start)/t_step)
    tot_conf=len(configurations)
    GM_tot=np.zeros((steps+1))
    RPO_tot=np.zeros((steps+1))
    
    k=0
    top_gauge,conf_read=analyzer.Count_index_gf(folder,configurations)
    for i in range(0, steps+1):
        t=t_start+t_step*i
        if int(t)==t: t=int(t)
        GM={}
        RPO={}
        tot_conf=0
        for conf in configurations:
            Topology_smooth=folder+"profile4dt2c"+str(conf)+"to.dat"
            density_smooth,sizes=Read.topology_1d(Topology_smooth)

            Topology_t=folder+"profile4dt"+str(t)+"c"+str(conf)+"to.dat"
            density_t,sizes=Read.topology_1d(Topology_t)

            GM[conf]=Geom_mean_1d(density_t,density_smooth)
            RPO[conf]=Relative_point(np.absolute(density_t),np.absolute(density_smooth),RPO_threshold)           
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


def Xi(densityA,densityB):
    Xi=0
    meanA=np.mean(densityA)
    meanB=np.mean(densityB)
    sizes=[len(densityA[0]),len(densityA[1]),len(densityA[2]),len(densityA[3])]
    volume=sizes[0]*sizes[1]*sizes[2]*sizes[3]
    meanA=np.mean(densityA)
    meanB=np.mean(densityB)
    for x in range(sizes[0]):
        for y in range(sizes[1]):
            for z in range(sizes[2]):
                for t in range(sizes[3]):
                    Xi+=((densityA[x,y,z,t]-meanA)*(densityB[x,y,z,t]-meanB))/volume
    return Xi

def Geom_mean(densityA,densityB):
    GM=(Xi(densityA,densityB)*Xi(densityA,densityB))/(Xi(densityA,densityA)*Xi(densityB,densityB))
    
    return GM

def Xi_1d(densityA,densityB):
    Xi=0
    meanA=np.mean(densityA)
    meanB=np.mean(densityB)
    sizes=[len(densityA)]
    volume=sizes[0]
    meanA=np.mean(densityA)
    meanB=np.mean(densityB)
    for t in range(sizes[0]):
        Xi+=((densityA[t]-meanA)*(densityB[t]-meanB))/volume
    return Xi

def Geom_mean_1d(densityA,densityB):
    GM=(Xi_1d(densityA,densityB)*Xi_1d(densityA,densityB))/(Xi_1d(densityA,densityA)*Xi_1d(densityB,densityB))
    
    return GM

def IPR(density):
    IPR=0
    for i in range(0,len(density)):
        IPR+=density[i]*density[i]
    return(len(density)*IPR)


def Relative_point(densityA,densityB,thresholdA):
    volA=0
    volB=0
    vol_union=0
    vol_inter=0
    for element in densityA:
        if element > thresholdA:
            volA+=1
    if volA==0:
        return(0)
    stop=0
    stop_max=10000
    thresholdB=thresholdA
    while (volA!=volB) and (stop<stop_max): 
        volB=0
        stop+=1
        for element in densityB:
            if element > thresholdB:
                volB+=1
        thresholdB-=0.01
    thresholdB+=0.01    
    if (stop<stop_max):
        for t in range(0,len(densityA)):
            if (densityA[t]>thresholdA) and (densityB[t]>thresholdB):
                vol_inter+=1
            if (densityA[t]>thresholdA) or (densityB[t]>thresholdB):
                vol_union+=1  
        return (vol_inter/vol_union)
    else:
        trehsoldB=thresholdA
        stop=0
        while (volA!=volB) and (stop<1000):
            
            volB=0
            stop+=1
            for element in densityB:
                if element > thresholdB:
                    volB+=1
            thresholdB+=0.01 
        thresholdB-=0.01    
        if (stop<stop_max):
            for t in range(0,len(densityA),1):
                if (densityA[t]>thresholdA) and (densityB[t]>thresholdB):
                    vol_inter+=1
                if (densityA[t]>thresholdA) or (densityB[t]>thresholdB):
                    vol_union+=1   
            return (vol_inter/vol_union)
        
        else:
            print("RPO volume did not converged")
            return (0)
