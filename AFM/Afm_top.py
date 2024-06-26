analyzer.py                                                                                         0000644 0001750 0001750 00000030530 14607720211 012131  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import Read
import math
import re
import Compare

plt.rcParams.update({'font.size': 12})        

# Reads the eigenvalue of the supersymmetric operator of each mode
def Real_eigenvalue(file,pattern):
    Measure_file = open(file, "r")
    dictionary={}
    for line in Measure_file:
        if re.search(pattern,line):
            line_split=line.split(":")
            if line_split[1] in dictionary:
                dictionary[line_split[1]].append(float(line_split[8]))
            else:
                dictionary.update({line_split[1]:[float(line_split[8])]})
    return dictionary

def GM_matrix(modes1, modes2):
    GM=np.zeros((len(modes1),len(modes2)))
    for i in range(0,len(modes1)):
        for j in range(0,len(modes2)):
            GM[i,j]=analyzer.Geom_mean_1d(modes1(i),modes2(j))
    return(GM)


def End_spectrum(folder,configurations):
    #Reads the last eigenvalue of each configuration and returns the smallest one
    #Overlap
    end_overlap=100
    highest_overlap_s0={}
    with open(folder+"sector_0/Measure.seq") as file:
        for line in file:
            match=re.search(":OverlapFilterModeC:",line)
            if match:
                string=line.split(":")
                if not string[1] in highest_overlap_s0:
                    highest_overlap_s0[string[1]]=abs(float(string[8]))
                else:
                    if (abs(float(string[8]))>highest_overlap_s0[string[1]]):
                        highest_overlap_s0[string[1]]=abs(float(string[8]))
                    
    #Overlap
    highest_overlap_s1={}
    with open(folder+"sector_1/Measure.seq") as file:
        for line in file:
            match=re.search(":OverlapFilterModeC:",line)
            if match:
                string=line.split(":")
                if not string[1] in highest_overlap_s1:
                    highest_overlap_s1[string[1]]=abs(float(string[8]))
                else:
                    if (abs(float(string[8]))>highest_overlap_s1[string[1]]):
                        highest_overlap_s1[string[1]]=abs(float(string[8]))
      
    for key in highest_overlap_s0:
        if highest_overlap_s0[key]<end_overlap:
            end_overlap=highest_overlap_s0[key]

    for key in highest_overlap_s1:
        if highest_overlap_s1[key]<end_overlap:
            end_overlap=highest_overlap_s1[key]
            
    #susy
    end_susy=100
    highest_susy_s0={}
    with open(folder+"sector_0/Measure.seq") as file:
        for line in file:
            match=re.search(":OverlapFilterModeR:",line)
            if match:
                string=line.split(":")
                for conf in configurations:
                    if string[1]==str(conf):
                        if not string[1] in highest_susy_s0:
                            highest_susy_s0[string[1]]=abs(float(string[8]))
                        else:
                            if (abs(float(string[8]))>highest_susy_s0[string[1]]):
                                highest_susy_s0[string[1]]=abs(float(string[8]))
                                
    #susy
    highest_susy_s1={}
    with open(folder+"sector_1/Measure.seq") as file:
        for line in file:
            match=re.search(":OverlapFilterModeR:",line)
            if match:
                string=line.split(":")
                for conf in configurations:
                    if string[1]==str(conf):
                        if not string[1] in highest_susy_s1:
                            highest_susy_s1[string[1]]=abs(float(string[8]))
                        else:
                            if (abs(float(string[8]))>highest_susy_s1[string[1]]):
                                highest_susy_s1[string[1]]=abs(float(string[8]))
                    
    for key in highest_susy_s1:
        if highest_susy_s1[key]<end_susy:
            end_susy=highest_susy_s1[key]
            conf_m=key

    for key in highest_susy_s0:
        if highest_susy_s0[key]<end_susy:
            end_susy=highest_susy_s0[key]
            conf_m=key
  
    return(end_susy,end_overlap,conf_m)

def Count_index(folder,measurement,threshold,configurations,max_modes):
    count={}
    for conf in configurations:
        count[str(conf)]=0

    with open(folder) as file:
        for line in file:
            match=re.search(measurement,line)
            if match:
                string=line.split(":")
                if (abs(float(string[8]))<threshold) and (string[1] in count):
                    if count[string[1]]<max_modes:
                        count[string[1]]+=1

    return(count)

def Count_index_cut(folder_modes,measure,threshold,conf_read,max_modes,pattern):
    
    dictionary_s1=Real_eigenvalue(folder_modes+"/sector_1/Measure.seq",pattern)
    dictionary_s0=Real_eigenvalue(folder_modes+"/sector_0/Measure.seq",pattern)
                    
    susy_read_s0=Count_index(folder_modes+measure+"/sector_0/Measure.seq",
                                      pattern,threshold,conf_read,max_modes)
    susy_read_s1=Count_index(folder_modes+measure+"/sector_1/Measure.seq",
                                      pattern,threshold,conf_read,max_modes)
    start_spectrum=False
    conf_min=0
    for conf in conf_read:
        if susy_read_s0[str(conf)]==0 and susy_read_s1[str(conf)]==0:
            start_spectrum=False
            conf_min=conf
            break
            #ev1=float(dictionary_s1[str(conf)][0])
            #ev2=float(dictionary_s0[str(conf)][0])
            #if ev1<ev2:
            #    susy_read_s1[str(conf)]=1
            #else:
            #    susy_read_s0[str(conf)]=1

    
    return(susy_read_s0,susy_read_s1,start_spectrum,conf_min)

def Count_index_gap(folder_modes,measure,gap,conf_read,max_modes,pattern):
    
    dictionary_s1=Real_eigenvalue(folder_modes+"/sector_1/Measure.seq",pattern)
    dictionary_s0=Real_eigenvalue(folder_modes+"/sector_0/Measure.seq",pattern)

    susy_read_s0={}
    susy_read_s1={}
    for conf in conf_read:
        i=len(dictionary_s1[conf])-1
        gap_find=False
        susy_read_s1[conf]=i
        while i>0 and not gap_find:
            if (dictionary_s1[conf][i]-dictionary_s1[conf][i-1])>gap:
                susy_read_s1[conf]=i
                gap_find=True
            i-=1
        if not gap_find:
            susy_read_s1[conf]=1

        i=len(dictionary_s0[conf])-1
        print(i)
        susy_read_s0[conf]=i
        gap_find=False
        while i>=0 and not gap_find:
            if (dictionary_s0[conf][i]-dictionary_s0[conf][i-1])>gap:
                susy_read_s0[conf]=i
                gap_find=True
            i-=1
        if not gap_find:
            susy_read_s0[conf]=1
    print(susy_read_s0,susy_read_s1)

    return(susy_read_s0,susy_read_s1)
    
def Count_index_impr(folder,measurement,threshold,modes_used):
    count={}
    with open(folder) as file:
        for line in file:
            match=re.search(measurement,line)
            if match:
                string=line.split(":")
                for conf in configurations:
                    if (abs(float(string[8]))<threshold) and modes_used[string[7]] and string[1]==str(conf):
                        if string[1] in count:
                            count[string[1]]+=1
                        else:
                            count[string[1]]=1
                        break
    return(count)

def Round_half_integers(number):
    return(round(number*2)/2)

def Count_index_gf(folder,configurations):
    count_gauge={}
    conf_read=[]
    with open(folder+"Measure.seq") as file:
        for line in file:
            match=re.search(":TopologicalCharge:t:4",line)
            if match:
                string=line.split(":")
                for conf in configurations: #check that the configuration is on the list of configurations to be read
                    if string[1]==str(conf):
                        count_gauge[string[1]]=round(float(string[7]))
                        break
                        
    for conf in configurations:
        Topology =folder+"profile4dt4c"+str(conf)+"to.dat"
        density_top,sizes=Read.topology_1d(Topology)
        for element in density_top:
            if element>0.1 or element<-0.1:
                conf_read.append(str(conf))
                break
    return(count_gauge,conf_read)

def Topology_dic(folder,threshold,configurations):
    
    #GF + configurations with non_trivial topological content
    count_gauge,conf_read=Count_index_gf(folder,configurations)
    
    count_s0_susy=Count_index(folder+"sector_0/Measure.seq", ":OverlapFilterModeR:",threshold,conf_read)
    count_s1_susy=Count_index(folder+"sector_1/Measure.seq", ":OverlapFilterModeR:",threshold,conf_read)

    #Overlap
    count_s0_ov=Count_index(folder+"sector_0/Measure.seq", ":OverlapFilterModeC:",threshold,conf_read)
    count_s1_ov=Count_index(folder+"sector_1/Measure.seq", ":OverlapFilterModeC:",threshold,conf_read)

    #for key in count_gauge: print((count_s1_ov[key] + count_s0_ov[key]))
    ov_top_dif={}
    susy_top_dif={}
    for conf in conf_read:
        ov_top_dif[conf]=(count_s1_ov[conf] - count_s0_ov[conf])/4. - count_gauge[conf]
        susy_top_dif[conf]=(count_s1_susy[conf] - count_s0_susy[conf])/2. - count_gauge [conf]

    return(ov_top_dif,susy_top_dif,conf_read)
    
def Topology_dic_impr(folder,threshold,modes_used_s0,modes_used_s1):

    #Overlap
    count_s0_ov=Count_index(folder+"sector_0/Measure.seq", ":OverlapFilterModeC:",configurations,threshold)
    count_s1_ov=Count_index(folder+"sector_1/Measure.seq", ":OverlapFilterModeC:",configurations,threshold)
    count_s0_susy=Count_index_impr(folder+"sector_0/Measure.seq", ":OverlapFilterModeR:",configurations,threshold,modes_used_s0)
    count_s1_susy=Count_index_impr(folder+"sector_1/Measure.seq", ":OverlapFilterModeR:",configurations,threshold,modes_used_s1)
    #GF
    count_gauge={}
    with open(folder+"../gf/Measure.seq") as file:
        for line in file:
            match=re.search(":TopologicalCharge:t:4",line)
            if match:
                string=line.split(":")
                count_gauge[string[1]]=round((float(string[7])))
    for key in count_gauge:
        if key not in count_s0_ov: count_s0_ov[key]=0
        if key not in count_s1_ov: count_s1_ov[key]=0
        if key not in count_s0_susy: count_s0_susy[key]=0
        if key not in count_s1_susy: count_s1_susy[key]=0
    conf_tot=len(count_gauge)
    #for key in count_gauge: print((count_s1_ov[key] + count_s0_ov[key]))
    ov_top_dif={key: (count_s1_ov[key] - count_s0_ov[key])/4. - count_gauge[key] for key in count_gauge}
    susy_top_dif={key: (count_s1_susy[key] - count_s0_susy[key])/2. - count_gauge[key] for key in count_gauge}
    
    return(ov_top_dif,susy_top_dif,conf_tot)

def susy_plot(folder_in,folder_out,sizes,colors,spin_length,max_modes,lambda_opt,configurations):
    
    dictionary_s1=Real_eigenvalue(folder_in+"./sector_1/Measure.seq")
    dictionary_s0=Real_eigenvalue(folder_in+"./sector_0/Measure.seq")
    
    top_gauge,conf_read=Count_index_gf(folder_in,configurations)
    susy_read_s0=Count_index(folder_in+"sector_0/Measure.seq",":OverlapFilterModeR:",lambda_opt,conf_read)
    susy_read_s1=Count_index(folder_in+"sector_1/Measure.seq",":OverlapFilterModeR:",lambda_opt,conf_read)
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
        density_susy=Compare.Construct_susy(folder_in,susy_read_s0[conf],susy_read_s1[conf],
                                        conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes)
        density_susy=density_susy*(normalization/np.sum(np.abs(density_susy)))

        #Plot the three densities
        plt.plot(density_top_1, label=r'top. density $\tau=0.5$')
        plt.plot(density_top_2, label=r'top. density $\tau=2$')
        plt.plot(density_top_3, label=r'top. density $\tau=4$')
        plt.plot(density_susy, label="AFM")
        plt.legend(loc="lower left", ncol=2)
        plt.savefig(folder_out+"./susy_mode_"+conf+"c.png",dpi=150, bbox_inches='tight')
        plt.close()      
        np.savetxt(folder_out+"./susy_mode_"+conf+"c.txt",density_susy)
        
    return()

                                                                                                                                                                        Compare.py                                                                                          0000644 0001750 0001750 00000056533 14607720202 011705  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import sys

import Read
import importlib
import analyzer
import re
import pickle
plt.rcParams.update({'font.size': 16})

#plt.rcParams['text.usetex'] = True
                
    
def Construct_susy_density(folder,susy_read_s0,susy_read_s1,conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes): #It assumes eigenvalues are ordered
    
    #Read supersymmetric modes up to a threshold of the eigenvalue
    density_susy=np.zeros(sizes[3])
    for mode in range(0,susy_read_s1):
        if mode<max_modes:
            Mode = folder+"sector_1/SusyMode_bin_"+str(mode)+"-"+str(conf)
            density_s1,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
            density_susy+=density_s1/2
        
    for mode in range(0,susy_read_s0):
        if mode<max_modes:
            Mode = folder+"sector_0/SusyMode_bin_"+str(mode)+"-"+str(conf)
            density_s0,sizes=Read.bin_mode_1d(Mode,sizes,colors,spin_length)
            density_susy-=density_s0/2
    return(density_susy)  

def Construct_susy_modes(folder,susy_read_s0,susy_read_s1,conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes): #It assumes eigenvalues are ordered
    
    #Read supersymmetric modes up to a threshold of the eigenvalue
    susy_mode_s1=np.zeros((colors, spin_length, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
    susy_mode_s0=np.zeros((colors, spin_length, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
    for mode in range(0,susy_read_s1):
        if mode<max_modes:
            Mode = folder+"sector_1/SusyMode_bin_"+str(mode)+"-"+str(conf)
            zmode,density,sizes=Read.bin_mode(Mode,sizes,colors,spin_length)
            susy_mode_s1+=zmode
        
    for mode in range(0,susy_read_s0):
        if mode<max_modes:
            Mode = folder+"sector_0/SusyMode_bin_"+str(mode)+"-"+str(conf)
            zmode,density,sizes=Read.bin_mode(Mode,sizes,colors,spin_length)
            susy_mode_s0+=zmode

    density= Read.mode_to_density(susy_mode_s1,colors,spin_length,sizes)- Read.mode_to_density(susy_mode_s0,colors,spin_length,sizes)
    density_susy=density.sum(axis=(0,1,2))/2
    return(density_susy)

def Construct_susy_overlap_density(folder,susy_read_s0,susy_read_s1,conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes):
     #Read supersymmetric modes up to a threshold of the eigenvalue
    susy_mode_s1=np.zeros((colors, spin_length, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
    susy_mode_s0=np.zeros((colors, spin_length, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
    
    for mode in range(0,susy_read_s1):
        if mode<max_modes:
            Mode_f = folder+"sector_1/OverlapMode"+str(conf)+"_"+str(mode)
            Mode,density,sizes=Read.ascii_mode(Mode_f)
            susy_mode_s1+=Mode/(2*np.sqrt(2))
            
    for mode in range(0,susy_read_s0):
        if mode<max_modes:
            Mode_f = folder+"sector_0/OverlapMode"+str(conf)+"_"+str(mode)
            Mode,density,sizes=Read.ascii_mode(Mode_f)
            susy_mode_s0+=Mode/(2*np.sqrt(2))

    density_s1=Read.mode_to_density(susy_mode_s1,colors,spin_length,sizes)    
    density_s0=Read.mode_to_density(susy_mode_s0,colors,spin_length,sizes)
    
    density_susy=density_s1-density_s0
    return(density_susy)

def Construct_susy_overlap_modes(folder,susy_read_s0,susy_read_s1,conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes):
     #Read supersymmetric modes up to a threshold of the eigenvalue
    density_susy=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])

    for mode in range(0,susy_read_s1):
        if mode<max_modes:
            Mode_f = folder+"sector_1/OverlapMode"+str(conf)+"_"+str(mode)
            Mode,density,sizes=Read.ascii_mode(Mode_f)
            density_susy+=density/4 #/2 normalizing to q=1/2, /2 because modes come in pair with the same density
            
    for mode in range(0,susy_read_s0):
        if mode<max_modes:
            Mode_f = folder+"sector_0/OverlapMode"+str(conf)+"_"+str(mode)
            Mode,density,sizes=Read.ascii_mode(Mode_f)
            density_susy-=density/4 #/2 normalizing to q=1/2, /2 because modes come in pair with the same density
    
    density_1d=density_susy.sum(axis=(0,1,2))
    return(density_1d)
            
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

def GM_doublers(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf_read,tao_compare,susy_read_s0,susy_read_s1,pattern,save=True):
    dictionary_s1=analyzer.Real_eigenvalue(folder_in+"./sector_1/Measure.seq",pattern)
    dictionary_s0=analyzer.Real_eigenvalue(folder_in+"./sector_0/Measure.seq",pattern)
    
    GM=[]
    if pattern=="OverlapFilterModeR":
        for conf in conf_read:
            for i in range(0,max_modes):
                Mode_i = folder_in+"sector_1/SusyMode_bin_"+str(i)+"-"+str(conf)
                density_s1, sizes=Read.bin_mode_1d(Mode_i,sizes,colors,spin_length)
                for j in range(0,max_modes):
                    Mode_j = folder_in+"sector_0/SusyMode_bin_"+str(j)+"-"+str(conf)
                    density_s0, sizes=Read.bin_mode_1d(Mode_j,sizes,colors,spin_length)
                    GM.append([int(conf), i, j, Geom_mean_1d(density_s1,density_s0)])
    else:
        for conf in conf_read:
            for i in range(0,max_modes):
                susy_mode_s1=np.zeros((spin_length, colors, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
                Mode_i = folder_in+"sector_1/OverlapMode"+str(conf)+"_"+str(i)
                Mode_s1,density,sizes=Read.ascii_mode(Mode_i)
                susy_mode_s1+=Mode_s1/(2*np.sqrt(2))

                i+=1
                Mode_i = folder_in+"sector_1/OverlapMode"+str(conf)+"_"+str(i)
                Mode_s1,density,sizes=Read.ascii_mode(Mode_i)
                susy_mode_s1+=Mode_s1/(2*np.sqrt(2))

                density_s1=Read.mode_to_density(susy_mode_s1,colors,spin_length,sizes)
                
                for j in range(0,max_modes):
                    susy_mode_s0=np.zeros((spin_length, colors, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
                    Mode_j = folder_in+"sector_0/OverlapMode"+str(conf)+"_"+str(j)
                    Mode_s0,density,sizes=Read.ascii_mode(Mode_j)
                    susy_mode_s0+=Mode_s0/(2*np.sqrt(2))
                    
                    j+=1
                    Mode_j = folder_in+"sector_0/OverlapMode"+str(conf)+"_"+str(j)
                    Mode_s0,density,sizes=Read.ascii_mode(Mode_j)
                    susy_mode_s0+=Mode_s0/(2*np.sqrt(2))

                    density_s0=Read.mode_to_density(susy_mode_s0,colors,spin_length,sizes)    
                    GM.append([int(conf), i, j, Geom_mean_1d(density_s1,density_s0)])
    
                    
    with open(folder_out+"GM_doublers.txt", 'w') as f:
        for element in GM:
            print(element, file=f)
    
    return()

def GM_doublers_cut(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf_read,tao_compare,susy_read_s0,susy_read_s1,pattern,save=True):
    dictionary_s1=analyzer.Real_eigenvalue(folder_in+"./sector_1/Measure.seq",pattern)
    dictionary_s0=analyzer.Real_eigenvalue(folder_in+"./sector_0/Measure.seq",pattern)
    GM=[]
    if pattern=="OverlapFilterModeR":
        for conf in conf_read:
            for i in range(0,susy_read_s1[str(conf)]):
                Mode_i = folder_in+"sector_1/SusyMode_bin_"+str(i)+"-"+str(conf)
                density_s1, sizes=Read.bin_mode_1d(Mode_i,sizes,colors,spin_length)
                for j in range(0,susy_read_s0[str(conf)]):
                    Mode_j = folder_in+"sector_0/SusyMode_bin_"+str(j)+"-"+str(conf)
                    density_s0, sizes=Read.bin_mode_1d(Mode_j,sizes,colors,spin_length)
                    GM.append([int(conf), i, j, Geom_mean_1d(density_s1,density_s0)])
    else:
        for conf in conf_read:
            for i in range(0,susy_read_s1[str(conf)]):
                susy_mode_s1=np.zeros((spin_length, colors, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
                Mode_i = folder_in+"sector_1/OverlapMode"+str(conf)+"_"+str(i)
                Mode_s1,density,sizes=Read.ascii_mode(Mode_i)
                susy_mode_s1+=Mode_s1/(2*np.sqrt(2))

                i+=1
                Mode_i = folder_in+"sector_1/OverlapMode"+str(conf)+"_"+str(i)
                Mode_s1,density,sizes=Read.ascii_mode(Mode_i)
                susy_mode_s1+=Mode_s1/(2*np.sqrt(2))

                density_s1=Read.mode_to_density(susy_mode_s1,colors,spin_length,sizes)
                
                for j in range(0,susy_read_s0[str(conf)]):
                    susy_mode_s0=np.zeros((spin_length, colors, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
                    Mode_j = folder_in+"sector_0/OverlapMode"+str(conf)+"_"+str(j)
                    Mode_s0,density,sizes=Read.ascii_mode(Mode_j)
                    susy_mode_s0+=Mode_s0/(2*np.sqrt(2))
                    
                    j+=1
                    Mode_j = folder_in+"sector_0/OverlapMode"+str(conf)+"_"+str(j)
                    Mode_s0,density,sizes=Read.ascii_mode(Mode_j)
                    susy_mode_s0+=Mode_s0/(2*np.sqrt(2))

                    density_s0=Read.mode_to_density(susy_mode_s0,colors,spin_length,sizes)    
                    GM.append([int(conf), i, j, Geom_mean_1d(density_s1,density_s0)])
    
    with open(folder_out+"GM_doublers_cut.txt", 'w') as f:
        for element in GM:
            print(element, file=f)
    
    return()


    

def GM_RPO(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf_read,tao_compare,susy_read_s0,susy_read_s1,save=True):

    dictionary_s1=analyzer.Real_eigenvalue(folder_in+"./sector_1/Measure.seq")
    dictionary_s0=analyzer.Real_eigenvalue(folder_in+"./sector_0/Measure.seq")

    GM={}
    for conf in conf_read:
        #Read GF
        Topology =folder_in+"../gf/profile4dt"+str(tao_compare)+"c"+str(conf)+"to.dat"
        density_top,sizes=Read.topology_1d(Topology)

        #Construct susy mode
        density_susy=Construct_susy(folder_in,susy_read_s0[conf],susy_read_s1[conf],conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes)

        #Compute the distance between the susy and the topological density
        GM[conf]=Geom_mean_1d(density_susy,density_top)
        #RPO[conf]=Relative_point(np.absolute(density_susy),np.absolute(density_top),RPO_threshold)

    if save:
        with open(folder_out+"GM_hist_0.txt", 'w') as f:
            for key in GM:
                print(key,GM[key], file=f)
    #Store the means        
    GM_mean=0
    #RPO_mean=0
    count_meas=0
    for key in GM:
        GM_mean+=GM[key]
        count_meas+=1
    GM_mean=GM_mean/count_meas
    if not count_meas==len(conf_read):
        print("GM measures left out")
        
    #Saving the GM and RPO
    np.savetxt(folder_out+"GM.txt",np.vstack((GM_mean,count_meas)))
    return()         


def GM_RPO_cut(folder_in,folder_out,sizes,max_modes,colors,spin_length,conf_read,lambdas,RPO_threshold,tao_compare,pattern,method,modes):
    
    #Param and definitions
    steps=len(lambdas)
    GM_tot=np.zeros((steps))
    RPO_tot=np.zeros((steps))      
    param=np.zeros((steps))
    
    #Parsing the eigenvalues
    dictionary_s1=analyzer.Real_eigenvalue(folder_in+"./sector_1/Measure.seq",pattern)
    dictionary_s0=analyzer.Real_eigenvalue(folder_in+"./sector_0/Measure.seq",pattern)
    
    ov_max=0
    susy_max=0 #To know in which step of lambda one runs out of values
    k=0
    for threshold in lambdas:
        #Check how many modes for each configuration we need to read )
        if method=="cut":
            susy_read_s0,susy_read_s1,start_spectrum = analyzer.Count_index_cut(folder_in,"",threshold,conf_read,max_modes,pattern)
        elif method=="gap":
            susy_read_s0,susy_read_s1 = analyzer.Count_index_gap(folder_in,"",threshold,conf_read,max_modes,pattern)
        else:
            susy_read_s0,susy_read_s1=[],[]
            print("method:", method, " not found")
            exit()
        #Saving the number of modes
        diff={}
        for key in susy_read_s1:
            diff[key]=susy_read_s1[key]-susy_read_s0[key]

        modes_file=open(folder_out+"modes_used_"+str(k)+".txt", "w")
        for conf in conf_read:
            print(conf, susy_read_s1[conf], susy_read_s0[conf], diff[conf], file=modes_file)
        modes_file.close()
        GM={}
        j=0
        for conf in conf_read:
            #Read GF
            Topology =folder_in+"../gf/profile4dt"+str(tao_compare)+"c"+str(conf)+"to.dat"
            density_top,sizes=Read.topology_1d(Topology)
            
            #Construct susy mode
            if "OverlapFilterModeC" in pattern:
                if modes == "0" :
                    density_susy=Construct_susy_overlap_modes(folder_in,susy_read_s0[str(conf)],susy_read_s1[str(conf)],conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes)
                    print(modes)
                else:
                    density_susy=Construct_susy_overlap_density(folder_in,susy_read_s0[str(conf)],susy_read_s1[str(conf)],conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes)
                    print(modes)
            else:
                if modes == "0":
                    density_susy=Construct_susy_modes(folder_in,susy_read_s0[str(conf)],susy_read_s1[str(conf)],conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes)
                else:
                    density_susy=Construct_susy_density(folder_in,susy_read_s0[str(conf)],susy_read_s1[str(conf)],conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes)
            #Compute the distance between the susy and the topological density
            if susy_read_s0[conf]==0 and susy_read_s1[conf]==0:
                GM[conf]=0
            else:
                GM[conf]=Geom_mean_1d(density_susy,density_top)
        
        #Save the GM for histogram
        with open(folder_out+"GM_hist_"+str(k)+".txt", 'w') as f:
            for key in GM:
                print(key,GM[key], file=f)
                

        #Store the mean        
        GM_mean=0
        count_meas=0
        for key in GM:
            GM_mean+=GM[key]
            count_meas+=1  
        GM_mean=GM_mean/count_meas
        GM_tot[k]=GM_mean
        print(k) #to see the progress
        print(len(conf_read))

        k+=1
        
    #To check when one runs out of modes
    end_susy,end_overlap=analyzer.End_spectrum(folder_in,conf_read) #Computes the last value of the spectrum

    #Check which in which lambda one crosses the last value of the spectrum
    i=0
    susy_max=[lambdas[len(lambdas)-1], len(lambdas)] #The position of the end of the spectrum initilized at the end of the cuts
    for threshold in lambdas:
        if "OverlapFilterModeR" in pattern:
            if threshold > end_susy:
                susy_max=[threshold,i]
                break
        else :
            if threshold > end_overlap:
                susy_max=[threshold,i]
                break
        i+=1             
        
    #Saving the GM and RPO
    np.savetxt(folder_out+"end_spectrum.txt", np.array(susy_max)) #Save the position of the cut where the spectrum ends
    np.savetxt(folder_out+"GM.txt",np.vstack((GM_tot,lambdas))) #Save all GM means
    
    
    #Check end spectrum
    i=0
    start_spectrum=False
    if method=="cut":
        for threshold in lambdas:
            susy_read_s0,susy_read_s1,start_spectrum = analyzer.Count_index_cut(folder_in,"",threshold,conf_read,max_modes,pattern)
            if start_spectrum:
                np.savetxt(folder_out+"start_spectrum.txt", np.array([threshold,i]))
                break
                
        if not start_spectrum:
            np.savetxt(folder_out+"start_spectrum.txt", np.array([lambdas[len(lambdas)-1],len(lambdas)]))
        
    #finding the optimal lambda
    max_xi=0
    lambda_opt=np.zeros((2))
    k=0
    for threshold in lambdas:
        if threshold<=susy_max[0] and GM_tot[k]>max_xi:
            max_xi=GM_tot[k]
            lambda_opt[0]=threshold
            lambda_opt[1]=k
        k+=1
    np.savetxt(folder_out+"lambda_opt.txt",lambda_opt) #Save optimal lambda

    #Saving the number of modes used = instanton numbers at the optimal threshold
    if method=="cut":
        susy_read_s0,susy_read_s1,start_spectrum = analyzer.Count_index_cut(folder_in,"",lambda_opt[0],conf_read,max_modes,pattern)
    if method=="gap":
        susy_read_s0,susy_read_s1 = analyzer.Count_index_gap(folder_in,"",lambda_opt[0],conf_read,max_modes,pattern)
    modes_file=open(folder_out+"modes_used_opt.txt", "w")
    for conf in conf_read:
        print(conf, susy_read_s1[conf], susy_read_s0[conf], file=modes_file)
    modes_file.close()

    return()         


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
    #np.savetxt(folder+"./Index_ov.txt", np.vstack((ov_dif_count/conf_tot,param)))
    #np.savetxt(folder+"./Top_ov.txt", np.vstack((ov_dif_mean,param)))
    #np.savetxt(folder+"./Index_susy.txt",  np.vstack((susy_dif_count/conf_tot,param)))
    #np.savetxt(folder+"./Top_susy.txt", np.vstack((susy_dif_mean,param)))
    #np.savetxt(folder+"./Top_susy.txt", np.vstack((susy_dif_mean,param)))
    
    with open('used_conf.pkl', 'wb') as f:
        pickle.dump(susy, f)
    return()

def GF_vs_GF(folder_in,folder_out,conf_read, t_start,t_end,t_step,RPO_threshold,tao_compare):
    
    steps=int((t_end-t_start)/t_step)
    tot_conf=len(conf_read)
    GM_tot=np.zeros((steps+1))
    RPO_tot=np.zeros((steps+1))
    
    k=0
    for i in range(0, steps+1):
        t=t_start+t_step*i
        if int(t)==t: t=int(t)
        GM={}
        RPO={}
        tot_conf=0
        for conf in conf_read:
            Topology_smooth=folder_in+"profile4dt"+str(tao_compare)+"c"+str(conf)+"to.dat"
            density_smooth,sizes=Read.topology_1d(Topology_smooth)

            Topology_t=folder_in+"profile4dt"+str(t)+"c"+str(conf)+"to.dat"
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
        #print(tot_conf)
    param=np.arange(t_start, t_end+t_step, t_step)
    np.savetxt(folder_out+"./GM_GF.txt",np.vstack((GM_tot,param)))
    np.savetxt(folder_out+"./RPO_GF.txt", np.vstack((RPO_tot,param)))  
    
    return(param,GM_tot)


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
                                                                                                                                                                     Fitting.py                                                                                          0000644 0001750 0001750 00000053660 14607720207 011726  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import Read
import Maxima_find
import Plotting
import copy
plt.rcParams.update({'font.size': 12})

def instanton_parameter(directory,gf_zm,name,pref,guess_number_maxima,t_init,t_final,t_step,X1,rho1,height1):   

    Cap_number_maxima=guess_number_maxima
    
    #Variables to store the parameters
    Xmax_final=[[[0 for k in range(4)] for j in range(Cap_number_maxima)] for i in range(t_init,t_final,t_step)]
    height_final=[[[0 for k in range(4)] for j in range(Cap_number_maxima)] for i in range(t_init,t_final,t_step)]
    rho_final=[[[0 for k in range(4)] for j in range(Cap_number_maxima)] for i in range(t_init,t_final,t_step)]
    Xmax_previous=[[]]
    Xmax_disordered=[[0 for k in range(4)] for j in range(Cap_number_maxima)]
    height_disordered=[[0 for k in range(4)] for j in range(Cap_number_maxima)]
    rho_disordered=[[0 for k in range(4)] for j in range(Cap_number_maxima)]
    # Just tries to count the time as it appears in the file naming.
    t=0
    l=0
    res1=[[]]
    for i in range(t_init,t_final,t_step):
        time=str(i/10)
        if not (i%10):
            time=str(int(i/10))
        if (gf_zm):
            density,sizes=Read.density(directory+name+time+pref)
        else:
            density,sizes=Read.susy_mode(directory+name+time+pref)
        if ((density.sum())<0):
            density=-density
        res=copy.deepcopy(Maxima_find.simple(density,sizes)) #find the maxima
        #print(len(res[0]))
        #print(res)
        if (len(res[0])<=Cap_number_maxima and len(res[0])>0): # Cap to smooth enough configuration by the number of maxima
            rho=[]
            Xmax=[]
            height=[]
            for j in range(0,len(res[0])): # For every maxima we do a fit.
                Xmax_find=np.array([round(res[0][j]),round(res[1][j]),round(res[2][j]),round(res[3][j])])
                rho_temp=np.array([10,0,0,0])
                Xmax_temp=np.array([10,0,0,0])
                height_temp=np.array([10,0,0,0])
                for d in (0,1,2,3): # the fit is done independently in each direction
                    rho_temp[d], Xmax_temp[d], height_temp[d]=fitting_instanton(density,d,sizes,Xmax_find,0)
                #print(j,Xmax_temp)
                #Xmax_disordered[j]=copy.deepcopy(Xmax_temp)
                Xmax_disordered[j]=np.copy(Xmax_find)
                height_disordered[j]=np.copy(height_temp)
                rho_disordered[j]=np.copy(rho_temp)
            #print(i,Xmax_disordered)
            if (l==0 and not X1): #if it is not the first configuration and we do not have a set of maxima, we define it
                Xmax_final[t]=np.copy(Xmax_disordered)
                height_final[t]=np.copy(height_disordered)
                rho_final[t]=np.copy(rho_disordered)
                res1=res
            elif (l==0):
                order_index=order_maxima(len(res[0]),  [Xmax_disordered,height_disordered,rho_disordered], [X1,height1,rho1], len(X1), sizes)
                #print(order_index)
                for p in range(0, len(res[0])):
                    #print(order_index[p])
                    Xmax_final[t][order_index[p]]=Xmax_disordered[p]
                    height_final[t][order_index[p]]=height_disordered[p]
                    rho_final[t][order_index[p]]=rho_disordered[p]
                Xmax_previous=np.copy(Xmax_final[t])
                res1=res
            else:
                order_index=order_maxima(len(res[0]),  [Xmax_disordered,height_disordered,rho_disordered], [Xmax_final[t-1],height_final[t-1],rho_final[t-1]], len(res1[0]), sizes)
                for p in range(0, len(res[0])):
                    #print(order_index[p])
                    Xmax_final[t][order_index[p]]=Xmax_disordered[p]
                    height_final[t][order_index[p]]=height_disordered[p]
                    rho_final[t][order_index[p]]=rho_disordered[p]
                Xmax_previous=np.copy(Xmax_final[t])
                res1=res
            l=l+1
        t=t+1
    #print(t)
    #print(l)
    return(Xmax_final, rho_final, height_final)

def periodic_ind(shift, size):
    ind=shift%size
    return ind

def fitting_instanton(density, d, sizes, Xmax, plot): #needs the density, the direction of the section d, the position of the maxima Xmax, returns the coefficient in front of x^2
    #1-D parabola
    def func(x, h, Xmax_extr, rho):
        return h - 1/(rho**2)*(x-Xmax_extr)**2
    
    # pick up one maxima
    #Xmax_int=[res[0][0]/5,res[1][0]/5,res[2][0]/5,res[3][0]/5]
    #Xmax=[round(Xmax_int[0]),round(Xmax_int[1]),round(Xmax_int[2]),round(Xmax_int[3])]
    
    # prepare the data: we will use 3 points, the Xmax and the nearest neighbours
    x=np.array([Xmax[d]-1,Xmax[d],Xmax[d]+1])
    y=np.array([Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4]])
    z=np.array([Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4]])
    t=np.array([Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4]])
    shift= { d:x, (d+1)%4:y, (d+2)%4:z, (d+3)%4:t} # dictionary from integer to axis
    data=density[periodic_ind(shift[0],sizes[0]),periodic_ind(shift[1],sizes[1]),
                     periodic_ind(shift[2],sizes[2]),periodic_ind(shift[3],sizes[3])]
    
    try:
        popt,pcov= curve_fit(func, x, data)
        
    except RuntimeError:
        print("error")
        return [0, 0, 0]
    rho=popt[2]
    Xmax_extr=popt[1]
    height=popt[0]
    
    #xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
    #plt.plot(xdata, func(xdata, *popt))
    #plt.plot(x, data)
    #plt.show()
    
    if plot:
        plt.plot(x, data)
        x=np.array([Xmax[d]-2,Xmax[d]-1, Xmax[d],Xmax[d]+1, Xmax[d]+2])
        y=np.array([Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4]])
        z=np.array([Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4]])
        t=np.array([Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4]])
        shift= { d:x, (d+1)%4:y, (d+2)%4:z, (d+3)%4:t} # dictionary from integer to axis
        data=density[periodic_ind(shift[0],sizes[0]),periodic_ind(shift[1],sizes[1]),
                     periodic_ind(shift[2],sizes[2]),periodic_ind(shift[3],sizes[3])]
        
        xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
        plt.plot(xdata, func(xdata, *popt))
        #plt.show()
        
    return popt

def fitting_instanton_impr(density, d, sizes, Xmax, plot): #needs the density, the direction of the section d, the position of the maxima Xmax, returns the coefficient in front of x^2
    #1-D parabola
    def func(x, h, Xmax_extr, rho):
        return h - 1/(rho**2)*(x-Xmax_extr)**2
    
    # prepare the data: we will use 5 points, the Xmax and the nearest neighbours

    x=np.array([Xmax[d]-2,Xmax[d]-1, Xmax[d],Xmax[d]+1, Xmax[d]+2])
    y=np.array([Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4]])
    z=np.array([Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4]])
    t=np.array([Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4]])
    shift= { d:x, (d+1)%4:y, (d+2)%4:z, (d+3)%4:t} # dictionary from integer to axis
    data=density[periodic_ind(shift[0],sizes[0]),periodic_ind(shift[1],sizes[1]),
                 periodic_ind(shift[2],sizes[2]),periodic_ind(shift[3],sizes[3])]

    xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
    #plt.plot(xdata, func(xdata, *popt))
    #plt.plot(x, data)
    #plt.show()
    try:
        popt, pcov = curve_fit(func, x, data)
        #popt = curve_fit(func, x, data)
    except RuntimeError:
        print("error")
        return [0, 0, 0], [0, 0, 0]
        #return 0,0,0,[0, 0, 0]
    rho=popt[2]
    Xmax_extr=popt[1]
    height= popt[0]
    
    #xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
    #plt.plot(xdata, func(xdata, *popt))
    #plt.plot(x, data)
    #plt.show()
    
    if plot:
        x=np.array([Xmax[d]-2,Xmax[d]-1, Xmax[d],Xmax[d]+1, Xmax[d]+2])
        y=np.array([Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4],Xmax[(d+1)%4]])
        z=np.array([Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4],Xmax[(d+2)%4]])
        t=np.array([Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4],Xmax[(d+3)%4]])
        shift= { d:x, (d+1)%4:y, (d+2)%4:z, (d+3)%4:t} # dictionary from integer to axis
        data=density[periodic_ind(shift[0],sizes[0]),periodic_ind(shift[1],sizes[1]),
                     periodic_ind(shift[2],sizes[2]),periodic_ind(shift[3],sizes[3])]
        
        xdata = np.linspace(Xmax[d]-2, Xmax[d]+2, 50)
        plt.plot(xdata, func(xdata, *popt))
        plt.plot(x, data)
        #plt.show()
    return (popt, pcov)

    #return (height, Xmax_extr, rho, pcov)

def instanton_parameter_np(directory,gf_zm,name,pref,guess_number_maxima,t_init,t_final,t_step,X1,rho1,height1):   

    Cap_number_maxima=guess_number_maxima*10
    #Variables to store the parameters
    Xmax_final=np.zeros((int((t_final-t_init)/t_step),Cap_number_maxima, 4))
    height_final=np.zeros((int((t_final-t_init)/t_step),Cap_number_maxima, 4))
    rho_final=np.zeros((int((t_final-t_init)/t_step),Cap_number_maxima, 4))
    Xmax_previous=np.zeros((Cap_number_maxima, 4))
    Xmax_disordered=np.zeros((Cap_number_maxima, 4))
    height_disordered=np.zeros((Cap_number_maxima, 4))
    rho_disordered=np.zeros((Cap_number_maxima, 4))
    # Just tries to count the time as it appears in the file naming.
    t=0
    res1=[[]]
    for i in range(t_init,t_final,t_step):
        time=str(i/10)
        if not (i%10):
            time=str(int(i/10))
        if (gf_zm):
            density,sizes=Read.density(directory+name+time+pref)
            density*=-density
        else:
            density,sizes=Read.susy_mode(directory+name+time+pre)
            density*=-density
        if ((density.sum())<0):
            density=-density
        res=Maxima_find.simple(density,sizes) #find the maxima
        #print(len(res[0]))
        #print(res)
        if (len(res[0])>0): # Just if there is some maxima otherwise it breaks
            rho=[]
            Xmax=[]
            height=[]
            for j in range(0,len(res[0])): # For every maxima we do a fit.
                Xmax_find=[round(res[0][j]),round(res[1][j]),round(res[2][j]),round(res[3][j])]
                #print(len(res[0]))
                rho_temp=[10,0,0,0]
                Xmax_temp=[10,0,0,0]
                height_temp=[10,0,0,0]
                for d in (0,1,2,3): # the fit is done independently in each direction
                    rho_temp[d], Xmax_temp[d], height_temp[d]=fitting_instanton(density,d,sizes,Xmax_find,0)
                Xmax_disordered[j]=np.copy(Xmax_temp)
                height_disordered[j]=np.copy(height_temp)
                rho_disordered[j]=np.copy(rho_temp)
            if (t==0 and not np.size(X1)):
                Xmax_final[t]=np.copy(Xmax_disordered)
                height_final[t]=np.copy(height_disordered)
                rho_final[t]=np.copy(rho_disordered)
                res1=np.copy(res)
            elif (t==0 and np.size(X1)):
                order_index=order_maxima(len(res[0]),  [Xmax_disordered,height_disordered,rho_disordered], [X1,height1,rho1], len(X1), sizes)
                #print(order_index)
                for p in range(0, len(res[0])):
                    #print(order_index[p])
                    Xmax_final[t][order_index[p]]=np.copy(Xmax_disordered[p])
                    height_final[t][order_index[p]]=np.copy(height_disordered[p])
                    rho_final[t][order_index[p]]=np.copy(rho_disordered[p])
                Xmax_previous=np.copy(Xmax_final[t])
                res1=res
            else:
                order_index=order_maxima(len(res[0]),  [Xmax_disordered,height_disordered,rho_disordered], [Xmax_final[t-1],height_final[t-1],rho_final[t-1]], len(res1[0]), sizes)
                for p in range(0, len(res[0])):
                    #print(order_index[p])
                    Xmax_final[t][order_index[p]]=np.copy(Xmax_disordered[p])
                    height_final[t][order_index[p]]=np.copy(height_disordered[p])
                    rho_final[t][order_index[p]]=np.copy(rho_disordered[p])
                Xmax_previous=np.copy(Xmax_final[t])
                res1=res
            t=t+1
        #print(t)
            #print(res)
            #print()
    return(Xmax_final, rho_final, height_final)
 

def distance_param(param_in, param_out,sizes):
    #print(param_in)
    return((np.abs(param_in[0]-param_out[0]))%(sizes[0]) + (np.abs(param_in[1]-param_out[1]))%(sizes[1]) + (np.abs(param_in[2]-param_out[2]))%(sizes[2]) + (np.abs(param_in[3]-param_out[3]))%(sizes[3]))

def distance_param2(param_in, param_out,sizes):
    return(np.abs(param_in[0]-param_out[0])+ np.abs(param_in[1]-param_out[1]) + np.abs(param_in[2]-param_out[2]) + np.abs(param_in[3]-param_out[3]))

def order_maxima(num_max_1, param_disordered, param_previous, num_max_2, sizes):
    #print(num_max_1, num_max_2)
    tmp_max_1=max(num_max_1,num_max_2)
    Max_range=np.arange((tmp_max_1))
    tmp_max_2=min(num_max_1,num_max_2)
    index=[[]]*tmp_max_2
    blocked_list=[]
    #print(len(index),np.size(Max_range))
    #print(num_max_1, num_max_2)
    for p in range(0,tmp_max_2): # We compare with the previous parameters to identify the maxima and see how the parameters change
        distance=1000000
        for l in range(0,tmp_max_1):
            distance_temp=distance_param(param_disordered[0][p],param_previous[0][l],sizes)#+distance_param2(param_disordered[1][p],param_previous[1][l],sizes)+distance_param2(param_disordered[2][p],param_previous[2][l],sizes)
            if (distance_temp<=distance and l not in blocked_list): #and Xmax_previous[l][0]>0 and Xmax_disordered[p][0]>0): # if the point is closer we substitute it
                index[p]=l #We store the position on the matrix 
                distance=distance_temp
        blocked_list.append(index[p])
    # We correct for the missing ones
    if (tmp_max_1==num_max_1):
        missing=list((set(Max_range).difference(index)))
        joined_index=index+missing
    else:
        joined_index=index
    #print(joined_index,'\n')
    #if (num_max_1>num_max_2):
    #    index=[0]*num_max_1
    #    print("bla")
    #if (num_max_1==0):  
    #    index=[0]
    #    print("ble")
    #print(index)
    return joined_index

def instanton_evolution_new():
    title=" "
    directory="./q0/gf_unbug/"
    name="profile4dt"
    legend4='wilson'
    pref="c100to.dat"
    max_guess=8
    gf_zm=True
    t_init=10
    t_final=14
    t_steps=4
    for t in range(t_init,t_final,t_steps):
        file=directory+name+str(t)+pref
        density,sizes=Read.density(file)
        density=np.abs(density)
        res=Maxima_find.improve(density,sizes)
        print(res)
        for i in range(0,len(res[0]),1):
            for d in range(0,3,1):
                Xmax_find=[round(res[0][i]),round(res[1][i]),round(res[2][i]),round(res[3][i])]
                rho,Xmax,height=fitting_instanton(density,d,sizes,Xmax_find,0)
    print("finished")
    return(rho,Xmax,height)

def instanton_evolution():
    title=" "
    
    directory="./notwist/wilson/"
    #directory="./heated/heated25/epsilon_1/"
    #directory="./new/epsilon_1/"
    name="profile4dt"
    legend4='wilson'
    pref="c100to.dat"
    max_guess=8
    gf_zm=True
    t_init=40
    t_final=360
    t_steps=4
    X4,rho4,height4=instanton_parameter_np(directory,gf_zm,name,pref,max_guess,t_init,t_final,t_steps,[],[],[])
    print("finished")
    
    #directory="./notwist/wilson/"
    directory="./heated/heated25/wilson/"
    #directory="./new/wilson/"
    legend1='wilson'
    name="profile4dt"
    pref="c100to.dat"

    X1,rho1,height1=instanton_parameter_np(directory,gf_zm,name,pref,max_guess,t_init,t_final,t_steps,[],[],[])
    print("finished")
    
    #directory="./notwist/symanzik/"
    directory="./heated/heated25/symanzik/"
    #directory="./new/symanzik/"
    name="profile4dt"
    legend2='symanzik'
    pref="c100to.dat"
    gf_zm=True
    X2,rho2,height2=instanton_parameter_np(directory,gf_zm,name,pref,max_guess,t_init,t_final,t_steps,[],[],[])
    
    #directory="./notwist/overimproved/"
    directory="./heated/heated25/overimproved/"
    #directory="./new/overimproved/"
    name="profile4dt"
    legend3='too_overimproved'
    pref="c100to.dat"
    gf_zm=True
    X3,rho3,height3=instanton_parameter_np(directory,gf_zm,name,pref,max_guess,t_init,t_final,t_steps,[],[],[])
    
    t_points=[]
    for i in range(t_init,t_final,t_steps):
        t_points.append(np.sqrt(8*i/100))
    for maxima_number in range(0,max_guess):
        f=plt.figure(figsize=(20,15))
        for direction in range(1,2):
            X1_plot=[]
            X2_plot=[]
            X3_plot=[]
            X4_plot=[]
            rho1_plot=[]
            rho2_plot=[]
            rho3_plot=[]
            rho4_plot=[]
            height1_plot=[]
            height2_plot=[]
            height3_plot=[]
            height4_plot=[]

            for t in range(0,len(X1)):
                try:
                    X1_plot.append(X1[t][maxima_number][direction])
                    X2_plot.append(X2[t][maxima_number][direction])
                    rho1_plot.append(rho1[t][maxima_number][direction])
                    rho2_plot.append(rho2[t][maxima_number][direction])
                    height1_plot.append(height1[t][maxima_number][direction])
                    height2_plot.append(height2[t][maxima_number][direction])
                    
                    X3_plot.append(X3[t][maxima_number][direction])
                    rho3_plot.append(rho3[t][maxima_number][direction])
                    height3_plot.append(height3[t][maxima_number][direction])
                    
                    X4_plot.append(X4[t][maxima_number][direction])
                    rho4_plot.append(rho4[t][maxima_number][direction])
                    height4_plot.append(height4[t][maxima_number][direction])
                except IndexError:
                    break 
                    
            
            #.plot_list(t_points,X1_plot,f, "", 
            #                   "flow_time","X",legend1,1,1,"blue",True,True,3,4,direction+1)
            #Plotting.plot_list(t_points,X2_plot,f, "", 
            #                   "flow_time","X",legend2,1,1,"red",True,True,3,4,direction+1)
            #Plotting.plot_list(t_points,X3_plot,f, "", 
            #                   "flow_time","X",legend3,1,1,"green",True,True,3,4,direction+1)
            plt.legend(loc="lower left", prop={'size': 20})
            Plotting.plot_list(t_points,rho1_plot,f, "", 
                               "","",legend1,1,1,"blue",True,True,2,1,2)
            Plotting.plot_list(t_points,rho2_plot,f, "", 
                               "","",legend2,1,1,"red",True,True,2,1,2)
            Plotting.plot_list(t_points,rho4_plot,f, "", 
                               "","",legend4,1,1,"green",True,True,2,1,2)

            #Plotting.plot_list(t_points,height1_plot,f, "",
            #                   "flow_time","height",legend1,1,1,"blue",True,True,3,4,direction+1+8)
            #Plotting.plot_list(t_points,height2_plot,f, "",
            #                   "flow_time","height",legend2,1,1,"red",True,True,3,4,direction+1+8)
            #Plotting.plot_list(t_points,height3_plot,f, "",
            #                   "flow_time","height",legend3,1,1,"green",True,True,3,4,direction+1+8)
        #plt.savefig("Maxima"+str(maxima_number))
        
            
            #Plotting.plot_list(t_points,rho1_plot,f,"", 
            #                   "flow_time","rho",legend1,1,1,"blue",True,True,3,1,2)
            #Plotting.plot_list(t_points,rho2_plot,f, "", 
            #Plotting.plot_list(t_points,rho3_plot,f, "", 
            #                   "flow_time","rho",legend3,1,1,"green",True,True,3,1,2)

            plt.legend(loc="lower left", prop={'size': 20})
            Plotting.plot_list(t_points,height1_plot,f, "",
                               "","",legend1,1,1,"blue",True,True,2,1,1)
            Plotting.plot_list(t_points,height2_plot,f, "",
                               "","",legend2,1,1,"red",True,True,2,1,1)
            #Plotting.plot_list(t_points,height3_plot,f, "",
            #                   "flow_time","height",legend3,1,1,"green",True,True,3,1,2)
            Plotting.plot_list(t_points,height4_plot,f, "",
                               "","",legend4,1,1,"purple",True,True,2,1,1)
            plt.legend(loc="lower left", prop={'size': 20})
            
            #Plotting.plot_list(t_points,X1_plot,f, title, 
            #                   "flow_time","X_0",legend1,1,1,"blue",True,True,3,1,1)
            #Plotting.plot_list(t_points,X2_plot,f, title, 
            #                   "flow_time","X_0",legend2,1,1,"red",True,True,3,1,1)
           # Plotting.plot_list(t_points,X3_plot,f, title, 
           #                    "flow_time","X_0",legend3,1,1,"green",True,True,3,1,1)
            #Plotting.plot_list(t_points,X4_plot,f, title, 
             #                  "flow_time","X_0",legend4,1,1,"purple",True,True,3,1,1)
            plt.legend(loc="lower left", prop={'size': 20})
            plt.savefig(title+"_Maxima"+str(maxima_number)+".png",dpi=100)
            plt.show()
            
    count1=max_count(height1)
    count2=max_count(height2)
    #count3=max_count(height3)
    count4=max_count(height4)
    
    Plotting.plot_list(t_points,count1,f, "",
                               "flow_time","#maxima",legend1,1,1,"blue",True,True,1,1,1)
    Plotting.plot_list(t_points,count2,f, "",
                               "flow_time","#maxima",legend2,1,1,"red",True,True,1,1,1)
   # Plotting.plot_list(t_points,count3,f, "",
    #                           "flow_time","#maxima",legend3,1,1,"green",True,True,1,1,1)
    Plotting.plot_list(t_points,count4,f, "",
                               "flow_time","#maxima",legend4,1,1,"purple",True,True,1,1,1)
    plt.rcParams.update({'font.size': 12})
    lt.legend(loc="lower left", prop={'size': 12})
    plt.savefig("Maxima_evolution.pdf")
    
    return(height1, height2, height3, height4)

def max_count(height):
    count=np.zeros((np.shape(height)[0]))
    for i in range(0, np.shape(height)[0]):
            for j in range(0, np.shape(height)[1]):
                if (height[i][j][0]>1e-12):
                    count[i]+=1
    return(count)                                                                                Loop_AFM.py                                                                                         0000644 0001750 0001750 00000001672 14513543615 011714  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import sys
import Plot_generator
plt.rcParams.update({'font.size': 16})
plt.rcParams['text.usetex'] = False

folder_in=str(sys.argv[1]) #./gf_afm_1p5t/ 
tau_compare=str(sys.argv[2]) # 1.5
folder_out=str(sys.argv[3]) #./compare_1p5t/gf_afm_1p5t/  

sizes=[4,4,4,32]
max_modes=8
colors=3
spin_length=4
max_modes=8

folder_gf="/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/gf/"
top_gauge,conf_read=analyzer.Count_index_gf(folder_gf,configurations)

tony_folder="../cases/"
tony_data=np.loadtxt(tony_folder+"tony_instantons.txt")
susy_read_s1={}
susy_read_s0={}
for element in tony_data:
    susy_read_s1[element[2]]=element[0]
    susy_read_s0[element[2]]=element[1]
Plot_generator.susy_plot("/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/"+str(sys.argv[1]),folder_out,sizes,colors,spin_length,max_modes,conf_read,
                         susy_read_s0,susy_read_s1,Load=False,Plot=True)

                                                                      main_fit.py                                                                                         0000644 0001750 0001750 00000004437 14607720210 012100  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   mport matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d, Axes3D
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interpn
import Projecting
import Plotting
import Read
import Maxima_find
import Fitting
import subprocess
import math
import re
from subprocess import call

plt.rcParams.update({'font.size': 12})


def fit_impr(folder,file,zm_top,normalization):
    
    if zm_top:
        density,sizes=Read.ascii_mode(file)
        title="susy"
    else:
        density,sizes=Read.topology(file)
        title="top"
    res=Maxima_find.simple(density,sizes)
    density=density*4*np.pi*np.pi*len(res[0])*normalization
    Plotting.plot_all_peaks(folder,density, sizes, title)

    param_all=[[]]
    param_all.pop(0)
    param_temp=[[]]
    param_temp.pop(0)

    err_all=[[]]
    err_all.pop(0)
    err_temp=[[]]
    err_temp.pop(0)
    for i in range(0,len(res[0])):
        param_temp=[[]]
        param_temp.pop(0)
        for d in range(0,4):
            Xmax=[res[0,i],res[1,i],res[2,i],res[3,i]]
            popt_temp,pcov_temp=Fitting.fitting_instanton_impr(density, d, sizes, Xmax, "0") #h,x,rho
            param_temp.append(popt_temp)
            err_temp.append(np.sqrt(np.diag(pcov_temp)))
        param_all.append(param_temp)
        err_all.append(err_temp)

    param=np.array(param_all)
    err=np.array(err_all)
 
    return density, param, err 

def fit_simple(folder,file,zm_top,normalization):
    if zm_top:
        density,sizes=Read.ascii_mode(file)
        res=Maxima_find.simple(density,sizes)
        density=density*4*np.pi*np.pi*len(res[0])
        title="susy"
    else:
        density,sizes=Read.topology(file)
        res=Maxima_find.simple(density,sizes)
        title="top"
    
    Plotting.plot_all_peaks(folder,density, sizes, title)

    param_all=[[]]
    param_all.pop(0)
    param_temp=[[]]
    param_temp.pop(0)

    for i in range(0,len(res[0])):
        param_temp=[[]]
        param_temp.pop(0)
        for d in range(0,4):
            Xmax=[res[0,i],res[1,i],res[2,i],res[3,i]]
            popt_temp=Fitting.fitting_instanton(density, d, sizes, Xmax, "0")
            param_temp.append(popt_temp)
        param_all.append(param_temp)
    param=np.array(param_all)
    return density, param                                                                                                                                                                                                                                 main_plot_polyakov.py                                                                               0000600 0001750 0001750 00000003706 14607720210 014206  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
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

                                                          main_plot.py                                                                                        0000600 0001750 0001750 00000007343 14607720165 012274  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import sys
import Plot_generator
plt.rcParams.update({'font.size': 16})
plt.rcParams['text.usetex'] = True

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
t_end=float(tau_compare)

t_step=0.25
print("GF vs AFM")
print(folder_out)
#Plot_generator.GF_vs_AFM(folder_out, folder_gf, folder_out, conf_read, t_start, t_end, t_step,
#                         RPO_threshold,tau_compare,measures,time_measures,observable)

'''for measure in measures:
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

    cut_i=0
    for cut in lambdas:
        print("susy modes at ", str(cut))
        lambda_opt,index_opt=float(cut), cut_i
        if method=="cut":
            susy_read_s0,susy_read_s1,start_spectrum = analyzer.Count_index_cut("./"+measure,"",lambda_opt,conf_read,max_modes,operator)
        if method=="gap":
            susy_read_s0,susy_read_s1 = analyzer.Count_index_gap("./"+measure,"",lambda_opt,conf_read,max_modes,operator)
            Plot_generator.susy_plot(folder_modes+measure,folder_out+measure,sizes,colors,
                                     spin_length,max_modes,conf_read,susy_read_s0,susy_read_s1,operator,str(cut_i))
        cut_i+=1'''

    
                                                                                                                                                                                                                                                                                             main_polyakov.py                                                                                    0000600 0001750 0001750 00000002052 14607720165 013152  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
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

conf_start=10
conf_end=1000
conf_step=10
configurations=np.arange(conf_start,conf_end,conf_step)
folder_gf="/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/gf/"

top_gauge,conf_read=analyzer.Count_index_gf(folder_gf,configurations)
tony_folder="/home/mi37fud/b2p44_new/cases/"
tony_data=np.loadtxt(tony_folder+"tony_instantons.txt",dtype=int)
susy_read_s1={}
susy_read_s0={}

for element in tony_data:
    susy_read_s1[str(element[2])]=element[0]
    susy_read_s0[str(element[2])]=element[1]

Compare.GM_RPO("/beegfs/mi37fud/4x4x4x32_su2/b2p44_new/"+str(sys.argv[1]),folder_out,sizes,max_modes,colors,spin_length,conf_read,tau_compare,susy_read_s0,susy_read_s1)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      main.py                                                                                             0000644 0001750 0001750 00000004230 14607720165 011236  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
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
                                                                                                                                                                                                                                                                                                                                                                        Maxima.py                                                                                           0000644 0001750 0001750 00000025540 14607720210 011524  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import numpy as np

def simple(density,sizes):
    maxima=np.zeros((4, ), dtype=int)
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            for k in range(0,sizes[2]):
                for l in range(0,sizes[3]):
                    #First the nearest neighbours
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],k,l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,k,(l-1)%sizes[3]])):
                        continue
                    maxima=np.column_stack((maxima, [i,j,k,l]))
    maxima=np.delete(maxima, 0, 1)
    return(maxima)

def simple_2d(density,sizes):
    maxima=np.zeros((2, ), dtype=int)
    for i in range(0,sizes[0]):
        for j in range(0,sizes[3]):
            #First the nearest neighbours
            if ((density[i,j] < density [(i+1)%sizes[0],j]) or (density[i,j] < density [(i-1)%sizes[0],j])):
                continue
            if ((density[i,j] < density [i,(j+1)%sizes[3]]) or (density[i,j] < density [i,(j-1)%sizes[3]])):
                continue
            maxima=np.column_stack((maxima, [i,j]))
    maxima=np.delete(maxima, 0, 1)
    return(maxima)


def simple_1d(density,sizes):
    maxima=np.zeros((1, ), dtype=int)
    for i in range(0,sizes):
            #First the nearest neighbours
            if ((density[i] < density [(i+1)%sizes]) or (density[i] < density [(i-1)%sizes])):
                continue
            maxima=np.column_stack((maxima, [i]))
    if (len(maxima)>0):
        maxima=np.delete(maxima, 0)
    return(maxima)

def improve_1d(density,sizes):
    maxima=np.zeros((1, ), dtype=int)
    for i in range(0,sizes):
            #First the nearest neighbours
            if ((density[i] < density [(i+1)%sizes]) or (density[i] < density [(i-1)%sizes]) or (density[i] < density [(i+2)%sizes]) or (density[i] < density [(i-2)%sizes])):
                continue
            maxima=np.column_stack((maxima, [i]))
    if (len(maxima)>0):
        maxima=np.delete(maxima, 0)
    return(maxima)

def improve_2d(density,sizes):
    maxima=np.zeros((2, ), dtype=int)
    for i in range(0,sizes[0]):
        for j in range(0,sizes[3]):
                    #First the nearest neighbours
                    if ((density[i,j] < density [(i+1)%sizes[0],j]) or (density[i,j] < density [(i-1)%sizes[0],j])):
                        continue
                    if ((density[i,j] < density [i,(j+1)%sizes[3]]) or (density[i,j] < density [i,(j-1)%sizes[3]])):
                        continue
                    #Now the diagonals of i j
                    if ((density[i,j] < density [(i+1)%sizes[0],(j+1)%sizes[1]]) or (density[i,j] < density [(i-1)%sizes[0],(j-1)%sizes[1]])) or (density[i,j] < density [(i+1)%sizes[0],(j-1)%sizes[1]]) or (density[i,j] < density [(i-1)%sizes[0],(j+1)%sizes[1]]):
                        continue
                    maxima=np.column_stack((maxima, [i,j]))

    maxima=np.delete(maxima, 0, 1)
    return(maxima)
                   
def improve(density,sizes):
    maxima=np.zeros((4, ), dtype=int)
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            for k in range(0,sizes[2]):
                for l in range(0,sizes[3]):
                    #First the nearest neighbours
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],k,l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,k,(l-1)%sizes[3]])):
                        continue
                    #Now the diagonals of i j
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],(j+1)%sizes[1],k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],(j-1)%sizes[1],k,l])) or (density[i,j,k,l] < density [(i+1)%sizes[0],(j-1)%sizes[1],k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],(j+1)%sizes[1],k,l]):
                        continue
                    #Now the diagonals of i k    
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,(k-1)%sizes[2],l])) or (density[i,j,k,l] < density [(i+1)%sizes[0],j,(k-1)%sizes[2],l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,(k+1)%sizes[2],l]):
                        continue
                    #Now the diagonals of i l           
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,(l-1)%sizes[3]])) or (density[i,j,k,l] < density [(i+1)%sizes[0],j,k,(l-1)%sizes[3]]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,(l+1)%sizes[3]]):
                        continue
                    #Now the diagonals of j k     
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],(k-1)%sizes[2],l])) or (density[i,j,k,l] < density [i,(j+1)%sizes[1],(k-1)%sizes[2],l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],(k+1)%sizes[2],l]):
                    #Now the diagonals j l
                        continue
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,(l-1)%sizes[3]])) or (density[i,j,k,l] < density [i,(j+1)%sizes[1],k,(l-1)%sizes[3]]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,(l+1)%sizes[3]]):
                        continue
                    #Now the diagonals k l
                    if ((density[i,j,k,l] < density [i,j,(k+1)%sizes[2],(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],(l-1)%sizes[3]])) or (density[i,j,k,l] < density [i,j,(k+1)%sizes[2],(l-1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],(l+1)%sizes[3]]):
                        continue
                    maxima=np.column_stack((maxima, [i,j,k,l]))

    maxima=np.delete(maxima, 0, 1)
    return(maxima)

def nnneighbours(density,sizes):
    maxima=np.zeros((4, ), dtype=int)
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            for k in range(0,sizes[2]):
                for l in range(0,sizes[3]):
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,l]) or (density[i,j,k,l] < density [(i+2)%sizes[0],j,k,l]) or (density[i,j,k,l] < density [(i-2)%sizes[0],j,k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],k,l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,l]) or (density[i,j,k,l] < density [i,(j+2)%sizes[1],k,l]) or (density[i,j,k,l] < density [i,(j-2)%sizes[1],k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],l]) or (density[i,j,k,l] < density [i,j,(k+2)%sizes[2],l]) or (density[i,j,k,l] < density [i,j,(k-2)%sizes[2],l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,k,(l-1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,k,(l+2)%sizes[3]]) or (density[i,j,k,l] < density [i,j,k,(l-2)%sizes[3]])):
                        continue
                    #Now the diagonals of i j
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],(j+1)%sizes[1],k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],(j-1)%sizes[1],k,l])) or (density[i,j,k,l] < density [(i+1)%sizes[0],(j-1)%sizes[1],k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],(j+1)%sizes[1],k,l]):
                        continue
                    #Now the diagonals of i k    
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,(k-1)%sizes[2],l])) or (density[i,j,k,l] < density [(i+1)%sizes[0],j,(k-1)%sizes[2],l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,(k+1)%sizes[2],l]):
                        continue
                    #Now the diagonals of i l           
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,(l-1)%sizes[3]])) or (density[i,j,k,l] < density [(i+1)%sizes[0],j,k,(l-1)%sizes[3]]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,(l+1)%sizes[3]]):
                        continue
                    #Now the diagonals of j k     
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],(k-1)%sizes[2],l])) or (density[i,j,k,l] < density [i,(j+1)%sizes[1],(k-1)%sizes[2],l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],(k+1)%sizes[2],l]):
                    #Now the diagonals j l
                        continue
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,(l-1)%sizes[3]])) or (density[i,j,k,l] < density [i,(j+1)%sizes[1],k,(l-1)%sizes[3]]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,(l+1)%sizes[3]]):
                        continue
                    #Now the diagonals k l
                    if ((density[i,j,k,l] < density [i,j,(k+1)%sizes[2],(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],(l-1)%sizes[3]])) or (density[i,j,k,l] < density [i,j,(k+1)%sizes[2],(l-1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],(l+1)%sizes[3]]):
                        continue
                    maxima=np.column_stack((maxima, [i,j,k,l]))
    #maxima=np.delete(maxima, 0, 1)
    return(maxima)

def repeated(res,sizes):
    distance=np.zeros([len(res[0]),len(res[0])])
    for i in range(0,len(res[0])):         
            for j in range(0,len(res[0])):
                for k in range(0,4):
                       distance[i][j]+=np.copy((np.abs(res[k][i]-res[k][j])%sizes[k])*(np.abs(res[k][i]-res[k][j])%sizes[k]))
                distance[i][j]=np.copy(np.sqrt(distance[i][j]))
    repeated=0
    for i in range(0,len(res[0])):         
        for j in range(0,len(res[0])):
            if (distance[i][j]<=2 and i<j):
                print(i,j)
                repeated+=1    
    print("repeated maxima = ",repeated)
    return()                                                                                                                                                                mode_plot_2d.py                                                                                     0000600 0001750 0001750 00000001773 14607720207 012657  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys
import Maxima
Sum_over=(1,2)
Mode = str(sys.argv[1])
Bin = int(sys.argv[2])
Chirality= int(sys.argv[3])
colors = 8
spin_length=4

if Bin:
    #sizes=[4,4,4,32]
    sizes=[40,6,6,40]
    density_mode,sizes=Read.bin_mode(Mode,sizes,colors,spin_length)
    #Mode=Mode.replace("./SusyMode_bin_","")
    #Mode=Mode.split("-")
    #file_name="SusyMode_"+Mode[1]+"_"+Mode[0]
    file_name=Mode
else:
    zero_mode,density_mode,sizes=Read.ascii_mode(Mode)  
    file_name=Mode

density_2d=Chirality*density_mode.sum(axis=(Sum_over))
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# Make data.
X = np.arange(0,sizes[0])
Y = np.arange(0,sizes[3])
X, Y = np.meshgrid(X, Y)

#ax.set_zlim3d(-0.05,0)
ax.set_box_aspect((1,1,1)) 

# Plot the surface.
surf = ax.plot_surface(X, Y, density_2d, cmap='viridis')

plt.savefig(file_name+".png",dpi=150, bbox_inches='tight')
maxima=Maxima.simple(Chirality*density_mode,sizes)
print(maxima)
 
     mode_plot_all_ov.py                                                                                 0000600 0001750 0001750 00000004767 14607720206 013633  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys

Sum_over=(0,1,2)
max_conf = int(sys.argv[1])
All_summed = int(sys.argv[2])
Chirality= int(sys.argv[3])
colors = 3
spin_length=4
max_modes=12

if All_summed:
    sizes=[4,4,4,32]
    t=np.arange(sizes[3])
    for conf in range(10,max_conf,10):
        susy_mode_1=np.zeros(sizes[3])
        file="./OverlapMode"+str(conf)+"_"
        i=0
        susy_mode=np.zeros((spin_length, colors, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
        while i<max_modes:
            mode,density,sizes=Read.ascii_mode(file+str(i))
            susy_mode+=mode/np.sqrt(2)
            i+=1
            
            mode,density,sizes=Read.ascii_mode(file+str(i))
            susy_mode+=mode/np.sqrt(2)
            i+=1
            
        density_susy=Read.mode_to_density(susy_mode,colors,spin_length,sizes)
        plt.title("Configuration:"+str(conf))
        plt.plot(t,Chirality*density_susy, label="susy_mode")
        Topology="../../gf/profile4dt2c"+str(conf)+"to.dat"
        density_top,sizes=Read.topology_1d(Topology)
        plt.plot(t,density_top, label="top charge")
        plt.legend(loc="upper right")
        plt.title("Configuration:"+str(conf))
        plt.savefig("./Susymode_"+str(conf)+".png", dpi=150, bbox_inches='tight')
        plt.close()
else:
    sizes=[4,4,4,32]
    t=np.arange(sizes[3])
    for conf in range(10,max_conf,10):
        susy_mode_1=np.zeros(sizes[3])
        file="./OverlapMode"+str(conf)+"_"
        i=0
        while i<max_modes:
            susy_mode=np.zeros((spin_length, colors, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
            mode,density,sizes=Read.ascii_mode(file+str(i))
            susy_mode+=mode/np.sqrt(2)
            i+=1
            
            mode,density,sizes=Read.ascii_mode(file+str(i))
            susy_mode+=mode/np.sqrt(2)
            
            
            density_susy=Read.mode_to_density(susy_mode,colors,spin_length,sizes)
            plt.title("Configuration:"+str(conf))
            plt.plot(t,Chirality*density_susy, label="susy_mode")
            Topology="../../gf/profile4dt2c"+str(conf)+"to.dat"
            density_top,sizes=Read.topology_1d(Topology)
            plt.plot(t,density_top, label="top charge")
            plt.legend(loc="upper right")
            plt.title("Configuration:"+str(conf))
            plt.savefig("./Overlapmode_"+str(conf)+"_"+str(i)+".png", dpi=150, bbox_inches='tight')
            plt.close()
            
            i+=1         mode_plot_all.py                                                                                    0000600 0001750 0001750 00000002061 14607720165 013114  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys

Sum_over=(0,1,2)
max_conf = int(sys.argv[1])
Bin = int(sys.argv[2])
Chirality= int(sys.argv[3])
colors = 3
spin_length=4
max_modes=12

if Bin:
    sizes=[8,8,8,64]
    t=np.arange(sizes[3])
    for conf in range(120,max_conf,20):
        susy_mode_1=np.zeros(sizes[3])
        for i in range(0,max_modes):
            file="./SusyMode_bin_"+str(i)+"-"+str(conf)
            density,sizes=Read.bin_mode_1d(file,sizes,colors,spin_length)
            susy_mode_1+=Chirality*density/2
        plt.title("Configuration:"+str(conf))
        plt.plot(t,susy_mode_1, label="susy_mode")
        Topology="../../gf/profile4dt2c"+str(conf)+"to.dat"
        density_top,sizes=Read.topology_1d(Topology)
        plt.plot(t,density_top, label="top charge")
        plt.legend(loc="upper right")
        plt.title("Configuration:"+str(conf))
        plt.savefig("./Susymode_"+str(conf)+".png", dpi=150, bbox_inches='tight')
        plt.close()
else:
    density_mode,sizes=Read.ascii_mode(Mode)  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                               mode_plot_ov.py                                                                                     0000600 0001750 0001750 00000001167 14607720211 012766  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys

Sum_over=(0,1,2)
Mode_name=str(sys.argv[1])
Mode_number=str(sys.argv[2])
Chirality= 1
colors = 3
spin_length=4
sizes=[4,4,4,32]

susy_mode=np.zeros((spin_length, colors, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
for i in range(0,int(Mode_number)):
    mode,density,sizes=Read.ascii_mode(Mode_name+str(i))
    susy_mode+=mode/np.sqrt(2)

density_susy=Read.mode_to_density(susy_mode,colors,spin_length,sizes)

t=np.arange(0,sizes[3])
plt.plot(t,density_susy/2)

plt.savefig(Mode_name+"density.png",dpi=150, bbox_inches='tight')
                                                                                                                                                                                                                                                                                                                                                                                                         mode_plot.py                                                                                        0000644 0001750 0001750 00000001324 14607720206 012271  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys

Sum_over=(0,1,2)
Mode = str(sys.argv[1])
Bin = int(sys.argv[2])
Chirality= int(sys.argv[3])
colors = 3
spin_length=4

if Bin:
    sizes=[4,4,4,32]
    #sizes=[40,6,6,40]
    density_mode,sizes=Read.bin_mode(Mode,sizes,colors,spin_length)
    #Mode=Mode.replace("./SusyMode_bin_","")
    #Mode=Mode.split("-")
    #file_name="SusyMode_"+Mode[1]+"_"+Mode[0]
    file_name=Mode
else:
    zero_mode,density_mode,sizes=Read.ascii_mode(Mode)  
    file_name=Mode

density_mode=Chirality*density_mode.sum(axis=Sum_over)

t=np.linspace(0,sizes[3],sizes[3])
plt.plot(t,density_mode)
plt.savefig(file_name+".png",dpi=150, bbox_inches='tight')
                                                                                                                                                                                                                                                                                                            mode_plot_susy.py                                                                                   0000600 0001750 0001750 00000002057 14607720201 013343  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import Read
import analyzer
import sys

Sum_over=(0,1,2)
conf = str(sys.argv[1])
colors = 2
spin_length=4
max_modes_1=int(sys.argv[2])
max_modes_2=int(sys.argv[3])

sizes=[8,8,8,8]
t=np.arange(sizes[3])
susy_mode_1=np.zeros(sizes[3])
for i in range(0,max_modes_1):
    file="./sector_1/SusyMode_bin_"+str(i)+"-"+str(conf)
    density,sizes=Read.bin_mode_1d(file,sizes,colors,spin_length)
    susy_mode_1+=density/2

susy_mode_0=np.zeros(sizes[3])
for i in range(0,max_modes_2):
    file="./sector_0/SusyMode_bin_"+str(i)+"-"+str(conf)
    density,sizes=Read.bin_mode_1d(file,sizes,colors,spin_length)
    susy_mode_0+=-1*density/2

plt.title("Configuration:"+str(conf))
plt.plot(t,susy_mode_1+susy_mode_0, label="susy_mode")

Topology="../gf/profile4dt2c"+str(conf)+"to.dat"
#density_top,sizes=Read.topology_1d(Topology)
#plt.plot(t,density_top, label="top charge")
plt.legend(loc="upper right")
plt.title("Configuration:"+str(conf))
plt.savefig("./Susymode_"+str(conf)+".png", dpi=150, bbox_inches='tight')
plt.close()
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Plot_generator.py                                                                                   0000644 0001750 0001750 00000024075 14607720211 013277  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import analyzer
import Compare
import Read
import Maxima
import re
import os
import sys
import pickle
plt.rcParams.update({"text.usetex": True, "font.size": 16})



def susy_plot(folder_in,folder_out,sizes,colors,spin_length,max_modes,conf_read,susy_read_s0, susy_read_s1,pattern, cut, Load=False,Plot=True,Polyakov=False):
    
    dictionary_s1=analyzer.Real_eigenvalue(folder_in+"./sector_1/Measure.seq",pattern)
    dictionary_s0=analyzer.Real_eigenvalue(folder_in+"./sector_0/Measure.seq",pattern)

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
            density_susy=np.loadtxt(folder+measure+"susy_mode_"+str(conf)+"c_"+cut_+"cut.txt")
        else:
            if "OverlapFilterModeC" in pattern:
                density_susy=Compare.Construct_susy_overlap_density(folder_in,susy_read_s0[str(conf)],susy_read_s1[str(conf)],conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes)
            else:
                density_susy=Compare.Construct_susy_density(folder_in,susy_read_s0[str(conf)],susy_read_s1[str(conf)],conf,sizes,colors,spin_length,dictionary_s1,dictionary_s0,max_modes)

            np.savetxt(folder_in+"susy_mode_"+str(conf)+"c.txt", density_susy)
            

        #Plot the three densities
        plt.plot(density_top_1, label='top. density t=0.5')
        plt.plot(density_top_2, label='top. density t=2')
        plt.plot(density_top_3, label='top. density t=4')
        plt.plot(density_susy, label="AFM")
        plt.legend(loc="lower left", ncol=2)
        plt.savefig(folder_out+"./susy_mode_"+str(conf)+"c_"+cut+"cut.png",dpi=150, bbox_inches='tight')
        plt.close()      
        np.savetxt(folder_out+"./susy_mode_"+str(conf)+"c_"+cut+"cut.txt",density_susy)
    print(folder_out)
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
                plt.xlabel(r'Configuration')
                plt.ylabel(r'$$ \mbox{\huge $ \Xi$}$$')
                #plt.xlabel('Configuration')
                #plt.ylabel('Xi')
                plt.ylim(0.0,1.1)
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
        string=measure.split("_")
        time=re.sub("t","",string[2])
        time=re.sub("p",".", time)
        time=re.sub("/","", time)
        if len(string)>3:
            kappa=re.sub("/","",string[3])
            kappa="0."+re.sub("k","",kappa)
            label_afm=r'$\tau$'+"="+time +r', $\kappa$'+"="+kappa
            if string[3]=="_ov":
                label_afm=r'$\tau$'+"="+time + "_ov"
        else:
            label_afm=r'$\tau$'+"="+time    
        data=np.loadtxt(folder_in+measure+observable+".txt")
        susy_max=np.loadtxt(folder_in+measure+"end_spectrum.txt")
        #if method=="cut":
        susy_min=np.loadtxt(folder_in+measure+"start_spectrum.txt")

        #Read the error
        error=[]
        for t in range(0,len(data[0])):
            with open(folder_in+measure+observable+"_error_"+str(t)+".txt", 'r') as f:
                error.append(float(f.readline()))
            f.close()

        color = next(ax._get_lines.prop_cycler)['color']
        plt.errorbar(data[1,0:(int(susy_max[1])+1)],data[0,0:(int(susy_max[1])+1)], yerr=error[0:(int(susy_max[1])+1)], label=label_afm, color=color)
        plt.fill_between(data[1,0:(int(susy_max[1])+1)], data[0,0:(int(susy_max[1])+1)]-error[0:(int(susy_max[1])+1)], data[0,0:(int(susy_max[1])+1)]+error[0:(int(susy_max[1])+1)], alpha=0.1, color=color)
        #plt.scatter(susy_min[0],data[0,(int(susy_min[1])-1)], marker="v", color=color)

    plt.xlabel(r'$$ \mbox{\huge $\lambda$}_{cut} $$')
    plt.ylabel(r'$$ \mbox{\huge $ \Xi$}$$')
    #plt.xlabel('lambda')
    #plt.ylabel('Xi')
    plt.legend(loc="upper right", ncol=1)
    plt.ylim([0,1.1])
    plt.xlim([0.003,0.0175])
    plt.xticks([0.004, 0.008,0.012,0.016])

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

def histogram(folder,file_name,conf_read,max_modes, susy_read_s0, susy_read_s1):
    
    GM_hist={}
    print(folder+file_name)
    with open(folder+file_name+".txt") as f:
        for line in f:
            conf=int(line.split(",")[0].replace("[","").replace("]",""))
            i=int(line.split(",")[1].replace("[","").replace("]",""))
            j=int(line.split(",")[2].replace("[","").replace("]",""))
            GM_hist[str(conf)+","+str(i)+","+str(j)]=float(line.split(",")[3].replace("[","").replace("]",""))

    max_xi_dic={}
    for conf in conf_read:
        for i in range(0,susy_read_s1[str(conf)]):
            max_xi_dic[str(conf)+","+str(i)]=0
            for j in range(0,susy_read_s0[str(conf)]):
                if max_xi_dic[str(conf)+","+str(i)]<GM_hist[str(conf)+","+str(i)+","+str(j)]:
                    max_xi_dic[str(conf)+","+str(i)]=GM_hist[str(conf)+","+str(i)+","+str(j)]
    print(max_xi_dic)
    max_xi_list=[]
    for key in max_xi_dic:
        max_xi_list.append(max_xi_dic[key])

    max_xi_np=np.array(max_xi_list)
    
    x_min=0
    x_max=10
    #bins=10
    #x_coord=np.range(x_min,x_max,bins)
    plt.hist(max_xi_np, bins=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], density=False)
    plt.xlabel(r'$\Xi$')
    plt.ylabel(r'Frequency')
    print(folder, file_name)
    plt.savefig(folder+file_name+".pdf",dpi=150, bbox_inches='tight')
    plt.close()
    return()
    
                                                                                                                                                                                                                                                                                                                                                                                                                                                                   plot_susy.py                                                                                        0000600 0001750 0001750 00000001607 14607720210 012337  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import Read
import Maxima
import re
import os
import sys
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import importlib
importlib.reload(Read)
importlib.reload(Maxima)

sizes=[4,4,4,32]
colors=3
spin=4

if sys.argv[2]:
    sizes=[16,8,8,8]
    density,sizes=Read.bin_mode(sys.argv[1],sizes,colors,spin)
    density_2d=density.sum(axis=(1,2))
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    X=np.arange(0,sizes[3])
    Y=np.arange(0,sizes[0])
    
    X, Y = np.meshgrid(X, Y)
    
    ax.plot_surface(X,Y,density_2d, cmap="viridis")
    plt.savefig("./"+sys.argv[1] + ".png",dpi=150, bbox_inches='tight')
    
else:
    density,sizes=Read.bin_mode_1d(sys.argv[1],sizes,colors,spin)
    t=np.linspace(0,sizes[3],32)
    plt.plot(t,density)
    plt.savefig("./"+sys.argv[1] + ".png",dpi=150, bbox_inches='tight')                                                                                                                         plot_top.py                                                                                         0000600 0001750 0001750 00000000637 14607720164 012150  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import matplotlib.pyplot as plt
import numpy as np
import Read
import Maxima
import re
import os
import sys
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import importlib
importlib.reload(Read)
importlib.reload(Maxima)

density,sizes=Read.topology_1d(sys.argv[1])
t=np.linspace(0,sizes[3],sizes[3])
plt.plot(t,density)
plt.savefig("./"+sys.argv[1] + ".png",dpi=150, bbox_inches='tight')
                                                                                                 Read.py                                                                                             0000644 0001750 0001750 00000027212 14607720205 011165  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   import numpy as np
import struct
import array
import os 

def ascii_mode(directory):
    sizes=[[]]*4
    count=0
    i=0
    j=0
    k=0
    l=0
    
    with open(directory,"r") as file:
        file.readline()
        file.readline()
        file.readline()
        info=file.readline()
        info=info.split(" ")

        for m in range (0,4):
            sizes[m]=int(info[m+2])
        range_color=int(info[0])
        range_spin=int(info[1])
        elements=sizes[0]*sizes[1]*sizes[2]*sizes[3]   
        zeromode=np.zeros((range_color, range_spin, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
        for l in range(0, sizes[3]):
            for k in range(0, sizes[2]):
                for j in range(0, sizes[1]):
                    for i in range(0,sizes[0]):
                        for spin in range(0,range_spin):
                            for color in range(0,range_color):
                                count+=1
                                num=file.readline()
                                num=num.replace("(","")
                                num=num.replace(")","")
                                num=num.replace("\n","")
                                split=num.split(",")
                                zeromode[color,spin,i,j,k,l]=complex(float(split[0]),float(split[1]))
    density=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])

    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for spin in range(0,range_spin):
                        for color in range(0,range_color):                    
                            density[i,j,k,l]+=np.copy(zeromode[color,spin,i,j,k,l].imag*zeromode[color,spin,i,j,k,l].imag + 
                                                      zeromode[color,spin,i,j,k,l].real*zeromode[color,spin,i,j,k,l].real)
    norm=density.sum()                       
    return(zeromode,density,sizes)

def mode_to_density(zeromode,range_color,range_spin,sizes):
    density=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for spin in range(0,range_spin):
                        for color in range(0,range_color):
                            density[i,j,k,l]+=np.copy(zeromode[color,spin, i,j,k,l].imag*zeromode[color,spin, i,j,k,l].imag + 
                                                      zeromode[color,spin, i,j,k,l].real*zeromode[color,spin, i,j,k,l].real)          
    #density_1d=density.sum(axis=(0,1,2))
    return(density)


def mode_real_density(zeromode,range_color,range_spin,sizes,chirality):
    density=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])
    if chirality==1:
        spin=0
    if chirality==0:
        spin=2        
    for l in range(0, sizes[3]):
        for k in range(0, sizes[2]):
            for j in range(0, sizes[1]):
                for i in range(0,sizes[0]):
                    for color in range(0,range_color):
                            zeromode[color,spin,i,j,k,l]=np.real(zeromode[color,spin,i,j,k,l])
    
    density=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])
    
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for spin in range(0,range_spin):
                        for color in range(0,range_color):
                            density[i,j,k,l]+=np.copy(zeromode[spin,color,i,j,k,l].imag*zeromode[spin,color,i,j,k,l].imag + 
                                                      zeromode[spin,color,i,j,k,l].real*zeromode[spin,color,i,j,k,l].real)          
    density_1d=density.sum(axis=(0,1,2))
    return(density_1d)

def bin_mode(file,sizes,range_color,range_spin):

    mode=np.memmap(file, dtype=np.double).byteswap()
    zeromode=np.zeros((range_color, range_spin, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
    density=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])
    index=0
    for l in range(0, sizes[3]):
        for k in range(0, sizes[2]):
            for j in range(0, sizes[1]):
                 for i in range(0,sizes[0]):
                    for spin in range(0,range_spin):
                        for color in range(0,range_color):
                            zeromode[color,spin,i,j,k,l]=complex(mode[index],mode[index+1])
                            density[i,j,k,l]+=mode[index]*mode[index]+mode[index+1]*mode[index+1]
                            index+=2    
    return zeromode,density,sizes

def table_fund(file,sizes,range_color,range_spin):
    size=sizes[0]*sizes[1]*sizes[2]*sizes[3]*range_color*range_spin
    index=size*[[]]
    with open(file,"r") as file:
        for line in file:
            #line=file.readline()
            sline=line.split(" ")
            index[int(sline[1])]=[int(sline[0]),int(sline[2]),int(sline[3]),int(sline[4]),int(sline[5]),int(sline[6]),int(sline[7])]
    return index

def CoordToSite(table,sizes):
    size=sizes[0]*sizes[1]*sizes[2]*sizes[3]
    coord=size*[[]]
    for i in range(0,len(table)):
        coord[table[i][0]]=[table[i][3],table[i][4],table[i][5],table[i][6]]
        
    return coord
    

def bin_space(file,sizes,range_color,range_spin,max_modes):

    eigenspace=np.memmap(file, dtype=np.csingle)
    modes={}
    eigenvalues=np.zeros([max_modes],dtype=np.csingle)
    
    volume=range_color*range_spin*sizes[0]*sizes[1]*sizes[2]*sizes[3]
    
    for m in range(0,max_modes):
        modes[m]=np.zeros([volume],dtype=np.csingle)
        for index in range(0,volume):
            modes[m][index]=eigenspace[index+m*(volume+1)] #+volume for next
            
        eigenvalues[m]=eigenspace[volume+(m*(volume+1))]

    #os.remove(temp_file)
    return modes,eigenvalues,eigenspace


def scalar_product(mode1, mode2,sizes,range_color,range_spin):
    scalar=0
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for color in range(0,range_color):
                        for spin in  range(0,range_spin):
                            scalar+=mode1[color,spin,i,j,k,l].conjugate()*mode2[color,spin,i,j,k,l]
    return np.sqrt(scalar),sizes

def scalar_g5_product(mode1, mode2,sizes,range_color,range_spin):
    scalar=0
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for color in range(0,range_color):
                        for spin in range(0,range_spin):
                            if spin==0 or spin==1:
                                g5=1
                            if spin==2 or spin==3:
                                g5=-1
                            scalar+=g5*mode1[color,spin, i,j,k,l].conjugate()*mode2[color,spin, i,j,k,l] 
    return np.sqrt(scalar),sizes

def real_cond(mode1,sizes,range_color,range_spin,chirality):
    scalar=0
    if chirality==0:
        spin=2
    if chirality==1:
        spin=0        
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for color in range(0,range_color):
                            scalar+=np.real(mode1[color,spin, i,j,k,l])*np.real(mode1[color,spin, i,j,k,l])
    return np.sqrt(scalar)

def imag_cond(mode1,sizes,range_color,range_spin,chirality):
    scalar=0
    if chirality==0:
        spin=2
    if chirality==1:
        spin=0        
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for color in range(0,range_color):
                            scalar+=np.imag(mode1[color,spin, i,j,k,l])*np.imag(mode1[color,spin, i,j,k,l])
    return np.sqrt(scalar)

def scalar_product(mode1, mode2,sizes,range_color,range_spin):
    scalar=0
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for color in range(0,range_color):
                        for spin in  range(0,range_spin):
                            scalar+=mode1[color,spin, i,j,k,l].conjugate()*mode2[color,spin, i,j,k,l] + mode1[color,spin, i,j,k,l].conjugate()*mode2[color,spin, i,j,k,l]
    return np.sqrt(scalar),sizes

def scalar_g5_product(mode1, mode2,sizes,range_color,range_spin):
    scalar=0
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for color in range(0,range_color):
                        for spin in range(0,range_spin):
                            if spin==0 or spin==1:
                                g5=1
                            if spin==2 or spin==3:
                                g5=-1
                            scalar+=g5*mode1[color,spin, i,j,k,l].conjugate()*mode2[color,spin, i,j,k,l] + g5*mode1[color,spin, i,j,k,l].conjugate()*mode2[color,spin, i,j,k,l]
    return np.sqrt(scalar),sizes

def ascii_mode_1d(namefile):
    density,sizes=ascii_mode(namefile)
    density_1d=density.sum(axis=(0,1,2))
    return(density_1d,sizes)

def bin_mode_1d(namefile,sizes,colors,spin):
    zeromode,density,sizes=bin_mode(namefile,sizes,colors,spin)
    density_1d=density.sum(axis=(0,1,2))
    return(density_1d,sizes)

def topology(file):
    sizes=[[]]*4
    precision=0

    count=0
    i=0
    j=0
    k=0
    l=0
    norm=0
    with open(file,"rb") as file:
        for i in range (0,4):
            sizes[i]=int(file.readline())
        precision=int(file.readline())
        colors=int(file.readline())
        elements=sizes[0]*sizes[1]*sizes[2]*sizes[3]
        density=np.zeros((sizes[0],sizes[1],sizes[2],sizes[3]))
        for i in range(0, sizes[0]):
            for j in range(0, sizes[1]):
                for k in range(0, sizes[2]):
                    for l in range(0,sizes[3]):
                        count+=1
                        num=struct.unpack('d', file.read(precision))
                        density[i,j,k,l]=num[0]
                        norm+=num[0]
    if (count==elements):
        return (density,sizes)

    if (count!=elements):
        print("Four dimensional array wrongly readed")
        return (density,sizes)


def topology_1d(file):
    density,sizes=topology(file)
    density_1d=density.sum(axis=(0,1,2))
    return(density_1d,sizes)

def shift_t(density,sizes,shift):
    temp=density
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            for k in range(0,sizes[2]):
                for l in range(0,sizes[3]):
                    temp[i,j,k,l]=density[(i+shift)%sizes[0],j,k,l]
    return temp

def shift_x(density,sizes,shift):
    temp=density
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            for k in range(0,sizes[2]):
                for l in range(0,sizes[3]):
                    temp[i,j,k,l]=density[i,(j+shift)%sizes[1],k,l]
    return temp

def shift_y(density,sizes,shift):
    temp=density
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            for k in range(0,sizes[2]):
                for l in range(0,sizes[3]):
                    temp[i,j,k,l]=density[i,j,(k+shift)%sizes[2],l]
    return temp

def shift_z(density,sizes,shift):
    temp=density
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            for k in range(0,sizes[2]):
                for l in range(0,sizes[3]):
                    temp[i,j,k,l]=density[i,j,k,(l+shift)%sizes[3]]
    return temp

                                                                                                                                                                                                                                                                                                                                                                                      Tools.py                                                                                            0000600 0001750 0001750 00000000223 14607720210 011367  0                                                                                                    ustar   ivsol                           ivan                                                                                                                                                                                                                   def charge_inv(mu,nu):
    
    
    return(index_mu,index_nu
           
def g_5(alpha):
    if alpha>
    
    return(index_mu,index_nu)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      