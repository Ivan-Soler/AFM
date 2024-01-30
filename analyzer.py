import matplotlib.pyplot as plt
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

    for key in highest_susy_s0:
        if highest_susy_s0[key]<end_susy:
            end_susy=highest_susy_s0[key]
  
    return(end_susy,end_overlap)

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
    start_spectrum=True

    for conf in conf_read:
        if susy_read_s0[str(conf)]==0 and susy_read_s1[str(conf)]==0:
            start_spectrum=False
            break
            #ev1=float(dictionary_s1[str(conf)][0])
            #ev2=float(dictionary_s0[str(conf)][0])
            #if ev1<ev2:
            #    susy_read_s1[str(conf)]=1
            #else:
            #    susy_read_s0[str(conf)]=1
    
    return(susy_read_s0,susy_read_s1,start_spectrum)

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

