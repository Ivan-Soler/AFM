import matplotlib.pyplot as plt
import numpy as np
import Read
import math
import re


plt.rcParams.update({'font.size': 12})        

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

def End_spectrum(folder):
    #Overlap
    end_overlap=100
    highest_overlap={}
    with open(folder+"sector_0/Measure.seq") as file:
        for line in file:
            match=re.search(":OverlapFilterModeC:",line)
            if match:
                string=line.split(":")
                if not string[1] in highest_overlap:
                    highest_overlap[string[1]]=abs(float(string[8]))
                else:
                    if (abs(float(string[8]))>highest_overlap[string[1]]):
                        highest_overlap[string[1]]=abs(float(string[8]))
                    
    #Overlap
    with open(folder+"sector_1/Measure.seq") as file:
        for line in file:
            match=re.search(":OverlapFilterModeC:",line)
            if match:
                string=line.split(":")
                if not string[1] in highest_overlap:
                    highest_overlap[string[1]]=abs(float(string[8]))
                else:
                    if (abs(float(string[8]))>highest_overlap[string[1]]):
                        highest_overlap[string[1]]=abs(float(string[8]))
                    
    for key in highest_overlap:
        if highest_overlap[key]<end_overlap:
            end_overlap=highest_overlap[key]
            
    #susy
    end_susy=100
    highest_susy={}
    with open(folder+"sector_0/Measure.seq") as file:
        for line in file:
            match=re.search(":OverlapFilterModeR:",line)
            if match:
                string=line.split(":")
                if not string[1] in highest_susy:
                    highest_susy[string[1]]=abs(float(string[8]))
                else:
                    if (abs(float(string[8]))>highest_susy[string[1]]):
                        highest_susy[string[1]]=abs(float(string[8]))
                    
    #susy
    with open(folder+"sector_1/Measure.seq") as file:
        for line in file:
            match=re.search(":OverlapFilterModeR:",line)
            if match:
                string=line.split(":")
                if not string[1] in highest_susy:
                    highest_susy[string[1]]=abs(float(string[8]))
                else:
                    if (abs(float(string[8]))>highest_susy[string[1]]):
                        highest_susy[string[1]]=abs(float(string[8]))
                    
    for key in highest_susy:
        if highest_susy[key]<end_susy:
            end_susy=highest_susy[key]
  
    return(end_overlap,end_susy)
 

def Count_index(folder,measurement,threshold,configurations):
    count={}
    for conf in configurations:
        count[str(conf)]=0

    with open(folder) as file:
        for line in file:
            match=re.search(measurement,line)
            if match:
                string=line.split(":")
                for conf in configurations:
                    if (abs(float(string[8]))<threshold) and string[1]==str(conf):
                        count[string[1]]+=1
                        continue
    return(count)

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
    with open(folder+"../gf/Measure.seq") as file:
        for line in file:
            match=re.search(":TopologicalCharge:t:4",line)
            if match:
                string=line.split(":")
                for conf in configurations: #check that the configuration is on the list of configurations to be read
                    if string[1]==str(conf):
                        count_gauge[string[1]]=round(float(string[7]))
                        break
    return(count_gauge)

def Topology_dic(folder,threshold,configurations):

    #Real
    count_s0_susy=Count_index(folder+"sector_0/Measure.seq", ":OverlapFilterModeR:",threshold,configurations)
    count_s1_susy=Count_index(folder+"sector_1/Measure.seq", ":OverlapFilterModeR:",threshold,configurations)

    #Create the dictionary with the configurations that were readed. Probably we have to change this such that one 
    #uses the GF and puts a cut of what is an instanton. For now we use the susy with a large threashold
    conf_read={}
    for key in count_s0_susy:
        conf_read[key]=True
    for key in count_s1_susy: 
        conf_read[key]=True

    #GF
    count_gauge=Count_index_gf(folder,conf_read)

    #Now we have to redo the count for the susy, this we can change in the future when we use GF to set the conf_read
    count_s0_susy={}
    count_s1_susy={}
    count_s0_susy=Count_index(folder+"sector_0/Measure.seq", ":OverlapFilterModeR:",threshold,conf_read)
    count_s1_susy=Count_index(folder+"sector_1/Measure.seq", ":OverlapFilterModeR:",threshold,conf_read)

    #Overlap
    count_s0_ov=Count_index(folder+"sector_0/Measure.seq", ":OverlapFilterModeC:",threshold,conf_read)
    count_s1_ov=Count_index(folder+"sector_1/Measure.seq", ":OverlapFilterModeC:",threshold,conf_read)

    # Need to fill all the read configurations

    #for key in count_gauge: print((count_s1_ov[key] + count_s0_ov[key]))
    ov_top_dif={key: (count_s1_ov[key] - count_s0_ov[key])/4. - count_gauge[key] for key in conf_read}
    susy_top_dif={key: (count_s1_susy[key] - count_s0_susy[key])/2. - count_gauge[key] for key in conf_read}

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
