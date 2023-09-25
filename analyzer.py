import matplotlib.pyplot as plt
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

    return density, param

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


def RPO(densityA,densityB,thresholdA):
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
 

def Count_index(folder,measurement,threshold):
    count={}
    with open(folder) as file:
        for line in file:
            match=re.search(measurement,line)
            if match:
                string=line.split(":")
                if (abs(float(string[8]))<threshold):
                    if string[1] in count:
                        count[string[1]]+=1
                    else:
                        count[string[1]]=1
    return(count)

def Count_index_impr(folder,measurement,threshold,modes_used):
    count={}
    with open(folder) as file:
        for line in file:
            match=re.search(measurement,line)
            if match:
                string=line.split(":")
                if (abs(float(string[8]))<threshold) and modes_used[string[7]]:
                    if string[1] in count:
                        count[string[1]]+=1
                    else:
                        count[string[1]]=1
    return(count)

def Topology_dic(folder,threshold):

    #Overlap
    count_s0_ov=Count_index(folder+"sector_0/Measure.seq", ":OverlapFilterModeC:",threshold)
    count_s1_ov=Count_index(folder+"sector_1/Measure.seq", ":OverlapFilterModeC:",threshold)
    count_s0_susy=Count_index(folder+"sector_0/Measure.seq", ":OverlapFilterModeR:",threshold,)
    count_s1_susy=Count_index(folder+"sector_1/Measure.seq", ":OverlapFilterModeR:",threshold)
    #GF
    top_gauge={}
    susy_top_tot={}
    with open(folder+"../gf/Measure.seq") as file:
        for line in file:
            match=re.search(":TopologicalCharge:t:4",line)
            if match:
                string=line.split(":")
                top_gauge[string[1]]=round((float(string[7])))
    for key in top_gauge:
        if key not in count_s0_ov: count_s0_ov[key]=0
        if key not in count_s1_ov: count_s1_ov[key]=0
        if key not in count_s0_susy: count_s0_susy[key]=0
        if key not in count_s1_susy: count_s1_susy[key]=0
    conf_tot=len(top_gauge)
    #for key in top_gauge: print((count_s1_ov[key] + count_s0_ov[key]))
    ov_top_dif={key: (count_s1_ov[key] - count_s0_ov[key])/4. - top_gauge[key] for key in top_gauge}
    susy_top_dif={key: (count_s1_susy[key] - count_s0_susy[key])/2. - top_gauge[key] for key in top_gauge}
    
    for key in top_gauge:
        if count_s0_susy[key]!=0 and count_s1_susy[key]!=0 :
            susy_top_tot[key]=True
        else :
            susy_top_tot[key]=False
    return(ov_top_dif,susy_top_dif,conf_tot,susy_top_tot)
    
def Topology_dic_impr(folder,threshold,modes_used_s0,modes_used_s1):

    #Overlap
    count_s0_ov=Count_index(folder+"sector_0/Measure.seq", ":OverlapFilterModeC:",threshold)
    count_s1_ov=Count_index(folder+"sector_1/Measure.seq", ":OverlapFilterModeC:",threshold)
    count_s0_susy=Count_index(folder+"sector_0/Measure.seq", ":OverlapFilterModeR:",threshold,modes_used_s0)
    count_s1_susy=Count_index(folder+"sector_1/Measure.seq", ":OverlapFilterModeR:",threshold,modes_used_s1)
    #GF
    top_gauge={}
    with open(folder+"../gf/Measure.seq") as file:
        for line in file:
            match=re.search(":TopologicalCharge:t:4",line)
            if match:
                string=line.split(":")
                top_gauge[string[1]]=round((float(string[7])))
    for key in top_gauge:
        if key not in count_s0_ov: count_s0_ov[key]=0
        if key not in count_s1_ov: count_s1_ov[key]=0
        if key not in count_s0_susy: count_s0_susy[key]=0
        if key not in count_s1_susy: count_s1_susy[key]=0
    conf_tot=len(top_gauge)
    #for key in top_gauge: print((count_s1_ov[key] + count_s0_ov[key]))
    ov_top_dif={key: (count_s1_ov[key] - count_s0_ov[key])/4. - top_gauge[key] for key in top_gauge}
    susy_top_dif={key: (count_s1_susy[key] - count_s0_susy[key])/2. - top_gauge[key] for key in top_gauge}
    
    return(ov_top_dif,susy_top_dif,conf_tot)