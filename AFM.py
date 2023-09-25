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
import Fit_1d
import subprocess
from subprocess import call
import analyzer
import re

threshold_ov=float(input("Enter threshold for Overlap zero modes: "))
threshold_susy=float(input("Enter threshold for Susy zero modes: "))
#Reading eigenvalues chirality 0 and applying the threshold
Measure_file=None
try:
    Measure_file = open("./sector_0/Measure.seq", "r")
except IOError:
    print("File couldn't open")
else:         
    pattern_ov=":OverlapFilterModeC:"
    pattern_susy=":OverlapFilterModeR:"
    eigven_susy_s0=[[]]
    eigven_susy_s0.pop(0)
    eigven_ov_s0=[[]]
    eigven_ov_s0.pop(0)
    for line in Measure_file:
        for line in Measure_file:
            if re.search(pattern_ov,line):
                line_split=line.split(":")
                if (float(line_split[8])<threshold_ov):
                    eigven_ov_s0.append([line_split[1],line_split[4],line_split[8]])
                    
            if re.search(pattern_susy,line):
                line_split=line.split(":")
                if (float(line_split[8])<threshold_susy):
                    eigven_susy_s0.append([line_split[1],line_split[4],line_split[8]])
finally:
    if Measure_file: Measure_file.close()

#Reading eigenvalues chirality 1 and applying the threshold
Measure_file=None
try:
    Measure_file = open("./sector_1/Measure.seq", "r")
except IOError:
    print("File couldn't open")
else:         
    pattern_ov=":OverlapFilterModeC:"
    pattern_susy=":OverlapFilterModeR:"
    eigven_susy_s1=[[]]
    eigven_susy_s1.pop(0)
    eigven_ov_s1=[[]]
    eigven_ov_s1.pop(0)
    for line in Measure_file:
        for line in Measure_file:
            if re.search(pattern_ov,line):
                line_split=line.split(":")
                if (float(line_split[8])<threshold_ov):
                    eigven_ov_s1.append([line_split[1],line_split[4],line_split[8]])
                    
            if re.search(pattern_susy,line):
                line_split=line.split(":")
                if (float(line_split[8])<threshold_susy):
                    eigven_susy_s1.append([line_split[1],line_split[4],line_split[8]])
finally:
    if Measure_file: Measure_file.close()


#Reading Overlap eigenvectors and summing over
density_overlap={}
count_ov_s0={}
for element in eigven_ov_s0:
    density_temp,sizes=Read.susy_mode_1d("./sector_0/OverlapMode"+element[0]+element[1])
    if element[0] in density_overlap:
        density_overlap[element[0]]-=density_temp
        count_ov_s0[element[0]]+=1
    else:
        density_overlap.update({element[0]:-density_temp})
        count_ov_s0.update({element[0]:1})

count_ov_s1={}        
for element in eigven_ov_s1:
    density_temp,sizes=Read.susy_mode_1d("./sector_1/OverlapMode"+element[0]+element[1])
    if element[0] in density_overlap:
        density_overlap[element[0]]+=density_temp
    else:
        density_overlap.update({element[0]:density_temp})
    if element[0] in count_ov_s1:
        count_ov_s1[element[0]]+=1
    else:
        count_ov_s1.update({element[0]:1})

t=np.arange(0,sizes[3])
for key in density_overlap:
    Plotting.plot_list(t,density_overlap[key],"./OverlapMode_sum_"+key+".png","black")
    Plotting.plot_list(t,Read.density_1d("../gf/profile4dt2c"+key+"to.dat")[0],"./OverlapMode_sum_"+key+".png","red")
    plt.close()

print("Eigenvectors with 0 chirality")
print(count_ov_s0)

print("Eigenvectors with 1 chirality")
print(count_ov_s1)
    
    
#Reading Susy eigenvectors and summing over
density_susy={}
count_susy_s0={}
for element in eigven_susy_s0:
    density_temp,sizes=Read.susy_mode_1d("./sector_0/SusyMode"+element[0]+element[1])
    if element[0] in density_susy:
        density_susy[element[0]]-=density_temp
        count_susy_s0[element[0]]+=1
    else:
        density_susy.update({element[0]:-density_temp})
        count_susy_s0.update({element[0]:1})

count_susy_s1={}        
for element in eigven_susy_s1:
    density_temp,sizes=Read.susy_mode_1d("./sector_1/SusyMode"+element[0]+element[1])
    if element[0] in density_susy:
        density_susy[element[0]]+=density_temp
    else:
        density_susy.update({element[0]:density_temp})
    if element[0] in count_susy_s1:
        count_susy_s1[element[0]]+=1
    else:
        count_susy_s1.update({element[0]:1})

t=np.arange(0,sizes[3])
for key in density_susy:
    Plotting.plot_list(t,density_susy[key],"./SusyMode_sum_"+key+".png","black")
    Plotting.plot_list(t,Read.density_1d("../gf/profile4dt2c"+key+"to.dat")[0],"./SusyMode_sum_"+key+".png","red")
    plt.close()

print("Susy vectors with 0 chirality")
print(count_susy_s0)

print("Susy vectors with 1 chirality")
print(count_susy_s1)
    