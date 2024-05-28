import os
import matplotlib.pyplot as plt
import sys
import re
import numpy as np
import tools
import pandas as pd
from collections import OrderedDict
import pickle
import copy

n_list=np.arange(10,80,2)

deltaq_dic={}
odd_dic={}

deltaq_dici={}
odd_dici={}

deltaq_dicd={}
odd_dicd={}

for n in n_list:
    deltaq_dic_temp=pickle.load(open('deltaq_frac/deltaq_'+str(n)+'.pkl', "rb" ))
    odd_dic_temp=pickle.load(open('deltaq_frac/odd.pkl_'+str(n)+'.pkl', "rb" ))
    
    for key in deltaq_dic_temp:
        if key not in deltaq_dic:
            deltaq_dic[key]=[]
        if key not in odd_dic:
            odd_dic[key]=[]
        if deltaq_dic_temp[key]:
            deltaq_dic[key].append(deltaq_dic_temp[key][0])
        if odd_dic_temp[key]:                   
            odd_dic[key].append(odd_dic_temp[key][0])

for n in n_list:
    deltaq_dic_temp=pickle.load(open('deltaq_inst/deltaq_'+str(n)+'.pkl', "rb" ))
    for key in deltaq_dic_temp:
        if key not in deltaq_dici:
            deltaq_dici[key]=[]
        if deltaq_dic_temp[key]:
            deltaq_dici[key].append(deltaq_dic_temp[key][0])

for n in n_list:
    deltaq_dic_temp=pickle.load(open('deltaq_dinst/deltaq_'+str(n)+'.pkl', "rb" ))
    for key in deltaq_dic_temp:
        if key not in deltaq_dicd:
            deltaq_dicd[key]=[]
        if deltaq_dic_temp[key]:
            deltaq_dicd[key].append(deltaq_dic_temp[key][0])
            
            
#print(n_list)
print(odd_dic[key])
for key in deltaq_dic:
    if deltaq_dic[key]:
        plt.plot(n_list,deltaq_dic[key],label="fractionals")
        plt.plot(n_list,deltaq_dici[key],label=" Q=1 inst")
        plt.plot(n_list,deltaq_dicd[key],label=" Q=2 inst")
        plt.legend(loc="upper right")
        plt.savefig("./deltaq/"+str(key)+".png")
        plt.close()
    #if odd_dic[key]:         
    #    plt.plot(n_list,odd_dic[key])
    #    plt.savefig("./odd/"+str(key)+".png")
    #    plt.close()

#print(deltaq_dic)
#print(odd_dic)