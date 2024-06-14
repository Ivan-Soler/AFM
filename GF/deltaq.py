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

deltaq_dic={}
odd_dic={}

deltaq_dici={}
odd_dici={}

deltaq_dicd={}
odd_dicd={}
file_list=os.listdir("./")

for file in  file_list:
  if ".pkl" in file:
    print(file)
    scale=re.search(r'_s(.*?(?=.pkl))',file).group().replace('_s','')
    deltaq_dic_temp=pickle.load(open(file, "rb" ))
    odd_dic_temp=pickle.load(open(file, "rb" ))
    
    for key in deltaq_dic_temp:
        if key not in deltaq_dic:
            deltaq_dic[key]=[]
        if key not in odd_dic:
            odd_dic[key]=[]
        deltaq_dic[key].append([float(scale),deltaq_dic_temp[key]['top-frac'][0],deltaq_dic_temp[key]['top-frac'][1],deltaq_dic_temp[key]['top-frac'][2]])                  
        odd_dic[key].append([float(scale),odd_dic_temp[key]["odds"]])

    if False:
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
            
for key in deltaq_dic:
  deltaq_dic[key]=np.array((deltaq_dic[key]))
  deltaq_dic[key]=np.array((deltaq_dic[key][deltaq_dic[key][:,0].argsort()]))
  odd_dic[key]=np.array((odd_dic[key]))
  odd_dic[key]=np.array((odd_dic[key][odd_dic[key][:,0].argsort()]))


for key in deltaq_dic:
      plt.plot(deltaq_dic[key][:,0]/10,deltaq_dic[key][:,1],label="fractionals")
      plt.plot(deltaq_dic[key][:,0]/10,deltaq_dic[key][:,2],label="+ instantons")
      plt.plot(deltaq_dic[key][:,0]/10,deltaq_dic[key][:,3],label="all fractionals")
      #plt.plot(n_list,deltaq_dicd[key],label=" Q=2 inst")
      plt.legend(loc="upper right")
      #plt.ylim(0,5)
      plt.ylabel("Q_top - Q_inst")
      plt.xlabel("delta norm")
      plt.savefig("./deltaq/"+str(key)+".png")
      plt.close()
  
      plt.plot(odd_dic[key][:,0],odd_dic[key][:,1],label="fractionals")
      plt.xlabel("delta norm")
      plt.ylabel("% of odd configurations")
      plt.savefig("./odd/"+str(key)+".png")
      plt.close()

#print(odd_dic)