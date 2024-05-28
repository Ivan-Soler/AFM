import os
import matplotlib.pyplot as plt
import sys
import re
import numpy as np
import tools
import pandas as pd

def extract_feature(line, param):
    splitted=line.split(" ")
    feature=[]
    for i in range(0,len(splitted)):
        if not (i-param-5) % 7:
            feature.append(float(splitted[i]))
     
    return(feature)

def plot_histo(data,bins,xlabel,figure_name,xrange=())
    plt.hist(frac,bins=100, density=True)
    title="# points:" +str(len(data))
    plt.ylabel("frequency")
    plt.xlabel(xlabel)
    plt.title(title)
    plt.xlim(xrange)
    plt.savefig(figure_name)
    plt.close()

    return 
        
        
#nt=sys.argv[1]
t=sys.argv[1]
hist=int(sys.argv[2])
norm_cut=float(sys.argv[3])

lis_dir=os.listdir("./")
folders=[]
for direct in lis_dir:
    if "runns" in direct:
        folders.append(direct)
        
	
a={
	"2.4":0.12130,
	"2.5":0.08194,
	"2.6":0.05938,
	"2.7":0.04337,
	"2.8":0.03238,
	"2.9":0.02456
}

lt_frac={
   "4":0.02,
   "5":0.02,
   "6":0.02,
   "7":0.01,
   "8":0.01,
   "9":0.005,
   "10":0.005,
   "11":0.004,
   "12":0.004,
   "13":0.004,
   "14":0.004,
   "15":0.00,
   "16":0.00,
   "17":0.00,
   "18":0.00,
   "19":0.00,
   "20":0.00
}


lt_inst={
   "4":0.1,
   "5":0.08,
   "6":0.05,
   "7":0.05,
   "8":0.05,
   "9":0.05,
   "10":0.05,
   "11":0.5,
   "12":0.5,
   "13":0.5,
   "14":0.5,
   "15":1,
   "16":1,
   "17":1,
   "18":1,
   "19":1,
   "20":1,
}


norm_min={
   "4":0.5,
   "5":0.5,
   "6":0.5,
   "7":0.5,
   "8":0.5,
   "9":0.5,
   "10":0.5,
   "11":0.5,
   "12":0.5,
   "13":0.5,
   "14":0.5,
   "15":5,
   "16":5,
   "17":5,
   "18":5,
   "19":5,
   "20":5
}

norm_max={
   "4":2.25,
   "5":2.25,
   "6":2.25,
   "7":2.25,
   "8":2.25,
   "9":2.25,
   "10":2.25,
   "11":25,
   "12":25,
   "13":25,
   "14":25,
   "15":25,
   "16":25,
   "17":25,
   "18":25,
   "19":25,
   "20":25,
}


#nt_list=["4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"]
nt_list=["4","5","6","7","8","9","10","11","12","13","14"]
nr_list=["32","45","104"]
table_ensembles={}
for nt in nt_list:
    for nr in nr_list:
        mean={}
        height_mean={}
        rho_mean={}
        dens={}
        error={}
        for folder in folders:
            x=os.listdir(folder)
            beta=re.search(r"b\S{3}",folder).group().replace("b","")
            ls=a[beta]*int(nt)
            Vol=a[beta]*a[beta]*int(nr)*int(nr)
            for file_name in x:
                if "identification" + str(nt)+"nt" in file_name and str(t) +"t" in file_name and str(nr)+ "nr"  in file_name:
                    table_ensemble[file_name]={'beta':beta, 'a':a[beta],'ls':ls,'lt':ls, 'nr':nr, 'nt':nt, 'vol':vol, 'configurations':0, 'top Charge':0, 'fractionals':0, 'anti_fractionals':0, 'odds':0, 'top-frac':0, 'means':[], 'erros':[] 'histo_dens':[]}
                    #print(file_name)
                    
                    f=open(folder+"/"+file_name,"r")
                    norm=[]
                    rho=[]
                    height=[]

                    for line in f:  
                        splited=line.split()
                        q_top=float(splited[5])
                        norm_temp=extract_feature(line, 6)
                        rho_temp=extract_feature(line, 5)
                        height_temp=extract_feature(line,7)
                        frac=0
                        afrac=0
                        
                        for i in range(0,len(height_temp)):
                            if rho[i]>0.025 and norm_temp[i]>0.75:
                            norm.append(norm_temp[i])
                            rho.append(rho_temp[i])
                            height.append(height_temp[i])
                                if height_temp[i]>0:
                                    frac+=1
                                elif height_temp[i]<0:
                                    afrac+=1                          
                        if (frac_temp+afrac_temp)%2:
                            table_ensemble[file_name]['odds']+=1  
                        table_ensemble[file_name]['configurations']+=1
                        table_ensemble[file_name]['top-frac']+=abs(q_top-(frac/2-afrac/2))
                        table_ensemble[file_name]['histo_dens'].append(frac+afrac)
                    
                    norm=np.array((norm))
                    rho=np.array((rho))
                    height=np.array((height))
                    table_ensemble[file_name]['histo_dens']=np.array((table_ensemble[file_name]['histo_dens']))
                    
                    table_ensembles['means']=np.array((np.mean(table_ensemble[file_name]['histo_dens']),np.mean(norm),np.mean(rho),np.mean(height)))
                    table_ensembles['errors']=np.array((np.sqrt(np.std(table_ensemble[file_name]['histo_dens'])),np.sqrt(np.std(norm)),np.sqrt(np.std(rho)),np.sqrt(np.std(height))))

                    if hist==0:
                        #print(nr)
                        figure="./density_hist/hist_nt"+str(nt)+"b"+str(beta)+"t"+str(t)+"nr"+str(nr)+".png"
                        xlabel="density"
                        plot_histo(frac,bins,title,xlabel,figure,xrange=())
                        
                        figure="./norm_hist/norm_hist_nt"+str(nt)+"b"+str(beta)+"t"+str(t)+"nr"+str(nr)+".png"
                        xlabel="norm"
                        plot_histo(norm,bins,title,xlabel,figure,xrange=())
                        
                        figure="./height_hist/height_hist_nt"+str(nt)+"b"+str(beta)+"t"+str(t)+"nr"+str(nr)+".png"
                        xlabel="height"
                        plot_histo(height,bins,title,xlabel,figure,xrange=())
                        
                        figure="./rho_hist/rho_hist_nt"+str(nt)+"b"+str(beta)+"t"+str(t)+"nr"+str(nr)+".png"
                        xlabel="rho"
                        plot_histo(rho,bins,title,xlabel,figure,xrange=())
                f.close()

with open('scaling_'+str(norm_cut)+'.pkl','wb') as fp:
    pickle.dump(table_ensembles,fp)
    
    