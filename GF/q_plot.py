import os
import matplotlib.pyplot as plt
import sys
import re
import numpy as np
import tools
import pandas as pd
import pickle
import jackknife

def extract_feature(line, param):
    splitted=line.split(" ")
    feature=[]
    for i in range(0,len(splitted)):
        if not (i-param-5) % 7:
            feature.append(float(splitted[i]))
     
    return(feature)

def plot_histo(data,bins,xlabel,figure_name,xrange=()):
    plt.hist(data,bins=100, density=True)
    title="# points:" +str(len(data))
    plt.ylabel("frequency")
    plt.xlabel(xlabel)
    plt.title(title)
    if xrange:
        plt.xlim(xrange)
    plt.savefig(figure_name)
    plt.close()

    return 
        
rho_min=float(sys.argv[1])
rho_max=float(sys.argv[2])
norm_min=float(sys.argv[3])
norm_max=float(sys.argv[4])

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

#nt_list=["4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"]
nt_list=["4","5","6","7","8","9","10","11","12","13","14"]
#nt_list=["15","16","17","18","19"]
nr_list=["32", "45", "104"]
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
            vol=a[beta]*a[beta]*int(nr)*int(nr)
            for file_name in x:
                if "identification" + str(nt)+"nt" in file_name and str(nr)+ "nr"  in file_name and float(beta)<2.75:
                    table_ensembles[file_name]={'beta':beta, 'a':a[beta],'ls':ls,'lt':ls, 'nr':nr, 'nt':nt, 'vol':vol, 'configurations':0, 'top Charge':0, 'fractionals':0, 'anti_fractionals':0, 'odds':0, 'top-frac':0, 'means':[0,0,0,0], 'erros':[0,0,0,0], 'histo_dens':[]}
          
                    f_dens=open("counting/count_b"+str(beta)+"nt"+str(nt)+".txt", "w")
                    f=open(folder+"/"+file_name,"r")
                    norm=[]
                    rho=[]
                    height=[]

                    for line in f:  
                        splited=line.split()
                        conf_number=int(splited[0])
                        q_top=float(splited[5])
                        norm_temp=extract_feature(line, 6)
                        rho_temp=extract_feature(line, 5)
                        height_temp=extract_feature(line,7)
                        frac=0
                        afrac=0
                        
                        for i in range(0,len(height_temp)):
                            if rho_temp[i]>rho_min and rho_temp[i]<rho_max and norm_temp[i]<norm_max and norm_temp[i]>norm_min:
                                norm.append(norm_temp[i])
                                rho.append(rho_temp[i])
                                height.append(height_temp[i])
                                if height_temp[i]>0:
                                    frac+=1
                                elif height_temp[i]<0:
                                    afrac+=1                          
                        if (frac+afrac)%2:
                            table_ensembles[file_name]['odds']+=1  
                        table_ensembles[file_name]['configurations']+=1
                        table_ensembles[file_name]['top-frac']+=abs(q_top-(frac/2-afrac/2))
                        table_ensembles[file_name]['histo_dens'].append([conf_number,frac+afrac])
                    
                    norm=np.array((norm))
                    rho=np.array((rho))
                    height=np.array((height))
                    table_ensembles[file_name]['histo_dens']=np.array((table_ensembles[file_name]['histo_dens']))
                    
                    #densities
                    dens, error=jackknife.jackknife_for_primary(np.array((table_ensembles[file_name]['histo_dens'][:,1])), int(1))
                    f_dens.write(str(beta)+" "+str(nt)+" " + str(ls) + "\n")
                    f_dens.write(str(dens)+" " +str(error) +"\n")
                    for element in table_ensembles[file_name]['histo_dens']:
                        f_dens.write(str(element[0])+" "+str(element[1])+"\n")
                    f_dens.close()
                    #dens=np.mean(np.array((table_ensembles[file_name]['histo_dens'])))
                    #error=np.sqrt(1/table_ensembles[file_name]['configurations']*np.std(np.array((table_ensembles[file_name]['histo_dens']))))
                    table_ensembles[file_name]['means']=np.array((dens,np.mean(norm),np.mean(rho),np.mean(height)))
                    table_ensembles[file_name]['errors']=np.array((error,np.sqrt(np.std(norm)/table_ensembles[file_name]['configurations']),np.sqrt(np.std(rho)/table_ensembles[file_name]['configurations']),np.sqrt(np.std(height)/table_ensembles[file_name]['configurations'])))

                    bins=100
                    figure="./density_hist/hist_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+".png"
                    xlabel="density"
                    plot_histo(table_ensembles[file_name]['histo_dens'],bins,xlabel,figure)

                    figure="./norm_hist/norm_hist_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+".png"
                    xlabel="norm"
                    xrange=()
                    plot_histo(norm,bins,xlabel,figure,xrange)

                    figure="./height_hist/height_hist_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+".png"
                    xlabel="height"
                    plot_histo(height,bins,xlabel,figure)

                    figure="./rho_hist/rho_hist_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+".png"
                    xlabel="rho"
                    plot_histo(rho,bins,xlabel,figure,xrange)

#with open('scaling_'+str(norm_cut)+'.pkl','wb') as fp:
with open('scaling.pkl','wb') as fp:
    pickle.dump(table_ensembles,fp)
    
    