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
    mode=0
    splitted=line.split(" ")
    feature=[]
    for i in range(0,len(splitted)):
        if not (i-(param+mode)-6) % 16 and i>6:
            if splitted[i] != ' ' and splitted[i] != '\n':
                feature.append(float(splitted[i]))
    return(np.array((feature)))

def plot_histo(data,bins,xlabel,figure_name,xrange=()):
    plt.rcParams.update({'font.size': 18})
  
    if xrange:
      plt.hist(data,bins=100, range=xrange, density=True)
    else:
      plt.hist(data,bins=100, density=True)
    
    xs = np.linspace(0.2,8,200)
    #plt.vlines(x=[0.5,0.75], ymin=0, ymax=8, colors='black', ls='--', lw=2)

    #title="# points:" +str(len(data))
    #plt.axis._axinfo['label']['space_factor'] = 2.8
    plt.ylabel("Frequency \n")
    plt.xlabel(xlabel)
    #plt.title(title)

    plt.savefig(figure_name,dpi=800, bbox_inches='tight')
    plt.close()

    return 

def check_folder(nt,nr,beta,folder):
    fbeta=re.search(r"b\S{3}",folder).group().replace("b","")
    fnr=re.search(r"nt(.*?(?=b))",folder).group().replace("nt","")
    fnt=re.search(r"runns(.*?(?=nt))",folder).group().replace("runns","")

    if fnr==nr and fnt==nt and fbeta==beta:
        return(True)
    else:
        return(False)



rho_min=float(sys.argv[1])
rho_max=float(sys.argv[2])
norm_min=float(sys.argv[3])
norm_max=float(sys.argv[4])
histograms=bool(int(sys.argv[5]))
lis_dir=os.listdir("./")
folders=[]

for direct in lis_dir:
    if "runns" in direct:
        folders.append(direct)
a={
	"2.4":0.11530,
	"2.5":0.08194,
	"2.6":0.05938,
	"2.7":0.04337,
	"2.8":0.03238,
	"2.9":0.02456
}

nt_list=["4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"]
nr_list=["45", "64", "104"]
beta_list=["2.4","2.5","2.6","2.7"]
t_list=["10","15","25","35","50"]
t_list=["3.52", "7.34", "15", "27.37", "47.84"]
table_ensembles={}

deltaq=[]
oddity=[]

#norm_dic=pickle.load(open('norm','rb'))
#print(norm_dic)

for nt in nt_list:
    for nr in nr_list:
      for t in t_list:
        for beta in beta_list:
            ls=a[beta]*int(nt)
            vol=a[beta]*a[beta]*int(nr)*int(nr)
            key=nt+"nt"+nr+"nr"+beta+"b"+t+"tau"
            table_ensembles[key]={'beta':beta, 'a':a[beta],'ls':ls,'lt':ls, "t":t,'nr':nr, 'nt':nt, 'vol':vol, 'configurations':0, 'top Charge':0, 'fractionals':0, 'anti_fractionals':0, 'odds':0, 'top-frac':np.array((0,0),dtype=float), 'means':[0,0,0,0], 'erros':[0,0,0,0],'histo_dens':[],'pos_frac':{}, 'pos_afrac':{}}
            #mean={}
            #height_mean={}
            #rho_mean={}
            #dens={}
            #error={}
            histo=[]
            for folder in folders:
                if check_folder(nt,nr,beta,folder):
                    x=os.listdir(folder)
                    for file_name in x:
                        if "identification" in file_name and "nr"+t+"t" in file_name:# and key in norm_dic:
                            f=open(folder+"/"+file_name,"r")
                            for line in f:  
                                splited=line.split()
                                conf_number=int(splited[0])
                                #conf_number=int(re.search(r"dt(.*?(?=c))",splited[0]).group().replace("dt",""))
                                table_ensembles[key]['pos_frac'][conf_number]=[]
                                table_ensembles[key]['pos_afrac'][conf_number]=[]
                                x=extract_feature(line, 1)
                                y=extract_feature(line, 2)
                                q_top=float(splited[5])
                                norm_temp=extract_feature(line, 6)
                                rho_temp=extract_feature(line, 5)
                                height_temp=extract_feature(line,8)
                                duality_temp=extract_feature(line,9)
                                cov_temp=extract_feature(line,7)
                                
                                
                                frac=0
                                afrac=0
                                inst=0
                                ainst=0
                                dinst=0
                                dainst=0

                                for i in range(0,len(height_temp)):
                                    if rho_temp[i]>rho_min and rho_temp[i]<rho_max and norm_temp[i]>norm_min and norm_temp[i]<norm_max and cov_temp[i]<1000:
                                        histo.append([norm_temp[i],rho_temp[i],abs(height_temp[i]),duality_temp[i],cov_temp[i]])
                                        if height_temp[i]>0:
                                            frac+=1
                                            table_ensembles[key]['pos_frac'][conf_number].append([x[i],y[i]])
                                        elif height_temp[i]<0:
                                            afrac+=1   
                                            table_ensembles[key]['pos_afrac'][conf_number].append([x[i],y[i]])
                                    elif False: #rho_temp[i]>rho_min and norm_temp[i]>3:
                                        if height_temp[i]>0:
                                            inst+=1
                                            #posi.append([x[i],y[i]])
                                        elif height_temp[i]<0:
                                            ainst+=1   
                                            #posai.append([x[i],y[i]])
                                            
                                if (frac+afrac)%2:
                                    table_ensembles[key]['odds']+=1  
                                table_ensembles[key]['configurations']+=1
                                table_ensembles[key]['top-frac']+=np.array((abs(q_top-(inst+frac/2-ainst-afrac/2)),abs(q_top-(frac/2-afrac/2))))
                                table_ensembles[key]['histo_dens'].append([conf_number,frac+afrac,inst+ainst])
                            #End reading line of each identification file 
                            f.close()
                        #End if for identification file
                    #End loop files in folder
                #End if check folder
            #End loop folders in directory  
            if table_ensembles[key]['configurations']<10:
              del table_ensembles[key]
            else:
                table_ensembles[key]['histo_dens']=np.array((table_ensembles[key]['histo_dens']))
                table_ensembles[key]['odds']=table_ensembles[key]['odds']/table_ensembles[key]['configurations']*100
                table_ensembles[key]['top-frac']=table_ensembles[key]['top-frac']/table_ensembles[key]['configurations']
                
                #densities
                oddity.append(table_ensembles[key]['odds'])
                dens, error=jackknife.jackknife_for_primary(np.array((table_ensembles[key]['histo_dens'][:,1])), int(1))
                f_dens=open("./counting/count_"+key+".txt","w")
                f_dens.write(str(beta)+" "+str(nt)+" " + str(ls) + " "+str(table_ensembles[key]['odds'])+"\n")
                f_dens.write(str(dens)+" " +str(error) +"\n")
                for element in table_ensembles[key]['histo_dens']:
                    f_dens.write(str(element[0])+" "+str(element[1])+"\n")
                f_dens.close()
            
                histo=np.array((histo))
                norm=np.array((histo[:,0]))
                rho=np.array((histo[:,1]))
                height=np.array((histo[:,2]))
                duality=np.array((histo[:,3]))
                cov=np.array((histo[:,4]))
                cov=np.nan_to_num(cov, copy=False, nan=0, posinf=0, neginf=0)
                table_ensembles[key]['means']=np.array((dens,np.mean(norm),np.mean(rho),np.mean(height)))
                table_ensembles[key]['errors']=np.array((error,np.std(norm)/np.sqrt(table_ensembles[key]['configurations']),np.std(rho)/np.sqrt(table_ensembles[key]['configurations']),np.std(height)/np.sqrt(table_ensembles[key]['configurations'])))
    
                if histograms:
                    bins=100
                    xrange=(-100,100)
                    figure="./fractionals/frac_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+"t"+str(t)+".png"
                    xlabel="density"
                    plot_histo(table_ensembles[key]['histo_dens'][1],bins,xlabel,figure)

                    figure="./instantons/inst_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+"t"+str(t)+".png"
                    xlabel="density"
                    plot_histo(table_ensembles[key]['histo_dens'][2],bins,xlabel,figure)

                    xrange=(0,2.0)
                    figure="./norm_hist/norm_hist_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+"t"+str(t)+".png"
                    xlabel="$Q_{fit}$"
                    plot_histo(norm,bins,xlabel,figure,xrange)

                    #xrange=(0,0.02)
                    figure="./height_hist/height_hist_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+"t"+str(t)+".png"
                    xlabel="height"
                    plot_histo(height,bins,xlabel,figure)
    
                    figure="./rho_hist/rho_hist_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+"t"+str(t)+".png"
                    xlabel="rho"
                    plot_histo(rho,bins,xlabel,figure)

                    figure="./duality/duality_his_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+"t"+str(t)+".png"
                    xlabel="duality"
                    plot_histo(duality,bins,xlabel,figure)

                    xrange=(0,3)
                    figure="./cov/cov_his_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+"t"+str(t)+".png"
                    xlabel="covariance"
                    plot_histo(cov,bins,xlabel,figure,xrange)
        
#with open('scaling_'+str(norm_cut)+'.pkl','wb') as fp:
#with open('scaling_s'+str(norm_scale)+'.pkl','wb') as fp:
#    pickle.dump(table_ensembles,fp)

norm_dic={}
for key in table_ensembles:
  norm_dic[key]=table_ensembles[key]['means'][1]
  print(key, norm_dic[key])
  
with open('norm','wb') as fp:
  pickle.dump(norm_dic,fp)
    
    