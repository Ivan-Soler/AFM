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
        if not (i-param-5) % 8 and i>5:
            if splitted[i] != ' ' and splitted[i] != '\n':
                feature.append(float(splitted[i]))
    return(np.array((feature)))

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
norm_scale=float(sys.argv[3])
print(norm_scale)
#norm_min=float(sys.argv[3])
#norm_max=int(sys.argv[4])
histograms=bool(int(sys.argv[4]))
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
table_ensembles={}

deltaq=[]
oddity=[]

#norm_dic=pickle.load(open('norm','rb'))
#print(norm_dic)

for nt in nt_list:
    for nr in nr_list:
        for beta in beta_list:
            ls=a[beta]*int(nt)
            vol=a[beta]*a[beta]*int(nr)*int(nr)
            key=nt+"nt"+nr+"nr"+beta+"b"
            table_ensembles[key]={'beta':beta, 'a':a[beta],'ls':ls,'lt':ls, 'nr':nr, 'nt':nt, 'vol':vol, 'configurations':0, 'top Charge':0, 'fractionals':0, 'anti_fractionals':0, 'odds':0, 'top-frac':np.array((0,0),dtype=float), 'means':[0,0,0,0], 'erros':[0,0,0,0],'histo_dens':[],'pos_frac':{}, 'pos_afrac':{}}
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
                        if "identification" in file_name:# and key in norm_dic:
                            f=open(folder+"/"+file_name,"r")
                            for line in f:  
                                splited=line.split()
                                conf_number=int(splited[0])
                                table_ensembles[key]['pos_frac'][conf_number]=[]
                                table_ensembles[key]['pos_afrac'][conf_number]=[]
                                x=extract_feature(line, 1)
                                y=extract_feature(line, 2)
                                q_top=float(splited[5])
                                norm_temp=extract_feature(line, 6)
                                rho_temp=extract_feature(line, 5)
                                height_temp=extract_feature(line,7)
                                duality_temp=np.pi*extract_feature(line,8)
                                
                                
                                frac=0
                                afrac=0
                                inst=0
                                ainst=0
                                dinst=0
                                dainst=0

                                for i in range(0,len(height_temp)):
                                    if True:# rho_temp[i]>rho_min and norm_temp[i] > 0.5 and norm_temp[i]<3:
                                        histo.append([norm_temp[i],rho_temp[i],abs(height_temp[i]),duality_temp[i]])
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
                                table_ensembles[key]['histo_dens'].append([conf_number,frac+afrac])
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
                table_ensembles[key]['means']=np.array((dens,np.mean(norm),np.mean(rho),np.mean(height)))
                table_ensembles[key]['errors']=np.array((error,np.std(norm)/np.sqrt(table_ensembles[key]['configurations']),np.std(rho)/np.sqrt(table_ensembles[key]['configurations']),np.std(height)/np.sqrt(table_ensembles[key]['configurations'])))
    
                if histograms:
                    bins=100
                    figure="./density_hist/hist_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+".png"
                    xlabel="density"
                    plot_histo(table_ensembles[key]['histo_dens'],bins,xlabel,figure)
    
                    figure="./norm_hist/norm_hist_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+".png"
                    xlabel="norm"
                    xrange=()
                    plot_histo(norm,bins,xlabel,figure)

                    #xrange=(0,0.02)
                    figure="./height_hist/height_hist_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+".png"
                    xlabel="height"
                    plot_histo(height,bins,xlabel,figure,xrange)
    
                    figure="./rho_hist/rho_hist_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+".png"
                    xlabel="rho"
                    plot_histo(rho,bins,xlabel,figure)

                    figure="./duality/duality_his_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+".png"
                    xlabel="duality"
                    plot_histo(duality,bins,xlabel,figure)

                    figure="./cov/cov_his_nt"+str(nt)+"b"+str(beta)+"nr"+str(nr)+".png"
                    xlabel="covariance"
                    plot_histo(cov,bins,xlabel,figure)
        
#with open('scaling_'+str(norm_cut)+'.pkl','wb') as fp:
#with open('scaling_s'+str(norm_scale)+'.pkl','wb') as fp:
#    pickle.dump(table_ensembles,fp)

norm_dic={}
for key in table_ensembles:
  norm_dic[key]=table_ensembles[key]['means'][3]
  print(key, norm_dic[key])
  
with open('norm','wb') as fp:
  pickle.dump(norm_dic,fp)
    
    