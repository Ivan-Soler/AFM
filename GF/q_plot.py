import os
import matplotlib.pyplot as plt
import sys
import re
import numpy as np

#nt=sys.argv[1]
#t=sys.argv[2]
lis_dir=os.listdir("./")
folders=[]
for direct in lis_dir:
    if "runns" in direct:
        folders.append(direct)

	
a={
	"2.4":0.12130,
	"2.5":0.08190,
	"2.6":0.05938,
	"2.7":0.04337,
	"2.8":0.03238,
	"2.9":0.02456
}


nt_list=["4","5","6","7","8","9","10"]
#nt_list=["4","6"]
t="20"
for nt in nt_list:
    mean={}
    dens={}
    error={}
    for folder in folders:
        x=os.listdir(folder)
        beta=re.search(r"b\S{3}",folder).group().replace("b","")
        L_s=a[beta]*int(nt)
        Vol=a[beta]*a[beta]*104*104
        for file_name in x:
            #print(file_name)
            if "identification" in file_name and str(nt)+"nt" in file_name and str(t) +"t" in file_name:
                f=open(folder+"/"+file_name,"r")
                conf=[]
                q_hist=[]
                mean_temp=0
                N=0
                for line in f:
                    splited=line.split()
                    conf.append(int(splited[0]))
                    q_hist.append(int(splited[1]) + int(splited[2]))
                    mean_temp+=int(splited[1]) + int(splited[2])
                    N+=1
                std2=0
                for i in range(0,len(q_hist)):
                    std2+=(q_hist[i]-mean_temp/N)**2
                std=np.sqrt(std2/N)
                if (beta=="2.4"or beta=="2.5" or beta=="2.6" or beta=="2.7"):
                    dens[L_s]=mean_temp/(N*Vol)
                    error[L_s]=std/(Vol)
                    print(beta, nt, N, mean_temp/N, error[L_s])
                #print(beta, nt, N)
                f.close()
                #plt.plot(conf,q_hist,label=file_name.replace("identification","").replace(".txt",""))
    #plt.legend(ncol=2, loc="upper right")
    #plt.savefig("qf_hist_nt"+str(nt)+"t"+str(t)+".pdf")
    #plt.close()

    x=[]
    y=[]
    z=[]
    sort_dens=dict(sorted(dens.items()))
    for key in sort_dens:
        x.append(key)
        y.append(dens[key])
        z.append(error[key])
    plt.errorbar(x,y, yerr=z,label="nt="+nt, marker="o")
    
#plt.ylim(0.2,1.2)
plt.xlabel("l_s(fm)")
plt.ylabel("rho_frac(1/fm^2)")   
plt.legend(ncol=2, loc="lower right")   
plt.savefig("scaling.pdf")
plt.close()
