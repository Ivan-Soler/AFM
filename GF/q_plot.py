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
        

#nt=sys.argv[1]
t=sys.argv[1]
hist=int(sys.argv[2])

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
            #nr=104
            L_s=a[beta]*int(nt)
            Vol=a[beta]*a[beta]*int(nr)*int(nr)
            for file_name in x:
                #print(file_name)
                if "identification" + str(nt)+"nt" in file_name and str(t) +"t" in file_name and str(nr)+ "nr"  in file_name:
                    #print(file_name)
                    f=open(folder+"/"+file_name,"r")
                    conf=[]
                    frac=[]
                    norm=[]
                    rho=[]
                    height=[]
                    deltaq=[]
                    inst=[]

                    N=0
                    odd=0

                    for line in f:  
                        frac_temp=0
                        afrac_temp=0
                        inst_temp=0
                        ainst_temp=0
                        norm_temp=[]
                        rho_temp=[]
                        height_temp=[]
                        
                        splited=line.split()
                        conf.append(int(splited[0]))
                        norm_temp=extract_feature(line, 6)
                        rho_temp=extract_feature(line, 5)
                        height_temp=extract_feature(line,7)
                        
                        
                        for i in range(0,len(height_temp)):
                            #if abs(height_temp[i])>0 and rho_temp[i]>2:
                            #if abs(height_temp[i])>0 and rho_temp[i]>1.5 and rho_temp[i]<12 : #for t=5 and nr=35,42
                            #if abs(height_temp[i])>0 and rho_temp[i]>0 and rho_temp[i]<12 and norm_temp[i]>0.75 and norm_temp[i]<2.25: #for t=20 and nr=104
                            if abs(height_temp[i])>lt_frac[nt] and abs(height_temp[i])<lt_inst[nt] and rho_temp[i]>2 and norm_temp[i]>norm_min[nt] and norm_temp[i]<norm_max[nt]:# for t=20 nr=104
                               #if abs(height_temp[i])>0.075 :
                                    #print(height_temp[i])
                                    #if height_temp[i]>lt_frac[nt]:
                                    #    inst_temp+=1
                                    #elif height_temp[i]<lt_frac[nt]:
                                    #    ainst_temp+=1   
                               #else:
                               norm.append(norm_temp[i])
                               rho.append(rho_temp[i])
                               height.append(height_temp[i])
                               if height_temp[i]>0:
                                    frac_temp+=1
                               elif height_temp[i]<0:
                                    afrac_temp+=1                           
                        if (frac_temp+afrac_temp)%2:
                            odd+=1         
                            #arclist=pd.read_pickle(directory+'dataprofile.pkl')
                            #flowtime=tau
                            #d=arclist[arclist['FlowTime'] == float(flowtime)]
                            #d.sort_values('ConfigNumber',inplace=True)
                            #try:
                            #   tar= tarfile.open(folder+"/"+"profile4dt"+str(t)+"c"+str(conf[i])+"to.dat")
                            #except PermissionError:
                            #   print("permission error for " + directory+fname)
                            #   continue   
                             
                        frac.append(frac_temp+afrac_temp)
                        inst.append(inst_temp+ainst_temp)
                        deltaq_temp=abs(float(splited[5])-(inst_temp+frac_temp/2-afrac_temp/2-ainst_temp))
                        #print(float(splited[5]), frac_temp, afrac_temp, inst_temp, ainst_temp)
                        #print("\n")
                        deltaq.append(deltaq_temp)
                        N+=1
                    #print(height)
                    mean_temp=sum(frac)
                    #print(mean_temp,sum(frac))
                    std2=0
                    if N>0:
                        for i in range(0,N):
                            std2+=(frac[i]-mean_temp/N)**2
                        std=np.sqrt(std2/N)
                        if (beta== "2.4" or beta=="2.5" or beta=="2.6" or beta=="2.7"):
                            height_mean[L_s]=0
                            count=0
                            for element in height_temp:
                                count+=1
                                height_mean[L_s]+=abs(element)
                            height_mean[L_s]=height_mean[L_s]/(count*a[beta]*a[beta])
                            rho_mean[L_s]=sum(rho_temp)*a[beta]/count
                            dens[L_s]=mean_temp/(N*Vol)
                            print(sum(frac),N,sum(frac)/N,mean_temp)
                            error[L_s]=std/(np.sqrt(N)*Vol)
                            #print(beta, nt, N, mean_temp/N, error[L_s])
                            if hist==0:
                                #print(nr)
                                plt.hist(frac,bins=100)
                                #plt.legend(loc="upper left")
                                plt.title("# conf="+str(N))
                                plt.savefig("./density_hist/dens_hist_nt"+str(nt)+"b"+str(beta)+"t"+str(t)+"nr"+str(nr)+".png")
                                plt.close()
                                
                                plt.hist(norm,bins=30,density=True)
                                #plt.hist(norm,bins=30,range=(0,2.75),density=True)
                                #plt.legend(loc="upper left")
                                plt.title("# structures="+str(np.array(frac).sum()))
                                plt.ylabel("frequency")
                                plt.xlabel("norm")
                                plt.savefig("./norm_hist/norm_hist_nt"+str(nt)+"b"+str(beta)+"t"+str(t)+"nr"+str(nr)+".png")
                                plt.close()
                                
                                plt.hist(height,bins=30)
                                #plt.hist(height,bins=60,range=(-0.1,0.1),density=True)
                                plt.title("# structures="+str(np.array(frac).sum()))
                                plt.ylabel("frequency")
                                plt.xlabel("height")
                                plt.savefig("./height_hist/height_hist_nt"+str(nt)+"b"+str(beta)+"t"+str(t)+"nr"+str(nr)+".png")
                                plt.close()
                                
                                #two dimensiona plot
                                #rhopy=np.array((rho))
                                #normpy=np.array((norm))
                                #data=np.dstack((rhopy,normpy))
                                #data=pd.DataFrame(data[0],columns=["rho","height"])
                                #tools.plot_with_seaborn(data,"./2d_hist/2d_hist_nt"+str(nt)+"b"+str(beta)+"t"+str(t)+"nr"+str(nr)+".png","# conf="+str(N))
                                
                                #plt.hist(rho,bins=30,density=True)
                                plt.hist(rho,bins=30,range=(0,3),density=False)
                                #plt.legend(loc="upper left")
                                plt.title("# conf="+str(N))
                                plt.ylabel("frequency")
                                plt.xlabel("rho")
                                plt.savefig("./rho_hist/rho_hist_nt"+str(nt)+"b"+str(beta)+"t"+str(t)+"nr"+str(nr)+".png")
                                plt.close()
                                
                                plt.hist(deltaq,density=True)
                                #plt.legend(loc="upper left")
                                plt.title("# conf="+str(N))
                                plt.savefig("./deltaq_hist/deltaq_hist_nt"+str(nt)+"b"+str(beta)+"t"+str(t)+"nr"+str(nr)+".png")
                                plt.close()
                        #print(beta, nt, N)
                            if beta=="2.4":
                                print(mean_temp,N,mean_temp/N,dens[L_s],nr,a[beta],Vol)
                    f.close()
                    #plt.plot(conf,frac,label=file_name.replace("identification","").replace(".txt",""))
        #plt.legend(ncol=2, loc="upper right")
        #plt.savefig("qf_hist_nt"+str(nt)+"t"+str(t)+".pdf")
        #plt.close()

        x=[]
        y=[]
        z=[]
        height_plot=[]
        #height_error_plot=[]
        rho_plot=[]
        #rho_error_plot=[]
        sort_dens=dict(sorted(dens.items()))
        for key in sort_dens:
            x.append(key)
            y.append(dens[key])
            height_plot.append(height_mean[key])
            rho_plot.append(rho_mean[key])
            z.append(error[key])
            #height_error_plot.append(height_error[key])
            #rho_error_plot.append(rho_error[key])
        if z and hist==1:
            plt.errorbar(x,y, yerr=z,label="nt="+nt +", nr="+nr, marker="o")
        if height_plot and hist==2:
            plt.errorbar(x,height_plot,label="nt="+nt +", nr="+nr, marker="o")
        if rho_plot and hist==3:
            plt.errorbar(x,np.sqrt(rho_plot),label="nt="+nt +", nr="+nr, marker="o")
            #plt.errorbar(x,y, yerr=z, marker="o")
    #plt.errorbar(x,y, yerr=z,marker="o")
#plt.ylim(0.2,1.2)
plt.legend(ncol=2, loc="lower right")
plt.xlabel("l_s(fm)")
plt.ylabel("size_frac**2(1/fm^2)")   
#plt.legend(ncol=2, loc="lower right")   
plt.savefig("scaling_size.pdf")
plt.close()
