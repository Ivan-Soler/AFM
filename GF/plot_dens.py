import pickle 
import matplotlib.pyplot as plt
import numpy as np
import tools
import sys
import tarfile
import os
import pandas as pd

directory=param[0]
param=re.search(r"runns(.+?(?=/))",directory).group().replace("runns","")
beta=re.search(r"b\S{3}",param).group().replace("b","")
nt=re.search(r"nt(.+?(?=b))",param).group().replace("nt","")
nr=re.search(r"nr(.+?(?=/)",param).group().replace("nr","")

arclist=pd.read_pickle(directory+'dataprofile.pkl')
flowtime=tau
d=arclist[arclist['FlowTime'] == float(flowtime)]
d.sort_values('ConfigNumber',inplace=True)

dcount="/data/datatopoym/iscalero/T2xR2/counting/"
for index,row in d.iterrows():
    fname=row['ArchiveName']
    filen=row['FileName']
    if "to.dat" in filen and "dt"+str(tau) in filen:
        try:
            tar= tarfile.open(directory+fname)
            fcount=open(dcount+"count_b"+str(beta)+"nt"+str(nt)+"nr"+str(nr))
        except PermissionError:
            error.write("permission error for " + directory+fname+"\n")
            continue
            
    conf=filen.replace("profile4dt"+str(tau)+"c", "")
    conf=conf.replace("to.dat", "")
    
    for line in fcount:
        splitted=line.split(" ")
        if splitted[0]==conf and int(splitted[1])%2:
            if "identification" + str(nt)+"nt" in file_name and str(nr)+ "nr"  in file_name and float(beta)<2.75:
                frac=np.loadtxt(dcount+"position_fractionals"+str(beta)+"nt"+str(nt)+"nr"+str(nr)+".txt")
                afrac=np.loadtxt(dcount+"position_anti_fractionals"+str(beta)+"nt"+str(nt)+"nr"+str(nr)+".txt")
                tar.extract(filen)
                top_density,sizes=tools.read_top(filen)
                density_2d_top,sizes_big,index_smal=tools.projection_2d(top_density,sizes)
                tools.plot_dens_2d(filen,density_2d_top,sizes_big, frac, afrac)
                continue
    fcount.close()