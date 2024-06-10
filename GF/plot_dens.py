import pickle 
import matplotlib.pyplot as plt
import numpy as np
import tools
import sys
import tarfile
import os
import pandas as pd
import re


config=sys.argv[1]
param=tools.parse_in_file(config)
directory=param[0]
tau=param[17]

param=re.search(r"runns.+?(?=/)",directory).group().replace("runns","")
beta=re.search(r"b\S{3}",param).group().replace("b","")
nr=re.search(r"nt(.+?(?=b))",param).group().replace("nt","")
nt=re.search(r"(.+?(?=nt))",param).group().replace("nt","")

arclist=pd.read_pickle(directory+'dataprofile.pkl')
flowtime=tau
d=arclist[arclist['FlowTime'] == float(flowtime)]
d.sort_values('ConfigNumber',inplace=True)

dcount="/data/datatopoym/iscalero/T2xR2/counting/"
table_ensembles=pickle.load(open('../scaling.pkl', "rb" ))
key=nt+"nt"+nr+"nr"+beta+"b"

for index,row in d.iterrows():
    fname=row['ArchiveName']
    filen=row['FileName']
    if "to.dat" in filen and "dt"+str(tau) in filen:
        try:
            tar= tarfile.open(directory+fname)
        except PermissionError:
            error.write("permission error for " + directory+fname+"\n")
            continue

        conf=filen.replace("profile4dt"+str(tau)+"c", "")
        conf=int(conf.replace("to.dat", ""))

        if conf in table_ensembles[key]['pos_frac']:
            frac=table_ensembles[key]['pos_frac'][conf]
            afrac=table_ensembles[key]['pos_afrac'][conf]
            
            tar.extract(filen)
            top_density,sizes=tools.read_top(filen)
            density_2d_top,sizes_big,index_smal=tools.projection_2d(top_density,sizes)
            tools.plot_dens_2d(filen,density_2d_top,sizes_big, frac, afrac)
            os.remove(filen)
