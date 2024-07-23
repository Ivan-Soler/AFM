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
tau=sys.argv[1]

param=re.search(r"runns.+?(?=/)",directory).group().replace("runns","")
beta=re.search(r"b\S{3}",param).group().replace("b","")
nr=re.search(r"nt(.+?(?=b))",param).group().replace("nt","")
nt=re.search(r"(.+?(?=nt))",param).group().replace("nt","")


tarfiles=os.listdir("profiles")
for directory in tarfiles:
      try:
          tar= tarfile.open("/data/datatopoym/iscalero/T2xR2/flow_13ns/runns13nt64b2.6nr64/profiles/"+directory)
      except Exception as e:
          error.write("permission error for " + directory+fname+"\n")
          continue
      list_files=tar.getmembers()
      for profile in list_files:
        filen=profile.name
        print(filen)
        if "to" in filen and "4dt"+str(tau)+"c" in filen:
            conf=filen.replace("profile4dt"+str(tau)+"c", "")
            conf=filen.replace("to.dat", "")

            en_file=filen.replace("to","en")
            tar.extract(filen)
            tar.extract(en_file)

            top_density,sizes=tools.read_top(filen)
            density_2d_top,sizes_big,index_smal=tools.projection_2d(top_density,sizes)
            tools.plot_dens_2d(filen,density_2d_top,sizes_big, frac, afrac)
            os.remove(filen)
