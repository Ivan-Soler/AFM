import pickle

directory="/data/datatopoym/gbergner/pure_gauge_runs_with_twist/T2xR2_scan5/m-serv2/runns4nt104b2.6nr104/gf2/"
with open(directory+"dataprofile.pkl", 'rb') as f:
    data = pickle.load(f)
print(data)
    
