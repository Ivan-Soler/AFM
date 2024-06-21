import pickle
import pandas as pd
import sys
file=sys.argv[1]

with open(file, 'rb') as f:
    data = pd.read_pickle(f)
print(data)
    
