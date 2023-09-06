# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio
from scipy.signal import find_peaks
from glob import glob


## loading in the data
#path ='C:/Users/mason/OneDrive - USNH/Documents/Research/Data Analysis/NREL instumented blade/Data_Analysis_pats_structure/Data_Files/MODAQ/Vlink/Processed_Data/*'
path = 'C:/Users/mb1536/OneDrive - USNH\Documents/Research/Data Analysis/NREL instumented blade/Data_Analysis_pats_structure/Data_Files/MODAQ/Vlink/Processed_Data/*'
files = sorted(glob(path))
v_8145 = pd.read_csv(files[0])
v_28175 = pd.read_csv(files[1])


## truncating the data to what we want
up_bound = 335729    # turbine stops rotating

v_8145 = v_8145.truncate(before=None,after=up_bound)
v_28175 = v_28175.truncate(before=None,after=up_bound)

## aligning the times
v_8145 = v_8145.drop(labels=None,axis=0,index=0)
v_28175 = v_28175.drop(labels=None,axis=0,index=up_bound)

## finding the peaks (troughs at sensor 3)
v_8145.plot.line(x='Unnamed: 0',y='2')
date=v_8145['Unnamed: 0'][:].to_numpy()
s3 = v_8145['2'][:].to_numpy()
plt.plot(date,s3)
peak_ind, _ = find_peaks(-v_8145['2'],distance=100,prominence=15)
v_8145.plot.scatter(x='Unnamed: 0',y='2')
pk_date = v_8145['Unnamed: 0'][peak_ind].to_numpy()
pk_st = v_8145['2'][peak_ind].to_numpy()
plt.plot(v_8145.iloc[:,0],v_8145.iloc[:,3])
plt.plot(v_8145['Unnamed: 0'],v_8145['2'])
plt.plot(v_8145['Unnamed: 0'][peak_ind],v_8145['2'][peak_ind],'r*')
plt.plot(pk_date,pk_st)
plt.figure()
plt.plot(v_8145.iloc[:,0],v_8145.iloc[:,3])
plt.show()
plt.plot(v_8145.iloc[peak_ind,0],v_8145.iloc[peak_ind,3],'r*')

# remove the means
for i in np.arange(0,8):
    if i < 4:
        v_8145.iloc[:,i+1] = v_8145.iloc[:,i+1]-v_8145.iloc[:,i+1].mean()
    else:
        v_28175.iloc[:,i-3] = v_28175.iloc[:,i-3]-v_28175.iloc[:,i-3].mean()
        

v_8145.plot.line(x='Unnamed: 0',y='2')

#does this even work