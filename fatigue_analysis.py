# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 06:50:41 2024

@author: mb1536
"""
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import pickle
import datetime as dt
import numpy as np
import matplotlib.dates as mdates
import rainflow
import pandas as pd

with open('strain_data.pkl', 'rb') as f:
    [strain] = pickle.load(f)
    
 
# datetimes = np.zeros_like(strain[0])
# for i in range(0,len(strain[0])-1):
#     datestr = strain[0][i+1]
#     datestr = datestr[:-3].strip()
#     datestr = dt.datetime.strptime(datestr,'%Y-%m-%d %H:%M:%S.%f')
#     datetimes[i] = datestr

#%% Convert strain to stress
E = 68.9e3    # for 6061 aluminum in MPa
sig = strain[3]*1e-6*E     #factor of 1e-6 because microstrain. 1e6-6 =/= 10e6. dumbass
sig = sig.to_numpy()
#%%
# peak_ind, _ = find_peaks(-sig,distance=100,prominence=15)
# pos_peak_ind, _ = find_peaks(sig,distance=100,prominence=15)
     
# mod = peak_ind + 1
# pos_mod = pos_peak_ind + 1
# pk_date = datetimes[peak_ind]
# pk_s = sig[mod]
# pos_pk_date = datetimes[pos_peak_ind]
# pos_pk_s = sig[pos_mod]

# plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
# plt.plot(datetimes,sig)
# plt.plot(pk_date,pk_s,'r*')
# plt.plot(pos_pk_date,pos_pk_s,'k*')
# plt.grid('on')

#%%
sig[-1] = sig[-2]


    
cycles = rainflow.count_cycles(sig,binsize = .125)
cycles = cycles[20:]
s = np.empty(1)
for i in range(0,len(cycles)):
    s = np.append(s,cycles[i][0])

s = s[1:]

# put it in a dataframe
df = pd.DataFrame(cycles, columns = ['Range [MPa]','Cycles'])
df.plot(kind = 'bar', x = 'Range [MPa]')
plt.ylabel('Counts')



#%% damage calculation
R = -0.8
s_ksi = 0.145038*s
s_eq = s_ksi*(1-R)**0.63
N = 10**(20.68 - 9.84*np.log10(s_eq))

D_i = N.copy()     # preallocation of damage
for i in range(0,len(cycles)):
    D_i[i] = cycles[i][1]/N[i]
    
D = sum(D_i)

life_yrs = 90/D/60/24/365
