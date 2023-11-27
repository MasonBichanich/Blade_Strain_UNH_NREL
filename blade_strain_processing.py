# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 13:33:00 2023

@author: mb1536
"""

#phase averaging done in matlab


#from nptdms import TdmsFile
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import mat73
import math
import scipy.io
import pandas as pd
from datetime import datetime
#from my_funcs import *
import pickle

# path = "../../../../V_link/8145_Vlink/Raw Data/"
# files = glob(path + "*.tdms")

# # allocation of the pre variety
# vlink_tdms = {}
# DATA = {}
# Vlink = [{} for _ in range(len(files))]
# # reading all .tdms files, putting JUST the data into Vec
# for i in range(0,len(files)): 
#     try:
#         vlink_tdms[i] = TdmsFile.read(files[i])
#         DATA[i] = vlink_tdms[i]._channel_data
#         indexed_keys = dict(enumerate(DATA[i].keys()))
#         vlink = {}
#         for j in range(0,len(DATA[i])):
#             vlink[j] = DATA[i][indexed_keys[j]].data
#         Vlink[i] = vlink
#     except:
#         print("Bad file")
# # concatenating the data so we have continuous time series      
# vlink_8145 = {}
# for k in Vlink[2].keys():
#   vlink_8145[k] = np.concatenate(list(vlink_8145[k] for vlink_8145 in Vlink))
#%%
file = 'Vector_UTC_PostQC.mat'

ADV = scipy.io.loadmat(file,squeeze_me = True)
Vector = ADV['Vector']

# %%
heading=105.05
INS_X_VEL_heading_correction = -Vector[:,1]*math.sin(np.deg2rad(heading-90))-Vector[:,2]*math.sin(np.deg2rad(180-heading))
INS_Y_VEL_heading_correction =  Vector[:,1]*math.cos(np.deg2rad(heading-90))-Vector[:,2]*math.cos(np.deg2rad(180-heading))

speed_complex=(INS_X_VEL_heading_correction +INS_Y_VEL_heading_correction*1j)
direc_rad = np.angle(speed_complex)
Speed_Direction = np.array([np.abs(speed_complex),np.rad2deg(direc_rad)])

Current_Magnitude = np.transpose(Speed_Direction[0,:])
Current_Direction = np.transpose(Speed_Direction[1,:])

flood_i = Current_Direction > 0
Current_Magnitude[flood_i] = -1*Current_Magnitude[flood_i]
#%%
path = '../../../Data_Files/MODAQ/Vlink/Processed_Data/*'
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

time = np.zeros_like(v_8145['Unnamed: 0'])

# for i in range(0,len(v_8145)):
#     date = v_8145['Unnamed: 0'].iloc[i]
#     date = date[:-3].strip()
#     date = datetime.fromisoformat(date).timestamp()
#     time[i] = date
for i in range(0,len(v_8145)):
    date = v_8145['Unnamed: 0'].iloc[i]
    date = date[:-3].strip()
    date = datetime.strptime(date,'%Y-%m-%d %H:%M:%S.%f')
    time[i] = date
    
for i in np.arange(0,8):
    if i < 4:
        v_8145.iloc[:,i+1] = v_8145.iloc[:,i+1]-v_8145.iloc[:,i+1].mean()
    else:
        v_28175.iloc[:,i-3] = v_28175.iloc[:,i-3]-v_28175.iloc[:,i-3].mean()
 
with open('strain_data.pkl','wb') as f:
    pickle.dump([v_8145,v_28175],f)
#%%
Vector_date = np.zeros_like(Vector[:,0],dtype=object)
for i in range(0,len(Vector[:,0])):
    Vector_date[i] = (datenum2datetime(Vector[i,0]))

#%%
spin_v = np.logical_and(Vector_date >= time[0],
                        Vector_date <= time[-1])
#%%

LCtime = scipy.io.loadmat('thrust_datenum.mat')
thrust_datenum = LCtime['LC_date']

LCforce = scipy.io.loadmat('thrust.mat')
thrust = LCforce['Thrust_Force']

LC_date = np.zeros_like(thrust_datenum,dtype=object)
for i in range(0,len(thrust_datenum)):
    LC_date[i] = datenum2datetime(thrust_datenum[i,0])


#%%
spin_thr = np.logical_and(LC_date >= time[0],
                        LC_date <= time[-1])                       
#%%
ax0 = plt.subplot(311)
ax1 = ax0.twinx()
ax2 = plt.subplot(312)
ax3 = ax2.twinx()
ax4 = plt.subplot(313)
ax5 = ax4.twinx()
ax1.set_ylabel('Microstrain')
ax3.set_ylabel('Microstrain')
ax5.set_ylabel('Microstrain')
ax1.get_shared_y_axes().join(ax1,ax3,ax5)
ax0.get_shared_x_axes().join(ax0,ax2,ax4)
ax0.set_ylabel('Current Magnitude [m/s]')
ax2.set_ylabel('Current Direction')
ax4.set_ylabel('Turbine Thrust [kN]')
ax1.set_ylim([-70,60])
ax0.plot(Vector_date[spin_v],Current_Magnitude[spin_v],'b') 
ax1.plot(time,v_28175['2'],'k')
ax2.plot(Vector_date[spin_v],Current_Direction[spin_v],'r')
ax3.plot(time,v_28175['2'],'k')
ax4.plot(LC_date[spin_thr],thrust[spin_thr])
ax5.plot(time,v_28175['2'],'k')


