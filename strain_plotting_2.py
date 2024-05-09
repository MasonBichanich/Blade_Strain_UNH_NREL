# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 16:31:52 2024

@author: mb1536
"""
import pickle
import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import my_funcs as my

with open('strain_plotting.pkl', 'rb') as f:
    [vec_date,mag,direction,str_date,strain,lc_date,thrust] = pickle.load(f)
    
amp = scipy.io.loadmat('amp.mat',squeeze_me = True)
amp=amp['amp']

amp_time=scipy.io.loadmat('amp_time.mat',squeeze_me = True)
amp_time=amp_time['amp_time']

amp_date = np.zeros_like(amp_time,dtype=object)
for i in range(0,len(amp_time)):
    amp_date[i] = my.matdatenum2datetime(amp_time[i])

for i in range(0,len(amp)):
    amp[i] = amp[i]/np.max(amp[i])
  
amp_n=np.zeros_like(amp_time,dtype=object)
for i in range(0,len(amp)):
    amp_n+=amp[i]
amp_n=amp_n/8

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
ax1.set_ylim([0,1])
# ax0.plot(vec_date,mag,'b') 
# ax1.plot(str_date,strain,'k')
# ax2.plot(vec_date,direction,'r')
# ax3.plot(str_date,strain,'k')
# ax4.plot(lc_date,thrust)
# ax5.plot(str_date,strain,'k')
ax0.plot(vec_date,mag,'b') 
ax1.plot(amp_date,amp_n,'k')
ax2.plot(vec_date,direction,'r')
ax3.plot(amp_date,amp_n,'k')
ax4.plot(lc_date,thrust)
ax5.plot(amp_date,amp_n,'k')

plt.figure()
ax = plt.subplot()
ax.plot(str_date,strain)
ax.set(ylabel='Microstrain', title='Strain Gauge 3')