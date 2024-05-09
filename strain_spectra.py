# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 17:38:39 2023

@author: mason
"""

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import my_funcs as my
import pickle


with open('strain_data.pkl', 'rb') as f:
    [strain] = pickle.load(f)
#%%
# replace = np.isfinite(ts) == 0
# ts[replace] = np.nanmean(ts)
# plt.plot(ts)
# k = 0
# for j in [v_8145,v_28175]:
#     for i in [1,2,3,4]:
#         S,f,DOF,CI = f_spectra(np.array(j.iloc[:,i]),1/64,10,2,0.05)
#         plt.figure()
#         plt.loglog(f,S,'k')
#         plt.xlim([10**-3,50])
#         plt.ylim([10**-4,10**4]) #microstrain square per Hz
#         if k == 0:
#             plt.title('Strain sensor '+str(i))
#         else:
#             plt.title('Strain sensor '+str(i+4))
#         plt.loglog(10,100,'.r')
#         plt.loglog([10,10],
#                      [CI[0]*100,CI[1]*100],
#                      'r')
#         plt.annotate("DOF = " + str(DOF),
#                      [10,200])
#         plt.grid()
#         plt.show()
#     k += 1


#%% fillling nans of time series
xi = np.arange(len(strain[8]))
mask = np.isfinite(strain[8])
filled = np.interp(xi,xi[mask],strain[8][mask])

#%% spectra

S,f,DOF,CI = my.f_spectra(filled,1/64,10,5,0.05)

fig,ax = plt.subplots()
ax.grid()
ax.loglog(f,S,'k')
ax.set(xlim = [10**-3,32],
       xlabel = 'Frequency [Hz]',
       ylim = [10**-4,10**4],
       ylabel = 'Energy density [ustrain^2/Hz]')
rot_f = 1.4922/4
rot_std = 0.1493
ax.annotate("Turbine rotation rate ="+str(rot_f),
              [.002,.01])
ax.axvspan(rot_f-rot_std,rot_f+rot_std,
            alpha = 0.2,
            color = 'red')
