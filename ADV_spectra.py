# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 09:08:32 2023

@author: mb1536
"""

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import my_funcs as my
import math


file = 'Vector_UTC_PostQC.mat'

ADV = scipy.io.loadmat(file,squeeze_me = True)
Vector = ADV['Vector']
xy_dir = np.arctan(Vector[:,1]/Vector[:,2])*180/3.1416
pos = xy_dir > 0 
xy_dir = xy_dir[pos]*-1

ts = Vector[142000:525000,1]
ts1 = xy_dir[142000:525000]
replace = np.isfinite(ts) == 0
replace1 = np.isfinite(ts1) == 0
ts[replace] = np.nanmean(ts)
ts1[replace1] = np.nanmean(ts1)
plt.plot(ts1)

S,f,DOF,CI = my.f_spectra(ts,1/64,10,10,0.05)
plt.loglog(f,S,'k')
plt.xlim([0,50])
plt.loglog(10,.1,'or')
plt.loglog([10,10],
            [CI[0]*0.1,CI[1]*0.1],
            'r')
plt.annotate("DOF = " + str(DOF),
              [10,.2])
plt.grid()

blade_f = 1.4922
blade_f_std = 0.1493
plt.annotate("Blade passage freq ="+str(blade_f),
              [10,0.3])
plt.axvspan(blade_f-2*blade_f_std,blade_f+2*blade_f_std,
            alpha = 0.2,
            color = 'red')
plt.axvspan(2*(blade_f-2*blade_f_std),2*(blade_f+2*blade_f_std),
            alpha = 0.2,
            color = 'green')


# S,f,DOF,CI = my.f_spectra(ts1,1/64,10,2,0.05)
# plt.loglog(f,S,'k')
# plt.xlim([0,50])
# plt.loglog(10,.1,'or')
# plt.loglog([10,10],
#            [CI[0]*0.1,CI[1]*0.1],
#            'r')
# plt.annotate("DOF = " + str(DOF),
#              [10,.2])
# plt.grid()