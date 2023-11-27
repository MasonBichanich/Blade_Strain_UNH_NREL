# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 17:38:39 2023

@author: mason
"""

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
from my_funcs import *
import pickle


with open('strain_data.pkl', 'rb') as f:
    [v_8145,v_28175] = pickle.load(f)
#%%
# replace = np.isfinite(ts) == 0
# ts[replace] = np.nanmean(ts)
# plt.plot(ts)
k = 0
for j in [v_8145,v_28175]:
    for i in [1,2,3,4]:
        S,f,DOF,CI = f_spectra(np.array(j.iloc[:,i]),1/64,10,2,0.05)
        plt.figure()
        plt.loglog(f,S,'k')
        plt.xlim([10**-3,50])
        plt.ylim([10**-4,10**4]) #microstrain square per Hz
        if k == 0:
            plt.title('Strain sensor '+str(i))
        else:
            plt.title('Strain sensor '+str(i+4))
        plt.loglog(10,100,'.r')
        plt.loglog([10,10],
                     [CI[0]*100,CI[1]*100],
                     'r')
        plt.annotate("DOF = " + str(DOF),
                     [10,200])
        plt.grid()
        plt.show()
    k += 1


