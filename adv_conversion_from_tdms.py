# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 20:25:09 2023

@author: mb1536
"""

from nptdms import TdmsFile
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

path = "../../../Data_Files/MODAQ/Fall_22_Blade_Strain/Vector_Velocity_Data_2/"
files = glob(path + "*.tdms")

# allocation of the pre variety
ADV_tdms = {}
DATA = {}
Vec = [{} for _ in range(len(files))]

# reading all .tdms files, putting JUST the data into Vec
for i in range(0,len(files)):    
    ADV_tdms[i] = TdmsFile.read(files[i])
    DATA[i] = ADV_tdms[i]._channel_data
    indexed_keys = dict(enumerate(DATA[i].keys()))
    vec = {}
    for j in range(0,len(DATA[i])):
        vec[j] = DATA[i][indexed_keys[j]].data
    Vec[i] = vec

# concatenating the data so we have continuous time series      
Vector = {}
for k in Vec[0].keys():
  Vector[k] = np.concatenate(list(Vector[k] for Vector in Vec))

# just a quick check of x-velocity
plt.plot(Vector[0],Vector[1],'.')