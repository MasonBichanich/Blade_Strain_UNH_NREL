# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 10:13:21 2023

@author: mason
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
from glob import glob
from PIL import Image
import pandas as pd
import matplotlib.dates as mdates
import datetime as dt
with open('strain_data.pkl', 'rb') as f:
    [strain] = pickle.load(f)
    
    
const = [0.005127,
         0.332171,
         -0.02996,
         4.621859,
         3.501752,
         0.199047,
         0.296224,
         1.685407]
coef = [-0.02675,
        -.03087,
        -.03112,
        -0.01909,
        -0.01851,
        -0.01874,
        -0.03016,
        -0.02539]

moment = pd.DataFrame(index=strain.index,columns=strain.columns)
time = strain.iloc[:,0]
for i in range(0,len(time)-1):
    temp = time[i+1][:-3].strip()
    moment.iloc[i,0] = dt.datetime.strptime(temp,'%Y-%m-%d %H:%M:%S.%f')
# for i in range(0,len(v_8145)):
#     date = v_8145['Unnamed: 0'].iloc[i]
#     date = date[:-3].strip()
#     date = datetime.strptime(date,'%Y-%m-%d %H:%M:%S.%f')
#     time[i] = date

#moment.iloc[:,0] = [dt.datetime.strptime(d,'%Y-%m-%d %H:%M:%S.%f').date() for d in moment.iloc[:,0]]

for i in range(0,8):
    moment.iloc[:,i+1] = np.array((strain.iloc[:,i+1]-const[i]) / coef[i])

plt.plot(moment.iloc[:,0],moment.iloc[:,7])

#%%

# x_loc = [4.114,9.114,14.114,27.5,33.5,39.5,52.886,62.886]
# for i in range(len(moment)):

#     fig,ax = plt.subplots()
#     ax = plt.plot(x_loc,moment[i,:])
#     #ax.set(xlim=(0,75))
#     plt.ylim([-5000,4500])
#     plt.grid()
#     #ax.set(ylabel = 'Internal moment [lb-in]',xlabel = 'Location along blade [in]')   
#     fig.savefig(f'moment_profiles/{i:00005}')
    
#     plt.close()

# #%%
# imgs = glob("moment_profiles/*.png")
# frames = []
# for i in imgs:
#     new_frame = Image.open(i)
#     frames.append(new_frame)
    
# frames[0].save('moment_profile.gif', format='GIF',
#                append_images=frames[1:],
#                save_all=True,
#                duration=20,
#                Loop=0)
