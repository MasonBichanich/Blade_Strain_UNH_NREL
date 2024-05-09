# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio
from scipy.signal import find_peaks
from glob import glob
from datetime import datetime


## loading in the data
#path ='C:/Users/mason/OneDrive - USNH/Documents/Research/Data Analysis/NREL instumented blade/Data_Analysis_pats_structure/Data_Files/MODAQ/Vlink/Processed_Data/*'
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

## finding the peaks (troughs at sensor 3)
date = v_8145['Unnamed: 0'].to_numpy()
s3 = v_8145['2'].to_numpy()

peak_ind, _ = find_peaks(-v_8145['2'],distance=100,prominence=15)
q=v_8145['Unnamed: 0'][1]
time = np.zeros_like(v_8145['Unnamed: 0'])

for i in range(0,len(v_8145)):
    date = v_8145['Unnamed: 0'].iloc[i]
    date = date[:-3].strip()
    date = datetime.fromisoformat(date).timestamp()
    time[i] = date

#plt.figure()
#plt.plot(time,v_8145['2'])

#

pk_date = v_8145['Unnamed: 0'].iloc[peak_ind].to_numpy()
pk_st = v_8145['2'].iloc[peak_ind].to_numpy()

#plt.plot(time[peak_ind],pk_st,'r*')
#plt.show()


# remove the means
for i in np.arange(0,8):
    if i < 4:
        v_8145.iloc[:,i+1] = v_8145.iloc[:,i+1]-v_8145.iloc[:,i+1].mean()
    else:
        v_28175.iloc[:,i-3] = v_28175.iloc[:,i-3]-v_28175.iloc[:,i-3].mean()
        
#plt.figure()
#plt.plot(time,v_8145['2'])
#plt.show()
# Preallocation
broken = pd.DataFrame()
broken_time = pd.DataFrame()

# Initializing
for i in range(8):
    if i < 4:
        broken.iloc[i][0] = v_28175.iloc[i+1][0:peak_ind[0]]
    else:
        broken.iloc[i][0] = v_8145.iloc[i-3][0:peak_ind[0]]
    broken_time.iloc[0] = time[0:peak_ind[0]]

# Breaking the time
for i in range(0, len(peak_ind)):
    broken_time[i] = time[peak_ind[i]:peak_ind[i+1]]

# Breaking the data
for j in range(8):
    for i in range(1, len(peak_ind)):
        if j < 4:
            broken[j][i] = v_28175[peak_ind[i - 1]:peak_ind[i], j + 1]
        else:
            broken[j][i] = v_8145[peak_ind[i - 1]:peak_ind[i], j - 3]

# Now let's average these phases. They are all different lengths. First try will be to make them all the same length
# and Nan the missing data
# create a new compass, 0-360, and the length of each cell and use outer join

# Create imaginary compass
comp = [pd.DataFrame(np.round(np.linspace(0, 359, len(pd.DataFrame(broken[0][i])))), columns=[0]) for i in
        range(len(broken[0]))]

# Combine imaginary compass with broken data
broken_comp = [[] for _ in range(8)]
for j in range(8):
    for i in range(len(broken[0])):
        broken_comp[j].append(pd.concat([comp[i], broken[j][i]], axis=1))

# Combining them all into a DataFrame (each degree has a bunch of NaNs and data points to average over)
# (this takes a long time)
import time
start_time = time.time()

joined = [broken_comp[0][0] for _ in range(8)]
for j in range(8):
    for i in range(1, len(broken_comp[0])):
        joined[j] = pd.merge(broken_comp[j][i], joined[j], left_on=0, right_on=0, how='outer')

for j in range(8):
    joined[j] = joined[j].to_numpy()

end_time = time.time()
print("Elapsed time for combining data: ", end_time - start_time, " seconds")

# The data should be ensemble averaged BEFORE being phase averaged
dg = 1
degree = str(dg)
ens_ave = [[] for _ in range(8)]
for k in range(8):
    for i in range(len(joined[0])):
        for j in range(int(len(joined[0][:, 0]) / dg)):
            ens_ave[k].append(np.nanmean(joined[k][dg * j:dg * (j + 1), 1:], axis=0))

for i in range(8):
    ens_ave[i] = np.vstack(ens_ave[i])

# Plotting phase-averaged data
import matplotlib.pyplot as plt
cmap = plt.get_cmap('jet')
cmap = cmap(np.linspace(0, 1, 8) / 1.2)
plt.figure()

for i in range(8):
    phase_averaged_i = np.column_stack([ens_ave[i][:, 0], np.nanmean(ens_ave[i][:, 1:], axis=1)])
    means_i = np.nanmean(phase_averaged_i[:, 1])
    pa_a_i = phase_averaged_i[:, 1] - means_i

    pa_end_i = pa_a_i[int(1 / 4.5 * len(phase_averaged_i)):].tolist()
    pa_start_i = pa_a_i[:int(1 / 4.5 * len(phase_averaged_i) - 1)].tolist()
    pa_shifted_i = pa_end_i + pa_start_i

    plt.plot(phase_averaged_i[:, 0], pa_shifted_i, linewidth=1, color=cmap[i])

plt.axhline(0, color='black')
plt.legend(['SG-1', 'SG-2', 'SG-3', 'SG-4', 'SG-5', 'SG-6', 'SG-7', 'SG-8', 'y = 0'])
plt.grid(True)
plt.xlim([0, 360])
plt.xlabel('Degree of rotation')
plt.ylabel('Microstrain')
plt.show()


# Averaging only the similar length phases
new_bc = []
for i in range(1966):
    new_bc.append(len(broken_comp[2][i]) > 175 and len(broken_comp[2][i]) < 300)

njoined = broken_comp[2][0]
for i in range(1, len(broken_comp[2])):
    if new_bc[i]:
        njoined = pd.merge(broken_comp[2][i], njoined, left_on=0, right_on=0, how='outer')

njoined = njoined.to_numpy()

# Ensemble averaging for similar length phases
dg = 1
degree = str(dg)
nens_ave = []
for i in range(183):
    nens_ave_i = []
    for j in range(int(len(njoined[:, 0]) / dg)):
        nens_ave_i.append(np.nanmean(njoined[dg * j:dg * (j + 1), 1:], axis=0))
    nens_ave.append(np.vstack(nens_ave_i))

nens_ave = np.array(nens_ave)

# Phase averaging for similar length phases
nphase_averaged = np.column_stack([nens_ave[:, 0], np.nanmean(nens_ave[:, 1:], axis=1)])
nmeans = np.nanmean(nphase_averaged[:, 1])
npa_a = nphase_averaged[:, 1] - nmeans
npa_end = npa_a[int(1 / 4.5 * len(nphase_averaged)):].tolist()
npa_start = npa_a[:int(1 / 4.5 * len(nphase_averaged) - 1)].tolist()
npa_shifted = npa_end + npa_start

plt.figure()
plt.plot(nphase_averaged[:, 0], npa_shifted, linewidth=1)
plt.xlabel('Degree of rotation')
plt.ylabel('Microstrain')
plt.title('Phase-averaged data for similar length phases')
plt.grid(True)
plt.xlim([0, 360])
plt.axhline(0, color='black')
plt.show()

# Histogram of lengths of phases
lens = [len(broken_time[i]) for i in range(len(broken_time))]
plt.figure()
plt.hist(np.array(lens) / 64, bins=25)
plt.xlabel('Phase Length (Degrees)')
plt.ylabel('Frequency')
plt.title('Histogram of Lengths of Phases')
plt.show()
