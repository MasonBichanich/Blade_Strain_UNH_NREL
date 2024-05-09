from nptdms import TdmsFile
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import datetime
import matplotlib.dates as mdate

path = "../../../Data_Files/MODAQ/V_link/8145_Vlink/Raw Data/"
files = glob(path + "*.tdms")

# allocation of the pre variety
ADV_tdms = {}
DATA = {}
Vec = [{} for _ in range(len(files))]

# reading all .tdms files, putting JUST the data into Vec
for i in range(0,len(files)): 
    try:
        ADV_tdms[i] = TdmsFile.read(files[i])
        DATA[i] = ADV_tdms[i]._channel_data
        indexed_keys = dict(enumerate(DATA[i].keys()))
        vec = {}
        for j in range(0,len(DATA[i])):
            vec[j] = DATA[i][indexed_keys[j]].data
        Vec[i] = vec
    except:
        print("Bad file at",files[i])
        
# concatenating the data so we have continuous time series      
Vector = {}

for k in Vec[0].keys():
    Vector[k] = np.concatenate(list(Vector[k] for Vector in Vec))
strain = Vector
# just a quick check of x-velocity
# plt.plot(strain[1],strain[2],'.')
# plt.show()
epoch_time = strain[1]/(10**9)
# time = [None]*len(epoch_time)
# for i in range(0,len(epoch_time)):
#     time[i] = datetime.datetime.fromtimestamp((epoch_time[i])).strftime('%Y-%m-%d %H:%M:%S')
time = mdate.epoch2num(epoch_time)

fig,ax = plt.subplots()
plt.plot_date(time,strain[3],'.k')
date_fmt = '%d-%m-%y %H:%M:%S'
date_formatter = mdate.DateFormatter(date_fmt)
ax.xaxis.set_major_formatter(date_formatter)
fig.autofmt_xdate()



with open('8145vlink_data.pkl','wb') as f:
    pkl.dump([strain],f)