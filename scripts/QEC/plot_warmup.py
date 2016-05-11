import numpy as np
import os,sys
from matplotlib import pyplot as plt
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
import time
import matplotlib.cm as cm

#Import data from previous warmups


def get_warmup_folders(number_of_warmups = 2):
    i=0
    folder_array=[]
    timestamp=toolbox.get_timestamp_from_now()
    while i < number_of_warmups:
        warmupdata_folder=toolbox.latest_data(contains = 'monitor_warmup', older_than= timestamp )
        folder_array.append(warmupdata_folder)
        [date,time]=toolbox.get_date_time_string_from_folder(warmupdata_folder)
        timestamp=  date+'_'+time
        i+=1
    return folder_array


def plt_warmup_data(folder):
    [date,time]=toolbox.get_date_time_string_from_folder(folder)
    data = np.genfromtxt(folder+ '/'+time+'_monitor_warmup.dat',names=['t','V','T'])
    cdcheck=np.sum(data['T'][0:20])
    if cdcheck < 200 and size(data['t']) > 10000:
        #find first value above 5K and its place
        Tstart=next(x for x in data['T'] if x > 4.4)
        Tstart_index=np.where(data['T']==Tstart)
        time_offset=Tstart_index[0]*(data['t'][2]-data['t'][1])       
        plt.plot(data['t']-time_offset+0.5,data['T'], label='Warmup'+date)
        return toolbox.datetime_from_timestamp(date+time)
    else:
        return 
  
warmup_folders = get_warmup_folders(number_of_warmups=15)
date_array=[]


fig = plt.figure(figsize=(12,6))
for folder in warmup_folders:
    date_array.append(plt_warmup_data(folder))
legend_array=[x for x in date_array if x is not None]


plt.legend(legend_array, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
plt.xlim([0,50])
fig.suptitle('Warmups shifted in time for increased overlap', fontsize=20)
plt.xlabel('Time in Hrs')
plt.ylabel('Temperature in K')
plt.savefig(os.path.join(warmup_folders[0], 'Comparison_of_warmups.png'), format='png',pad_inches=0.1,bbox_inches='tight')

plt.show()

fig = plt.figure(figsize=(12,6))
for folder in warmup_folders:
    date_array.append(plt_warmup_data(folder))
legend_array=[x for x in date_array if x is not None]


plt.legend(legend_array, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
plt.xlim([0,10])
plt.ylim([4,25])
fig.suptitle('Warmups shifted in time for increased overlap Zoomed', fontsize=20)
plt.xlabel('Time in Hrs')
plt.ylabel('Temperature in K')
plt.savefig(os.path.join(warmup_folders[0], 'Comparison_of_warmups_zoom.png'), format='png',pad_inches=0.1,bbox_inches='tight')

plt.show()
