import sys, os,time
import h5py
import numpy as np
import pylab as plt
'''
File that is made to test importing and plotting hdf5 data
'''
#Import the data
h5filepath='/Users/Adriaan/Documents/teamdiamond/data_for_analysis/20140312/172721_ElectronRabi_Hans_sil1_Rabi-1/172721_ElectronRabi_Hans_sil1_Rabi-1.hdf5'
f = h5py.File(h5filepath,'r')
name = f.keys()[0]
g = f[name]
adwingrpname = g.keys()[1]
adwingrp = g[adwingrpname]

reps = adwingrp['completed_reps'].value
sweep_pts = adwingrp.attrs['sweep_pts'] #in ns
RO_data = adwingrp['RO_data'].value
normalized_RO_data = RO_data/(float(reps/len(sweep_pts)))


plt.plot(sweep_pts,normalized_RO_data, 'ro')
plt.xlabel('time ns')
plt.ylabel('F(|0>)')
plt.title(r'')
plt.grid(True)

plt.show()
