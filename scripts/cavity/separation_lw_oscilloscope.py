import numpy as np 
import scipy.constants
import scipy.interpolate
from matplotlib import pyplot as plt
from math import pi, sqrt, sin, cos, fabs
import csv
from itertools import *
from scipy import signal 
import itertools
from scipy.signal import argrelextrema
import pandas as pd
import os
from peakdetect import peakdet
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.stats import threshold

order_peak_detection = 700

indir="D:\measuring\data/20160426\separation_lw" 
data = pd.read_csv(os.path.join(indir,"LONG000.csv"), skiprows=16, names = ["time","ramp","res","nothing"],usecols=[0,2]) #creating a dataframe in pandas and importing the data

res_array = np.asarray(data['res']) # argrelextrema only takes an array
indices = argrelextrema(res_array, np.greater, order=order_peak_detection) # the number 100 is somewhat arbitrary, but seems to work. 

peak_time= [] #connect indices with the values for the wavelength, creating arrays
peak_res=[]

for i in indices[0]:
    peak_time = np.append(peak_time,data.loc[i,'time'])
    peak_res = np.append(peak_res,data.loc[i,'res'])

sep = []
for i in xrange(0,len(peak_time)-2,2):
    jump = fabs(peak_time[i]-peak_time[i+2])
    sep = np.append(sep,jump)

th_sep = threshold(sep,0.0015)
th_list=th_sep.tolist()
th_list[:] = (value for value in th_list if value != 0.0)
#th_[th_sep != 0]

# sep1=[]
# for i in sep:
#     if i >= 0.0015:
#         sep1 = np.append(sep1,sep)

# print sep1

# best fit of data
(mu, sigma) = norm.fit(th_list)
bins = np.linspace(0.0005,0.0025,100)


plt.hist(th_list, bins=bins,facecolor='green')
#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--', linewidth=2)

plt.show()
#make histogram of sep


