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

order_peak_detection = 500
peak_thresh = 0.1

indir="D:\measuring\data/20160426\separation_lw" 
data = pd.read_csv(os.path.join(indir,"1s_scan_tubeon.csv"), skiprows=16, names = ["time","ramp","res","nothing"],usecols=[0,2]) #creating a dataframe in pandas and importing the data

xmin_index = 0
xmax_index = 100000
res_array = np.asarray(data['res'])[xmin_index:xmax_index] # argrelextrema only takes an array
time_array = np.asarray(data['time'])[xmin_index:xmax_index]

fig,ax = plt.subplots(figsize=(15,4))

indices = argrelextrema(res_array, np.greater, order=order_peak_detection) # the number 100 is somewhat arbitrary, but seems to work. 
# peaks, min_peaks = peakdet(res_array[0:40000],2)

# indices_2 = peaks[:,0].astype(int)
# print indices_2
peak_time= [] #connect indices with the values for the wavelength, creating arrays
peak_res=[]
# peak_time_2= [] #connect indices with the values for the wavelength, creating arrays
# peak_res_2=[]

for i in indices[0]:
    if data.loc[i,'res']<peak_thresh:
        pass
    else:
        peak_time = np.append(peak_time,data.loc[i,'time'])
        peak_res = np.append(peak_res,data.loc[i,'res'])

print len(peak_time)
# for i in indices_2:
#     peak_time_2 = np.append(peak_time_2,data.loc[i,'time'])
#     peak_res_2 = np.append(peak_res_2,data.loc[i,'res'])


ax.plot(time_array,res_array,'o')
# print peak_time
# print peak_time_2
ax.plot(peak_time,peak_res,'ro')
# ax.plot(peak_time_2,peak_res_2,'co')

ax.set_xlim((time_array[0],time_array[-1]))
plt.show()

sep = []
for i in xrange(0,len(peak_time)-2,2):
    jump = fabs(peak_time[i]-peak_time[i+2])
    sep = np.append(sep,jump)


th_sep = threshold(sep,0.0)
th_list=th_sep.tolist()
th_list[:] = (value for value in th_list if value != 0.0)
th_list = np.multiply(th_list,1000) 

#th_[th_sep != 0]

# sep1=[]
# for i in sep:
#     if i >= 0.0015:
#         sep1 = np.append(sep1,sep)

# print sep1

# best fit of data
# (mu, sigma) = norm.fit(th_list)
# print(mu, sigma)
bins = np.linspace(0,10,150)

def f(x, a, b, c):
    return a * py.exp(-(x - b)**2.0 / (2 * c**2))

# Generate data from bins as a set of points 
x = [0.5 * (th_list[1][i] + th_list[1][i+1]) for i in xrange(len(th_list[1])-1)]
y = th_list[0]

popt, pcov = optimize.curve_fit(f, x, y)

x_fit = py.linspace(x[0], x[-1], 100)
y_fit = f(x_fit, *popt)

plot(x_fit, y_fit, lw=4, color="r")

fig2,ax2 = plt.subplots(figsize=(10,4))


ax2.hist(th_list, bins=bins,facecolor='green')
# y = mlab.normpdf( bins, mu, sigma)
# l = plt.plot(bins, y, 'r--', linewidth=)
ax2.set_xlabel('dt (ms)')
plt.show()
#make histogram of sep


