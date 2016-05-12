#################################
# L.C. Coenen
# coenen.lisanne@gmail.com
#################################
# Data analysis script of cavity characteristics and modes


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

"""
This script is written to import data from the spectrometer in the CSV format to Python and characterize the cavity. 
"""

### This is the only part where values should be inserted.#### 
n_xticks = 9
n_yticks = 8
peak_detect_H_order_modes = False
peak_detect_TEM00 = True
order_peak_detection = 100
order_peak_detect_higher_modes = 30
 
# Open file and create dataframes in pandas

indir="D:\measuring\data/20160430\ON_diamond" 
data = pd.read_csv(os.path.join(indir,"2016-04-30 19_34_17 cal_finepiezo.csv"), usecols=[2,4,5])

max_WL = data['Wavelength'].max()
min_WL = data['Wavelength'].min()
max_I = data['Intensity'].max()
min_I = data['Intensity'].min()
data_mean=data.groupby('Column').agg([np.mean])

# Peak detection method for higher order modes 

if peak_detect_H_order_modes == True:

	I_array = np.asarray(data_mean['Intensity']) # write intensity dataframe as array
	peak_1=peakdet(I_array, order_peak_detect_higher_modes) # use peakdet definition to detect all peaks from intensity
	peak_WL_1=np.transpose(np.asarray(data_mean.loc[peak_1[0][:,0],'Wavelength']))[0] #get corresponding wavelengths
	peak_I_1=peak_1[0][:,1]
	peak_freq_1 = 3.e8/(peak_WL_1*1.e-9)
else:
	print 'There are no higher order modes in this plot.'

# Peak detection method for TEM00 modes

if peak_detect_TEM00 == True:

	I_avg_array = np.asarray(data_mean['Intensity']) # argrelextrema only takes an array
	indices = argrelextrema(I_avg_array, np.greater, order=order_peak_detection) # the number 100 is somewhat arbitrary, but seems to work. 

	peak_WL= [] #connect indices with the values for the wavelength, creating arrays
	peak_I=[]
	peak_freq=[]

	for i in indices[0]:
		peak_WL = np.append(peak_WL,data_mean.loc[i,'Wavelength'])
		peak_I = np.append(peak_I,data_mean.loc[i,'Intensity'])
		peak_f = 3.e8/(data_mean.loc[i,'Wavelength']*1.e-9)
		peak_freq = np.append(peak_freq,peak_f)

	FSR=fabs(peak_WL[-2]-peak_WL[-1]) #calculate the free spectral range
	print 'The FSR is', FSR,'nm.'

	FSR_freq=fabs(peak_freq[-2]-peak_freq[-1]) # calculating the free spectral range in frequency in Hz
	print 'The FSR is', round(FSR_freq*1.e-12,2), 'THz.'

	peak_1_WL = peak_WL[0]
	print 'The wavelength of the 1st peak is', peak_1_WL
	peak_2_WL = peak_WL[1]
	print 'The wavelength of the 2nd peak is', peak_2_WL

	L = 3.e8/(2*FSR_freq) # Calculating the length of the cavity in micrometer
	print 'The Cavity Length is', round(L*1.e6,2), 'um.'

else:
	print 'No TEM00 peaks are detected for this data, so no FSR and cavity length is determined.'
	
#Plotting the data and setting the axes 

ax = data_mean.plot(x='Wavelength',y='Intensity', legend=False) 

if peak_detect_TEM00 == True:
	plt.plot(peak_WL, peak_I, 'ro') # for assigning the peaks
else:
	print 'No TEM00 peaks assigned in plot'
if peak_detect_H_order_modes == True:
	plt.plot(peak_WL_1, peak_I_1, 'r+') # for assigning the higher order peaks
else:
	print 'No peaks for higher order modes assigned in plot.'

ax.set_xlabel("Wavelength (nm)", fontsize = 14)
ax.set_ylabel("Intensity (a.u.)", fontsize = 14)
ax.tick_params(which = 'both', direction = 'out')
#ax.set_aspect('equal')
xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],n_xticks)
xticklabels = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],n_xticks)

#yticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],n_yticks)
#ytickslabels = np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],n_yticks)

xticklabels_round=[]
for j in xticklabels:
	round_ = round(j,0)
	xticklabels_round = np.append(xticklabels_round,round_)

#ytickslabels_round=[]
#for i in ytickslabels:
#	round_ = round(i,0)
#	ytickslabels_round = np.append(ytickslabels_round,round_)

ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels_round)

#ax.set_yticks(yticks)
#ax.set_yticklabels(ytickslabels_round)

plt.show()






