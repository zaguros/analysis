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
This script is written to import data from the oscilloscope in the CSV format to Python and characterize the cavity. 
"""

### This is the only part where values should be inserted.#### 
n_xticks = 1.e-7
n_yticks = 8
peak_detect_TEM00 = False
order_peak_detection = 300
 
# Open file and create dataframes in pandas

indir="D:\DATA/20160419 Linewidth" #folder where the data is
data = pd.read_csv(os.path.join(indir,"FRI013_FSR_80nm_RFON_20DB.csv")) #creating a dataframe in pandas and importing the data

#Could be used for x/y ticks and labels
max_X = data['X'].max()
min_X = data['X'].min()
max_I = data['Y'].max()
min_I = data['Y'].min()

# Peak detection method for TEM00 modes doesn't work so wel for this though

if peak_detect_TEM00 == True:

	I_avg_array = np.asarray(data['Y']) # argrelextrema only takes an array
	indices = argrelextrema(I_avg_array, np.greater, order=order_peak_detection) # the number 100 is somewhat arbitrary, but seems to work. 

	peak_X= [] #connect indices with the values for the wavelength, creating arrays
	peak_I=[]
	#peak_freq=[]

	for i in indices[0]:
		peak_X = np.append(peak_X,data.loc[i,'X'])
		peak_I = np.append(peak_I,data.loc[i,'Y'])
		#peak_f = 3.e8/(data_mean.loc[i,'X']*1.e-9)
		#peak_freq = np.append(peak_freq,peak_f)

	# FSR=fabs(peak_WL[-2]-peak_WL[-1]) #calculate the free spectral range
	# print 'The FSR is', FSR,'nm.'

	# FSR_freq=fabs(peak_freq[-2]-peak_freq[-1]) # calculating the free spectral range in frequency in Hz
	# print 'The FSR is', round(FSR_freq*1.e-12,2), 'THz.'

	# L = 3.e8/(2*FSR_freq) # Calculating the length of the cavity in micrometer
	# print 'The Cavity Length is', round(L*1.e6,2), 'um.'

else:
	print 'No TEM00 peaks are detected for this data, so no FSR and cavity length is determined.'
	

#Plotting the data and setting the axes 

ax = data.Y.plot(color = 'r') #super simple plot of the data 

if peak_detect_TEM00 == True:
	plt.plot(indices[0], peak_I, 'bo') # for assigning the peaks
else:
	print 'No TEM00 peaks assigned in plot'

ax.set_xlabel("Time (us)", fontsize = 14)
ax.set_ylabel("Intensity (a.u.)", fontsize = 14)

### you can use this for setting axes ###

#xticks = np.linspace(min_X,max_X,n_xticks)
#xticklabels = np.linspace(min_X,max_X,n_xticks)

#yticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],n_yticks)
#ytickslabels = np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],n_yticks)

# xticklabels_round=[]
# for j in xticklabels:
# 	round_ = round(j,0)
# 	xticklabels_round = np.append(xticklabels_round,round_)

#ytickslabels_round=[]
#for i in ytickslabels:
#	round_ = round(i,0)
#	ytickslabels_round = np.append(ytickslabels_round,round_)

#ax.set_xticks(xticks)
#ax.set_xticklabels(xticklabels_round)

#ax.set_yticks(yticks)
#ax.set_yticklabels(ytickslabels_round)

plt.show()











