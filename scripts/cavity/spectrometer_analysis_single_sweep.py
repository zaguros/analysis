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

# values to input 

rows = 20 #number of rows in spectrometer 
refr= 1 #refractive index 

# Open file and create dataframes in pandas

indir="D:\measuring\data/20160418" 
data = pd.read_csv(os.path.join(indir,"Z_884_FP_8V 2016 april 18 16_07_51.csv"), usecols=[2,4,5])
data_mean=data.groupby('Column').agg([np.mean])


# # Peak detection method for higher order modes 

# I_array = np.asarray(data_mean['Intensity'])
# peak_1=peakdet(I_array, 1000)
# peak_WL_1=np.transpose(np.asarray(data_mean.loc[peak_1[0][:,0],'Wavelength']))[0]
# peak_I_1=peak_1[0][:,1]
# peak_freq_1 = 3.e8/(peak_WL_1*1.e-9)

# Peak detection method for FSR

I_avg_array = np.asarray(data_mean['Intensity']) # argrelextrema only takes an array
indices = argrelextrema(I_avg_array, np.greater, order=100) # the number 100 is somewhat arbitrary, but seems to work. 

peak_WL= [] #connect indices with the values for the wavelength, creating arrays
peak_I=[]
peak_freq=[]

for i in indices[0]:
    peak_WL = np.append(peak_WL,data_mean.loc[i,'Wavelength'])
    peak_I = np.append(peak_I,data_mean.loc[i,'Intensity'])
    peak_f = 3.e8/(data_mean.loc[i,'Wavelength']*1.e-9)
    peak_freq = np.append(peak_freq,peak_f)

# Calculating the free spectral range in wavelength

FSR=fabs(peak_WL[-2]-peak_WL[-1]) #calculate the free spectral range
print 'The FSR is', FSR,'nm.'

# Calculating the free spectral range in frequency in Hz

FSR_freq=fabs(peak_freq[-2]-peak_freq[-1])

print 'The average FSR in frequency is', round(FSR_freq*1.e-12,2), 'THz.'

# Calculating the length of the cavity in micrometer

L = ((peak_WL[-3]*1.e-9)**2)/(refr*2*FSR*1.e-9) #in meter
print 'The Cavity Length is', round(L*1.e6,2), 'um.'

L2 = 3.e8/(2*refr*FSR_freq)
print 'The Cavity Length is', round(L2*1.e6,2), 'um.'

# The optical length of the cavity 

l = 2 * L 
print 'The optical length of the cavity is', round(l*1.e6,2), 'um.'
#Calculating the Finesse of the cavity with an assumption for the linewidth

dv = 1.3*1.e9 
F = FSR_freq/dv

print 'The Finesse is', round(F,0), ',under the assumption that the linewidth is', round(dv*1.e-9,2), 'GHz.'

# #Calculating the effective radius of curvature

# df_trans=fabs(peak_freq_1[-1]-peak_freq_1[-2]) # value from data: 2 nm (not conclusive yet!)
# T=(df_trans/FSR_freq)*pi
# ROC =l*(1/(1-(cos(T)**2)))
# print 'The ROC is', round(ROC*1.e6,2),'um.' 

# #Calculating the beam waist and mode volume 

# w_0= sqrt((peak_WL[-2]*1.e-9)/pi)*(l*(ROC-l))**(1/4)
# V_0 = (pi*(w_0**2)*l)/4

# print 'The beam waist is', round(w_0*1.e6,2)
# print 'The mode volume is', round(V_0*1.e12,2)

#Plotting the data  

data_mean.plot(x='Wavelength',y='Intensity',title='FSR bare cavity', legend=False) #ylim=[0,14000], xticks=range(600,700,10), yticks=range(0,14000,2000))
plt.savefig("test.png", format="png")
plt.plot(peak_WL, peak_I, 'ro') # for assigning the peaks
#plt.plot(peak_WL_1, peak_I_1, 'r+') # for assigning the higher order peaks
plt.show()







