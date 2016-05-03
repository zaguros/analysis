#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from stat import S_ISREG, ST_MTIME, ST_MODE
import os, sys, time
import matplotlib.image as mpimg
from scipy.signal import argrelextrema


indir="C:\Users\suzannevandam\Documents\localdata\spectrometer_data/20160330/CSV_2Dplot_diamond_sweep_l" 
outdir="H:\My Documents\Cavities\SpectrumData"

# parameters to vary per measurement Note: you might have to change the vmin and vmax of the colorbar inside the script! 
V_min = 8
V_max = 2
n_xticks= 7 #how many ticks you want on the x-axis
n_yticks = 11 #how many ticks you want on the y-axis
peak_detect_TEM00 = False
order_peak_detection = 100


# load the files in the dataframe, only pick out Intensity and Column. All files in folder. 

files = (os.path.join(indir, fn) for fn in os.listdir(indir))
files = ((os.stat(path), path) for path in files)
files = ((stat[ST_MTIME], path)
           for stat, path in files if S_ISREG(stat[ST_MODE]))

# for cdate, path in sorted(files)[:-1]:
#     print time.ctime(cdate), os.path.basename(path)

dataframes = [pd.read_csv(path, usecols = [2,5]) for cdate,path in sorted(files)]

# dataframes=np.array([])
# print dataframes
# for cdate, path in sorted(files):
#   dataframes=np.append(dataframes,pd.read_csv(path, usecols=[2,5]))
#   print pd.read_csv(path, usecols=[2,5])
#   print path
# print dataframes

#load the file in the dataframe with a certain name in a certain folder
#number of files
#dataframes = [pd.read_csv(os.path.join(indir,"160311_FP_%s_Par1.csv") % i, usecols=[2,5]) for i in xrange(0,11)]

# Calculating maximum and minimum values for wavelength and intensity. Intensity we have to use a for loop because we want the overall maximum from all files. 

max_WL = dataframes[0]['Wavelength'].max()
min_WL = dataframes[0]['Wavelength'].min()

max_I_all=[]
for i in dataframes:
    max_I_single = i['Intensity'].max()
    max_I_all = np.append(max_I_all,max_I_single)

max_I = np.amax(max_I_all)

min_I_all=[]
for i in dataframes:
    min_I_single = i['Intensity'].min()
    min_I_all = np.append(min_I_all,min_I_single)

min_I = np.amin(min_I_all)

# Group the Intensity data by column and compute the mean
mean_val1=[]
for i in dataframes:
	merged=i.groupby('Wavelength')
	mean=merged.agg([np.mean])
	mean_val1.append(mean)

if peak_detect_TEM00 == True:
	
    for data in mean_val:
        I_avg_array = np.asarray(data['Intensity']) # argrelextrema only takes an array
        indices = argrelextrema(I_avg_array, np.greater, order=order_peak_detection) # the number 100 is somewhat arbitrary, but seems to work. 
				
        peak_WL= [] #connect indices with the values for the wavelength, creating arrays
        peak_I=[]
        peak_freq=[]

		for i in indices:
			peak_WL = np.append(peak_WL,data.ix['Wavelength'])
			peak_I = np.append(peak_I,data.loc[i,'Intensity'])
			peak_f = 3.e8/(data.loc[i,'Wavelength']*1.e-9)
			peak_freq = np.append(peak_freq,peak_f)

		FSR = []
		FSR=fabs(peak_WL[-2]-peak_WL[-1]) #calculate the free spectral range
		
		print 'The FSR is', FSR,'nm.'

        FSR_freq=fabs(peak_freq[-2]-peak_freq[-1]) # calculating the free spectral range in frequency in Hz
        print 'The FSR is', round(FSR_freq*1.e-12,2), 'THz.'

        L = 3.e8/(2*FSR_freq) # Calculating the length of the cavity in micrometer
        print 'The Cavity Length is', round(L*1.e6,2), 'um.'

        range_WL = linspace(max_WL,min_WL,n_yticks)
        plt.plot(range_WL,L)
else:
    print 'No TEM00 peaks are detected for this data, so no FSR and cavity length is determined.'


# Concatenate in one dataframe
mean_val=pd.concat(mean_val1,axis=1)


# Plot using seaborn function and set all the axes
ax=sns.heatmap(mean_val, vmin = 0, vmax= 5000, cmap='YlGnBu')


ax.set_xlabel("Voltage (V)", fontsize = 14)
ax.set_ylabel("Wavelength (nm)", fontsize = 14)

ax.set_xlabel("Voltage (V)", fontsize = 20)
ax.set_ylabel("Wavelength (nm)", fontsize = 20)

ax.tick_params(which = 'both', direction = 'out')
xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],n_xticks)
xticklabels = np.linspace(V_min,V_max,n_xticks)

yticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],n_yticks)
ytickslabels = np.linspace(max_WL,min_WL,n_yticks)


xticklabels_round = np.round(xticklabels,1)

ytickslabels_round=np.round(ytickslabels,0).astype(int)

ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels_round, rotation=0, fontsize=20)

ax.set_yticks(yticks)
ax.set_yticklabels(ytickslabels_round, fontsize=20)

# plt.show()
try:
    plt.savefig(os.path.join(outdir,'2Dplot_diamond_cavity.jpg'))
    plt.savefig(os.path.join(outdir,'2Dplot_diamond_cavity.eps'))
except:
    print('could not save figure')
#mpimg.imsave("out.png", fig)
#plt.savefig(os.path.join(outdir, "2Dplot.eps"), format="png")


# # plot the maximum value of the intensity for the single piezo sweeps

# xaxis = np.linspace(V_min,V_max,21)
# plt.plot(xaxis, max_I_all,'ro', )
# plt.xlabel("Voltage (V)", fontsize = 14, fontweight='bold')
# plt.ylabel("Intensity (a.u.)", fontsize = 14, fontweight = 'bold')
# plt.xlim(V_min-0.1, V_max+0.1)
# plt.show()


















