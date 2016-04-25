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
from analysis.lib import fitting
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot

reload(common)


"""
This script is written to import data from the oscilloscope in the CSV format to Python and characterize the cavity. 
"""

### This is the only part where values should be inserted.#### 
n_xticks = 5
n_yticks = 8
peak_detect_TEM00 = False
order_peak_detection = 700
EOM_ON = True
EOM_freq = 6  #GHz

# fit guess parameters

g_a1 = 0.01
g_A1 = 0.05
g_x01 = 0
g_gamma1 = 0.15
g_A2 = 0.09
g_x02 = 0.7
g_gamma2 = 0.2
g_A3 = 0.3
g_x03 = 0.4
g_gamma3 = 0.15

 
# Open file and create dataframes in pandas

indir="D:\measuring\data/20160419" #folder where the data is
data = pd.read_csv(os.path.join(indir,"FRI020_FSR_50nm_RFON_20DB.csv"), skiprows=16, names = ["X","Y","Z"],usecols=[0,1]) #creating a dataframe in pandas and importing the data

#Could be used for x/y ticks and labels
max_X = data['X'].max()
min_X = data['X'].min()
max_I = data['Y'].max()
min_I = data['Y'].min()

x = 1.e5*np.asarray(data['X'])
y = np.asarray(data['Y'])

fixed=[]

if EOM_ON == True:

	p0, fitfunc, fitfunc_str = common.fit_3lorentz(g_a1, g_A1, g_x01, g_gamma1, g_A2, g_x02, g_gamma2, g_A3, g_x03, g_gamma3)
	fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)

	x01 = fit_result['params_dict']['x01']
	x02 = fit_result['params_dict']['x02']
	x03 = fit_result['params_dict']['x03']
	gamma1 = fit_result['params_dict']['gamma1']
	gamma2 = fit_result['params_dict']['gamma2']
	gamma3 = fit_result['params_dict']['gamma3']

	h = x02 - x03
	j = x03 - x01
	scaling = EOM_freq/h
	lw = gamma2*scaling


#	Calculating


	print 'The cavity linewidth is:', round(lw,2), 'GHz'


	#Plotting


	fig,ax = plt.subplots(figsize=(8,4))


	plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],len(x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)
	ax.set_xlabel("Frequency (GHz)", fontsize = 14)
	ax.set_ylabel("Intensity (a.u.)", fontsize = 14)
	# ax.set_xlim(x[0],x[-1])
	# xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],n_xticks)
	# #rescaling for x-axis in GHz
	# X_min_freq = ax.get_xlim()[0]*scaling
	# X_max_freq = ax.get_xlim()[-1]*scaling

	# xticklabels = np.linspace(X_min_freq,X_max_freq,n_xticks)

	# xticklabels_round=[]
	# for j in xticklabels:
	# 	round_ = round(j,0)
	# 	xticklabels_round = np.append(xticklabels_round,round_)

	# ax.set_xticks(xticks)
	# ax.set_xticklabels(xticklabels_round)

else:


	p0, fitfunc, fitfunc_str = common.fit_lorentz(g_a1, g_A1, g_x01, g_gamma1)
	fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)
	fig,ax = plt.subplots(figsize=(8,4))
	plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],len(x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)
	ax.set_xlim(x[0],x[-1])
	ax.set_xlabel("Time (a.u.)", fontsize = 14)
	ax.set_ylabel("Intensity (a.u.)", fontsize = 14)



#Plotting the data and setting the axes 


## you can use this for setting axes ###


# freq_range = ((max_X1-min_X1)*EOM_freq)/h
# print (max_X1-min_X1), EOM_freq, h

# xticks = np.linspace(min_X1,max_X1,n_xticks)
# xticklabels = np.linspace(min_X1,max_X1,n_xticks)

# yticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],n_yticks)
# ytickslabels = np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],n_yticks)

# xticklabels_round=[]
# for j in xticklabels:
# 	round_ = round(j,0)
# 	xticklabels_round = np.append(xticklabels_round,round_)

# ytickslabels_round=[]
# for i in ytickslabels:
# 	round_ = round(i,0)
# 	ytickslabels_round = np.append(ytickslabels_round,round_)



# ax.set_yticks(yticks)
# ax.set_yticklabels(ytickslabels_round)





plt.show()






