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
from analysis.scripts.cavity.peakdetect import peakdet
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
EOM_ON = True
EOM_freq = 6  #GHz

# fit guess parameters

# # Parameters for FSR 80 nm
# g_a1 = 0.02
# g_A1 = 0.22
# g_x01 = 0.4
# g_gamma1 = 0.2
# g_dx = 0.3
# g_A2 = 0.08
# g_gamma2 = 0.2
# g_A3 = 0.08
# g_gamma3 = 0.2


# Parameters for
g_a1 = 0.02
g_A1 = 0.02
g_x01 = 2.14
g_gamma1 = 0.1
g_A2 = 0.02
#g_x02 = 0.02
g_gamma2 = 0.1
g_dx = 0.24

 
# Open file and create dataframes in pandas

indir="D:\measuring\data/20160524\ON_diamond\Linewidth" #folder where the data is
data = pd.read_csv(os.path.join(indir,"OND_MODON000.csv"), skiprows=16, names = ["X","Y","None"],usecols=[0,1]) #creating a dataframe in pandas and importing the data

print 'found data'
print data

#Could be used for x/y ticks and labels
max_X = data['X'].max()
min_X = data['X'].min()
max_I = data['Y'].max()
min_I = data['Y'].min()

x = 1.e3*np.asarray(data['X'])#[3000:] # multiplied by factor 1.e3  to get ms  
y = np.asarray(data['Y'])#[3000:]

fixed=[]

print 'plotting data'

fig,ax = plt.subplots(figsize=(6,4))
ax.plot(x,y)
plt.show()

x=x[3000:]
y=y[3000:]

if EOM_ON == True:

    print 'fitting data to 3 lorentzians'
    p0, fitfunc, fitfunc_str = common.fit_3lorentz_symmetric(g_a1, g_A1, g_x01, g_gamma1, g_dx, g_A2, g_gamma2)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)

    x01 = fit_result['params_dict']['x01']
    dx = fit_result['params_dict']['dx']
    gamma1 = fit_result['params_dict']['gamma1']
    gamma2 = fit_result['params_dict']['gamma2']
    u_gamma1 = fit_result['error_dict']['gamma1']
    scaling = EOM_freq/dx #scale the known EOM freq with the separation here.
    lw = gamma1*scaling #scale the linewidth to get linewidht in frequency
    lw_error = u_gamma1*scaling
    print 'gamma = ', round(lw,2), '+-',round(lw_error,3),'GHz'

    #Plotting

    fig,ax = plt.subplots(figsize=(8,4))
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],len(x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)
    ax.set_xlabel("Time (ms)", fontsize = 14)
    ax.set_ylabel("Intensity (a.u.)", fontsize = 14)
    ax.set_xlim(x[0],x[-1])
    xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],n_xticks)
    #rescaling for x-axis in GHz
    X_min_freq = ax.get_xlim()[0]*scaling
    X_max_freq = ax.get_xlim()[-1]*scaling

    #xticklabels = np.linspace(X_min_freq,X_max_freq,n_xticks)

    # xticklabels_round=[]
    # for j in xticklabels:
    #   round_ = round(j,0)
    #   xticklabels_round = np.append(xticklabels_round,round_)

    # ax.set_xticks(xticks)
    # ax.set_xticklabels(xticklabels_round)

else:

    p0, fitfunc, fitfunc_str = common.fit_lorentz(g_a1, g_A1, g_x01, g_gamma1)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)
    fig,ax = plt.subplots(figsize=(8,4))
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],len(x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)
    ax.set_xlim(x[0],x[-1])
    ax.set_xlabel("Time (ms)", fontsize = 14)
    ax.set_ylabel("Intensity (a.u.)", fontsize = 14)


## you can use this for setting axes ###


# freq_range = ((max_X1-min_X1)*EOM_freq)/h
# print (max_X1-min_X1), EOM_freq, h

# xticks = np.linspace(min_X1,max_X1,n_xticks)
# xticklabels = np.linspace(min_X1,max_X1,n_xticks)

# yticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],n_yticks)
# ytickslabels = np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],n_yticks)

# xticklabels_round=[]
# for j in xticklabels:
#   round_ = round(j,0)
#   xticklabels_round = np.append(xticklabels_round,round_)

# ytickslabels_round=[]
# for i in ytickslabels:
#   round_ = round(i,0)
#   ytickslabels_round = np.append(ytickslabels_round,round_)

# ax.set_yticks(yticks)
# ax.set_yticklabels(ytickslabels_round)

plt.show()






