# This script is written in order to analyse the data obtained by the cavity.
# In the first part of the script the information about the to-be-analysed data should be entered
# This script makes use of some formulas, embedded in the script fitting_formula
# 2015



import os
import h5py
import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot; reload(plot)
import pylab as plt
from analysis.lib import fitting
from analysis.lib.fitting import esr, common
from analysis.lib.fitting import fit,esr, common
import matplotlib
from matplotlib import rc, cm
from scipy import interpolate
from scipy.optimize import curve_fit
import operator
import sys


import analysis.scripts.cavity.analysis_linewidth_fsr as fitting; reload(fitting)
#get_FSR, get_piezo_linewidth, get_laser_linewidth, fit_laser_linewidth, show_piezo_plots, show_laser_plots, show_fine_laser_plots, fit_fine_laser_plots, fit_fine_laser_plots2, bla_fine_laser_plots, fit_piezo_plots
reload(toolbox)


'''
This is the only part where values should be entered (if the script works as intended). The folder with the data can be 
selected, the circumstances from the measurement can be appended in 'measurement_type', so that it will be added to plot titles.
Furthermore for the fitting method either 'gauss', or 'lorentz' can be selected. It might be nice to see the graphs first before
fitting them, then fitting should be False (or zero). In order to start the fitting, fitting should be True (or 1).
By means of looking at the plots, the first and last timestamp should be provided of the desired range. Finally the min and
max value for the Voltage can be chosen in order to make the plots clearer.
'''


''' Enter your data '''
date = '20160126'
save_data = True

''' Select the time stamps you want to include or remove'''
time_start = 123331    #161322#172207#231234
time_stop = 123334#161341#145900#150629#172727#231357
remove = []#['225951','225948']

#pick a threshold for the peak detection
threshold = 0.01#0.008#0.015


''' Switch for if you want to set the x-limits of the piezo-scans. 
	NB: these are very sensitive, if you take it too narrow, it will not fit anymore '''
set_range_V = False
V_min = 2.13
V_max = 2.17

''' Switch for if you want to set the x-limits of the laser scans.
	NB: these are very sensitive, if you take it too narrow, it will not fit anymore '''
set_range_f = False
f_min = -30
f_max = 30


###########################################################################################################################################################
###########################################################################################################################################################

'''
Select all the folders in the preferred folder and sort them. Thereafter, a directory is created where the plots will be saved.
For all the timestamps in the folder, only the timestamp will be selected from the title. And thereafter only the preferred
timestamps are taken into account.
'''

folder = 'D:/measuring/data/20160126/'
all_folders = [f for f in os.listdir(folder) if f[0] in '1233314740']
all_folders.sort()

all_timestamp = [item[0:6] for item in all_folders]
timestamp = [all_timestamp[i] for i in range(len(all_timestamp)) if int(all_timestamp[i]) >= time_start and int(all_timestamp[i]) <= time_stop and all_timestamp[i] not in remove]


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################



#fitting.get_piezo_linewidth(date = date, timestamp = timestamp, folder = folder, measurement_type = measurement_type, V_min = V_min, V_max = V_max, set_range_V = set_range_V, save_data = save_data, newpath = newpath)
#fitting.get_laser_linewidth(date = date, timestamp = timestamp, folder = folder, measurement_type = measurement_type, threshold = threshold, f_max = f_max, f_min = f_min, set_range_f = set_range_f, save_data = save_data, newpath = newpath)
#fitting.fit_laser_linewidth(date = date, timestamp = timestamp, folder = folder, measurement_type = measurement_type, threshold = threshold, f_max = f_max, f_min = f_min, set_range_f = set_range_f, newpath = newpath)

# fitting.get_FSR(date=date, timestamp=timestamp, folder=folder, newpath=newpath, V_min=V_min, V_max=V_max, set_range_V=set_range_V, labda = 637, save_data = False)
#fitting.show_piezo_plots(date = date, timestamp = timestamp, folder = folder, measurement_type = measurement_type, V_min = V_min, V_max = V_max, set_range_V = set_range_V)
#fitting.show_laser_plots(date = date, timestamp = timestamp, folder = folder, measurement_type = measurement_type, f_min = f_min, f_max = f_max, set_range_f = set_range_f)

#fitting.show_fine_laser_plots(date = date, timestamp = timestamp, folder = folder, measurement_type = measurement_type, f_min = f_min, f_max = f_max, set_range_f = set_range_f)

#fitting.fit_fine_laser_plots2(date = date, timestamp = timestamp, folder = folder, measurement_type = measurement_type, f_min = f_min, f_max = f_max, set_range_f = False, save_data = False, threshold = threshold)

# fitting.bla_fine_laser_plots(date = date, timestamp = timestamp, folder = folder, measurement_type = measurement_type, f_min = f_min, f_max = f_max, set_range_f = False, save_data = False)




#fitting.fit_fine_laser_plots(date = date, timestamp = timestamp, folder = folder,  f_min = f_min, f_max = f_max, set_range_f = set_range_f, save_data = save_data, threshold=threshold)
fitting.fit_piezo_plots(date = date, timestamp = timestamp, folder = folder, V_min = V_min, V_max = V_max, 
    set_range_V = set_range_V, save_data = save_data, threshold = threshold,show_plots=False)

# fitting.fit_birefringence(date = date, timestamp = timestamp, folder = folder, V_min = V_min, 
#     V_max = V_max, set_range_V = set_range_V, save_data = save_data, threshold = threshold,show_plots=True)
