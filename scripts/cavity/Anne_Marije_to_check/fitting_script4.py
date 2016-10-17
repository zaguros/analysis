import os
import h5py
import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
import pylab as plt
from analysis.lib import fitting
from analysis.lib.fitting import esr, common
from analysis.lib.fitting import fit,esr, common
from matplotlib import rc, cm
from scipy import interpolate
from scipy.optimize import curve_fit
import operator



reload(toolbox)


'''
These two functions do the following:
the first function, written by Cristian, picks out all the folders in the folder of the day that you indicated with the preferred 
timestamp, which should be just one folder. Thereafter it picks the 'hdf5' file from this folder (i.e. the file that contains the 
data) and selects the x and y data (which is different for laser scans and piezo scans).

The second function loads a piezo scan specifically and specifies x and y.
'''

def load_data (folder, timestamp, scan_type):

	all_folders = [f for f in os.listdir(folder) if timestamp in f]
	#print all_folders
	curr_fold = os.path.join(folder, all_folders[0])
	all_files =[ f for f in os.listdir(curr_fold) if os.path.isfile(os.path.join(curr_fold,f)) ]
	file_name = [f for f in all_files if '.hdf5' in f]

	if (scan_type =='piezo'):
		grp_name = '/piezo_scan'
		x_name = 'piezo_voltage'
	elif (scan_type == 'lr_scan'):
		grp_name = '/lr_scan'
		x_name = 'frequency_GHz'
	
	f = h5py.File(os.path.join(curr_fold, file_name[0]),'r')
	ls_grp = f[grp_name]
	x = ls_grp[x_name].value
	y = ls_grp['PD_signal'].value
	return x, y


def load_piezo_scan (folder, timestamp):
	V, y = load_data (folder = folder, timestamp = timestamp, scan_type = 'piezo')
	V = np.ravel(V)
	y = np.ravel(y)

	# ind = np.where(abs(V)>5000)
	# V = np.delete(V, ind)
	# y = np.delete(y, ind)
	return V, y

def load_lr_scan (folder, timestamp):
	f, y = load_data (folder = folder, timestamp = timestamp, scan_type = 'lr_scan')
	f = np.ravel(f)
	y = np.ravel(y)

	ind = np.where(abs(f)>1000)
	f = np.delete(f, ind)
	y = np.delete(y, ind)
	return f, y

##########################################################################################################################################################
###########################################################################################################################################################
##########################################################################################################################################################
###########################################################################################################################################################
'''
This is the only part where values should be enetered (if the script works as intended). The folder with the data can be 
selected, the circumstances from the measurement can be appended in 'measurement_type', so that it will be added to plot titles.
Furthermore for the fitting method either 'gauss', or 'lorentz' can be selected. It might be nice to see the graphs first before
fitting them, then fitting should be False (or zero). In order to start the fitting, fitting should be True (or 1).
By means of looking at the plots, the first and last timestamp should be provided of the desired range. Finally the min and
max value for the Voltage can be chosen in order to make the plots clearer.
'''


# folder = 'H:/My Documents/TU/Afstuderen/Data/20150630/RT_tubeoff_Nin/piezo_scans'

''' Enter your data '''
date = '20150728'
measurement_type = 'LT_tubelow_newL(2u)_syncoff'
scan_type = ''
# 'piezo_scans'

''' Fitting yes/no '''
fitting = 1
fitting_method = 'gauss'
fixed = [0]

''' Select the time stamps you want to include or remove'''
time_start = 100000
time_stop = 235900
# remove = ['201245','201351','201507','201949','202023','202450', '203635', '203718', '203632', '203630','203634', '203644', '203659', '212152', '212206']
# remove = ['185132','192231','191401','201010','201015','201022','201035','201040','201245','201351','201428','201507','201949','202023','202450', '203635', '203644', '203718', '204314', '212152', '212157', '212206', '212543', '212618', '212652', '212618', '212229', '212543', '212322', '213108','184254', '213513','213626', '183919','184335','185004', '185059']
# remove = ['145047','145540','150244','150905','151009','151112', '145540', '145646', '151129', '150944', '151246', '145244','150347', '150122', '120429', '212322', '120517', '121146','122150', '122259', '121059']
# remove = ['175318','164643', '164647', '164657', '164659', '164703', '164715', '164717', '164718', '164720','155755', '160818', '161605', '162219','173210', '173212','173214', '173215', '173841', '154408', '154413', '154420', '154430', '154440', '154723']
remove = []


''' Switch for if you want to set the x-limits of the piezo-scans. 
	NB: these are very sensitive, if you take it too narrow, it will not fit anymore '''
set_range_V = False
V_min = 6.4
V_max = 6.42

''' Switch for if you want to set the x-limits of the laser scans.
	NB: these are very sensitive, if you take it too narrow, it will not fit anymore '''
set_range_f = 1
f_min = -600
f_max = -400

threshold = 1


###########################################################################################################################################################
###########################################################################################################################################################

'''
Select all the folders in the preferred folder and sort them. Thereafter, a directory is created where the plots will be saved.
For all the timestamps in the folder, only the timestamp will be selected from the title. And thereafter only the preferred
timestamps are taken into account.
'''

folder = 'K:/ns/qt/Diamond/Projects/Cavities/data/' + date + '/' + measurement_type
 # + '/' + scan_type
all_folders = [f for f in os.listdir(folder)]
all_folders.sort()

newpath =  'K:/ns/qt/Diamond/Projects/Cavities/data_fitting/vibration measurements/' + date + '/' + measurement_type
 # + '/' + scan_type
if not os.path.exists(newpath):
	os.makedirs(newpath)

all_timestamp = [item[0:6] for item in all_folders if item[0] in '1234567890']
timestamp = [all_timestamp[i] for i in range(len(all_timestamp)) if int(all_timestamp[i]) >= time_start and int(all_timestamp[i]) <= time_stop and all_timestamp[i] not in remove]


print '### YOU\'RE ANALYSING PIEZO SCANS ### \n'
LW = []
Quality = []
t = []
drift = []
drift_error = []
frequency = []
standard_dev = []





for i in range(len(timestamp)):
	time = timestamp[i]
	V,y = load_piezo_scan(folder = folder, timestamp = time)
	# V,y = zip(*sorted(zip(V, y)))
	if set_range_V == True:
		ind_min = np.where(V < V_min)
		ind_max = np.where(V > V_max)
		# print len(V)
		# print len(ind_min[0])
		V = np.delete(V, ind_min)
		V = np.delete(V, ind_max)
		y = np.delete(y, ind_min)
		y = np.delete(y, ind_max)
		# print len(V)

	max_index, max_value = max(enumerate(y), key = operator.itemgetter(1))

	
	'''
	Determine mean and standard deviation of the peaks
	'''
	# norm = y/sum(y)
	# mean = sum([V[i]*norm[i] for i in range(len(norm))])
	# variance = sum([norm[i]*(V[i]-mean)**2 for i in range(len(norm))])
	# st_dev = np.sqrt(variance)
	



	peak = np.where(y > 0.5)
	X = [V[peak[0][i]] for i in range(len(peak[0]))]
	mean = sum(X)/len(X)
	variance = sum([(x-mean)**2 for x in X])/len(X)
	st_dev = np.sqrt(variance)
	st_dev_3 = 3*st_dev
	LW.append(st_dev_3)

	X1 = mean - (st_dev_3/2)
	X2 = mean + (st_dev_3/2)




	# V1 = V[peak[0][0]]
	# V2 = V[peak[0][-1]]

	
	
	

	Q1 = 637./st_dev
	freq_st_dev = (((3e8/637e-9)/Q1)/1e9)
	standard_dev.append(st_dev)







	fig,ax = plt.subplots(figsize=(6,4.7))
	# ax.plot(V,norm, 'b', label = 'piezo scan')
	ax.plot(V,y,'r')
# 		plot.plot_fit1d(fit_result, np.linspace(V[0],V[-1],len(V)), ax=ax, show_guess=True, plot_data=False)
	# ax.plot(V,function,'g', label = 'testplot')
	ax.legend()
	plt.plot(X1,0, 'go')
	plt.plot(X2,0, 'go')
# 		plt.plot(V[Half_Max_r],A/2, 'bo')
# 		plt.plot(V[Half_Max_l], A/2, 'bo')
	plt.xlabel('Voltage [V]')
	plt.ylabel('Transmission')
	# ax.set_xlim(V_min,V_max)
	plt.title(time + '\n' + measurement_type + '\n Linewidth is %s V' %(round(st_dev_3,5)))
	fig.savefig(newpath +'/' + time + ".piezoplot2.png")

	
''' 
When all the plots are fitted, the drift can be determined and fitted by means of the fitting parameter x0, 
which is expected to represent the middle of the peak. These drift values are plotted over time, whereby the
error found in the fitting is considered as well.
'''




''' Determine linewidth by means of sigma
'''
average_sigma = sum(LW)/float(len(LW))
variance_average_sigma = sum([(x-average_sigma)**2 for x in LW])/len(LW)
st_dev_average_sigma = np.sqrt(variance_average_sigma)

# print sum(standard_dev)/len(standard_dev)



print 'AVERAGE SIGMA IS ', round(average_sigma,5), 'WITH A ST DEV OF ', round(st_dev_average_sigma,5)

	


# 	print 'The average linewidth for timestamps %s through %s is %s nm with a standard deviation of %s.' %(time_start, time_stop, linewidth_average, linewidth_stdev) 
# 	# print 'The average quality of the cavity is %s, with a standard deviation of %s.' %(time_start, time_stop, Quality_average, Q_stdev)

