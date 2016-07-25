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
date = '20150716'
measurement_type = 'RT_tubeoff_Nin'
scan_type = 'piezo_scans'
# 'piezo_scans'

''' Fitting yes/no '''
fitting = 1
fitting_method = 'gauss'
fixed = [0]

''' Select the time stamps you want to include or remove'''
time_start = 100000
time_stop = 235900
# remove = ['201245','201351','201507','201949','202023','202450', '203635', '203718', '203632', '203630','203634', '203644', '203659', '212152', '212206']
remove = ['185132','192231','191401','201010','201015','201022','201035','201040','201245','201351','201428','201507','201949','202023','202450', '203635', '203644', '203718', '204314', '212152', '212157', '212206', '212543', '212618', '212652', '212618', '212229', '212543', '212322', '213108','184254', '213513','213626', '183919','184335','185004', '185059']
# remove = ['145047','145540','150244','150905','151009','151112', '145540', '145646', '151129', '150944', '151246', '145244','150347', '150122', '120429', '212322', '120517', '121146','122150', '122259', '121059']
# remove = ['175318','164643', '164647', '164657', '164659', '164703', '164715', '164717', '164718', '164720','155755', '160818', '161605', '162219','173210', '173212','173214', '173215', '173841', '154408', '154413', '154420', '154430', '154440', '154723']



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

folder = 'K:/ns/qt/Diamond/Projects/Cavities/data/' + date + '/' + measurement_type + '/' + scan_type
all_folders = [f for f in os.listdir(folder)]
all_folders.sort()

newpath =  'K:/ns/qt/Diamond/Projects/Cavities/data_fitting/vibration measurements/' + date + '/' + measurement_type + '/' + scan_type
if not os.path.exists(newpath):
	os.makedirs(newpath)

all_timestamp = [item[0:6] for item in all_folders if item[0] in '1234567890']
timestamp = [all_timestamp[i] for i in range(len(all_timestamp)) if int(all_timestamp[i]) >= time_start and int(all_timestamp[i]) <= time_stop and all_timestamp[i] not in remove]

if scan_type == 'piezo_scans':
	print '### YOU\'RE ANALYSING PIEZO SCANS ### \n'
	LW = []
	Quality = []
	t = []
	drift = []
	drift_error = []
	frequency = []
	standard_dev = []



	if fitting == True:

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
			norm = y/sum(y)
			mean = sum([V[i]*norm[i] for i in range(len(norm))])
			variance = sum([norm[i]*(V[i]-mean)**2 for i in range(len(norm))])
			st_dev = np.sqrt(variance)
			
			X1 = mean - (st_dev)
			X2 = mean + (st_dev)



			peak = np.where(y > 0.2)
			# V1 = V[peak[0][0]]
			# V2 = V[peak[0][-1]]

			X = [V[peak[0][i]] for i in range(len(peak[0]))]
			
			

			Q1 = 637./st_dev
			freq_st_dev = (((3e8/637e-9)/Q1)/1e9)
			standard_dev.append(st_dev)


	
	


			# V = V[0:len(V):10]
			# y = y[0:len(y):10]

			'''
			Start fitting
			'''
			

			offset = 0.00
			amplitude = max_value
			x0 = V[max_index]
			exponent = 1
			decay_constant = 3
			sigma = 0.0004
			gamma = 0.5


			if fitting_method == 'gauss':
				p0, fitfunc, fitfunc_str = common.fit_gauss(offset, amplitude, x0, sigma)
			elif fitting_method == 'exp':
				p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, x0, decay_constant,exponent)
			elif fitting_method == 'lorentz':
				p0, fitfunc, fitfunc_str = common.fit_lorentz(offset, amplitude, x0, gamma)

			fit_result = fit.fit1d(V,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

			# # 	# chi2 = fit_result['chisq']
			# # 	# print "chi squared is %s" %chi2


			''' 
				In order to calculate the linewidth a help-function is created which resembles the fitting function (by means of
				the fitting parameters). Of this function all indices where the function exceeds the half maximum are stored in
				Half_Max. Of this list the very first and very last values are considered to by the indices at the half maximum and
				the linewidth (in nm) is then calculated by subtracting the voltages of those points and multiplying that by 25 nm/V, 
				a spec that is provided by JPE. Subsequently Q can be calculated by the well-known formula Q=omega/(delta omega), for 
				which omega is the wavelength of the laser and delta omega the linewidth
			'''

			A = fit_result['params_dict']['A']
			a = offset
			x0 = fit_result['params_dict']['x0']
			sigma = fit_result['params_dict']['sigma']

			function = a + A*np.exp(-(V-x0)**2/(2*sigma**2))

				# 25 nm might not be the right value, this should be looked into more thoroughly
			Half_Max = np.where(function > A/2)
			Half_Max_l = Half_Max[0][0]
			Half_Max_r = Half_Max[0][-1]
			Linewidth = (abs(V[Half_Max_r]-V[Half_Max_l])*25)
			Q = 637./Linewidth
			freq_LW = (((3e8/637e-9)/Q)/1e9)
			LW.append(Linewidth)
			Quality.append(Q)
			frequency.append(freq_LW)

			print 'BELANGRIJK!!! ST_DEV, DAN LINEWIDTH ', st_dev*25, Linewidth, '\n'


			''' Determine the drift and the corresponding times'''
			t0 = int(timestamp[0][0:2])*3600 + int(timestamp[0][2:4])*60 + int(timestamp[0][4:6])
			t_nu = (int(time[0:2])*3600 + int(time[2:4])*60 + int(time[4:6]))
			t.append(abs(t_nu-t0))
			# drift.append(mean)
			drift.append(fit_result['params_dict']['x0'])
			drift_error.append(fit_result['error_dict']['x0'] if fit_result['error_dict']['x0'] < 1e-4 else 0)




			fig,ax = plt.subplots(figsize=(6,4.7))
			ax.plot(V,norm, 'b', label = 'piezo scan')
			ax.plot(V,y,'r')
	# 		plot.plot_fit1d(fit_result, np.linspace(V[0],V[-1],len(V)), ax=ax, show_guess=True, plot_data=False)
			# ax.plot(V,function,'g', label = 'testplot')
			ax.legend()
			plt.plot(mean,0, 'go')
			# plt.plot(X2,0, 'go')
	# 		plt.plot(V[Half_Max_r],A/2, 'bo')
	# 		plt.plot(V[Half_Max_l], A/2, 'bo')
			plt.xlabel('Voltage [V]')
			plt.ylabel('Transmission')
			ax.set_xlim(V_min,V_max)
			plt.title(time + '\n' + measurement_type + '\n Linewidth is %s V' %(round(st_dev,5)))
			fig.savefig(newpath +'/' + time + ".piezoplot2.png")

			
		''' 
		When all the plots are fitted, the drift can be determined and fitted by means of the fitting parameter x0, 
		which is expected to represent the middle of the peak. These drift values are plotted over time, whereby the
		error found in the fitting is considered as well.
		'''

		drift_slope = (t[-1] - t[0])/(drift[-1] - drift[0])
		drift_offset = t[0]-(drift_slope*drift[0])
		fixed = [0]

		drift = np.array(drift)
		t = np.array([float(x) for x in t])

		p0, fitfunc, fitfunc_str = common.fit_line(drift_offset, drift_slope)
		drift_fit_result = fit.fit1d(drift,t, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

		fig, ax = plt.subplots(figsize=(6,4.7))
		plt.errorbar(drift, t, xerr = drift_error, ecolor = 'r')
		plot.plot_fit1d(drift_fit_result, drift, ax=ax, show_guess=True, plot_data=False)
		plt.xlabel('Volt [V]')
		plt.ylabel('Time [s]')
		plt.title('Drift for %s from timestamp %s through %s' %(measurement_type, time_start, time_stop))
		fig.savefig(newpath +'/' + str(time_start) + ".### DRIFTPLOT ###.png")



















	# 	''' Finally, the linewidth of the data of selected timestamps is averaged. '''
	# 	print 'LW = ', LW
	# 	print 'Q = ', Quality
	# 	linewidth_average = sum(LW)/float(len(LW))
	# 	Quality_average = sum(Quality)/float(len(Quality))
	# 	LW_freq_av = sum(frequency)/float(len(frequency))

	# 	var_LW = []
	# 	var_Q = []
	# 	for i in range(len(LW)):
	# 		var_LW.append((LW[i]-linewidth_average)**2)
	# 		# var_Q.append((Q[i]-Quality_average)**2)

	# 	linewidth_stdev = np.sqrt(sum(var_LW)/float(len(LW)))
	# 	# Q_stdev = np.sqrt(sum(var_Q)/float(len(Q)))

	# 	linewidth_av_freq = ((3e8/637e-9)/(637./linewidth_average))/1e9
	# 	linewidth_av_freq2 = ((3e8/637e-9)/Quality_average)/1e9

	# 	print 'LW_freq', linewidth_av_freq, linewidth_av_freq2, LW_freq_av

	# 	# LW = 0.084
	# 	# Q = 637./LW
	# 	# Freq = 3e8/(637e-9)
	# 	# Linewidth = Freq/Q
	# 	# print Linewidth/1e9

		''' Determine linewidth by means of sigma
		'''
		average_sigma = sum(standard_dev)/float(len(standard_dev))
		variance_average_sigma = sum([(x-average_sigma)**2 for x in standard_dev])/len(standard_dev)
		st_dev_average_sigma = np.sqrt(variance_average_sigma)

		print sum(standard_dev)/len(standard_dev)



		print 'AVERAGE SIGMA IS ', round(average_sigma,5), 'WITH A ST DEV OF ', round(st_dev_average_sigma,5)

			


	# 	print 'The average linewidth for timestamps %s through %s is %s nm with a standard deviation of %s.' %(time_start, time_stop, linewidth_average, linewidth_stdev) 
	# 	# print 'The average quality of the cavity is %s, with a standard deviation of %s.' %(time_start, time_stop, Quality_average, Q_stdev)


	else:
		for i in range(len(all_timestamp)):
			time = all_timestamp[i]
			V,y = load_piezo_scan(folder = folder, timestamp = time)
			fig,ax = plt.subplots(figsize=(6,4.7))
			ax.plot(V,y)
			plt.xlabel('Volt [V]')
			plt.ylabel('Transmission')
			plt.title(time + '\n Plot without fit')

	#################################################################################################################################################################
	#################################################################################################################################################################

elif scan_type == 'laser_scans':
	print '### YOU\'RE ANALYSING LASER SCANS ### \n'
	LW = []
	Quality = []
	t = []
	drift = []
	drift_error = []
	delete = []
	standard_dev = []


	if fitting == True:

		for i in range(len(timestamp)):
			time = timestamp[i]
			f,y = load_lr_scan (folder=folder, timestamp = time)
			# f,y = zip(*sorted(zip(f, y)))

			if set_range_f == True:
				ind_min = np.where(f < f_min)
				ind_max = np.where(f > f_max)
				# print len(ind_min[0])
				f = np.delete(f, ind_min)
				f = np.delete(f, ind_max)
				y = np.delete(y, ind_min)
				y = np.delete(y, ind_max)
				# print len(f)

			max_index, max_value = max(enumerate(y), key = operator.itemgetter(1))
			if max_value < 0.02:
				continue


			'''
			Determine mean and standard deviation of the peaks
			'''
			peak = np.where(y > 0.2)
			X = [f[peak[0][i]] for i in range(len(peak[0]))]
			mean = sum(X)/len(X) # in Volt
			variance = sum([(x - mean)**2 for x in X])/len(X)
			st_dev = np.sqrt(variance) # in nm
			
			X1 = mean - (st_dev)
			X2 = mean + (st_dev)

			standard_dev.append(st_dev)


			f,y = zip(*sorted(zip(f, y)))

			# f_min = max(a-200,min(f))
			# f_max = min(a+200,max(f))

			# f = f[0:len(f):10]
			# y = y[0:len(y):10]

			offset = 0.00
			amplitude = max_value
			x0 = f[max_index]
			exponent = 1
			decay_constant = 3
			sigma = 2.413625
			gamma = 20

			if fitting_method == 'gauss':
				p0, fitfunc, fitfunc_str = common.fit_gauss(offset, amplitude, x0, sigma)
			elif fitting_method == 'exp':
				p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, x0, decay_constant,exponent)
			elif fitting_method == 'lorentz':
				p0, fitfunc, fitfunc_str = common.fit_lorentz(offset, amplitude, x0, gamma)

			fit_result = fit.fit1d(f,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)


			# # 	# chi2 = fit_result['chisq']
			# # 	# print "chi squared is %s" %chi2

			''' 
				In order to calculate the linewidth a help-function is created which resembles the fitting function (by means of
				the fitting parameters). Of this function all indices where the function exceeds the half maximum are stored in
				Half_Max. Of this list the very first and very last values are considered to by the indices at the half maximum and
				the linewidth (in nm) is then calculated by subtracting the voltages of those points and multiplying that by 25 nm/V, 
				a spec that is provided by JPE. Subsequently Q can be calculated by the well-known formula Q=omega/(delta omega), for 
				which omega is the wavelength of the laser and delta omega the linewidth
			'''
			
			A = fit_result['params_dict']['A']
			a = offset
			x0 = fit_result['params_dict']['x0']
			sigma = fit_result['params_dict']['sigma']

			function = a + A*np.exp(-(f-x0)**2/(2*sigma**2))

			Half_Max = np.where(function > A/2)
			Half_Max_l = Half_Max[0][0]
			Half_Max_r = Half_Max[0][-1]

			# print 'Half_Max_l', f[Half_Max_l], 'Half_Max_r', f[Half_Max_r]

			Linewidth = abs(f[Half_Max_r]-f[Half_Max_l])
			LW.append(Linewidth)


			''' Determine the drift and the corresponding times'''
			t0 = int(timestamp[0][0:2])*3600 + int(timestamp[0][2:4])*60 + int(timestamp[0][4:6])
			t_nu = (int(time[0:2])*3600 + int(time[2:4])*60 + int(time[4:6]))
			t.append(abs(t_nu-t0))
			drift.append(fit_result['params_dict']['x0'])
			drift_error.append(fit_result['error_dict']['x0'] if fit_result['error_dict']['x0'] < 50 else 0)

			fig,ax = plt.subplots(figsize=(6,4.7))
			ax.plot(f,y, 'r', label = 'laser scan')
			# plot.plot_fit1d(fit_result, np.linspace(f[0],f[-1],len(f)), ax=ax,label='Fit',show_guess=True, plot_data=False)
			# ax.plot(f,function, 'g', label='testplot')
			ax.legend()
			# plt.plot(f[Half_Max_r], A/2, 'bo')
			# plt.plot(f[Half_Max_l], A/2, 'bo')
			# plt.plot(X1,0, 'go')
			# plt.plot(X2,0, 'go')
			plt.xlabel('frequency [GHz]')
			plt.ylabel('Transmission')
			# plt.xlim(f_min,f_max)
			plt.title(time + '\n' + measurement_type + '\n Linewidth is %s GHz' %st_dev)
			fig.savefig(newpath +'/' + time + ".laserplot.png")

		
		''' 
		When all the plots are fitted, the drift can be determined and fitted by means of the fitting parameter x0, 
		which is expected to represent the middle of the peak. These drift values are plotted over time, whereby the
		error found in the fitting is considered as well.
		'''

		drift_slope = (t[-1] - t[0])/(drift[-1] - drift[0])
		drift_offset = t[0]-(drift_slope*drift[0])
		fixed = [0]

		drift = np.array(drift)
		t = np.array([float(x) for x in t])

		p0, fitfunc, fitfunc_str = common.fit_line(drift_offset, drift_slope)
		drift_fit_result = fit.fit1d(drift,t, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

		fig, ax = plt.subplots(figsize=(6,4.7))
		plt.errorbar(drift, t, xerr = drift_error, ecolor = 'r')
		plot.plot_fit1d(drift_fit_result, drift, ax=ax, show_guess=True, plot_data=False)
		plt.xlabel('frequency [GHz]')
		plt.ylabel('Time [s]')
		plt.title('Drift for ' + measurement_type)
		fig.savefig(newpath +'/' + str(time_start) + ".### DRIFTPLOT ###.png")


		# ''' Finally, the linewidth of the data of selected timestamps is averaged. '''
		# print 'LW = ', LW
		# linewidth_average = sum(LW)/float(len(LW))
		# print 'The average linewidth for timestamps %s through %s is %s GHz' %(time_start, time_stop, linewidth_average)

		''' 
		Determine linewidth by means of sigma
		'''
		average_sigma = sum(standard_dev)/float(len(standard_dev))
		variance_average_sigma = sum([(x-average_sigma)**2 for x in standard_dev])/len(standard_dev)
		st_dev_average_sigma = np.sqrt(variance_average_sigma)

		print 'AVERAGE SIGMA IS ', round(average_sigma,2), 'WITH A ST DEV OF ', round(st_dev_average_sigma,2)


	else:
		for i in range(len(all_timestamp)):
			time = all_timestamp[i]


			f,y = load_lr_scan (folder=folder, timestamp = time)

			if set_range_f == True:
				ind_min = np.where(f < f_min)
				ind_max = np.where(f > f_max)
				# print len(ind_min[0])
				f = np.delete(f, ind_min)
				f = np.delete(f, ind_max)
				y = np.delete(y, ind_min)
				y = np.delete(y, ind_max)
				# print len(f)
			f,y = zip(*sorted(zip(f, y)))
			fig,ax = plt.subplots(figsize=(6,4.7))
			ax.plot(f,y)
			plt.xlabel('frequency [GHz]')
			plt.ylabel('Transmission')
			plt.title(time + '\n Plot without fit')
