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
import csv


reload(toolbox)


'''
Uitleg van deze functie
# '''
# with open('K:/ns/qt/Diamond/Projects/Cavities/data/20150630/RT_tubeoff_Nin/scope_traces/U1A000','rb') as file:
# 	contents = csv.reader(file)
# 	[x for x in contents]


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




def load_scope_traces(folder, sample):
	all_files = [ f for f in os.listdir(folder) if sample in f]
	current_file = os.path.join(folder, all_files[0])
	# file_name = [f for f in current_file if 'U1' in f]

	y = np.loadtxt(current_file, delimiter=',', skiprows=16, usecols=(1,))
	# t = np.loadtxt(dtype = str, fname = current_file, delimiter = ',', usecols = (1,))

	y = np.ravel(y)
	T_sampling = 1.6e-6
	# float(t[8]) 

	return y, T_sampling
##########################################################################################################################################################################
##########################################################################################################################################################################

# date = '20150716'
# measurement_type = 'RT_tubehigh_Nin'
# scan_type = 'scope_traces'
remove = []

steps = ['9999', '2000','500','500','5000','5000_higher_res']

folder = 'K:/ns/qt/Diamond/Projects/Cavities/data/piezo_linearity_(forAM)'
# # folder = 'K:/ns/qt/Diamond/Projects/Cavities/data/' + date + '/' + measurement_type + '/' + scan_type
all_files = [f for f in os.listdir(folder)]
all_files.sort()
all_samples = [item[0:6] for item in all_files if item[0] == 'U']
samples = [all_sample[i] for i in range(len(all_samples)) if all_sample[i] not in remove]




for i in range(len(samples)):
	sample = samples[i]

	print len(steps), 'steps', len(samples)
	print steps[i]
	
	y, T_sampling = load_scope_traces(folder = folder, sample = sample)
	F_sampling = 1./T_sampling
	n = len(y)
	# print F_sampling
	t = np.linspace(0, n*T_sampling, n)

	fig = plt.figure()

	# ax1 = fig.add_subplot(211)
	plt.plot(t,y, 'b-')
	plt.xlabel('time')
	plt.ylabel('Voltage output adwin')
	plt.title('scope_traces for ' + steps[i] + 'steps')
	plt.show()



# newpath =  'K:/ns/qt/Diamond/Projects/Cavities/data_fitting/vibration measurements/' + date + '/' + measurement_type + '/' + scan_type
# if not os.path.exists(newpath):
# 	os.makedirs(newpath)


# for i in range(len(sample)):
# 	s = sample[i]
# 	print s

# 	if len(s) > 12:
# 		print '### FROM HERE NEW KIND OF MEASUREMENTS ###'

# 	y, T_sampling = load_scope_traces(folder = folder, sample = s)
# 	# print y
# 	# print T_sampling

# 	F_sampling = 1./T_sampling
# 	n = len(y)
# 	print F_sampling
# 	t = np.linspace(0, n*T_sampling, n)

	

# 	Y = np.fft.fft(y)/n
# 	Y = Y[range(n/2)]
# 	Y = abs(Y)
# 	X = np.fft.fftfreq(n, T_sampling)
# 	X = X[range(n/2)]
# 	X = np.delete(X,0)
# 	Y = np.delete(Y,0)

# 	ind = np.where(Y > 0.005)
# 	place = X[ind]

# 	''' Add all scope traces for 1 measurement type '''
# 	if i == 0:
# 		C = Y
# 	else:
# 		C = C + Y




# 	frequencies = [g for g in place if g > 200]

# 	frequencies2 = []
# 	print len(frequencies) - 1
# 	for g in range((len(frequencies)-1)):
# 		if frequencies[g+1] - frequencies[g] > 100:
# 			frequencies2.append(frequencies[g])

# 	if len(frequencies) > 0:
# 		frequencies2.append(frequencies[-1])
# 	print frequencies2



# 	fig = plt.figure()

# 	ax1 = fig.add_subplot(211)
# 	ax1.plot(t,y, 'k-')
# 	plt.xlabel('time')
# 	plt.ylabel('amplitude [V]')
# 	plt.title('scope_traces \n' + date + ' ' + measurement_type + ' ' + s)

# 	ax2 = fig.add_subplot(212)
# 	ax2.plot(X,Y, 'b-')
# 	plt.xlabel('freq(Hz)')
# 	plt.ylabel('|Y(freq)|')
# 	plt.title(date + '\n Cavity vibrations for ' + measurement_type + ' at ' + s)
# 	plt.xlim(0,30000)
# 	plt.ylim(0,0.05)

# 	fig.subplots_adjust(hspace=0.8)
# 	plt.show()
# 	fig.savefig(newpath +'/' + s + ".scopeplot.png")

# fig1, ax = plt.subplots(figsize=(6,4.7))
# plt.plot(X,C,'r-')
# plt.xlabel('time')
# plt.ylabel('amplitude [V]')
# plt.xlim(0,20000)
# plt.show()
# fig1.savefig(newpath +'/' + measurement_type + ".scopeplotTOTAL.png")


# ind = np.where(C > 0.02)
# place = X[ind]

# frequencies3 = [g for g in place if g > 200]
# print frequencies3


# # frequencies4 = []
# # print len(frequencies3) - 1
# # for g in range((len(frequencies3)-1)):
# # 	if frequencies3[g+1] - frequencies3[g] > 100:
# # 		frequencies4.append(frequencies[g])

# # # if len(frequencies) > 0:
# # # 	frequencies2.append(frequencies[-1])
# # print frequencies4




