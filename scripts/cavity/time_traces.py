from analysis.lib.tools.oscilloscope_cvs_trace import CVSTrace as csvTrace
from analysis.lib.cavity import cavity_tools as cavTools
import numpy as np
import pylab as plt
from matplotlib import rc, cm
import csv

reload(cavTools)

matplotlib.rc ('xtick', labelsize=15)
matplotlib.rc ('ytick', labelsize=15)

def time_trace(folder, file_name, do_limits = False, do_save = False):

	#folder = r'D:/Research/cavity/low_temperature/20150605/oscilloscope_traces/take_1_pulsetubeON/'
	#fname = "NewFile"+str(idx)+".csv"

	t = csvTrace()
	t.load_trace_yoko(filename=folder+file_name)
	x = t.x_axis
	y = t.trace

	f1 = plt.figure (figsize=(12, 4))
	plt.plot (x*1000, y, '.b')
	plt.plot (x*1000,y,'r')
	plt.xlabel ('time [ms]', fontsize=15)
	plt.ylabel ('photodiode signal', fontsize=15)
	if do_limits:
		plt.xlim ([-10, 50])
	if do_save:
		f1.savefig (folder+'time_trace_'+str(idx)+'.png')
	plt.show()

	yF = np.fft.fftshift(np.fft.fft(y, n=1000000))
	f = 0.001*np.linspace (-1, 1,  len(yF))*(1./(x[2]-x[1]))
	yF = yF[500500:]
	f = f[500500:]

	f2 = plt.figure (figsize=(12, 4))
	plt.plot (f, np.abs(yF)**2)
	plt.axis("tight")
	plt.xlim([0.001, 20])

	if do_limits:
		plt.ylim ([1, 200])
	plt.xlabel ('freq [kHz]', fontsize=15)
	if do_save:
		f2.savefig (folder+'fft_'+str(idx)+'.png')
	plt.show()

	return x, y, f, yF

def frequency_scan():
	folder = r'D:/Research/cavity/low_temperature/20150605/'
	folder = r'D:/Research/cavity/room_temperature/20150601/'
	f, y = cavTools.load_data (folder=folder, timestamp='101830', scan_type='lr_scan')
	f = np.ravel(f)
	y = np.ravel(y)
	ind = np.argsort(f)
	f = f[ind]
	y=y[ind]
	do_save = False
	f4 = plt.figure (figsize=(8, 5))
	plt.plot (f, y, 'ob')
	plt.plot (f, y, 'r')
	plt.axis("tight")
	plt.xlim([140, 200])

	#if do_limits:
	#	plt.ylim ([1, 200])
	plt.xlabel ('freq [GHz]', fontsize=15)
	if do_save:
		f3.savefig (folder+'sum_fft_'+str(idx)+'.png')
	plt.show()

def piezo_scan():
	folder = r'D:/Research/cavity/low_temperature/20150605/piezo_scan_7/'
	fnames, times = cavTools.get_files_in_folder (contains='piezo_scan', folder=folder)
	#cavTools.combine_piezoscans_1D (folder=folder, folder_array=fnames, time_array=times)#, do_save=True)
	cavTools.track_single_peak (folder=folder, folder_array=fnames, time_array=times)#, do_save=True)

piezo_scan()
#frequency_scan()

'''
folder = r'D:/Research/cavity/low_temperature/20150605/piezo_scan_7/'
plt.figure (figsize=(15,5))
x, y = cavTools.load_data (folder=folder, timestamp='165056', scan_type='piezo')
plt.plot (x,y, 'ob')
plt.plot (x,y, 'r')

plt.xlim ([6., 6.07])
plt.xlabel ('piezo voltage [V]', fontsize=15)
plt.show()
'''
#cavTools.piezoscan_movie()
'''
folder = r'D:/Research/cavity/low_temperature/20150605/high_res_scope/'
fname = 'U1A004.csv'

t = CSVTrace()
t.load_trace(filename=folder+fname, scope_type = 'yoko')
#x = t.x_axis
y = t.trace[10:]
x = np.arange(len(y))*t.h_res

f1 = plt.figure (figsize=(18, 4))
plt.plot (x*1000, y, '.b')
plt.plot (x*1000, y,'r')
plt.xlim([15, 20])
plt.xlabel ('time [s]', fontsize=15)
plt.ylabel ('photodiode signal', fontsize=15)
plt.show()


yF = np.fft.fftshift(np.fft.fft(y, n=5000000))
f = 0.001*np.linspace (-1, 1,  len(yF))*(1./(x[2]-x[1]))
yF = yF[2500500:]
f = f[2500500:]

f2 = plt.figure (figsize=(18, 4))
plt.plot (f, np.abs(yF)**2)
plt.axis("tight")
plt.xlim([0., 100])
plt.xlabel ('f [KHz]', fontsize=15)
plt.show()
'''