
import h5py
from analysis.lib.tools import toolbox
import pylab as plt
from analysis.lib import fitting
from analysis.lib.fitting import fit,esr, common
from matplotlib import rc, cm
from scipy import interpolate
from analysis.lib.tools import oscilloscope_cvs_trace as scope_trace

reload(toolbox)
reload (scope_trace)


def analize_time_trace(folder, file_name, limits = None, do_save = False, correlate = False):

	#folder = r'D:/Research/cavity/low_temperature/20150605/oscilloscope_traces/take_1_pulsetubeON/'
	#fname = "NewFile"+str(idx)+".csv"

	t = scope_trace.CVSTrace()
	t.load_trace_yoko(filename=folder+file_name)
	x = t.x_axis
	y = t.trace

	f1 = plt.figure (figsize=(25, 4))
	plt.plot  (x*1000, y, '.b')
	plt.plot (x*1000, y,'r')
	plt.xlabel ('time [ms]', fontsize=15)
	plt.ylabel ('photodiode signal', fontsize=15)
	#plt.xlim([80000,82000])
	if limits:
		plt.xlim (limits)
	if do_save:
		f1.savefig (folder+'time_trace_'+str(idx)+'.png')
	plt.show()

	yF = np.fft.fftshift(np.fft.fft(y, n=1000000))
	f = 0.001*np.linspace (-1, 1,  len(yF))*(1./(x[2]-x[1]))
	#yF = yF[500500:]
	#f = f[500500:]

	f2 = plt.figure (figsize=(12, 4))
	plt.plot (f, np.abs(yF)**2)
	plt.axis("tight")
	plt.xlim([0.001, 20])

	#if do_limits:
	plt.ylim ([1, 1000000])
	plt.xlabel ('freq [kHz]', fontsize=15)
	if do_save:
		f2.savefig (folder+'fft_'+str(idx)+'.png')
	plt.show()

	if correlate:
		y = np.asarray(y)
		yC = np.correlate (y, y, 'full')
		xC = np.arange (len(yC))*t.h_res*1000000
		plt.figure (figsize=(25,5))
		n = len(yC)/2
		plt.plot (xC[n+1:]-xC[n],yC[n+1:])
		plt.xlim([0,300])
		plt.xlabel ('time [us]', fontsize =15)
		plt.show()


	return x, y, f, yF

def frequency_scan(folder, timestamp, do_save=False, xlimit=None):

	f, y = cavTools.load_data (folder=folder, timestamp=timestamp, scan_type='lr_scan')
	f = np.ravel(f)
	y = np.ravel(y)
	ind = np.argsort(f)
	f = f[ind]
	y=y[ind]
	do_save = False
	f4 = plt.figure (figsize=(15, 5))
	plt.plot (f, y, 'ob')
	plt.plot (f, y, 'r')
	plt.axis("tight")
	#plt.xlim([-2000, -1500])

	if xlimit:
		plt.xlim (xlimit)
	plt.xlabel ('freq [GHz]', fontsize=15)
	if do_save:
		f3.savefig (folder+'sum_fft_'+str(idx)+'.png')
	plt.show()

def track_piezo_scan(folder, contains, threshold, selection = None):
	fnames, times = cavTools.get_files_in_folder (contains=contains, folder=folder)
	#cavTools.combine_piezoscans_1D (folder=folder, folder_array=fnames, time_array=times)#, do_save=True)
	if selection:
		cavTools.track_single_peak (folder=folder, folder_array=fnames[selection[0]:selection[1]], time_array=times[selection[0]:selection[1]], conv_factor = 1, threshold=threshold)#, do_save=True)
	else:
		cavTools.track_single_peak (folder=folder, folder_array=fnames[:], time_array=times[:], conv_factor = 1, threshold=threshold)#, do_save=True)


def plot_single_piezoscan (folder, timestamp, xlimit=None):
	plt.figure (figsize=(20,5))
	x, y = cavTools.load_data (folder=folder, timestamp=timestamp, scan_type='piezo')
	plt.plot (x,y, 'ob')
	plt.plot (x,y, 'r')

	if xlimit:
		plt.xlim (xlimit)
	plt.xlabel ('piezo voltage [V]', fontsize=15)
	plt.show()

def plot_2D_piezo_laser (folder, timestamp):
	cavTools.process_2D_scan (folder=folder, timestamp=timestamp)

def plot_piezoscan_in_time (folder, contains, threshold=None, xlimit=None, selection = None):
	fnames, times = cavTools.get_files_in_folder (contains=contains, folder=folder)
	if selection:
		cavTools.combine_piezoscans_1D (folder = folder, folder_array = fnames[selection[0]:selection[1]], time_array = times[selection[0]:selection[1]], threshold = threshold, xlimit=xlimit)
	else:
		cavTools.combine_piezoscans_1D (folder = folder, folder_array = fnames[:], time_array = times[:], threshold = threshold, xlimit=xlimit)

#frequency_scan()

#folder = r'D:/Research/cavity/low_temperature/20150617_emptyCav/'
#track_piezo_scan(folder=folder, contains = '_switchOFFcompr', threshold = 0.6, selection =[1, 40])

#plot_piezoscan_in_time (folder = folder, contains = 'N3_z=270_zoom')#, xlimit = [6,7], selection =[1,40])
#plot_single_piezoscan (folder = folder, timestamp = '123740')#, xlimit=[6.1, 6.3])
#folder = r'M:/tnw/ns/qt/Diamond/Projects/Cavities/data/20150611/'
#frequency_scan (folder=folder, timestamp='165011', xlimit=[-1800, -500])
#plot_2D_piezo_laser (folder=folder, timestamp='120548')


folder = r'M:/tnw/ns/qt/Diamond/Projects/Cavities/data/piezo_linearity_(forAM)/'

fname = 'U1A004.csv'
analize_time_trace (folder=folder, file_name=fname, limits = None, correlate = False)#[40, 60])




'''
folder = r'O:/'
fname = 'U1A005.csv'
t = CVSTrace()
t.load_trace_yoko(filename=folder+fname)
x = t.x_axis[100:50000]
y = t.trace[100:50000]

f1 = plt.figure (figsize=(25, 4))
plt.plot  (x*1000,y, '.b')
plt.plot (x*1000,y,'r')
plt.xlabel ('time [ms]', fontsize=15)
plt.ylabel ('photodiode signal', fontsize=15)
plt.show()
'''