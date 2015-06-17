#from analysis.lib.tools.oscilloscope_cvs_trace import CVSTrace as csvTrace
from analysis.lib.cavity import cavity_tools as cavTools
import numpy as np
import pylab as plt
from matplotlib import rc, cm
import csv

reload(cavTools)

matplotlib.rc ('xtick', labelsize=15)
matplotlib.rc ('ytick', labelsize=15)


class CVSTrace():

	def __init__(self):
		self.trace = None
		self.x_axis = None

	def load_trace(self, filename, delimiter = ',', quotechar='|', scope_type = None):

		with open(filename, 'rb') as csvfile:
			row_list = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)
			ind = 0
			x = []
			y = []
			for row in row_list:
				ind = ind + 1
				try:
					if (scope_type=='yoko'):
						y.append (float(row[0]))
						print row[0], row[1]
					else:
						x.append(float(row[0]))
						y.append(float(row[1]))
				except:
					pass
					#print 'row ', ind, ' is not a number'

		self.trace = np.asarray(y)
		self.x_axis = np.asarray (x)

	def load_trace_yoko (self, filename, delimiter = ',', quotechar='|'):
		print 'Loading YOKO data...', filename
		with open(filename, 'rb') as csvfile:
			row_list = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)
			ind = 0
			x = []
			y = []
			for row in row_list:
				try:
					if (ind>15):
						y.append( float(row[1]))
				except:
					print 'no value'
				#if (ind<15):
				#	print ind, row
				if (ind==8):
					self.h_res = float(row[1])
				if (ind==7):
					self.sampling_rate = float(row[1])
				ind = ind+1
		self.trace = np.asarray(y)
		self.x_axis = np.arange(len(self.trace))*self.h_res
		print 'resolution: ', self.h_res*1e8, '[ns]'

	def plot_trace (self):
		plt.figure(figsize = (10,6))
		plt.plot (self.x_axis, self.trace, '.b')
		plt.plot (self.x_axis, self.trace, 'r')
		plt.show()



def time_trace(folder, file_name, limits = None, do_save = False, correlate = False):

	#folder = r'D:/Research/cavity/low_temperature/20150605/oscilloscope_traces/take_1_pulsetubeON/'
	#fname = "NewFile"+str(idx)+".csv"

	t = CVSTrace()
	t.load_trace_yoko(filename=folder+file_name)
	x = t.x_axis
	y = t.trace

	f1 = plt.figure (figsize=(25, 4))
	plt.plot  (x*1000,y, '.b')
	plt.plot (x*1000,y,'r')
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

def frequency_scan(folder, timestamp, do_save=False):

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

	#if do_limits:
	#	plt.ylim ([1, 200])
	plt.xlabel ('freq [GHz]', fontsize=15)
	if do_save:
		f3.savefig (folder+'sum_fft_'+str(idx)+'.png')
	plt.show()

def piezo_scan():
	folder = r'M:/tnw/ns/qt/Diamond/Projects/Cavities/data/20150611/'
	fnames, times = cavTools.get_files_in_folder (contains='piezo_scan', folder=folder)
	#cavTools.combine_piezoscans_1D (folder=folder, folder_array=fnames, time_array=times)#, do_save=True)
	cavTools.track_single_peak (folder=folder, folder_array=fnames, time_array=times)#, do_save=True)

def plot_single_piezoscan (folder, timestamp):
	plt.figure (figsize=(20,5))
	x, y = cavTools.load_data (folder=folder, timestamp=timestamp, scan_type='piezo')
	plt.plot (x,y, 'ob')
	plt.plot (x,y, 'r')

	#plt.xlim ([-1.7, -1.65])
	plt.xlabel ('piezo voltage [V]', fontsize=15)
	plt.show()

def plot_2D_piezo_laser (folder, timestamp):
	cavTools.process_2D_scan (folder=folder, timestamp=timestamp)
#piezo_scan()
#frequency_scan()

folder = r'D:/Research/cavity/low_temperature/20150617/'
#plot_single_piezoscan (folder = folder, timestamp = '142333')
#folder = r'M:/tnw/ns/qt/Diamond/Projects/Cavities/data/20150611/'
frequency_scan (folder=folder, timestamp='145041')
#plot_2D_piezo_laser (folder=folder, timestamp='120548')
'''
folder = r'D:/Research/cavity/low_temperature/20150617/scope_traces/'

fname = 'U1A027.csv'
time_trace (folder=folder, file_name=fname, limits = None, correlate = True)#[40, 60])
'''