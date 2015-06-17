
import h5py
from analysis.lib.tools import toolbox
import pylab as plt
from analysis.lib import fitting
from analysis.lib.fitting import fit,esr, common
from matplotlib import rc, cm
from scipy import interpolate

reload(toolbox)


def get_files (newer_than, older_than, contains, folder=None):
	a, b, folder_names = toolbox.latest_data (contains=contains, newer_than=newer_than, older_than=older_than, return_all=True, folder=folder)
	print folder_names
	times = []
	for f in folder_names:
		times.append (int(f[9:11])*3600 + int(f[11:13])*60 + int(f[13:15]))
	return folder_names, times

def get_files_in_folder (contains, folder):
	folder_names = [f for f in os.listdir(folder) if contains in f]
	times = []
	for f in folder_names:
		times.append (int(f[0:2])*3600 + int(f[2:4])*60 + int(f[4:6]))
	times = np.array(times)
	times = times-times[0]
	return folder_names, times


def load_data (folder, timestamp, scan_type):

	all_folders = [f for f in os.listdir(folder) if timestamp in f]
	print all_folders
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

def load_calibration (timestamp=''):
	if (timestamp==''):
		folder_name = toolbox.latest_data (contains='laser_calib')
	else:
		folder_name = toolbox.latest_data (contains=timestamp)
	
	all_files = [f for f in os.listdir(folder_name) if timestamp in f]  
	file_name = [f for f in all_files if '.hdf5' in f]
	
	f = h5py.File(os.path.join(folder_name, file_name[0]),'r')
	ls_grp = f['/calibration']
	V = ls_grp['V'].value
	freq = ls_grp['freq_GHz'].value
	return V, freq


def fit_calibration(V, freq):
	guess_b = (freq[-1]-freq[0])/(V[-1]-V[0])
	guess_a = freq[0]
	print 'Initial guess: ', guess_a, guess_b

	a = fit.Parameter(guess_a, 'a')
	b = fit.Parameter(guess_b, 'b')
	
	p0 = [a, b]
	fitfunc_str = ''

	def fitfunc(x):
		return a()+b()*x


	fit_result = fit.fit1d(V,freq, None, p0=p0, fitfunc=fitfunc, fixed=[],
        	do_print=False, ret=True)
	a_fit = fit_result['params_dict']['a']
	b_fit = fit_result['params_dict']['b']
	print 'a= ',a_fit
	print 'b=',b_fit
	return a_fit, b_fit

def combine_piezoscans (timestamp_array):

	voltages = np.array([])
	PD_signal = np.array([])

	for k in timestamp_array:
		V, Y = load_data (k)
		voltages = np.hstack ([voltages, V])
		PD_signal = np.hstack ([PD_signal, Y])

	f = plt.figure(figsize=(12,6))
	ax = f.add_subplot (1,1,1)
	ax.plot (voltages, PD_signal, 'b', linewidth=2)
	for tick in ax.xaxis.get_major_ticks():
	    tick.label.set_fontsize(16) 
	for tick in ax.yaxis.get_major_ticks():
	    tick.label.set_fontsize(16) 
	ax.set_xlabel ('piezo voltage [V]', fontsize= 16)
	ax.set_ylabel ('transmission [a.u.]', fontsize=16)
	plt.show()

	f = plt.figure(figsize=(12,6))
	ax = f.add_subplot (1,1,1)
	ax.plot (voltages[1300:1800], PD_signal[1300:1800], 'b', linewidth=2)
	for tick in ax.xaxis.get_major_ticks():
	    tick.label.set_fontsize(16) 
	for tick in ax.yaxis.get_major_ticks():
	    tick.label.set_fontsize(16) 
	ax.set_xlabel ('piezo voltage [V]', fontsize= 16)
	ax.set_ylabel ('transmission [a.u.]', fontsize=16)
	plt.show()

def combine_piezoscans_2D (folder, folder_array, y_axis = None, V_lims=None):

	voltages = np.array([])
	PD_signal = np.array([])
	i = 0
	for k in folder_array:
		V, Y = load_data (folder=folder, timestamp = k, scan_type='piezo')
		if (i==0):
			PD_signal = Y
		else:
			PD_signal = np.vstack ([PD_signal, Y])
		i = i + 1

	if y_axis == None:
		y_axis = np.arange(i)
	A, B = np.meshgrid (V, y_axis)

	f = plt.figure(figsize=(12,6))
	ax = f.add_subplot (1,1,1)
	ax.pcolor (A, B, PD_signal)
	for tick in ax.xaxis.get_major_ticks():
	    tick.label.set_fontsize(16) 
	for tick in ax.yaxis.get_major_ticks():
	    tick.label.set_fontsize(16) 
	if V_lims:
		ax.set_xlim(V_lims)
	ax.set_xlabel ('piezo voltage [V]', fontsize= 16)
	ax.set_ylabel ('time [s]', fontsize=16)
	plt.show()


def combine_piezoscans_1D (folder, folder_array, time_array):
	matplotlib.rc ('xtick', labelsize=20)
	matplotlib.rc ('ytick', labelsize=20)
	colori = cm.gist_earth(np.linspace(0,0.75, len(folder_array)))
	f = plt.figure(figsize=(8, 40))
	ax = f.add_subplot(1,1,1)
	i = 0
	for k in folder_array:
		#print k
		V, Y = load_data (folder=folder, timestamp = k, scan_type='piezo')
		ax.plot (V[5:500], Y[5:500]+time_array[i], '.', color=colori[i])
		ax.plot (V[5:500], Y[5:500]+time_array[i], color=colori[i])
		i= i+1
	ax.set_xlabel ('piezo voltage [V]', fontsize=20)
	ax.set_ylabel ('time [min]', fontsize=20)
	f.savefig ('D:/measuring/low_temp_cavity.png', dpi=300)
	plt.show()

def process_2D_scan (folder):
	file_name = [f for f in os.listdir(folder) if '.hdf5' in f]
	idx = 0
	done = False
	datafile = h5py.File(os.path.join(folder, file_name[0]),'r')
	nr_points = 300
	plt.figure(figsize = (12, 8))

	scan_data = np.zeros((10, nr_points))
	piezo_voltage = np.zeros(10)

	while idx<10:
		#try:
		curr_grp = datafile[str(idx)]
		f = curr_grp['frq'].value
		y = curr_grp['sgnl'].value
		piezo_voltage[idx] = curr_grp.attrs['pzV']
		#clean-up wavemeter errors 
		ind = find (abs(f)>5000)
		f = np.delete(f, ind)
		y = np.delete(y, ind)
		if (idx==0):
			frq = np.linspace (-300, 300, nr_points)
		interp_funct = interpolate.interp1d(f, y)
		y_new = interp_funct (frq)
		scan_data [idx, :] = y_new
		plt.plot (frq, y_new+0.1*idx, label = str(piezo_voltage[idx])+'V')
		idx = idx + 1
		#except:
		#	done = True
	datafile.close()
	plt.xlabel (' frq [GHz]', fontsize = 15)
	plt.legend()
	plt.show()

	FF, PZ = meshgrid (frq, piezo_voltage)
	plt.figure (figsize = (12, int(idx/5.)))
	plt.pcolor (FF, PZ, scan_data)
	plt.xlabel (' frq [GHz]', fontsize = 15)
	plt.ylabel ('piezo [V]', fontsize = 15)
	plt.show()


#t_array, times = get_files (contains='laser', older_than = '20150422_215600', newer_than = '20150422_215138')

#times = np.array(times)
#times = times - min(times)
#combine_piezoscans_2D (timestamp_array = t_array, y_axis = times)

#a, b = load_data (folder='D:/measuring/data/20150528/roomT_drift_mechExcit/', timestamp='105524', scan_type='piezo')

#f_array, t_array = get_files_in_folder (contains='piezo', folder = r'D:\measuring\data\20150529\PT_on')
#combine_piezoscans_2D (folder = r'D:\measuring\data\20150529\PT_on', folder_array = f_array, y_axis = t_array, V_lims = None)
#combine_piezoscans_1D (folder = r'D:\measuring\data\20150529\PT_on', folder_array = f_array, time_array = t_array/60.)
def load_lr_scan (folder, timestamp):
	f, y = load_data (folder = folder, timestamp = timestamp, scan_type = 'lr_scan')
	f = np.ravel(f)
	y = np.ravel(y)

	ind = find (abs(f)>5000)
	f = np.delete(f, ind)
	y = np.delete(y, ind)
	return f, y

folder = 'D:/measuring/data/20150615/'
f, y = load_lr_scan (folder=folder, timestamp = '124917')
plt.figure(figsize=(12,5))
plt.plot (f,y)
plt.xlim ([-1550, -1450])
plt.show()

#process_2D_scan (r'D:\measuring\data\201506152\124917_lr_scan_6360_6400')

'''
V, Y = load_data (timestamp = '112232')
Vf, freq = load_calibration ()
freq = freq-min(freq)

a_fit, b_fit = fit_calibration (V=Vf, freq=freq)
x_fit = np.linspace(V[0], V[-1], 1000)
y_fit = a_fit + b_fit*x_fit

plt.figure()
plt.plot (Vf, freq, 'b', linewidth=2)
plt.plot (x_fit, y_fit, 'r')
plt.show()

X = V*b_fit
guess_c = 10
guess_b = 5
guess_a = 50
a = fit.Parameter(guess_a, 'a')
b = fit.Parameter(guess_b, 'b')
c = fit.Parameter(guess_c, 'c')

p0 = [a, b, c]
fitfunc_str = ''

def fitfunc(x):
	return a()/((x-b())**2 + c()**2)


fit_result = fit.fit1d(V*b_fit,Y, None, p0=p0, fitfunc=fitfunc, fixed=[],
    	do_print=False, ret=True)
a_fit = fit_result['params_dict']['a']
b_fit = fit_result['params_dict']['b']
c_fit = fit_result['params_dict']['c']
c_err = fit_result['error_dict']['c']
print 'a= ',a_fit
print 'b=',b_fit
print 'c=', c_fit, ' -- error: ', c_err

xxx = np.linspace(X[0], X[-1], 1000)
yyy = a_fit/((xxx-b_fit)**2 + c_fit**2)

f = plt.figure()
ax = f.add_subplot (1,1,1)
ax.plot (X, Y, '.b', linewidth=2)
ax.plot (xxx, yyy, 'r')
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
ax.set_xlabel ('freq [GHz]', fontsize= 14)
ax.set_ylabel ('transmission [a.u.]', fontsize=14)

plt.show()
'''