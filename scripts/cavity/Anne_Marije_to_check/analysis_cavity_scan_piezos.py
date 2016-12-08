
import h5py
from analysis.lib.tools import toolbox
import pylab as plt
from analysis.lib import fitting
from analysis.lib.fitting import fit,esr, common
reload(toolbox)


def get_files (newer_than, older_than, contains):
	a, b, folder_names = toolbox.latest_data (contains=contains, newer_than=newer_than, older_than=older_than, return_all=True)
	times = []
	for f in folder_names:
		times.append (int(f[9:11])*3600 + int(f[11:13])*60 + int(f[13:15]))
	return folder_names, times



def load_data (timestamp):
	folder_name = toolbox.latest_data (contains=timestamp)
	all_files = [f for f in os.listdir(folder_name) if timestamp in f]  
	file_name = [f for f in all_files if '.hdf5' in f]
	
	f = h5py.File(os.path.join(folder_name, file_name[0]),'r')
	ls_grp = f['/laserscan']
	V = ls_grp['V'].value
	PD_signal = ls_grp['PD_signal'].value
	return V, PD_signal

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

def combine_piezoscans_2D (timestamp_array, y_axis = None):

	voltages = np.array([])
	PD_signal = np.array([])
	i = 0
	for k in timestamp_array:
		V, Y = load_data (k)
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
	ax.set_xlabel ('piezo voltage [V]', fontsize= 16)
	ax.set_ylabel ('time [s]', fontsize=16)

	plt.show()


def combine_piezoscans_1D (timestamp_array):

	colori = cm.gist_earth(np.linspace(0,0.75, len(timestamp_array)))
	f = plt.figure(figsize=(12,8))
	ax = f.add_subplot(1,1,1)
	i = 0
	for k in timestamp_array:
		print k
		V, Y = load_data (k)
		ax.plot (V[5:500], Y[5:500]+i*0.1, '.', color=colori[i])
		ax.plot (V[5:500], Y[5:500]+i*0.1, color=colori[i])
		i= i+1
	plt.show()

t_array, times = get_files (contains='laser', older_than = '20150422_215600', newer_than = '20150422_215138')

#times = np.array(times)
#times = times - min(times)
#combine_piezoscans_2D (timestamp_array = t_array, y_axis = times)



V, Y = load_data (timestamp = '112004')
#Vf, freq = load_calibration ()
#freq = freq-min(freq)

#a_fit, b_fit = fit_calibration (V=Vf, freq=freq)
#x_fit = np.linspace(V[0], V[-1], 1000)
#y_fit = a_fit + b_fit*x_fit

#plt.figure()
#plt.plot (Vf, freq, 'b', linewidth=2)
#plt.plot (x_fit, y_fit, 'r')
#plt.show()

#X = V*b_fit
guess_c = 0.002
guess_b = 1.145
guess_a = 50
a = fit.Parameter(guess_a, 'a')
b = fit.Parameter(guess_b, 'b')
c = fit.Parameter(guess_c, 'c')

p0 = [a, b, c]
fitfunc_str = ''

def fitfunc(x):
	return a()/((x-b())**2 + c()**2)


fit_result = fit.fit1d(V,Y, None, p0=p0, fitfunc=fitfunc, fixed=[],
    	do_print=False, ret=True)
a_fit = fit_result['params_dict']['a']
b_fit = fit_result['params_dict']['b']
c_fit = fit_result['params_dict']['c']
c_err = fit_result['error_dict']['c']
print 'a= ',a_fit
print 'b=',b_fit
print 'c=', c_fit, ' -- error: ', c_err

xxx = np.linspace(V[0], V[-1], 1000)
yyy = a_fit/((xxx-b_fit)**2 + c_fit**2)

f = plt.figure()
ax = f.add_subplot (1,1,1)
ax.plot (V, Y, '.b', linewidth=2)
ax.plot (xxx, yyy, 'r')
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
ax.set_xlabel ('Voltage [V]', fontsize= 14)
ax.set_ylabel ('transmission [a.u.]', fontsize=14)

plt.show()
