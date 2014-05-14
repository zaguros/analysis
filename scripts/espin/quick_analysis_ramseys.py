
for uniquevar in [var for var in globals().copy() if var[0] != "_" and var != 'clearall']:
	del globals()[uniquevar]

import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt


from analysis.lib import fitting
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, ramsey
from analysis.lib.tools import plot
from analysis.lib.math import error


def bin_data (data = [], bin_size = 50):
	binned_data = []
	N = len(data)/bin_size
	for n in np.arange(N):
		v = np.sum (data[n*bin_size:(n+1)*bin_size-1])/(bin_size+1e-9)
		binned_data.append(v)
	return binned_data

def clean_data (folder = '', sub_string = ''):
	#takes all files in the folder with a given sub-strin gin the filename,
	#bins the data and checks at which point the data becomes zero due to ionization.
	#Only datasets with a number of valid data greater then 'threshold_data' are kept
	
	b_size = 50
	threshold_data = 1200
	thr_ind = threshold_data/b_size
	ind = 0
	all_results = {}
	for d in os.listdir(folder):
		if (d.find (sub_string)>-1):
			print d
			folder0 = folder+d
			a = sequence.SequenceAnalysis(folder0)
			a.get_sweep_pts()
			a.get_readout_results(name='ssro')
			results = a.ssro_results
			binned =  bin_data (data = results, bin_size=b_size)
			binned.append(0)
			binned = np.array(binned)
			indice = np.nonzero(binned<0.001) [0] [0] 
			if (indice>thr_ind):
				results = results [:(indice-1)*b_size]
				all_results [str(ind)] = results
				ind = ind + 1
			a.finish()


	all_results ['max_ind'] = ind
	print 'Nr of good datasets: '+str(ind)
	return all_results

def correlation (results = [], bin_size = 100, N =500):
	wait_time =0.
	print 'N = ',N
	nr_datasets = results ['max_ind']
	ind = 0
	indice_massimo = []
	decay = []
	err_decay = []
	for counter in np.arange(nr_datasets):
		dati = results [str(counter)]
		if (len(dati)>2*N+1):
			if (bin_size>1 ):
				b = bin_data (data = dati, bin_size = bin_size)
			else:
				b = dati
							
			t = np.arange(len(b))*bin_size*(20e-6+wait_time*1e-6)	
			phase = b
			mu = np.mean(phase)
			sigma = np.std(phase)
			corr = np.correlate (phase-mu, phase-mu, 'full')/(np.correlate(phase-mu, phase-mu)+0.)
			t_corr = (np.arange (len(corr))-len(corr)/2.)*(wait_time+20.)*1e-6*bin_size

			nn = len(corr)
			corr2 = corr [nn/2-N:nn/2+N]
			t_corr2 = t_corr [nn/2-N:nn/2+N]
			ind_max = corr2.argmax()
			indice_massimo.append (ind_max)	
			
			if (ind == 0):
				avg_corr = corr2			
			else:
				avg_corr = avg_corr+corr2
			ind = ind + 1
	
	avg_corr[N] = 0
	avg_corr = avg_corr/max(avg_corr)			
	return t_corr2, avg_corr
	
def plot_time_series (results = [], bin_size=1, nr = []):

	nr_datasets = results ['max_ind']
	if (nr==[]):
		nr = np.arange(nr_datasets)
	for counter in nr:
		if (counter<nr_datasets):
			dati = results [str(counter)]
			if (bin_size>1 ):
				b = bin_data (data = dati, bin_size = bin_size)
			else:
				b = data
		
		t = np.arange(len(b))*(20e-3)*bin_size
		plt.figure()
		plt.plot (t, b, 'b')
		plt.plot (t, b, '.r')
		plt.xlabel ('time [msec]')
		plt.title ('dataset nr. '+str(counter)) 
		plt.show()
		
	
	
def fit_exp (x = [], y = [], plot_fits=False):

	guess_b = 10.
	guess_c = 0.
	b = fit.Parameter(guess_b, 'b')

	def fitfunc(x):
		return np.exp(-np.abs((x)/b()))

	xx = np.linspace(0, x[-1], 1000)
	fit_result = fit.fit1d(x,y, None, p0=[b], fitfunc=fitfunc, do_print=False, ret=True)

	b_fit = fit_result ['params_dict'] ['b']
	err_b = fit_result ['error_dict'] ['b']	

	if plot_fits:
		plt.figure()
		plt.plot (x*1000, y, '.b')
		plt.plot (xx*1000, np.exp(-np.abs((xx)/b_fit)), 'r')
		plt.xlabel ('time [msec]')
		plt.show()

	return b_fit, err_b


def plot_data_28april ():
	bin_size = 1
	N = 200/bin_size
	f0 = r'/home/cristian/Work/Research/adaptive magnetometry/dati/20140430/'
	#save_folder = r'/home/cristian/Work/Research/adaptive magnetometry/dati/20140403_analisi/spin_flips/'

	stringa1 = '_45V'
	results1 = clean_data (folder = f0, sub_string = stringa1)
	print "Calculating correlations..."
	t, corr45 = plot_correlation (results= results1, b_size = bin_size, N=N)

	stringa1 = '_35V'
	results1 = clean_data (folder = f0, sub_string = stringa1)
	print "Calculating correlations..."
	t, corr35 = plot_correlation (results= results1, b_size = bin_size, N=N)


	stringa1 = '_28V'
	results1 = clean_data (folder = f0, sub_string = stringa1)
	print "Calculating correlations..."
	t, corr28 = plot_correlation (results= results1, b_size = bin_size, N=N)


	stringa1 = '_30V'
	results1 = clean_data (folder = f0, sub_string = stringa1)
	print "Calculating correlations..."
	t, corr30 = plot_correlation (results= results1, b_size = bin_size, N=N)


	stringa1 = '_20V'
	results1 = clean_data (folder = f0, sub_string = stringa1)
	print "Calculating correlations..."
	t, corr20 = plot_correlation (results= results1, b_size = bin_size, N=N)


	stringa1 = '_25V'
	results1 = clean_data (folder = f0, sub_string = stringa1)
	print "Calculating correlations..."
	t, corr25 = plot_correlation (results= results1, b_size = bin_size, N=N)


	stringa1 = 'Rep'
	results1 = clean_data (folder = f0+'data_0V/', sub_string = stringa1)
	print "Calculating correlations..."
	t, corr0 = plot_correlation (results= results1, b_size = bin_size, N=N)

	stringa1 = '_18V'
	results1 = clean_data (folder = f0, sub_string = stringa1)
	print "Calculating correlations..."
	t, corr18 = plot_correlation (results= results1, b_size = bin_size, N=N)

	stringa1 = '_5V'
	results1 = clean_data (folder = f0, sub_string = stringa1)
	print "Calculating correlations..."
	t, corr5 = plot_correlation (results= results1, b_size = bin_size, N=N)

	stringa1 = '_7V'
	results1 = clean_data (folder = f0, sub_string = stringa1)
	print "Calculating correlations..."
	t, corr7 = plot_correlation (results= results1, b_size = bin_size, N=N)

	stringa1 = '_10V'
	results1 = clean_data (folder = f0, sub_string = stringa1)
	print "Calculating correlations..."
	t, corr10 = plot_correlation (results= results1, b_size = bin_size, N=N)

	stringa1 = '_36V'
	results1 = clean_data (folder = f0, sub_string = stringa1)
	print "Calculating correlations..."
	t, corr36 = plot_correlation (results= results1, b_size = bin_size, N=N)


	t = np.arange(len(corr0))-N
	plt.plot (t, corr0, '.b')
	plt.plot (t, corr5, '.r')
	plt.plot (t, corr25, '.k')
	plt.plot (t, corr36, '.g')
	plt.xlabel ('nr of Ramseys')
	plt.show()

	i = N/2
	b0, e0 = fit_exp(x = t[:], y=corr0[:], plot_fits=True)
	b5, e5 = fit_exp(x = t[:], y=corr5[:], plot_fits=True)
	b7, e7 = fit_exp(x = t[:], y=corr7[:], plot_fits=False)
	b10, e10 = fit_exp(x = t[:], y=corr10[:], plot_fits=False)
	b25, e25 = fit_exp(x = t[:], y=corr25[:], plot_fits=False)
	b18, e18 = fit_exp(x = t[:], y=corr18[:], plot_fits=False)
	b20, e20 = fit_exp(x = t[:], y=corr20[:], plot_fits=False)
	b36, e36 = fit_exp(x = t[:], y=corr36[:], plot_fits=False)
	b30, e30 = fit_exp(x = t[:], y=corr30[:], plot_fits=False)
	b28, e28 = fit_exp(x = t[:], y=corr28[:], plot_fits=False)
	b35, e35 = fit_exp(x = t[:], y=corr35[:], plot_fits=False)
	b45, e45 = fit_exp(x = t[:], y=corr45[:], plot_fits=False)

	bb = np.array([0,5,7,10,25, 20, 30, 28, 35])
	aa = np.array([11.1, 12.6,13.17, 13.,13.02, 12.82, 13.56, 13.34, 13.74, 13.63, 13.84, 14.34])
	b_fits = np.abs(np.array([b36, b18, b0,b5,b7,b10,b25, b20, b30, b28, b35, b45]))
	e_fits = np.abs(np.array([e36, e18,e0,e5,e7,e10,e25, e20, e30, e28, e35, e45]))

	plt.figure()
	plt.errorbar (aa, b_fits, yerr=e_fits, fmt='ob')
	plt.xlabel ('Freq diff E''-Ex [GHz]')
	plt.show()


f0 = r'Z:/Diamond/Projects/Magnetometry with adaptive measurements/Data/20140505/'
data = clean_data (folder=f0, sub_string = '45_hf_')
plot_time_series (results = data, bin_size=25)
t, corr = correlation (results=data, bin_size=1)

x = np.arange (len(corr))-len(corr)/2
plt.plot (x, corr, linewidth =2)
plt.xlim ([0, 150])
plt.ylim ([0,1])
plt.xlabel ('nr of ramseys')
plt.show()


