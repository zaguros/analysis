import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import  sequence, mbi #sequence_ssro,
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
reload(rabi)
reload(sequence)

from analysis.lib.tools import plot


def analyze_data(timestamp, output = False):   #,	guess_frq = 1 / 20e9, 	guess_amp = 0.5):

	# o = fit.Parameter(guess_of, 'o')
	# f = fit.Parameter(guess_frq, 'f')
	# A = fit.Parameter(guess_amp, 'A')
	# phi = fit.Parameter(guess_phi, 'phi')
	# k = fit.Parameter(guess_k, 'k')
	# p0 = [f, A, phi, o, k]
	# fitfunc_str = ''

	if timestamp != None:
		folder = toolbox.data_from_time(timestamp)
	else:
		folder = toolbox.latest_data('Pi_Calibration')

	print folder


	a = sequence.SequenceAnalysis(folder)
	a.get_sweep_pts()
	a.get_readout_results('ssro')
	a.get_electron_ROC()
	


	ax = a.plot_result_vs_sweepparam(ret='ax')

	# plt.ylim([0,0.052])
	x = a.sweep_pts
	y = a.p0

	if output:
		return folder, x, y

def fit_sine(x, y, 
	guess_frq, 
	guess_amp,
	guess_phi,
	guess_of,
	guess_k,
	fixed,
	folder,
	show_fit = True,
	save = True,):
	

	o = fit.Parameter(guess_of, 'o')
	f = fit.Parameter(guess_frq, 'f')
	A = fit.Parameter(guess_amp, 'A')
	phi = fit.Parameter(guess_phi, 'phi')
	k = fit.Parameter(guess_k, 'k')
	p0 = [f, A, phi, o, k]
	fitfunc_str = ''

	fitfunc_str = 'o - A + A*e^(-(kx)**2)*cos(2pi (fx-phi))'

	def fitfunc(x):
	    return (o()-A()) + A() * np.exp(-(k()*x)**2) * np.cos(2*np.pi*(f()*x - phi()))

	fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=fixed,
	        do_print=True, ret=True)
	
	if show_fit:
		plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201),
	       plot_data=True)
	# print "pi pulse = {:.5f} ".format(1/f()/2.) + a.sweep_name

	# ax.set_title(a.timestamp+'\n'+a.measurementstring)
	if save:

		plt.savefig(os.path.join(folder, 'Calibrate_Pi_analysis_fit.png'))

def fit_parabola(x, y, 
	guess_x0, 
	guess_amp,
	guess_of,
	fixed,
	show_fit = True,
	save = True):
	
	x0 = fit.Parameter(guess_x0, 'x0')
	A = fit.Parameter(guess_amp, 'A')
	o = fit.Parameter(guess_of, 'o')
	p0 = [x0, A, o]
	fitfunc_str = ''

	fitfunc_str = 'o + A*(x - x0)^2'

	def fitfunc(x):
	    return o() + A() * (x - x0) ** 2.

	fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=fixed,
	        do_print=True, ret=True)

	if show_fit:
		plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201),
	       plot_data=True)
	# print "pi pulse = {:.5f} ".format(1/f()/2.) + a.sweep_name

	# ax.set_title(a.timestamp+'\n'+a.measurementstring)
	if save:
		plt.savefig(os.path.join(folder, 'Calibrate_Pi_analysis_parabolafit.png'))