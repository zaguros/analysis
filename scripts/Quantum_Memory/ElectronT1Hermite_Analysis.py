import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import  sequence, mbi #sequence_ssro,
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi, common
reload(rabi)
reload(sequence)

from analysis.lib.tools import plot


def analyze_data(timestamp, output = False):   #,	guess_frq = 1 / 20e9, 	guess_amp = 0.5):

	if timestamp != None:
		folder = toolbox.data_from_time(timestamp)
	else:
		folder = toolbox.latest_data('ElectronT1')

	print folder

	a = sequence.SequenceAnalysis(folder)
	a.get_sweep_pts()
	a.get_readout_results('ssro')
	a.get_electron_ROC()
	
	# print a.p0

	ax = a.plot_result_vs_sweepparam(ret='ax')

	# plt.ylim([0,0.052])
	x = a.sweep_pts
	y = a.p0

	if output:
		return x, y, ax
	else:
		return ax

def fit_exp(x, y,
	guess_a,
	guess_A,
	guess_tau,
	fixed = [],
	show_fit = True,
	save = True,
	timestamp = '',
	**kw):
	'''
	Fit single exponential to T1 data
	'''

	p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(guess_a,guess_A, guess_tau)

	fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=fixed,
		        do_print=True, ret=True)
	
	# Plot fit and format if requested
	keys = sorted(kw.keys())

	if show_fit:
		if len(keys) == 0:
			plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201),
	       plot_data=True)
		elif len(keys) != 0 and 'ax' in keys:
			ax = kw['ax']
			plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201),
	       plot_data=False, ax = ax)
		elif len(keys) != 0 and 'ax' not in keys:
			print 'length of keyword arguments =', len(keys)
			raise Exception("Your keywords for plot formatting didn't contain a pyplot axis with keyword 'ax'. Please provide it.")

	if 'xlabel' in keys:
		ax.set_xlabel(kw['xlabel'])
	if 'ylabel' in keys:
		ax.set_ylabel(kw['ylabel'])


	if timestamp != '':
		if timestamp == None:
			folder = toolbox.latest_data('ElectronT1')
		else:
			folder = toolbox.data_from_time(timestamp)

	if save:
		plt.savefig(os.path.join(folder, 'Calibrate_Pi_analysis_fit.png'))
	return fit_result	