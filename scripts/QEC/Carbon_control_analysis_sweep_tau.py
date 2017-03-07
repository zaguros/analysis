import os, sys
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi; reload(mbi)
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot
from analysis.lib.math import error

import analysis.scripts.pulse_calibration.calibration_funcs as funcs


def carbon_control_sweep_tau(
timestamp =None,
name = 'Decoupling',#'20150721182946' #20140505090504' #'154236'#'175555' #'160434'
guess_x0 = 6.62,
guess_of = 0.1,
guess_a = 0,
do_fit = False):

	### script
	if timestamp != None:
	    folder = toolbox.data_from_time(timestamp)
	else:
	    folder = toolbox.latest_data(name)


	a = mbi.MBIAnalysis(folder)
	a.get_sweep_pts()
	a.get_readout_results(name='adwindata')
	a.get_electron_ROC()
	ax = a.plot_results_vs_sweepparam(name='adwindata', ret='ax', fmt = '-o')

	x = a.sweep_pts.reshape(-1)[:]
	y = a.p0.reshape(-1)[:]

	print 'minimum'
	print x[np.argmin(y)]

	if do_fit:
		res = funcs.calibrate_pulse_amplitude(x, y, ax, guess_x0, guess_of, guess_a)

	plt.savefig(os.path.join(folder, 'sweep_tau.pdf'),
	        format='pdf')
	plt.savefig(os.path.join(folder, 'sweep_tau.png'),
	        format='png')

# fig = a.default_fig(figsize=(6,4))
# ax = a.default_ax(fig)
# ax.plot(a.sweep_pts, a.p0, '-ob', lw=1)
