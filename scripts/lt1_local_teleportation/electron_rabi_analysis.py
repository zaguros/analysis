import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.ssro_m2 import sequence_ssro
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot

timestamp = None # '20130107231602'
guess_frq = 1./4.
guess_amp = 0.05
guess_yof = 0.83
guess_tau = 5
guess_slope = 0.
guess_xof = 0.

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('Rabi')

a = sequence_ssro.ElectronRabiAnalysis(folder)
x = a.get_sweep_pts() / 1e3
y = a.get_ssro_results()

fit_result = fit.fit1d(x, y, rabi.fit_rabi_damped_exp_with_offset_on_linslope,
        guess_frq, guess_amp, guess_yof, guess_tau, guess_xof, guess_slope,
        fixed=[4,5],
        do_print=True, ret=True)
ax = plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ret='ax')

ax.set_xlabel(r'MW pulse length ($\mu$s)')
ax.set_ylabel(r'uncorrected fidelity $F(|0\rangle)$')
ax.set_title(a.timestamp+'\n'+a.measurementstring)

plt.savefig(os.path.join(folder, 'electronrabi_analysis.pdf'), 
        format='pdf')
