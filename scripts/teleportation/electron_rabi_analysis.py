import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import sequence_ssro
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
reload(rabi)

from analysis.lib.tools import plot

timestamp = None # '20130107231602'
guess_frq = 1./300.
guess_amp = 1
guess_yof = 0.6
guess_tau = 10000
guess_slope = 0.
guess_xof = 0.

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('Rabi')

a = sequence.SequenceAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results('ssro')
a.get_electron_ROC()

x = a.sweep_pts
y = a.p0

fit_result = fit.fit1d(x, y, rabi.fit_rabi_multiple_detunings,
        guess_amp, guess_yof, guess_frq, guess_tau, (0, 0), (-2.193e-3, 0), (2.193e-3, 0), fixed=[],
        do_print=True, ret=True)
ax = plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ret='ax')

ax.set_xlabel(r'MW pulse length (ns)')
ax.set_ylabel(r'uncorrected fidelity $F(|0\rangle)$')
ax.set_title(a.timestamp+'\n'+a.measurementstring)

plt.savefig(os.path.join(folder, 'electronrabi_analysis.png'))

# (-2.193e-3, 0), (2.193e-3, 0)
