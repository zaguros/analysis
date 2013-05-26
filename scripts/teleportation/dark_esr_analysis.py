import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import sequence
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot

### settings
timestamp = None # 
guess_offset = 0.78
guess_ctr = 2830.13
guess_splitB = 30.
guess_splitN = 2.189
guess_splitC = .377 #12.78
guess_width = 0.38
guess_amplitude = 0.05


### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('DarkESR')
print folder

a = sequence.SequenceAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results('ssro')

x = a.sweep_pts*1e3 # convert to MHz
y = a.normalized_ssro

# try fitting
fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
        guess_amplitude, guess_width, guess_ctr,
        # (2, guess_splitN), 
        # (2, guess_splitC),
        # (2, guess_splitB), 
        (3, guess_splitN),
        do_print=True, ret=True, fixed=[4])
ax = plot.plot_fit1d(fit_result, x, ret='ax', plot_data=True)

ax.set_xlabel('MW frq (MHz)')
ax.set_ylabel(r'uncorrected fidelity $F(|0\rangle)$')
ax.set_title(a.timestamp+'\n'+a.measurementstring)

plt.savefig(os.path.join(folder, 'darkesr_analysis.png'), 
        format='png')






