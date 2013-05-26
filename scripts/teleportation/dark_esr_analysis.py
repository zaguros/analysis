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
guess_offset = 1
guess_ctr = 2.827
guess_splitB = 30.
guess_splitN = 2.187e-3
guess_splitC = .377e-3 #12.78
guess_width = 0.03e-3
guess_amplitude = 0.2


### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('DarkESR')
print folder

a = sequence.SequenceAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results('ssro')
a.get_electron_ROC()


x = a.sweep_pts # convert to MHz
y = a.p0.reshape(-1)[:]
ax = a.plot_result_vs_sweepparam(ret='ax', name='ssro')
ax.set_ylim(0.75,1.05)

# try fitting
fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
        guess_amplitude, guess_width, guess_ctr,
        # (2, guess_splitN), 
        # (2, guess_splitC),
        # (2, guess_splitB),
        (3, guess_splitN), 
        do_print=True, ret=True, fixed=[])
plot.plot_fit1d(fit_result, x, ax=ax, plot_data=False)

ax.set_xlabel('MW frq (GHz)')
ax.set_ylabel(r'uncorrected fidelity $F(|0\rangle)$')
ax.set_title(a.timestamp+'\n'+a.measurementstring)

plt.savefig(os.path.join(folder, 'darkesr_analysis.png'), 
        format='png')






