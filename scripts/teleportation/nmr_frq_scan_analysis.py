import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib.m2.ssro import ssro, mbi
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot

timestamp = None # '20130107231602'
guess_a = 1
guess_A = -0.8
guess_F = 0.005
guess_x0 = 7.138

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('NMR')

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results()
a.get_N_ROC(0.99, 0.02, 0.94, 0.005)
ax = a.plot_results_vs_sweepparam(ret='ax', )

fit_result=fit.fit1d(a.sweep_pts, a.p0.reshape(-1), 
        rabi.fit_population_vs_detuning,
        guess_a, guess_A, guess_F, guess_x0,
        fixed=[0], do_print=True, ret=True)

plot.plot_fit1d(fit_result, np.linspace(a.sweep_pts[0],a.sweep_pts[-1],201), 
        ax=ax, plot_data=False)

plt.savefig(os.path.join(folder, 'nmr_analysis.pdf'),
        format='pdf')
