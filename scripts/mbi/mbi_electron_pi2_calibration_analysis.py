import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.math import error


timestamp = None#'20141128_130516' #211609'#None#'174154' 
guess_A = [0.5, 0.5]


### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('Ramsey')

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
ax = a.plot_results_vs_sweepparam(ret='ax', name='adwindata' )
ax.set_ylim([0.4,0.6])
fit_result = fit.fit1d(a.sweep_pts, a.p0.reshape(-1), common.fit_poly,
        guess_A,
        do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(0,a.sweep_pts[-1],201), ax=ax,
        plot_data=False)

plt.savefig(os.path.join(folder, 'electronramsey_analysis.pdf'),
        format='pdf')
plt.savefig(os.path.join(folder, 'electronramsey_analysis.png'),
        format='png')

