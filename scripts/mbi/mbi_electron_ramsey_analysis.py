import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, ramsey
from analysis.lib.tools import plot
from analysis.lib.math import error


timestamp = '20141203_180041'#'211609'#None#'174154' 
guess_f1 = 0.6e-3 #in GHz
guess_A1 = 0.5
guess_phi1 = 0.

guess_f2 =0.40e-3 #in GHz
guess_A2 = 0.5
guess_phi2 = 0.

guess_f3 = 2.188e-3 #in GHz
guess_A3 = 0.05
guess_phi3= 180.

guess_tau = 4600
guess_a = 0.5


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
ax.set_ylim([-0.4,1.05])
fit_result = fit.fit1d(a.sweep_pts, a.p0.reshape(-1), ramsey.fit_ramsey_gaussian_decay,
        guess_tau, guess_a, (guess_f1, guess_A1, guess_phi1),
        (guess_f2, guess_A2, guess_phi2),
        #(guess_f3, guess_A3, guess_phi3),
         fixed=[1,4,7],
        do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(0,a.sweep_pts[-1],201), ax=ax,
        plot_data=False)

plt.savefig(os.path.join(folder, 'electronramsey_analysis.pdf'),
        format='pdf')
plt.savefig(os.path.join(folder, 'electronramsey_analysis.png'),
        format='png')


### FFT
# p0_fft = np.fft.fft(a.p0.reshape(-1), n=32)
# frq = np.fft.fftfreq(32, d=(a.sweep_pts[1]-a.sweep_pts[0])/1e3)
#  
# fig = a.default_fig()
# ax = a.default_ax(fig)
# ax.plot(frq, p0_fft, 'o-')
