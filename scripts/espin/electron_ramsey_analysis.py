import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, ramsey
from analysis.lib.tools import plot
from analysis.lib.math import error

N_HF_splt = 2.17e-3
detuning = 3e-3

timestamp = None#20130530183305'#None #
guess_f1 = 0.8e-3#detuning - N_HF_splt
guess_A1 = -0.15
guess_phi1 = 0.
guess_f2 = 3e-3 #detuning
guess_A2 = -0.15
guess_phi2 = 0.
guess_f3 = 5.2e-3 #detuning + N_HF_splt
guess_A3 = -0.15
guess_phi3 = 0.


guess_tau = 6000
guess_a = 0.5


### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('ElectronRamsey')

a = sequence.SequenceAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='ssro')
a.get_electron_ROC()
ax = a.plot_result_vs_sweepparam(ret='ax', name='ssro')

p0, fitfunc_0, fitfunc_str_0 = ramsey.fit_ramsey_gaussian_decay(guess_tau, guess_a, (guess_f1, guess_A1, guess_phi1), (guess_f2,guess_A2, guess_phi2), (guess_f3,guess_A3, guess_phi3))
fit_xvals_0=np.linspace(0,a.sweep_pts[-1],1000)
ax.plot(fit_xvals_0,fitfunc_0(fit_xvals_0), 'r--', lw=1)
#fit_xvals=np.linspace(res['x'][0],res['x'][-1],fit_num_points)


fit_result = fit.fit1d(a.sweep_pts, a.p0.reshape(-1), ramsey.fit_ramsey_gaussian_decay,
        guess_tau, guess_a, (guess_f1, guess_A1, guess_phi1), (guess_f2,guess_A2, guess_phi2), (guess_f3,guess_A3, guess_phi3), fixed=[1,4,7,10],
        do_print=True, ret=True)
#fit_result = False
if fit_result != False :
	plot.plot_fit1d(fit_result, np.linspace(0,a.sweep_pts[-1],201), ax=ax,
        plot_data=False)



plt.savefig(os.path.join(folder, 'electronramsey_analysis.pdf'),
        format='pdf')




### FFT
# p0_fft = np.fft.fft(a.p0.reshape(-1), n=32)
# frq = np.fft.fftfreq(32, d=(a.sweep_pts[1]-a.sweep_pts[0])/1e3)
#
# fig = a.default_fig()
# ax = a.default_ax(fig)
# ax.plot(frq, p0_fft, 'o-')
