import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import  sequence, mbi #sequence_ssro,
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
reload(rabi)
reload(sequence)

from analysis.lib.tools import plot

pulse_shape = 'Square'


timestamp = None #'20150407144111' #'20150403163852'#None#'20140408125318'
guess_frq = 1./4000.

guess_amp = 0.5
guess_of = 0.5
# guess_slope = 0.
guess_phi = 0.
guess_k = 1/2000.
fixed = [2]

mbi_analysis = False

# ylim = (-0.05, 0.1)
ylim = (-.05, 1.05)


o = fit.Parameter(guess_of, 'o')
f = fit.Parameter(guess_frq, 'f')
A = fit.Parameter(guess_amp, 'A')
phi = fit.Parameter(guess_phi, 'phi')
k = fit.Parameter(guess_k, 'k')
p0 = [f, A, phi, o, k]
fitfunc_str = ''

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    if pulse_shape == 'Hermite':
        folder = toolbox.latest_data('ElectronRabiHermite')
    elif pulse_shape == 'Square':
        folder = toolbox.latest_data('ElectronRabi')

if mbi_analysis:
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('adwindata')
    a.get_electron_ROC()
    ax = a.plot_results_vs_sweepparam(ret='ax', name = 'adwindata' )

else:
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()
    ax = a.plot_result_vs_sweepparam(ret='ax', ylim = ylim)

x = a.sweep_pts
y = a.p0


# fit_result = fit.fit1d(x, y, rabi.fit_rabi_multiple_detunings,
#         guess_amp, guess_yof, guess_frq, guess_tau, (0, 0), (-2.193e-3, 0), (2.193e-3, 0), fixed=[],
#         do_print=True, ret=True)

fitfunc_str = 'o - A + A*e^(-(kx)**2)*cos(2pi (fx-phi))'

def fitfunc(x):
    return (o()-A()) + A() * np.exp(-(k()*x)**2) * np.cos(2*np.pi*(f()*x - phi()))

fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=fixed,
        do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax,
       plot_data=False)
ax.set_ylim([0.,1.05])
print "pi pulse = {:.5f} ".format(1/f()/2.) + a.sweep_name

# ax.set_title(a.timestamp+'\n'+a.measurementstring)
plt.savefig(os.path.join(folder, 'electronrabi_analysis_fit.png'))

