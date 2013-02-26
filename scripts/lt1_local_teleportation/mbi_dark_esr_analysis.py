import os, sys
import numpy as np
import h5py
import logging
import sympy

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import mbi
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, esr
from analysis.lib.tools import plot
from analysis.lib.math import error

### settings
timestamp = None # 

guess_offset = 0.91
guess_A_min1 = 0.
guess_A_plus1 = 0.
guess_A_0 = 0.4
guess_x0 = 2863.08 - 2820
guess_sigma = 0.45

splitting = 2.184


### fitfunction
A_min1 = fit.Parameter(guess_A_min1, 'A_min1')
A_plus1 = fit.Parameter(guess_A_plus1, 'A_plus1')
A_0 = fit.Parameter(guess_A_0, 'A_0')
o = fit.Parameter(guess_offset, 'o')
x0 = fit.Parameter(guess_x0, 'x0')
sigma = fit.Parameter(guess_sigma, 'sigma')

def fitfunc(x):
    return o() - A_min1()*np.exp(-((x-(x0()-splitting))/sigma())**2) \
            - A_plus1()*np.exp(-((x-(x0()+splitting))/sigma())**2) \
            - A_0()*np.exp(-((x-x0())/sigma())**2) 

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('PostInitDarkESR')

a = mbi.PostInitDarkESRAnalysis(folder)
x = a.get_sweep_pts()
y = a.get_readout_results().reshape(-1)

# try fitting
fit_result = fit.fit1d(x, y, None, p0 = [A_min1, A_plus1, A_0, o, x0, sigma],
        fitfunc = fitfunc, do_print=True, ret=True, fixed=[])
ax = plot.plot_fit1d(fit_result, x, ret='ax', plot_data=True)

ax.set_xlabel('MW frq (MHz)')
ax.set_ylabel(r'uncorrected fidelity $F(|0\rangle)$')
ax.set_title(a.timestamp+'\n'+a.measurementstring)

plt.savefig(os.path.join(folder, 'mbi_darkesr_analysis.pdf'), 
        format='pdf')

pol = error.Formula()
a0, am1, ap1 = sympy.symbols('a0, am1, ap1')
pol.formula = a0 / (a0 + ap1 + am1)
pol.values[a0] = A_0()
pol.values[am1] = A_min1()
pol.values[ap1] = A_plus1()
pol.uncertainties[a0] = fit_result['error_dict']['A_0']
pol.uncertainties[am1] = fit_result['error_dict']['A_min1']
pol.uncertainties[ap1] = fit_result['error_dict']['A_plus1']

print 'Spin polarization = %.3f +/- %.3f' \
        % (float(pol.value()), float(pol.uncertainty()))
