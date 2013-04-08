import os, sys
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.fitting import fit
from analysis.lib.tools import plot

timestamp = '20130303155729' # '20130107231602'

g_f = 1./360
g_A = 0.1
g_o = 0.5
g_phi = 0.

f = fit.Parameter(g_f, 'f')
A = fit.Parameter(g_A, 'A')
o = fit.Parameter(g_o, 'o')
phi = fit.Parameter(g_phi, 'phi')

def fitfunc(x):
    return o() + A() * np.cos(2*pi*(f()*x + phi()/360.))

if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results()
a.get_N_ROC(0.99, 0.02, 0.94, 0.005)
ax = a.plot_results_vs_sweepparam(ret='ax', )

fit_result = fit.fit1d(a.sweep_pts, a.p0.reshape(-1), None,
        fitfunc=fitfunc, p0=[f,A,o,phi], fixed=[], do_print=True, ret=True)

plot.plot_fit1d(fit_result, np.linspace(a.sweep_pts[0],a.sweep_pts[-1],201), ax=ax,
        plot_data=False)





