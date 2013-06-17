import os, sys
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.fitting import fit
from analysis.lib.tools import plot

timestamp = None

g_f = 1./360
g_A = 0.5
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
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
ax = a.plot_results_vs_sweepparam(ret='ax', )

x = a.sweep_pts.reshape(-1)[:] - 1
y = a.p0.reshape(-1)[:]

fit_result = fit.fit1d(x, y, None,
        fitfunc=fitfunc, p0=[A,o,phi], fixed=[], do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(a.sweep_pts[0],a.sweep_pts[-1],201), ax=ax,
        plot_data=False)





