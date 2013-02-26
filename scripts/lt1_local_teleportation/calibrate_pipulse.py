import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import mbi

if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('calib')

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results()
a.get_electron_ROC()
ax = a.plot_results_vs_sweepparam(ret='ax', )

x = a.sweep_pts.reshape(-1)
y = a.p0.reshape(-1)

x0 = fit.Parameter(0.15, 'x0')
of = fit.Parameter(1., 'of')
a = fit.Parameter(0., 'a')
fitfunc_str = '(1-of) + a (x-x0)**2'

def fitfunc(x):
    return (1.-of()) + a() * (x-x0())**2

fit_result = fit.fit1d(x,y, None, p0=[x0,of,a], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False)


