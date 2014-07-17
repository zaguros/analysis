import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import  sequence, mbi #sequence_ssro,
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit
from analysis.lib.tools import plot

timestamp = None

if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('Pi')


a = sequence.SequenceAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results('ssro')
a.get_electron_ROC()    
x = a.sweep_pts
y = a.p0
u_y = a.u_p0
n = a.sweep_name
a.finish()

x2 = x[::2]
print x2
y2 = y[1::2] - y[::2]
u_y2 = np.sqrt(  u_y[1::2]**2 + u_y[::2]**2 )    

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,4), sharex=True)
ax1.errorbar(x2, y2, yerr=u_y2, fmt='o')
ax1.set_xlabel(n)
ax1.set_title('Difference btw. Pi/2-Pi and Pi/2')
ax1.set_ylabel('Difference')

m = fit.Parameter((y[-1]-y[0])/(x[-1]-x[0]), 'm')
x0 = fit.Parameter(x2.mean(), 'x0')
p0 = [m, x0]

def ff(x):
    return m() * (x-x0())
fitfunc_str = 'm * (x - x0)'

fit_result = fit.fit1d(x2, y2, None, p0=p0, fitfunc=ff,
    fitfunc_str=fitfunc_str, do_print=True, ret=True)    
    
ax2.errorbar(x2, y[0::2], yerr=u_y[0::2], fmt='o',
             label='Pi/2')
ax2.errorbar(x2, y[1::2], yerr=u_y[1::2], fmt='o',
             label='Pi/2 - Pi')
ax2.legend(frameon=True, framealpha=0.5)
ax2.set_ylabel('P(0)')
ax2.set_xlabel(n)

if fit_result != False  :
    plot.plot_fit1d(fit_result, np.linspace(x2[0],x2[-1],201), ax=ax1,
        plot_data=False, print_info=True)
    if  a.sweep_pts[0] < x0() < a.sweep_pts[-1] :
        ax2.axvline(x0(), c='k', lw=2)
        ax2.axhline(0.5, c='k', lw=2)
        ax2.set_title('X marks the spot')

fig.savefig(os.path.join(folder, 'pi2_calibration.png'))

