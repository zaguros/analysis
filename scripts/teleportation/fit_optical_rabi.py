import os, sys
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot

from analysis.lib.tools import toolbox

timestamp = None #'20130801132824'

x_i = 0
x_f = 1200-1
pts = 1200

g_f = 1./360.#8.593e-3#2.19290
g_A = 0.5
g_o = 0.5
g_x0 = 0
g_tau = 200.

f = fit.Parameter(g_f, 'f')
A = fit.Parameter(g_A, 'A')
o = fit.Parameter(g_o, 'o')
x0 = fit.Parameter(g_x0, 'x0')
tau = fit.Parameter(g_tau, 'c')


def fitfunc(x):
    return o() + A() * np.cos(2*pi*(f()*(x - x0()))) 

if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()

datafile = open('Optical_Rabi.dat', 'r')
data = []
for row in datafile:
    data.append(row.strip().split(','))

fclose(datafile)



x = np.linspace(x_i,x_f,pts)
y = data[:,1]

fit_result = fit.fit1d(x, y, rabi.fit_rabi_damped_exp_with_offset, 
		g_f, g_A, g_o, g_tau, g_x0, fixed=[], do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False)

ax.title
