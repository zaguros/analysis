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

from analysis.lib.tools import plot

y = np.loadtxt(fname='D:/measuring/data/20140924/scope_trace.txt')
folder = r'D:/measuring/data/20140924/scope_trace/'
y = y[1000:]
x_ns = 1e9*np.arange(len(y))*8e-12

guess_frq = 1./30
guess_amp = 0.25
guess_of = 0.0
# guess_slope = 0.
guess_phi = 0.
guess_k = 0.

mbi_analysis = False

o = fit.Parameter(guess_of, 'o')
f = fit.Parameter(guess_frq, 'f')
A = fit.Parameter(guess_amp, 'A')
phi = fit.Parameter(guess_phi, 'phi')
k = fit.Parameter(guess_k, 'k')
p0 = [f, A, phi, o, k]
fitfunc_str = ''

plt.plot (x_ns, y)
plt.xlabel ('ns')


fitfunc_str = 'o - A + A*e^(-kx)*cos(2pi (fx-phi))'

def fitfunc(x):
    return o() + A() * np.exp(-k()*x) * np.cos(2*np.pi*(f()*x - phi()))

fit_result = fit.fit1d(x_ns,y, None, p0=p0, fitfunc=fitfunc, fixed=[4],do_print=True, ret=True)

plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax,
        plot_data=False)

#print "pi pulse = {:.2f} ".format(1/f()/2.) + a.sweep_name
A_fit = fit_result['params_dict']['A']
o_fit = fit_result['params_dict']['o']
f_fit = fit_result['params_dict']['f']
phi_fit = fit_result['params_dict']['phi']

xx = np.linspace(0, 250, 10000)
yy = o_fit +A_fit * np.cos(2*np.pi*(f_fit*xx - phi_fit))

plt.plot (xx,yy,'r')
plt.show()
# ax.set_title(a.timestamp+'\n'+a.measurementstring)
plt.savefig(os.path.join(folder, 'oscilloscope_trace_fit.png'))

