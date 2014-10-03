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

fName = 
y = np.loadtxt(fname='D:/measuring/data/20140924/scope_trace.txt')
x_ns = 1e9*np.arange(len(a))*1.25e-12

guess_frq = 1./15000
guess_amp = 0.2
guess_of = 0.1
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
plt.show()

'''
fitfunc_str = 'o - A + A*e^(-kx)*cos(2pi (fx-phi))'

def fitfunc(x):
    return (o()-A()) + A() * np.exp(-k()*x) * np.cos(2*np.pi*(f()*x - phi()))

fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=[2],
        do_print=True, ret=True)

plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax,
        plot_data=False)

print "pi pulse = {:.2f} ".format(1/f()/2.) + a.sweep_name

# ax.set_title(a.timestamp+'\n'+a.measurementstring)
plt.savefig(os.path.join(folder, 'electronrabi_analysis_fit.png'))
'''
