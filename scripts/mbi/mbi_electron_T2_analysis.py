import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot
from analysis.lib.math import error

# fit_startup = False

timestamp = None#'20150404123117'# '20150403184301'#one#'20150403165700' #'20130907183620' # None
guess_frq = 0#/200.
guess_amp = 0.5
guess_Tdec = 1/0.18
guess_phi = 180.
guess_o = 1.
plot_total_free_evolution_time = False
### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
print a.sweep_pts
a.get_readout_results(name='adwindata')
a.get_electron_ROC()

if plot_total_free_evolution_time:
    a.sweep_name='Total free evolution time (ms)'
    #print a.g.attrs['Number_of_pulses'][0]
    a.sweep_pts=a.g.attrs['Number_of_pulses'][0]*2*a.sweep_pts*1e-3
ax = a.plot_results_vs_sweepparam(ret='ax', )

t_shift = 0
# t_start=0
# t_0=10

x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]
print x
print y
o = fit.Parameter(guess_o, 'o')
f = fit.Parameter(guess_frq, 'f')
A = fit.Parameter(guess_amp, 'A')
phi = fit.Parameter(guess_phi, 'phi')
Tdec = fit.Parameter(guess_Tdec, 'Tdec')
p0 = [f, A,o,phi,Tdec]
fitfunc_str = ''

def fitfunc(x) : 
    return (o()-A()) + A() * exp(-(x/float(Tdec()))**1) * cos(2*pi*(f()*x - phi()))

fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True,fixed=[0])

#plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
#        plot_data=False)

plt.savefig(os.path.join(folder, 'mbi_erabi_analysis.pdf'),
        format='pdf')
plt.savefig(os.path.join(folder, 'mbi_erabi_analysis.png'),
        format='png')




# print 'The pulse length shift is:' + str(t_shift)

# else:
    # fit_result = fit.fit1d(a.sweep_pts[:], a.p0.reshape(-1)[:], rabi.fit_rabi_fixed_upper,
        # guess_frq, guess_amp, guess_phi, guess_k, fixed=[],
        # do_print=True, ret=True)
    # plot.plot_fit1d(fit_result, np.linspace(0,a.sweep_pts[-1], 201), ax=ax,
        # plot_data=False)

    # plt.savefig(os.path.join(folder, 'electronrabi_analysis.pdf'),
        # format='pdf')

        
        
### FFT
# p0_fft = np.fft.fft(a.p0.reshape(-1), n=32)
# frq = np.fft.fftfreq(32, d=(a.sweep_pts[1]-a.sweep_pts[0])/1e3)
#  
# fig = a.default_fig()
# ax = a.default_ax(fig)
# ax.plot(frq, p0_fft, 'o-')
