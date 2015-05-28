import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot
from analysis.lib.math import error
reload(plot)
# fit_startup = False

timestamp ='20150202_162810'
guess_frq = 1/4000.
guess_amp = 1
guess_k = 0.
guess_phi = 0.
guess_o = 1.

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('NuclearRFRabi')


a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
ax = a.plot_results_vs_sweepparam(ret='ax', )

t_shift = 0
# t_start=0
# t_0=10

x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

if x[-1] < x[0]:
    x = x[::-1]
    y = y[::-1]

o = fit.Parameter(guess_o, 'o')
f = fit.Parameter(guess_frq, 'f')
A = fit.Parameter(guess_amp, 'A')
phi = fit.Parameter(guess_phi, 'phi')
k = fit.Parameter(guess_k, 'k')
p0 = [A,o,f]
fitfunc_str = ''

def fitfunc(x) :
    return (o()-A()) + A() * exp(-k()*x) * cos(2*pi*(f()*x - phi()))

fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)

plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, color='0.25',add_txt=False,
        plot_data=False)


ax.set_ylim([-0.05,1.05])
ax.set_xlim([4.688e2,4.6925e2])
plt.xticks(np.linspace(0,2000,5))
plt.yticks(np.linspace(0, 1, 3))
mpl.rcParams['axes.linewidth'] = 2
ax.tick_params(axis='x', which='major', labelsize=15)
ax.tick_params(axis='y', which='major', labelsize=15)
ax.set_xlabel('Pulse Length (us)',fontsize =15)
plt.tight_layout()
ax.tick_params('both', length=4, width=2, which='major')
ax.set_ylabel(r'$F(|0\rangle)$', fontsize=15)
ax.set_title('')

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
