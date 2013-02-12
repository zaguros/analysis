import numpy as np
from matplotlib import pyplot as plt
from analysis.lib.m2.ssro import mbi
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot


timestamp = '20130208172813' # None
guess_frq = 1./2.8
guess_amp = 0.15
guess_yof = 0.7
guess_tau = 10
guess_xof = 0.

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('MBIElectronRabi')

a = mbi.ElectronRabiAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results()

x = a.sweep_pts
y = a.normalized_ssro[:,0]

fit_result = fit.fit1d(x, y, rabi.fit_rabi_damped_exp_with_offset,
        guess_frq, guess_amp, guess_yof, guess_tau, guess_xof, fixed=[4],
        do_print=True, ret=True)

if fit_result == None:
    print 'cannot do fit'

elif fit_result['success'] == False:
    print 'cannot do fit'

else:
    
    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111)
    ax.errorbar(x,y, fmt='o', yerr=a.u_normalized_ssro[:,0])

    plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax)

    ax.set_xlabel(r'MW pulse length ($\mu$s)')
    ax.set_ylabel(r'uncorrected fidelity $F(|0\rangle)$')
    ax.set_title(a.timestamp+'\n'+a.measurementstring)
