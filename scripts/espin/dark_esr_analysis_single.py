import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot

### settings
timestamp = None# '125821' #' #'114103_PulsarD' #YYYYmmddHHMMSS
guess_offset = 1
guess_ctr = 2.8280
guess_splitB = 30.
guess_width = 0.2e-3
guess_amplitude = 0.3

def analyze_dark_esr(folder,center_guess = False, ax=None, ret=None,min_dip_depth = 0.85 , **kw):

    if ax == None:
        fig, ax = plt.subplots(1,1)

    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()

    x = a.sweep_pts # convert to MHz
    y = a.p0.reshape(-1)[:]
    a.plot_result_vs_sweepparam(ret=ret, name='ssro', ax=ax)
    ax.set_ylim(0.6,1.05)


    if center_guess == True:
        guess_ctr = float(raw_input('Center guess?'))
    else:
        j=0
        while y[j]>min_dip_depth and j < len(y)-2:  #y[j]>0.93*y[j+1]: # such that we account for noise
            k = j
            j += 1
        #j = len(y)-2
        if k > len(y)-3:
            print 'Could not find dip'
            return
        else:
            guess_ctr = x[k] #convert to GHz and go to middle dip
            print 'guess_ctr= '+str(guess_ctr)

    # try fitting
    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            do_print=True, ret=True, fixed=[])
    plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), ax=ax, plot_data=False, **kw)

    ax.set_xlabel('MW frq (GHz)')
    ax.set_ylabel(r'fidelity wrt. $|0\rangle$')
    ax.set_title(a.timestamp+'\n'+a.measurementstring)

    plt.savefig(os.path.join(folder, 'darkesr_analysis.png'),
            format='png')
    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']

        ax.text(f0, 0.8, '$f_0$ = ({:.3f} +/- {:.3f})'.format(
            (f0-2.8)*1e3, u_f0*1e3), ha='center')

        return (f0-2.8)*1e3, u_f0*1e3

### script
if __name__ == '__main__':
    if timestamp != None:
        folder = toolbox.latest_data(timestamp)
    else:
        folder = toolbox.latest_data('DarkESR')

    analyze_dark_esr(folder)







