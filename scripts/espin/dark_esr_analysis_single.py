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
timestamp = None # '125821' #' #'114103_PulsarD' #YYYYmmddHHMMSS
guess_offset = 1
guess_ctr = 2.8280
guess_splitB = 30.
guess_width = 0.2e-3
guess_amplitude = 0.3
guess_splitN = 181e-6

def analyze_dark_esr_single(timestamp = None,center_guess = False, ax=None, ret=None,min_dip_depth = 0.85 , **kw):

    if ax == None:
        fig, ax = plt.subplots(1,1)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('DarkESR')


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
        guess_ctr = x[y.argmin()]
        print 'guess_ctr = '+str(guess_ctr)
    
    # try fitting
    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            do_print=False, ret=True, fixed=[])
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



def analyze_dark_esr_double(timestamp = None,center_guess = False, ax=None, ret=None,min_dip_depth = 0.85 , **kw):

    if ax == None:
        fig, ax = plt.subplots(1,1)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('DarkESR')


    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()

    x = a.sweep_pts # convert to MHz
    # y = a.get_readout_results('ssro')
    y = a.p0.reshape(-1)[:]
    a.plot_result_vs_sweepparam(ret=ret, name='ssro', ax=ax)
    ax.set_ylim(0.6,1.05)


   # try fitting

    guess_offset = 1.0
    guess_A_min = 0.3
    guess_A_plus = 0.3
    guess_x0 = 1.74666
    guess_sigma = 0.100e-3
    guess_Csplit = 0.210e-3/2

    if center_guess == True:
        guess_x0 = float(raw_input('Center guess?'))
    else:
        guess_x0 = x[y.argmin()]
        guess_x0 = x[len(x)/2.]
        print 'guess_ctr = '+str(guess_x0)
    


    ### fitfunction
    A_min = fit.Parameter(guess_A_min, 'A_min')
    A_plus = fit.Parameter(guess_A_plus, 'A_plus')
    o = fit.Parameter(guess_offset, 'o')
    x0 = fit.Parameter(guess_x0, 'x0')
    sigma = fit.Parameter(guess_sigma, 'sigma')
    Csplit = fit.Parameter(guess_Csplit, 'Csplit')

    def fitfunc(x):
        return o() - A_min()*np.exp(-((x-(x0()-Csplit()))/sigma())**2) \
                - A_plus()*np.exp(-((x-(x0()+Csplit()))/sigma())**2) \


    fit_result = fit.fit1d(x, y, None, p0 = [A_min, A_plus, o, x0, sigma, Csplit],
            fitfunc = fitfunc, do_print=True, ret=True, fixed=[])
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
# if __name__ == '__main__':









