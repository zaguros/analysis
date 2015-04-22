import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import sequence, mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot

reload(plot)
### settings
timestamp = '20150203_144510' # '125821' #' #'114103_PulsarD' #YYYYmmddHHMMSS
guess_offset = 0.9
#guess_ctr = 469.023e3 #4.08306e8
#guess_splitB = 30.
guess_width = 0.5
guess_amplitude = 0.8
# guess_splitN = 181e-6

def analyze_nmr_single(timestamp = None,center_guess = False, ax=None, ret=None,min_dip_depth = 0.85 , **kw):

    if ax == None:
        fig, ax = plt.subplots(1,1)
        
    x = np.zeros((0))
    y = np.zeros((0))
    ysig = np.zeros((0))
    # rdts = np.zeros((0))


    for i in range(15):
        timestamp, folder = toolbox.latest_data(contains = 'NuclearRFRabi_111_1_sil18Rabi_C5_el1_positive_'+str(i),return_timestamp = True)
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC()    
        x = np.append(x,a.sweep_pts[:]) 
        y = np.append(y,a.p0[:])
        ysig = np.append(ysig,a.u_p0[:])
    print dir(a)

    for i in range(len(x)-1):
        print x[i], x[i+1]-x[i]



    if center_guess == True:
        guess_ctr = float(raw_input('Center guess?'))
    else:
        guess_ctr = x[y.argmin()]
        print 'guess_ctr = '+str(guess_ctr)
    
    # try fitting
    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            do_print=False, ret=True, fixed=[])
    fit_result['yerr'] = ysig
    plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), ax=ax, plot_data=True, linewidth=2, color='0.25', figsize=(6,4.7), add_txt=False,**kw)
    plt.tight_layout()
    # plt.xticks(np.linspace(468.9,469.2,4),[str(468.9),str(469.0),str(469.1),str(469.2)])
    # ax.set_ylim([0.5,1])
    # ax.set_xlim([4.6885e2,4.6921e2])
    # plt.yticks([0.5, 0.75, 1])
    mpl.rcParams['axes.linewidth'] = 2
    ax.tick_params(axis='x', which='major', labelsize=15)
    ax.tick_params(axis='y', which='major', labelsize=15)
    ax.set_xlabel('RF frq (kHz)',fontsize =15)
    
    ax.tick_params('both', length=4, width=2, which='major')
    ax.set_ylabel(r'fidelity wrt. $|0\rangle$', fontsize=15)
    # ax.set_title(a.timestamp+'\n'+a.measurementstring)

    plt.savefig(os.path.join(folder, 'nmr_analysis.pdf'),
            format='pdf')
    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']

        ax.text(f0, 0.8, '$f_0$ = ({:.3f} +/- {:.3f})'.format(
            (f0-2.8)*1e3, u_f0*1e3), ha='center')

        return (f0-2.8)*1e3, u_f0*1e3


'''OLD'''
def analyze_dark_esr_double(timestamp = None,center_guess = False, ax=None, ret=None,min_dip_depth = 0.85 ,do_ROC = True, **kw):

    if ax == None:
        fig, ax = plt.subplots(1,1)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('DarkESR')


    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    # a.get_electron_ROC()

    x = a.sweep_pts # convert to MHz
    # y = a.get_readout_results('ssro')

    if do_ROC == True:
        a.get_electron_ROC()
        y = a.p0.reshape(-1)[:]
    else:
        y = a.get_readout_results('ssro')
    # y = a.p0.reshape(-1)[:]
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
if __name__ == '__main__':
    analyze_nmr_single()




