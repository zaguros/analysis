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
from analysis.lib.m2.ssro import sequence
import analysis.scripts.pulse_calibration.calibration_funcs as funcs; reload(funcs)


def analyse_pulse_calibration(name='Pi', timestamp=None, guess_x0 = None):
    ### parameters
    if guess_x0 == None:
        guess_x0 = 0.5

    msmt_type = 'sequence'

    guess_of = 0.973
    guess_a = 0.

    ### script
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data(name)

    if msmt_type == 'sequence':
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        a.get_electron_ROC()
        ax = a.plot_result_vs_sweepparam(ret='ax')

        x = a.sweep_pts
        y = a.p0

    elif msmt_type == 'mbi':
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC()
        ax = a.plot_results_vs_sweepparam(name='adwindata', ret='ax')
    
    else:
        raise Exception('Unknown msmt type')

    ax.set_ylim(-0.1,1)

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]


    res = funcs.calibrate_pulse_amplitude(x, y, ax, guess_x0, guess_of, guess_a)
    plt.savefig(os.path.join(folder, 'pulse_calibration.pdf'),
            format='pdf')
    plt.savefig(os.path.join(folder, 'pulse_calibration.png'),
            format='png')



def analyse_pi2_calibration(name='Pi2', timestamp=None, guess_x0 = None):
   

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data(name)


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
    ax1.set_title(a.timestamp+'\nDifference btw. Pi/2-Pi and Pi/2')
    ax1.set_ylabel('Difference')
    ax2.set_title(a.timestamp)
    m = fit.Parameter((y[-1]-y[0])/(x[-1]-x[0]), 'm')
    x0 = fit.Parameter(x2.mean(), 'x0')
    p0 = [m, x0]

    def ff(x):
        return m() * (x-x0())
    fitfunc_str = 'm * (x - x0)'

    fit_result = fit.fit1d(x2, y2, None, p0=p0, fitfunc=ff,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)    
        
    ax2.errorbar(x2, y[0::2], yerr=u_y[0::2], fmt='o',
                 label='Pi/2 - Pi')
    ax2.errorbar(x2, y[1::2], yerr=u_y[1::2], fmt='o',
                 label='Pi/2')
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


def bloch_vector_length(name = '',multiplicity = '01',timestamp=None, guess_x0 = None):
    ### parameters
    if guess_x0 == None:
        guess_x0=0.5

    guess_of = 0.973
    guess_a = 0.

    ### script
    x_dict = {}
    y_dict = {}
    y_err_dict = {}

    for RO in ['X','-X','Y','-Y','Z','-Z']:
        if timestamp != None:
            folder = toolbox.data_from_time(timestamp)
        else:
            folder = toolbox.latest_data(name+RO+multiplicity)
            print folder
        
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        a.get_electron_ROC()
        # ax = a.plot_result_vs_sweepparam(ret='ax')
        x_dict[RO] ={}
        y_dict[RO] ={}
        y_err_dict[RO] ={}

        x_dict[RO] = a.sweep_pts
        y_dict[RO] = ((a.p0.reshape(-1)[:])-0.5)*2
        y_err_dict[RO] = 2*a.u_p0.reshape(-1)[:] 

    y_f = (((y_dict['X']-y_dict['-X'])/2.)**2+((y_dict['Y']-y_dict['-Y'])/2.)**2+((y_dict['Z']-y_dict['-Z'])/2.)**2)**(1/2.)


    print y_f
    

    fig,ax = plt.subplots()
    ax.plot(x_dict['X'],y_f,'o',color = 'k')
    ax.set_ylim(-0.1,1.1)
    # ax.set_xlim(x_dict['X'][0]-0.01,x_dict['X'][:]+0.01)
    ax.set_title(folder)
    # ax.hlines([-1,0],x_dict['X'][0]-1,x_dict['X'][-1]+1,linestyles='dotted')
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Bloch vector length')

    res = funcs.calibrate_pulse_amplitude(x_dict['X'], y_f, ax, guess_x0, guess_of, guess_a)

    plt.savefig(os.path.join(folder, 'bloch_length.pdf'),
            format='pdf')
    plt.savefig(os.path.join(folder, 'bloch_length.png'),
            format='png')



if __name__ == '__main__':
    analyse_pulse_calibration(angle='_pi_1')



