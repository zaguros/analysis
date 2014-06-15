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

reload(sequence)

### settings
timestamp = None#'20140603_134433' #' #'114103_PulsarD' #YYYYmmddHHMMSS

def adaptive_analysis (folder):

    #if ax == None:
    #    fig, ax = plt.subplots(1,1)
    ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration')
    print ssro_calib_folder
    a = sequence.SequenceAnalysis(folder)
  
    RO_data, phase = a.get_magnetometry_results(name='adwindata')
    print RO_data
    print phase
    phase_values = np.unique(phase)
    ssro_results = np.zeros(len(phase_values))
    ind = 0
    for j in phase_values:
        multiplicity = len(np.where(phase==j)[0])
        ssro_results [ind] = np.sum(RO_data[np.where(phase==j)])/(multiplicity+0.)
        ind = ind+1
        print 'phase value ', j, ' is present ', multiplicity, ' times'


    plt.plot(phase_values, ssro_results)
    plt.show()    
    return ssro_results

def fpga_calibration_analysis (folder):

    ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration')
    print ssro_calib_folder
    a = sequence.SequenceAnalysis(folder)
    a.get_magnetometry_phase_calibration(name='adwindata')
    a.get_electron_ROC()
    y = a.p0
    uy = a.u_p0

    a.get_sweep_pts()
    x = a.sweep_pts
    nn = a.sweep_name
    ax = a.plot_result_vs_sweepparam(ret='ax')

    do_fit = True

    if do_fit:
        guess_frq = 1./360
        guess_amp = 0.3
        guess_of = 1
        # guess_slope = 0.
        guess_phi = 0.
        guess_k = 0.


        o = fit.Parameter(guess_of, 'o')
        f = fit.Parameter(guess_frq, 'f')
        A = fit.Parameter(guess_amp, 'A')
        phi = fit.Parameter(guess_phi, 'phi')
        k = fit.Parameter(guess_k, 'k')
        p0 = [f, A, phi, o, k]
        fitfunc_str = ''


        fitfunc_str = 'o - A + A*e^(-kx)*cos(2pi (fx-phi))'

        def fitfunc(x):
            return (o()-A()) + A() * np.exp(-k()*x) * np.cos(2*np.pi*(f()*x - phi()))

        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=[0,4],
                do_print=True, ret=True)
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
                plot_data=False)

        # ax.set_title(a.timestamp+'\n'+a.measurementstring)
        plt.savefig(os.path.join(folder, 'fpga_pulse_analysis_fit.png'))




    '''

    plt.errorbar(x, y, yerr=uy, fmt = 'o')
    plt.ylim([0,1.05])
    plt.xlabel(nn)
    plt.show()    
    return x,y, uy
    '''

if __name__ == '__main__':
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
        if folder==None:
            folder = toolbox.latest_data(timestamp) 
    else:
        folder = toolbox.latest_data('adptv_estimation')
    print folder
    fpga_calibration_analysis (folder)
    #adaptive_analysis(folder)







