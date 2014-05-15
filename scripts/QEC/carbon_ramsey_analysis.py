import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)

reload(common)

def Carbon_Ramsey(timestamp=None, measurement_name = ['adwindata'],frequency = 1, offset = 0.5, amplitude = 0.5,  decay_constant = 200,phase =0, plot_fit = False, do_print = False, show_guess = True):
    ''' Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    '''

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('CarbonRamsey')

    fit_results = []
    for k in range(0,len(measurement_name)):
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC()
        ax = a.plot_results_vs_sweepparam(ret='ax')

        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]




        ax.plot(x,y)
        p0, fitfunc, fitfunc_str = common.fit_decaying_cos(frequency,offset, amplitude, phase, decay_constant)

        #plot the initial guess
        if show_guess:
            ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)

        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[3])


        ## plot data and fit as function of total time
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)

        fit_results.append(fit_result)

        plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder, 'analyzed_result.png'),
        format='png')

    return fit_results