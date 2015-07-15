import numpy as np
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import plot
from analysis.lib.tools import toolbox

#from matplotlib import pyplot as plt #Not being used
from analysis.lib.fitting import fit, SE
reload(SE)


def electron_T2_anal(timestamp=None, measurement_name = ['ms0'],Amplitude = 1/2., T2 = 10e3, offset = 1./2, k=4, do_print = False, ylim=(0.,1.)):
    ''' Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    Based on electron_T1_anal, modified by Adriaan Rol
    '''

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('Decoupling')

    fit_results = []
    for k in range(0,len(measurement_name)):
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(measurement_name[k])
        a.get_electron_ROC()
        ax = a.plot_result_vs_sweepparam(ret='ax')
        ax.set_ylim(ylim)
        x = a.sweep_pts
        p0, fitfunc, fitfunc_str = SE.fit_echo(T2, Amplitude, offset, k)
        y = a.p0




        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True, fixed=[3])
        #fit_result = False
        if fit_result != False :
            plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)
        fit_results.append(fit_result)

        plot.ylim=([0,0.2])

    return fit_results
