import numpy as np
import os, sys
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import plot
from analysis.lib.tools import toolbox

from matplotlib import pyplot as plt
from analysis.lib.fitting import fit, common
from matplotlib import pyplot as plt


def electron_T1_anal(timestamp=None, measurement_name = ['ms0'],Amplitude = 0.1, T1 = 1000, offset =1,do_print = False):
    ''' Function to analyze T1 measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    '''

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('T1')

    fit_results = []
    for k in range(0,len(measurement_name)):
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        #a.get_readout_results(measurement_name[k]) #commented out since newer measurements are no longer grouped which renders the whole for-loop obsolete. NK 20141110
        a.get_electron_ROC()
        ax = a.plot_result_vs_sweepparam(ret='ax')

        x = a.sweep_pts
        p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(offset, Amplitude, T1)
        y = a.p0
        print y
        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True)
        plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)
        fit_results.append(fit_result)

        plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder, 'analyzed_result.png'),
        format='png')


    return fit_results

def get_T1_data(folder):
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()
    x = a.sweep_pts
    y = a.p0
    y_err = a.u_p0
    return x,y,y_err





def electron_T1_anal_mul(older_than='20150227_184920',newer_than='20150227_101541',mode='contrast', Amplitude=1, Offset=1, T1=1e7):
    ''' Function to analyze T1 measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    older_than = older than timestamp
    newer_than = newer than timestamp
    mode = init_0 init_1 contrast
    '''

    Folder_list = toolbox.latest_data(contains='T1',older_than='20150227_184920',newer_than='20150227_101541',return_all=True)
    Directory = Folder_list[0]
    Date = Folder_list[1]
    init0 = [k for k in Folder_list[2] if 'init_0' in k]
    init1 = [k for k in Folder_list[2] if 'init_1' in k]


    def get_T1_data(folder): #Assures that the data is in the right order
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        a.get_electron_ROC()
        x = a.sweep_pts
        y = a.p0
        y_var = (a.u_p0)**2
        minloc = -np.where(x == min(x))[0][0]
        
        x = np.roll(x,minloc)
        y = np.roll(y,minloc)
        y_var = np.roll(y_var,minloc)

        return x,y,y_var

    if mode == 'init0' or mode == 'contrast':
        x0_tot = np.zeros((6,))
        y0_tot = np.zeros((6,))
        y0_var_tot = np.zeros((6,))

        for i in range(len(init0)):
            Folder = Directory + "\\" + Date + "\\" + init0[i]
            x,y,y_var = get_T1_data(Folder)
            y0_tot += y
            y0_var_tot += y_var

        y0_tot /= len(init0)
        y0_var_tot /= len(init0)


    if mode == 'init1' or mode == 'contrast':
        x1_tot = np.zeros((6,))
        y1_tot = np.zeros((6,))
        y1_var_tot = np.zeros((6),)

        for i in range(len(init1)):
            Folder = Directory + "\\" + Date + "\\" + init1[i]
            x,y,y_var = get_T1_data(Folder)
            y1_tot += y
            y1_var_tot += y_var
        
        y1_tot /= len(init0)
        y1_var_tot /= (len(init0)) 

    if mode == 'init1':
        y_tot = y1_tot
        y_var_tot = y1_var_tot
    elif mode == 'init0':
        y_tot = y1_tot
        y_var_tot = y1_var_tot
    elif mode == 'contrast':
        y_diff = y0_tot - y1_tot
        var_diff = (y0_var_tot + y1_var_tot) / 2
    else:
        raise Exception('Mode not specified')


    p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(Amplitude, Offset, T1)


    a = sequence.SequenceAnalysis(Directory + "\\" + Date + "\\" + init0[0])
    a.get_sweep_pts()
    a.sweep_pts = x / 1e6
    a.sweep_name = 'Times (sec)'
    a.get_readout_results('ssro')
    a.get_electron_ROC()

    a.p0 = y_diff
    a.u_p0 = var_diff**0.5
    ax = a.plot_result_vs_sweepparam(ret='ax')

    fit_result = fit.fit1d(x / 1e6,y_diff, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(0,x[-1]/1e6,201), ax=ax, plot_data=False)