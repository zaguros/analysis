import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(fit)
reload(common)
reload(mbi)

def exp_sin(contains = '',timestamp=None, measurement_name = ['adwindata'],ssro_folder = None,
            offset=[0], amplitude = [0.5], center = [0], decay_constant = [200], exp_power = [0],
            frequency = [1], phase =[0],
            fixed = [], ylim = [-0.5, 1.05],ssro_tstamp ='',base_folder = None,
            plot_fit = False, do_print = False, show_guess = True,correct_ionization = True):
    ''' Function to fit mbi-type data with exponential and sinusoidal functions or combinations thereof.
    timestamp       : format yyyymmdd_hhmmss or hhmmss or None. None takes the last data.
    measurement_name: list of measurement names
    '''
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp,folder = base_folder)
    else:
        folder = toolbox.latest_data(contains,folder = base_folder)
    print folder
    if ssro_folder == None:
        if ssro_tstamp == '':
            ssro_folder = ssro_tstamp

        else:
            ssro_folder = toolbox.data_from_time(ssro_tstamp)
        
    fit_results = []

    for k in range(0,len(measurement_name)):
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name=measurement_name[k],CR_after_check = correct_ionization)
        a.get_electron_ROC(ssro_folder)
        ax = a.plot_results_vs_sweepparam(ret='ax')

        if ylim != None:
            ax.set_ylim(ylim[0],ylim[1])

        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]

        ### fit depending on the number of the frequencies
        if len(frequency) == 1:
            print exp_power[0]
            p0, fitfunc, fitfunc_str = common.fit_exp_cos(offset[0],
                    amplitude[0], center[0], decay_constant[0], exp_power[0],
                    frequency[0], phase[0])
            if show_guess:
                ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(0,x[-1],201)), lw=2)
            print 'starting fit.fit1d'
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
            fit_result['y_u'] = a.u_p0.reshape(-1)[:]
        elif len(frequency) == 2:
            p0, fitfunc, fitfunc_str = common.fit_gaussian_decaying_2cos(offset[0],amplitude[0],decay_constant[0],amplitude[0],
                frequency[0],  phase[0], amplitude[1], frequency[1],  phase[1])
            if show_guess:
                ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
            fit_result['y_u'] = a.u_p0.reshape(-1)[:]
            ### Also plot the individual curves

            # p0_0, fitfunc_0, fitfunc_str_0 = common.fit_double_decaying_cos(fit_result['params'][0], 0, phase[0], fit_result['params'][2], fit_result['params'][3], fit_result['params'][4], phase[1], fit_result['params'][5])
            # ax.plot(np.linspace(0,x[-1],201), fitfunc_0(np.linspace(0,x[-1],201)), 'b-', lw=1)
            # p0_1, fitfunc_1, fitfunc_str_1 = common.fit_double_decaying_cos(fit_result['params'][0], fit_result['params'][1], phase[0], fit_result['params'][2], fit_result['params'][3],0, phase[1], fit_result['params'][5])
            # ax.plot(np.linspace(0,x[-1],201), fitfunc_1(np.linspace(0,x[-1],201)), 'm-', lw=1)

        ## plot fit
        if plot_fit == True:

            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax, plot_data=False)

        fit_results.append(fit_result)
        print folder
        plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder, 'analyzed_result.png'),
        format='png')

        plt.show()
    return fit_results


def general_exponential(contains = '',timestamp=None, measurement_name = ['adwindata'],ssro_folder = None,
            offset=[0], amplitude = [0.5], center = [0], decay_constant = [200], exp_power = [0],
            fixed = [], ylim = [-0.5, 1.05],ssro_tstamp ='',base_folder = None,
            plot_fit = False, do_print = False, show_guess = True):
    
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp,folder = base_folder)
    else:
        folder = toolbox.latest_data(contains,folder = base_folder)

    if ssro_folder == None:
        if ssro_tstamp == '':
            ssro_folder = ssro_tstamp

        else:
            ssro_folder = toolbox.data_from_time(ssro_tstamp)
        
    fit_results = []

    for k in range(0,len(measurement_name)):
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name=measurement_name[k])
        a.get_electron_ROC(ssro_folder)
        ax = a.plot_results_vs_sweepparam(ret='ax')

        if ylim != None:
            ax.set_ylim(ylim[0],ylim[1])

        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]


        p0, fitfunc, fitfunc_str = common.fit_general_exponential_fixed_offset(offset,
                amplitude, center, decay_constant, exp_power)

        if show_guess:
            ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(0,x[-1],201)), lw=2)

        print 'starting fit.fit1d'
        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
        fit_result['y_u'] = a.u_p0.reshape(-1)[:]

        ## plot fit
        if plot_fit == True:

            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax, plot_data=False)

        fit_results.append(fit_result)
        print folder
        plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder, 'analyzed_result.png'),
        format='png')

        plt.show()
    return fit_result

def get_CR_histos(contains='',timestamp = None,):
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data(contains)

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()

    before,after = a.get_CR_before_after()

    fig = plt.figure()
    ax = plt.subplot()
    ax.hist(before,abs(max(before)-min(before)+1),normed=True,label = 'before')
    ax.hist(after,abs(max(after)-min(after)+1),normed=True,label = 'after')
    ax.set_title(a.default_plot_title)
    ax.set_xlabel('counts during CR check')
    ax.set_ylabel('probability')
    plt.legend()

    plt.show()
    plt.close('all')
