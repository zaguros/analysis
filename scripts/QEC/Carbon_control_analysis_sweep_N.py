import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(plot)

def Carbon_control_sweep_N(timestamp='20140507125650', measurement_name = ['adwindata'], 
            frequency = [1], amplitude = [0.5],  decay_constant = [200],phase =[0], offset = 0.5,
            fitfunc_type = 'single', plot_fit = False, do_print = False, show_guess = True,
            yaxis = [-0.5,1.05], fixed=[2,6]):
    ''' Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    '''
    # offset = 0.5
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('Decoupling')

    fit_results = []
    for k in range(0,len(measurement_name)):
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC()
        ax = a.plot_results_vs_sweepparam(ret='ax')
        ax.set_ylim(yaxis)

        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]

        if len(frequency) == 1:
            
            p0, fitfunc, fitfunc_str = common.fit_decaying_cos(frequency[0], offset, amplitude[0], phase[0], decay_constant[0])
            #plot the initial guess
            if show_guess:
                ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[3])
        elif len(frequency) == 2:
            p0, fitfunc, fitfunc_str = common.fit_double_decaying_cos(frequency[0], amplitude[0], phase[0], 
                            decay_constant[0], frequency[1], amplitude[1], phase[1], decay_constant[1],offset)
            if show_guess:
                ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
           
            ## Also plot the individual curves
            
            p0_0, fitfunc_0, fitfunc_str_0 = common.fit_double_decaying_cos(fit_result['params'][0], 0, phase[0], fit_result['params'][2], fit_result['params'][3], fit_result['params'][4], phase[1], fit_result['params'][5], fit_result['params'][6])
            ax.plot(np.linspace(0,x[-1],201), fitfunc_0(np.linspace(0,x[-1],201)), 'r-', lw=1,alpha=0.5)
            p0_1, fitfunc_1, fitfunc_str_1 = common.fit_double_decaying_cos(fit_result['params'][0], fit_result['params'][1], phase[0], fit_result['params'][2], fit_result['params'][3],0, phase[1], fit_result['params'][5], fit_result['params'][6])
            ax.plot(np.linspace(0,x[-1],201), fitfunc_1(np.linspace(0,x[-1],201)), 'm-', lw=1,alpha=0.5)

        ## plot fit
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False,print_info = True)

        fit_results.append(fit_result)
        print folder
        plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder, 'analyzed_result.png'),
        format='png')

        # freq = fit_results[0]['params_dict']['f1']
        # period = 1/freq 
        # print 'Period is %s pulses ' %(period)
        # # N_pi = round(period*.5/2)*2.0 
        # N_pi2 = round(period/2*.25)*2.0
        # # print 'Pi pulse: %s pulses' %N_pi
        # print 'Pi2 pulse: %s pulses' %N_pi2
    return fit_results


def Carbon_control_sweep_N_zoom(timestamp=None, measurement_name = ['adwindata'], 
            A = [0.5, 0.5],
            fitfunc_type = 'single', plot_fit = False, do_print = False, show_guess = True,
            yaxis = [-0.05,1.05]):
    ''' Function to analyze data for optimization of the number of pulses for a controlled C13 gate. 
    '''
    
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('Decoupling')

    fit_results = []
    for k in range(0,len(measurement_name)):
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC()
        ax = a.plot_results_vs_sweepparam(ret='ax')
        ax.set_ylim(yaxis)
        ax.axhspan(0,0.5,fill=False,ls='dotted')

        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]
        y_u = a.u_p0.reshape(-1)[:]

        ax.set_xlim(x[0]-1,x[-1]+1)

        p0, fitfunc, fitfunc_str = common.fit_poly(A)
        
        if show_guess:
            ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)
        
        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[3])
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax, plot_data=False,print_info = True)

        fit_results.append(fit_result)
        print folder
        plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder, 'analyzed_result.png'),
        format='png')

        diff = np.abs(y - 0.5)
        print diff
        print 'Optimum number of pulses N = ' + str(x[np.argmin(diff)])
        print 'with y-0.5 = ' + str(y[np.argmin(diff)]-0.5) + ' +/- ' + str(y_u[np.argmin(diff)])

        # freq = fit_results[0]['params_dict']['f1']
        # period = 1/freq 
        # print 'Period is %s pulses ' %(period)
        # # N_pi = round(period*.5/2)*2.0 
        # N_pi2 = round(period/2*.25)*2.0
        # # print 'Pi pulse: %s pulses' %N_pi
        # print 'Pi2 pulse: %s pulses' %N_pi2
    return fit_results

