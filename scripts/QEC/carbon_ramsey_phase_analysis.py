import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(plot)


def Carbon_Ramsey(timestamp=None, measurement_name = ['adwindata'], 
            frequency = [1], amplitude = [0.5],  decay_constant = [200],phase =[0], 
            fitfunc_type = 'single', plot_fit = False, do_print = False, show_guess = True,fitexp = None):
    ''' Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    '''
    offset = 0.5
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
        a.sweep_pts = a.sweep_pts*1e6
        ax = a.plot_results_vs_sweepparam(ret='ax')
        ax.set_ylim(-0.05,1.05)

        ax.set_xlim(-0.05,40)

        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]
        print fitexp
        print
        print 
        if len(frequency) == 1:
            if fitexp == 'Gaussian':
                p0, fitfunc, fitfunc_str = common.fit_gaussian_decaying_cos(frequency[0], offset, amplitude[0], phase[0], decay_constant[0])
                print 'ok'
            else:
                p0, fitfunc, fitfunc_str = common.fit_decaying_cos(frequency[0], offset, amplitude[0], phase[0], decay_constant[0])
            #plot the initial guess
            print 'no'
            if show_guess:
                ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[4])
        elif len(frequency) == 2:
            print 'yes'
            p0, fitfunc, fitfunc_str = common.fit_double_decaying_cos(frequency[0], amplitude[0], phase[0], decay_constant[0], frequency[1], amplitude[1], phase[1], decay_constant[1])
            if show_guess:
                ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[2,3,6,7])
           
            ### Also plot the individual curves
            
            # p0_0, fitfunc_0, fitfunc_str_0 = common.fit_double_decaying_cos(fit_result['params'][0], 0, phase[0], fit_result['params'][2], fit_result['params'][3], fit_result['params'][4], phase[1], fit_result['params'][5])
            # ax.plot(np.linspace(0,x[-1],201), fitfunc_0(np.linspace(0,x[-1],201)), 'b-', lw=1)
            # p0_1, fitfunc_1, fitfunc_str_1 = common.fit_double_decaying_cos(fit_result['params'][0], fit_result['params'][1], phase[0], fit_result['params'][2], fit_result['params'][3],0, phase[1], fit_result['params'][5])
            # ax.plot(np.linspace(0,x[-1],201), fitfunc_1(np.linspace(0,x[-1],201)), 'm-', lw=1)



        ## plot fit
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False,print_info = False)

        fit_results.append(fit_result)
        ax.set_xlabel('Free evolution time (us)')

        print folder
        plt.savefig(os.path.join(folder, 'analyzed_result2.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder, 'analyzed_result2.png'),
        format='png')

        # freq = fit_results[0]['params_dict']['f1']
        # period = 1/freq 
        # print 'Period is %s pulses ' %(period)
        # N_pi = round(period*.5/2)*2.0 
        # N_pi2 = round(period/2*.25)*2.0
        # print 'Pi pulse: %s pulses' %N_pi
        # print 'Pi2 pulse: %s pulses' %N_pi2
    return fit_results

######################### RAMSEY T2* ##############################################

# Carbon_Ramsey(timestamp='20140519135216', measurement_name = ['adwindata'], 
#             frequency = [575], amplitude = [0.5],  decay_constant = [0.009],phase =[0], 
#             fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False,fitexp='Gaussian')

# Carbon_Ramsey(timestamp='20140519183105', measurement_name = ['adwindata'], 
#             frequency = [575], amplitude = [0.5],  decay_constant = [0.007],phase =[0], 
#             fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False,fitexp='Gaussian')

# Carbon_Ramsey(timestamp='20140520134826', measurement_name = ['adwindata'], 
#             frequency = [205], amplitude = [0.5],  decay_constant = [0.007],phase =[0], 
#             fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False,fitexp='Gaussian')

######################## RAMSEY 2 freq ##############################################

Carbon_Ramsey(timestamp='20140521164658', measurement_name = ['adwindata'], 
            frequency = [338e-3,14e-3], amplitude = [1,1],  decay_constant = [0.009e6,0.009e6],phase =[0,0], 
            fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False)

Carbon_Ramsey(timestamp='20140521165249', measurement_name = ['adwindata'], 
            frequency = [305e-3,16e-3], amplitude = [1,1],  decay_constant = [0.009e6,0.009e6],phase =[0,0], 
            fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False)

Carbon_Ramsey(timestamp='20140521170735', measurement_name = ['adwindata'], 
            frequency = [350e-3,27e-3], amplitude = [1,1],  decay_constant = [0.009e6,0.009e6],phase =[0,0], 
            fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False)

Carbon_Ramsey(timestamp='20140521162939', measurement_name = ['adwindata'], 
            frequency = [344e-3,19e-3], amplitude = [1,1],  decay_constant = [0.009e6,0.009e6],phase =[0,0], 
            fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False)

# ######################### RAMSEY simple ##############################################

# Carbon_Ramsey(timestamp='20140507114220', measurement_name = ['adwindata'], 
#             frequency = [350e-3], amplitude = [-0.5],  decay_constant = [900e6],phase =[0], 
#             fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False,fitexp='Gaussian')

# Carbon_Ramsey(timestamp='20140507124543', measurement_name = ['adwindata'], 
#             frequency = [350e-3], amplitude = [-0.5],  decay_constant = [900e6],phase =[0], 
#             fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False,fitexp='Gaussian')

# Carbon_Ramsey(timestamp='20140507151219', measurement_name = ['adwindata'], 
#             frequency = [350e-3], amplitude = [-0.5],  decay_constant = [900e6],phase =[0], 
#             fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False,fitexp='Gaussian')

# Carbon_Ramsey(timestamp='20140429141113', measurement_name = ['adwindata'], 
#             frequency = [350e-3], amplitude = [-0.5],  decay_constant = [900e6],phase =[0], 
#             fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False,fitexp='Gaussian')