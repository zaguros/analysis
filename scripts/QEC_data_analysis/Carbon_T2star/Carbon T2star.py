import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)





def Carbon_Ramsey(timestamp=None, measurement_name = ['adwindata'], ssro_folder =None,
            frequency = 1, 
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            decay_constant = 200, 
            phase =0, 
            exponent = 2, 
            plot_fit = False, do_print = False, fixed = [2], show_guess = True,
            return_phase = False,
            return_freq = False,
            return_A = False,
            return_results = True,
            close_plot = False,
            title = 'Carbon'
            ):
    ''' 
    Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data(title)

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')

    else:
        ssro_calib_folder = ssro_folder
    
    fit_results = []
    for k in range(0,len(measurement_name)):
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        ax = a.plot_results_vs_sweepparam(ret='ax')

        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]
        y_err = u_p0(-1)[:]

        ax.plot(x,y)
        p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset, amplitude, 
                x0, decay_constant,exponent,frequency ,phase )

        #plot the initial guess
        if show_guess:
            ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

        print 'fitfunction: '+fitfunc_str

        ## plot data and fit as function of total time
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

        fit_results.append(fit_result)
        if title == None:
            title = 'analyzed_result'
        plt.savefig(os.path.join(folder, title + '.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder, title + '.png'),
        format='png')
        if close_plot == True:
            plt.close()

        if return_freq == True:
            f0 = fit_result['params_dict']['f']
            u_f0 = fit_result['error_dict']['f']
            return f0, u_f0

        if return_phase == True:
            phi0 = fit_result['params_dict']['phi']
            u_phi0 = fit_result['error_dict']['phi']
            return phi0, u_phi0

        if return_phase == True and return_A == True:
            phi0 = fit_result['params_dict']['phi']
            u_phi0 = fit_result['error_dict']['phi']
            A = fit_result['params_dict']['A']
            u_A = fit_result['error_dict']['A']
            return phi0, u_phi0, A, u_A
            print 'ok'

    if return_results == True:
        return fit_results

    return x, y, y_err, fit_results


#################### Plot data ########################
folder = r'D:\measuring\data\Analyzed figures\Nuclear Ramseys'
ssro_folder = r'D:\measuring\data\20141021\081901_AdwinSSRO_SSROCalibration_111_1_sil18'
ssro_timestamp = '20141021_081901'
# C1, ms = 0

timestamp 	= '20141021_104203'
figure_name = 'Ramsey_C1_ms0'

x, y, y_err, fit_results =  Carbon_Ramsey(timestamp=timestamp,ssro_folder =ssro_folder)

fig = self.default_fig(figsize=figsize)
ax = self.default_ax(fig)

ax.errorbar(x,y,y_err = y_err, color = 'r',ls = '', marker = 'o')
plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

plt.savefig(os.path.join(folder, figure_name + '.pdf'),
format='pdf')