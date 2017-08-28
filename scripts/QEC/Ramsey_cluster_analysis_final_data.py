import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(mbi)


def Carbon_Ramsey(timestamp=None, carbon=None, transition=None, measurement_name = ['adwindata'], ssro_calib_timestamp =None,
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
            return_amp = False,
            return_results = True,
            close_plot = False,
            color='b',
            title = 'Carbon'):
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
    elif carbon != None:
        folder = toolbox.latest_data(contains='C'+str(carbon)+'_ms'+str(transition))
    else:
        folder = toolbox.latest_data(title)

    print folder

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSROCalibration')

    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
        print ssro_calib_folder


    fit_results = []
    for k in range(0,len(measurement_name)):
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        
        x = 1.0e3*a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]
        e= a.u_p0
        
        fig = plt.figure(1,figsize=(8,6))
        ax2 = fig.add_subplot(111)
        ax2.errorbar(x.flatten(),y.flatten(),yerr=e,fmt='o',label='',color=color)
        ax2.set_xlabel('Free evolution time (ms)')
        ax2.set_ylabel('State Fidelity')
        p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset, amplitude, 
                x0, decay_constant,exponent,frequency ,phase )

        #plot the initial guess
        if show_guess:
            ax2.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

        print 'fitfunction: '+fitfunc_str

        ## plot data and fit as function of total time
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax2, add_txt = False, color=color, plot_data=False)

        fit_results.append(fit_result)
        if title == None:
            title = 'analyzed_result'
        

        plt.savefig(os.path.join(folder, title + '.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder, title + '.png'),
        format='png')
        if close_plot == True:
            plt.close()

        # for item in fit_result['params_dict']:
        #     print item
        # return fit_result
        print folder

        if return_freq == True:
            f0 = fit_result['params_dict']['f']
            u_f0 = fit_result['error_dict']['f']
            return f0, u_f0

        if return_phase == True and return_amp == False:
            phi0 = fit_result['params_dict']['phi']
            u_phi0 = fit_result['error_dict']['phi']
            return phi0, u_phi0

        if return_phase == True and return_amp == True:
            phi0 = fit_result['params_dict']['phi']
            u_phi0 = fit_result['error_dict']['phi']
            A = fit_result['params_dict']['A']
            u_A = fit_result['error_dict']['A']
            return phi0, u_phi0, A, u_A

        if return_amp == True:
            A = fit_result['params_dict']['A']
            u_A = fit_result['error_dict']['A']
            return A, u_A

    if return_results == True:
        return fit_results