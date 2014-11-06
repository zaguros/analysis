"""
Script to analyze the C13 T1 data
"""

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)


def Carbon_T1_analysis(measurement_name = ['adwindata'], ssro_calib_timestamp =None, 
            offset = 0.5, 
            amplitude = 0.5,  
            decay_constant = 0.2, 
            exponent = 1, 
            Addressed_carbon=1,
            plot_fit = True, do_print = True, show_guess = False):
    ''' 
    Function to gather and analyze T1 measurements of a specific carbon.

    Addressed_carbon: selects the used timestamps and therefore the analyzed carbon

    measurement_name: list of measurement names

    Possible inputs for initial guesses: offset, amplitude, decay_constant,exponent
    '''
    
    ######################################
    # el_RO = positive & carbon_init= up #
    ######################################


    if Addressed_carbon == 1:
        timestamp=['20141029_183047','20141029_215155','20141029_185921','20141029_205756','20141029_175201','20141029_180826']
    elif Addressed_carbon == 5:
        timestamp=['20141030_073729','20141030_012052']

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print ssro_calib_folder


    ##accumulate data##

    cum_pts = 0

    for kk in range(len(timestamp)):
        folder = toolbox.data_from_time(timestamp[kk])

        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        cum_pts += a.pts

        if kk == 0:
            cum_sweep_pts = a.sweep_pts
            cum_p0 = a.p0
            cum_u_p0 = a.u_p0
        else:
            cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
            cum_p0 = np.concatenate((cum_p0, a.p0))
            cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

    a.pts   = cum_pts
    a.sweep_pts = cum_sweep_pts
    a.p0    = cum_p0
    a.u_p0  = cum_u_p0


    ## generate plot of the raw data ##

    fig = a.default_fig(figsize=(6,4))
    ax = a.default_ax(fig)
    ax.plot(np.log10(a.sweep_pts), a.p0, '.b', lw=1)

    ## fit to a general exponential##

    fit_results = []

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]

    ax.plot(x,y)
    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, x0, decay_constant,exponent)

    #plot the initial guess

    if show_guess:
        ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[2])


    ## plot data and fit as function of total time

    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    fit_results.append(fit_result)

    filename= 'C13_T1_analysis_up_positive_C'+str(Addressed_carbon)

    print 'plots are saved in ' + folder
    plt.savefig(os.path.join(folder, filename+'.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, filename+'analyzed_result.png'),
    format='png')


if __main__ == '__main__':

    Carbon_T1_analysis(Addressed_carbon=1, ssro_calib_timestamp =None)
    Carbon_T1_analysis(Addressed_carbon=5, ssro_calib_timestamp =None)

#TODO:  fill in timestamps for ssro calibration
#       gather all timestamps for T1 measurements. (done so far: up_positive)