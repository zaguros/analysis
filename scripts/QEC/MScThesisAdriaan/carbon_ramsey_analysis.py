import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)


def Carbon_Ramsey(timestamp=None, measurement_name = ['adwindata'], ssro_calib_timestamp =None,
        frequency = 1, offset = 0.5, x0 = 0,  amplitude = 0.5,  decay_constant = 200,phase =0, exponent = 2,
        plot_fit = False, do_print = False, fixed = [2], show_guess = True, print_info = True,
        figsize = (3,2), linewidth = 2, markersize = 2, fontsize =10,
        title ='' ,savename ='Gen_exp_fit'):
    '''
    Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed')
    offset, amplitude, decay_constant,exponent,frequency ,phase
    '''
    plt.rc('font', size=fontsize)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('Carbon')

    if ssro_calib_timestamp == None:
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print ssro_calib_folder


    fit_results = []
    for k in range(0,len(measurement_name)):
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        ax = a.plot_results_vs_sweepparam(ret='ax')

        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]

        p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset, amplitude,
                x0, decay_constant,exponent,frequency ,phase )



        #plot the initial guess
        if show_guess:
            ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=fixed)

        print 'fitfunction: '+fitfunc_str

        ## plot data and fit as function of total time
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False, print_info =print_info, lw = linewidth)

        ### Beautify plot
        ax.set_xlim(0.0,x[-1])
        ax.set_ylim(-0.05,1.05)
        ax.yaxis.set_ticks( [0,0.25,0.5,0.75,1])

        locs,labels = plt.xticks()
        labels = []
        for l in locs :
            labels.append(str(l*1e3))
        plt.xticks(locs,labels)

        ax.set_title(title)
        ax.set_xlabel(r'Free evolution time (ms)',fontsize =fontsize)
        ax.set_ylabel(r'$F$ $\left( |0\rangle \right)$', fontsize = fontsize)



        fit_results.append(fit_result)

        plt.savefig(os.path.join(folder, savename+'.pdf'),
        format='pdf',bbox_inches='tight')


    return fit_results
