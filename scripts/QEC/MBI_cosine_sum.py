import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(fit) 
reload(toolbox)


def CosineSum_MBI_data(timestamp=None, measurement_name = ['adwindata'], ssro_calib_timestamp =None,
        frequency = [1,1], offset =0.5, amplitude =[ 0.5,0.5],  phase =[0,0], 
        fixed = [], 
        plot_fit = False, do_print = False, show_guess = True):
    ''' 
    Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    [freq, offset, Amplitude, phase] 
    '''


    if timestamp == None:
        timestamp, folder   = toolbox.latest_data('CarbonR',return_timestamp =True)
    else: 
        folder = toolbox.data_from_time(timestamp) 

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'

    fit_result = [None]

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder)
    ax = a.plot_results_vs_sweepparam(ret='ax')
    x = a.sweep_pts.reshape(-1)[:]
    y= a.p0.reshape(-1)[:]

    p0, fitfunc, fitfunc_str = common.fit_sum_2cos(offset,amplitude[0],frequency[0],phase[0],amplitude[1],frequency[1],phase[1]) 
    if show_guess:
        ax.plot(np.linspace(0e-6,x[-1],201), fitfunc(np.linspace(0e-6,x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=fixed)
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(0e-6,x[-1],201), ax=ax, 
                plot_data=False,print_info = True)
    
    fit.write_to_file(fit_result,folder,fitname = 'Sum of cosine fit') 


    ## plot data and fit as function of total time

    plt.savefig(os.path.join(folder, 'CosineSumFit.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'CosineSumFit.png'),
    format='png')
    return fitfunc

    # return fit_result
