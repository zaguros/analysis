import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(plot)

def CosineSum_MBI_data(timestamp=None, measurement_name = ['adwindata'], ssro_calib_timestamp =None,
        frequency = [1,1], offset =0.5, amplitude =[ 0.5,0.5],  phase =[0,0],
        fixed = [],
        plot_fit = False, do_print = False, show_guess = True, print_info = True,
        figsize = (3,2), linewidth = 2, markersize = 2, fontsize =10,
        title ='' ,savename ='Cosine_sum2'):
    '''
    Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed')
    [freq, offset, Amplitude, phase]
    '''
    plt.rc('font', size=fontsize)


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




    # ax = a.plot_results_vs_sweepparam(ret='ax',figsize=figsize,markersize=0.1)
    fig = a.default_fig(figsize=figsize)
    ax = a.default_ax(fig)
    # ax.set_xlim(xlims)


    for i in range(a.readouts):
        ax.errorbar(a.sweep_pts, a.p0[:,i], fmt='o',
            yerr=a.u_p0[:,i],markersize=markersize)
    start, end = ax.get_xlim()
    print start
    print end
    # ax.xaxis.set_ticks(np.arange(start, end, .25))
    ax.xaxis.set_ticks( np.arange(start, end+1e-6, 15e-6))
    ax.set_ylim(-0.05,1.05)
    ax.yaxis.set_ticks( [0,0.5,1])




    #ax.set_markersize = 0.1
    x = a.sweep_pts.reshape(-1)[:]
    y= a.p0.reshape(-1)[:]

    p0, fitfunc, fitfunc_str = common.fit_sum_2cos(offset,amplitude[0],frequency[0],phase[0],amplitude[1],frequency[1],phase[1])
    if show_guess:
        ax.plot(np.linspace(0e-6,x[-1],201), fitfunc(np.linspace(0e-6,x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=fixed)
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(0e-6,x[-1],201), ax=ax,
                plot_data=False,print_info = print_info, lw = linewidth)


    fit.write_to_file(fit_result,folder,fitname = 'Sum of cosine fit')



    ## plot data and fit as function of total time
    locs,labels = plt.xticks()
    labels = []
    for l in locs :
        labels.append(str(l*1e6))
    plt.xticks(locs,labels)

    ax.set_title(title)
    ax.set_xlabel(r'Free evolution time ($\mu$s)',fontsize =fontsize)
    ax.set_ylabel(r'$F$ $\left( |0\rangle \right)$', fontsize = fontsize)



    plt.savefig(os.path.join(folder, savename+'.pdf'),
    format='pdf',bbox_inches='tight')
    plt.savefig(os.path.join(folder, savename+'png'),
    format='png',bbox_inches='tight')
    return fig



#CR.Carbon_Ramsey(timestamp='20140521164658', measurement_name = ['adwindata'],
#            frequency = [338e-3,14e-3], amplitude = [1,1],  decay_constant = [0.009e6,0.009e6],phase =[0,0],
#            fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False)

#CR.Carbon_Ramsey(timestamp='20140521165249', measurement_name = ['adwindata'],
#            frequency = [305e-3,16e-3], amplitude = [1,1],  decay_constant = [0.009e6,0.009e6],phase =[0,0],
#            fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False)

#CR.Carbon_Ramsey(timestamp='20140521170735', measurement_name = ['adwindata'],
#            frequency = [350e-3,27e-3], amplitude = [1,1],  decay_constant = [0.009e6,0.009e6],phase =[0,0],
#            fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False)

#CR.Carbon_Ramsey(timestamp='20140521162939', measurement_name = ['adwindata'],
#            frequency = [344e-3,19e-3], amplitude = [1,1],  decay_constant = [0.009e6,0.009e6],phase =[0,0],
#            fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = False)
