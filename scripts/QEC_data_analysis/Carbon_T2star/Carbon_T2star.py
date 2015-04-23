import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot; reload(plot)
from analysis.lib.fitting import fit, common; reload(common)
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt

def Carbon_Ramsey(timestamp=None, measurement_name = ['adwindata'], ssro_folder =None,
            frequency = 1, 
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            decay_constant = 0.015, 
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
    '''

    folder = toolbox.data_from_time(timestamp)
    ssro_calib_folder = ssro_folder

    fit_results = []

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder)
    ax = a.plot_results_vs_sweepparam(ret='ax')

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    y_err = a.u_p0.reshape(-1)[:]

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

    return x, y, y_err, fit_result


#################### Plot data ########################
def plot_Ramsyes(folder = r'D:\measuring\data\Analyzed figures\Nuclear Ramseys'):

    ### Data locations
    timestamp_list  = ['20141021_104203','20141021_174316','20141024_162620','20141024_163936','20141021_180647','20141021_182509']
    
    SSRO_timestamp, ssro_folder_C1_C5 = toolbox.latest_data(contains = 'AdwinSSRO', older_than = '20141021_081902',return_timestamp = True)
    SSRO_timestamp, ssro_folder_C2    = toolbox.latest_data(contains = 'AdwinSSRO', older_than = '20141024_082135',return_timestamp = True)
    SSRO_folder_list = [ssro_folder_C1_C5,ssro_folder_C1_C5,ssro_folder_C2,ssro_folder_C2,ssro_folder_C1_C5,ssro_folder_C1_C5]
    
    # C1, ms = 0
    figure_name_list = ['Ramsey_C1_ms0','Ramsey_C1_ms1','Ramsey_C2_ms0','Ramsey_C2_ms1','Ramsey_C5_ms0','Ramsey_C5_ms1']
    color_list  = ['r','r','b','b','g','g']
    x_max_list  = [14e-3,16e-3,16e-3,16e-3,32e-3,32e-3]
    f_list      = [250,200,288,150,100,110]
    decay_list  = [9e-3,9e-3,12e-3,12e-3,21e-3,21e-3]
    spacing     = [4e-3]*4+[5e-3]*2

    for ii, timestamp in enumerate(timestamp_list):
        figure_name = figure_name_list[ii]
        ssro_folder = SSRO_folder_list[ii]
        color = color_list[ii]

        x, y, y_err, fit_result =  Carbon_Ramsey(timestamp=timestamp,ssro_folder =ssro_folder,frequency = f_list[ii],decay_constant = decay_list[ii])

        xticks = np.arange(0,x_max_list[ii]+1e-3,spacing[ii])
        yticks = np.linspace(0,1,3)
        xtick_rescaled = (xticks*1e3).astype(int)
        fig,ax = plt.subplots()


        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False,add_txt = False,linestyle = '-',linewidth = 3,color = '0.25')

        errlines = ax.errorbar(x,y,yerr = y_err, color = color,ls = '', marker = 'o',markersize = 8,capsize=5, lw = 2)
    	# set(errlines, linewidth=6)

        ax.set_ylim([0,1])
        ax.set_xlim([0-0.1e-3,x_max_list[ii]+0.1e-3])
        plt.xticks(xticks)
        plt.yticks(yticks)
        ax.set_xticklabels(xtick_rescaled)
        ax.set_xlabel('Free evolution time (ms)',fontsize = 25)
        ax.set_ylabel(r'F$(|$x$\rangle)$',fontsize = 30)
        # ax.hlines([0.5],x[0]-1,x[-1]+1,linestyles='dotted',linewidth = 1.5)
        ax.set_ylim([0,1])

        ax.tick_params(axis='x', which='major', labelsize=25)
        ax.tick_params(axis='y', which='major', labelsize=25)
        plt.rcParams['axes.linewidth'] = 2
        ax.tick_params('both', length=4, width=2, which='major')
        # plt.savefig(os.path.join(folder, figure_name+'_'+timestamp_list[ii] + '.pdf'),format='pdf',bbox_inches='tight')
        # plt.savefig(os.path.join(folder, figure_name+'_'+timestamp_list[ii] + '.png'),format='png',bbox_inches='tight')
