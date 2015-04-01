import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import pickle
reload(common)
reload(plot)

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

    if ssro_folder == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')

    else:
        ssro_calib_folder = ssro_folder
    
    fit_results = []
    for k in range(0,len(measurement_name)):
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        # ax = a.plot_results_vs_sweepparam(ret='ax')

        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]
        y_err = a.u_p0.reshape(-1)[:]

        # ax.plot(x,y)
        p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset, amplitude, 
                x0, decay_constant,exponent,frequency ,phase )

        #plot the initial guess
        # if show_guess:
        #     ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

        print 'fitfunction: '+fitfunc_str

        ## plot data and fit as function of total time
        # if plot_fit == True:
        #     plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

        fit_results.append(fit_result)
        if title == None:
            title = 'analyzed_result'
        # plt.savefig(os.path.join(folder, title + '.pdf'),
        # format='pdf')
        # plt.savefig(os.path.join(folder, title + '.png'),
        # format='png')
        # if close_plot == True:
        #     plt.close()

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

    # if return_results == True:
    #     return fit_results

    return x, y, y_err, fit_result


#################### Plot data ########################

folder = r'D:\measuring\data\Analyzed figures\Nuclear Ramseys'
folder = r'D:\measuring\data\QEC_data\figs\final figures'

ssro_folder_C1_C5 = r'D:\measuring\data\20141021\081901_AdwinSSRO_SSROCalibration_111_1_sil18'
ssro_folder_C2 = r'D:\measuring\data\20141024\082134_AdwinSSRO_SSROCalibration_111_1_sil18'
ssro_timestamp = '20141021_081901'
# C1, ms = 0

timestamp_list 	= ['20141021_104203','20141021_174316','20141024_162620','20141024_163936','20141021_180647','20141021_182509']
SSRO_folder_list = [ssro_folder_C1_C5,ssro_folder_C1_C5,ssro_folder_C2,ssro_folder_C2,ssro_folder_C1_C5,ssro_folder_C1_C5]
figure_name_list = ['Ramsey_C1_ms0','Ramsey_C1_ms1','Ramsey_C2_ms0','Ramsey_C2_ms1','Ramsey_C5_ms0','Ramsey_C5_ms1']

c_green = (9/255.,232/255.,94/255.)
c_grey = (64/255.,78/255.,77/255.)#(240/255.,242/255.,166/255.)
c_blue = (68/255.,204/255.,255/255.)
c_red = (150/255.,52/255.,132/255.)
c_orange = (242/255.,129/255.,35/255.)
c_orange_2 = (242/255.,129/255.,35/255.)

# color_list = ['r','r','b','b','g','g']
color_list = [c_orange_2,c_orange_2,c_blue,c_blue,c_green,c_green]

x_max_list = [14e-3,16e-3,16e-3,16e-3,32e-3,32e-3]
f_list = [250,200,288,150,100,110]
decay_list = [9e-3,9e-3,12e-3,12e-3,21e-3,21e-3]
spacing = [5e-3]*4+[5e-3]*2

for ii, timestamp in enumerate(timestamp_list):
    figure_name = figure_name_list[ii]
    ssro_folder = SSRO_folder_list[ii]
    color = color_list[ii]

    x, y, y_err, fit_result =  Carbon_Ramsey(timestamp=timestamp,ssro_folder =ssro_folder,frequency = f_list[ii],decay_constant = decay_list[ii])

    xticks = np.arange(0,x_max_list[ii]+1e-3,spacing[ii])
    yticks = np.linspace(0,1,3)
    xtick_rescaled = (xticks*1e3).astype(int)
    # fig,ax = plt.subplots()


    # plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False,add_txt = False,linestyle = '-',linewidth = 1,color = color)
    # fig, ax = plt.subplots(figsize = (10,3))
    fig, ax = plt.subplots(figsize = (13,3))
    res = fit_result
    fit_xvals=np.linspace(res['x'][0],res['x'][-1],1001)
    ax.plot(fit_xvals, res['fitfunc'](fit_xvals), linestyle = '-',color = color, linewidth = 2 )

    errlines = ax.errorbar(x,y,yerr = y_err, color = color,ls = '', marker = 'o',markersize = 7,markeredgecolor = color,capsize=5)
	# set(errlines, linewidth=6)

    ax.set_ylim([0,1])
    ax.set_xlim([0-0.1e-3,x_max_list[ii]+0.1e-3])
    plt.xticks(xticks)
    plt.yticks(yticks)
    ax.set_yticks(np.arange(-1,1.1,0.2), minor = True)
    ax.set_xticklabels(xtick_rescaled)
    ax.set_xlabel('Free evolution time (ms)',fontsize = 30)
    ax.set_ylabel('State fidelity',fontsize = 30)
    # ax.hlines([0.5],x[0]-1,x[-1]+1,linestyles='dotted',linewidth = 1.5)
    ax.set_ylim([0,1])

    ax.tick_params(axis='x', which='major', labelsize=30)
    ax.tick_params(axis='y', which='major', labelsize=30)
    mpl.rcParams['axes.linewidth'] = 1
    ax.tick_params('both', length=4, width=1, which='major')
    plt.savefig(os.path.join(folder, figure_name+'_'+timestamp_list[ii] + 'small.pdf'),format='pdf',bbox_inches='tight')
    plt.savefig(os.path.join(folder, figure_name+'_'+timestamp_list[ii] + 'small.png'),format='png',bbox_inches='tight')

#     if ii%2 ==0:
#         fitfunc_dict['x'] = fit_xvals
#         fitfunc_dict['y'+figure_name] =  res['fitfunc'](fit_xvals)
#         print fit_result['params_dict']
#         p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(fit_result['params_dict']['a'], fit_result['params_dict']['A'], 
#                                                     0, fit_result['params_dict']['T'],fit_result['params_dict']['n'],0 ,0)
#         # fitfunc_dict['y_dec'+figure_name] =  2.*fitfunc(fit_xvals)-1

#         ax.plot(fit_xvals,fitfunc(fit_xvals),'m')
#         # fitfunc_dict['y_dec'+figure_name] =  2.*fitfunc(fit_xvals)-1
# # pickle.dump(fitfunc_dict, open( "ramseys.p", "wb" ) )
    plt.show()