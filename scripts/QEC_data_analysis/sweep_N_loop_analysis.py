'''
Script to analyze the dynamical decoupling data
'''

import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import analysis.lib.QEC.nuclear_spin_characterisation as SC 
import matplotlib.cm as cm
import os
from matplotlib import pyplot as plt
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
reload(fit)
def sweep_N_analysis(tau,N_steps,plot_fit = False, do_print = True, show_guess = False):


    ###################
    # Add simulated spins #
    ###################

    total_pts = 320/N_steps+1
    pts_per_run=11

    ## Data location ##
    ssro_calib_folder = 'd:\\measuring\\data\\20140419\\111949_AdwinSSRO_SSROCalibration_Hans_sil1'

    a, folder = load_mult_dat(tau, number_of_msmts = total_pts/pts_per_run,N_steps=N_steps, ssro_calib_folder=ssro_calib_folder)
    a.p0 = a.p0*2-1

    ############
    ## Plotting ###
    ############


    ax = a.plot_results_vs_sweepparam(ret='ax')
    # ax.set_xlim(23.4,25)
    # ax.set_xlim(0,73)
    # start, end = ax.get_xlim()
    # ax.xaxis.set_ticks(np.arange(start, end, 0.5))

    ax.set_ylim(-1.05,1.05)
    
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    fit_results = []
    p0, fitfunc, fitfunc_str = common.fit_decaying_cos(1/45.,0, 1, 0, 40000)
    print fitfunc
    #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)
    ax.plot(x,y)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[3])
    # print fitfunc(np.linspace(0,x[-1],201))
    ## plot data and fit as function of total time
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)

    fit_results.append(fit_result)

    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')

def load_mult_dat(tau, number_of_msmts, N_steps=4,ssro_calib_folder=''):
  cum_pts = 0
  for kk in range(number_of_msmts):
    folder = toolbox.latest_data(contains=str(tau), older_than='140424')
    # print folder
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='measurement'+ str(kk))
    a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)
    cum_pts += a.pts
       
    pts         = 11
    N_start    = (kk+0)     * (pts-1)*N_steps 
    N_end      = (kk+1+0) * (pts-1)*N_steps
    N_list     = np.linspace(N_start, N_end, pts)

    if kk == 0:
        cum_sweep_pts = np.linspace(N_start, N_end, 11)#a.sweep_pts
        cum_p0 = a.p0
        cum_u_p0 = a.u_p0
    else:
        cum_sweep_pts = np.concatenate((cum_sweep_pts, np.linspace(N_start,N_end,11)))
        cum_p0 = np.concatenate((cum_p0, a.p0))
        cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))


  a.pts   = cum_pts
  a.sweep_pts = cum_sweep_pts
  a.p0    = cum_p0
  a.u_p0  = cum_u_p0

  return a, folder

# sweep_N_analysis(8640,4)
tau_list = [6.522,18.102,31.138,54.294,6.620,9.560,8.088,12.500,11.010,8.640,8.828,9.558,12.038,12.086,18.576,6.456,
                15.066,6.726,9.712,12.706,7.068,10.214,13.352,14.918,21.192,21.211,24.328,24.366,13.528,9.824,12.848,
                17.378,9.854,11.370,12.888,15.920,22.570,24.132,30.354,11.701,14.820,24.172,27.294,30.412,16.500]
N_steps_list = [4,4,4,4,8,8,8,8,8,4,4,4,4,4,4,8,8,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,8,8,8,8,16,16,16,16,16,16,16,16,16]
    
for jj in range(len(tau_list)):
    print int(tau_list[jj]*1e3)
    sweep_N_analysis(int(tau_list[jj]*1e3),N_steps_list[jj])







