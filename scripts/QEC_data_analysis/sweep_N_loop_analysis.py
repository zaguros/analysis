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

import analysis.scripts.QEC_data_analysis.hyperfine_params as hyperfine_params; reload(hyperfine_params) 
hf = hyperfine_params.hyperfine_params

def sweep_N_analysis(tau,N_steps):


    ###################
    # Add simulated spins #
    ###################

    total_pts = 320/N_steps+1
    pts_per_run=11

    ## Data location ##
    ssro_calib_folder = 'd:\\measuring\\data\\20140419\\111949_AdwinSSRO_SSROCalibration_Hans_sil1'

    a, folder = load_mult_dat(tau, number_of_msmts = total_pts/pts_per_run,N_steps=N_steps, ssro_calib_folder=ssro_calib_folder)


    ############
    ## Plotting ###
    ############

    fig = a.default_fig(figsize=(14,8))
    ax = a.default_ax(fig)
    # ax.set_xlim(23.4,25)
    # ax.set_xlim(0,73)
    # start, end = ax.get_xlim()
    # ax.xaxis.set_ticks(np.arange(start, end, 0.5))

    # ax.set_ylim(-0.05,1.05)
   
    ax.plot(a.sweep_pts, a.p0, '.-k', lw=0.4,label = 'data') #N = 16

    # plt.legend(loc=4)

    print folder
    plt.savefig(os.path.join(folder, 'sweep_N_analysis.pdf'),
        format='pdf')
    plt.savefig(os.path.join(folder, 'sweep_N_analysis.png'),
        format='png')


def load_mult_dat(tau, number_of_msmts, N_steps=4,ssro_calib_folder=''):
  cum_pts = 0
  print number_of_msmts
  for kk in range(number_of_msmts):
    print kk
    print str(tau)
    folder = toolbox.latest_data(contains=str(tau), older_than=None)
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



sweep_N_analysis(30412,16)







