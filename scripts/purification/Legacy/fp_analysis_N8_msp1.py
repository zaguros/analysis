'''
Module to analyze the dynamical decoupling data for N=8, ms=+1 and ~300 Gauss
'''

import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
# import analysis.lib.QEC.nuclear_spin_characterisation as SC
import matplotlib.cm as cm
import os
from matplotlib import pyplot as plt

# import analysis.scripts.Purification.hyperfine_params as module_hyperfine_params; reload(module_hyperfine_params) 
import hyperfine_params as module_hyperfine_params; reload(module_hyperfine_params) 
hf = module_hyperfine_params.hyperfine_params_hans_SIL1_msm1

import fp_purification_funcs as fp_funcs; reload(fp_funcs)
import nuclear_spin_char as SC;  reload(SC)


def fingerprint(disp_sim_spin = True):
    

    ###################
    ## Data location ##
    ###################
    
    timestamp ='20160110_121238'
    ssro_calib_folder = 'd:\\measuring\\data\\20160107\\172632_AdwinSSRO_SSROCalibration_Pippin_SIL1'
    a, folder = fp_funcs.load_mult_dat(timestamp, 
              number_of_msmts = 150,
              x_axis_step     = 0.1,
              x_axis_start    = 3.5,
              x_axis_pts_per_msmnt= 51,
              ssro_calib_folder=ssro_calib_folder)

    #######################
    # Add simulated spins #
    #######################

    if disp_sim_spin == True:
        print 'Starting Simulation'  
        HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = 'plus', NV = 'Pippin')
        print 'HF_perp = ' + str(HF_perp)
        print 'HF_par = ' + str(HF_par)
        B_Field = 417.268
        tau_lst = np.linspace(0, 20e-6, 5000)
        Mt16 = SC.dyn_dec_signal(HFs_par = HF_par, HFs_orth = HF_perp,
          B_field = B_Field, N = 8, tau = tau_lst)
        FP_signal16 = ((Mt16+1)/2)

    
    ############
    ## Plotting ###
    ############

    fig = a.default_fig(figsize=(35,5))
    ax = a.default_ax(fig)
    ax.set_xlim(a.sweep_pts[0],a.sweep_pts[-1])
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 0.5))

    ax.set_ylim(-0.05,1.05)
   
    ax.plot(a.sweep_pts, a.p0, '.-k', lw=0.4,label = 'data') #N = 16
    
    if disp_sim_spin == True:
      colors = cm.rainbow(np.linspace(0, 1, len(HF_par)))
      for tt in range(len(HF_par)):
        ax.plot(tau_lst*1e6, FP_signal16[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])
    if False:
        tot_signal = np.ones(len(tau_lst))
        for tt in range(len(HF_par)):
          tot_signal = tot_signal * Mt16[tt,:]
        fin_signal = (tot_signal+1)/2.0   
        ax.plot(tau_lst*1e6, fin_signal,':g',lw=.8,label = 'tot')
    

    plt.legend(loc=4)

    print folder
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint_short.pdf'),
        format='pdf')
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint_short.png'),
        format='png')

    