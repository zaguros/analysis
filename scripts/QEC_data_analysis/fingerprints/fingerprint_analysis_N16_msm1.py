'''
Script to analyze the dynamical decoupling data
'''
import numpy as np
import os
import analysis.lib.QEC.nuclear_spin_characterisation as SC #used for simulating FP response of spins
import matplotlib.cm as cm
from matplotlib import pyplot as plt


import fingerprint_funcs as fp_funcs; reload(fp_funcs)

def fingerprint(disp_sim_spin = True, RO = 'x'):
    
    ###################
    ## Data location ##
    ###################
    
    if RO == '-x':
        timestamp ='20140730_140911' # for the -x msmt
    elif RO == 'x':
        timestamp ='20140730_134956' # for the +x msmt
        timestamp = '20140730_184039'
    
    ssro_calib_folder = 'D:\\measuring\data\\20140730\\115839_AdwinSSRO_SSROCalibration_Hans_sil1'
    a, folder = fp_funcs.load_mult_dat(timestamp, 
                number_of_msmts = 80, 
                x_axis_step     = 0.5,
                x_axis_pts_per_msmnt= 51,
                ssro_calib_folder=ssro_calib_folder)

    #######################
    # Add simulated spins #
    #######################

    if disp_sim_spin == True:
            
            HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = 'min')
            B_Field = 304.49 
            tau_lst = np.linspace(0, 72e-6, 10000)
            Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,16,tau_lst)
            FP_signal16 = ((Mt16+1)/2)
 
    ###############
    ## Plotting ###
    ###############

    fig = a.default_fig(figsize=(35,5))
    ax = a.default_ax(fig)
    ax.set_xlim(0,40)
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
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.pdf'),
        format='pdf')
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.png'),
        format='png')

