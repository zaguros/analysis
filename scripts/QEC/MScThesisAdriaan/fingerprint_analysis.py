'''
Script to analyze the dynamical decoupling data
Copy to be self contained and used for generating figure of MSc thesis of M.A.Rol

'''
import numpy as np
import os
import analysis.lib.QEC.nuclear_spin_characterisation as SC #used for simulating FP response of spins
import matplotlib.cm as cm
from matplotlib import pyplot as plt

import fingerprint_funcs_MAR as fp_funcs; reload(fp_funcs)

def fingerprint(disp_sim_spin = True, n_spins_to_disp = 'all' ,N = 16, xlims = [2,10],B_Field =304.12,disp_total_sig=False,figsize= (25,5),fontsize = 10,showlegend = True,title ='Fingerprint of Hans Sil 01 at 304.12G for N =32 pulses'):

    ###################
    ## Data location ##
    ###################

    if N == 16:
        timestamp ='20140419_005744' # for the -x msmt
        ssro_calib_folder = '//Users//Adriaan//Documents//teamdiamond//data//20140419//111949_AdwinSSRO_SSROCalibration_Hans_sil1'
        a, folder = fp_funcs.load_mult_dat(timestamp,
                    number_of_msmts = 140,
                    x_axis_step     = 0.5,
                    x_axis_pts_per_msmnt= 51,
                    ssro_calib_folder=ssro_calib_folder)

    elif N == 32:
        timestamps =['20140418_185913', '20140419_233953']
        ssro_calib_folders = ['//Users//Adriaan//Documents//teamdiamond//data//20140419//111949_AdwinSSRO_SSROCalibration_Hans_sil1',
                '//Users//Adriaan//Documents//teamdiamond//data//20140419//123556_AdwinSSRO_SSROCalibration_Hans_sil1']
        a, folder = fp_funcs.load_mult_dat_2(timestamps = timestamps, number_of_msmts = [90,90], ssro_calib_folders =ssro_calib_folders)

    else:
        print 'No data known for N = %s'  %N
        return


    #######################
    # Add simulated spins #
    #######################

    HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = 'plus')
    tau_lst = np.linspace(0, 72e-6, 10000)
    Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,N,tau_lst)
    FP_signal16 = ((Mt16+1)/2)

    ###############
    ## Plotting ###
    ###############

    fig = a.default_fig(figsize=figsize)
    ax = a.default_ax(fig)
    ax.set_xlim(xlims)
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 0.5))
    ax.set_ylim(-0.05,1.05)
    ax.plot(a.sweep_pts, a.p0, 'k', lw=0.4,label = 'data') #N = 16


    if disp_sim_spin == True:
      # colors = cm.rainbow(np.linspace(0, 1, len(HF_par)))
      if n_spins_to_disp =='all':
        n_spins_to_disp = len(HF_par)
      for tt in range(n_spins_to_disp):
        ax.plot(tau_lst*1e6, FP_signal16[tt,:] ,'-',lw=.2,label = 'spin' + str(tt+1))#, color = colors[tt])
    if disp_total_sig==True:
        tot_signal = np.ones(len(tau_lst))
        for tt in range(len(HF_par)):
          tot_signal = tot_signal * Mt16[tt,:]
        fin_signal = (tot_signal+1)/2.0
        ax.plot(tau_lst*1e6, fin_signal,':b',lw=.8,label = 'tot')

    ax.set_title(title)
    ax.set_xlabel(r'$\tau$ ($\mu$s) ',fontsize =fontsize)
    ax.set_ylabel(r'$F$',fontsize =fontsize) #Maybe I want to have contrast instead <X> from -1 to +1

    if showlegend == True:
        plt.legend(loc=4)




    print 'Figures saved in: %s' %folder
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.pdf',),
        format='pdf',bbox_inches='tight')
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.png'),
        format='png',bbox_inches='tight')


return a





