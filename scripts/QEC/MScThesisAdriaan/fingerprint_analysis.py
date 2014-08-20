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

def fingerprint(disp_sim_spin = True, n_spins_to_disp = 'all' ,N = 16, xlims = [2,10],B_Field =304.12,disp_total_sig=False,
        figsize= (25,5), fontsize = 10, showlegend = True, plot_contrast = False,
        title ='Fingerprint',ms = 0.1, tau_larmor_axis = False):
    plt.rc('font', size=fontsize)

    ###################
    ## Data location ##
    ###################
    gamma_C = 1.0705e3
    tau_larmor = 1/(gamma_C*B_Field)

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
    if tau_larmor_axis == False:
        ax.xaxis.set_ticks(np.arange(start, end+1e-6, 1))
        print np.arange(start, end+1e-6, 1)
        ax.set_xlabel(r'$\tau$ ($\mu$s) ',fontsize =fontsize)
    else:

        xticks = np.arange((tau_larmor*1e6),3*(tau_larmor*1e6)+1e-6,tau_larmor*1e6*.5)
        ax.xaxis.set_ticks(xticks)
        xticklabels = [1,1.5,2,2.5,3]
        print plt.gca()
        ax.set_xlabel(r'$\tau /\tau _L$ ',fontsize =fontsize)
        ax.set_xticklabels(xticklabels)
        print xticks

    if plot_contrast == True:
        a.p0= ((a.p0.reshape(-1))-0.5)*2
        for tt in range(n_spins_to_disp):
            FP_signal16[tt,:] = ((FP_signal16[tt,:].reshape(-1))-0.5)*2
        ax.set_ylim(-1.05,1.05)
        ax.yaxis.set_ticks( [-1,-0.5,0,0.5,1])
        ax.set_ylabel(r'$\langle X \rangle $',fontsize =fontsize)
    else:
        ax.set_ylim(-0.05,1.05)
        ax.yaxis.set_ticks( [0,0.5,1])
        ax.set_ylabel(r'$F$',fontsize =fontsize)


    ax.plot(a.sweep_pts, a.p0, 'k.-',ms = ms, lw=0.4,label = 'data') #N = 16


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

    if showlegend == True:
        plt.legend(loc=4)




    print 'Figures saved in: %s' %folder
    plt.savefig(os.path.join(folder, 'fingerprint'+str(N)+'.pdf',),
        format='pdf',bbox_inches='tight')









