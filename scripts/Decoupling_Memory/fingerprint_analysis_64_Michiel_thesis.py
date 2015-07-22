'''
Script to analyze the dynamical decoupling data
'''

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import analysis.lib.QEC.nuclear_spin_characterisation as SC #used for simulating FP response of spins
# import magnettools as mt # Does not work atm because of qt lab being imported in MT
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from analysis.lib.tools import toolbox

import fingerprint_funcs as fp_funcs; reload(fp_funcs)

def fingerprint(disp_sim_spin = True,n_sim_spins = 5,N=32,xrange = [0,20],carbons=[1,2,3,4,5]):



    ###################
    # Add simulated spins #
    ###################

    if disp_sim_spin == True:
            HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = 'min',carbon_spins=carbons)
            #msmp1_f from hdf5 file
            # msm1 from hdf5 file
            # ZFG g_factor from hdf5file
            B_Field = 403.555 # use magnet tools  Bz = (msp1_f**2 - msm1_f**2)/(4.*ZFS*g_factor)
            tau_L = 1/(1.0705e3*B_Field)*1e6
            print tau_L

            tau_lst = np.linspace(0,72e-6,10000)
            Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,N,tau_lst)
            FP_signal16 = ((Mt16+1)/2)

    ## Data location ##
    save_folder_str = 'D:/Dropbox/QEC LT/Decoupling memory/Fingerprints/N=' + str(N)
    if N==32:
      timestamp ='20141016_205842'
      ssro_calib_folder = 'D:\\measuring\\data\\20141016\\150451_AdwinSSRO_SSROCalibration_111_1_sil18'
      a, folder = load_mult_dat(timestamp, number_of_msmts = 40, ssro_calib_folder=ssro_calib_folder)

    elif N==64:
      timestamp ='20141016_234534'
      ssro_calib_folder = 'D:\\measuring\\data\\20141016\\150451_AdwinSSRO_SSROCalibration_111_1_sil18'
      a, folder = load_mult_dat(timestamp, number_of_msmts = 50, ssro_calib_folder=ssro_calib_folder)

    else:
      raise Exception('Only N=32 or N=64')
    print tau_L

    carbon_relabel = {}
    carbon_relabel['1'] = '2'
    carbon_relabel['2'] = '1'
    carbon_relabel['3'] = '4'
    carbon_relabel['5'] = '3'
    carbon_relabel['6'] = '5'
    ###############
    ## Plotting ###
    ###############
    first_tau = int(np.ceil(xrange[0]/tau_L))*tau_L
    last_tau = int(np.ceil(xrange[1]/tau_L))*tau_L
    last_tau_nr = int(np.ceil(xrange[1]/tau_L))
    fig = a.default_fig(figsize=(16,4))
    ax = a.default_ax(fig)
    # ax.set_xlim(15.0,15.5)
    ax.set_title('')
    ax.set_xlabel(r'$\tau$ $(\mu s)$',fontsize = 20)
    ax.set_ylabel(r'$\langle$X$\rangle$',fontsize = 20)
    ax.set_xlim(xrange)
    #ax.hlines(0.,xrange[0],xrange[1],linestyles='dotted',linewidth = 1.5)
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(np.ceil(xrange[0]), xrange[1], 1))
    ax.vlines(first_tau,0.1,2,linestyles='dotted',linewidth=1.5)
    for tau_value in np.arange(first_tau+tau_L, last_tau, tau_L):
      ax.vlines(tau_value,-2,2,linestyles='dotted',linewidth=1.5)
    ax.tick_params(axis='x', which='major', labelsize=15)
    ax.tick_params(axis='y', which='major', labelsize=15)
    #ax.set_xticklabels(range(1,last_tau_nr+1))
    ax.set_ylim(-1.05,1.05)
    handles = []
    ax.plot(a.sweep_pts, a.p0*2-1, '.-k', lw=0.4) #N = 16
    if disp_sim_spin == True:
      colors = cm.rainbow(np.linspace(0, 1, n_sim_spins))
      for tt, carbon_nr in enumerate(carbons):
        h = ax.plot(tau_lst*1e6, FP_signal16[tt,:]*2-1 ,'-',lw=.8,label = 'C' + carbon_relabel[str(carbon_nr)], color = colors[tt])
        handles.append(h)
    if False:
        tot_signal = np.ones(len(tau_lst))
        for tt in range(n_sim_spins):
          tot_signal = tot_signal * Mt16[tt,:]
        fin_signal = (tot_signal+1)/2.0
        ax.plot(tau_lst*1e6, fin_signal,':g',lw=.8,label = 'tot')

    #handles, labels = ax.get_legend_handles_labels()
    #order= [i[0] for i in sorted(enumerate(labels), key=lambda x:x[1])]
    #handles = handles[order]
    #labels = labels[order]
    # sort both labels and handles by labels
    
    #labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    #ax.legend(handles, labels)
    plt.legend(loc='lower left',frameon=False,fontsize=15,columnspacing=0.2,handletextpad=0.0)

    print folder
    plt.savefig(save_folder_str+'fingerprint.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(save_folder_str+'fingerprint.png', format='png')

def load_mult_dat(timestamp, number_of_msmts, ssro_calib_folder=''):
   cum_pts = 0
   for kk in range(number_of_msmts):
       folder = toolbox.data_from_time(timestamp,folder = 'd:\measuring\data')
       a = mbi.MBIAnalysis(folder)
       a.get_sweep_pts()
       a.get_readout_results(name='measurement' + str(kk))
       a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)
       cum_pts += a.pts

       if kk == 0:
           cum_sweep_pts = np.linspace(2.0, 2.2, 51)
           cum_p0 = a.p0
           cum_u_p0 = a.u_p0
       else:
           cum_sweep_pts = np.concatenate((cum_sweep_pts, np.linspace(2.0+kk*(50)*4e-3, 2.2+kk*(50)*4e-3, 51)))
           cum_p0 = np.concatenate((cum_p0, a.p0))
           cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

   a.pts   = cum_pts
   a.sweep_pts = cum_sweep_pts
   a.p0    = cum_p0
   a.u_p0  = cum_u_p0

   return a, folder
