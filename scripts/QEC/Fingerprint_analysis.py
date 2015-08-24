'''
Script to analyze the dynamical decoupling data of 21 and 22 /07
N = 4,8,16,32 and 64, contrast measurement
no Nitrogen MBI, Hermite pulses
'''
import pickle
import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import analysis.lib.QEC.nuclear_spin_characterisation as SC #used for simulating FP response of spins
# import magnettools as mt # Does not work atm because of qt lab being imported in MT
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from analysis.lib.tools import toolbox
reload(toolbox)

from analysis.scripts.QEC_data_analysis.fingerprints import fingerprint_funcs as fp_funcs; reload(fp_funcs)


def fingerprint_single(disp_sim_spin = True,n_sim_spins = 2,xrange = [0,20],tag = '', step_size = 10e-3,start = 0, tot = 5,pts = 51,start_tau = 0.52,
                   older_than = None,return_data = False, ms = 'plus'):


    ###################
    # Add simulated spins #
    ###################

    if disp_sim_spin == True:
            
            HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = ms)
            #msmp1_f from hdf5 file
            # msm1 from hdf5 file
            # ZFG g_factor from hdf5file
            B_Field = 403.555 # use magnet tools  Bz = (msp1_f**2 - msm1_f**2)/(4.*ZFS*g_factor)

            tau_lst = np.linspace(0,72e-6,10000)
            Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,32,tau_lst)
            FP_signal16 = ((Mt16+1)/2)

    ## Data location ##
    print older_than
    timestamp, ssro_calib_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = older_than,return_timestamp = True)
    print ssro_calib_folder
    start_x = start_tau + (start)* (pts-1)*step_size
    a, folder = load_mult_dat_tag(tag,older_than, number_of_msmts = tot, ssro_calib_folder=ssro_calib_folder,start = start_x, pts = pts, step_size = step_size)

    ###############
    ## Plotting ###
    ###############

    fig = a.default_fig(figsize=(10,5))
    ax = a.default_ax(fig)
    ax.set_xlim(4.9,5.1)
    ax.set_xlim(xrange)
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(xrange[0], xrange[1], (xrange[1]- xrange[0])/10.))

    ax.set_ylim(-0.05,1.05)

    y = a.p0

    ax.plot(a.sweep_pts, y, '.-k', lw=0.4,label = 'data') #N = 16
    if disp_sim_spin == True:
      colors = cm.rainbow(np.linspace(0, 1, n_sim_spins))
      for tt in range(n_sim_spins):
        ax.plot(tau_lst*1e6, FP_signal16[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])
    if False:
        tot_signal = np.ones(len(tau_lst))
        for tt in range(n_sim_spins):
          tot_signal = tot_signal * Mt16[tt,:]
        fin_signal = (tot_signal+1)/2.0
        ax.plot(tau_lst*1e6, fin_signal,':g',lw=.8,label = 'tot')


    plt.legend(loc=4)
    # plt.show(block = False)

    print folder
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.pdf'),
        format='pdf')
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.png'),
        format='png')

    if return_data == True:
      return a.sweep_pts, a.p0

def load_mult_dat_tag(tag,older_than, number_of_msmts, ssro_calib_folder='',start = 3.0, pts = 51, step_size = 10e-3):
   cum_pts = 0
   for kk in range(number_of_msmts):
       folder = toolbox.latest_data(contains = tag, older_than = older_than,folder = 'd:\measuring\data')
       a = mbi.MBIAnalysis(folder)
       a.get_sweep_pts()
       a.get_readout_results(name='measurement' + str(kk))
       a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)
       cum_pts += a.pts

       if kk == 0:
           cum_sweep_pts = np.linspace(start, start+(pts-1)*step_size,pts)
           cum_p0 = a.p0
           cum_u_p0 = a.u_p0
       else:
           cum_sweep_pts = np.concatenate((cum_sweep_pts, np.linspace(start+kk*(pts-1)*step_size, start+(pts-1)*step_size+kk*(pts-1)*step_size, pts)))
           cum_p0 = np.concatenate((cum_p0, a.p0))
           cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

   a.pts   = cum_pts
   a.sweep_pts = cum_sweep_pts
   a.p0    = cum_p0
   a.u_p0  = cum_u_p0

   return a, folder


