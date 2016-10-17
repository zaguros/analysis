'''
Script to analyze the dynamical decoupling data of 21 and 22 /07
N = 4,8,16,32 and 64, contrast measurement
no Nitrogen MBI, Hermite pulses
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
reload(toolbox)

from analysis.scripts.QEC_data_analysis.fingerprints import fingerprint_funcs as fp_funcs; reload(fp_funcs)


def fingerprint_contrast(disp_sim_spin = True,n_sim_spins = 2,pts = 51,xrange = [0,20],tag_p = '',tag_n = '', older_than = None,return_data = False):


    ###################
    # Add simulated spins #
    ###################

    if disp_sim_spin == True:
            
            HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = 'min')
            #msmp1_f from hdf5 file
            # msm1 from hdf5 file
            # ZFG g_factor from hdf5file
            B_Field = 403.555 # use magnet tools  Bz = (msp1_f**2 - msm1_f**2)/(4.*ZFS*g_factor)

            tau_lst = np.linspace(0,72e-6,10000)
            Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,32,tau_lst)
            FP_signal16 = ((Mt16+1)/2)

    ## Data location ##
    
    ssro_calib_folder = folder = toolbox.latest_data(contains = 'AdwinSSRO_SSROCalibration_111_1_sil18', older_than = older_than,folder = 'd:\measuring\data')
    a, folder = load_mult_dat_contrast(tag_p,older_than, number_of_msmts = 50,pts = pts, ssro_calib_folder=ssro_calib_folder)
    b, folder_b = load_mult_dat_contrast(tag_n,older_than, number_of_msmts = 50,pts = pts, ssro_calib_folder=ssro_calib_folder)

    ###############
    ## Plotting ###
    ###############

    fig = a.default_fig(figsize=(35,5))
    ax = a.default_ax(fig)
    ax.set_xlim(4.9,5.1)
    ax.set_xlim(xrange)
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(xrange[0], xrange[1], (xrange[1]- xrange[0])/10.))

    ax.set_ylim(-1.05,1.05)

    y = (a.p0-b.p0)

    print a.sweep_pts

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
    plt.show(block = False)

    print folder
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint_contrast.pdf'),
        format='pdf')
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint_contrast.png'),
        format='png')

    if return_data == True:
      return a.sweep_pts, (y+1)/2.

def fingerprint_single(disp_sim_spin = True,n_sim_spins = 2,xrange = [0,20],tag = '', step_size = 10e-3,
                   older_than = None,return_data = False):


    ###################
    # Add simulated spins #
    ###################

    if disp_sim_spin == True:
            
            HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = 'min')
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
    a, folder = load_mult_dat_tag(tag,older_than, number_of_msmts = 18, ssro_calib_folder=ssro_calib_folder)

    ###############
    ## Plotting ###
    ###############

    fig = a.default_fig(figsize=(35,5))
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
    plt.show(block = False)

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



def fingerprint(disp_sim_spin = True,n_sim_spins = 2,step_size= 10e-3, xrange = [0,20],return_data = False,    timestamp ='20141016_205842',
    ssro_calib_folder = 'D:\\measuring\\data\\20141016\\150451_AdwinSSRO_SSROCalibration_111_1_sil18'):


    ###################
    # Add simulated spins #
    ###################

    if disp_sim_spin == True:
            
            HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = 'min')
            #msmp1_f from hdf5 file
            # msm1 from hdf5 file
            # ZFG g_factor from hdf5file
            B_Field = 403.555 # use magnet tools  Bz = (msp1_f**2 - msm1_f**2)/(4.*ZFS*g_factor)

            tau_lst = np.linspace(0,72e-6,10000)
            Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,32,tau_lst)
            FP_signal16 = ((Mt16+1)/2)

    ## Data location ##

    a, folder = load_mult_dat(timestamp,step_size = step_size, number_of_msmts = 40, ssro_calib_folder=ssro_calib_folder)

    ###############
    ## Plotting ###
    ###############

    fig = a.default_fig(figsize=(35,5))
    ax = a.default_ax(fig)
    ax.set_xlim(4.9,5.1)
    ax.set_xlim(xrange)
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(xrange[0], xrange[1], (xrange[1]- xrange[0])/10.))

    ax.set_ylim(-0.05,1.05)
    ax.plot(a.sweep_pts, a.p0, '.-k', lw=0.4,label = 'data') #N = 16
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

    print folder
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.pdf'),
        format='pdf')
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.png'),
        format='png')

    if return_data == True:
      return a.sweep_pts, a.p0

def load_mult_dat(timestamp, number_of_msmts,older_than = '', ssro_calib_folder='',start = 2.0, pts = 51, step_size = 10e-3):
   cum_pts = 0
   for kk in range(number_of_msmts):
       folder = toolbox.data_from_time(timestamp,folder = 'd:\measuring\data')
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



tag_p_list = [
'mx111_1_sil18_64pulses',
            'mx111_1_sil18_32pulses',
            'mx111_1_sil18_16pulses',
            'mx111_1_sil18_8pulses',
            'mx111_1_sil18_4pulses'
                  ]
tag_n_list = [
            'px111_1_sil18_64pulses',
            'px111_1_sil18_32pulses',
            'px111_1_sil18_16pulses',
            'px111_1_sil18_8pulses',
            'px111_1_sil18_4pulses'
                  ]                  

older_than = '20150724_090000'

# for ii, tag_p in enumerate(tag_p_list):
#   fingerprint_contrast(disp_sim_spin = True,n_sim_spins = 0,xrange = [2,27],tag_p =tag_p,tag_n = tag_n_list[ii],older_than = older_than,return_data = False)
  # fingerprint_single(disp_sim_spin = True,n_sim_spins = 0,xrange = [2,27],tag =tag_p,older_than = older_than,return_data = False)

# fingerprint_contrast(disp_sim_spin = True,n_sim_spins = 0,xrange = [2,14],tag_p ='mx111_1_sil18_128pulses',
#                     tag_n = 'px111_1_sil18_128pulses',pts = 21,older_than = older_than,return_data = False)




H_x, H_y = fingerprint_single(disp_sim_spin = True,n_sim_spins = 0,xrange = [2,27],tag ='px111_1_sil18_32pulses',
                          older_than = '20150722_170000',return_data = True)

S_x, S_y = fingerprint_single(disp_sim_spin = True,n_sim_spins = 0,xrange = [2,27],tag ='px111_1_sil18_32pulses',
                          older_than = '20150723_020000',return_data = True)


fig, ax = plt.subplots(figsize=(20,5))
# ax.set_xlim(4.9,5.1)
ax.set_xlim([0,15])
ax.plot(H_x,H_y,'r',label = 'Hermite')
ax.plot(S_x,S_y,'b',label = 'Square')
plt.title('32 pulses')
plt.legend()
plt.show(block = False)

H_x, H_y = fingerprint_single(disp_sim_spin = True,n_sim_spins = 0,xrange = [2,27],tag ='px111_1_sil18_64pulses',
                          older_than = '20150722_170000',return_data = True)

S_x, S_y = fingerprint_single(disp_sim_spin = True,n_sim_spins = 0,xrange = [2,27],tag ='px111_1_sil18_64pulses',
                          older_than = '20150723_020000',return_data = True)


fig, ax = plt.subplots(figsize=(20,5))
# ax.set_xlim(4.9,5.1)
ax.set_xlim([0,15])
ax.plot(H_x,H_y,'r',label = 'Hermite')
ax.plot(S_x,S_y,'b',label = 'Square')
plt.title('64 pulses')
plt.legend()
plt.show(block = False)

# #### old data
# # N = 16
# timestamp_16 ='20141016_182539'
# ssro_calib_folder_16 = 'D:\\measuring\\data\\20141016\\150451_AdwinSSRO_SSROCalibration_111_1_sil18'

# x_old, y_old = fingerprint(disp_sim_spin = True,n_sim_spins = 0,xrange = [0,20],timestamp = timestamp_16,ssro_calib_folder = ssro_calib_folder_16,return_data = True)
# x_new, y_new = fingerprint_contrast(disp_sim_spin = True,n_sim_spins = 0,xrange = [0,20],tag_p = 'mx111_1_sil18_16pulses',tag_n = 'px111_1_sil18_16pulses', 
#                 older_than = None,return_data = True)


# fig, ax = plt.subplots(figsize=(20,5))
# ax.set_xlim(4.9,5.1)
# ax.set_xlim([0,27])
# ax.plot(x_old,y_old,'r',label = 'old')
# ax.plot(x_new,y_new,'b',label = 'new')
# plt.legend()
# plt.show(block = False)

# timestamp_32 ='20141016_205842'
# ssro_calib_folder_32 = 'D:\\measuring\\data\\20141016\\150451_AdwinSSRO_SSROCalibration_111_1_sil18'

# # fingerprint()

# x_old, y_old = fingerprint(disp_sim_spin = True,n_sim_spins = 0,xrange = [0,20],timestamp = timestamp_32,ssro_calib_folder = ssro_calib_folder_32,return_data = True)
# x_new, y_new = fingerprint_contrast(disp_sim_spin = True,n_sim_spins = 0,xrange = [0,20],tag_p = 'mx111_1_sil18_32pulses',tag_n = 'px111_1_sil18_32pulses', 
#                 older_than = None,return_data = True)


# fig, ax = plt.subplots(figsize=(20,5))
# ax.set_xlim(4.9,5.1)
# ax.set_xlim([0,27])
# ax.plot(x_old,y_old,'r',label = 'old')
# ax.plot(x_new,y_new,'b',label = 'new')
# plt.legend()
# plt.show(block = False)


# timestamp_64 ='20141016_234534'
# ssro_calib_folder_64 = 'D:\\measuring\\data\\20141016\\150451_AdwinSSRO_SSROCalibration_111_1_sil18'
    

# x_old, y_old = fingerprint(disp_sim_spin = True,n_sim_spins = 0,xrange = [0,20],timestamp = timestamp_64,ssro_calib_folder = ssro_calib_folder_64,return_data = True,step_size= 4e-3)
# x_new, y_new = fingerprint_contrast(disp_sim_spin = True,n_sim_spins = 0,xrange = [0,20],tag_p = 'mx111_1_sil18_64pulses',tag_n = 'px111_1_sil18_64pulses', 
#                 older_than = None,return_data = True)


# fig, ax = plt.subplots(figsize=(20,5))
# ax.set_xlim(4.9,5.1)
# ax.set_xlim([0,27])
# ax.plot(x_old,y_old,'r',label = 'old')
# ax.plot(x_new,y_new,'b',label = 'new')
# plt.legend()
# plt.show(block = False)
