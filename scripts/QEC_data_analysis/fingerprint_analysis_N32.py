'''
Script to analyze the dynamical decoupling data
'''

import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import analysis.lib.QEC.nuclear_spin_characterisation as SC #used for simulating FP response of spins
# import magnettools as mt # Does not work atm because of qt lab being imported in MT
import matplotlib.cm as cm
from matplotlib import pyplot as plt

import analysis.scripts.QEC_data_analysis.hyperfine_params as hyperfine_params; reload(hyperfine_params) 
hf = hyperfine_params.hyperfine_params
hf = hyperfine_params.hyperfine_params

def fingerprint(disp_sim_spin = True):


    ###################
    # Add simulated spins #
    ###################

    if disp_sim_spin == True:
            HF_par =   [hf['C4']['par'], hf['C5']['par'], hf['C6']['par'], hf['C7']['par'], hf['C8']['par'], hf['C9']['par'], hf['C10']['par'],   hf['C11']['par']]#hf['C1']['par'],hf['C2']['par'],hf['C3']['par'],
            HF_perp =   [hf['C4']['perp'], hf['C5']['perp'], hf['C6']['perp'], hf['C7']['perp'], hf['C8']['perp'], hf['C9']['perp'], hf['C10']['perp'],   hf['C11']['perp']]#hf['C1']['perp'],hf['C2']['perp'],hf['C3']['perp'], 
            #msmp1_f from hdf5 file
            # msm1 from hdf5 file
            # ZFG g_factor from hdf5file
            B_Field = 304.12 # use magnet tools  Bz = (msp1_f**2 - msm1_f**2)/(4.*ZFS*g_factor)

            tau_lst = np.linspace(0,72e-6,10000)
            Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,32,tau_lst)
            FP_signal16 = ((Mt16+1)/2)

    ## Data location ##
    timestamps =['20140418_185913', '20140419_233953'] 
    ssro_calib_folders = ['d:\\measuring\\data\\20140419\\111949_AdwinSSRO_SSROCalibration_Hans_sil1',
                        'd:\\measuring\\data\\20140419\\123556_AdwinSSRO_SSROCalibration_Hans_sil1']
    a, folder = load_mult_dat(timestamps, number_of_msmts = [90,90], ssro_calib_folders =ssro_calib_folders)

    ############
    ## Plotting ###
    ############

    fig = a.default_fig(figsize=(200,5))
    ax = a.default_ax(fig)
    ax.set_xlim(23.4,25)
    ax.set_xlim(0,73)
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 0.5))

    ax.set_ylim(-0.05,1.05)
   
    ax.plot(a.sweep_pts, a.p0, '.-k', lw=0.4,label = 'data') #N = 16
    if disp_sim_spin == True:
      colors = cm.rainbow(np.linspace(0, 1, len(HF_par)))
      for tt in range(len(HF_par)):
        ax.plot(tau_lst*1e6, FP_signal16[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])
    plt.legend(loc=4)

    print folder
    plt.savefig(os.path.join(folder, 'fingerprint.pdf'),
        format='pdf')
    plt.savefig(os.path.join(folder, 'fingerprint.png'),
        format='png')


def load_mult_dat(timestamps, number_of_msmts, ssro_calib_folders=['']):
   cum_pts = 0
   for tt in range(len(timestamps)):
    for kk in range(number_of_msmts[tt]):
      folder = toolbox.data_from_time(timestamps[tt])
      a = mbi.MBIAnalysis(folder)
      a.get_sweep_pts()
      a.get_readout_results(name='measurement' + str(kk))
      a.get_electron_ROC(ssro_calib_folder=ssro_calib_folders[tt])
      cum_pts += a.pts

      if kk == 0 and tt==0:
        cum_sweep_pts = np.linspace(2.0, 2.5, 51)#a.sweep_pts
        cum_p0 = a.p0
        cum_u_p0 = a.u_p0
      elif tt==0:
        cum_sweep_pts = np.concatenate((cum_sweep_pts, np.linspace(2.0+kk*(50)*10e-3, 2.5+kk*(50)*10e-3, 51)))
        cum_p0 = np.concatenate((cum_p0, a.p0))
        cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))
      elif tt==1:  
        cum_sweep_pts = np.concatenate((cum_sweep_pts, np.linspace(2.0+(kk+90)*(50)*10e-3, 2.5+(kk+90)*(50)*10e-3, 51)))
        cum_p0 = np.concatenate((cum_p0, a.p0))
        cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

   a.pts   = cum_pts
   a.sweep_pts = cum_sweep_pts
   a.p0    = cum_p0
   a.u_p0  = cum_u_p0

   return a, folder


fingerprint()







