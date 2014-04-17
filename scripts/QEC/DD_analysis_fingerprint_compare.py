'''
Script to analyze the dynamical decoupling data
'''

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
import analysis.lib.qec.nuclear_spin_characterisation as SC #used for simulating FP response of spins
# import magnettools as mt # Does not work atm because of qt lab being imported in MT

from matplotlib import pyplot as plt

reload(common)

def fingerprint_compare(disp_sim_spin = True):

    ###################
    # Add simulated spins #
    ###################
    if disp_sim_spin == True:
            HF_par = [30e3, 27e3,-62.5e3]
            HF_orth =[80e3,28.5e3,132e3]
            #msmp1_f from hdf5 file
            # msm1 from hdf5 file
            # ZFG g_factor from hdf5file
            B_Field = 304.12 # use magnet tools  Bz = (msp1_f**2 - msm1_f**2)/(4.*ZFS*g_factor)

            tau_lst = np.linspace(0,30e-6,10000)
            Mt16 = SC.dyn_dec_signal(HF_par,HF_orth,B_Field,16,tau_lst)
            FP_signal16 = ((Mt16+1)/2)
            Mt32 = SC.dyn_dec_signal(HF_par,HF_orth,B_Field,32,tau_lst)
            FP_signal32 = ((Mt32+1)/2)





    ## Data location ##
    measurement_name = ['adwindata']

    timestamps_16 = ['20140407_223450', '20140407_223801', '20140407_224119',
    '20140407_224446', '20140407_231158','20140407_225215', '20140407_225614',
     '20140407_230023', '20140407_231645', '20140407_232118', '20140407_232603',
     '20140407_233057', '20140407_233605', '20140407_234654', '20140408_111324',
     '20140408_111725', '20140408_112126', '20140408_112529', '20140408_112930',
     '20140408_113614', '20140408_114015', '20140408_114416', '20140408_114818',
     '20140408_115220', '20140408_115622', '20140408_120024', '20140408_120426',
     '20140408_120825', '20140408_130753']#,

    timestamps_32 = ['20140417_110229','20140417_110529','20140417_110827','20140417_111123','20140417_111421','20140417_111716',
     '20140417_112017','20140417_112315','20140417_112613','20140417_112914','20140417_113214','20140417_113517',
     '20140417_113817','20140417_114122','20140417_114424','20140417_114727','20140417_115034','20140417_115338',
     '20140417_115644','20140417_115953','20140417_120300','20140417_120607','20140417_120915','20140417_121224',
     '20140417_121535','20140417_121847','20140417_122159','20140417_122511','20140417_122511','20140417_122824',
     '20140417_123137','20140417_123453','20140417_123808','20140417_124125','20140417_124444','20140417_124803',
     '20140417_125124','20140417_125444','20140417_125807','20140417_130129'] #with 32 pulses



    a = load_mult_dat(timestamps_16)
    b = load_mult_dat(timestamps_32)

    ############
    ## Plotting ###
    ############

    # N = 16
    fig = a.default_fig(figsize=(18,4))
    ax = a.default_ax(fig)
    ax.set_xlim(0,23)
    ax.plot(a.sweep_pts, a.p0, '.-b', lw=1,label = 'data') #N = 16
    ax.plot(tau_lst*1e6, FP_signal16[0,:] ,'--',lw=.5,label = 'spin1')
    ax.plot(tau_lst*1e6, FP_signal16[1,:] ,'--',lw=.5,label = 'spin2')
    ax.plot(tau_lst*1e6, FP_signal16[2,:] ,'--',lw=.5,label = 'spin3')
    plt.legend(loc=4)

    fig = b.default_fig(figsize=(18,4))
    ax = b.default_ax(fig)
    ax.set_xlim(0,23)
    ax.plot(b.sweep_pts, b.p0, '.-b', lw=1,label = 'data' )
    ax.plot(tau_lst*1e6, FP_signal32[0,:] ,'--',lw=.5,label = 'spin1')
    ax.plot(tau_lst*1e6, FP_signal32[1,:] ,'--',lw=.5,label = 'spin2')
    ax.plot(tau_lst*1e6, FP_signal32[2,:] ,'--',lw=.5,label = 'spin3')
    plt.legend(loc=4)

    fig = b.default_fig(figsize=(18,4))
    ax = b.default_ax(fig)
    ax.set_xlim(0,23)
    ax.plot(a.sweep_pts, a.p0, '.-b', lw=1,label = 'N=16')
    ax.plot(b.sweep_pts, b.p0, '.-r', lw=1,label = 'N=32')
    plt.legend(loc=4)

    # plt.savefig(os.path.join(folder, 'combined_result.pdf'),
    # format='pdf')
    # plt.savefig(os.path.join(folder, 'combined_result.png'),
    # format='png')


def load_mult_dat(timestamp):
   cum_pts = 0
   for kk in range(len(timestamp)):
       folder = toolbox.data_from_time(timestamp[kk])
       a = mbi.MBIAnalysis(folder)
       a.get_sweep_pts()
       a.get_readout_results(name='adwindata')
       a.get_electron_ROC()
       cum_pts += a.pts

       if kk == 0:
           cum_sweep_pts = a.sweep_pts
           cum_p0 = a.p0
           cum_u_p0 = a.u_p0
       else:
           cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
           cum_p0 = np.concatenate((cum_p0, a.p0))
           cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

   a.pts   = cum_pts
   a.sweep_pts = cum_sweep_pts
   a.p0    = cum_p0
   a.u_p0  = cum_u_p0

   return a



fingerprint_compare()







