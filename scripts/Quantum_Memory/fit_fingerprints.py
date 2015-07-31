
import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import analysis.lib.qec.nuclear_spin_characterisation as SC #used for simulating FP response of spins
import matplotlib.cm as cm
from matplotlib import pyplot as plt

import analysis.lib.qec.hyperfine_params as hyperfine_params; reload(hyperfine_params) 

hf = hyperfine_params

ZFS                 = 2.87180e9
g_factor            = 2.8025e6

def plot_sim_vs_Bx(spin_list=['C1'],Bx_list = [0],B_Field = 12, N =2):

    for ii in range(len(spin_list)):
        print spin_list
        fig = figure(ii+1,figsize=(25,6))
        ax = fig.add_subplot(111)
        # ax.title('Vary Bx for Spin '+ spin_list[ii])
        start, end = ax.get_xlim()
        # ax.xaxis.set_ticks(np.arange(start, end, 0.1))

        for b in range(len(Bx_list)):
            Bx = Bx_list[b]
            #print spin_list[ii]

            HF_par=[hf[spin_list[ii]]['par']  - hf[spin_list[ii]]['perp']*Bx/B_Field]#HF_par 	= [10e3]
            HF_perp=[hf[spin_list[ii]]['perp'] + hf[spin_list[ii]]['par']*Bx/B_Field]#HF_perp = [100e3]
            '''
            print Bx/B_Field

            print hf[spin_list[ii]]['par']
            print hf[spin_list[ii]]['perp']

            print HF_par
            print HF_perp
            '''
            print 'ii =' ,ii
            
            tau_lst = np.linspace(2e-6,20e-6,621)
            if ii == 0:
                Mt = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,N,tau_lst)
            else:
                Mt=Mt*SC.dyn_dec_signal(HF_par,HF_perp,B_Field,N,tau_lst)
            #if ii == 0:
            FP_signal = ((Mt+1)/2)
            #else:
            #    FP_signal = FP_signal*((Mt+1)/2)
            ax.plot(tau_lst*1e6, FP_signal[0,:],'.-',lw=.8,label = 'spin_'+spin_list[ii]+'_Bx_'+str(Bx))
    ax.set_xlabel('tau (us)')
#     ax.set_ylim([0.5,1.05])
    plt.legend(loc=4)

plot_sim_vs_Bx(spin_list=['C1'],Bx_list = [21],B_Field = 22, N =8)
plot_sim_vs_Bx(spin_list=['C1'],Bx_list = [7],B_Field = 22, N =8)