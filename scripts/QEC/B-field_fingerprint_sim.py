
import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import analysis.lib.QEC.nuclear_spin_characterisation as SC #used for simulating FP response of spins
# import magnettools as mt # Does not work atm because of qt lab being imported in MT
import matplotlib.cm as cm
from matplotlib import pyplot as plt

import analysis.scripts.QEC_data_analysis.hyperfine_params as hyperfine_params; reload(hyperfine_params) 
hf = hyperfine_params.hyperfine_params

ZFS                 = 2.877480e9
g_factor            = 2.8025e6

def plot_sim_vs_Bx(spin_list=['C1'],Bx_list = [-0.7,0,0.7],B_Field = 304.12, N =32):


	for ii in range(len(spin_list)):
		print spin_list
		fig = figure(ii+1,figsize=(10,5))
		ax = fig.add_subplot(111)
		# ax.title('Vary Bx for Spin '+ spin_list[ii])
		start, end = ax.get_xlim()
		# ax.xaxis.set_ticks(np.arange(start, end, 0.1))

		for b in range(len(Bx_list)):
			Bx = Bx_list[b]
			print spin_list[ii]

			HF_par 	= [hf[spin_list[ii]]['par']  - hf[spin_list[ii]]['perp']*Bx/B_Field]
			HF_perp = [hf[spin_list[ii]]['perp'] + hf[spin_list[ii]]['par']*Bx/B_Field]

			print Bx/B_Field

			print hf[spin_list[ii]]['par']
			print hf[spin_list[ii]]['perp']

			print HF_par
			print HF_perp

			tau_lst = np.linspace(11.8e-6,12.3e-6,5000)
			Mt = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,N,tau_lst)
			FP_signal = ((Mt+1)/2)
			ax.plot(tau_lst*1e6, FP_signal[0,:],'-',lw=.8,label = 'spin_'+spin_list[ii]+'_Bx_'+str(Bx))
		ax.set_xlabel('tau (us)')
		plt.legend(loc=4)

plot_sim_vs_Bx(spin_list=['C3'],Bx_list = [-5,-1.5,0,1.5,5],B_Field = 304.12, N =32)
# ['C1','C3','C10']