'''
Script to analyze the dynamical decoupling data from January 2016

Data covered:

N = 8, 16, 32, 64
Electron transition 0 to -1 (min) and 0 to +1 (plus), this flips the sign on the parallel hyperfine component

'''

import numpy as np
import os
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy.misc import imread
from analysis.lib.tools import toolbox

import fp_purification_funcs as fp_funcs; reload(fp_funcs)

# import analysis.lib.QEC.nuclear_spin_characterisation as SC #used for simulating FP response of spins
import nuclear_spin_char as SC;  reload(SC)


def fingerprint(disp_sim_spin = True, N = [8], el_trans = 'min', HF_perp = None, HF_par = None):

	# allowed params:
	# el_trans = ['min', 'plus']
	# N = 8, 16, 32, 64

	###################
	## Data location ##
	###################

	timestamps = {}
	timestamps['min'] = {'N8' : '20160229_114914', 
						'N16' : '20160229_133036',
						'N32' : '20160229_152201',
						'N64' : '20160229_174524'}
	timestamps['plus'] = {'N8' : '20160110_121238', 'N16' : '20160110_143232','N32' : '20160110_170758','N64' : '20160110_202511'}
	ssro_calib_folder = 'd:\\measuring\\data\\20160107\\172632_AdwinSSRO_SSROCalibration_Pippin_SIL1'	


	### Load hyperfine params	
	if disp_sim_spin == True:	
		if (HF_perp == None) & (HF_par== None):
			HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = el_trans, NV = 'Pippin')
		elif el_trans == 'min':
			HF_par =  [x * (-1) for x in HF_par]

		# security check could be removed
		if len(HF_perp) == len(HF_par):
			pass
		else:
			print 'Unequal amount of Parallel and Perpendicular HF parameters'

		print 'HF_perp = ' + str(HF_perp)
		print 'HF_par = ' + str(HF_par)

	else:
		print 'No HF simulation'

	a = []
	folder = []
	print N
	for i in range(len(N)):
		print i
		# load data removing append for memory error!
		print 'loading data'
		a, folder = fp_funcs.load_mult_dat(timestamps[el_trans]['N'+str(N[i])], 
              number_of_msmts = 100,
              x_axis_step     = 0.1,
              x_axis_start    = 3.5,
              x_axis_pts_per_msmnt= 51,
              ssro_calib_folder=ssro_calib_folder)

		print 'data N' +str(N[i]) + ' loaded'
		##########################
		### 	 plot data ######
		#########################
		
		fig = a.default_fig(figsize=(35,5))
		ax = a.default_ax(fig)
		ax.set_xlim(3.5,13.5)
		start, end = ax.get_xlim()
		ax.xaxis.set_ticks(np.arange(start, end, 0.5))
		ax.set_ylim(-0.05,1.05)
		ax.plot(a.sweep_pts, a.p0, '.-k', lw=0.4,label = 'data')


		#######################
		# Add simulated spins #
		#######################

		if disp_sim_spin == True:
			print 'Starting Simulation for N = ' + str(N[i]) + ' on transtion ' + str(el_trans) 
			B_Field = 417.22
			tau_lst = np.linspace(0, 72e-6, 10000)
			Mt16 = SC.dyn_dec_signal(HFs_par = HF_par, HFs_orth = HF_perp,
				B_field = B_Field, N = N[i], tau = tau_lst)
			FP_signal16 = ((Mt16+1)/2)
			
			colors = cm.rainbow(np.linspace(0, 1, len(HF_par)))
			for tt in range(len(HF_par)):
			  ax.plot(tau_lst*1e6, FP_signal16[tt,:] ,'-',lw=1.5,label = str(tt + 1) + ': HF_par = ' +str(HF_par[tt]) + '; HF_perp = ' +str(HF_perp[tt]), color = colors[tt])
			plt.legend(loc=4)

			print folder
			plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.pdf'),
			    format='pdf')
			plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.png'),
			    format='png')

