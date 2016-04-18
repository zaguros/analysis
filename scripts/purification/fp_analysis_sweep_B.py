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


###################
## Data location ##
###################
timestamps = {}
#uber stamps!
timestamps['min'] = {'N8' : '20160107_173313', 'N16' : '20160107_201229','N32' : '20160107_222009','N64' : '20160108_004236'}
# timestamps['plus'] = {'N8' : '20160112_174142', 'N16' : '20160112_192510','N32' : '20160112_211902','N64' : '20160112_234557'}
timestamps['plus'] = {'N8' : '20160116_140951', 'N16' : '20160116_161812','N32' : '20160116_181552','N64' : '20160116_202552'}

#crappy stamps
# timestamps['min'] = {'N8' : '20160107_173313', 'N16' : '20160107_201229','N32' : '20160107_222009','N64' : '20160108_004236'}
# timestamps['plus'] = {'N8' : '20160110_121238', 'N16' : '20160110_143232','N32' : '20160110_170758','N64' : '20160110_202511'}
ssro_calib_folder = 'd:\\measuring\\data\\20160107\\172632_AdwinSSRO_SSROCalibration_Pippin_SIL1'	

def load_data(N = [8], el_trans = 'min'):
	a = []
	folder = []
	for i in range(len(N)):
		print i
		# load data removing append for memory error!
		print 'loading data'
		a_temp, folder_temp = fp_funcs.load_mult_dat(timestamps[el_trans]['N'+str(N[i])], 
				number_of_msmts = 100,
				x_axis_step     = 0.1,
				x_axis_start    = 3.5,
				x_axis_pts_per_msmnt= 51,
				ssro_calib_folder=ssro_calib_folder)
		a.append(a_temp)
		folder.append(folder_temp)
		print 'data N' +str(N[i]) + 'for el_trans ' + el_trans + ' loaded'
	print 'All data for the specified N loaded via the timestamps in fp_ls'
	return a, folder

def fingerprint(a = None, folder = None, disp_sim_spin = True, N = [8], el_trans = 'min', 
	HF_perp = None, HF_par = None, B_list = [418]):

	# allowed params:
	# el_trans = ['min', 'plus']
	# N = 8, 16, 32, 64
	if (a == None) | (folder == None):
		print 'Folder path or data (a) missing fool!'
		return

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

	
	print 'N = ' + str(N)

	for ii in range(len(N)):

		# bad way of dealing with complete loaded data
		if  N[ii] == 8:
			i = 0
		elif N[ii] == 16:
			i = 1
		elif N[ii] == 32:
			i = 2
		elif N[ii] == 64:
			i = 3
		else:
			print 'N is not a standard value (8, 16, 32, 64)'
			return

		print i
		
		# print 'loading data'
		# a, folder = fp_funcs.load_mult_dat(timestamps[el_trans]['N'+str(N[i])], 
		#             number_of_msmts = 100,
		#             x_axis_step     = 0.1,
		#             x_axis_start    = 3.5,
		#             x_axis_pts_per_msmnt= 51,
		#             ssro_calib_folder=ssro_calib_folder)

		# print 'data N' +str(N[i]) + ' loaded'
		##########################
		### 	 plot data ######
		#########################
		
		#lw default 0.4
		fig = a[i].default_fig(figsize=(35,5))
		ax = a[i].default_ax(fig)
		ax.set_xlim(3.5,13.5)
		start, end = ax.get_xlim()
		ax.xaxis.set_ticks(np.arange(start, end, 0.5))
		ax.set_ylim(-0.05,1.05)
		ax.plot(a[i].sweep_pts, a[i].p0, '.-k', lw=0.4,label = 'data')


		#######################
		# Add simulated spins #
		#######################

		if disp_sim_spin == True:
			print 'Starting Simulation for N = ' + str(N[ii]) + ' on transition ' + str(el_trans) 
			colors = cm.rainbow(np.linspace(0, 1, len(HF_par)+len(B_list)))
			for jj,B in enumerate(B_list):	
				tau_lst = np.linspace(3.5e-6, 22.5e-6, 5000)
				Mt16 = SC.dyn_dec_signal(HFs_par = HF_par, HFs_orth = HF_perp,
					B_field = B, N = N[ii], tau = tau_lst)
				FP_signal16 = ((Mt16+1)/2)
				
				# colors = ['m', 'b', 'r', 'g', 'c']
				for tt in range(len(HF_par)):
				  ax.plot(tau_lst*1e6, FP_signal16[tt,:] ,'-',lw=1,label =  str(tt + 1) + ': B = ' + str(B) + 
				  	': HF_par = ' +str(HF_par[tt]) + '; HF_perp = ' +str(HF_perp[tt]), color = colors[tt+jj])
				

				print folder[i]
			plt.legend(loc=4)

			# plt.savefig(os.path.join(folder[i], str(disp_sim_spin)+'fingerprint.pdf'),
			#     format='pdf')
			# plt.savefig(os.path.join(folder[i], str(disp_sim_spin)+'fingerprint.png'),
			#     format='png')

