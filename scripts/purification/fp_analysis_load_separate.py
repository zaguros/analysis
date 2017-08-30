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



#### Pippin sil 3 timestamps

# timestamps['min'] = 	{'N8' : ['20170228_184913'],
# 					 	'N16' : ['20170228_193016'],
# 					 	'N32' : ['20170228_201540'],
# 						'N64' : ['20170228_211001']}
# timestamps['plus'] = 	{'N8' : ['20170228_231052'],
# 					 	'N16' : ['20170301_001356'],
# 					 	'N32' : ['20170301_012432'],
# 						'N64' : ['20170301_025042']}
# ssro_calib_folder = 'd:\\measuring\\data\\20170228\\172427_AdwinSSRO_SSROCalibration_Pippin_SIL3'


####################################
#### Pippin SIL 2 timestamps
####################################


# timestamps['min'] = 	{'N8' : ['20160420_183630'],
# 					 	'N16' : ['20160420_200748'],
# 					 	'N32' : ['20160420_222506'],
# 						'N64' : ['20160421_004107']}

# timestamps['plus'] = 	{'N8'  : ['20160324_181353'],
# 						 'N16' : ['20160324_223312'],
# 						 'N32' : ['20160325_030316'],
# 						 'N64' : ['20160325_160138']}

# ssro_calib_folder = 'd:\\measuring\\data\\20160413\\162401_AdwinSSRO_SSROCalibration_111no2_SIL2'


####################################
#### 111no2 SIL 2 timestamps
####################################


timestamps['min'] = 	{'N8' : ['20170809_132903'],#['20160420_183630'],
					 	'N16' : ['20170809_145256'],#['20160420_200748'],
					 	'N32' : ['20170809_162715'],#['20160420_222506'],
						'N64' : ['20170809_183126']}#['20160421_004107']}

timestamps['plus'] = 	{'N8'  : ['20170804_083153'],#'N8'  : ['20160309_173550'],
						 'N16' : ['20170804_090407'],#'N16' : ['20160309_193753'],
						 'N32' : ['20170804_090407'],# 32 has not been taken in the latest data set'N32' : ['20160309_220938'],
						 'N64' : ['20170804_094039']}#'N64' : ['20160310_012834']}

ssro_calib_folder = 'd:\\measuring\\data\\20170704\\074219_SSRO_calib_MWInit_111no2_SIL2_SSROCalibration_MWInit'


data_folder = None




def load_data(N = [8], el_trans = 'min'):
	a = {}
	folder = {}
	for i in range(len(N)):
		# load data removing append for memory error!
		print 'loading data'
		for ii,tstamp in enumerate(timestamps[el_trans]['N'+str(N[i])]):
			a_temp, folder_temp = fp_funcs.load_mult_dat(tstamp, 
					ssro_calib_folder=ssro_calib_folder, data_folder=data_folder)

			if ii == 0:
				sweep_pts = a_temp.sweep_pts
				pts = a_temp.pts
				p0 = a_temp.p0
				u_p0 = a_temp.u_p0

			else:
				sweep_pts = np.concatenate((sweep_pts,a_temp.sweep_pts))
				pts+=a_temp.pts
				p0 = np.concatenate((p0,a_temp.p0))
				u_p0 = np.concatenate((u_p0,a_temp.u_p0))

		a_temp.pts = pts
		a_temp.p0 = p0
		a_temp.u_p0 = u_p0
		a_temp.sweep_pts = sweep_pts
		a['N'+str(N[i])] = a_temp
		folder['N'+str(N[i])] = folder_temp

		print 'data N' +str(N[i]) + 'for el_trans ' + el_trans + ' loaded'
	print 'All data for the specified N loaded via the timestamps in fp_ls'
	return a, folder

def fingerprint(a = None, folder = None, disp_sim_spin = True, N = [8], 
	el_trans = 'min', HF_perp = None, HF_par = None,xlim=None,xticks=None,**kw):

	# allowed params:
	# el_trans = ['min', 'plus']
	# N = 8, 16, 32, 64
	if (a == None) | (folder == None):
		print 'Folder path or data (a) missing fool!'
		return

	### Load hyperfine params	
	if disp_sim_spin == True:	
		if (HF_perp == None) & (HF_par == None):
			HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = el_trans, NV = 'Pippin_SIL1')
		elif el_trans == 'min':
			# print 'I did a thing!@'
			# needs to be flipped for simulation
			HF_par =  [x * (-1) for x in HF_par]

		# security check could be removed
		if len(HF_perp) == len(HF_par):
			pass
		else:
			print 'Unequal amount of Parallel and Perpendicular HF parameters'

		# print 'HF_perp = ' + str(HF_perp)
		# print 'HF_par = ' + str(HF_par)

	else:
		print 'No HF simulation'

	
	# print 'N = ' + str(N)

	N_keys = ['N'+str(pulse_no) for pulse_no in N]
	# print N_keys
	for pulses,data,datafolder in zip(N,[a[x] for x in N_keys],[folder[x] for x in N_keys]):

		##########################
		### 	 plot data ######
		#########################
		
		if xlim == None:
			fig = data.default_fig(figsize=(35,5))
			ax = data.default_ax(fig)
			ax.set_xlim(3.5,23.5)
			xlim = [3.5,23.5]
		else:
			fig = data.default_fig(figsize=(5+3*(xlim[1]-xlim[0]),2))
			ax = data.default_ax(fig)
			ax.set_xlim(xlim)

		start, end = ax.get_xlim()
		if xticks == None:
			ax.xaxis.set_ticks(np.arange(start, end, 0.5))
		else:
			ax.xaxis.set_ticks(np.arange(start, end, xticks))
		ax.set_ylim(-0.05,1.05)
		ax.plot(data.sweep_pts, data.p0, '.-k', lw=0.6,label = 'data',zorder=10)
		# print 'these are the sweep_pts',data.sweep_pts
		plt.ylabel(r'F($|0 \rangle$)',fontsize = 12)
		plt.xlabel(r'$\tau$ ($\mu s$)',fontsize = 12)
		#######################
		# Add simulated spins #
		#######################
		if disp_sim_spin == True:
			# print 'Starting Simulation for N = ' + str(pulses) + ' on transition ' + str(el_trans) 
			B_Field = kw.pop('B_field',414.2) # in gauss
			# print B_Field
			tau_lst = np.linspace(xlim[0]*1e-6, xlim[1]*1e-6, 2000)
			Mt16 = SC.dyn_dec_signal(HFs_par = HF_par, HFs_orth = HF_perp,
				B_field = B_Field, N = pulses, tau = tau_lst)
			FP_signal16 = ((Mt16+1)/2)
			

			# plot simulated results
			# colors = ['m', 'b', 'r', 'g', 'c']
			colors = cm.rainbow(np.linspace(0, 1, len(HF_par)))

			
			if el_trans == 'min':
				# flip sign back after simulation for correct graph legend
				HF_par_legend =  [x * (-1) for x in HF_par]

			else:
				HF_par_legend = HF_par
			for tt in range(len(HF_par)):
			  ax.plot(tau_lst*1e6, FP_signal16[tt,:] ,'-',lw=1.5,label = str(tt + 1) + ': HF_par = ' +str(HF_par_legend[tt]) + '; HF_perp = ' +str(HF_perp[tt]), color = colors[tt])
			plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

		print datafolder
		plt.savefig(os.path.join(datafolder, str(disp_sim_spin)+'fingerprint.pdf'),bbox_inches = 'tight',
			format='pdf')
		plt.savefig(os.path.join(datafolder, str(disp_sim_spin)+'fingerprint.png'),bbox_inches = 'tight',
			format='png')

def fingerprint_v2(a = None, folder = None, disp_sim_spin = True, N = [8], 
	el_trans = 'min', HF_perp = None, HF_par = None,xlim=None,xticks=None):

	# allowed params:
	# el_trans = ['min', 'plus']
	# N = 8, 16, 32, 64
	if (a == None) | (folder == None):
		print 'Folder path or data (a) missing fool!'
		return

	### Load hyperfine params	
	if disp_sim_spin == True:	
		if (HF_perp == None) & (HF_par == None):
			HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = el_trans, NV = 'Pippin_SIL3')
		elif el_trans == 'min':
			# needs to be flipped for simulation
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

	N_keys = ['N'+str(pulse_no) for pulse_no in N]

	for pulses,data,datafolder in zip(N,[a[x] for x in N_keys],[folder[x] for x in N_keys]):

		##########################
		### 	 plot data ######
		#########################
		
		if xlim == None:
			fig = data.default_fig(figsize=(35,5))
			ax = data.default_ax(fig)
			ax.set_xlim(3.5,23.5)
			xlim = [3.5,23.5]
		else:
			# 5+30*(xlim[1]-xlim[0])
			# 5+5*(xlim[1]-xlim[0])
			fig = data.default_fig(figsize=(5+2*(xlim[1]-xlim[0]),5))
			ax = data.default_ax(fig)
			ax.set_xlim(xlim)

		start, end = ax.get_xlim()
		if xticks == None:
			ax.xaxis.set_ticks(np.arange(start, end, 0.5))
		else:
			ax.xaxis.set_ticks(np.arange(start, end, xticks))
		ax.set_ylim(-0.05,1.05)
		ax.plot(data.sweep_pts, data.p0, '.-k', lw=0.4)#,label = 'Data')

		#######################
		# Add simulated spins #
		#######################
		if disp_sim_spin == True:
			# print 'Starting Simulation for N = ' + str(pulses) + ' on transition ' + str(el_trans) 
			B_Field = kw.pop('B_field',414.2)
			# print B_Field
			tau_lst = np.linspace(xlim[0]*1e-6, xlim[1]*1e-6, 2000)
			Mt16 = SC.dyn_dec_signal(HFs_par = HF_par, HFs_orth = HF_perp,
				B_field = B_Field, N = pulses, tau = tau_lst)
			FP_signal16 = ((Mt16+1)/2)
			

			# plot simulated results
			# colors = ['m', 'b', 'r', 'g', 'c']
			colors = cm.rainbow(np.linspace(0, 1, len(HF_par)))

			
			if el_trans == 'min':
				# flip sign back after simulation for correct graph legend
				HF_par =  [x * (-1) for x in HF_par]
			for tt in range(len(HF_par)):
				# ax.text(tau_c,h_c,str(tt+1),color = colors[tt])# better in illustrator
				ax.plot(tau_lst*1e6, FP_signal16[tt,:] ,'-',lw=1,label = 'C' + str(tt + 1), color = colors[tt])
			plt.legend(loc=3, borderaxespad=0.,frameon = False)

		plt.title('Fingerprint for N = ' +str(pulses) + ' pulses')
		plt.ylabel(ur'$\langle X_e \rangle$',fontsize = 20)
		plt.xlabel(ur'$\tau (\mu s)$',fontsize = 20)


		print datafolder
		plt.savefig(os.path.join(datafolder, str(disp_sim_spin)+'fingerprint.pdf'),
		    format='pdf')
		plt.savefig(os.path.join(datafolder, str(disp_sim_spin)+'fingerprint.png'),
		    format='png')
