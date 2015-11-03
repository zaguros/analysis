import pickle
import numpy as np
import os
from matplotlib import pyplot as plt
import matplotlib.cm as cm

import analysis.scripts.QEC_data_analysis.fingerprints.nr1_sil18.fitting_fingerprints.hyperfine_params as module_hyperfine_params; reload(module_hyperfine_params) 
import analysis.scripts.QEC_data_analysis.fingerprints.nr1_sil18.fitting_fingerprints.nuclear_spin_characterisation as SC; reload(SC)

hf = module_hyperfine_params.hyperfine_params

def get_hyperfine_params(ms = 'plus', carbon_spins = 'all',NV = None):
    
    HF_perp = []
    HF_par 	= []
    # hf = module_hyperfine_params.hyperfine_params
    if carbon_spins == 'all':
      for kk in range(len(hf)):
        carbon_string     = 'C' + str(kk+1) 
        HF_perp.append(hf[carbon_string]['perp'])
        if ms == 'plus':
          HF_par.append(hf[carbon_string]['par'])
        elif ms == 'min':
          HF_par.append(-1*hf[carbon_string]['par'])

    else:
      for carbon in carbon_spins:
        carbon_string     = 'C' + str(carbon) 
        HF_perp.append(hf[carbon_string]['perp'])
        if ms == 'plus':
          HF_par.append(hf[carbon_string]['par'])
        elif ms == 'min':
          HF_par.append(-1*hf[carbon_string]['par'])


    return HF_perp, HF_par

def fingerprint_plot(disp_sim_spin = True,n_sim_spins = 2,
        xrange = [0,20],figsize=(35,5), Nr_of_pulses = None, ms = 'plus'):


	###################
	# Add simulated spins #
	###################

	if disp_sim_spin == True:
	    
	    HF_perp, HF_par = get_hyperfine_params(ms = ms)

	    B_Field = 403.555 

	    tau_lst = np.linspace(0,72e-6,10000)
	    Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,Nr_of_pulses,tau_lst)
	    FP_signal16 = ((Mt16+1)/2)



	# data = pickle.load( open( 'sil18_fingerprint_ms_'+ ms +'_N'+str(Nr_of_pulses)+'.p', 'rb' ) )

	# x = data['x']
	# y = data['y']
	# x,y = np.loadtxt('sil18_fingerprint_ms_plus_N'+str(Nr_of_pulses)+'.txt')


	fig,ax = plt.subplots(figsize=figsize)
	ax.set_xlim(4.9,5.1)
	ax.set_xlim(xrange)
	start, end = ax.get_xlim()
	ax.xaxis.set_ticks(np.arange(xrange[0], xrange[1]+0.1, (xrange[1]- xrange[0])/10.))

	ax.set_ylim(-1.05,1.05)

	# ax.plot(x, y, '.-k', lw=0.4,label = 'data') #N = 16
	ax.set_xlabel('tau (us)')
	ax.set_ylabel('Contrast')

	if disp_sim_spin == True:
		colors = cm.rainbow(np.linspace(0, 1, n_sim_spins))
		for tt in range(n_sim_spins):
			ax.plot(tau_lst*1e6, Mt16[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])


	lgd = plt.legend(loc=4)
	plt.show(block = False)


# for ms in ['plus','min']:
# 	for Nr_of_pulses in [4,8,16,32,64]:

# 		fingerprint_plot(disp_sim_spin = True,n_sim_spins = 8,
# 		xrange = [8,15],figsize=(35,5), Nr_of_pulses = Nr_of_pulses, ms = ms)


fingerprint_plot(disp_sim_spin = True,n_sim_spins = 2,
		xrange = [8,15],figsize=(35,5), Nr_of_pulses = 32, ms = 'min')		