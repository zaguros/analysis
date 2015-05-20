
import numpy as np
from analysis.lib.fitting import fit, ramsey, common
from analysis.lib.tools import plot
import random
from matplotlib import rc, cm
import os, sys
import h5py
import logging, time

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import  sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot
from matplotlib import rc, cm
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
from analysis.lib.magnetometry import adwin_debug_magnetometry as adwin_mgnt

reload(magnetometry)
reload(adwin_mgnt)

def B_vs_time_single (label):

	f = toolbox.latest_data(contains=label)
	s = magnetometry.RamseySequence_Exp (folder = f)
	s.set_exp_pars (T2=96e-6, fid0=0.87, fid1=1-.975)
	print f
	s.load_exp_data()
	s.convert_to_dict()
	B_dict, index_dict = s.B_vs_index()
	print 'Array before cleaning...', B_dict, index_dict
	s.CR_after_postselection()
	s.convert_to_dict()
	B_dict, index_dict = s.B_vs_index()
	print 'Array after cleaning...', B_dict, index_dict

	return B_dict, index_dict, s.CR_after, s.discarded_elements

def B_vs_time (nr, label):

	index  = []
	B_field = []
	for i in np.arange(nr)+10:
		try:
			label0=label+'_%d'%i
			B_dict, index_dict, CR, disc = B_vs_time_single (label=label0)
			last = len(index)
			for k in B_dict:
				for idx in index_dict[k]:
					index.append (last+idx)
					B_field.append (B_dict[k])
		except:
			print 'Excluded...', i
	B_field = np.array(B_field)
	index=np.array(index)				
	#plt.plot (index, B_field, 'r')
	plt.plot (index, B_field, 'ob')
	plt.xlabel('repetition')
	plt.show()

	error=[]
	for i in np.arange(nr)+10:
		try:
			label0=label+'_%d'%i
			beta_exp, p_exp, ave_exp,err_exp, mB, sB=analyze_single_instance(label=label0, compare_to_simulations=False)
			error.append(sB)
			print i, sB
		except:
			print 'Excluded...', i
	plt.figure()
	plt.plot (np.array(error))
	plt.plot (np.array(error),'ob')
	plt.ylabel ('std [MHz]')
	plt.xlabel ('instance nr')
	return error

def single_B_field(nr, label):
	index = []
	B_field = []
	label0=label+'_%d'%nr
	B_dict, index_dict, CR, discarded = B_vs_time_single (label=label0)
	print B_dict, index_dict
	for k in B_dict:
		for idx in index_dict[k]:
			index.append (idx)
			B_field.append (B_dict[k])
	for i in discarded:
		plt.axvline (x=i, ymin=0, ymax = 2, color='k')
	B_field = np.array(B_field)
	index=np.array(index)				
	#plt.plot (index, B_field, 'r')
	plt.plot (index, B_field, 'ob')
	plt.xlabel('repetition')
	plt.show()




def analyze_single_instance(label='adptv_estimation_det', compare_to_simulations=True):

	
	f = toolbox.latest_data(contains=label)
	s = magnetometry.RamseySequence_Exp (folder = f)
	s.set_exp_pars (T2=96e-6, fid0=0.87, fid1=1-.975)
	print f
	s.load_exp_data()
	s.CR_after_postselection()
	s.convert_to_dict()
	#print 'phases dict'
	#print s.phases_dict
	s.print_results()
	
	B_dict, index_dict = s.B_vs_index()

	#beta, prob, err, mB, sB = s.mean_square_error(do_plot=True, save_plot=True)
	if compare_to_simulations:
		beta_sim, p_sim, ave_exp,err_sim, a, b=s.compare_to_simulations(show_plot = True, verbose=True,do_save=True,plot_log=True)
	else:
		beta_sim, p_sim, ave_exp,err_sim, a, b=s.mean_square_error(show_plot = True, save_plot=True, do_plot=True)

	#s.analyse_ramsey()
	return beta_sim, p_sim, ave_exp,err_sim, a, b

'''
def temporal_evolution_B(label, nr):

	all_msmnts = magnetometry.RamseySequence (folder = f)
	all_msmnts.set_exp_pars (T2=96e-6, fid0=0.87, fid1=1-.975)
	s.load_exp_data()


	for i in np.arange(nr):
		f = toolbox.latest_data(contains=label)
		s = magnetometry.RamseySequence_Exp (folder = f)
		s.set_exp_pars (T2=96e-6, fid0=0.87, fid1=1-.975)
		s.load_exp_data()
'''

def single_repetition_evolution(label='adptv_estimation_det'):
	f = toolbox.latest_data(contains=label)
	s = magnetometry.RamseySequence_Exp (folder = f)
	s.set_exp_pars (T2=96e-6, fid0=0.87, fid1=1-.975)
	#print f
	s.load_exp_data()
	s.convert_to_dict()
	s.CR_after_postselection()
	#beta, prob = s.analysis_dict (phase = curr_phase, msmnt_results = curr_msmnt, times = np.rint(self.msmnt_times))
	return s
def analyze_sweep_field(G=1,F=2,nr_periods=1,phase_update=False,older_than=None,newer_than=None,CR_after_threshold=2):

	mgnt_exp = magnetometry.AdaptiveMagnetometry(N=14, tau0=20e-9)
	mgnt_exp.set_protocol (G=G,K=13,F=F)
	mgnt_exp.error_bars=True
	mgnt_exp.set_sweep_params (nr_periods = nr_periods, nr_points_per_period=7)
	mgnt_exp.set_exp_params (T2=96e-6, fid0=0.87, fid1=1-.975)
	mgnt_exp.phase_update=phase_update
	#mgnt_exp.load_sweep_field_data (N=1,older_than=older_than,newer_than=newer_than)
	mgnt_exp.load_sweep_field_data (N=2,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	mgnt_exp.load_sweep_field_data (N=3,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	mgnt_exp.load_sweep_field_data (N=4,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	mgnt_exp.load_sweep_field_data (N=5,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	mgnt_exp.load_sweep_field_data (N=6,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	'''
	mgnt_exp.load_sweep_field_data (N=7,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	mgnt_exp.load_sweep_field_data (N=8,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	mgnt_exp.load_sweep_field_data (N=9,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	mgnt_exp.load_sweep_field_data (N=10,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	mgnt_exp.load_sweep_field_data (N=11,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	mgnt_exp.load_sweep_field_data (N=12,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	mgnt_exp.load_sweep_field_data (N=13,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	mgnt_exp.load_sweep_field_data (N=14,older_than=older_than,newer_than=newer_than,CR_after_threshold=CR_after_threshold)
	'''
	#print 'mgnt_exp repetitions', mgnt_exp.repetitions
	plt.figure()

	mgnt_exp.plot_msqe_dictionary(y_log=True)

	mgnt_exp.plot_sensitivity_scaling(do_fit = True)
	mgnt_exp.save()
	#mgnt_exp.load_analysis(timestamp='20141014')

#analyze_sweep_field()
#analyze_single_instance(label='det=-12.5MHz_N=1', compare_to_simulations=True)
#error= B_vs_time(nr=90, label = 'CR40b_manyreps')
#single_B_field (nr = 69, label = 'CR40b_manyreps')

def check_adwin_realtime(label, newer_than=False, print_details=True):

	dir0, daydir, m_dirs = toolbox.latest_data (contains = label, return_all=True, newer_than = newer_than)

	for i in m_dirs:
		f = os.path.join(dir0, daydir, i)
		print '################## - current folder: ', i
		exp = magnetometry.RamseySequence_Exp (folder = f)
		exp.load_exp_data()
		#print 'Timer: '+str('{0:.2f}'.format(exp.timer))+' ms'
		exp.check_realtime_phases()
		p_tn_exp = exp.p_tn
		p_2tn_exp = exp.p_2tn

		s = adwin_mgnt.RamseySequence_Adwin (N_msmnts = exp.N, reps=1, tau0=20e-9)
		s.M = exp.M
		s.renorm_ssro = False
		s.verbose = False
		s.maj_reps = 1
		s.maj_thr = 0	
		phases_th, p_tn_th, p_2tn_th = s.sim_cappellaro_track_coefficients (msmnt_result = exp.msmnt_results[0,:])
		#print np.round(phases_th*180/np.pi)


		if print_details:
			for n in np.arange (exp.N)+1:
				err_real = abs((np.real(p_2tn_th[n-1]) - np.real(p_2tn_exp[n-1])))
				err_imag = abs((np.imag(p_2tn_th[n-1]) - np.imag(p_2tn_exp[n-1])))
				abs_value = (np.real(p_2tn_th[n-1])**2+np.imag(p_2tn_th[n-1])**2)**0.5
				error = 100*(err_real + err_imag)/abs_value

				print ' ---- N = ', n, '  ---  t_n = ', 2**(exp.N-n), '  ---- opt_phase = '+str('{0:.1f}'.format(phases_th[n-1]*180/np.pi))
				#print '		p[tn] = ',p_tn_th[n-1], ' (th)  ----- ', p_tn_exp[n-1],' (exp)'
				print '		p[2tn] = ',p_2tn_th[n-1], ' (th)  ----- ', p_2tn_exp[n-1],' (exp)'
				print '     err: real part = ', error, '%'
		s.adwin_optimal (ext_outcomes = exp.msmnt_results[0,:])

def check_adwin_realtime_record_pk(label, newer_than=False):


	f, axarr = plt.subplots(2, sharex=True, figsize=(10,10))
	for n in [1]:
		for m in [1,2,3,4]:
			dir0, daydir, m_dirs = toolbox.latest_data (contains = label+'_test_pk_(n='+str(n)+'_m='+str(m)+')', return_all=True, newer_than = newer_than)

			for i in m_dirs:
				f = os.path.join(dir0, daydir, i)
				exp = magnetometry.RamseySequence_Exp (folder = f)
				exp.load_exp_data()
				#exp.check_realtime_phases()
				x = np.arange(2**exp.N+1)
				axarr[0].plot (exp.real_pk_adwin, ':k')
				axarr[1].plot (exp.imag_pk_adwin, ':k')
				axarr[0].plot (exp.real_pk_adwin, 'o', label=str(n)+'_'+str(m))
				axarr[1].plot (exp.imag_pk_adwin, 'o', label=str(n)+'_'+str(m))

	axarr[0].set_title('real part')
	axarr[1].set_title('imaginary part')
	plt.xlim ([0, 17])
	plt.legend()
	plt.show()


def check_adwin_realtime_plots (N, M, outcomes = [], do_plot=True, do_print = False, do_plot_adwin = True, newer_than = None):

	s = adwin_mgnt.RamseySequence_Adwin (N_msmnts = N, reps=1, tau0=20e-9)
	s.M = M
	s.renorm_ssro = False
	s.verbose = False
	s.maj_reps = 1
	s.maj_thr = 0	
	s.fid0=1
	s.fid1=0
	s.T2= 96e-6
	print 'Outcomes: ', outcomes
	s.outcomes = outcomes
	s.compare_adwin_python (do_plot=do_plot, newer_than = newer_than, use_fid_bayesian_update = True)


#result = '4021'
#check_adwin_realtime (label = result+'_test_pk_(n=4_m=1)', newer_than = '102000')
#check_adwin_realtime_record_pk(label = result, newer_than = '102000')
#analyze_single_instance(label='105854',compare_to_simulations=True)
#l=['N = 2','N = 3','N = 4','N = 5','N = 6','N = 7']


#for n,label in enumerate(l):
#	print label
#analyze_single_instance(label='154340',compare_to_simulations=False)
#check_adwin_realtime (label = 'rtAdwin', newer_than = '124600', print_details=False)
#check_adwin_realtime_plots (N=4, M=5, outcomes = [5,0,2,1,5,5,0,2], newer_than='145500', do_plot=True)


