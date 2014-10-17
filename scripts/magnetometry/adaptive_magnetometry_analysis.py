
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

reload(magnetometry)

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
			beta_exp, p_exp, err_exp, mB, sB=analyze_single_instance(label=label0, compare_to_simulations=False)
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
	s.convert_to_dict()
	s.print_results()
	s.CR_after_postselection()
	B_dict, index_dict = s.B_vs_index()

	#beta, prob, err, mB, sB = s.mean_square_error(do_plot=True, save_plot=True)
	if compare_to_simulations:
		beta_exp, p_exp, err_exp, mB, sB=s.compare_to_simulations(show_plot = False, do_save=True,plot_log=True)
	else:
		beta_exp, p_exp, err_exp, mB, sB=s.mean_square_error(show_plot = False, save_plot=True, do_plot=False)

	#s.analyse_ramsey()
	return beta_exp, p_exp, err_exp, mB, sB

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
def analyze_sweep_field():

	mgnt_exp = magnetometry.AdaptiveMagnetometry(N=7, tau0=20e-9)
	mgnt_exp.set_protocol (M=3, maj_reps = 7, maj_thr = 2)
	mgnt_exp.set_sweep_params (nr_periods = 2, nr_points_per_period=7)
	mgnt_exp.load_sweep_field_data (N=2)
	mgnt_exp.load_sweep_field_data (N=3)
	mgnt_exp.load_sweep_field_data (N=4)
	mgnt_exp.load_sweep_field_data (N=5)
	mgnt_exp.load_sweep_field_data (N=6)
	mgnt_exp.load_sweep_field_data (N=7)
	plt.figure()

	mgnt_exp.plot_msqe_dictionary(y_log=True)

	mgnt_exp.plot_sensitivity_scaling(do_fit = True)
	mgnt_exp.save()
	#mgnt_exp.load_analysis(timestamp='20141014')

analyze_sweep_field()
#analyze_single_instance(label='CR40gr_manyreps_2', compare_to_simulations=True)
#error= B_vs_time(nr=90, label = 'CR40b_manyreps')
#single_B_field (nr = 69, label = 'CR40b_manyreps')