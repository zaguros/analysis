
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
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot
from matplotlib import rc, cm
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
#from analysis.lib.magnetometry import adwin_debug_magnetometry as adwin_mgnt

reload(magnetometry)
#reload(adwin_mgnt)

def simulate_cappellaro ():

	F = 1
	G = 100
	K = 0

	set_magnetic_field =12.5e6/2.
	s = magnetometry.RamseySequence_Simulation (N_msmnts = 1, reps=201, tau0=20e-9,)

	s.setup_simulation (magnetic_field_hz = set_magnetic_field, G=G,F=F,K=K)
	s.T2 = 96000e-6
	s.fid0 = 1.-0.025
	s.fid1 = 0.025

	#s.table_based_simulation()
	#s.sim_cappellaro_majority()
	s.sim_cappellaro_variable_M()

def simulate_sql (B=3*12.5e6/4.):

	F = 0
	G = 100
	K = 0

	set_magnetic_field = B
	print set_magnetic_field
	s = magnetometry.RamseySequence_Simulation (N_msmnts = 1, reps=101, tau0=10e-9,)

	s.setup_simulation (magnetic_field_hz = set_magnetic_field, G=G,F=F,K=K)
	s.T2 = 96000e-6
	s.fid0 = 0.87
	s.fid1 = 0.025
	s.B_max = 1/(4*s.t0)
	s.n_points = 1000

	print 'Simulating data...\n'
	s.sim_SQL(phase_deg=90)
	print 'Done!'
	print 'Plotting normal results...'
	print 'Size results matrix: ', np.shape(s.msmnt_results)
	s.convert_to_dict()
	beta_py, p_py, av_exp_py,H_py, m_py, s_py = s.mean_square_error(set_value=set_magnetic_field, do_plot=True)
	plt.show()

	for i in [200, 500, 1000]:
		print 'Reshaping msmnt result into G = ', i
		print 'Plotting...'
		s.reshape_SQL(G=i)
		print 'Size results matrix: ', np.shape(s.msmnt_results)
		print 'Size phase matrix: ', np.shape(s.msmnt_phases)
		print 'Size time array: ', np.shape(s.msmnt_times)
		s.convert_to_dict()
		beta_py, p_py, av_exp_py,H_py, m_py, s_py = s.mean_square_error(set_value=set_magnetic_field, do_plot=True)
		plt.show()


def simulate_berry (do_adaptive):

	F = 1
	G = 2
	K = 3
	set_magnetic_field = 12.5e6/2.
	s = magnetometry.RamseySequence_Simulation (N_msmnts = 4, reps=1000, tau0=20e-9)

	s.setup_simulation (magnetic_field_hz = set_magnetic_field, G=G,F=F,K=K)
	s.T2 = 96000e-6
	s.fid0 = .87
	s.fid1 = 0.02

	s.sim_berry_protocol(do_adaptive=do_adaptive)
	s.convert_to_dict()
	beta_py, p_py, av_exp_py,H_py, m_py, s_py = s.mean_square_error(set_value=set_magnetic_field, do_plot=True)
	plt.show()


def test_swarm_opt ():

	F = 2
	G = 5
	K = 8
	set_magnetic_field = 12.5e6/2.
	s = magnetometry.RamseySequence_Simulation (N_msmnts = K+1, reps=10, tau0=20e-9)

	s.setup_simulation (magnetic_field_hz = set_magnetic_field, G=G,F=F,K=K)
	s.T2 = 96e-6
	s.fid0 = .87
	s.fid1 = 0.02

	s.sim_swarm_optim()
	s.convert_to_dict()
	beta_py, p_py, av_exp_py,H_py, m_py, s_py = s.mean_square_error(set_value=set_magnetic_field, do_plot=True)
	plt.show()




def simulate_cappellaro_debug_adwin ():
	maj_reps = 1
	M = 10

	set_magnetic_field =25/(20e-9*2**8)
	s = magnetometry.RamseySequence_Adwin (N_msmnts = 8, reps=100, tau0=20e-9)

	s.setup_simulation (magnetic_field_hz = set_magnetic_field, M=M)
	s.T2 = 96e-6
	s.fid0 = 0.9
	s.fid1 = 0.02
	s.renorm_ssro = True
	s.maj_reps = maj_reps
	s.maj_thr = 0

	#s.table_based_simulation()
	#s.sim_cappellaro_majority_with_plots(N_sim = [1,2, 3,4])
	s.adwin_ultrafast_M(nr_coeff=51)
	s.convert_to_dict()
	s.print_results()
	#s.print_table_positions()
		
	beta, p, err,a, b = s.mean_square_error(set_value=set_magnetic_field, do_plot=True)
	#plt.plot (beta, p)
	#plt.show()

	#s.sim_cappellaro_majority()
	#s.convert_to_dict()
	s.print_results()
	#beta, p, err,a,b = s.mean_square_error(set_value=set_magnetic_field, do_plot=True)



def simulate_nonadaptive ():
	set_magnetic_field = 4e6 
	s = magnetometry.RamseySequence_Simulation (N_msmnts = 7, reps=100, tau0=20e-9)

	s.setup_simulation (magnetic_field_hz = set_magnetic_field,  G=G,F=F,K=K)
	s.T2 = 96e-6
	s.fid0 = 0.9
	s.fid1 = 0.02

	#s.table_based_simulation()
	s.sim_cappellaro_majority()
	s.convert_to_dict()
	s.print_results()
		
	beta, p, err,a,b = s.mean_square_error(set_value=set_magnetic_field, do_plot=True)


def simulate_sweep_field_variable_M(G,F,K,fid0, protocol, fid1=0.02,print_results=False,reps=101, error_bars = True, specific_B=False):
#def simulate_sweep_field_variable_M(G,F,K,fid0, do_adaptive, fid1=0.02,print_results=False,reps=101, phase_update=False, error_bars = True, always_recalculate_phase=False,specific_B=False):

	#try:
	print '############### Simulate #####################'
	N=K+1
	mgnt_exp = magnetometry.AdaptiveMagnetometry(N=N, tau0=20e-9)
	mgnt_exp.set_protocol (G=G,K=K,F=F)
	mgnt_exp.set_sweep_params (reps =reps, nr_periods = 21, nr_points_per_period=15)
	mgnt_exp.set_exp_params( T2 = 5000e-6, fid0 = fid0, fid1 = fid1)
	mgnt_exp.error_bars = error_bars
	for n in np.arange(N-1)+2:
		mgnt_exp.set_protocol (G=G,K=n-1,F=F)
		mgnt_exp.verbose=True
		#mgnt_exp.sweep_field_simulation (N=n,do_adaptive=do_adaptive,print_results=print_results, phase_update=phase_update, always_recalculate_phase=always_recalculate_phase,specific_B=specific_B)
		print 'T2* = ', mgnt_exp.T2
		mgnt_exp.sweep_field_simulation (N=n, protocol = protocol ,print_results=print_results, specific_B=specific_B)

		plt.figure()
		
		mgnt_exp.plot_msqe_dictionary(y_log=True)
	mgnt_exp.plot_sensitivity_scaling()
	mgnt_exp.save(name='_noT2')


def simulate_sweep_field_SQL (fid0, fid1=0.02,print_results=False,reps=501, error_bars = True, specific_B=False):
#def simulate_sweep_field_variable_M(G,F,K,fid0, do_adaptive, fid1=0.02,print_results=False,reps=101, phase_update=False, error_bars = True, always_recalculate_phase=False,specific_B=False):

	#try:
	print '############### Simulate #####################'
	N=1
	mgnt_exp = magnetometry.magnetometrySQL(tau0=0.5*20e-9, max_G = 3000)
	mgnt_exp.set_sweep_params (reps =reps, nr_periods = 1, nr_points_per_period=3)
	mgnt_exp.set_exp_params( T2 = 96e-6, fid0 = fid0, fid1 = fid1)
	mgnt_exp.error_bars = error_bars
	mgnt_exp.G = 100
	mgnt_exp.verbose=True
	mgnt_exp.sweep_field_simulation_SQL (print_results=print_results, specific_B=specific_B)
	plt.figure()
		
	mgnt_exp.plot_msqe_dictionary(y_log=True)
	mgnt_exp.plot_sensitivity_scaling()
	mgnt_exp.save()



def analyze_saved_simulations (timestamp,error_bars=False):
	mgnt_exp = magnetometry.AdaptiveMagnetometry(N=14, tau0=20e-9)
	mgnt_exp.error_bars=error_bars
	mgnt_exp.load_analysis (timestamp=timestamp)
	#mgnt_exp.plot_msqe_dictionary(y_log=True, save_plot=False)
	#mgnt_exp.G=G
	#mgnt_exp.F=F
	#mgnt_exp.K=K
	mgnt_exp.error_bars=error_bars
	mgnt_exp.plot_sensitivity_scaling(save_plot=False)
	return mgnt_exp

def simulate_adwin ():
	N = 7
	pp = int(np.random.randint(1, 2**(N-1)))
	set_magnetic_field = pp*(1/20e-9)/(2**N) 
	print 'B-field:', set_magnetic_field/1e6
	s = magnetometry.RamseySequence_Adwin (N_msmnts = N, reps=20, tau0=20e-9)

	s.setup_simulation (magnetic_field_hz = set_magnetic_field, M=5)
	s.T2 = 96e-6
	s.fid0 = 1.00
	s.fid1 = 0.0
	s.renorm_ssro = False
	s.maj_reps = 1
	s.maj_thr = 0
	
	s.adwin_ultrafast()
	s.convert_to_dict()
	s.print_results()
	beta, p, err,a, b = s.mean_square_error(set_value=set_magnetic_field, do_plot=False)

	plt.figure()
	#plt.plot (beta_adwin, p_adwin, 'b')
	plt.plot (beta*1e-6, p, 'r')
	plt.show()


def check_simulated_adwin_phases ():
	N = 7
	pp = int(np.random.randint(1, 2**(N-1)))
	set_magnetic_field = pp*(1/20e-9)/(2**N) 
	print 'B-field:', set_magnetic_field/1e6

	s = magnetometry.RamseySequence_Adwin (N_msmnts = N, reps=20, tau0=20e-9)
	s.setup_simulation (magnetic_field_hz = set_magnetic_field, M=5)
	s.T2 = 96e-6
	s.fid0 = 1.00
	s.fid1 = 0.0
	s.renorm_ssro = False
	s.maj_reps = 1
	s.maj_thr = 0
	
	s.check_adwin_phases(nr_coeff = 101)



def check_adwin_code (N,M, outcomes):
	s = magnetometry.RamseySequence_Adwin (N_msmnts = N, reps=1, tau0=20e-9)
	s.setup_simulation (magnetic_field_hz = 0, M=M)
	s.T2 = 96e-6
	s.fid0 = 1.00
	s.fid1 = 0.0
	s.renorm_ssro = False
	s.maj_reps = 1
	s.maj_thr = 0
	s.adwin_algorithm (outcomes = outcomes)


def benchmark_exec_speed ():

	nr_N=10
	tempo_brute = np.zeros(nr_N)
	tempo_positive  = np.zeros(nr_N)
	tempo_nonzero = np.zeros(nr_N)
	tempo_ultrafast = np.zeros(nr_N)

	ind = 0
	for n in np.arange (nr_N)+2:
		print '##### ',n
		set_magnetic_field = 2*(1/20e-9)/(2**n) 
		s = magnetometry.RamseySequence_Adwin (N_msmnts = n, reps=10, tau0=20e-9)

		s.setup_simulation (magnetic_field_hz = set_magnetic_field, M=1)
		s.T2 = 96e-6
		s.fid0 = 1.00
		s.fid1 = 0.00
		s.renorm_ssro = False
		s.maj_reps = 1
		s.maj_thr = 0
		ex1 = s.basic_adwin_algorithm(exec_speed = True)
		tempo_brute[ind] = ex1
		ex2 = s.adwin_only_positive(exec_speed = True)
		tempo_positive[ind] = ex2
		ex3 = s.adwin_positive_nonzero(exec_speed = True)
		tempo_nonzero[ind] = ex3
		ex4 = s.adwin_ultrafast(exec_speed = True)
		tempo_ultrafast[ind] = ex4
		ind +=1

	plt.semilogy (np.arange (nr_N)+2, tempo_brute, linewidth=2, label = 'brute force')
	plt.semilogy (np.arange (nr_N)+2, tempo_positive, linewidth=2, label = 'only k>0')
	plt.semilogy (np.arange (nr_N)+2, tempo_nonzero, linewidth=2, label = 'k>0, p[k] not 0')
	plt.semilogy (np.arange (nr_N)+2, tempo_ultrafast, linewidth=2, label = 'ultrafast')
	plt.xlabel ('N')
	plt.ylabel ('execution time [s]')
	plt.legend()
	plt.show()

def adwin_phase_angle (real_part, imag_part):

		c_nr = real_part+1j*imag_part

		if (real_part>0):
			if (imag_part>0):
				th = (180/np.pi)*np.arctan (imag_part/real_part)
			else:
				th = (180/np.pi)*np.arctan (imag_part/real_part)+360
		else:
			if (imag_part>0):
				th = 180-(180/np.pi)*np.arctan (-imag_part/real_part)
			else:
				th = 180+(180/np.pi)*np.arctan (imag_part/real_part)

		print 'Number :', c_nr
		print 'Using python angle: ', np.mod((180/np.pi)*np.angle(c_nr), 360)
		print 'Using arctan: ', th


def test_adwin_sims(N, M, outcomes = [], do_plot = False, do_print=False):

	a = adwin_mgnt.RamseySequence_Adwin (N_msmnts = N, reps=1, tau0=20e-9)
	a.renorm_ssro = False
	a.verbose = False
	a.maj_reps = 1
	a.maj_thr = 0	
	a.M = M
	phase_adwin, phase_python, diff, p_2tn_adwin, p_2tn_python = a.compare_algorithms(outcomes=outcomes, do_plot = do_plot, do_print=do_print)

	print '-----Phases:'
	print '** adwin: ', np.round(phase_adwin*180/np.pi)
	print '** python: ', np.round(phase_python*180/np.pi)

	diff_real = np.abs(np.real(p_2tn_adwin)-np.real(p_2tn_python))
	diff_imag = np.abs(np.imag(p_2tn_adwin)-np.imag(p_2tn_python))
	avg_phase_error = np.sum(np.abs(phase_adwin-phase_python)*180/np.pi)/float(a.N)

	if do_plot:
		f, axarr = plt.subplots(2, sharex=True, figsize=(10,10))
		axarr[0].plot (np.abs(phase_adwin-phase_python)*180/np.pi, ':g')
		axarr[0].plot (np.abs(phase_adwin-phase_python)*180/np.pi, 'og')
		axarr[1].plot (diff_real, 'ob', label = 'real_part')
		axarr[1].plot (diff_imag, 'or', label = 'imag_part')
		
		axarr[0].set_title('difference in set phases')
		axarr[1].set_title('dipperences in p[2*tn]')
		axarr[0].legend()
		axarr[1].legend()
		plt.show()

	print 'Test adwin. Msmsnt outcome: ', outcomes,' --- avg_phase_err = ', avg_phase_error, ' deg'
	return avg_phase_error


def simulate_adwin (N,M):
	a = adwin_mgnt.RamseySequence_Adwin (N_msmnts = N, reps=2000, tau0=20e-9)
	field = 100/(a.t0*2**N)
	a.renorm_ssro = False
	a.verbose = False
	a.maj_reps = 1
	a.maj_thr = 0	
	a.setup_simulation (magnetic_field_hz = field, M=M)

	a.T2 = 96e-6
	a.fid0 = 0.88
	a.fid1 = 0.02
	a.adwin_optimal(use_fid_bayesian_update=True)

	a.convert_to_dict()
	a.print_results()
		
	beta, p, err,a, b = a.mean_square_error(set_value=field, do_plot=True)
	plt.show()






#simulate_sweep_field (N=1, M=3, maj_reps=5, maj_thr=1, fid0=0.95)
'''
simulate_sweep_field (N=9, M=4, maj_reps=7, maj_thr=2, fid0=0.87)
simulate_sweep_field (N=10, M=3, maj_reps=6, maj_thr=2, fid0=0.87)
simulate_sweep_field (N=9, M=4, maj_reps=6, maj_thr=2, fid0=0.87)
simulate_sweep_field (N=10, M=3, maj_reps=5, maj_thr=1, fid0=0.87)
simulate_sweep_field (N=9, M=4, maj_reps=5, maj_thr=1, fid0=0.87)

#check_adwin_code(N=4, M=5, outcomes = [5,3,0,4])

#simulate_cappellaro_debug_adwin()
#check_simulated_adwin_phases ()



mean_error = []
m_list = np.arange(30)+1
for m in m_list:
	rep = 30
	err = []
	for i in np.arange(rep):
		a = np.array(np.random.randint(m+1, size=6))
		err.append(1e-18+test_adwin_sims(N=6, M=m, outcomes=a, do_plot=False))
	mean_error.append(np.mean(np.array(err)))
plt.plot (m_list, np.array(mean_error), ':k')
plt.plot (m_list, np.array(mean_error), 'o')
plt.ylabel ('avg phase error [deg]')
plt.xlabel('M')
plt.legend()
plt.show()
'''

#simulate_adwin(N=8, M= 20)
#test_adwin_sims(N=7, M=5, outcomes=[3,0,4,4,0,4,4], do_plot=False, do_print = True)
#simulate_cappellaro()

'''
def overnight_simulations_15dec2014():
	fid0=0.87
	fid1=0.02
	reps=11
	for g in [3,4,5,10, 20]:
		for f in [0,1,2,3,4,5,7, 10]:
			simulate_sweep_field_variable_M (G=g,K=9,F=f , fid0=fid0,fid1=fid1,print_results=False,reps=reps, phase_update=False, error_bars = True, do_adaptive=True, always_recalculate_phase= True, N1_sweep=True) 
			simulate_sweep_field_variable_M (G=g,K=9,F=f , fid0=fid0,fid1=fid1,print_results=False,reps=reps, phase_update=False, error_bars = True, do_adaptive=True, always_recalculate_phase= True, N1_sweep=False)
			simulate_sweep_field_variable_M (G=g,K=9,F=f , fid0=fid0,fid1=fid1,print_results=False,reps=reps, phase_update=True, error_bars = True, do_adaptive=False, always_recalculate_phase= False, N1_sweep=False)
'''


def suppl_info_simulations_adptv (G, fid0, sweep_f=[0,1,2,3,4,5]):
	fid1=0.02
	reps=31
	for fff in sweep_f:
		simulate_sweep_field_variable_M (protocol = 'modified_cappellaro',G=G,K=9,F=fff, fid0=fid0,fid1=fid1,reps=reps, error_bars = True)
		simulate_sweep_field_variable_M (protocol = 'non_adaptive',G=G,K=9,F=fff, fid0=fid0,fid1=fid1,reps=reps, error_bars = True)

def suppl_info_simulations_nn_adptv (G, fid0, sweep_f=[0,1,2,3,4,5]):
	fid1=0.02
	reps=31
	for fff in sweep_f:
		simulate_sweep_field_variable_M (protocol = 'non_adaptive',G=G,K=9,F=fff, fid0=fid0,fid1=fid1,reps=reps, error_bars = True)


#mgnt_MNp1_WRONG_lessreps=analyze_saved_simulations('20141105_112326',G=2,F=1,K=7)

#simulate_berry(do_adaptive=False)

#test_swarm_opt()
#simulate_sweep_field_variable_M(G=5,F=2,K=9,fid0=0.87, protocol='swarm_optimization', fid1=0.02,print_results=False,reps=21, error_bars = True, specific_B=False)

#simulate_sweep_field_SQL (fid0=0.87, fid1=0.02,print_results=False,reps=501, error_bars = True, specific_B=False)
#simulate_sql()

suppl_info_simulations_adptv (G=5, fid0=1.0, sweep_f = [0,1])
suppl_info_simulations_nn_adptv (G=5, fid0=1.0, sweep_f = [1,2,3])
suppl_info_simulations_adptv (G=5, fid0=0.95, sweep_f = [0,1,2])
suppl_info_simulations_nn_adptv (G=5, fid0=0.95, sweep_f = [1,2,3])
suppl_info_simulations_adptv (G=5, fid0=0.9, sweep_f = [1,2,3])
suppl_info_simulations_nn_adptv (G=5, fid0=0.9, sweep_f = [1,2,3])
suppl_info_simulations_adptv (G=5, fid0=0.85, sweep_f = [1,2,3])
suppl_info_simulations_nn_adptv (G=5, fid0=0.85, sweep_f = [3,4,5])
suppl_info_simulations_adptv (G=5, fid0=0.8, sweep_f = [1,2,3])
suppl_info_simulations_nn_adptv (G=5, fid0=0.8, sweep_f = [4,5,6])
suppl_info_simulations_adptv (G=5, fid0=0.75, sweep_f = [2,3,4])
suppl_info_simulations_nn_adptv (G=5, fid0=0.75, sweep_f = [5,6,7])
suppl_info_simulations_adptv (G=5, fid0=0.7, sweep_f = [3,4,5])
suppl_info_simulations_nn_adptv (G=5, fid0=0.7, sweep_f = [6,7,8])


#suppl_info_simulations (G=5, fid0=0.75)
