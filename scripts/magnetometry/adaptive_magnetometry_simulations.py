
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
#from analysis.lib.m2.ssro import sequence
#from analysis.lib.tools import toolbox
##from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot
from matplotlib import rc, cm
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry

reload(magnetometry)

def simulate_cappellaro ():
	maj_reps = 5
	M = 5

	set_magnetic_field =-22.1875e6
	s = magnetometry.RamseySequence_Simulation (N_msmnts = 7, reps=200, tau0=20e-9)

	s.setup_simulation (magnetic_field_hz = set_magnetic_field, M=M)
	s.T2 = 96e-6
	s.fid0 = 0.868
	s.fid1 = 1-0.978
	s.renorm_ssro = True
	s.maj_reps = maj_reps
	s.maj_thr = 1

	s.table_based_simulation()
	#s.sim_cappellaro_majority()
	s.convert_to_dict()
	s.print_table_positions()
		
	beta, p, err,a, b = s.mean_square_error(set_value=set_magnetic_field, do_plot=True)

	s.sim_cappellaro_majority()
	s.convert_to_dict()
	s.print_results()
	beta, p, err,a,b = s.mean_square_error(set_value=set_magnetic_field, do_plot=True)



def simulate_nonadaptive ():
	set_magnetic_field = 4e6 
	s = magnetometry.RamseySequence_Simulation (N_msmnts = 7, reps=100, tau0=20e-9)

	s.setup_simulation (magnetic_field_hz = set_magnetic_field, M=M)
	s.T2 = 96e-6
	s.fid0 = 0.9
	s.fid1 = 0.02
	s.renorm_ssro = False
	s.maj_reps = maj_reps
	s.maj_thr = 1
	#s.table_based_simulation()
	s.sim_cappellaro_majority()
	s.convert_to_dict()
	s.print_results()
		
	beta, p, err,a,b = s.mean_square_error(set_value=set_magnetic_field, do_plot=True)


def simulate_sweep_field(M, maj_reps, maj_thr):
	print '############### Simulate #####################'
	mgnt_exp = magnetometry.AdaptiveMagnetometry(N=6, tau0=20e-9)
	mgnt_exp.set_protocol (M=M, maj_reps = maj_reps, maj_thr = maj_thr)
	mgnt_exp.set_sweep_params (reps =250, nr_periods = 5, nr_points_per_period=11)
	mgnt_exp.set_exp_params( T2 = 96e-6, fid0 = 0.87, fid1 = 0.02)
	mgnt_exp.sweep_field_simulation (N=2)
	mgnt_exp.sweep_field_simulation (N=3)
	mgnt_exp.sweep_field_simulation (N=4)	
	mgnt_exp.sweep_field_simulation (N=5)
	mgnt_exp.sweep_field_simulation (N=6)	
	#mgnt_exp.sweep_field_simulation (N=7)
	plt.figure()
	mgnt_exp.plot_msqe_dictionary()
	mgnt_exp.plot_scaling()
	mgnt_exp.save()


def analyze_saved_simulations (timestamp):
	mgnt_exp = magnetometry.AdaptiveMagnetometry(N=6, tau0=20e-9)
	mgnt_exp.load_analysis (timestamp=timestamp)
	mgnt_exp.plot_msqe_dictionary(y_log=True)
	mgnt_exp.plot_sensitivity_scaling()

analyze_saved_simulations (timestamp='120405')

#simulate_sweep_field (M=7, maj_reps=5, maj_thr=1)
