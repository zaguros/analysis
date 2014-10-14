
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

reload(magnetometry)




def analyze_single_instance():

	#f = toolbox.latest_data(contains='adptv_estimation_det')
	f = toolbox.latest_data(contains='153422')
	s = magnetometry.RamseySequence_Exp (folder = f)
	s.set_exp_pars (T2=96e-6, fid0=0.87, fid1=1-.975)
	print f
	s.load_exp_data()
	s.convert_to_dict()
	s.print_results()
	#beta, prob, err, mB, sB = s.mean_square_error(do_plot=True, save_plot=True)
	s.compare_to_simulations(show_plot = False, do_save=True)
	#s.analyse_ramsey()

def analyze_sweep_field():


	mgnt_exp = magnetometry.AdaptiveMagnetometry(N=3, tau0=20e-9)
	mgnt_exp.set_protocol (M=7, maj_reps = 5, maj_thr = 1)
	mgnt_exp.set_sweep_params (nr_periods = 2, nr_points_per_period=11)
	mgnt_exp.load_sweep_field_data (N=2)
	mgnt_exp.load_sweep_field_data (N=3)
	mgnt_exp.load_sweep_field_data (N=5)
	#mgnt_exp.load_sweep_field_data (N=6)
	mgnt_exp.plot_msqe_dictionary()
	mgnt_exp.plot_scaling()
	mgnt_exp.save()
	#mgnt_exp.load_analysis(timestamp='20141014')

#analyze_single_instance()
analyze_sweep_field()









