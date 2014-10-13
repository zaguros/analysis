
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

	f = toolbox.latest_data(contains='adptv_estimation_det')
	#f = toolbox.latest_data(contains='170245')
	s = magnetometry.RamseySequence_Exp (folder = f)
	s.set_exp_pars (T2=96e-6, fid0=0.876, fid1=1-.964)
	print f
	s.load_exp_data()
	s.convert_to_dict()
	s.print_results()
	beta, prob, err, mB, sB = s.mean_square_error(do_plot=True, save_plot=True)
	#s.compare_to_simulations()
	#s.analyse_ramsey()

def analyze_sweep_field():

	mgnt_exp = magnetometry.AdaptiveMagnetometry(N=6, tau0=20e-9)




analyze_single_instance()





