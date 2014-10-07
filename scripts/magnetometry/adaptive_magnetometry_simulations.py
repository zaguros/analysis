
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


maj_reps = 5
M = 5
set_magnetic_field = 4e6 
s = magnetometry.RamseySequence_Simulation (N_msmnts = 7, reps=100, tau0=20e-9)
#s.B_max = 600e6
std = []
thresholds = []

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
	
beta, p, err = s.mean_square_error(set_value=set_magnetic_field, do_plot=True)





