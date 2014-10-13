
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



def analysis_exp():
	f = toolbox.latest_data(contains='220832')
	s = RamseySequence_Exp (folder = f)
	s.set_exp_pars (T2=4000e-9, fid0=0.85, fid1=0.015)
	#s.M=1
	#s.set_ideal()
	s.load_data()
	p = s.plot_avg_phase_distribution()
	#s.print_phases()
	p = s.phase_distribution_scaling (do_plot=True)
	p.plot_scaling()

	


'''
for n in np.arange(10)+1:
	for m in np.arange(7)+1:
		
		print '##### N = '+str(n)+' --- M = '+str(m)
		
		t = AdaptiveTable (N=6,M=3)
		t.verbose = False
		t.generate()
		t.save_table()
'''


'''
bbb = -0*4*93.285e6
maj_reps = 3
s = RamseySequence_Simulation (N_msmnts = 6, reps=100, tau0=1e-9)
s.B_max = 600e6
std = []
thresholds = []
#for thr in [0, 1,2, 3]:
s.setup_simulation (magnetic_field_hz = bbb, M=3)
s.T2 = 4000e-9
s.fid0 = 1#0.85
s.fid1 = 0#0.02
s.renorm_ssro = False
s.maj_reps = maj_reps
s.maj_thr = 0
s.table_based_simulation()
s.convert_to_dict()
	
beta, p, err = s.mean_square_error(set_value=bbb)
plt.plot (beta*1e-6, p, linewidth =1.5)
#media = np.sum (beta*p)
#v = np.sum(beta*beta*p)-media**2
#	std.append((v**0.5)*1e-6)
plt.xlabel ('magnetic field [MHz]')
plt.legend()
plt.show()
'''


'''
plt.plot (thresholds, std)
plt.plot (thresholds, std, 'ob')
	thresholds.append(thr)
plt.title ('majority reps: '+str(maj_reps))
plt.xlabel ('threshold')
plt.ylabel ('std [MHz]')
plt.show()
'''


f = toolbox.latest_data(contains='adptv_estimation_det')
#f = toolbox.latest_data(contains='152727')
s = magnetometry.RamseySequence_Exp (folder = f)
s.set_exp_pars (T2=96e-6, fid0=0.85, fid1=0.015)
print f
s.load_exp_data()
s.convert_to_dict()
s.print_results()
beta, prob, err = s.mean_square_error(do_plot=True, save_plot=True)
s.analyse_ramsey()






