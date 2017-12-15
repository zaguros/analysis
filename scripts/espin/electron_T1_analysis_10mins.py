import numpy as np
import os, sys
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import plot
from analysis.lib.tools import toolbox

from matplotlib import pyplot as plt
from analysis.lib.fitting import fit, common
from matplotlib import pyplot as plt

import matplotlib as matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

def get_T1_data_uncorrected(folder):
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    #a.get_electron_ROC()
    #a.get_cr_results('ssro')
    #a.plot_cr_vs_sweep()
    x = a.sweep_pts
    #y = a.p0
    y=a.ssro_results
    reps=a.reps
    print reps
    print y
    #y_err = a.u_p0
    #y_err=a.u_normalized_ssro
    #ax = a.plot_result_vs_sweepparam(ret='ax')
    return x,y,reps


def electron_T1_mul_3_uncorrected(older_than='20161111_091500',newer_than='20161110_224400', Amplitude=1, offset=1, T1=1e9, do_print = False,contains='T1'):

	Folder_list = toolbox.latest_data(contains=contains,older_than=older_than,newer_than=newer_than,return_all=True)
	x_tot=np.zeros(3)
	y_tot=np.zeros(3)
	reps_tot=np.zeros(3)
	y_var_tot=np.zeros(3)

	for i in range(len(Folder_list)):
	    print Folder_list[len(Folder_list)-i-1]
	    Folder = Folder_list[len(Folder_list)-i-1]
	    x,y,reps = get_T1_data_uncorrected(Folder)
	    y_tot+=y
	    reps_tot+=reps

	a = sequence.SequenceAnalysis(Folder)
	a.get_sweep_pts()
	a.get_readout_results('ssro')

	print reps_tot
	a.ssro_results=y_tot
	a.reps=reps_tot
	print a.reps


	a.normalized_ssro = a.ssro_results/((a.reps)/len(a.sweep_pts))
	a.u_normalized_ssro = \
	(a.normalized_ssro*(1.-a.normalized_ssro)/((a.reps)/len(a.sweep_pts)))**0.5  #this is quite ugly, maybe replace?

	a.get_electron_ROC()
	y_tot=a.p0
	y_var_tot=a.u_p0

	return x_tot,y_tot,y_var_tot

#################################################################
# x1,y1,y_var1 = electron_T1_mul_3(older_than='20161119_205000',newer_than='20161118_110000',contains='10_min')

# x_tot.append(x1)
# y_tot.append(y1)
# y_var_tot.append(y_var1)


# x_tot=[]
# y_tot=[]
# y_var_tot=[]