import os
import numpy as np
from matplotlib import pyplot as plt
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import ssro, sequence
from analysis.lib.fitting import fit, ramsey
reload(ramsey)

def analyse_Ramsey(folder='', T2 = 3e3, Ampl = -1./3, detuning = 3e-3,hf_N = 2.17e-3, *arg):

	timestamp = None

	guess_tau = T2
	guess_a = 0.5
	guess_A = Ampl

	guess_hf_N = hf_N
	guess_det = detuning
	guess_hf_C = hf_C
	

	if timestamp != None:
	    folder = toolbox.data_from_time(timestamp)
	elif folder !='':
	    folder = toolbox.latest_data('Ramsey')

	a = sequence.SequenceAnalysis(folder)
	a.get_sweep_pts()
	a.get_readout_results(name='ssro')
	a.get_electron_ROC()
	x= a.sweep_pts
	y=a.p0


	ax = a.plot_result_vs_sweepparam(ret='ax', name='ssro')

	params_0, fitfunc_0, fitfunc_str_0 = ramsey.fit_ramsey_14N_fixed_13C_opt(guess_tau, guess_A, guess_a, guess_det, guess_hf_N)
	x_0=np.linspace(0,a.sweep_pts[-1],1000)
	ax.plot(x_0,fitfunc_0(x_0), 'r--', lw=1)
	#fit_xvals=np.linspace(res['x'][0],res['x'][-1],fit_num_points)


	fit_result = fit.fit1d(x, y, ramsey.fit_ramsey_14N_fixed_13C_opt,
	        guess_tau, guess_A, guess_a, guess_det, guess_hf_N,  fixed=[],
	        do_print=True, ret=True)
	#fit_result = False
	if fit_result != False :
		plot.plot_fit1d(fit_result, np.linspace(0,a.sweep_pts[-1],201), ax=ax,
	        plot_data=False)



	plt.savefig(os.path.join(folder, 'electronramsey_analysis.pdf'),
	        format='pdf')
