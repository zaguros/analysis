import os
import numpy as np
from matplotlib import pyplot as plt
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import ssro, sequence
from analysis.lib.fitting import fit, ramsey
reload(ramsey)

from analysis.lib.tools import plot



def analyse_Rabi(guess_frq = 2., guess_amp = 0.2, guess_of = 0.1, **kw) :

    timestamp    = kw.pop('timestamp', None)
    guess_phi    = kw.pop('guess_phi', 0.)
    guess_k      = kw.pop('guess_k', 0.)
    mbi_analysis = kw.pop('mbi_analysis', False)
    do_print     = kw.pop('do_print', False)

    o = fit.Parameter(guess_of, 'o')
    f = fit.Parameter(guess_frq, 'f')
    A = fit.Parameter(guess_amp, 'A')
    phi = fit.Parameter(guess_phi, 'phi')
    k = fit.Parameter(guess_k, 'k')
    p0 = [f, A, phi, o, k]
    fitfunc_str = ''


    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else :
        folder = toolbox.latest_data('ElectronRabi')

    if mbi_analysis:
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('adwindata')
        a.get_electron_ROC()
        ax = a.plot_results_vs_sweepparam(ret='ax', name = 'adwindata')

    else:
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        a.get_electron_ROC()
        ax = a.plot_result_vs_sweepparam(ret='ax')

    x = a.sweep_pts
    y = a.p0

    fitfunc_str = 'o - A + A*e^(-kx)*cos(2pi (fx-phi))'

    def fitfunc(x):
    	return (o()-A()) + A() * np.exp(-k()*x) * np.cos(2*np.pi*(f()*x - phi()))

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=[2],
        	do_print=do_print, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax,
        	plot_data=False)

    print "\npi pulse at {:.3f} for .\n".format(1/f()/2.) + a.sweep_name

    # ax.set_title(a.timestamp+'\n'+a.measurementstring)
    plt.savefig(os.path.join(folder, 'electronrabi_analysis_fit.png'))






def analyse_Ramsey(folder='', T2 = 3e3, Ampl = -1./3, detuning = 3e-3,hf_N = 2.17e-3, *arg):

	timestamp = kw.pop(timestamp, None)

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
