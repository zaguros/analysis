import os
import numpy as np
from matplotlib import pyplot as plt
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import ssro, sequence
from analysis.lib.fitting import fit, ramsey, esr
reload(ramsey)

from analysis.lib.tools import plot
from analysis.scripts.bell import calibration_tools
reload(calibration_tools)


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

    return fit_result


def analyse_dark_esr(**kw):

    timestamp     = kw.pop('timestamp', None)
    ax            = kw.pop('ax', None)
    ret           = kw.pop('ret', None)
    min_dip_depth = kw.pop('min_dip_depth', 0.85)
    do_print      = kw.pop('do_print', False)
    plot_initial_guess = kw.pop('plot_initial_guess', False)

    guess_amplitude = kw.pop('guess_amplitude', 0.3)
    guess_width     = kw.pop('guess_width', 0.2e-3)
    guess_offset    = kw.pop('guess_offset', 1)
    guess_x0        = kw.pop('guess_x0', 2.805)
    guess_Nsplit    = kw.pop('guess_Nsplit', 2.196e-3)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
        if folder==None:
            folder = toolbox.latest_data(timestamp) 
    else:
        folder = toolbox.latest_data('DarkESR')
    print 'Folder for the Dark ESR analysis : ', folder


    if ax == None:
        fig, ax = plt.subplots(1,1)
    ssro_calib_folder = toolbox.latest_data(contains='130113_AdwinSSRO_SSROCalibration')
    print ssro_calib_folder

    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

    x = a.sweep_pts # convert to MHz
    y = a.p0.reshape(-1)[:]
    a.plot_result_vs_sweepparam(ret=ret, name='ssro', ax=ax)
    #ax.set_ylim(0.1,1.05)

    
    j=0
    while y[j]>min_dip_depth and j < len(y)-2:  #y[j]>0.93*y[j+1]: # such that we account for noise
        k = j
        j += 1
    #j = len(y)-2
    if k > len(y)-5:
        print 'Could not find dip'
        return
    else:
        guess_ctr = x[k]+ guess_Nsplit #convert to GHz and go to middle dip
        print 'guess_ctr= '+str(guess_ctr)


    
    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            # (2, guess_splitN),
            # (2, guess_splitC),
            # (2, guess_splitB),
            (3, guess_Nsplit),
            do_print=do_print, ret=True, fixed=[])
        
    plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), ax=ax, plot_data=False, **kw)


    plot_initial_guess = False
    if plot_initial_guess :
        params_0, fitfunc_0, fitfunc_str = esr.fit_ESR_gauss (guess_offset,
            guess_amplitude, guess_width, guess_ctr,(2, guess_splitC), (3, guess_splitN) )

    ax.set_xlabel('MW frq (GHz)')
    ax.set_ylabel(r'fidelity wrt. $|0\rangle$')
    ax.set_title(a.timestamp+'\n'+a.measurementstring)

    plt.savefig(os.path.join(folder, 'darkesr_analysis.png'),
            format='png')

    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']

        ax.text(f0, 0.8, '$f_0$ = ({:.3f} +/- {:.3f})'.format(
            (f0-2.8)*1e3, u_f0*1e3), ha='center')

        return (f0-2.8)*1e3, u_f0*1e3
    
    return fit_result


def analyse_pi_pulse(**kw):
    '''
    Perform a parabolic fit of the data to extract the optimal voltage for the pi pulse.
    '''
    timestamp     = kw.pop('timestamp', None)
    do_print      = kw.pop('do_print', False)

    x0_guess      = kw.pop('x0_guess', 0.6)
    a_guess       = kw.pop('a_guess', 12)
    of_guess      = kw.pop('of_guess', 0.9)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
        if folder==None:
            folder = toolbox.latest_data(timestamp) 
    else:
        folder = toolbox.latest_data('Pi')


    fig, ax = plt.subplots(1,1, figsize=(4.5,4))
    fit_result=calibration_tools.fit_parabolic(folder, x0_guess=x0_guess,a_guess=a_guess,of_guess=of_guess, ax=ax, do_print = do_print)#, info_xy=(x0_guess,0.5))
    calibration_tools.plot_result(folder, ax=ax, ret=True)
    #ax.set_ylim(0.,0.4)

    plt.savefig(os.path.join(folder, 'pi_pulse_calibration_analysis_fit.png'))

    return fit_result


def analyse_pi2_pulse(**kw):
    
    timestamp     = kw.pop('timestamp', None)
    do_print      = kw.pop('do_print', False)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
        if folder==None:
            folder = toolbox.latest_data(timestamp) 
    else:
        folder = toolbox.latest_data('Pi2')
    
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()    
    x = a.sweep_pts
    y = a.p0
    u_y = a.u_p0
    n = a.sweep_name
    a.finish()
    
    x2 = x[::2]
    y2 = y[1::2] - y[::2]
    u_y2 = np.sqrt(  u_y[1::2]**2 + u_y[::2]**2 )    
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,4), sharex=True)
    ax1.errorbar(x2, y2, yerr=u_y2, fmt='o')
    ax1.set_xlabel(n)
    ax1.set_title('Difference btw. Pi/2-Pi and Pi/2')
    ax1.set_ylabel('Difference')
    ax2.set_title(a.timestamp)
    m = fit.Parameter((y[-1]-y[0])/(x[-1]-x[0]), 'm')
    x0 = fit.Parameter(x2.mean(), 'x0')
    p0 = [m, x0]
    
    def ff(x):
        return m() * (x-x0())
    fitfunc_str = 'm * (x - x0)'
    
    fit_result = fit.fit1d(x2, y2, None, p0=p0, fitfunc=ff,
        fitfunc_str=fitfunc_str, do_print=do_print, ret=True)    
        
    ax2.errorbar(x2, y[0::2], yerr=u_y[0::2], fmt='o',
                 label='Pi/2 - Pi')
    ax2.errorbar(x2, y[1::2], yerr=u_y[1::2], fmt='o',
                 label='Pi/2')
    ax2.legend(frameon=True, framealpha=0.5)
    ax2.set_ylabel('P(0)')
    ax2.set_xlabel(n)
    

    if fit_result != False  :
        plot.plot_fit1d(fit_result, np.linspace(x2[0],x2[-1],201), ax=ax1,
            plot_data=False, print_info=do_print)
        if  a.sweep_pts[0] < x0() < a.sweep_pts[-1] :
            ax2.axvline(x0(), c='k', lw=2)
            ax2.axhline(0.5, c='k', lw=2)
            ax2.set_title('X marks the spot')
    
    fig.savefig(os.path.join(folder, 'pi2_calibration.png'))

    return fit_result


def analyse_pi4_pulse(**kw):
    
    timestamp     = kw.pop('timestamp', None)
    do_print      = kw.pop('do_print', False)
    pi4_calib     = kw.pop('pi4_calib', '1')

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
        if folder==None:
            folder = toolbox.latest_data(timestamp) 
    else:
        folder = toolbox.latest_data('Pi4')
    
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()    
    x = a.sweep_pts
    y = a.p0
    u_y = a.u_p0
    n = a.sweep_name
    a.finish()
    
    x2 = x[::2]
    y2 = y[1::2] - y[::2]
    u_y2 = np.sqrt(  u_y[1::2]**2 + u_y[::2]**2 )    
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,4), sharex=True)
    ax1.errorbar(x2, y2, yerr=u_y2, fmt='o')
    ax1.set_xlabel(n)
    if pi4_calib == '1' :
        ax1.set_title('Difference btw. Pi/4-Pi-Pi/4 and Pi/2')
    elif pi4_calib == '2' :
        ax1.set_title('Difference btw. Pi/2-Pi-Pi/4 and Pi/4')
    ax1.set_ylabel('Difference')
    ax2.set_title(a.timestamp)
    m = fit.Parameter((y[-1]-y[0])/(x[-1]-x[0]), 'm')
    x0 = fit.Parameter(x2.mean(), 'x0')
    p0 = [m, x0]
    
    def ff(x):
        return m() * (x-x0())
    fitfunc_str = 'm * (x - x0)'
    
    fit_result = fit.fit1d(x2, y2, None, p0=p0, fitfunc=ff,
        fitfunc_str=fitfunc_str, do_print=do_print, ret=True)    
    
    if pi4_calib == '1':    
        ax2.errorbar(x2, y[0::2], yerr=u_y[0::2], fmt='o',
                     label='Pi/4 - Pi - Pi/4')
        ax2.errorbar(x2, y[1::2], yerr=u_y[1::2], fmt='o',
                     label='Pi/2')
    elif pi4_calib == '2':
        ax2.errorbar(x2, y[0::2], yerr=u_y[0::2], fmt='o',
                     label='Pi/2 - Pi - Pi/4')
        ax2.errorbar(x2, y[1::2], yerr=u_y[1::2], fmt='o',
                     label='Pi/4')

    ax2.legend(frameon=True, framealpha=0.5)
    ax2.set_ylabel('P(0)')
    ax2.set_xlabel(n)
    

    if fit_result != False  :
        plot.plot_fit1d(fit_result, np.linspace(x2[0],x2[-1],201), ax=ax1,
            plot_data=False, print_info=do_print)
        if  a.sweep_pts[0] < x0() < a.sweep_pts[-1] :
            ax2.axvline(x0(), c='k', lw=2)
            if pi4_calib == '1':
            	ax2.axhline(0.5, c='k', lw=2)
            elif pi4_calib == '2':
	            ax2.axhline(np.sqrt(np.sin(np.pi/4.)), c='k', lw=2)
            ax2.set_title('X marks the spot')
    
    if pi4_calib == '1':
        fig.savefig(os.path.join(folder, 'pi4_calibration_1.png'))
    elif pi4_calib == '2':
        fig.savefig(os.path.join(folder, 'pi4_calibration_2.png'))

    return fit_result





# BUG TO FIX !!!
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
