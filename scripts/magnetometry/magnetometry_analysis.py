import os, sys
import numpy as np
import h5py
import logging
import csv

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot

reload(sequence)

### settings
timestamp = None#'20140603_134433' #' #'114103_PulsarD' #YYYYmmddHHMMSS

def adaptive_analysis (folder):

    #if ax == None:
    #    fig, ax = plt.subplots(1,1)
    ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration')
    print ssro_calib_folder
    a = sequence.SequenceAnalysis(folder)
  
    RO_data, phase = a.get_magnetometry_results(name='adwindata')
    print RO_data
    print phase
    phase_values = np.unique(phase)
    ssro_results = np.zeros(len(phase_values))
    ind = 0
    for j in phase_values:
        multiplicity = len(np.where(phase==j)[0])
        ssro_results [ind] = np.sum(RO_data[np.where(phase==j)])/(multiplicity+0.)
        ind = ind+1
        print 'phase value ', j, ' is present ', multiplicity, ' times'


    plt.plot(phase_values, ssro_results)
    plt.show()    
    return ssro_results

def fpga_calibration_analysis (folder):

    ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration')
    print ssro_calib_folder
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    t=a.get_magnetometry_data(name='adwindata')
    print 'a clicks', np.sum(a.clicks, axis=0)
    a.get_electron_ROC()
    y = a.p0
    uy = a.u_p0

    a.get_sweep_pts()
    x = a.sweep_pts
    print 'x ',x
    print 'y', y
    print 'a',a
    #print y
    #print x
    #print a.reps
    #print len(a.sweep_pts)
    #print sum(y)
    nn = a.sweep_name
    ax = a.plot_result_vs_sweepparam(ret='ax')

    do_fit = True
    do_fit2 = False

    if do_fit:
        guess_frq = 1./3500
        guess_amp = 0.3
        guess_of = 1
        # guess_slope = 0.
        guess_phi = 0.
        guess_k = 0.


        o = fit.Parameter(guess_of, 'o')
        f = fit.Parameter(guess_frq, 'f')
        A = fit.Parameter(guess_amp, 'A')
        phi = fit.Parameter(guess_phi, 'phi')
        k = fit.Parameter(guess_k, 'k')
        p0 = [f, A, phi, o, k]
        fitfunc_str = ''


        fitfunc_str = 'o - A + A*e^(-kx)*cos(2pi (fx-phi))'

        def fitfunc(x):
            return (o()-A()) + A() * np.exp(-k()*x) * np.cos(2*np.pi*(f()*x - phi()))

        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=[],
                do_print=True, ret=True)
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
                plot_data=False)

        # ax.set_title(a.timestamp+'\n'+a.measurementstring)
        plt.savefig(os.path.join(folder, 'fpga_pulse_analysis_fit.png'))

    if do_fit2:
        guess_frq = 1./80
        guess_amp = 0.3
        guess_of = 0.5
        guess_phi = 0.
        guess_frq2 = 1./10
        guess_amp2 = 0.05
        guess_of2 = 1
        guess_phi2 = 0.




        o = fit.Parameter(guess_of, 'o')
        f = fit.Parameter(guess_frq, 'f')
        A = fit.Parameter(guess_amp, 'A')
        phi = fit.Parameter(guess_phi, 'phi')
        f2 = fit.Parameter(guess_frq, 'f')
        A2 = fit.Parameter(guess_amp, 'A')
        phi2 = fit.Parameter(guess_phi, 'phi')
        p0 = [f, A, phi, o, f2, A2, phi2]

        fitfunc_str = 'o - A + A*e^(-kx)*cos(2pi (fx-phi))'

        def fitfunc(x):
            return o() + A() * np.exp(np.cos(2*np.pi*(f()*x - phi()))) + A2() * np.exp(np.cos(2*np.pi*(f2()*x - phi2()))) 

        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=[],
                do_print=True, ret=True)
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
                plot_data=False)

        # ax.set_title(a.timestamp+'\n'+a.measurementstring)
        plt.savefig(os.path.join(folder, 'fpga_pulse_analysis_fit.png'))


    '''

    plt.errorbar(x, y, yerr=uy, fmt = 'o')
    plt.ylim([0,1.05])
    plt.xlabel(nn)
    plt.show()    
    return x,y, uy
    '''

def load_scope_trace(fName=''):
    with open(fName, 'rb') as f:
        reader = csv.reader(f)
        ind = 0
        v = []
        for row in reader:
            if (ind>20):
                v.append(float(row[1]))
            ind = ind + 1 
        v = np.squeeze(np.array(v))
        return v

'''

def p_to_prob (p=[], points = 0):
    points = (len (p)-1)/2
    g = np.zeros(np.size(p))+1j*np.zeros(np.size(p))
    g[0:points] = p[points+1:]
    g[points:] = p[0:points+1]
    xx = np.abs(np.fft.fftshift(np.fft.fft(g)))
    norm = np.sum(xx)
    time_domain_prob = xx/(norm+0.)
    return time_domain_prob


def adaptive_1_msmnt (N = 5, beta0 = 0., theta = 0., debug = False, adaptive = False, c = 1., corr_time=20, sigma=0.3):

    points = 2**(N+2)
    discr_steps = 2*points+1
    
    x = np.linspace (0, 2*np.pi, discr_steps)
    m = np.zeros (N+1)
    t = np.zeros (N+1)
    th = np.zeros(N+1)
    p = np.zeros ((N+1, discr_steps))+1j*np.zeros ((N+1, discr_steps))
    p [0, points] = 1/(2.*np.pi)
    t_total = 2**(N+1)+1
    
    
    for n in np.arange(N)+1:
    
        t[n] = 2**(N-n)
        ttt = -2**(N-n+1)
        th[n] = #load from file
        m[n] = #load from file
        
        #update rule:
        for k in np.arange(-t_total-1, t_total+1):
            p [n, k+points] = 0.5*p[n-1, k+points] + 0.25*(np.exp(1j*(m[n]*np.pi+th[n]))*p [n-1, k-t[n]+points] + 
                    np.exp(-1j*(m[n]*np.pi+th[n]))*p [n-1, k+t[n]+points])
        p [n, :] = (p[n, :]/(0.+p[n,points]))*(1./6.28)
    
    
    return t, th, m, p, np.array(field)
'''



if __name__ == '__main__':
    #timestamp = '124812'
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
        if folder==None:
            folder = toolbox.latest_data(timestamp) 
    else:
        folder = toolbox.latest_data('adptv_estimation')
    print folder

    #v = load_scope_trace(fName='M:/tnw/ns/qt/Diamond/Projects/Magnetometry with adaptive measurements/oscilloscope_traces/003.csv')    
    #plt.plot (v)
    #plt.show()

    fpga_calibration_analysis (folder)
    #adaptive_analysis(folder)







