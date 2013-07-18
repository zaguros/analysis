import os, sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pprint

from analysis.lib.fitting import fit
from analysis.lib.tools import plot
from measurement.lib.tools import toolbox

from analysis.lib.m2.ssro import mbi
from analysis.lib.m2.ssro import ssro

import dark_esr_analysis

reload(ssro)
reload(dark_esr_analysis)

# adapt
name = 'sil2'
pi2_4mhz_value = 1. - 0.473



def stage_1_calibrations():
    # ssro first
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,4))

    print 80*'='
    print 'Dark ESR'
    print 80*'='
    f0,u_f0 = dark_esr_analysis.analyze_dark_esr(
        toolbox.latest_data('DarkESR'), ax=ax1, ret='f0',
        print_info=False)
    
    print 80*'='
    print 'Slow pi pulse for MBI'
    print 80*'='
    slow_pi_length = slow_pi(ax2)


def stage_2_calibrations():
    fig, ax = plt.subplots(1,1)

    print 80*'='
    print '4 MHz electron Rabi'
    print 80*'='
    CORPSE_freq = rabi_4mhz(ax)

def stage_3_calibrations():
    fig, ([ax1, ax2],[ax3, ax4], [ax5, ax6]) = plt.subplots(3,2, figsize=(10,15))
    
    print 80*'='
    print '4 MHz Pi'
    print 80*'='
    pi_4mhz_amp = pi_4mhz(ax1)
    
    #print 80*'='
    #print '4 MHz Pi/2'
    #print 80*'='
    #pi2_4mhz_amp = pi2_4mhz(ax2)

    print 80*'='
    print 'CORPSE pi'
    print 80*'='
    CORPSE_pi_amp = CORPSE_pi(ax3)

    print 80*'='
    print 'pi2pi'
    print 80*'='
    pi2pi_pi_amp = pi_pi2pi(ax4)

    print 80*'='
    print 'pi2pi mI=0'
    print 80*'='
    pi2pi_pi_mI0_amp = pi_pi2pi_mI0(ax5)

    ##print 80*'='
    #print 'Hard pi'
    #print 80*'='
    #hard_pi_amp = pi_hard(ax6)



### generic calibration functions, don't use those directly
def calibrate_epulse_amplitude(folder, ax, *args):
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    a.plot_results_vs_sweepparam(ax=ax, name='adwindata')
    
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    
    x0 = fit.Parameter(args[0], 'x0')
    of = fit.Parameter(args[1], 'of')
    a = fit.Parameter(args[2], 'a')
    fitfunc_str = '(1-of) + a (x-x0)**2'
    
    def fitfunc_parabolic(x):
        return (1.-of()) + a() * (x-x0())**2
    
    fit_result = fit.fit1d(x,y, None, p0=[of, a, x0], fitfunc=fitfunc_parabolic,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

def fit_linear(folder, ax, value, *args):
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    a.plot_results_vs_sweepparam(ax=ax, name='adwindata')
    
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]

    a = fit.Parameter(args[0], 'a')
    b = fit.Parameter(args[0], 'b')
    fitfunc_str = 'a x + b'

    def fitfunc(x):
        return a() * x + b()

    fit_result = fit.fit1d(x,y, None, p0=[a, b], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result


def calibrate_epulse_rabi(folder, ax, *args, **kws):
    fit_phi = kws.pop('fit_phi', True)
    fit_k = kws.pop('fit_k', True)

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    a.plot_results_vs_sweepparam(ax=ax, name='adwindata')
    
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]

    f = fit.Parameter(args[0], 'f')
    A = fit.Parameter(args[1], 'A')
    phi = fit.Parameter(0, 'phi')
    k = fit.Parameter(0, 'k')
    p0 = [f, A]
    if fit_phi:
        p0.append(phi)
    if fit_k:
        p0.append(k)
    fitfunc_str = '(1 - A) + A * exp(-kx) * cos(2pi f x - 2pi phi)'

    def fitfunc(x) : 
        return (1.-A()) + A() * np.exp(-k()*x) * \
            np.cos(2*np.pi*(f()*x - phi()))

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

def epulse_fidelity(folder, ax, *args):
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    a.plot_results_vs_sweepparam(ax=ax, name='adwindata')
    
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    
    of = fit.Parameter(args[0], 'of')
    fitfunc_str = '(1-of)'
    
    def fitfunc_fid(x):
        return (1.-of())
    
    fit_result = fit.fit1d(x,y, None, p0=[of], fitfunc=fitfunc_fid,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

### end generic calibration functions

### actual functions to call, for specific calibrations    
def slow_pi(ax=None):
    folder = toolbox.latest_data('cal_slow_pi_'+name)
    if ax==None:
        fig,ax = plt.subplots(1,1)
    fit_result = calibrate_epulse_rabi(folder, ax, 1./5000, 0.5)

    f = fit_result['params_dict']['f']
    u_f = fit_result['error_dict']['f']
    ax.text(500, 0.4, 'pi = (%.0f +/- %.0f) ns' % (0.5/f, 0.5/f**2 * u_f),
        va='bottom', ha='left')
    
    return (0.5/f, 0.5/f**2 * u_f)


def rabi_4mhz(ax=None):
    folder = toolbox.latest_data('cal_4mhz_rabi'+name)
    if ax==None:
        fig,ax = plt.subplots(1,1)
    fit_result = calibrate_epulse_rabi(folder, ax, 1./250, 0.5, fit_k=False)

    f = fit_result['params_dict']['f']
    u_f = fit_result['error_dict']['f']
    ax.text(100, 0.9, '$f_r$ = (%.3f +/- %.3f) MHz' % (f*1e3, u_f*1e3),
        va='bottom', ha='left')

    return (f*1e3, u_f*1e3)

def pi2_4mhz(ax=None):
    folder = toolbox.latest_data('cal_4MHz_pi_over_2_'+name)
    if ax==None:
        fig,ax = plt.subplots(1,1)

    #fit_result = calibrate_epulse_amplitude(folder, ax, 0.68,  1, 0)
    #A = fit_result['params_dict']['x0']
    #u_A = fit_result['error_dict']['x0']
    #ax.text(0.6, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))
    fit_result = fit_linear(folder, ax, -0.1/0.15 , 1)
    a = fit_result['params_dict']['a']
    u_a = fit_result['error_dict']['a']
    b = fit_result['params_dict']['b']
    u_b = fit_result['error_dict']['b']

    v = pi2_4mhz_value
    A = (v - b)/a
    u_A = np.sqrt( ((b-v)/a**2 * u_a)**2 + (1/a * u_b)**2 )
    ax.text(0.6, 0.2, 'pi/2 = ({:.3f} +/- {:.3f}) V'.format(A, u_A) )

    return A, u_A 
    
def pi_4mhz(ax=None):
    folder = toolbox.latest_data('cal_4MHz_pi_'+name)
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.68,  1, 0)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    ax.text(0.6, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A 

def CORPSE_pi(ax=None):
    folder = toolbox.latest_data('CORPSEPiCalibration') # _'+name)
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.68,  1, 0)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    ax.text(0.6, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A 

def pi_pi2pi(ax=None):
    folder = toolbox.latest_data('cal_pi2pi_pi_'+name)
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.165, 1, 0)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    ax.text(0.16, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A

def pi_pi2pi_mI0(ax=None):
    folder = toolbox.latest_data('cal_pi2pi_pi_mI0')
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.165, 1, 0)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    ax.text(0.16, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A

def pi_hard(ax=None):
    folder = toolbox.latest_data('cal_hard_pi_'+name)
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.85,  1, 0)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    ax.text(0.81, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A 

#### Not fully finished yet

def pi_397ns(ax=None):
    folder = toolbox.latest_data('cal_397ns_pi_'+name)
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.16,  1, 0)

def hard_pi_duration(ax=None):
    folder = toolbox.latest_data('cal_hard_pi_dur_'+name)
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 120,  1, 0)

def hard_pi_amplitude(ax=None):
    folder = toolbox.latest_data('cal_hard_pi_amp_'+name)
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.85,  1, 0)
          
    
def CORPSE(ax=None):
    folder = toolbox.latest_data('Test_CORPSE')
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.68,  1, 0)
    
def CORPSE_60(ax=None):
    folder = toolbox.latest_data('cal_CORPSE60')
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.68,  1, 0)

def CORPSE_60_duration(ax=None):
    folder = toolbox.latest_data('CORPSE_vs_60length')
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 52,  1, 0)
     
    
def CORPSE_m300(ax=None):
    folder = toolbox.latest_data('cal_CORPSEm300')
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.68,  1, 0)
    
def CORPSE_420(ax=None):
    folder = toolbox.latest_data('cal_CORPSE420')
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.68,  1, 0)
    
def CORPSE_fidelity(ax=None):
    folder = toolbox.latest_data('CORPSE')
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = epulse_fidelity(folder, ax, 1)
  
                                
if __name__ == '__main__':
    # calibrate_all()
    stage_1_calibrations()