import os, sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pprint

from analysis.lib.fitting import fit
from analysis.lib.tools import plot
from measurement.lib.tools import toolbox

from analysis.lib.m2.ssro import sequence
from analysis.lib.m2.ssro import ssro

from analysis.scripts.espin import dark_esr_analysis

reload(ssro)
reload(dark_esr_analysis)

# adapt
name = 'sil4'
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
    print '8 MHz electron Rabi'
    print 80*'='
    CORPSE_freq = rabi_8mhz(ax2)


def stage_2_calibrations(sil=4):
    fig, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2,2, figsize = (10,8))

    print 80*'='
    print 'CORPSE pi'
    print 80*'='
    CORPSE_pi_amp = CORPSE_pi(ax1)

    print 80*'='
    print 'CORPSE pi/2'
    print 80*'='
    CORPSE_pi2_amp = CORPSE_pi2(sil,ax=ax2)

    #print 80*'='
    #print 'CORPSE pi/2 alternative 2 x'
    #print 80*'='
    #CORPSE_pi2_amp = CORPSE_pi2_alt(sil,M=2,ax=ax3)
    #
    #print 80*'='
    #print 'CORPSE pi/2 alternative 4 x'
    #print 80*'='
    #CORPSE_pi2_amp = CORPSE_pi2_alt(sil,M=4,ax=ax4)


def stage_3_calibrations():

    fig, ax = plt.subplots(1,1, figsize = (5,4))

    print 80*'='
    print 'Delta t'
    print 80*'='
    delta_t = dd_delta_t(ax)


def CORPSE_fidelity_cal():
    fig, ax = plt.subplots(1,1)

    print 80*'='
    print 'CORPSE pi fidelity'
    print 80*'='
    CORPSE_pi_fid = CORPSE_fidelity(ax)

def CORPSE_pi2_fidelity_cal():
    fig, ax = plt.subplots(1,1)

    print 80*'='
    print 'CORPSE pi/2 fidelity'
    print 80*'='
    CORPSE_pi2_fid = CORPSE_pi2_fidelity(ax)

def dd_calibrate_fid(r=4):
    fig, ax = plt.subplots(1,1)

    print 80*'='
    print 'dd sequence fidelity'
    print 80*'='
    dd_sequence_fid = dd_calibrate_fidelity(r,ax)

### generic calibration functions, don't use those directly
def calibrate_epulse_amplitude(folder, ax, *args):
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='ssro')
    a.get_electron_ROC()
    a.plot_result_vs_sweepparam(ax=ax, name='ssro')
    
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
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='ssro')
    a.get_electron_ROC()
    a.plot_result_vs_sweepparam(ax=ax, name='ssro')
    
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

    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='ssro')
    a.get_electron_ROC()
    a.plot_result_vs_sweepparam(ax=ax, name='ssro')
    
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]

    f = fit.Parameter(args[0], 'f')
    A = fit.Parameter(args[1], 'A')
    x0 = fit.Parameter(0, 'x0')
    k = fit.Parameter(0, 'k')
    p0 = [f, A]
    if fit_phi:
        p0.append(x0)
    if fit_k:
        p0.append(k)
    fitfunc_str = '(1 - A) + A * exp(-kx) * cos(2pi f (x - x0))'

    def fitfunc(x) : 
        return (1.-A()) + A() * np.exp(-k()*x) * \
            np.cos(2*np.pi*(f()*(x - x0())))

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

def epulse_fidelity(folder, ax, *args):

    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='ssro')
    a.get_electron_ROC()
    a.plot_result_vs_sweepparam(ax=ax, name='ssro')
    
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    
    of = fit.Parameter(args[0], 'of')
    of2 = fit.Parameter(0, 'of2')
    of3 = fit.Parameter(0, 'of3')
    fitfunc_str = '(1-of)'
    
    def fitfunc_fid(x) :
        return (1.-of())

    fit_result = fit.fit1d(x,y, None, p0=[of], fixed = [], fitfunc=fitfunc_fid,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    #plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
    #    plot_data=False, print_info=False)
       
    return fit_result

### end generic calibration functions

### actual functions to call, for specific calibrations    

def rabi_8mhz(ax=None):
    folder = toolbox.latest_data('cal_8mhz_rabi')
    if ax==None:
        fig,ax = plt.subplots(1,1)
    fit_result = calibrate_epulse_rabi(folder, ax, 1./125, 0.5, fit_k=False)

    f = fit_result['params_dict']['f']
    u_f = fit_result['error_dict']['f']
    ax.text(50, 0.85, '$f_r$ = (%.3f +/- %.3f) MHz' % (f*1e3, u_f*1e3),
        va='bottom', ha='left')

    return (f*1e3, u_f*1e3)


def CORPSE_pi(ax=None):
    folder = toolbox.latest_data('CORPSEPiCalibration') 
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.75,  1, 0 )
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    ax.text(0.42, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A 


def CORPSE_pi2(sil,ax=None):
    folder = toolbox.latest_data('CORPSEPi2Calibration') 
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = fit_linear(folder, ax, -1,  1)
    a = fit_result['params_dict']['a']
    b = fit_result['params_dict']['b']
    u_a = fit_result['error_dict']['a']
    u_b = fit_result['error_dict']['b']
    A = (0.5 - b) / a
    u_A = np.sqrt(A**2) * np.sqrt ( (u_a/a)**2 + (u_b/b)**2 )
    ax.text(0.42, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A 
 
def CORPSE_pi2_alt(sil,M,ax=None):
    folder = toolbox.latest_data('CORPSEPi2Calibration_sil'+str(sil)+str(M)) 
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.46,  0, 0 )
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    ax.text(0.42, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A 

def CORPSE_fidelity(ax=None):
    folder = toolbox.latest_data('CORPSEPiCalibration') 
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = epulse_fidelity(folder, ax, 1)
    fid = fit_result['params_dict']['of']
    u_fid = fit_result['error_dict']['of']
    ax.text(0.42, 0.5, 'of = (%.3f +/- %.3f)' % (fid, u_fid))

    return fid, u_fid 

def CORPSE_pi2_fidelity(ax=None):
    folder = toolbox.latest_data('CORPSEPi2Calibration') 
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = epulse_fidelity(folder, ax, 0.5)
    fid = fit_result['params_dict']['of']
    u_fid = fit_result['error_dict']['of']
    ax.text(0.6, 0.2, 'of = (%.3f +/- %.3f)' % (fid, u_fid))

    return fid, u_fid 
                                

def dd_delta_t(ax=None):
    folder = toolbox.latest_data('DynamicalDecoupling')

    if ax==None:
        fig,ax = plt.subplots(1,1)
           
    fit_result = calibrate_epulse_amplitude(folder, ax, -0.24,  1, 0 )
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    ax.text(-0.24, 0.5, 'dt = (%.3f +/- %.3f) ns' % (A, u_A))

    return A, u_A 

def dd_calibrate_fidelity(r, ax=None):
    folder = toolbox.latest_data('DynamicalDecoupling')

    if ax==None:
        fig,ax = plt.subplots(1,1)
           
    fit_result = calibrate_epulse_amplitude(folder, ax, r*108,  1, 0 )
    A = fit_result['params_dict']['of']
    u_A = fit_result['error_dict']['of']
    ax.text(r*108, 0.5, 'F = (%.3f +/- %.3f)' % (A, u_A))

    return A, u_A 

if __name__ == '__main__':
    # calibrate_all()
    stage_1_calibrations()
    stage_2_calibrations()