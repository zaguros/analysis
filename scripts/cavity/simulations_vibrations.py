"""
This script provides the length distribtuions of a cavity based on measured vibrations,
needed to calculate the effect of vibrations on the 
"""
import numpy as np
import scipy.constants
import math
from matplotlib import pyplot as plt
from analysis.lib.fitting import common
from analysis.lib.fitting import fit
from analysis.lib.tools import plot
c = scipy.constants.c

n_diamond = 2.41
############# parameters during slow scan:
optical_cavity_length_slow_scans = 53.*(c/470.8e12)/2
dnu_vibr_slow_scans = 13.3e9 #gaussian fit in Hz

#####

def gaussian(x, mu, FWHM):
    sig = FWHM/(2.*np.sqrt(2*np.log(2)))
    return 1./(math.sqrt(2*sig*math.pi))*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def avg_p_ZPL_to_zero(ZPL_to_zero,p_cav_length):
    return np.sum(ZPL_to_zero*p_cav_length)/np.sum(p_cav_length)

def dnu_vibr_to_dL_vibr(dnu_vibr=dnu_vibr_slow_scans,cavity_length=optical_cavity_length_slow_scans,res_freq=c/637.e-9):
    return dnu_vibr/(res_freq)*(cavity_length)# in m


###################################
dL_vibr_slow_scans = dnu_vibr_to_dL_vibr(dnu_vibr_slow_scans,cavity_length=optical_cavity_length_slow_scans)
Pseudo_finesse_slow_scans = (c/(2*optical_cavity_length_slow_scans))/dnu_vibr_slow_scans
####################################




###################
#Functions to generate a lifetime curve, taking into account vibrations
###################

def weights_for_lifetime_curves(lifetime,into_ZPL,p_cav_length,det_eff=0.2,fs_det_eff=0.03):
    return (det_eff*into_ZPL+fs_det_eff*(1-into_ZPL))*1./lifetime*p_cav_length/np.sum(p_cav_length)#in Hz

def lifetime_under_vibrations(lifetime,intoZPL, p_cav_length, ts,det_eff=0.2,fs_det_eff=0.03):
    weights =  weights_for_lifetime_curves(lifetime,intoZPL,p_cav_length,det_eff,fs_det_eff)
    curves = np.zeros((len(ts),len(lifetime)))
    for i,(A,tau) in enumerate(zip(weights,lifetime)):
        curves[:,i] = A*np.exp(-ts/tau)
    return np.sum(curves,axis=1)

def fit_lifetime_curve_from_dL(p_cav_length,dL,lt,b,g_tau=2.,det_eff=0.2,fs_det_eff=0.03,ts=np.linspace(0,30e-9,31)):
    lifetime_curve = lifetime_under_vibrations(lt,b,p_cav_length,ts,det_eff=det_eff,fs_det_eff=fs_det_eff) 
    lifetime_curve_norm = lifetime_curve/lifetime_curve[0]
    plt.plot(lifetime_curve)
    plt.show()
    
    g_a = 0
    g_A = 1.
    fixed = [0]
    p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(g_a, g_A, g_tau)
    fit_result = fit.fit1d(ts*1e9,lifetime_curve_norm, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
    ax=plot.plot_fit1d(fit_result,ts*1.e9, add_txt=True,label='Fit',show_guess=False, plot_data=True,ret='ax')
    return lifetime_curve_norm,fit_result

