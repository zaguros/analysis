'''
Script to analyze the dynamical decoupling data
'''

import os, sys
import numpy as np
import h5py
import logging

from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr


def analyze_dark_esr(x,y,guess_ctr, guess_splitN,
guess_offset = 1,
guess_width = 0.2e-3,
guess_amplitude = 0.3,
ret='f0',
**kw):

    j=0
    min_dip_depth = 0.9
    print 'j = '+str(j)

    print y[21]
    k = len(y)
    while y[j]>min_dip_depth and j < len(y)-2:  #y[j]>0.93*y[j+1]: # such that we account for noise
        k = j
        j += 1
    #j = len(y)-2
    if k > len(y)-3:
        print 'Could not find dip'
        return
    else:
        guess_ctr = x[k]+ guess_splitN #convert to GHz and go to middle dip
        print 'guess_ctr= '+str(guess_ctr)
        print 'k'+str(k)
    ## I added this to be more robust for SSRO calibration.Please monitor if this is better - Machiel may-2014

    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            (3, guess_splitN),
            do_print=True, ret=True, fixed=[4])

    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']
        return f0, u_f0

def analyze_dark_esr_single(x,y,guess_ctr, 
guess_offset = 1,
guess_width = 0.2e-3,
guess_amplitude = 0.3,
ret='f0',
**kw):

    guess_ctr = x[y.argmin()]
    print 'guess_ctr = '+str(guess_ctr)
    ## I added this to be more robust for SSRO calibration.Please monitor if this is better - Machiel may-2014
    guess_offset=np.average(y)
    dip_threshold=guess_offset-1.5*np.std(y)
    print guess_offset
    print dip_threshold
    print min(y)
    if min(y) > dip_threshold:
        print 'Could not find dip'
        return

    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            do_print=True, ret=True, fixed=[])


    f0 = fit_result['params_dict']['x0']
    u_f0 = fit_result['error_dict']['x0']
    return f0, u_f0


def plot_fake_esr(centre=2.025,pts_course = 81, range_course = 5e-3,pts_fine = 51, range_fine =0.25):

    a= 0.963;    A = 0.182;sigma = 0.14e-3;x0 = centre+0.0003;    s0=2.179e-3

    xc = np.linspace(centre -range_course,centre+range_course,pts_course)

    yc = a-abs(A)*(0.6*np.exp(-((xc-x0-s0)/sigma)**2)+0.94*np.exp(-((xc-x0)/sigma)**2)+1.3*np.exp(-((xc-x0+s0)/sigma)**2)) 
    noise_c = np.random.normal(0,0.01,len(xc))

    yc = yc+noise_c
    yc_e = 0.01*np.ones(pts_course)

    xf = np.linspace(centre-range_fine,centre+range_fine,pts_fine)

    yf = a-abs(A)*(np.exp(-((xf-x0-s0)/sigma)**2)+np.exp(-((xf-x0)/sigma)**2)+np.exp(-((xf-x0+s0)/sigma)**2)) 
    noise_f = np.random.normal(0,0.005,len(xf))
    yf = yf+noise_f
    yf_e = 0.005*np.ones(pts_fine)

    f0c,uf0c = analyze_dark_esr(xc,yc,None,2.18e-3)
    a= analyze_dark_esr_single(xf,yf,None,2.18e-3)

    fig = figure(1,figsize=(10,5))
    ax = fig.add_subplot(111)
    ax.errorbar(xc,yc,yc_e)

    fig = figure(2,figsize=(10,5))
    ax = fig.add_subplot(111)
    ax.errorbar(xf,yf,yf_e)
    # ax.set_xlabel('msmt #')
    # ax.set_ylabel('relative course center freq (kHz) (offset 2.87748 GHz)')
    # plt.savefig('freq_vs_time_course',format='png')

    print f0c, uf0c
    print 'fine'
    print a




plot_fake_esr(centre = 2.025,pts_course = 51, range_course = 5e-3,pts_fine = 51, range_fine =0.25e-3)