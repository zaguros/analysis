import os, sys
import numpy as np
import h5py
import logging


from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr



def analyze_dark_esr(guess_ctr, guess_splitN,
guess_offset = 1,
guess_width = 0.2e-3,
guess_amplitude = 0.3,
timestamp = None,
ret='f0',
**kw):

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('DarkESR')

    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()

    x = a.sweep_pts # convert to MHz
    y = a.p0.reshape(-1)[:]

    # here we find the first of the three dips
    j=0
    while y[j]>0.85 and j < len(y)-2: # such that we account for noise
        k = j
        j = j+1
    #j = len(y)-2
    if k > len(y)-3:
        print 'Could not find dip'
        return
    else:
        guess_ctr = x[k]+ guess_splitN #convert to GHz and go to middle dip
        print 'guess_ctr= '+str(guess_ctr)

    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            (3, guess_splitN),
            do_print=True, ret=True, fixed=[4])


    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']

        return f0, u_f0








