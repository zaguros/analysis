import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib.m2 import m2
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.scripts.espin import dark_esr_analysis_p7889 as dea


def analyze_DESR_roomtemp(folder, 
    guess_ctr,
    ROI_start = 1670,
    ROI_end = 2200,
    guess_splitN = 2e6,
    guess_offset = 1,
    guess_width = 0.2e-3,
    guess_amplitude = 0.3,
    ret = 'f0'
    ):
    """
    Fit Gaussians to Dark ESR result.
    """

    # Fetch & postprocess data
    # a = dea.Analyse_DarkESR(folder)
    a = dea.DarkESRanalysis(folder)
    a.get_data('raw_data')
    a.eval_data(ROI_start, ROI_end)
    x = a.sweep_pts
    y = a.summation
    y_normalized = y / np.amax(y)

    # Fit Gaussians
    j = 0
    while y_normalized[j]>0.97 and j < len(y)-2: # find indices for 'non-dips'
        k = j
        j = j+1

    if k > len(y)-5:
        print 'Could not find dip'
        return
    else:
        guess_ctr = x[k]+ guess_splitN #convert to GHz and go to middle dip
        print 'guess_ctr = '+str(guess_ctr)

    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            (3, guess_splitN),
            do_print=True, ret=True, fixed=[4])


    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']

        return f0, u_f0


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
    print 'j = '+str(j)
    while y[j]>0.85 and j < len(y)-2: # such that we account for noise
        k = j
        j = j+1

    if k > len(y)-5:
        print 'Could not find dip'
        return
    else:
        guess_ctr = x[k]+ guess_splitN #convert to GHz and go to middle dip
        print 'guess_ctr = '+str(guess_ctr)

    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            (3, guess_splitN),
            do_print=True, ret=True, fixed=[4])


    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']

        return f0, u_f0








