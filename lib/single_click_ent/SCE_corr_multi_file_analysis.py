"""
provides functions to analyze spin-spin correlators between two NV centres
Functions should be executed on the computer that stores the PQ data, i.e. LT4. (othewise use pq_folder = 'xxxx' when calling instances of purify_pq)
Based on the analysis class purify_pq and some functionalities from purify_analysis.py
"""

import os
import numpy as np
from analysis.lib.lde import sscorr ### two qubit SSRO correction
from analysis.lib.purification import purify_pq as ppq; reload(ppq)
import Analysis_params_SCE as analysis_params; reload(analysis_params)
from analysis.lib.pq import pq_tools,pq_plots
from analysis.lib.tools import plot; reload(plot)
import analysis.lib.purification.purify_analysis as purify_analysis
from analysis.lib.tools import toolbox as tb
from analysis.lib.m2.ssro import ssro
from matplotlib import pyplot as plt
import h5py
from analysis.lib.fitting import fit, common



from SpCorr_ZPL_theta_sweep import temporal_filtering ### note that this function uses the same analysis parameters as SPCORRS!!!
import SpCorr_ZPL_theta_sweep
reload(SpCorr_ZPL_theta_sweep)

def get_data_objects(expm_folder,**kw):
        base_folder_lt3 = analysis_params['data_settings']['base_folder_lt3']
        lt3_folder = os.path.join(base_folder_lt3,expm_folder)
        lt3_ssro_folder = os.path.join(base_folder_lt3,'SSROs')
        base_folder_lt4 = analysis_params['data_settings']['base_folder_lt4']
        lt4_folder = os.path.join(base_folder_lt3,expm_folder)
        lt4_ssro_folder = os.path.join(base_folder_lt3,'SSROs')

        filename_str = analysis_params['data_settings']['filenames_for_expms'][expm_folder]

        b_list=tb.latest_data(filename_str,folder= lt3_folder,return_all = True)
        a_list=tb.latest_data(filename_str,folder =lt4_folder,return_all = True)

        if len(b_list) != len(a_list):
            raise(Exception('Different number of files for lt3 and lt4!'))

        ssro_b  = tb.latest_data('SSROCalib', folder = lt3_ssro_folder)
        ssro_a  = tb.latest_data('SSROCalib',  folder = lt4_ssro_folder)


    return a_list,b_list,ssro_a,ssro_b


def run_analysis(contains, **kw):

    a_list,b_list,ssro_a,ssro_b =  get_data_objects(contains,**kw)

    sca_list = []
    for i, folder_a,folder_b in enumerate(zip(a_list,b_list)):
        print 'Processing file {0}'.format(i)
        sca = singleClickAnalysis(folder_a,folder_b)
        sca.process_correlations(**kw)
        sca_list.append(sca)

    analyze_spspcorrs(sca_list,ssro_a,ssro_b,**kw)