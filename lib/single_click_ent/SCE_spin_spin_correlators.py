"""
provides functions to analyze spin-spin correlators between two NV centres
Functions should be executed on the computer that stores the PQ data, i.e. LT4. (othewise use pq_folder = 'xxxx' when calling instances of purify_pq)
Based on the analysis class purify_pq and some functionalities from purify_analysis.py
"""


import numpy as np
from analysis.lib.lde import sscorr; reload(sscorr) ### two qubit SSRO correction
from analysis.lib.purification import purify_pq as ppq; reload(ppq)
import Analysis_params_SCE as analysis_params; reload(analysis_params)
from analysis.lib.pq import pq_tools,pq_plots; reload(pq_plots)
import analysis.lib.purification.purify_analysis as purify_analysis
from analysis.lib.tools import toolbox as tb
from analysis.lib.m2.ssro import ssro
from matplotlib import pyplot as plt
import os

from SpCorr_ZPL_theta_sweep import temporal_filtering ### note that this function uses the same analysis parameters as SPCORRS!!!


def get_data_objects(contains_lt3, contains_lt4,**kw):
    folder=tb.latest_data(contains_lt4,**kw)
    a = ppq.purifyPQAnalysis(folder, hdf5_mode='r')
    lt3_folder=tb.latest_data(contains_lt3,folder =r'Z:\data',**kw)
    b = ppq.purifyPQAnalysis(lt3_folder, hdf5_mode='r')
    ssro_b  = tb.latest_data('SSROCalib', folder =r'Z:\data')
    ### lt3 should use SSRO with MW init. We therefore also supply the microwave transition.
    if 'p' in a.g.attrs['electron_transition']:
        trans_b = 'msp1'
    else:
        trans_b = 'msm1'
    ssro_a  = tb.latest_data('SSROCalib')


    return a,b,ssro_a,ssro_b,trans_b,folder

def analyze_spspcorrs(contains,**kw):
    """
    TO DO: ability to pass analysis_params via kwargs
    Plotting
    """

    if (isinstance(contains, list)):
        contains_lt3 = contains[0]
        contains_lt4 = contains[1]
    else:
        contains_lt3 = contains
        contains_lt4 = contains

    #### kws
    plot_filter = kw.pop('plot_filter',False)
    do_plot = kw.pop('do_plot',False)
    plot_raw_correlators = kw.pop('plot_raw_correlators',False)

    #### get files
    a_lt4,a_lt3,ssro_f_lt4,ssro_f_lt3,trans_lt3,folder_lt4 = get_data_objects(contains_lt3, contains_lt4,**kw)

    ### temporal filtering
    sn_lt,st_fltr_c0,st_fltr_c1 = temporal_filtering(a_lt4,plot_filter = plot_filter)
    st_fltr = np.logical_or(st_fltr_c0,st_fltr_c1)
    st_fltr_c0 = st_fltr_c0[st_fltr]
    st_fltr_c1 = st_fltr_c1[st_fltr]
    
    ### prepare filtered sync numbers
    for a in [a_lt3]:#,a_lt4]: for later functionatlity this could be extended to both files. Not sure if necessary yet.
        a.sn_filtered = sn_lt[st_fltr]
        adwin_filter,adwin_syncs = a.filter_adwin_data_from_pq_syncs(a.sn_filtered) ## adwin syncs within window
        results = a.agrp['ssro_results'].value
        adwin_filt_bool_psi0 = np.in1d(range(len(results)),adwin_filter[st_fltr_c0]) ### convert the adwin filter to boolean
        adwin_filt_bool_psi1 = np.in1d(range(len(results)),adwin_filter[st_fltr_c1]) ### convert the adwin filter to boolean
        if (a.g.attrs['PLU_during_LDE'] != 1):
            adwin_filt_bool= np.array([1]*len(results))
        
        ####################################
        ##### electron RO correlations #####
        ####################################
    correlators_per_sweep_pt_psi0 = get_time_filtered_correlations(a_lt3,a_lt4,adwin_filt_bool_psi0)
    correlators_per_sweep_pt_psi1 = get_time_filtered_correlations(a_lt3,a_lt4,adwin_filt_bool_psi1)
    print np.sum(st_fltr)
    
    a_lt4.get_sweep_pts()
    ### do ROC
    norm_correlators_psi0, norm_correlators_psi0_u = RO_correction_of_correlators(correlators_per_sweep_pt_psi0,a_lt3,a_lt4,ssro_f_lt3,ssro_f_lt4,**kw)
    norm_correlators_psi1, norm_correlators_psi1_u = RO_correction_of_correlators(correlators_per_sweep_pt_psi1,a_lt3,a_lt4,ssro_f_lt3,ssro_f_lt4,**kw)

    ### extract spin-spin expectation value from correlators
    exp_values_psi0,exp_values_psi0_u = get_exp_value_from_spin_spin_corr(norm_correlators_psi0,norm_correlators_psi0_u)
    exp_values_psi1,exp_values_psi1_u = get_exp_value_from_spin_spin_corr(norm_correlators_psi1,norm_correlators_psi1_u)

    if plot_raw_correlators:
        ### transpose the functions to be plotted.

        exp_values_psi0_trans = map(list, zip(*norm_correlators_psi0))
        exp_values_psi0_u_trans= map(list, zip(*norm_correlators_psi0_u))
        exp_values_psi1_trans = map(list, zip(*norm_correlators_psi1))
        exp_values_psi1_u_trans= map(list, zip(*norm_correlators_psi1_u))
        
        labels  = ['11','10','01','00']
        fig = plt.figure()
        ax = plt.subplot()
        for e,e_u,l in zip(exp_values_psi0_trans,exp_values_psi0_u_trans,labels):
            ax.errorbar(a_lt4.sweep_pts,e,e_u,fmt='o',label=l)
        plt.legend()
        ax.set_xlabel(a_lt4.sweep_name)
        ax.set_ylabel('Probability')
        ax.set_ylim([0,1])
        ax.set_title(a_lt4.timestamp+'\n'+a_lt4.measurementstring+ '\n' + 'psi0')
        plt.show()
        fig.savefig(
                    os.path.join(a_lt4.folder, '_probabilty_psi0.png'),
                    format='png')

        labels  = ['11','10','01','00']
        fig = plt.figure()
        ax = plt.subplot()
        for e,e_u,l in zip(exp_values_psi1_trans,exp_values_psi1_u_trans,labels):
            ax.errorbar(a_lt4.sweep_pts,e,e_u,fmt='o',label=l)
        plt.legend()
        ax.set_xlabel(a_lt4.sweep_name)
        ax.set_ylabel('Probability')
        ax.set_ylim([0,1])
        ax.set_title(a_lt4.timestamp+'\n'+a_lt4.measurementstring + '\n' + 'psi1')
        plt.show()
        fig.savefig(
                    os.path.join(a_lt4.folder, '_probabilty_psi1.png'),
                    format='png')

    if do_plot:
        #### plot exp value
        # fig = plt.figure()
        # ax = plt.subplot()
        a_lt4.result_corrected = True
        
        a_lt4.p0 = (np.array(exp_values_psi0)+1)/2.
        a_lt4.u_p0 = exp_values_psi0_u
        ax = a_lt4.plot_results_vs_sweepparam(ylabel = 'Expectation value',ret = 'ax', labels = 'psi0', save=False,**kw)
        a_lt4.readouts = 2
        a_lt4.p0 = (np.array(exp_values_psi1)+1)/2.
        a_lt4.u_p0 = exp_values_psi1_u
        a_lt4.plot_results_vs_sweepparam(ylabel = 'Expectation value', ax = ax, labels = 'psi1', save=True,**kw)
    else:
        return exp_values_psi0,exp_values_psi0_u



def sweep_analysis_parameter():
    """
    Needs to be programmed!!!
    """
    pass

def get_time_filtered_correlations(a_lt3,a_lt4,adwin_filt_bool,**kw):
    """
    Note that the raw data is handled in terms of photon numbers! I.e. 1 == ms=0 and 0 == ms=+-1

    Input: this function takes two data objects and a boolen list (adwin_filt_bool)

    Output: A 2D list with shape (nr_of_sweep_points,4)
            Each entry in the list contains the correlator: (LT3 ms=1 & LT4 ms = 1 // LT3 ms=1 & LT4 ms=0 // LT3 ms=0 & LT4 ms=1 // LT3 ms=0 & LT4 ms=)
            Or in short [11;10;01;00].
            One can then use these entries for two qubit RO Correction.
    """
    verbose = kw.pop('verbose',False)
    ### prepare RO results and sort them according to sweep point
    for a in [a_lt3,a_lt4]:
        a.pts = a.g.attrs['sweep_length']
        a.ssros = a.agrp['ssro_results'].value
        a.readouts = a.g.attrs['nr_of_ROsequences']
        # a.sorted_results = a_ssros.reshape((-1,a.pts,a.readouts))


    ### correlate the ROs with each other by making a boolean filter:
    ### variables here are described in terms of spin states!
    m00 = (a_lt3.ssros == 1)*(a_lt4.ssros == 1)
    m10 = (a_lt3.ssros == 1)*(a_lt4.ssros == 0)
    m01 = (a_lt3.ssros == 0)*(a_lt4.ssros == 1)
    m11 = (a_lt3.ssros == 0)*(a_lt4.ssros == 0)
    
    ### now define unique identifiers for each Ro correlation and recast the correlations into a single array.
    ### As identifieres I choose 1 = index 0 in the output list, i.e. 11; 2 = index 1 in the output list ... and so forth
    RO_correlators =  np.array(len(a_lt3.ssros)*[1])*m11 \
                    + np.array(len(a_lt3.ssros)*[2])*m10 \
                    + np.array(len(a_lt3.ssros)*[3])*m01 \
                    + np.array(len(a_lt3.ssros)*[4])*m00 
    print len(RO_correlators)
    ### now sort the correlators and the adwin fltr according to the sweep pts
    sorted_RO_correlators = RO_correlators.reshape((-1,a_lt3.pts,a_lt3.readouts))
    sorted_adwin_fltr = adwin_filt_bool.reshape((-1,a_lt3.pts,a_lt3.readouts))

    ### from now on: no numpy magic anymore. from here it is brutforce 'for-looping'
    ### (all conceived arrays will have different lengths due to temporal filtering. this break most np methods)
    ### although vstack and hstack would probably work...
    
    return_list = range(a_lt3.pts) ## all of these pts will be substituted with the correlator occurence
    for i in range(a_lt3.pts): 
        correlators_at_sweep_pt = [0,0,0,0]
        for j in [1,2,3,4]: ### loop over the correlator identifiers
            correlators_at_sweep_pt[j-1] = np.sum(np.logical_and(sorted_adwin_fltr[:,i,:],sorted_RO_correlators[:,i,:]==j)) ## exclude adwin filter and do a logical and with the correlator identifier. Then sum over the number of occurences


        return_list[i] = correlators_at_sweep_pt

    return return_list


def RO_correction_of_correlators(correlators_per_sweep_pt,a_lt3,a_lt4,ssro_f_lt3,ssro_f_lt4,**kw):

    verbose = kw.pop('verbose',False)

    ### to store the estimated statistical uncertainty
    norm_correlators_u = range(len(correlators_per_sweep_pt))
    norm_correlators = range(len(correlators_per_sweep_pt))
    ### get ssro_ROC for LT3 --> corresponds to setup B
    F0_LT3,F1_LT3 = get_RO_fildeities(ssro_f_lt3,a_lt3.g.attrs['electron_transition'],a_lt3.g.attrs['E_RO_durations'][0])
    ### get ssro_ROC for LT4 --> corresponds to setup A
    F0_LT4,F1_LT4 = get_RO_fildeities(ssro_f_lt4,a_lt4.g.attrs['electron_transition'],a_lt4.g.attrs['E_RO_durations'][0])

    for i in range(len(correlators_per_sweep_pt)):
        #### note that the function below assumes an error of 1% on the SSRO fidelities!
        norm_correlator,norm_correlator_u = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(correlators_per_sweep_pt[i]),F0_LT4,F0_LT3,F1_LT4,F1_LT3,
                                                                                                verbose = verbose,return_error_bars = True)
        norm_correlators[i] = norm_correlator
        norm_correlators_u[i] = norm_correlator_u

    return norm_correlators,norm_correlators_u

def get_RO_fildeities(ssro_folder,mw_transition,E_RO_duration):

    if 'MWInit' in ssro_folder:
        F_0,u_F0,F_1,u_F1 = ssro.get_SSRO_MWInit_calibration(ssro_folder,E_RO_duration,mw_transition)
    else:
        F_0,u_F0,F_1,u_F1 = ssro.get_SSRO_calibration(ssro_folder,E_RO_duration)

    return F_0,F_1 ### excludiung uncertainties for now.

def get_exp_value_from_spin_spin_corr(norm_correlators,norm_correlators_u):
    """
    Assumes that the input list of correlators has the form
    [
    sweep pt 1: [11,10,01,00],
    sweep pt 2: [11,10,01,00]
    and so forth....
    ]

    Expectation values are then estimated by evaluating the odd parity.
    """


    exp_vals = range(len(norm_correlators))
    exp_vals_u = range(len(norm_correlators))

    for i in range(len(norm_correlators)):
        exp_vals[i]     = 1-2*(norm_correlators[i][1]+norm_correlators[i][2])
        exp_vals_u[i]   = 2*np.sqrt((norm_correlators_u[i][1])**2+(norm_correlators_u[i][2])**2)

    return exp_vals,exp_vals_u