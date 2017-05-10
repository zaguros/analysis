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


def get_data_objects(contains_lt3, contains_lt4,**kw):
    folder=tb.latest_data(contains_lt4,**kw)
    a = ppq.purifyPQAnalysis(folder, hdf5_mode='r')
    lt3_folder=tb.latest_data(contains_lt3,folder =r'Z:\data',**kw)
    b = ppq.purifyPQAnalysis(lt3_folder, hdf5_mode='r')
    ssro_b  = tb.latest_data('103055', folder =r'Z:\data')
    ### lt3 should use SSRO with MW init. We therefore also supply the microwave transition.
    if 'p' in a.g.attrs['electron_transition']:
        trans_b = 'msp1'
    else:
        trans_b = 'msm1'
    ssro_a  = tb.latest_data('SSROCalib')


    return a,b,ssro_a,ssro_b,trans_b,folder

def analyze_spspcorrs(contains,**kw):
    """
    TO DO: ability to pass analysis_params via kwargs. e.g. for temporal filtering
    """

    if (isinstance(contains, list)):
        contains_lt3 = contains[0]
        contains_lt4 = contains[1]
    else:
        contains_lt3 = contains
        contains_lt4 = contains

    #### kws
    plot_filter          = kw.pop('plot_filter',True)
    do_plot              = kw.pop('do_plot',True)
    plot_raw_correlators = kw.pop('plot_raw_correlators',True)
    do_sine_fit          = kw.pop('do_sine_fit',False)
    ret                  = kw.pop('ret',False)
    save_corrs           = kw.pop('save_corrs',False)
    #### get files
    a_lt4,a_lt3,ssro_f_lt4,ssro_f_lt3,trans_lt3,folder_lt4 = get_data_objects(contains_lt3, contains_lt4,**kw)

    ### temporal filtering
    sn_lt,st_fltr_c0,st_fltr_c1 = temporal_filtering(a_lt4,plot_filter = plot_filter)
    st_fltr = np.logical_or(st_fltr_c0,st_fltr_c1)
    st_fltr_c0 = st_fltr_c0[st_fltr]
    st_fltr_c1 = st_fltr_c1[st_fltr]
    
    print 'Time filtered events ', np.sum(st_fltr_c0), np.sum(st_fltr_c1)
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
    adwin_record_filter = filter_on_adwin_parameters(a_lt3 = a_lt3,a_lt4 = a_lt4)
    adwin_filt_bool_psi0 = np.logical_and(adwin_record_filter,adwin_filt_bool_psi0)
    adwin_filt_bool_psi1 = np.logical_and(adwin_record_filter,adwin_filt_bool_psi1)

    print 'Adwin filtered events ', np.sum(adwin_filt_bool_psi0), np.sum(adwin_filt_bool_psi1)

    p0 = []
    p0_u = []
    counts_per_pt = []
    a_lt4.get_sweep_pts()
    for i, adwin_filt_bool in zip([0,1],[adwin_filt_bool_psi0,adwin_filt_bool_psi1]):

        correlators_per_sweep_pt = get_time_filtered_correlations(a_lt3,a_lt4,adwin_filt_bool)
        counts_per_pt.append(np.sum(correlators_per_sweep_pt, axis=1))
        ### do ROC
        norm_correlators, norm_correlators_u = RO_correction_of_correlators(correlators_per_sweep_pt,
                                                                    a_lt3,a_lt4,ssro_f_lt3,ssro_f_lt4,**kw)

        ### extract spin-spin expectation value from correlators
        exp_values,exp_values_u = get_exp_value_from_spin_spin_corr(norm_correlators,norm_correlators_u)

        p0.append(exp_values)
        p0_u.append(exp_values_u)

        if plot_raw_correlators:
            ### transpose the functions to be plotted.
            correlators_trans = map(list, zip(*norm_correlators))
            correlators_u_trans= map(list, zip(*norm_correlators_u))

            labels  = ['11','10','01','00']
            fig = plt.figure()
            ax = plt.subplot()
            for e,e_u,l in zip(correlators_trans,correlators_u_trans,labels):
                ax.errorbar(a_lt4.sweep_pts,e,e_u,fmt='o',label=l)
            plt.legend()
            ax.set_xlabel(a_lt4.sweep_name)
            ax.set_ylabel('Probability')
            ax.set_ylim([0,1])
            ax.set_title(a_lt4.timestamp+'\n'+a_lt4.measurementstring+ '\n' + 'psi'+str(i))
            plt.show()
            fig.savefig(
                        os.path.join(a_lt4.folder, '_probabilty_psi_' +str(i) +'.png'),format='png')
    
    a_lt4.result_corrected = True
    a_lt4.result_correlation_corrected = True
    a_lt4.p_correlations = np.transpose(np.array(p0)) ##XXXXXX CHECK IF THIS REALLY DOES THE THING without inversion of x
    a_lt4.u_p_correlations = np.transpose(np.array(p0_u))
    a_lt4.readouts = len(p0)
    if do_plot:
        #### plot exp value
        a_lt4.correlation_names = ['psi0','psi1']
        fig = a_lt4.plot_results_vs_sweepparam(   ylabel = 'Expectation value', 
                                            mode = 'correlations',
                                            ylim = [-1.05,1.05],
                                            save=False,
                                            ret = 'fig',
                                            **kw)
        if do_sine_fit:
            ax = fig.gca()
            phi = []
            for jj,p in zip(range(len(p0)),np.array(p0)):
                y = np.array(p)
                x = a_lt4.sweep_pts

                info_x = ax.get_xlim()[1] + (ax.get_xlim()[-1]-ax.get_xlim()[0])*0.02
                info_y = ax.get_ylim()[0] + (ax.get_ylim()[-1]-ax.get_ylim()[0])*0.02 + 0.5*jj
                info_xy = [info_x,info_y]
                ### estimate guess:
                phi.append(fit_and_plot_sine(x,y,fig.gca(),info_xy = info_xy,**kw))
            phi[0] = np.mod(-phi[0],360)
            phi[1] = np.mod(-phi[1]+180,360)
            avg_phi = np.mean(phi)
            print 'Phi angle ', phi
            print 'Avg phi ', avg_phi
        fig.savefig(
            os.path.join(a_lt4.folder, '{}_vs_sweepparam.'.format('correlations') + a_lt4.plot_format),
            format=a_lt4.plot_format)

    if save_corrs:

        name = os.path.join(a_lt4.folder, 'correlations.h5')
        with h5py.File(name, 'w') as hf:
            
            hf.create_dataset('correlations', data=a_lt4.p_correlations)
            hf.create_dataset('correlations_u', data=a_lt4.u_p_correlations)
            hf.create_dataset('norm_correlators', data=norm_correlators)
            hf.create_dataset('norm_correlators_u', data=norm_correlators_u)
            hf.create_dataset('counts_per_pt', data=counts_per_pt)
            

        
    if ret:
        return a_lt4.p_correlations,a_lt4.u_p_correlations



def sweep_analysis_parameter():
    """
    TODO. See SPCorr, Sweep theta.
    """
    pass

def fit_and_plot_sine(x,y,ax,info_xy,**kw):

    do_print_fit_result = kw.pop('do_print_fit_result', False)
    ### estimate some guess from the data:
    g_a = 0.0#np.sum(y)/len(y)
    g_A = np.amax(y)-g_a
    g_phi = x[np.argmax(y)] 
    ### frequency guess is hardcoded as this is mainly intended for specific entangling oscillations.
    g_f = 1/360.
    p0, fitfunc,fitfunc_str = common.fit_cos(g_f,g_a,g_A,g_phi)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print_fit_result,
                     ret=True,fixed=[0,1])

    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax, 
                    plot_data=False,print_info = True,info_xy = info_xy)

    return fit_result['params_dict']['phi']


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
    
    RO_correlators = RO_correlators[:(a.g.attrs['sweep_length']*(len(RO_correlators)/a.g.attrs['sweep_length']))]
    adwin_filt_bool =  adwin_filt_bool[:(a.g.attrs['sweep_length']*(len(RO_correlators)/a.g.attrs['sweep_length']))]

    
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

def filter_on_adwin_parameters(a_lt3,a_lt4):
    """
    creates a filter on the adwin RO array which is based on filters as specified in the analysis parameters
    Each filter parameters should exist in the hdf5 file as an adwin group and should have the same length as the 
    number of ssros. Prime example: CR_after and CR_before
    """
    fltr = np.array([True]*len(a_lt3.agrp['ssro_results'].value)) ### initially everything true


    for a,suffix in zip([a_lt3,a_lt4],['_lt3','_lt4']): ### loop over both files

        for key in analysis_params.SPSP_fltr_adwin_settings['list_of_adwin_params']: ### loop over the list of filter parameters
            [minimum,maximum] = analysis_params.SPSP_fltr_adwin_settings['fltr_dict'+suffix][key] 
            values = a.agrp[key].value
            fltr = np.logical_and(fltr,(values > minimum) & ( values < maximum)) ### update filter

    return fltr

def RO_correction_of_correlators(correlators_per_sweep_pt,a_lt3,a_lt4,ssro_f_lt3,ssro_f_lt4,**kw):

    verbose = kw.pop('verbose',False)
    do_ROC = kw.pop('do_ROC',True)
    ### to store the estimated statistical uncertainty
    norm_correlators_u = range(len(correlators_per_sweep_pt))
    norm_correlators = range(len(correlators_per_sweep_pt))
    ### get ssro_ROC for LT3 --> corresponds to setup B
    if do_ROC:
        F0_LT3,F1_LT3 = get_RO_fildeities(ssro_f_lt3,a_lt3.g.attrs['electron_transition'],a_lt3.g.attrs['E_RO_durations'][0])
        ### get ssro_ROC for LT4 --> corresponds to setup A
        F0_LT4,F1_LT4 = get_RO_fildeities(ssro_f_lt4,a_lt4.g.attrs['electron_transition'],a_lt4.g.attrs['E_RO_durations'][0])

        for i in range(len(correlators_per_sweep_pt)):
            #### note that the function below assumes an error of 1% on the SSRO fidelities!
            norm_correlator,norm_correlator_u = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(correlators_per_sweep_pt[i]),F0_LT4,F0_LT3,F1_LT4,F1_LT3,
                                                                                                    verbose = verbose,return_error_bars = True)
            norm_correlators[i] = np.squeeze(norm_correlator)
            norm_correlators_u[i] = norm_correlator_u

    else: 
        for i in range(len(correlators_per_sweep_pt)):
            norm_correlators[i] = np.array(correlators_per_sweep_pt[i],dtype = np.float64)/np.sum(np.array(correlators_per_sweep_pt[i]))
            norm_correlators_u[i] = np.array(np.sqrt(correlators_per_sweep_pt[i]),dtype = np.float64)/np.sum(np.array(correlators_per_sweep_pt[i]))

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