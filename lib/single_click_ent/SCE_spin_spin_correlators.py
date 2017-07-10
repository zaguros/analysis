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
from analysis.lib.pq import pq_tools,pq_plots; reload(pq_plots)
from analysis.lib.tools import plot; reload(plot)
import analysis.lib.purification.purify_analysis as purify_analysis
from analysis.lib.tools import toolbox as tb; reload(tb)
from analysis.lib.m2.ssro import ssro
from matplotlib import pyplot as plt
import h5py
from analysis.lib.fitting import fit, common



from SpCorr_ZPL_theta_sweep import temporal_filtering ### note that this function uses the same analysis parameters as SPCORRS!!!
import SpCorr_ZPL_theta_sweep
reload(SpCorr_ZPL_theta_sweep)

class twoSetupAnalysis:
    # Skeleton class combining the files from two setups into one analysis

    def __init__(self,folder_a,folder_b):

        self.a = ppq.purifyPQAnalysis(folder_a, hdf5_mode='r') 
        self.b = ppq.purifyPQAnalysis(folder_b, hdf5_mode='r')

class singleClickAnalysis(twoSetupAnalysis):
    # Each run of the single click entanglement expm is analyzed in this class. It pulls in the results from each setup
    # along with PQ events etc.

    def process_correlations(self,**kw):

        verbose = kw.get('verbose',False)
        ignore_HH = kw.pop('ignore_HH',False)

        # print self.a.g
        # print self.a.g.attrs['sweep_length']
        # print self.b.g.attrs['sweep_length']
        if self.a.g.attrs['sweep_length'] != self.b.g.attrs['sweep_length']:
            raise(Exception('Different sweep lengths for lt3 and lt4!'))
        if self.a.agrp['completed_reps'].value  != self.b.agrp['completed_reps'].value :
            print 'Different reps for lt3 and lt4!'

        if np.any(self.a.agrp['counted_awg_reps'].value != self.b.agrp['counted_awg_reps'].value):
            print 'Fix your syncs you fool'

        self.completed_reps = np.min([self.a.agrp['completed_reps'].value ,self.b.agrp['completed_reps'].value ])
            
        self.sweep_pts = self.a.g.attrs['sweep_pts'] # Specified in lt4
        self.sweep_length = self.a.g.attrs['sweep_length']
        
        ### temporal filtering of HH data
        sn_lt,st_fltr_c0,st_fltr_c1 = temporal_filtering(self.a,**kw)

        st_fltr = np.logical_or(st_fltr_c0,st_fltr_c1)
        self.st_fltr_c0 = st_fltr_c0[st_fltr]
        self.st_fltr_c1 = st_fltr_c1[st_fltr]
        
        ### prepare filtered sync numbers
        self.sn_filtered = sn_lt[st_fltr]

        if verbose: print 'Time filtered events ', np.sum(st_fltr_c0), np.sum(st_fltr_c1)

        ####################################
        ##### electron RO correlations #####
        ####################################
        # Apply adwin filter to data and combine with earlier temporal fitering
        self.filter_on_adwin_parameters(**kw) 
        adwin_filter,adwin_syncs = self.a.filter_adwin_data_from_pq_syncs(self.sn_filtered) ## adwin syncs within window
        # if ignore_HH:
        #     adwin_filter = np.array(range(self.completed_reps))
        #     print np.shape(self.st_fltr_c0)
        #     print self.st_fltr_c0
        self.combined_filt = np.array([np.logical_and(self.adwin_fltr,np.in1d(range(self.completed_reps),adwin_filter[st_fltr])) for st_fltr in [self.st_fltr_c0, self.st_fltr_c1]]) ### convert the adwin filter to boolean
        self.failed_events = np.logical_and(self.adwin_fltr,np.logical_not(np.logical_or(self.combined_filt[0],self.combined_filt[1]))) # where didnt get a plu signal

        if verbose: print 'Num failed events ', np.sum(self.failed_events)

        if ignore_HH:
            psi0_or_psi1_failed = np.random.rand(len(self.combined_filt[0])) > 0.5
            self.combined_filt[0] = np.logical_or(np.logical_and(self.failed_events,psi0_or_psi1_failed == 0),self.combined_filt[0])
            self.combined_filt[1] = np.logical_or(np.logical_and(self.failed_events,psi0_or_psi1_failed == 1),self.combined_filt[1])
            # self.combined_filt[0] = self.failed_events

        if verbose: print 'Adwin filtered events ', np.sum(self.combined_filt[0]), np.sum(self.combined_filt[1])

        self.get_correlations()
        self.counts_per_pt = (np.sum(self.correlators_per_sweep_pt, axis=2))
        self.get_tail_counts()


    def filter_on_adwin_parameters(self,**kw):
        """
        creates a filter on the adwin RO array which is based on filters as specified in the analysis parameters
        Each filter parameters should exist in the hdf5 file as an adwin group and should have the same length as the 
        number of ssros. Prime example: CR_after and CR_before

        Can manually override the dict by passing the arguement adwin_filter_params in the form {'setup_key (e.g. lt3)' : {'param (e.g. CR_after)' : [enabled, min, max]}}
        """
        
        filter_params = kw.pop('adwin_filter_params',{})
        if len(filter_params):
            old_params = analysis_params.SPSP_fltr_adwin_settings
        
        for setup_key,setup_dict in filter_params.iteritems():
            for key,params in setup_dict.iteritems():
                analysis_params.SPSP_fltr_adwin_settings['fltr_dict_'+setup_key][key] = params

        self.adwin_fltr = np.array([True]*self.completed_reps) ### initially everything true

        for a,suffix in zip([self.a,self.b],['lt4','lt3']): ### loop over both files
            for key,val in analysis_params.SPSP_fltr_adwin_settings['fltr_dict_'+suffix].iteritems(): ### loop over the list of filter parameters
                [filter_on,minimum,maximum] = val

                if filter_on:
                    if key == 'repetition_number':
                        values = np.array([i for i in range(self.completed_reps/self.sweep_length) for _ in range(self.sweep_length)]) ### Make an array of values corresponding to the current rep
                    elif key == 'LDE_attempts':
                        reps = a.agrp['counted_awg_reps'].value # I feel like a dirty person
                        reps = reps - np.append([0],reps[:-1]) ### Calculate reps since last click
                        values = reps % a.g.attrs['LDE_attempts'] ### Take mod with 
                    elif key == 'pst_msmt_phase':
                        v_1 = a.agrp['sampling_counts_1'].value
                        v_2 = a.agrp['sampling_counts_2'].value
                        g_0 = a.agrp.attrs['Phase_Msmt_g_0']
                        visibility = a.agrp.attrs['Phase_Msmt_Vis']
                        
                        cosvals = [2*(float(n0)/(float(n0)+float(n1)*g_0)-0.5)*visibility for n0,n1 in zip(v_1,v_2)]
                        cosvals = [cosval if np.abs(cosval) < 1 else (1.0 * np.sign(cosval)) for cosval in cosvals]
                        values = 180*np.arccos(cosvals)/np.pi
                    
                    else:
                        values = a.agrp[key].value
                    values = values[:self.completed_reps]
                    self.adwin_fltr = np.logical_and(self.adwin_fltr,(values >= minimum) & ( values <= maximum)) ### update filter

        if len(filter_params):
            analysis_params.SPSP_fltr_adwin_settings = old_params
        


    def get_correlations(self,**kw):
        """
        Note that the raw data is handled in terms of photon numbers! I.e. 1 == ms=0 and 0 == ms=+-1

        Output: A 2D list with shape (nr_of_sweep_points,4)
                Each entry in the list contains the correlator: (LT3 ms=1 & LT4 ms = 1 // LT3 ms=1 & LT4 ms=0 // LT3 ms=0 & LT4 ms=1 // LT3 ms=0 & LT4 ms=)
                Or in short [11;10;01;00].
                One can then use these entries for two qubit RO Correction.
        """
        verbose = kw.pop('verbose',False)
        ### prepare RO results and sort them according to sweep point
        for a in [self.a,self.b]:
            a.ssros = a.agrp['ssro_results'].value
            a.ssros = a.ssros[:self.completed_reps]
            
        ### correlate the ROs with each other by making a boolean filter:
        ### variables here are described in terms of spin states!
        m00 = (self.b.ssros == 1)*(self.a.ssros == 1)
        m10 = (self.b.ssros == 1)*(self.a.ssros == 0)
        m01 = (self.b.ssros == 0)*(self.a.ssros == 1)
        m11 = (self.b.ssros == 0)*(self.a.ssros == 0)
        
        ### now define unique identifiers for each Ro correlation and recast the correlations into a single array.
        ### As identifieres I choose 1 = index 0 in the output list, i.e. 11; 2 = index 1 in the output list ... and so forth
        RO_correlators =  np.array(self.completed_reps*[1])*m11 \
                        + np.array(self.completed_reps*[2])*m10 \
                        + np.array(self.completed_reps*[3])*m01 \
                        + np.array(self.completed_reps*[4])*m00 
        ### PH - added to make sure that has a full set of repetitions
        completed_points = (self.sweep_length*(self.completed_reps/self.sweep_length))
        RO_correlators = RO_correlators[:completed_points]
  
        self.combined_filt =  self.combined_filt[:,:completed_points]
        
        ### now sort the correlators and the adwin fltr according to the sweep pts
        sorted_RO_correlators = RO_correlators.reshape((-1,self.sweep_length))
        
        self.correlators_per_sweep_pt = np.zeros([2,self.sweep_length,4])

        for k, combined_filt in zip([0,1],self.combined_filt):
            sorted_combined_filt = combined_filt.reshape((-1,self.sweep_length))
            ### from now on: no numpy magic anymore. from here it is brutforce 'for-looping'
            ### (all conceived arrays will have different lengths due to temporal filtering. this break most np methods)
            ### although vstack and hstack would probably work...
            
            for i in range(self.sweep_length): 
                for j in [1,2,3,4]: ### loop over the correlator identifiers
                    self.correlators_per_sweep_pt[k][i][j-1] = np.sum(np.logical_and(sorted_combined_filt[:,i],sorted_RO_correlators[:,i]==j)) ## exclude adwin filter and do a logical and with the correlator identifier. Then sum over the number of occurences

    def get_tail_counts(self):

        counted_awg_reps = self.a.agrp['counted_awg_reps'].value

        awg_reps_per_attempt = counted_awg_reps - np.append([0],counted_awg_reps[:-1])
        ### PH - added to make sure that has a full set of repetitions
        awg_reps_per_attempt = awg_reps_per_attempt[: (self.sweep_length*(self.completed_reps/self.sweep_length))]
        repsPerClick = np.sum(np.reshape(awg_reps_per_attempt,[-1,self.sweep_length]),axis=0)
        counts_per_pt = (np.sum(self.counts_per_pt,axis=0)).astype(np.float)
        self.tail_per_pt = 10**4 * counts_per_pt/repsPerClick
        self.tail_per_pt_u = 10**4 * np.sqrt(counts_per_pt)/repsPerClick
        
def run_analysis(contains, **kw):
    ### This is a simple helper function for when only analysing one file

    a_list,b_list,ssro_a,ssro_b =  get_data_objects(contains,**kw)

    sca_list = []
    for i, (folder_a,folder_b) in enumerate(zip(a_list,b_list)):
        sca = singleClickAnalysis(folder_a,folder_b)
        sca.process_correlations(**kw)
        sca_list.append(sca)

    analyze_spspcorrs(sca_list,ssro_a,ssro_b,**kw)


def get_data_objects(contains,**kw):
    ### Pull data objects when analyzing only one file at a time

    analysis_computer = kw.pop('analysis_computer', 'lt4')

    if (isinstance(contains, list)):
        if len(contains) == 2:
            contains_lt3 = contains[0]
            contains_lt4 = contains[1]
            contains_lt3_ssro = 'SSROCalib'
            contains_lt4_ssro = 'SSROCalib'
        elif len(contains) == 4:
            contains_lt3 = contains[0]
            contains_lt4 = contains[1]
            contains_lt3_ssro = contains[2]
            contains_lt4_ssro = contains[3]
    else:
        contains_lt3 = contains
        contains_lt4 = contains
        contains_lt3_ssro = 'SSROCalib'
        contains_lt4_ssro = 'SSROCalib'

    if analysis_computer == 'lt4':
        folder_a =tb.latest_data(contains_lt4,**kw)
        folder_b =tb.latest_data(contains_lt3,folder =r'Z:\data',**kw)

        folder_a = folder_a if isinstance(folder_a, list) else [folder_a]
        folder_b = folder_b if isinstance(folder_b, list) else [folder_b]
        
        ssro_b  = tb.latest_data(contains_lt3_ssro, folder =r'Z:\data')
        ssro_a  = tb.latest_data(contains_lt4_ssro)

    if analysis_computer == 'lt3_analysis':
        folder_a =tb.latest_data(contains_lt4,folder= r'X:\data',**kw)
        folder_b =tb.latest_data(contains_lt3,folder =r'Y:\data',**kw)
        ssro_b  = tb.latest_data(contains_lt3_ssro, folder =r'Y:\data')
        ssro_a  = tb.latest_data(contains_lt4_ssro,  folder =r'X:\data')

    return folder_a,folder_b,ssro_a,ssro_b

def analyze_spspcorrs(singleClickAnalyses,ssro_a,ssro_b,**kw):
    """
    This function takes the results from each run (or can be only one run) and
    combines them. It then analyses them in every which way and does ROC etc.
    """

    #### kws
    plot_temporal_filter = kw.pop('plot_temporal_filter', False)
    plot_correlations    = kw.pop('plot_correlations',True)
    plot_raw_correlators = kw.pop('plot_raw_correlators',False)
    plot_tail            = kw.pop('plot_tail',False)
    verbose              = kw.pop('verbose',False)
    do_sine_fit          = kw.pop('do_sine_fit',False)
    combine_correlation_data = kw.pop('combine_correlation_data', False)
    flip_psi1            = kw.pop('flip_psi1',False)
    abs_corrs            = kw.pop('abs_corrs', False)
    ret                  = kw.pop('ret',False)
    save_corrs           = kw.pop('save_corrs',False)
    save_figs            = kw.pop('save_figs',False)

    if isinstance(singleClickAnalyses,singleClickAnalysis): # Make sure in an array

        sweep_name = singleClickAnalyses.a.g.attrs['sweep_name']
        timestamp = singleClickAnalyses.a.timestamp
        measurementstring = singleClickAnalyses.a.measurementstring
        save_folder = singleClickAnalyses.a.folder

        singleClickAnalyses = [singleClickAnalyses]
    else:
        print 'Gotta fix this!'

        sweep_name = singleClickAnalyses[0].a.g.attrs['sweep_name']
        timestamp = singleClickAnalyses[0].a.timestamp
        measurementstring = singleClickAnalyses[0].a.measurementstring
        save_folder = singleClickAnalyses[0].a.folder
        
    
    # Combine files together
    num_files = len(singleClickAnalyses)
    sweep_pts = singleClickAnalyses[0].sweep_pts
    sweep_length = singleClickAnalyses[0].sweep_length
    gen_sweep_pts = singleClickAnalyses[0].a.g.attrs['general_sweep_pts']
    
    electron_transitions = [singleClickAnalyses[0].a.g.attrs['electron_transition'], singleClickAnalyses[0].b.g.attrs['electron_transition']]
    E_RO_durations = [singleClickAnalyses[0].a.g.attrs['E_RO_durations'][0], singleClickAnalyses[0].b.g.attrs['E_RO_durations'][0]]

    correlators_per_sweep_pt = np.zeros([2,sweep_length,4])
    counts_per_pt = np.zeros([2,sweep_length])
    tail_per_pt = np.zeros(sweep_length)
    tail_per_pt_u = np.zeros(sweep_length)

    for sca in singleClickAnalyses:
    #     if (sca.sweep_pts != sweep_pts).any():
    #         raise(Exception('Mismatched sweep pts between runs!'))
        
        correlators_per_sweep_pt += sca.correlators_per_sweep_pt
        counts_per_pt += sca.counts_per_pt
        tail_per_pt += sca.tail_per_pt
        tail_per_pt_u += sca.tail_per_pt_u**2

    tail_per_pt = tail_per_pt/num_files
    tail_per_pt_u = np.sqrt(tail_per_pt_u)/num_files

    if plot_temporal_filter:
        plot_ph_hist_and_fltr(singleClickAnalyses)

    if plot_tail:
        fig = singleClickAnalyses[0].a.default_fig(figsize=(6,4))
        ax = singleClickAnalyses[0].a.default_ax(fig)
        plt.errorbar(sweep_pts, tail_per_pt,
                         fmt='o', yerr=tail_per_pt_u,markersize=6,capsize=3)
        xlims = plt.xlim()
        newxlims = [xlims[0] - 0.1*(xlims[1] - xlims[0]), xlims[1] + 0.1*(xlims[1] - xlims[0])]
        plt.xlim(newxlims)

    ### do ROC
    norm_correlators, norm_correlators_u = RO_correction_of_correlators(correlators_per_sweep_pt,
                                                                    electron_transitions,E_RO_durations,ssro_a,ssro_b,**kw)

    p0 = []
    p0_u = []

    for i,corrs,corrs_u in zip([0,1],norm_correlators,norm_correlators_u):

        ### extract spin-spin expectation value from correlators
        exp_values,exp_values_u = get_exp_value_from_spin_spin_corr(corrs,corrs_u)

        if flip_psi1 and i == 1:
            exp_values = -1.0*np.array(exp_values)

        p0.append(exp_values)
        p0_u.append(exp_values_u)

        if plot_raw_correlators:
            ### transpose the functions to be plotted.
            correlators_trans = map(list, zip(*corrs))
            correlators_u_trans= map(list, zip(*corrs_u))

            labels  = ['11','10','01','00']
            fig = plt.figure()
            ax = plt.subplot()
            for e,e_u,l in zip(correlators_trans,correlators_u_trans,labels):
                ax.errorbar(sweep_pts,e,e_u,fmt='o',label=l)
            plt.legend()
            ax.set_xlabel(sweep_name)
            ax.set_ylabel('Probability')
            ax.set_ylim([0,1])
            ax.set_title(timestamp+'\n'+measurementstring+ '\n' + 'psi'+str(i))
            xlims = plt.xlim()
            newxlims = [xlims[0] - 0.1*(xlims[1] - xlims[0]), xlims[1] + 0.1*(xlims[1] - xlims[0])]
            plt.xlim(newxlims)
            plt.show()
            if save_figs:
                fig.savefig(
                        os.path.join(save_folder, 'probabilty_psi_' +str(i) +'.png'),format='png')
    
    p0 = np.array(p0)
    p0_u = np.array(p0_u)



    if combine_correlation_data:
        plot_p0 = np.array([(p0[0] + p0[1])/2.0]) # Need to use with flip_psi1 usually 
        plot_p0_u = np.array([np.sqrt(p0_u[0]**2 + p0_u[1]**2)/2.0])
        labels = ['combined data']
    else:
        plot_p0 = p0
        plot_p0_u = p0_u
        labels = ['psi0','psi1']

    if plot_correlations:

        
        
        if len(gen_sweep_pts) == 0:
            fig = singleClickAnalyses[0].a.default_fig(figsize=(6,8))
            ax = []
            for n in range(len(singleClickAnalyses[0].a.g.attrs['general_sweep_pts2'])):
                ax.append(fig.add_subplot(len(singleClickAnalyses[0].a.g.attrs['general_sweep_pts2']),1,n+1))
            ax[0].set_title(singleClickAnalyses[0].a.timestamp+'\n'+singleClickAnalyses[0].a.measurementstring)
        else:
            fig = singleClickAnalyses[0].a.default_fig(figsize=(6,4))
            ax = singleClickAnalyses[0].a.default_ax(fig)


    phi = []
    phi_u = []

    for jj,(p,p_u) in enumerate(zip(plot_p0,plot_p0_u)):
        x = sweep_pts

        if plot_correlations:
            if abs_corrs:
                p_plot = np.abs(p)
                lims =  [0.,1.05]
            else:
                p_plot = p
                lims =  [-1.05,1.05]

            # Cleverness for if is a 2d sweep
            if len(gen_sweep_pts) == 0:
                general_sweep_pts1 = singleClickAnalyses[0].a.g.attrs['general_sweep_pts1']
                general_sweep_pts2 = singleClickAnalyses[0].a.g.attrs['general_sweep_pts2']
                x  =general_sweep_pts1

                p_plot = np.reshape(p_plot, [len(general_sweep_pts1),len(general_sweep_pts2)])
                p_u = np.reshape(p_u, [len(general_sweep_pts1),len(general_sweep_pts2)])
                for z, (lab, pp, pu) in enumerate(zip(general_sweep_pts2, p_plot.T, p_u.T)):
                    ax[z].errorbar(x, pp,
                         fmt='o', yerr=pu,markersize=6,capsize=3)
                    ax[z].set_ylabel('Correlations')
                    ax[z].set_xlabel(sweep_name)
                    ax[z].axhspan(0,1,fill=False,ls='dotted')
                print (1+np.sum(np.abs(p_plot),axis =1))/4
                print np.sqrt(np.sum(p_u**2,axis = 1))/4
                print p_plot
            else:
                ax.errorbar(x, p_plot,
                         fmt='o', yerr=p_u,markersize=6,capsize=3,label= labels[jj])
                ax.set_ylabel('Correlations')
                ax.set_xlabel(sweep_name)
                ax.axhspan(0,1,fill=False,ls='dotted')
                plt.ylim(lims)
                xlims = plt.xlim()
                newxlims = [xlims[0] - 0.1*(xlims[1] - xlims[0]), xlims[1] + 0.1*(xlims[1] - xlims[0])]
                plt.xlim(newxlims)
                print 'fidelity', (1+np.sum(np.abs(p_plot)))/4,  np.sqrt(np.sum(p_u**2))/4
                plt.legend()

        if do_sine_fit:    
            g_a = 0.0
            g_A = np.amax(p)-g_a
            g_phi = x[np.argmax(p)] 
            ### frequency guess is hardcoded as this is mainly intended for specific entangling oscillations.
            g_f = 1/360.
            p0s, fitfunc,fitfunc_str = common.fit_cos(g_f,g_a,g_A,g_phi)

            fit_result = fit.fit1d(x,p, None, p0=p0s, fitfunc=fitfunc,
                 ret=True,fixed=[0,1])

            if plot_correlations:
                plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax, 
                    plot_data=False,print_info = False)
                print 'A,phi ', fit_result['params_dict']['A'], np.mod(-fit_result['params_dict']['phi'],360)
            
            phi_t = fit_result['params_dict']['phi']
            if fit_result['params_dict']['A'] < 0:
                phi_t = phi_t + 180

            phi.append(np.mod(-phi_t,360))
            phi_u.append(fit_result['error_dict']['phi'])

    if do_sine_fit:
        if not(combine_correlation_data):
            phi[1] = np.mod(-phi[1]+180,360)
        
        phi = np.mean(phi)
        phi_u = np.sqrt(np.mean(np.array(phi_u)**2))
        print 'Phi angle ', phi, phi_u
            
    if plot_correlations and save_figs:
        fig.savefig(
            os.path.join(save_folder, '{}_vs_sweepparam.'.format('correlations') + a.plot_format),
            format=a.plot_format)

    # from analysis.lib.single_click_ent import Phase_stab_and_meas as psm; reload(psm)
    # psm.analyze_phase(contains = 'Xsweep', mode  = 'only_stab', plot_zoomed  = [0,20],start_rep_no = 3)

    if save_corrs:

        name = os.path.join(save_folder, 'correlations.h5')
        with h5py.File(name, 'w') as hf:
            hf.create_dataset('sweep_pts', data=sweep_pts)
            hf.create_dataset('correlations', data=p0)
            hf.create_dataset('correlations_u', data=p0_u)
            hf.create_dataset('norm_correlators', data=norm_correlators)
            hf.create_dataset('norm_correlators_u', data=norm_correlators_u)
            hf.create_dataset('counts_per_pt', data=counts_per_pt)
            hf.create_dataset('tail_counts', data=tail_per_pt)
            hf.create_dataset('tail_counts_u', data=tail_per_pt_u)
  
    if ret:
        return sweep_pts,p0,p0_u,norm_correlators,norm_correlators_u,counts_per_pt,tail_per_pt, tail_per_pt_u, phi, phi_u







def RO_correction_of_correlators(correlators_per_sweep_pt,electron_transitions,E_RO_durations,ssro_a,ssro_b,**kw):

    verbose = kw.pop('verbose',False)
    do_ROC = kw.pop('do_ROC',True)

    ### to store the estimated statistical uncertainty
    norm_correlators_u = np.zeros(np.shape(correlators_per_sweep_pt))
    norm_correlators = np.zeros(np.shape(correlators_per_sweep_pt))

    
    if do_ROC:
        ### get ssro_ROC for LT3 --> corresponds to setup B
        F0_LT3,F1_LT3 = get_RO_fildeities(ssro_b,electron_transitions[1],E_RO_durations[1])
        ### get ssro_ROC for LT4 --> corresponds to setup A
        F0_LT4,F1_LT4 = get_RO_fildeities(ssro_a,electron_transitions[0],E_RO_durations[0])

        for j, corrs in zip([0,1],correlators_per_sweep_pt):
            for i, corr in enumerate(corrs):
                #### note that the function below assumes an error of 1% on the SSRO fidelities!
                norm_correlator,norm_correlator_u = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(corr),F0_LT4,F0_LT3,F1_LT4,F1_LT3,
                                                                                                        verbose = verbose,return_error_bars = True)
                norm_correlators[j,i] = np.squeeze(norm_correlator)
                norm_correlators_u[j,i] = norm_correlator_u

    else: 
        for j, corrs in zip([0,1],correlators_per_sweep_pt):
            for i, corr in enumerate(corrs):
                norm_correlators[j,i] = np.array(corr,dtype = np.float64)/np.sum(np.array(corr))
                norm_correlators_u[j,i] = np.array(np.sqrt(corr),dtype = np.float64)/np.sum(np.array(corr))

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


def get_sweep_analysis_results(singleClick,ssro_a,ssro_b,parameter_name,parameter_range,parameter_kind,**kw):
    
    
    min_or_max = kw.pop('min_or_max','min')

    y = parameter_range
    z = []
    u_z = []
    tail = []
    u_tail = []

    x = singleClick.a.g.attrs['sweep_pts']
    xlabel = singleClick.a.g.attrs['sweep_name']
    
    for p in parameter_range:
        if parameter_kind == 'SPCorr':
            analysis_params.SPCorr_settings[parameter_name] = p
        elif parameter_kind == 'fltr_dict_lt3' or 'fltr_dict_lt4':
            if min_or_max == 'min':
                analysis_params.SPSP_fltr_adwin_settings[parameter_kind][parameter_name] = [1,p,analysis_params.SPSP_fltr_adwin_settings[parameter_kind][parameter_name][2]]
            elif min_or_max == 'max':
                analysis_params.SPSP_fltr_adwin_settings[parameter_kind][parameter_name] = [1,analysis_params.SPSP_fltr_adwin_settings[parameter_kind][parameter_name][1],p]
            elif min_or_max == 'both':
                analysis_params.SPSP_fltr_adwin_settings[parameter_kind][parameter_name] = [1,p[0],p[1]]
        
        singleClick.process_correlations(**kw)
        sweep_pts,p_correlations,u_p_correlations,norm_correlators,norm_correlators_u,counts_per_pt,tail_cts, tail_cts_u, phi = analyze_spspcorrs(singleClick,ssro_a,ssro_b,plot_correlations = False,ret = True,**kw)

        ### store sweep results
        z.append(p_correlations)
        u_z.append(u_p_correlations)
        tail.append(tail_cts)
        u_tail.append(tail_cts_u)

    #restore analysis parameters
    reload(analysis_params)
    return xlabel,x,y,z,u_z,tail,u_tail

### helper functions



def sweep_analysis_parameter(contains,parameter_name, parameter_range, **kw):
    
    parameter_kind = kw.pop('parameter_kind','SPCorr')

    plot_tail = kw.pop('plot_tail', True)
    
    folder_a,folder_b,ssro_a,ssro_b =  get_data_objects(contains,**kw)
    singleClick = singleClickAnalysis(folder_a,folder_b)

    xlabel,x,y,z_list,u_z_list,tail_list,u_tail_list= get_sweep_analysis_results(singleClick,ssro_a,ssro_b,parameter_name,parameter_range,parameter_kind, **kw)
   
    ylim = (-1,1)#kw.get('ylim', (-0.05, 1.05))
    fig1,ax1 = plt.subplots(1,1)
    
    ax1.set_title(input_data_files[0].timestamp+'\n'+input_data_files[0].measurementstring + '\n sweep: ' + parameter_kind + ' ' + parameter_name)

    if np.shape(z_list)[1] == 1:
        z_list = np.squeeze(z_list)
        u_z_list = np.squeeze(u_z_list)

        ax1.errorbar(y,z_list[:,0],u_z_list[:,0],label = 'psi0 ',fmt = '-o')
        ax1.errorbar(y,z_list[:,1],u_z_list[:,1],label = 'psi1 ',fmt = '-o')
        ax1.set_xlabel(parameter_name)          
    else:
        for i,z,u_z,yparam in zip(range(len(z_list)),z_list,u_z_list,y):

            ax1.errorbar(x,z[:,0],u_z[:,0],label = 'psi0 ' + str(yparam),fmt = '-o',color = (float(i)/len(z_list),0.1,0.1))
            ax1.errorbar(x,z[:,1],u_z[:,1],label = 'psi1 ' + str(yparam),fmt = '--o',color = (float(i)/len(z_list),0.1,0.1))
                    
        ax1.set_xlabel(xlabel)

    ## formatting
    ax1.set_ylabel(r'Correlations')
    ax1.set_ylim(ylim)

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if plot_tail:
        fig3,ax3 = plt.subplots(1,1)

        if np.shape(tail_list)[1] == 1:
            ax3.errorbar(y,tail_list,u_tail_list,fmt = '-o')
            ax3.set_xlabel(parameter_name)
        else:
            for i,t,u_t,yparam in zip(range(len(z_list)),tail_list,u_tail_list,y):
                ax3.errorbar(x,t,u_t,label = str(yparam),fmt = '-o',color = (float(i)/len(z_list),0.1,0.1))
            ax3.set_xlabel(xlabel)
            ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


        ## formatting
        ax3.set_ylabel(r'Tail counts *10^4')
        ax3.set_title(input_data_files[0].timestamp+'\n'+input_data_files[0].measurementstring + '\n' + 'parameter_name: ' + parameter_name)

       

def get_multiple_files(expm_name,**kw):

    base_folder_lt3 = analysis_params.data_settings['base_folder_lt3']
    lt3_folder = os.path.join(base_folder_lt3,expm_name)
    lt3_ssro_folder = os.path.join(base_folder_lt3,'SSROs')
    base_folder_lt4 = analysis_params.data_settings['base_folder_lt4']
    lt4_folder = os.path.join(base_folder_lt4,expm_name)
    lt4_ssro_folder = os.path.join(base_folder_lt4,'SSROs')

    filename_str = kw.pop('filename_str', analysis_params.data_settings['filenames_for_expms'][expm_name])
    print lt3_folder
    b_list=tb.latest_data(contains = filename_str,folder= lt3_folder,return_all = True,**kw)
    a_list=tb.latest_data(contains = filename_str,folder =lt4_folder,return_all = True,**kw)

    if len(b_list) != len(a_list):
        raise(Exception('Different number of files for lt3 and lt4!'))

    ssro_b  = tb.latest_data(contains = 'SSROCalib', folder = lt3_ssro_folder)
    ssro_a  = tb.latest_data(contains = 'SSROCalib',  folder = lt4_ssro_folder)


    return a_list,b_list,ssro_a,ssro_b


def run_multi_file_analysis(expm_name, **kw):
    combine_files = kw.pop('combine_files', True)
    save_corrs = kw.pop('save_corrs', False)

    a_list,b_list,ssro_a,ssro_b =  get_multiple_files(expm_name,**kw)

    sca_list = []
    for i, (folder_a,folder_b) in enumerate(zip(a_list,b_list)):
        sca = singleClickAnalysis(folder_a,folder_b)
        sca.process_correlations(**kw)
        sca_list.append(sca)

        if not(combine_files):
            analyze_spspcorrs(sca,ssro_a,ssro_b,**kw)
    
    if combine_files:
        sweep_pts,p0,p0_u,norm_correlators,norm_correlators_u,counts_per_pt,tail_per_pt, tail_per_pt_u, phi = analyze_spspcorrs(sca_list,ssro_a,ssro_b,ret=True,**kw)
        print p0
        print (1+np.sum(np.abs(p0),axis=1))/4
        if save_corrs:
            base_folder_lt4 = analysis_params.data_settings['base_folder_lt4']
            lt4_folder = os.path.join(base_folder_lt4,expm_name)    
            name = os.path.join(lt4_folder, 'correlations.h5')

            with h5py.File(name, 'w') as hf:
                hf.create_dataset('sweep_pts', data=sweep_pts)
                hf.create_dataset('correlations', data=p0)
                hf.create_dataset('correlations_u', data=p0_u)
                hf.create_dataset('norm_correlators', data=norm_correlators)
                hf.create_dataset('norm_correlators_u', data=norm_correlators_u)
                hf.create_dataset('counts_per_pt', data=counts_per_pt)
                hf.create_dataset('tail_counts', data=tail_per_pt)
                hf.create_dataset('tail_counts_u', data=tail_per_pt_u)
 

def calc_MW_phases(expm_name,single_file=True,save = False,plot_corrs = False,**kw):
    ### Helper function to get the phases from MW phase angle sweeps

    if single_file:
           
        folder_a,folder_b,ssro_a,ssro_b =  get_data_objects(expm_name,**kw)
        singleClick = singleClickAnalysis(folder_a,folder_b)
        singleClick.process_correlations(**kw)
        _,_,_,_,_,_,_,_,phi,phi_u = analyze_spspcorrs(sca,ssro_a,ssro_b,combine_correlation_data = True,flip_psi1 = True,do_sine_fit = True, ret= True, plot_correlations = plot_corrs, **kw)
        
        if save:
            name = os.path.join(folder_a, 'measured_calib_phase.h5')
            with h5py.File(name, 'w') as hf:
                hf.create_dataset('phi', data=phi) 

        return phi, phi_u  
    else:

        a_list,b_list,ssro_a,ssro_b =  get_multiple_files(expm_name,filename_str = 'EntangleXsweepXY',**kw)

        timestamps = []
        phis = []
        for i, (folder_a,folder_b) in enumerate(zip(a_list,b_list)):
            timestamps.append(os.path.split(folder_a)[1][:6])
            print timestamps[-1]
            sca = singleClickAnalysis(folder_a,folder_b)
            sca.process_correlations(**kw)
            _,_,_,_,_,_,_,_,phi = analyze_spspcorrs(sca,ssro_a,ssro_b,combine_correlation_data = True,flip_psi1 = True,do_sine_fit = True, ret= True, plot_correlations = True, **kw)
            phis.append(phi)

        base_folder_lt4 = analysis_params.data_settings['base_folder_lt4']
        lt4_folder = os.path.join(base_folder_lt4,expm_name)
        
        if save:
            name = os.path.join(lt4_folder, 'measured_calib_phases.h5')
            with h5py.File(name, 'w') as hf:
                
                hf.create_dataset('timestamps', data=timestamps)
                hf.create_dataset('phi', data=phis)
           

def extract_pqf_from_sca_list(sca_list):
    # Helper function to grab the pq files from a bunch of sca instances
    return [sca.a.pqf for sca in sca_list]

def plot_ph_hist_and_fltr(sca_list):
    pqf_list = extract_pqf_from_sca_list(sca_list)

    st_start = analysis_params.SPCorr_settings['st_start']
    st_len = analysis_params.SPCorr_settings['st_len']
    ch1_offset = analysis_params.SPCorr_settings['ch1_offset']
    f,(ax0,ax1) = plot_marker_filter_comparison(pqf_list,
                          mrkr_chan = 1,
                          start = st_start-20e3,
                          length= st_len+100e3,
                          hist_binsize = 1e2,save = False,log=True,ret=True)
    ax0.vlines(np.array([st_start,st_start+st_len])/1e3,0,1000,color='r',lw=2)
    ax1.vlines(np.array([st_start+ch1_offset,st_start+st_len+ch1_offset])/1e3,0,1000,color= 'r',lw=2)
    # ax1.set_xlim([(st_start-20e3+ch1_offset)*1e-3,(st_start+st_len+40e3+ch1_offset)*1e-3])
    

def plot_marker_filter_comparison(pqf_list,mrkr_chan = 2,ret=False,**kw):

    is_ph_with_PLU_mrkr = []
    for pqf in pqf_list:
        # get the PLU marked photons first
        is_ph_ch0, is_ph_ch1 = pq_tools.get_photons(pqf)
        is_ph = is_ph_ch0 | is_ph_ch1
        is_ph_with_PLU_mrkr.append(is_ph & pq_tools.filter_marker(pqf, mrkr_chan))

    if ret:
        return pq_plots.plot_photon_hist_filter_comparison(pqf_list,fltr =is_ph_with_PLU_mrkr,**kw)
    else:
        pq_plots.plot_photon_hist_filter_comparison(pqf_list,fltr =is_ph_with_PLU_mrkr,**kw)


