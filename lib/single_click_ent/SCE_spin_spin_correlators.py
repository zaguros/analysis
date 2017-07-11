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

############################# ANALYSIS CLASSES #############################################

class twoSetupMsmt:
    # Skeleton class combining the files from two setups into one analysis

    def __init__(self,folder_a,folder_b):

        self.a = ppq.purifyPQAnalysis(folder_a, hdf5_mode='r') 
        self.b = ppq.purifyPQAnalysis(folder_b, hdf5_mode='r')

class singleClickMsmt(twoSetupMsmt):
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

class SingleClickAnalysis:
    """
    This function takes the results from each run (or can be only one run) and
    analyses them in every which way and does ROC etc.
    """
    def __init__(self,sca_folders,ssro_a,ssro_b,**kw):
        # Pull together the different single click msmts and analyze them

        self.ssro_a = ssro_a
        self.ssro_b = ssro_b
        self.sca_folders = sca_folders

        self.process_single_click_msmts(**kw)
        self.calculate_corrs(**kw)

    def process_single_click_msmts(self,**kw):

        # Start by processing each file individually
        self.singleClickMsmts = []
        for i, (folder_a,folder_b) in enumerate(self.sca_folders):
            scm = singleClickMsmt(folder_a,folder_b)
            scm.process_correlations(**kw)
            self.singleClickMsmts.append(scm)
        
        # Now aggregate the results together
        self.save_folder  = kw.pop('save_folder',self.singleClickMsmts[0].a.folder)
        self.sweep_name = self.singleClickMsmts[0].a.g.attrs['sweep_name']
        self.timestamp = self.singleClickMsmts[0].a.timestamp
        self.measurementstring = self.singleClickMsmts[0].a.measurementstring
        
        self.num_files = len(self.singleClickMsmts)
        self.sweep_pts = self.singleClickMsmts[0].sweep_pts
        self.sweep_length = self.singleClickMsmts[0].sweep_length
        self.gen_sweep_pts = self.singleClickMsmts[0].a.g.attrs['general_sweep_pts']

        self.sweep_2d = len(self.gen_sweep_pts) == 0

        if self.sweep_2d:
            self.gen_sweep_pts1 = self.singleClickMsmts[0].a.g.attrs['general_sweep_pts1']
            self.gen_sweep_pts2 = self.singleClickMsmts[0].a.g.attrs['general_sweep_pts2']

        self.electron_transitions = [self.singleClickMsmts[0].a.g.attrs['electron_transition'], self.singleClickMsmts[0].b.g.attrs['electron_transition']]
        self.E_RO_durations = [self.singleClickMsmts[0].a.g.attrs['E_RO_durations'][0], self.singleClickMsmts[0].b.g.attrs['E_RO_durations'][0]]

        self.default_fig = self.singleClickMsmts[0].a.default_fig # Steal the functions
        self.default_ax = self.singleClickMsmts[0].a.default_ax

        self.correlators_per_sweep_pt = np.zeros([2,self.sweep_length,4])
        self.counts_per_pt = np.zeros([2,self.sweep_length])
        self.tail_per_pt = np.zeros(self.sweep_length)
        self.tail_per_pt_u = np.zeros(self.sweep_length)

        for scm in self.singleClickMsmts:
            
            self.correlators_per_sweep_pt += scm.correlators_per_sweep_pt
            self.counts_per_pt += scm.counts_per_pt
            self.tail_per_pt += scm.tail_per_pt
            self.tail_per_pt_u += scm.tail_per_pt_u**2

        self.tail_per_pt = self.tail_per_pt/self.num_files
        self.tail_per_pt_u = np.sqrt(self.tail_per_pt_u)/self.num_files

    def calculate_corrs(self,**kw):

        ### do ROC
        self.norm_correlators, self.norm_correlators_u = RO_correction_of_correlators(self.correlators_per_sweep_pt,
                                                                        self.electron_transitions,self.E_RO_durations,self.ssro_a,self.ssro_b,**kw)

        self.p0 = []
        self.p0_u = []

        for i,corrs,corrs_u in zip([0,1],self.norm_correlators,self.norm_correlators_u):

            ### extract spin-spin expectation value from correlators
            exp_values,exp_values_u = get_exp_value_from_spin_spin_corr(corrs,corrs_u)

            self.p0.append(exp_values)
            self.p0_u.append(exp_values_u)

        self.p0 = np.array(self.p0)
        self.p0_u = np.array(self.p0_u)

         # Cleverness for if is a 2d sweep
        if self.sweep_2d:

            self.p0 = np.reshape(self.p0, [len(self.p0),len(self.gen_sweep_pts1),len(self.gen_sweep_pts2)])
            self.p0_u = np.reshape(self.p0_u, [len(self.p0),len(self.gen_sweep_pts1),len(self.gen_sweep_pts2)])
            self.norm_correlators = np.reshape(self.norm_correlators, [len(self.p0),len(self.gen_sweep_pts1),len(self.gen_sweep_pts2),4])
            self.norm_correlators_u = np.reshape(self.norm_correlators_u, [len(self.p0),len(self.gen_sweep_pts1),len(self.gen_sweep_pts2),4])
            self.tail_per_pt = np.reshape(self.tail_per_pt, [len(self.gen_sweep_pts1),len(self.gen_sweep_pts2)])
            self.tail_per_pt_u = np.reshape(self.tail_per_pt_u, [len(self.gen_sweep_pts1),len(self.gen_sweep_pts2)])

    def save_corrs(self):

        if self.sweep_2d:
            # Useful to distinguish files so that obvious that a different kind of data
            name = os.path.join(self.save_folder, 'correlations_2d.h5')
           
        else:    
            name = os.path.join(self.save_folder, 'correlations.h5')

        with h5py.File(name, 'w') as hf:
            if self.sweep_2d:
                hf.create_dataset('sweep_pts1', data=self.gen_sweep_pts1)
                hf.create_dataset('sweep_pts2', data=self.gen_sweep_pts2)
            else:
                hf.create_dataset('sweep_pts', data=self.gen_sweep_pts)

            hf.create_dataset('correlations', data=self.p0)
            hf.create_dataset('correlations_u', data=self.p0_u)
            hf.create_dataset('norm_correlators', data=self.norm_correlators)
            hf.create_dataset('norm_correlators_u', data=self.norm_correlators_u)
            hf.create_dataset('counts_per_pt', data=self.counts_per_pt)
            hf.create_dataset('tail_counts', data=self.tail_per_pt)
            hf.create_dataset('tail_counts_u', data=self.tail_per_pt_u)
  
    def plot_temporal_filter(self):
        plot_ph_hist_and_fltr(self.singleClickMsmts)

    def plot_tail(self):
        fig = self.default_fig(figsize=(6,8))
        ax = self.default_ax(fig)

        plt.errorbar(self.sweep_pts, self.tail_per_pt,
                         fmt='o', yerr=self.tail_per_pt_u,markersize=6,capsize=3)
        xlims = plt.xlim()
        newxlims = [xlims[0] - 0.1*(xlims[1] - xlims[0]), xlims[1] + 0.1*(xlims[1] - xlims[0])]
        plt.xlim(newxlims)

    def plot_raw_correlators(self,**kw):

        save_figs            = kw.pop('save_figs',False)
        
        psi_labels = ['Psi0','Psi1']

        for ii,corrs,corrs_u in zip([0,1],self.norm_correlators,self.norm_correlators_u):
            
            # Check if a 2d sweep
            if self.sweep_2d:

                fig = self.default_fig(figsize=(6,8))
                fig.subplots_adjust(hspace=.3)
                for n, pt in enumerate(self.gen_sweep_pts2):
                    ax = fig.add_subplot(len(self.gen_sweep_pts2),1,n+1)
                    if n == 0:
                        plot_title = self.timestamp+'\n'+self.measurementstring + '\n' + psi_labels[ii] +  ', Sweep pt: ' + pt
                    else:
                        plot_title = psi_labels[ii] +  ', Sweep pt: ' + pt
                    ax.set_title(plot_title)

                    ### transpose the functions to be plotted.
                    correlators_trans = np.transpose(corrs[:,n])
                    correlators_u_trans = np.transpose(corrs_u[:,n])
                    labels  = ['11','10','01','00']
                    for e,e_u,l in zip(correlators_trans,correlators_u_trans,labels):
                        ax.errorbar(self.gen_sweep_pts1,e,e_u,fmt='o',label=l)
                    plt.legend()
                    ax.set_xlabel(self.sweep_name[0])
                    ax.set_ylabel('Probability')
                    ax.set_ylim([0,1])
                    xlims = plt.xlim()
                    newxlims = [xlims[0] - 0.1*(xlims[1] - xlims[0]), xlims[1] + 0.1*(xlims[1] - xlims[0])]
                    plt.xlim(newxlims)
                plt.show()
            else:

                ### transpose the functions to be plotted.
                correlators_trans = np.transpose(corrs)
                correlators_u_trans= np.transpose(corrs_u)

                labels  = ['11','10','01','00']
                fig = plt.figure()
                ax = plt.subplot()
                for e,e_u,l in zip(correlators_trans,correlators_u_trans,labels):
                    ax.errorbar(self.sweep_pts,e,e_u,fmt='o',label=l)
                plt.legend()
                ax.set_xlabel(self.sweep_name)
                ax.set_ylabel('Probability')
                ax.set_ylim([0,1])
                ax.set_title(self.timestamp+'\n'+self.measurementstring+ '\n' + 'psi'+str(ii))
                xlims = plt.xlim()
                newxlims = [xlims[0] - 0.1*(xlims[1] - xlims[0]), xlims[1] + 0.1*(xlims[1] - xlims[0])]
                plt.xlim(newxlims)
                plt.show()
            if save_figs:
                fig.savefig(
                        os.path.join(self.save_folder, 'correlators_vs_sweepparam_psi_' +str(ii) +'.pdf'),format='pdf')

    def plot_correlations(self,**kw): 
        print_fids           = kw.pop('print_fids', False)
        do_sine_fit          = kw.pop('do_sine_fit',False)
        combine_correlation_data = kw.pop('combine_correlation_data', False)
        flip_psi1            = kw.pop('flip_psi1',False)
        abs_corrs            = kw.pop('abs_corrs', False)
        save_figs            = kw.pop('save_figs',False)

        plot_p0 = self.p0
        plot_p0_u = self.p0_u

        if flip_psi1 :
            plot_p0[1] = -plot_p0[1]

        if abs_corrs:
            plot_p0 = np.abs(plot_p0)
            lims =  [0.,1.05]
        else:
            lims =  [-1.05,1.05]

        if combine_correlation_data:
            plot_p0 = np.array([(plot_p0[0] + plot_p0[1])/2.0]) # Need to use with flip_psi1 usually 
            plot_p0_u = np.array([np.sqrt(plot_p0_u[0]**2 + plot_p0_u[1]**2)/2.0])
            labels = ['combined data']
        else:
            labels = ['psi0','psi1']
        
        phis = []
        phis_u = []

        # Check if a 2d sweep
        if self.sweep_2d:
            fig = self.default_fig(figsize=(6,8))
            ax = []
            for n, pt in enumerate(self.gen_sweep_pts2):
                ax.append(fig.add_subplot(len(self.gen_sweep_pts2),1,n+1))
                if n == 0:
                    plot_title = self.timestamp+'\n'+self.measurementstring + '\n' + 'Sweep pt: ' + pt
                else:
                    plot_title = 'Sweep pt: ' + pt
                ax[n].set_title(plot_title)

            fig.subplots_adjust(hspace=.3)

            for jj,(p,p_u) in enumerate(zip(plot_p0,plot_p0_u)):
                
                x = self.gen_sweep_pts1
                
                if print_fids: print 'Fidelity ', (1+np.sum(np.abs(p),axis=1))/4,  np.sqrt(np.sum(p_u**2,axis=1))/4

                # Iterate over second set of sweep pts
                for z, (lab, pp, pu) in enumerate(zip(self.gen_sweep_pts2, p.T, p_u.T)):
                    plt.sca(ax[z])
                    plt.errorbar(x, pp,
                         fmt='o', yerr=pu,markersize=6,capsize=3)

                    set_corr_plot_properties(ax[z],lims)
                    ax[z].set_xlabel(self.sweep_name[0])
                    
                    
                    if do_sine_fit:    
                       phi,phi_u = fit_sin_and_calc_phi(x,pp,ax = ax[z])
                       phis.append(phi)
                       phis_u.append(phi_u)
        else:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)

            x = self.sweep_pts
            for jj,(p,p_u) in enumerate(zip(plot_p0,plot_p0_u)):
                ax.errorbar(x, p,
                         fmt='o', yerr=p_u,markersize=6,capsize=3,label= labels[jj])
                set_corr_plot_properties(ax,lims)
                ax.set_xlabel(self.sweep_name)
                plt.legend()

                if print_fids: print 'Fidelity ', (1+np.sum(np.abs(plot_p0)))/4,  np.sqrt(np.sum(plot_p0_u**2))/4
                

                if do_sine_fit:    
                    phi,phi_u = fit_sin_and_calc_phi(x,p,ax = ax)
                    phis.append(phi)
                    phis_u.append(phi_u)

        if do_sine_fit:
            
            if np.size(phi) > 1:
                phi = np.mean(phis)
                phi_u = quadrature_sum(phis_u)
                print 'Avg. phi angle ', phi, phi_u
                
        if save_figs:
            fig.savefig(
                os.path.join(self.save_folder, 'correlations_vs_sweepparam.pdf'),
                format='pdf')

############################# HELPER FUNCTIONS #############################################

def set_corr_plot_properties(ax,ylims):

    ax.axhspan(0,1,fill=False,ls='dotted')
    ax.axhspan(-1,0,fill=False,ls='dotted')
    plt.ylim(ylims)
    xlims = plt.xlim()
    newxlims = [xlims[0] - 0.02*(xlims[1] - xlims[0]), xlims[1] + 0.02*(xlims[1] - xlims[0])]
    plt.xlim(newxlims) 
    ax.set_ylabel('Correlations')
    

def fit_sin_and_calc_phi(x,y,ax = False):

    g_a = 0.0
    g_A = np.amax(y)
    g_phi = x[np.argmax(y)] 
    ### frequency guess is hardcoded as this is mainly intended for specific entangling oscillations.
    g_f = 1/360.
    p0s, fitfunc,fitfunc_str = common.fit_cos(g_f,g_a,g_A,g_phi)

    fit_result = fit.fit1d(x,y, None, p0=p0s, fitfunc=fitfunc,
         ret=True,fixed=[0,1])

    if ax:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax, 
            plot_data=False,print_info = False)
        A = fit_result['params_dict']['A']
        phi = fit_result['params_dict']['phi']
        if A < 0:
            phi = phi + 180
            A = -A
        phi = np.mod(-phi,360)
        phi_u = fit_result['error_dict']['phi']
        print 'A,phi,phi_u ', A, phi, phi_u
    

    return phi, phi_u  # Phi and Phi_u

def quadrature_sum(errors):
    return np.sqrt(np.sum(np.array(errors)**2))/len(errors)

def RO_correction_of_correlators(correlators_per_sweep_pt,electron_transitions,E_RO_durations,ssro_a,ssro_b,**kw):

    verbose = kw.pop('verbose',False)
    do_ROC = kw.pop('do_ROC',True)

    ### to store the estimated statistical uncertainty
    norm_correlators_u = np.zeros(np.shape(correlators_per_sweep_pt))
    norm_correlators = np.zeros(np.shape(correlators_per_sweep_pt))

    if do_ROC:
        ### get ssro_ROC for LT3 --> corresponds to setup B
        F0_LT3,F1_LT3 = get_RO_fidelities(ssro_b,electron_transitions[1],E_RO_durations[1])
        ### get ssro_ROC for LT4 --> corresponds to setup A
        F0_LT4,F1_LT4 = get_RO_fidelities(ssro_a,electron_transitions[0],E_RO_durations[0])

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

def get_RO_fidelities(ssro_folder,mw_transition,E_RO_duration):

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


def get_data_objects(contains,**kw):
    # Helper function to grab the desired file objects

    use_file_library = kw.pop('use_file_library', False)
    analysis_computer = kw.pop('analysis_computer', 'lt4')

    if not(use_file_library):

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
            a_list =tb.latest_data(contains_lt4,**kw)
            b_list =tb.latest_data(contains_lt3,folder =r'Z:\data',**kw)
            ssro_b  = tb.latest_data(contains_lt3_ssro, folder =r'Z:\data')
            ssro_a  = tb.latest_data(contains_lt4_ssro)

        elif analysis_computer == 'lt3_analysis':
            a_list =tb.latest_data(contains_lt4,folder= r'X:\data',**kw)
            b_list =tb.latest_data(contains_lt3,folder =r'Y:\data',**kw)
            ssro_b  = tb.latest_data(contains_lt3_ssro, folder =r'Y:\data')
            ssro_a  = tb.latest_data(contains_lt4_ssro,  folder =r'X:\data')

    else:

        # We have a library of final data on the lt3_analysis computer, the details of which are stored in Analysis params sce.

        base_folder_lt3 = analysis_params.data_settings['base_folder_lt3']
        lt3_folder = os.path.join(base_folder_lt3,contains)
        lt3_ssro_folder = os.path.join(base_folder_lt3,'SSROs')
        base_folder_lt4 = analysis_params.data_settings['base_folder_lt4']
        lt4_folder = os.path.join(base_folder_lt4,contains)
        lt4_ssro_folder = os.path.join(base_folder_lt4,'SSROs')
        filename_str = kw.pop('filename_str', analysis_params.data_settings['filenames_for_expms'][contains])

        b_list=tb.latest_data(contains = filename_str,folder= lt3_folder,return_all = True,**kw)
        a_list=tb.latest_data(contains = filename_str,folder =lt4_folder,return_all = True,**kw)
        ssro_b  = tb.latest_data(contains = 'SSROCalib', folder = lt3_ssro_folder)
        ssro_a  = tb.latest_data(contains = 'SSROCalib',  folder = lt4_ssro_folder)

    a_list = a_list if isinstance(a_list, list) else [a_list] # Can use kw flag to return multiple files from tb.latest_data even if not using library, which is useful.
    b_list = b_list if isinstance(b_list, list) else [b_list]
    
    if len(b_list) != len(a_list):
        raise(Exception('Different number of files for lt3 and lt4!'))

    sca_folders = zip(a_list,b_list)
    return sca_folders,ssro_a,ssro_b


def get_sweep_analysis_results(sca,parameter_name,parameter_range,parameter_kind,**kw):
    
    min_or_max = kw.pop('min_or_max','min')

    y = parameter_range
    z = []
    u_z = []
    tail = []
    u_tail = []

    x = sca.sweep_pts
    xlabel = sca.sweep_name
    
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
        
        sca.process_single_click_msmts(**kw) # Have to reprocess with new filter values
        sca.calculate_corrs(**kw)

        ### store sweep results
        z.append(sca.p0)
        u_z.append(sca.p0_u)
        tail.append(sca.tail_cts)
        u_tail.append(sca.tail_cts_u)

    #restore analysis parameters
    reload(analysis_params)
    return xlabel,x,y,z,u_z,tail,u_tail

def sweep_analysis_parameter(contains,parameter_name, parameter_range, **kw):
    # Note that not tested since extensive rewrite

    parameter_kind = kw.pop('parameter_kind','SPCorr')
    plot_tail = kw.pop('plot_tail', True)
    
    sca_folders,ssro_a,ssro_b =  get_data_objects(contains,**kw)
    sca = SingleClickAnalysis(sca_folders,ssro_a,ssro_b)
    xlabel,x,y,z_list,u_z_list,tail_list,u_tail_list= get_sweep_analysis_results(sca,parameter_name,parameter_range,parameter_kind, **kw)
   
    ylim = (-1,1)
    fig1,ax1 = plt.subplots(1,1)
    
    ax1.set_title(sca.timestamp+'\n'+sca.measurementstring + '\n sweep: ' + parameter_kind + ' ' + parameter_name)

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

      

############################# HIGH LEVEL ANALYSIS FUNCTIONS #############################################

def run_analysis(contains, **kw):
    # This function can do all the things.

    use_file_library = kw.get('use_file_library', False)
    save_corrs = kw.pop('save_corrs', False)

    plot_temporal_filter = kw.pop('plot_temporal_filter',False)
    plot_tail = kw.pop('plot_tail',False)
    plot_raw_correlators = kw.pop('plot_raw_correlators',False)
    plot_correlations = kw.pop('plot_correlations',False)

    sca_folders,ssro_a,ssro_b =  get_data_objects(contains,**kw)

    if use_file_library:    
            base_folder_lt4 = analysis_params.data_settings['base_folder_lt4']
            save_folder = os.path.join(base_folder_lt4,contains)
    else:
        save_folder = False

    sca = SingleClickAnalysis(sca_folders,ssro_a,ssro_b,save_folder = save_folder,**kw)

    if plot_temporal_filter: sca.plot_temporal_filter()
    if plot_tail: sca.plot_tail()
    if plot_raw_correlators: sca.plot_raw_correlators(**kw)
    if plot_correlations: sca.plot_correlations(**kw)
    
    if save_corrs:
        sca.save_corrs() 
        
def check_phase_calibration(**kw):
    run_analysis(kw.pop('contains', 'XsweepY'),plot_correlations =True,do_ROC = True,do_sine_fit = True,flip_psi1=True,combine_correlation_data=True, **kw)

def calc_MW_phases(expm_name,single_file=True,save = False,plot_corrs = False,**kw):
    ### Helper function to get the phases from MW phase angle sweeps NOTE NOT TESTED SINCE REWRITE

    if single_file:
           
        folder_a,folder_b,ssro_a,ssro_b =  get_data_objects(expm_name,**kw)
        singleClick = singleClickMsmt(folder_a,folder_b)
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
            sca = singleClickMsmt(folder_a,folder_b)
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
           
