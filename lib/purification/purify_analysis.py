"""
Provides analysis functions for purification measurements

NK 2016
"""


import purify_pq as ppq
import numpy as np
import os,h5py
from analysis.lib.tools import toolbox as tb; reload(tb)
from analysis.lib.pq import pq_tools
import copy as cp
from matplotlib import pyplot as plt
from analysis.lib.m2.ssro import ssro
from analysis.lib.m2 import m2
from analysis.lib.lde import sscorr; reload(sscorr)
import purify_analysis_params as analysis_params;reload(analysis_params)
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

#### standard parameters

save_plot_folder = r'D:\measuring\data\purification_plots'


class purify_analysis(object):
    """
    general class that stores prefiltered data and serves as analysis suite for non-local correlation measurements
    """

    def __init__(self,name,lt3_folder,lt4_folder,ROC_lt3_tstamp,ROC_lt4_tstamp,**kw):
        """
        data is in general stored in dictionaries associated with one of the two setups. each dictionary entry contains a list of np.arrays (each entry corresponding to one given timestamp)
        we additionally store the full purify_pq object for further processing of the raw data
        """
        self.name = name

        self.key_list_pq =  ['/PQ_sync_number-1','/PQ_channel-1','/PQ_sync_time-1','/PQ_special-1','/PQ_time-1']
        self.key_lst_adwin_data = ['ssro_results','electron_readout_result','Phase_correction_repetitions','CR_after','CR_before','carbon_readout_result','attempts_first','attempts_second']
        self.key_list_misc_to_save = ['counted_awg_reps_w1','counted_awg_reps_w2','total_elapsed_time','total_syncs']
        self.key_list_misc = ['tstamp','data_attrs','joint_attrs']

        self.lt4_dict = {}
        self.lt3_dict = {}

        ### initialize the data dictionaries
        for key in self.key_list_pq:
            self.lt3_dict.update({key:[]})
            self.lt4_dict.update({key:[]})

        for key in self.key_lst_adwin_data:
            self.lt3_dict.update({key:[]})
            self.lt4_dict.update({key:[]})
        
        for key in self.key_list_misc:
            self.lt3_dict.update({key:[]})
            self.lt4_dict.update({key:[]})

        for key in self.key_list_misc_to_save:
            self.lt3_dict.update({key:[]})
            self.lt4_dict.update({key:[]})

        self.lt3_folder = lt3_folder
        self.lt4_folder = lt4_folder

        self.ROC_lt3_tstamp = ROC_lt3_tstamp
        self.ROC_lt4_tstamp = ROC_lt4_tstamp

    def load_raw_data(self,lt3_timestamps,lt4_timestamps,force_evaluation = False,save=True,verbose = False):
        """
        this script takes a list of timestamps for both setups and prefilters them according to adwin filters
        creates a list of arrays associated with each time stamp for adwin_ssro,syncs, time, special, channel, marker and attaches is it to the data object for further processing.
        """
        length = len(lt3_timestamps)
        print 'loading the data, total number of files ', length

        ### initialize the data dictionaries
        for key in self.key_list_pq:
            self.lt3_dict.update({key:[]})
            self.lt4_dict.update({key:[]})

        for key in self.key_lst_adwin_data:
            self.lt3_dict.update({key:[]})
            self.lt4_dict.update({key:[]})
        
        for key in self.key_list_misc:
            self.lt3_dict.update({key:[]})
            self.lt4_dict.update({key:[]})

        for key in self.key_list_misc_to_save:
            self.lt3_dict.update({key:[]})
            self.lt4_dict.update({key:[]})

        i = 0

        for t_lt3,t_lt4 in zip(lt3_timestamps,lt4_timestamps):
            
            lt3_folder = tb.latest_data(t_lt3,folder = self.lt3_folder)
            lt4_folder = tb.latest_data(t_lt4,folder = self.lt4_folder)

            cleaned_data_lt4 = os.path.join(lt4_folder,'cleaned_data.hdf5')

            a_lt3 = ppq.purifyPQAnalysis(lt3_folder,hdf5_mode ='r')
            a_lt4 = ppq.purifyPQAnalysis(lt4_folder,hdf5_mode ='r')

            if not(force_evaluation) and check_if_file_exists(cleaned_data_lt4):

                for key in self.key_list_pq:
                    vals, attrs = tb.get_analysis_data(cleaned_data_lt4,str(key[1:])) # Need to drop first slash!
                    self.lt4_dict[key].append(vals)   

                for key in self.key_list_misc_to_save:
                    vals, attrs = tb.get_analysis_data(cleaned_data_lt4,str(key))
                    self.lt4_dict[key].append(vals)

            else:
    
                if save:
                    try:
                        f = h5py.File(cleaned_data_lt4, 'r')
                    except:
                        f = h5py.File(cleaned_data_lt4,'w-')
                    f.close() 

                ### filter the timeharp data according to adwin events / syncs

                ### need to create two sync filters here.!
                sync_filter = a_lt4.filter_pq_data_from_adwin_syncs() ### this syncfilter erases all data where from the PQ data where the adwin did NOT read out

                if len(sync_filter) == 0: # empty list --> no read outs.
                    # print 'file empty, skipping these time stamps:',t_lt3,t_lt4
                    # print
                    continue

                ### store relevant adwin results

                if -1 in a_lt3.agrp['attempts_second'].value: ### this indicates that the measurement ran without LDE2
                    syncs_w1 = (np.array(a_lt3.agrp['counted_awg_reps'].value))
                    print 'Ran without LDE2!'
                else:
                    syncs_w1 = (np.array(a_lt3.agrp['counted_awg_reps'].value)-np.array(a_lt3.agrp['attempts_second'].value))
                syncs_w2 = np.array(a_lt3.agrp['counted_awg_reps'].value)


                if len(syncs_w1) == 0:
                    print 'syncs_w1 empty, skipping these time stamps:',t_lt3,t_lt4
                    print
                    continue

                self.lt4_dict['counted_awg_reps_w2'].append(syncs_w2)
                self.lt4_dict['counted_awg_reps_w1'].append(syncs_w1)

                ### need to create two sync filters here. This is done by appending them to each other and sorting the relatively short array. (should never be more than 500 events??)
                ### this syncfilter erases all data from the PQ data where the adwin did NOT read out. Drastically reduces the data amount we have to handle.
                sync_filter = a_lt4.filter_pq_data_from_adwin_syncs(np.sort(np.append(syncs_w1,syncs_w2)))
                if len(sync_filter) == 0: # empty list --> no read outs.
                    print 'file empty, skipping these time stamps:',t_lt3,t_lt4
                    print
                    continue

                self.lt4_dict['total_elapsed_time'].append((a_lt4.pqf['/PQ_time-1'].value[-1]-a_lt4.pqf['/PQ_time-1'].value[0])*1e-12)
                self.lt4_dict['total_syncs'].append(a_lt4.pqf['/PQ_sync_number-1'].value[-1])

                for key in self.key_list_pq:
                    vals = np.array(a_lt4.pqf[key].value[sync_filter])
                    self.lt4_dict[key].append(vals)

                    if save:
                        tb.set_analysis_data(cleaned_data_lt4,unicode(key[1:]), vals, []) # Need to drop first slash!

                for key in self.key_list_misc_to_save:
                    if save:
                        tb.set_analysis_data(cleaned_data_lt4,unicode(key), self.lt4_dict[key][i], [])

            self.lt3_dict['tstamp'].append(t_lt3)
            self.lt4_dict['tstamp'].append(t_lt4)
            self.lt3_dict['data_attrs'].append(convert_attrs_to_dict(a_lt3.g.attrs.items()))
            self.lt4_dict['data_attrs'].append(convert_attrs_to_dict(a_lt4.g.attrs.items()))
            self.lt3_dict['joint_attrs'].append(convert_attrs_to_dict(a_lt3.joint_grp.attrs.items()))
            self.lt4_dict['joint_attrs'].append(convert_attrs_to_dict(a_lt4.joint_grp.attrs.items()))
            
            for key in self.key_lst_adwin_data:
                self.lt4_dict[key].append(np.array(a_lt4.agrp[key].value))
            for key in self.key_lst_adwin_data:
                self.lt3_dict[key].append(np.array(a_lt3.agrp[key].value))

            a_lt3.finish()
            a_lt4.finish()

            #### calculate the duty cycle for that specific file.
            # print 'lde length',a_lt3.joint_grp.attrs['LDE_element_length']
            # print'first and last time',a_lt4.pqf['/PQ_time-1'].value[0],a_lt4.pqf['/PQ_time-1'][-1]
            # print 'last elapsed time in sequence vs total elapsed time',  

            total_elapsed_time = self.lt4_dict['total_elapsed_time'][i]

            time_in_LDE_sequence = self.lt3_dict['joint_attrs'][i]['LDE_element_length']*self.lt4_dict['total_syncs'][i]
           
            if verbose:
                print 'file no ', i+1 , ' with duty cycle of', round(100*time_in_LDE_sequence/total_elapsed_time,1), ' %'
            i += 1


    def correct_pq_times(self,offsets = [],offsets_ch1 = []): # Some of the data has different timings, since we changed one of the ZPL apds. Here we can fix this!

        if np.size(offsets):
            for i,offset in enumerate(offsets):
                self.lt4_dict['/PQ_sync_time-1'][i] = self.lt4_dict['/PQ_sync_time-1'][i] + offset

        if np.size(offsets_ch1):
            for i,offset_ch1 in enumerate(offsets_ch1):
                ch1_fltr = self.lt4_dict['/PQ_channel-1'][i]==1
                self.lt4_dict['/PQ_sync_time-1'][i][ch1_fltr] = self.lt4_dict['/PQ_sync_time-1'][i][ch1_fltr] + offset_ch1   
   

    def check_tail_w1_w2(self,st_start = 2000e3,st_len = 50e3):
        """
        goes through the raw_data and selects only files according to the applied sync filter: returns the measured tail in each window
        for the purification experiment there will only be one window
        """

        i = 0
        tails_w1 = []

        for t_lt4 in self.lt4_dict['tstamp']:
                  
            a_lt4 = ppq.purifyPQAnalysis(tb.latest_data(t_lt4,folder = self.lt4_folder),hdf5_mode ='r')

            ### analysis for channel 0  && window 1
            w1_ch0 = self.get_total_number_of_clicks_in_window(a_lt4,0,st_start,st_len)
            w1_ch1 = self.get_total_number_of_clicks_in_window(a_lt4,1,st_start,st_len)
            last_sync = a_lt4.pqf['/PQ_sync_number-1'][-1]

            tail_w1 = round(1e4*(w1_ch0+w1_ch1)/last_sync,2)

            # print 'tail in w1 / w2 (1e-4)    ', tail_w1, ' / ', tail_w2
            tails_w1.append(tail_w1);

        f,ax = self.create_plot(ylabel = 'Tail (10e-4)', title = 'Tail counts vs run number; w1_start ' +  str(round(st_start/1e3,0)))
        self.plot_data(range(len(tails_w1)),tails_w1,label = 'w1')
        plt.legend()
        self.save_and_close_plot()
        # return np.array(tails_w1),np.array(tails_w2)


    def get_total_number_of_clicks_in_window(self,a,channel,st_start,st_len):
            

        is_ph = pq_tools.get_photons(a.pqf)[channel]
        bins = np.arange(st_start-.5,st_start+st_len,1e3)
        y,x=np.histogram(a.pqf['/PQ_sync_time-1'].value[np.where(is_ph)], bins=bins)
        x=x[:-1]
        # print 'Total clicks:', np.sum(y)

        return np.sum(y)

    def apply_temporal_filters_to_prefiltered_data(self,st_start = None,st_len = None,st_len_w2 = None,verbose = True):
        """
        applies temporal filters to all the registered photon detections in the relevant window
        """
        #####
        self.st_fltr_w1 = []
        self.st_fltr_w2 = []

        ### one can also apply manual filters if one wants to deviate from the prescribed parameter dictionary
        if st_start == None:
            st_start = analysis_params.filter_settings['st_start']
        if st_len == None:
            st_len = analysis_params.filter_settings['st_len']
        if st_len_w2 == None:
            st_len_w2 = analysis_params.filter_settings['st_len_w2']


        for st_filtered,sp_filtered in zip(self.lt4_dict['/PQ_sync_time-1'],self.lt4_dict['/PQ_special-1']):
            
            st_fltr_w1 = (st_filtered > st_start)  & (st_filtered < (st_start  + st_len)) & (sp_filtered == 0)
            st_fltr_w2 = (st_filtered > st_start)  & (st_filtered < (st_start + st_len_w2)) & (sp_filtered == 0)
            self.st_fltr_w1.append(st_fltr_w1)
            self.st_fltr_w2.append(st_fltr_w2)

            no_w1 = np.sum(st_fltr_w1)
            if verbose:
                print 'number of total filtered detection events : ', no_w1
        return


    def apply_sync_filter_w1_w2(self,verbose = True,max_w2 = None):
        """
        checks if associated sync number corresponds to a click in entanglement generation 1 or 2 and alters self.st_fltr_w1 and self.st_fltr_w2
        is applied after temporal filtering.
        Also has the ability to filter for decoherence of the nuclear spin state based on min and max attempts for LDE 1 and 2.
        """
        
        temp_fltr_w1 = cp.deepcopy(self.st_fltr_w1) ## store
        temp_fltr_w2 = cp.deepcopy(self.st_fltr_w2) ## store

        self.st_fltr_w1, self.st_fltr_w2 = [],[] ### reinitialize


        ### get deocherence filter
        max_w1 = analysis_params.filter_settings['max_reps_w1']
        min_w1 = analysis_params.filter_settings['min_reps_w1']
        if max_w2 == None:
            max_w2 = analysis_params.filter_settings['max_reps_w2']
        min_w2 = analysis_params.filter_settings['min_reps_w2']


        loop_array = zip(temp_fltr_w1,temp_fltr_w2,self.lt4_dict['/PQ_sync_number-1'],self.lt4_dict['counted_awg_reps_w1'],self.lt4_dict['counted_awg_reps_w2'],self.lt4_dict['attempts_first'],self.lt4_dict['attempts_second'])
        
        no_w1, no_w2 = 0, 0

        for fltr_w1,fltr_w2,sync_nrs,adwin_nrs_w1,adwin_nrs_w2,attempts_first,attempts_second in loop_array:

            decoherence_filter_w1 = (attempts_first >= min_w1) & (attempts_first <= max_w1)
            decoherence_filter_w2 = (attempts_second >= min_w2) & (attempts_second <= max_w2)


            dec_and_w1_filter = np.zeros(len(sync_nrs))
            # Start by finding when get right sync number       
            w1_sync_filtered = np.in1d(sync_nrs,adwin_nrs_w1)
            # Gotta pull the right element from the decoherence filter list. Messy
            dec_and_w1_filter[w1_sync_filtered] = decoherence_filter_w1[self.filter_adwin_data_from_pq_syncs(sync_nrs[w1_sync_filtered],adwin_nrs_w1)]


            dec_and_w2_filter = np.zeros(len(sync_nrs))
            # Start by finding when get right sync number       
            w2_sync_filtered = np.in1d(sync_nrs,adwin_nrs_w2)
            # Gotta pull the right element from the decoherence filter list. Messy
            dec_and_w2_filter[w2_sync_filtered] = decoherence_filter_w2[self.filter_adwin_data_from_pq_syncs(sync_nrs[w2_sync_filtered],adwin_nrs_w2)]

            # Combine with original filter
            st_fltr_w1 = np.logical_and(fltr_w1, dec_and_w1_filter)
            st_fltr_w2 = np.logical_and(fltr_w2, dec_and_w2_filter)

            self.st_fltr_w1.append(st_fltr_w1)
            self.st_fltr_w2.append(st_fltr_w2)

            no_w1 += np.sum(st_fltr_w1)
            no_w2 += np.sum(st_fltr_w2)
        
        if verbose:
            print 'number of filtered detection events in each window w1 / w2: ', no_w1, ' / ', no_w2


    def apply_is_purified_filter(self,signature = '11',verbose = True):
        """
        correlates the electron RO signature after the purification gate to "ms0 & ms0"
        Returns a filter for the adwin_ssro results.
        Input: the desired purification signature
        Output: none. it only modifies the attributes self.st_fltr_w1 and self.st_fltr_w2
        """

        temp_fltr_w1 = cp.deepcopy(self.st_fltr_w1) ## store
        temp_fltr_w2 = cp.deepcopy(self.st_fltr_w2) ## store

        self.st_fltr_w1, self.st_fltr_w2 = [],[] ### reinitialize

        loop_arrays = zip (temp_fltr_w1,temp_fltr_w2,self.lt4_dict['/PQ_sync_number-1'],self.lt4_dict['counted_awg_reps_w1'],self.lt4_dict['counted_awg_reps_w2'],self.lt3_dict['electron_readout_result'],self.lt4_dict['electron_readout_result'])
        
        no_w1, no_w2 = 0, 0
        counts = np.zeros(4)

        for fltr_w1,fltr_w2,sync_nrs,adwin_nrs_w1,adwin_nrs_w2,e_ro_lt3,e_ro_lt4 in loop_arrays:

            adwin_indices_w1  = self.filter_adwin_data_from_pq_syncs(sync_nrs[fltr_w1],adwin_nrs_w1)
            adwin_indices_w2  = self.filter_adwin_data_from_pq_syncs(sync_nrs[fltr_w2],adwin_nrs_w2)

            if verbose:
                sigs = ['00','01','10','11']
                
                for z, sig in enumerate(sigs):

                    purified_filter_w1 = np.core.defchararray.add((e_ro_lt3[adwin_indices_w1]).astype('str'),(e_ro_lt4[adwin_indices_w1]).astype('str')) == sig 
                    counts[z] = np.sum(purified_filter_w1)

            purified_filter_w1 = np.core.defchararray.add((e_ro_lt3[adwin_indices_w1]).astype('str'),(e_ro_lt4[adwin_indices_w1]).astype('str')) == signature 
            fltr_w1[fltr_w1] = purified_filter_w1

            purified_filter_w2 = np.core.defchararray.add((e_ro_lt3[adwin_indices_w2]).astype('str'),(e_ro_lt4[adwin_indices_w2]).astype('str')) == signature 
            fltr_w2[fltr_w2] = purified_filter_w2


            self.st_fltr_w1.append(fltr_w1)
            self.st_fltr_w2.append(fltr_w2)
            
            no_w1 += np.sum(fltr_w1)
            no_w2 += np.sum(fltr_w2)
        
        if verbose:

            probs = counts/float(np.sum(counts))

            row_format ="     {:>.2f} " * (len(sigs))
            headline_format = "{:>9} " * len(sigs)
            print '\nProbability of each purification outcome:'
            print headline_format.format(*sigs)

            print "-"*(50)
            print row_format.format(*probs)
            row_format =" "+"    ({:>d}) " * (len(sigs))
            print row_format.format(*(counts).astype('int'))

            print 'number of filtered detection events in each window w1 / w2: ', no_w1, ' / ', no_w2
            print

    def apply_CR_after_filter(self,min_cr_lt3 = None, min_cr_lt4 = None, verbose = True):
        '''
        Checks self.st_fltr_w1 and self.st_fltr_w2 if the CR check after the event was below a certain treshold.
        '''

        ### one can also apply manual filters if one wants to deviate from the prescribed parameter dictionary
        if min_cr_lt3 == None:
            min_cr_lt3 = analysis_params.filter_settings['min_cr_lt3_after']
        if min_cr_lt4 == None:
            min_cr_lt4 = analysis_params.filter_settings['min_cr_lt4_after']


        temp_fltr_w1 = cp.deepcopy(self.st_fltr_w1) ## store
        temp_fltr_w2 = cp.deepcopy(self.st_fltr_w2) ## store

        self.st_fltr_w1, self.st_fltr_w2 = [],[] ### reinitialize

        loop_arrays = zip (temp_fltr_w1,temp_fltr_w2,self.lt4_dict['/PQ_sync_number-1'],self.lt4_dict['counted_awg_reps_w1'],self.lt4_dict['counted_awg_reps_w2'],self.lt3_dict['CR_after'],self.lt4_dict['CR_after'])
        
        no_w1, no_w2 = 0, 0

        for fltr_w1,fltr_w2,sync_nrs,adwin_nrs_w1,adwin_nrs_w2,cr_lt3,cr_lt4 in loop_arrays:

            adwin_indices_w1  = self.filter_adwin_data_from_pq_syncs(sync_nrs[fltr_w1],adwin_nrs_w1)
            adwin_indices_w2  = self.filter_adwin_data_from_pq_syncs(sync_nrs[fltr_w2],adwin_nrs_w2)
            
            cr_filter_w1 = np.logical_and(cr_lt3[adwin_indices_w1] > min_cr_lt3, cr_lt4[adwin_indices_w1] > min_cr_lt3)
            cr_filter_w2 = np.logical_and(cr_lt3[adwin_indices_w2] > min_cr_lt3, cr_lt4[adwin_indices_w2] > min_cr_lt4)
            
            fltr_w1[fltr_w1] = cr_filter_w1
            fltr_w2[fltr_w2] = cr_filter_w2

            self.st_fltr_w1.append(fltr_w1)
            self.st_fltr_w2.append(fltr_w2)
            
            no_w1 += np.sum(fltr_w1)
            no_w2 += np.sum(fltr_w2)
        
        if verbose:

            print 'number of filtered detection events in each window w1 / w2: ', no_w1, ' / ', no_w2
            print

    def apply_CR_before_filter(self,min_cr_lt3 = None, min_cr_lt4 = None, verbose = True):
        '''
        Checks self.st_fltr_w1 and self.st_fltr_w2 if the CR check after the event was below a certain treshold.
        '''

        ### one can also apply manual filters if one wants to deviate from the prescribed parameter dictionary
        if min_cr_lt3 == None:
            min_cr_lt3 = analysis_params.filter_settings['min_cr_lt3_before']
        if min_cr_lt4 == None:
            min_cr_lt4 = analysis_params.filter_settings['min_cr_lt4_before']


        temp_fltr_w1 = cp.deepcopy(self.st_fltr_w1) ## store
        temp_fltr_w2 = cp.deepcopy(self.st_fltr_w2) ## store

        self.st_fltr_w1, self.st_fltr_w2 = [],[] ### reinitialize

        loop_arrays = zip (temp_fltr_w1,temp_fltr_w2,self.lt4_dict['/PQ_sync_number-1'],self.lt4_dict['counted_awg_reps_w1'],self.lt4_dict['counted_awg_reps_w2'],self.lt3_dict['CR_before'],self.lt4_dict['CR_before'])
        
        no_w1, no_w2 = 0, 0

        for fltr_w1,fltr_w2,sync_nrs,adwin_nrs_w1,adwin_nrs_w2,cr_lt3,cr_lt4 in loop_arrays:

            adwin_indices_w1  = self.filter_adwin_data_from_pq_syncs(sync_nrs[fltr_w1],adwin_nrs_w1)
            adwin_indices_w2  = self.filter_adwin_data_from_pq_syncs(sync_nrs[fltr_w2],adwin_nrs_w2)
            
            cr_filter_w1 = np.logical_and(cr_lt3[adwin_indices_w1] > min_cr_lt3, cr_lt4[adwin_indices_w1] > min_cr_lt3)
            cr_filter_w2 = np.logical_and(cr_lt3[adwin_indices_w2] > min_cr_lt3, cr_lt4[adwin_indices_w2] > min_cr_lt4)
            
            fltr_w1[fltr_w1] = cr_filter_w1
            fltr_w2[fltr_w2] = cr_filter_w2

            self.st_fltr_w1.append(fltr_w1)
            self.st_fltr_w2.append(fltr_w2)
            
            no_w1 += np.sum(fltr_w1)
            no_w2 += np.sum(fltr_w2)
        
        if verbose:

            print 'number of filtered detection events in each window w1 / w2: ', no_w1, ' / ', no_w2
            print


    def attach_state_filtered_syncs(self,apply_dt_filter = True, max_dt = None, only_one_outcome = False, ZeroOrOne = 0, verbose = True):
        """
        checks for the signatures of psi_minus or psi_plus and returns a list of numpy arrays where each numpy array corresponds the correpsonding sync number of the event
        also has the ability to filter by the dt between the two events
        """

        self.st_fltr_w1_ch1     = []
        self.st_fltr_w1_ch0     = []
        self.st_fltr_w2_ch1     = []
        self.st_fltr_w2_ch0     = []
        self.st_fltr_psi_plus   = []
        self.st_fltr_psi_minus  = []

        self.HH_sync_psi_plus    = []
        self.HH_sync_psi_minus   = []

        i = 0

        for fltr_w1,fltr_w2,channels,HH_sync,HH_time in zip(self.st_fltr_w1,self.st_fltr_w2,self.lt4_dict['/PQ_channel-1'],self.lt4_dict['/PQ_sync_number-1'],self.lt4_dict['/PQ_sync_time-1']):

            st_fltr_w1_ch1 = fltr_w1 & (channels == 1)
            st_fltr_w1_ch0 = fltr_w1 & (channels == 0)
            st_fltr_w2_ch1 = fltr_w2 & (channels == 1)
            st_fltr_w2_ch0 = fltr_w2 & (channels == 0)

            ### store filters in object for potential later processing.
            self.st_fltr_w1_ch1.append(cp.deepcopy(st_fltr_w1_ch1))
            self.st_fltr_w1_ch0.append(cp.deepcopy(st_fltr_w1_ch0))
            self.st_fltr_w2_ch1.append(cp.deepcopy(st_fltr_w2_ch1))
            self.st_fltr_w2_ch0.append(cp.deepcopy(st_fltr_w2_ch0))


            #### formulate filters according to the relevant state
            #### for psi_plus: the same detector has to click within one sync
            #### for psi_minus: different detectors have to click
            #### the filter w1 is shifted onto the filter of w2 by inserting two boolean Falses at the beginning (clicks must be consecutive accross windows)
            #### Need extra False to get past the special in between the two valid clicks
            st_fltr_w2_ch1 = np.append(st_fltr_w2_ch1,[False,False])[2:]
            st_fltr_w2_ch0 = np.append(st_fltr_w2_ch0,[False,False])[2:]

            st_fltr_psi_plus = np.logical_or(np.logical_and(st_fltr_w1_ch1,st_fltr_w2_ch1),np.logical_and(st_fltr_w1_ch0,st_fltr_w2_ch0))
            st_fltr_psi_minus = np.logical_or(np.logical_and(st_fltr_w1_ch0,st_fltr_w2_ch1),np.logical_and(st_fltr_w1_ch1,st_fltr_w2_ch0))

            st_fltr_psi_plus_1 = np.logical_and(st_fltr_w1_ch1,st_fltr_w2_ch1)
            st_fltr_psi_plus_0 = np.logical_and(st_fltr_w1_ch0,st_fltr_w2_ch0)
            st_fltr_psi_minus_0 = np.logical_and(st_fltr_w1_ch0,st_fltr_w2_ch1)
            st_fltr_psi_minus_1 = np.logical_and(st_fltr_w1_ch1,st_fltr_w2_ch0)
            
            if only_one_outcome:
                st_fltr_psi_plus = st_fltr_psi_plus_0 if ZeroOrOne == 0 else st_fltr_psi_plus_1
                st_fltr_psi_minus = st_fltr_psi_minus_0 if ZeroOrOne == 0 else st_fltr_psi_minus_1

            if apply_dt_filter:

                if max_dt == None:
                    max_dt = analysis_params.filter_settings['max_dt']


                shift_fltr_psi_plus = np.insert(st_fltr_psi_plus,0,[False,False])[:-2]
                shift_fltr_psi_minus = np.insert(st_fltr_psi_minus,0,[False,False])[:-2]
                
                st_fltr_psi_plus[st_fltr_psi_plus] = np.absolute(HH_time[st_fltr_psi_plus]-HH_time[shift_fltr_psi_plus]) <= max_dt
                st_fltr_psi_minus[st_fltr_psi_minus] = np.absolute(HH_time[st_fltr_psi_minus]-HH_time[shift_fltr_psi_minus]) <= max_dt
                
            self.st_fltr_psi_plus.append(st_fltr_psi_plus)
            self.st_fltr_psi_minus.append(st_fltr_psi_minus)

            if verbose:
                print 'Run ', i+1
                print 'total events w1 / w2', np.sum(st_fltr_w1_ch0)+np.sum(st_fltr_w1_ch1),'/',np.sum(st_fltr_w2_ch0)+np.sum(st_fltr_w2_ch1)
                print 'psi_minus events', np.sum(st_fltr_psi_minus)
                print 'psi_plus events', np.sum(st_fltr_psi_plus)
                print

            ### the required HH syncs are later used to get the corresponding adwin ROs
            self.HH_sync_psi_plus.append(HH_sync[st_fltr_psi_plus])
            self.HH_sync_psi_minus.append(HH_sync[st_fltr_psi_minus])

            i += 1

    

    def correlate_RO_results(self,apply_ROC = False, verbose = True,return_value = False):

        self.RO_data_LT3_plus = []
        self.RO_data_LT3_minus = []
        self.RO_data_LT4_plus = []
        self.RO_data_LT4_minus = []

        for HH_s_psi_p,HH_s_psi_m,adwin_nrs_w1,adwin_ro_lt3,adwin_ro_lt4 in zip(self.HH_sync_psi_plus,self.HH_sync_psi_minus,self.lt4_dict['counted_awg_reps_w1'],self.lt3_dict['ssro_results'],self.lt4_dict['ssro_results']):


            fltr_plus  = self.filter_adwin_data_from_pq_syncs(HH_s_psi_p,adwin_nrs_w1)
            fltr_minus = self.filter_adwin_data_from_pq_syncs(HH_s_psi_m,adwin_nrs_w1)


            self.RO_data_LT3_plus.append(adwin_ro_lt3[fltr_plus])
            self.RO_data_LT4_plus.append(adwin_ro_lt4[fltr_plus])
            self.RO_data_LT3_minus.append(adwin_ro_lt3[fltr_minus]) 
            self.RO_data_LT4_minus.append(adwin_ro_lt4[fltr_minus])

        all_m_lt3,all_m_lt4,all_p_lt3,all_p_lt4 = np.array([]),np.array([]),np.array([]),np.array([])
        for m_lt3,m_lt4,p_lt3,p_lt4 in zip(self.RO_data_LT3_minus,self.RO_data_LT4_minus,self.RO_data_LT3_plus,self.RO_data_LT4_plus):
            all_m_lt3 = np.append(all_m_lt3,m_lt3)
            all_m_lt4 = np.append(all_m_lt4,m_lt4)
            all_p_lt3 = np.append(all_p_lt3,p_lt3)
            all_p_lt4 = np.append(all_p_lt4,p_lt4)


        ### get overall events psi minus:
        m_correlations = [0,0,0,0]
        m_correlations[0], m_correlations[1], m_correlations[2] , m_correlations[3] = correlate_arrays(all_m_lt3,al_m_lt4)

        ### get overall events psi plus:
        p_correlations = [0,0,0,0]
        p_correlations[0], p_correlations[1], p_correlations[2], p_correlations[3] = correlate_arrays(all_p_lt3,all_p_lt4)

        if verbose:
            print 'The occurence of each event after filtering'
            print

            x = ['ms0 & ms0', 'ms0 & ms1', 'ms1 & ms0 ','ms1 & ms1 ',]
            row_format ="{:>12}" * (len(x) + 1)
            headline_format = "{:>15}"+"{:>12}" * len(x)
            print headline_format.format("", *x)


            

            for state, row in zip(['psi_minus'], [m_correlations]):
                print "-"*(12*4+15)
                print row_format.format(state+' |', *row)

            print


            for state, row in zip(['psi_plus'], [p_correlations]):
                print "-"*(12*4+15)
                print row_format.format(state+' |', *row)

        if apply_ROC: 
            ### get ssro_ROC for LT3 --> corresponds to setup B
            F0_LT3,F1_LT3 = self.find_RO_fidelities(self.ROC_lt3_tstamp,self.lt3_dict['data_attrs'][0],folder = self.lt3_folder)
            ### get ssro_ROC for LT4 --> corresponds to setup A
            F0_LT4,F1_LT4 = self.find_RO_fidelities(self.ROC_lt4_tstamp,self.lt4_dict['data_attrs'][0],folder = self.lt4_folder)

            ### apply ROC to the results --> input arrays for this function have to be reversed!
            corrected_psi_minus,u_minus = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(m_correlations[::-1]),F0_LT4,F0_LT3,F1_LT4,F1_LT3,verbose = verbose,return_error_bars = True)
            corrected_psi_plus, u_plus  = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(p_correlations[::-1]),F0_LT4,F0_LT3,F1_LT4,F1_LT3,verbose = verbose,return_error_bars = True)

            ### take care of the total number of correlations
            no_m = np.sum(m_correlations)
            raw_m_correlations = m_correlations
            raw_p_correlations = p_correlations
            no_p = np.sum(p_correlations)

            ### sscorr returns the array an inverted order ms 11, ms 10, ms 01, ms 00

            m_correlations = np.array(list(reversed(corrected_psi_minus.reshape(-1)))) ### so much recasting!
            p_correlations = np.array(list(reversed(corrected_psi_plus.reshape(-1)))) ### so much recasting!
            u_minus = np.array(list(reversed(u_minus.reshape(-1))))
            u_plus = np.array(list(reversed(u_plus.reshape(-1)))) ### so much recasting!

            if return_value:
                return m_correlations,u_minus,p_correlations,u_plus,no_m,no_p

            else:
                self.plot_corrected_RO_results(m_correlations,u_minus,raw_m_correlations,state = 'minus')
                print 'Results_psi_minus'
                print m_correlations
                print u_minus
                self.plot_corrected_RO_results(p_correlations,u_plus,raw_p_correlations,state = 'plus')
                print 'Results_psi_plus'
                print p_correlations
                print u_plus
        if return_value:
            return m_correlations,[],p_correlations,[],np.sum(m_correlations),np.sum(p_correlations)

    def find_RO_fidelities(self,timestamp,data_attrs,folder = ''):

        ssro_folder = tb.data_from_time(timestamp,folder = folder) 
        # print ssro_folder
        if 'MWInit' in ssro_folder:
            e_trans_string = data_attrs['electron_transition']
            e_trans_string = 'ms'+e_trans_string[1:]
            F_0,u_F0,F_1,u_F1 = ssro.get_SSRO_MWInit_calibration(ssro_folder,data_attrs['E_RO_durations'][0],e_trans_string)

        else:
            F_0,u_F0,F_1,u_F1 = ssro.get_SSRO_calibration(ssro_folder,data_attrs['E_RO_durations'][0]) #we onlky have one RO time in E_RO_durations

        return F_0,F_1



    def sweep_filter_parameter_vs_correlations(self,parameter_name,parameter_range,apply_ROC = False,do_plot = True):


        ### initialize results lists
        no_of_minus_events,no_of_plus_events = [],[]
        minus_correlation,plus_correlation   = [],[]
        minus_correlation_u,plus_correlation_u = [],[]
        minus_full_corrs, plus_full_corrs = [], []
        minus_full_corrs_u, plus_full_corrs_u = [], []

        last_x = 0

        for x in parameter_range:

            if parameter_name == 'bin_w2':
                analysis_params.filter_settings['max_reps_w2'] = x
                analysis_params.filter_settings['min_reps_w2'] = last_x + 1
                last_x = x
                
            else:
                analysis_params.filter_settings[parameter_name] = x ### commence sweep

            self.apply_temporal_filters_to_prefiltered_data(verbose = False)
            self.apply_sync_filter_w1_w2(verbose = False)
            self.apply_is_purified_filter(verbose = False)
            self.apply_CR_before_filter(verbose=False)
            self.apply_CR_after_filter(verbose=False)
            self.attach_state_filtered_syncs(verbose = False)
            psi_m_corrs,minus_u, psi_p_corrs,plus_u,no_m,no_p = self.correlate_RO_results(verbose=False,return_value = True,apply_ROC = apply_ROC)

            no_of_minus_events.append(no_m)
            no_of_plus_events.append(no_p)

            plus_full_corrs.append(psi_p_corrs)
            minus_full_corrs.append(psi_m_corrs)  
            plus_full_corrs_u.append(plus_u) 
            minus_full_corrs_u.append(minus_u) 

            no_anti_correlations_m = float(psi_m_corrs[1]+psi_m_corrs[2])
            no_correlations_m = float(psi_m_corrs[0]+psi_m_corrs[3])
            no_anti_correlations_p = float(psi_p_corrs[1]+psi_p_corrs[2])
            no_correlations_p = float(psi_p_corrs[0]+psi_p_corrs[3])

            
            Tomo = self.get_tomography_bases()


            if np.sum(psi_m_corrs) == 0:
                minus_correlation.append(0)
                minus_correlation_u.append(0)

            elif not apply_ROC:
                minus_correlation_u.append(np.sqrt((no_anti_correlations_m*(no_correlations_m**2)+no_correlations_m*(no_anti_correlations_m**2)))/((no_correlations_m + no_anti_correlations_m)**2))
                minus_correlation.append(float(no_anti_correlations_m)/np.sum(psi_m_corrs))
            else:
                if (Tomo == 'XX') or (Tomo == 'YY'):
                    #### we do ROC and return expectation values directly!!!
                    no_anti_correlations_m = (no_anti_correlations_m-0.5)*2
                    
                    no_anti_correlations_m,m_u = do_carbon_ROC(no_anti_correlations_m,np.sqrt(minus_u[1]**2+minus_u[2]**2)*2)
                    minus_correlation.append(no_anti_correlations_m)
                    minus_correlation_u.append(m_u)
                else:
                    no_correlations_m = (no_correlations_m-0.5)*2
                    no_correlations_m,m_u = do_carbon_ROC(no_correlations_m,np.sqrt(minus_u[0]**2+minus_u[3]**2)*2)
                    minus_correlation.append(no_correlations_m)
                    minus_correlation_u.append(m_u)

            if np.sum(psi_p_corrs) == 0:
                plus_correlation.append(0)
                plus_correlation_u.append(0)
            elif not apply_ROC:
                plus_correlation.append(float(no_correlations_p)/np.sum(psi_p_corrs))
                plus_correlation_u.append(np.sqrt((no_anti_correlations_p*(no_correlations_p**2)+no_correlations_p*(no_anti_correlations_p**2)))/((no_correlations_p + no_anti_correlations_p)**2))
            else:
                #### we do ROC and return expectation values directly!!!
                if (Tomo == 'XX') or (Tomo == 'ZZ'):
                    no_anti_correlations_p = (no_anti_correlations_p-0.5)*2
                    no_anti_correlations_p,p_u = do_carbon_ROC(no_anti_correlations_p,np.sqrt(plus_u[1]**2+plus_u[2]**2)*2)
                    plus_correlation.append(no_anti_correlations_p)
                    plus_correlation_u.append(p_u) 

                else:
                    no_correlations_p = (no_correlations_p-0.5)*2
                    no_correlations_p,p_u = do_carbon_ROC(no_correlations_p,np.sqrt(plus_u[0]**2+plus_u[3]**2)*2)
                    plus_correlation.append(no_correlations_p)
                    plus_correlation_u.append(p_u) 

            ### error bars are based on poissonian statistics for correlated events vs. uncorrelated
                        
        ### commence plotting

        ### check if we are dealing with timing or something else (like CR after...)
        if max(parameter_range) > 1000:
            x = parameter_range/1000
            x_units = ' (ns)'
        else:
            x = parameter_range
            x_units = ''

        reload(analysis_params) ### to negate the sweep changes

        if do_plot:
            self.create_plot(title = 'Number of events within filter', xlabel = parameter_name + x_units,ylabel = 'Occurences')
            self.plot_data(x,no_of_minus_events,label = '-')
            self.plot_data(x,no_of_plus_events,label = '+')
            plt.legend()
            self.save_and_close_plot()


            ### should really have an error calculation for the minus / plus correlation
            self.create_plot(title = 'Fraction of correct correlations', xlabel = parameter_name + x_units, ylabel = 'Expectation value')
            self.plot_data(x,minus_correlation,y_u = minus_correlation_u,label = '-')
            self.plot_data(x,plus_correlation,y_u = plus_correlation_u,label = '+')
            plt.legend(loc = 2)
            self.save_and_close_plot()

        else: # if we do not plot we assume that values are supposed to be returned
            return [np.array(no_of_minus_events),np.array(minus_correlation),np.array(minus_correlation_u)],[np.array(no_of_plus_events),np.array(plus_correlation),np.array(plus_correlation_u)],[np.transpose(np.array(minus_full_corrs)),np.transpose(np.array(minus_full_corrs_u))],[np.transpose(np.array(plus_full_corrs)),np.transpose(np.array(plus_full_corrs_u))]

    def filter_adwin_data_from_pq_syncs(self,filtered_sn,counted_awg_reps):
        """
        takes the filtered pq syncs as input and returns a boolean array.
        This array serves as filter for the adwin RO results

        TODO: generalize for arbitrary PQ data size (introduce 'index')
        """

        if np.sum(np.logical_not(np.in1d(filtered_sn,counted_awg_reps))) != 0:
                    print 'Connecting pq syncs to adwin data seems to be going wrong!'

        # print 'elen', len(filtered_sn)
        insert_pos = np.searchsorted(counted_awg_reps,filtered_sn)
        #insert_pos = np.searchsorted(filtered_sn,adwin_syncs)
        return insert_pos


    def correlate_RO_results_no_purification(self,apply_ROC = False, verbose = True,return_value = False):
        """
        very similar to to correlate_RO_results:
        Ignores the existence of psi_p and psi_m and only takes st_fltr_w1 into account.
        Is mainly used for population measurements such as eta_vs_theta
        """

        self.RO_data_LT3 = []
        self.RO_data_LT4 = []


        for sync_nrs,fltr_w1,adwin_nrs_w1,adwin_ro_lt3,adwin_ro_lt4 in zip(self.lt4_dict['/PQ_sync_number-1'],\
                                                                    self.st_fltr_w1,self.lt4_dict['counted_awg_reps_w1'],\
                                                                    self.lt3_dict['ssro_results'],\
                                                                    self.lt4_dict['ssro_results']):

            fltr  = self.filter_adwin_data_from_pq_syncs(sync_nrs[fltr_w1],adwin_nrs_w1)
            self.RO_data_LT3.append(adwin_ro_lt3[fltr])
            self.RO_data_LT4.append(adwin_ro_lt4[fltr])


        all_lt3,all_lt4 = np.array([]),np.array([])
        for r_lt3,r_lt4 in zip(self.RO_data_LT3,self.RO_data_LT4):

            all_lt3 = np.append(all_lt3,r_lt3)
            all_lt4 = np.append(all_lt4,r_lt4)

        #### print correlation matrix for RO results
        ### 


        ### get overall events psi minus:
        correlations = [0,0,0,0]
        correlations[0],correlations[1],correlations[2],correlations[3] = correlate_arrays(all_lt3,all_lt4)
        # print m_correlations
        # print headline_format.format("", *x)


        if verbose:
            print 'The occurence of each event after filtering'
            print

            x = ['ms0 & ms0', 'ms0 & ms1', 'ms1 & ms0 ','ms1 & ms1 ',]
            row_format ="{:>12}" * (len(x) + 1)
            headline_format = "{:>15}"+"{:>12}" * len(x)
            print headline_format.format("", *x)


            

            for state, row in zip(['ZZ'], [correlations]):
                print "-"*(12*4+15)
                print row_format.format(state+' |', *row)


        if apply_ROC: 
            ### get ssro_ROC for LT3 --> corresponds to setup B
            F0_LT3,F1_LT3 = self.find_RO_fidelities(self.ROC_lt3_tstamp,self.lt3_dict['data_attrs'][0],folder = self.lt3_folder)
            ### get ssro_ROC for LT4 --> corresponds to setup A
            F0_LT4,F1_LT4 = self.find_RO_fidelities(self.ROC_lt4_tstamp,self.lt4_dict['data_attrs'][0],folder = self.lt4_folder)

            ### apply ROC to the results --> input arrays for this function have to be reversed!
            corrected_corrs,u_corrected_corrs = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(correlations[::-1]),F0_LT4,F0_LT3,F1_LT4,F1_LT3,verbose = verbose,return_error_bars = True)

            ### take care of the total number of correlations
            no = np.sum(correlations)
            raw_correlations = correlations
            ### sscorr returns the array an inverted order ms 11, ms 10, ms 01, ms 00

            corrected_corrs = np.array(list(reversed(corrected_corrs.reshape(-1)))) ### so much recasting!
            u_corrected_corrs = np.array(list(reversed(u_corrected_corrs.reshape(-1))))

            if return_value:
                return corrected_corrs,u_corrected_corrs,no

            else:
                self.plot_corrected_RO_results(corrected_corrs,u_corrected_corrs,raw_correlations,state='one_click')
                print 'these are the results'
                print corrected_corrs
                print u_corrected_corrs
        if return_value:
            return m_correlations,[],np.sum(m_correlations)

    #############################################
    #### density matrix and state estimation ####
    #############################################

    def get_DM_correlations(self,verbose = True,apply_ROC = False):
        """
        Generates a dictionary of all correlations (and associated errors) for plus and minus signature.
        Dictionary keys are tomography bases: e.g. XY
        If apply_ROC: besides the probabilities to detect each correlation (four columns) we also obtain a fifth column that gives the total amount of counts (comes in handy for any ML estimation of the state)

        TODO: IMPLEMENT CARBON ROC
        """

        tomos = ['I','X','Y','Z']
        ### init correlation dictionary
        self.correlation_dict_m,self.correlation_dict_m_u  = {},{}
        self.correlation_dict_p,self.correlation_dict_p_u  = {},{}
        for t1 in tomos:
            for t2 in tomos:
                if 'I' in (t1+t2):
                    tomo_length = 2 ## single qubit values
                else:
                    tomo_length = 4
                self.correlation_dict_m.update({t1+t2:np.zeros(tomo_length+apply_ROC)})
                self.correlation_dict_p.update({t1+t2:np.zeros(tomo_length+apply_ROC)})
                self.correlation_dict_m_u.update({t1+t2:np.zeros(tomo_length)})
                self.correlation_dict_p_u.update({t1+t2:np.zeros(tomo_length)})

        ### init other lists
        self.RO_data_LT3_plus,self.RO_data_LT3_minus = [],[]
        self.RO_data_LT4_plus,self.RO_data_LT4_minus = [],[]
        self.LT3_tomo,self.LT4_tomo = [],[]

        ### filter RO data from HH syncs
        for HH_s_psi_p,HH_s_psi_m,adwin_nrs_w1,adwin_ro_lt3,adwin_ro_lt4,attrs_LT3,attrs_LT4 in zip(self.HH_sync_psi_plus,self.HH_sync_psi_minus,self.lt4_dict['counted_awg_reps_w1'], \
                                                                                self.lt3_dict['ssro_results'],self.lt4_dict['ssro_results'], \
                                                                                self.lt3_dict['data_attrs'],self.lt4_dict['data_attrs']):


            fltr_plus  = self.filter_adwin_data_from_pq_syncs(HH_s_psi_p,adwin_nrs_w1)
            fltr_minus = self.filter_adwin_data_from_pq_syncs(HH_s_psi_m,adwin_nrs_w1)

            self.RO_data_LT3_plus.append(adwin_ro_lt3[fltr_plus])
            self.RO_data_LT4_plus.append(adwin_ro_lt4[fltr_plus])
            self.RO_data_LT3_minus.append(adwin_ro_lt3[fltr_minus]) 
            self.RO_data_LT4_minus.append(adwin_ro_lt4[fltr_minus])
            self.LT3_tomo.append(attrs_LT3['Tomography_bases'][0])
            self.LT4_tomo.append(attrs_LT4['Tomography_bases'][0])


        for m_lt3,m_lt4,p_lt3,p_lt4,t_LT3,t_LT4 in zip(self.RO_data_LT3_minus, self.RO_data_LT4_minus, self.RO_data_LT3_plus,\
                                                        self.RO_data_LT4_plus, self.LT3_tomo, self.LT4_tomo):

            basis = t_LT3,t_LT4
            m11,m10,m01,m00 = correlate_arrays(m_lt3,m_lt4)

            self.correlation_dict_m[basis] += np.array([m11,m10,m01,m00,np.sum([m11,m10,m01,m00])]) 

            p11,p10,p01,p00 = correlate_arrays(p_lt3,p_lt4)

            self.correlation_dict_p[basis] += np.array([p11,p10,p01,p00,np.sum([p11,p10,p01,p00])]) 

        
        if apply_ROC:
            ### get ssro_ROC for LT3 --> corresponds to setup B
            F0_LT3,F1_LT3 = self.find_RO_fidelities(self.ROC_lt3_tstamp,self.lt3_dict['data_attrs'][0],folder = self.lt3_folder)
            ### get ssro_ROC for LT4 --> corresponds to setup A
            F0_LT4,F1_LT4 = self.find_RO_fidelities(self.ROC_lt4_tstamp,self.lt4_dict['data_attrs'][0],folder = self.lt4_folder)

            ### apply ROC to the results --> input arrays for this function have to be reversed!
            for k in self.correlation_dict_m.keys():
                corrected_psi_minus,u_minus = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(self.correlation_dict_m[k][3::-1]),F0_LT4,F0_LT3,F1_LT4,F1_LT3,verbose = verbose,return_error_bars = True)
                corrected_psi_plus, u_plus  = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(self.correlation_dict_p[k][3::-1]),F0_LT4,F0_LT3,F1_LT4,F1_LT3,verbose = verbose,return_error_bars = True)
        

                self.correlation_dict_m[k] = np.array(list(reversed(corrected_psi_minus.reshape(-1)))+[self.correlation_dict_m[k][-1]])
                self.correlation_dict_p[k] = np.array(list(reversed(corrected_psi_plus.reshape(-1)))+[self.corelation_dict_p[k][-1]]) 
                self.correlation_dict_m_u[k] = np.array(list(reversed(u_minus.reshape(-1))))
                self.correlation_dict_p_u[k] = np.array(list(reversed(u_plus.reshape(-1))))
        else:
            print 'NO ROC is not implemented for the density matrix yet!'

        #### TODO
        #### IMPLEMENT CARBON ROC FOR THE GIVEN CASE


        #### after getting the correlated results we need to determine contributions such as IX and ZI (single qubit contributions)
        ### trace out LT4 first and get expectation values for LT3
        for t in ['X','Y','Z']:
            cumulative_prob_m,cumulative_prob_p = np.array([0,0,0,0,0]),np.array([0,0,0,0,0])
            cumulative_prob_m_u,cumulative_prob_p_u = np.array([0,0,0,0]),np.array([0,0,0,0])
            for t2 in ['X','Y','Z']:
                cumulative_prob_m   += self.correlation_dict_m[t+t2]
                cumulative_prob_m_u += self.corelation_dict_m_u[t+t2]**2/9.
                cumulative_prob_p   += self.corelation_dict_p[t+t2]
                cumulative_prob_p_u += self.corelation_dict_p_u[t+t2]**2/9.

            cumulative_prob_m_u = np.sqrt(cumulative_prob_m_u) ## error propagation on the individual errorbars
            cumulative_prob_p_u = np.sqrt(cumulative_prob_p_u) ## error propagation on the individual errorbars

            self.correlation_dict_m[t+'I']   = np.array([(cumulative_prob_m[0]+cumulative_prob_m[1])/3.,(cumulative_prob_m[2]+cumulative_prob_m[3])/3.,cumulative_prob_m[4]])
            self.correlation_dict_p[t+'I']   = np.array([(cumulative_prob_p[0]+cumulative_prob_p[1])/3.,(cumulative_prob_p[2]+cumulative_prob_p[3])/3.,cumulative_prob_p[4]])
            self.correlation_dict_m_u[t+'I'] = np.array([np.sqrt(cumulative_prob_m_u[0]**2+cumulative_prob_m_u[1]**2),(cumulative_prob_m_u[2]**2+cumulative_prob_m_u[3]**2)])
            self.correlation_dict_p_u[t+'I'] = np.array([np.sqrt(cumulative_prob_p_u[0]**2+cumulative_prob_p_u[1]**2),(cumulative_prob_p_u[2]**2+cumulative_prob_p_u[3]**2)])

        #### now trace out LT3 and get expectation values for LT3
        for t in ['X','Y','Z']:
            cumulative_prob_m,cumulative_prob_p = np.array([0,0,0,0,0]),np.array([0,0,0,0,0])
            cumulative_prob_m_u,cumulative_prob_p_u = np.array([0,0,0,0]),np.array([0,0,0,0])
            for t2 in ['X','Y','Z']:
                cumulative_prob_m   += self.correlation_dict_m[t2+t]
                cumulative_prob_m_u += self.corelation_dict_m_u[t2+t]**2/9.
                cumulative_prob_p   += self.corelation_dict_p[t2+t]
                cumulative_prob_p_u += self.corelation_dict_p_u[t2+t]**2/9.

            cumulative_prob_m_u = np.sqrt(cumulative_prob_m_u) ## error propagation on the individual errorbars
            cumulative_prob_p_u = np.sqrt(cumulative_prob_p_u) ## error propagation on the individual errorbars

            self.correlation_dict_m['I'+t]   = np.array([(cumulative_prob_m[0]+cumulative_prob_m[2])/3.,(cumulative_prob_m[1]+cumulative_prob_m[3])/3.,cumulative_prob_m[4]])
            self.correlation_dict_p['I'+t]   = np.array([(cumulative_prob_p[0]+cumulative_prob_p[2])/3.,(cumulative_prob_p[1]+cumulative_prob_p[3])/3.,cumulative_prob_p[4]])
            self.correlation_dict_m_u['I'+t] = np.array([np.sqrt(cumulative_prob_m_u[0]**2+cumulative_prob_m_u[2]**2),(cumulative_prob_m_u[1]**2+cumulative_prob_m_u[3]**2)])
            self.correlation_dict_p_u['I'+t] = np.array([np.sqrt(cumulative_prob_p_u[0]**2+cumulative_prob_p_u[2]**2),(cumulative_prob_p_u[1]**2+cumulative_prob_p_u[3]**2)])


    def reconstruct_DMs(self,verbose = True,max_likelihood = False,aply_ROC = False):
        """
        reconstructs the density matrices for the plus and minus signature. 
        Options include maximum likelihood estimation of the DM from the original number of counts
        This step also applies the carbon RO correction.

        TODO: implement max_likelihood estimation
        """
        paulis = self.generate_pauli_matrices()


        if max_likelihood:
            print 'Not implemented yet!'
        else:
            dm_p,dm_m = np.kron(paulis[0],paulis[0])/4.,np.kron(paulis[0],paulis[0])/4.
            dm_p_u,dm_m_u = np.zeros((4,4)),np.zeros((4,4))

            t_dict = {'I' : 0, 'X':1, 'Y':2, 'Z':3}
            for t in ['I','X','Y','Z']:
                for t2 in ['I','X','Y','Z']:
                    
                    if t+t2 == 'II':
                        continue
                        
                    sigma_kron = np.kron(paulis[t_dict[t]],paulis[t_dict[t2]])

                    ### single or two qubit expectation value?
                    if 'I' in t+t2:
                        exp_p,exp_p_u = do_carbon_ROC(get_1q_expectation_val(self.correlation_dict_p[t+t2]),self.correlation_dict_p_u[t+t2])) 
                        exp_m,exp_m_u = do_carbon_ROC(get_1q_expectation_val(self.correlation_dict_m[t+t2]),self.correlation_dict_m_u[t+t2]))
                    else:
                        exp_p,exp_p_u = do_carbon_ROC(get_2q_expectation_val(self.correlation_dict_p[t+t2]), self.correlation_dict_p_u[t+t2])) 
                        exp_m,exp_m_u = do_carbon_ROC(get_2q_expectation_val(self.correlation_dict_m[t+t2]), self.correlation_dict_m_u[t+t2])) 

                    dm_p +=exp_p*sigma_kron/4.
                    dm_p_u += (exp_p_u**2)*sigma_kron/16.
                    dm_m +=exp_m*sigma_kron/4.
                    dm_m_u += (exp_m_u**2)*sigma_kron/16.


        dm_p_u = np.sqrt(dm_p_u.real)+np.sqrt(dm_p_u.imaginary)*1j
        dm_m_u = np.sqrt(dm_p_u.real)+np.sqrt(dm_p_u.imaginary)*1j

        if verbose:
            print 'Density matrix for the positive signature'
            print dm_p
            print 'Error matrix'
            print dm_p_u
            print 'Eigenvalues'
            print np.linalg.eigh(dm_p)
            print 'Density matrix for the negative signature'
            print dm_m
            print 'Error matrix'
            print dm_m_u
            print 'Eigenvalues'
            print np.linalg.eigh(dm_m)

        else:
            return dm_p,dm_p_u,dm_m,dm_m_u


    ##########################
    #### timing and rates ####
    ##########################

    def get_total_time(self):
        """
        loops over all raw lt4 data to extract the total time the experiment had been running
        Input: ///
        Output: total_time (in seconds)
        """

        return np.sum(self.lt4_dict['total_elapsed_time'])


    def estimate_sequence_time(self):
        """
        Approximates the time spent in the AWG sequence:
        instead of calculating the exact time spent in the awg sequence one can make a few basic assumptions to obtain a duration 
        a) the number of trials give the number of LDE attempts
        b) per 1000 attempts we do the full 5 carbon gates and two RO triggers of 90 us duration (the number 1000 is pure speculatioN!!!)
        c) keep in mind that one could get the exact duration spent in the sequence by looping over the syncnumbers and plu markers to extract what happened
        Output: total_time (in seconds)
        """

        No_of_pulses = self.lt4_dict['data_attrs'][0]['C4_Ren_N_m1'][0]
        tau = self.lt4_dict['data_attrs'][0]['C4_Ren_tau_m1'][0]

        total_syncs = np.sum(self.lt4_dict['total_syncs'])
        total_duration = total_syncs*self.lt4_dict['joint_attrs'][0]['LDE_element_length']

        return total_duration

    def calculate_sequence_time(self,lt4_timestamps,st_start=None,st_len=None,max_w2 = None):
        """
        Approximates the time spent in the AWG sequence:
        instead of calculating the exact time spent in the awg sequence one can make a few basic assumptions to obtain a duration 
        a) the number of trials give the number of LDE attempts
        b) per 1000 attempts we do the full 5 carbon gates and two RO triggers of 90 us duration (the number 1000 is pure speculatioN!!!)
        c) keep in mind that one could get the exact duration spent in the sequence by looping over the syncnumbers and plu markers to extract what happened
        Output: total_time (in seconds)
        """

        num_raw_successes, num_successes, total_attempts_to_first_click, num_first_clicks, est_resets, total_time, total_syncs = 0,0,0,0,0,0,0

        #for each file
        for t_lt4 in lt4_timestamps:
                   
            a_lt4 = ppq.purifyPQAnalysis(tb.latest_data(t_lt4,folder = self.lt4_folder),hdf5_mode ='r')
            
            if 'LDE1_attempts' in a_lt4.g.attrs:
                LDE1_attempts = a_lt4.g.attrs['LDE1_attempts']
            else:
                LDE1_attempts = 1000

            attempts_first = a_lt4.agrp['attempts_first'].value
            attempts_second = a_lt4.agrp['attempts_second'].value

            awg_reps_w1 = (np.array(a_lt4.agrp['counted_awg_reps'].value)-np.array(a_lt4.agrp['attempts_second'].value))
            awg_reps_w2 = np.array(a_lt4.agrp['counted_awg_reps'].value)

            syncs = np.array(a_lt4.pqf['/PQ_sync_number-1'].value)
            spec = np.array(a_lt4.pqf['/PQ_special-1'].value)
            time = np.array(a_lt4.pqf['/PQ_time-1'].value)
            sync_time = np.array(a_lt4.pqf['/PQ_sync_time-1'].value)

            plu_clicks = np.append((spec == 1),False)[1:] # Plu click is followed by a special value
            plu_syncs = syncs[plu_clicks]
            plu_sync_time = sync_time[plu_clicks]

            ### one can also apply manual filters if one wants to deviate from the prescribed parameter dictionary
            if st_start == None:
                st_start = analysis_params.filter_settings['st_start']
            if st_len == None:
                st_len = analysis_params.filter_settings['st_len']

            st_fltr = (plu_sync_time > st_start)  & (plu_sync_time < (st_start  + st_len))
            plu_syncs = plu_syncs[st_fltr]

            successful_pur_w1_sync = np.in1d(plu_syncs, awg_reps_w1) # Check if led to successful purification
            successful_pur_w2_sync = np.in1d(plu_syncs, awg_reps_w2) # Check if led to successful purification

            first_clicks = plu_syncs[np.logical_not(successful_pur_w2_sync)]

            # To estimate this, we only use data where we know we succeeded, since then no ambiguity about what happened.
            next_sync = plu_syncs[np.insert(successful_pur_w2_sync,0,[False])[:-1]] 
            elapsed_syncs = next_sync - plu_syncs[successful_pur_w2_sync][:len(next_sync)]


            # Finally, for rates, add in filtering on the number of allowed attempts
            max_w1 = analysis_params.filter_settings['max_reps_w1']
            min_w2 = analysis_params.filter_settings['min_reps_w2']
            if max_w2 == None:
                max_w2 = analysis_params.filter_settings['max_reps_w2']

            num_successes += np.sum(np.in1d(plu_syncs,awg_reps_w1[np.logical_and((attempts_first <= max_w1),(attempts_second <= max_w2),(attempts_second >= min_w2))]))

            # Add it all up:
            num_raw_successes += np.sum(successful_pur_w1_sync)
            total_attempts_to_first_click += (np.sum(elapsed_syncs))
            num_first_clicks += len(first_clicks)
            total_time += time[-1]/1e12
            total_syncs += syncs[-1]

        No_of_pulses = a_lt4.g.attrs['C4_Ren_N_m1'][0]
        tau = a_lt4.g.attrs['C4_Ren_tau_m1'][0]
        C13_manipulation_duration = 2*No_of_pulses*tau+1*90e-6
        reset_duration = 2*C13_manipulation_duration
        msmt_duration = 2*C13_manipulation_duration
        store_duration = 2*C13_manipulation_duration
        LDE_elem_length = a_lt4.joint_grp.attrs['LDE_element_length'] 

        avg_attempts_per_first_click = float(total_attempts_to_first_click)/num_raw_successes
        resets_per_first_click = (float(avg_attempts_per_first_click)/LDE1_attempts + 1)
        
        crude_click_prob = num_first_clicks/float(total_syncs)
        first_click_prob = 1/float(avg_attempts_per_first_click)

        # Times per success
        est_succ_prob_given_first_click = (1-(1-first_click_prob)**max_w2)
        meas_succ_prob_given_first_click = num_successes/float(num_first_clicks)

        succ_prob_given_first_click = meas_succ_prob_given_first_click

        entanglement_time = (1/first_click_prob)*LDE_elem_length * (1/succ_prob_given_first_click)
        reset_time = resets_per_first_click*reset_duration * (1/succ_prob_given_first_click)
        store_time = store_duration * (1/succ_prob_given_first_click)
        msmt_time = msmt_duration
        est_time_per_success = entanglement_time + reset_time + store_time + msmt_time

        # B and K comparisons
        tail_lt3 = 7.2 #9.1975 for SIL 3
        tail_lt4 = 5.1514
        click_prob_lt4 = first_click_prob*tail_lt4/(tail_lt4 + tail_lt4)
        click_prob_lt3 = first_click_prob*tail_lt3/(tail_lt4 + tail_lt4)

        prob_successive_clicks =  (0.5*(click_prob_lt4*click_prob_lt3) + 0.25*click_prob_lt4**2 + 0.25*click_prob_lt3**2)
        expected_disting_clicks = prob_successive_clicks * total_syncs

        print 'num_raw_successes: ', num_raw_successes
        print 'num filtered successes: ', num_successes
        print 'num_first_clicks: ', num_first_clicks
        print 'num_syncs: ', total_syncs

        print '\n'
        print 'crude click_prob (first_clicks/syncs): ', crude_click_prob
        print 'meas. first click prob: ', first_click_prob
        print 'avg_attempts_per_first_click: ', avg_attempts_per_first_click
        print 'est_resets_per_first_click: ', resets_per_first_click
        print 'est. success prob given first click: ', est_succ_prob_given_first_click
        print 'meas. success prob given first click: ', meas_succ_prob_given_first_click
        print 'inferred storage success rate (per side): ', np.sqrt(meas_succ_prob_given_first_click/est_succ_prob_given_first_click)
        print '\n'
        print 'Times per successful event:'
        print 'entanglement_time: ',entanglement_time 
        print 'reset_time: ',reset_time
        print 'store_time: ',store_time
        print 'msmt_time: ',msmt_time
        print 'est_total_time: ', num_successes * est_time_per_success
        print 'elapsed_total_time: ', total_time
        
        print '\n'
        print 'est. rate (inc. purific. succ.): ',(1/8.0)*(1/est_time_per_success)
        print 'crude rate from number of syncs (inc. purific. succ.): ',(1/8.0)*num_successes/(total_syncs*LDE_elem_length)
        print 'actual rate (inc. purific. succ.): ',(1/8.0)*float(num_successes)/total_time
        print 'est. duty cycle: ',(est_time_per_success*float(num_successes))/(total_time)
    
        print '\n'
        print 'TPQI expected disting clicks: ', expected_disting_clicks

    def get_sequence_time(self):
        """
        TO BE WRITTEN: this function gets the time we spent in a data file as accurate as possible
        """
    #######################################
    #### helper functions for plotting ####
    #######################################
    def generate_pauli_matrices(self):
        return [np.matrix([[1,0],[0,1]],dtype=complex),np.matrix([[0,1],[1,0]],dtype=complex),np.matrix([[0,-1j],[1j,0]],dtype=complex),np.matrix([[1,0],[0,-1]],dtype=complex)]
    

    def get_tomography_bases(self):
        """
        makes two strings for the tomography bases as chosen by both setups
        """


        try: ### this only work if we really read-out a nuclear spin at the very end and not an electron spin.
            LT3 = self.lt3_dict['data_attrs'][0]['Tomography_bases'][0]
            LT4 = self.lt4_dict['data_attrs'][0]['Tomography_bases'][0]
        except:
            print 'e spin experiment detected: using the MW amplitudes to determine tomography bases'
            LT3 = self.lt3_dict['data_attrs'][0]['LDE_final_mw_amplitude']
            LT4 = self.lt4_dict['data_attrs'][0]['LDE_final_mw_amplitude']

            if LT3 > 0:
                LT3 = 'X'
            else:
                LT3 = 'Z'
            if LT4 > 0:
                LT4 = 'X'
            else: 
                LT4 = 'Z'

        return LT3+LT4



    def num2str(self, num, precision): 
        return "%0.*f" % (precision, num)

    def create_plot(self,**kw):
        ylabel = kw.pop('ylabel',None)
        xlabel = kw.pop('xlabel',None)
        title = kw.pop('title',None)

        fig = plt.figure()
        ax = plt.subplot()

        if xlabel != None:
            plt.xlabel(xlabel)
        else:
            plt.xlabel('Run (#)')

        if ylabel != None:
            plt.ylabel(ylabel)
        else:
            plt.ylabel('Contrast')

        if title != None:
            plt.title(title)

        return fig,ax
        
    def plot_corrected_RO_results(self,corr,err=np.array([0,0,0,0]),events = np.array([0,0,0,0]),state = ''):

        tomo = self.get_tomography_bases()
        fig,ax = self.create_plot(title = '',xlabel = '',ylabel = 'Probability')
        # plt.title('Corrected')
        lw,fontsize = self.get_lw_and_fontsize()

        plt.bar(np.arange(0,4), np.squeeze(np.asarray(corr)), yerr=err, align = 'center',
                color = '#436CDE',linewidth = lw,error_kw=dict(lw =lw,capthick=lw,ecolor = 'black'))
        for k in range(4):
            plt.text(k,0.01, self.num2str(np.squeeze(np.asarray(corr))[k],2), ha = 'center')
            plt.text(k,err[k]+corr[k]+0.02, self.num2str(np.squeeze(np.asarray(events))[k],0), ha = 'center')

        if np.amax(corr)+np.amax(err) > 0.5:
            plt.ylim([0,np.amax(corr)+np.amax(err)])
        else:
            plt.ylim([0,0.5])
        x_tick1 = tomo[0]
        x_tick2 = tomo[1]
        plt.xticks([0,1,2,3], [ '|'+x_tick1+','+x_tick2+'>', '|'+x_tick1+','+'-'+x_tick2+'>', '|'+'-'+x_tick1+','+x_tick2+'>','|'+'-'+x_tick1+','+'-'+x_tick2+'>'])
        self.format_plot(fig,ax)
        ax.tick_params(top = 'off',bottom = 'off')
        self.save_and_close_plot(save = True,name = 'Purify_correlations_'+tomo+'_'+state)

    def save_and_close_plot(self,f = save_plot_folder, save = False, name = None):

        if name == None:
            name = 'Results'

        if save:
            plt.savefig(os.path.join(f,name+'.pdf'),format='pdf',bbox_inches='tight')
            plt.savefig(os.path.join(f,name+'.png'),format='png',bbox_inches='tight')

        plt.show()

        plt.close('all')

    def plot_data(self,x,y,**kw):
        label = kw.pop('label',None)
        y_u = kw.pop('y_u',None)
        if y_u != None:
            plt.errorbar(x,y,y_u,fmt = 'x',label = label,**kw)
        else: plt.plot(x,y,'x',label = label)


    def format_plot(self,fig,ax):
        """
        all formatting of standard plots goes here.
        """
        linewidths, textsize = self.get_lw_and_fontsize()

        ax.tick_params(labelsize = textsize,width =linewidths)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(linewidths)
        ax.yaxis.label.set_size(textsize)
        ax.xaxis.label.set_size(textsize)


    def get_lw_and_fontsize(self):
        """
        returns linewidths and font sizes
        this function is used for consistency reasons!
        """
        linewidth = 1.5
        textsize = 14
        return linewidth,textsize
    def tstamps_for_both_setups(self,day_string,contains = 'XX',newest_tstamp = '235959'):
        """
        takes a date as input and scans lt3 and lt4 for appropriate timestamps
        will throw an error if both setups have run the experiment an unequal amount of times!
        --> then you have to clean up the data folder of the setup with a longer list of timestamps
        input: day_string, e.g. '20160607'
        output: lt3_t_list,lt4_t_list
        """

        lt3_t_list = self.find_tstamps_of_day([],day_string,contains=contains,analysis_folder = self.lt3_folder ,newest_tstamp = newest_tstamp)
        lt4_t_list = self.find_tstamps_of_day([],day_string,contains=contains,analysis_folder = self.lt4_folder, newest_tstamp = newest_tstamp)

        return self.verify_tstamp_lists(lt3_t_list,lt4_t_list,day_string)


    def find_tstamps_of_day(self,ts_list,day_string,contains='XX',analysis_folder = 'throw exception',newest_tstamp = '235959'):

        latest_t = day_string + newest_tstamp # where in the day do you want to begin? 235959 mean: take the whole day
        newer_than = day_string+'_000000'

        while tb.latest_data(contains,older_than=latest_t,folder=analysis_folder,newer_than=newer_than,return_timestamp = True,raise_exc = False) != False:

            latest_t,f = tb.latest_data(contains,older_than=latest_t,folder=analysis_folder,newer_than=newer_than,return_timestamp = True,raise_exc=False)
            
            ### debug statement that prints the full timestamp and the relevant identifier.
            # print latest_t[8:],latest_t

            ### append found timestamp to list of timestamps
            ts_list.append(latest_t[8:]) 

        return ts_list

    def verify_tstamp_lists(self,lt3_t_list,lt4_t_list,date):

        # print len(lt3_t_list),len(lt4_t_list)
        print 'verifying timestamps of',date
        if len(lt3_t_list) != len(lt4_t_list):
            print 't_lt3 , t_lt4'
            for lt3_t,lt4_t in zip(lt3_t_list,lt4_t_list):
                print lt3_t,lt4_t
            raise Exception('The length of the time stamp lists is unequal. Clean out the data folders on each computer t3_0,t4_0: ',lt3_t_list[0],lt4_t_list[0]) 

        clean_t_list_lt3,clean_t_list_lt4 = [],[]

        ### check for contents
        newer_than = date+'_000000'
        for t_lt3,t_lt4 in zip(lt3_t_list,lt4_t_list):
            f_lt3 = tb.latest_data(t_lt3,folder = self.lt3_folder,newer_than = newer_than)
            f_lt4 = tb.latest_data(t_lt4,folder = self.lt4_folder,newer_than = newer_than)

            ## get file path and open file
            Datafile_lt3 = h5py.File(f_lt3+f_lt3[len(self.lt3_folder)+9:] + '.hdf5','r') 
            Datafile_lt4 = h5py.File(f_lt3+f_lt3[len(self.lt4_folder)+9::] + '.hdf5','r') 

            if (not u'PQ_hist' in Datafile_lt3.keys()) or (not u'PQ_hist' in Datafile_lt4.keys()): ### did we actually detect any clicks what so ever?
                continue

            ### we exploit the fact that .keys is ordered according to what is saved last.
            ### --> adwin data is saved last: therefore if the last key contains pq --> no adwin data

            if ('PQ' in Datafile_lt3.keys()[-1]) or ('PQ' in Datafile_lt4.keys()[-1]):
                continue

            ### check if there are adwin ssro events and if they have the same length for both files
            ssros_lt3 = Datafile_lt4[Datafile_lt4.keys()[-1]]['adwindata']['ssro_results'].value
            ssros_lt4 = Datafile_lt4[Datafile_lt4.keys()[-1]]['adwindata']['ssro_results'].value

            if (len(ssros_lt4) == len(ssros_lt3)) and (len(ssros_lt3) !=0):
                ### timestamp contains valuable information, add to clean lists
                clean_t_list_lt3.append(t_lt3)
                clean_t_list_lt4.append(t_lt4)
    
        ### return clean timestamp lists

        return clean_t_list_lt3,clean_t_list_lt4

    def plot_quantity_for_raw_data(self):
        """
        To be written. 
        loops over all elements in raw data (lt3_lt4)
        Input a function that operates on raw data objects (such as pq type objects) and the parameters associated with that function
        The function should return a 1D array of values associated with the raw data.
        plot_quantity_for_raw_data then groups the returned arrays element-wise and plots them
        """
        pass

    def plot_quantity_for_filtered_data(self):
        """
        TO BE WRITTEN.
        see above but takes entries of the filtered dictionary which are fed into the function alongside parameters.
        one example for correlations would be: parameters: windows
        dictionary entries: 'adwin_ssro', 'counted_awg_reps', 'all PQ data per entry'
        """
        pass


    # def format_figure(fig,scaling):
    #     ## writes a bunch of parameters to matplotlib config.
    #     ## input: effective figure scaling factor
    #     mpl.rcParams['axes.linewidth'] = 1.8*figscaling
    #     mpl.rcParams['xtick.major.width'] = 1.8*figscaling
    #     mpl.rcParams['ytick.major.width'] = 1.8*figscaling
    #     mpl.rcParams['font.size'] = (22-8*(figscaling-0.66*0.72)/(1-0.66*0.72))*figscaling
    #     mpl.rcParams['axes.titlesize'] = mpl.rcParams['font.size']
    #     mpl.rcParams['legend.fontsize'] = mpl.rcParams['font.size']
    #     mpl.rcParams['legend.labelspacing'] = 0.5*figscaling
    #     mpl.rcParams['legend.columnspacing'] = 1.5*figscaling
    #     mpl.rcParams['legend.handletextpad'] = 0.3*figscaling
    #     mpl.rcParams['legend.handlelength'] = 1.*figscaling
    #     mpl.rcParams['legend.borderpad'] = 0.2+0.2*(figscaling-0.66*0.72)/(1-0.66*0.72)
    #     mpl.rcParams['lines.markersize'] = 13.5*figscaling#7*figscaling 
    #     mpl.rcParams['figure.figsize'] = (6*figscaling,4*figscaling)
    #     mpl.rcParams['lines.markeredgewidth'] = 0.3/figscaling



def do_carbon_ROC(expectation_value,uncertainty):
    """
    takes a two-partite expectation value and applies carbon read-out correction 
    Returns the corrected expectation value and a new uncertainty via error propagation
    """
    #### these values were measured on 15-08-2016
    ### see onenote carbon control LT3/LT4 for details
    ROC = 0.951786
    uROC = 0.0033490

    u = np.sqrt((ROC*uncertainty)**2+(uROC*expectation_value)**2)/ROC**2

    return expectation_value/ROC,u

def calculate_ebits(parity_ZZ,parity_YY,correlations_XX):


    c0110, p00, p01, p10, p11 = sp.symbols('c0110 p00 p01 p10 p11')

    del_p00,del_p01,del_p10,del_p11, del_c0110 = sp.symbols('del_p00,del_p01,del_p10,del_p11, del_c0110')

    logNegSimp = parse_expr('log(p01+p10+sqrt(2 *c0110**2+p11**2+(-1+p01+p10+p11)**2+(-1+p01+p10)* sqrt(4 *c0110**2+(-1+p01+p10+2 *p11)**2))/sqrt(2)+sqrt(2 *c0110**2+p00**2+p11**2+sqrt((-1+p01+p10)**2 *(4* c0110**2+(-1+p01+p10+2 *p11)**2)))/sqrt(2))/log(2)')

    error_p00 = parse_expr('p00/(sqrt(2) * sqrt(2 *c0110**2+p00**2+p11**2+sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2 *p11)**2)))* (p01+p10+sqrt(2 *c0110**2+p11**2+(-1+p01+p10+p11)**2+(-1+p01+p10) *sqrt(4 *c0110**2+(-1+p01+p10+2 *p11)**2))/sqrt(2)+sqrt(2 *c0110**2+p00**2+p11**2+sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2 *p11)**2)))/sqrt(2))* log(2))')
    error_p01 = parse_expr('(8+(2 *sqrt(2)* (-1+p01+p10)* (4 *c0110**2+(-1+p01+p10)* (-1+p01+p10+2* p11)+(-1+p01+p10+2* p11)**2))/(sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2* p11)**2))* sqrt(2* c0110**2+p00**2+p11**2+sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2* p11)**2))))+(4 *sqrt(2)* (2 *c0110**2+(-1+p01+p10+p11)* (-1+p01+p10+2 *p11+sqrt(4 *c0110**2+(-1+p01+p10+2 *p11)**2))))/(sqrt(4 *c0110**2+(-1+p01+p10+2 *p11)**2) * sqrt(2 *c0110**2+p11**2+(-1+p01+p10+p11)**2+(-1+p01+p10) *sqrt(4* c0110**2+(-1+p01+p10+2* p11)**2))))/(8* (p01+p10+sqrt(2* c0110**2+p11**2+(-1+p01+p10+p11)**2+(-1+p01+p10)* sqrt(4* c0110**2+(-1+p01+p10+2* p11)**2))/sqrt(2)+sqrt(2* c0110**2+p00**2+p11**2+sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2 *p11)**2)))/sqrt(2)) *log(2))')
    error_p10 = parse_expr('(8+(2 *sqrt(2)* (-1+p01+p10)* (4 *c0110**2+(-1+p01+p10)* (-1+p01+p10+2* p11)+(-1+p01+p10+2* p11)**2))/(sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2* p11)**2))* sqrt(2* c0110**2+p00**2+p11**2+sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2* p11)**2))))+(4 *sqrt(2)* (2 *c0110**2+(-1+p01+p10+p11)* (-1+p01+p10+2 *p11+sqrt(4 *c0110**2+(-1+p01+p10+2 *p11)**2))))/(sqrt(4 *c0110**2+(-1+p01+p10+2 *p11)**2) * sqrt(2 *c0110**2+p11**2+(-1+p01+p10+p11)**2+(-1+p01+p10) *sqrt(4* c0110**2+(-1+p01+p10+2* p11)**2))))/(8* (p01+p10+sqrt(2* c0110**2+p11**2+(-1+p01+p10+p11)**2+(-1+p01+p10)* sqrt(4* c0110**2+(-1+p01+p10+2* p11)**2))/sqrt(2)+sqrt(2* c0110**2+p00**2+p11**2+sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2 *p11)**2)))/sqrt(2)) *log(2))')
    error_p11 = parse_expr('((2 *(-1+p01+p10+2 *p11+((-1+p01+p10)* (-1+p01+p10+2 *p11))/sqrt(4 *c0110**2+(-1+p01+p10+2* p11)**2)))/sqrt(2 *c0110**2+p11**2+(-1+p01+p10+p11)**2+(-1+p01+p10)* sqrt(4* c0110**2+(-1+p01+p10+2* p11)**2))+(2* (p11+((-1+p01+p10)**2 *(-1+p01+p10+2* p11))/sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2 *p11)**2))))/sqrt(2* c0110**2+p00**2+p11**2+sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2* p11)**2))))/(2 * sqrt(2) *(p01+p10+sqrt(2 *c0110**2+p11**2+(-1+p01+p10+p11)**2+(-1+p01+p10) *sqrt(4* c0110**2+(-1+p01+p10+2 *p11)**2))/sqrt(2)+sqrt(2* c0110**2+p00**2+p11**2+sqrt((-1+p01+p10)**2 *(4* c0110**2+(-1+p01+p10+2* p11)**2)))/sqrt(2)) *log(2))')
    error_c0110 = parse_expr('(sqrt(2) *c0110 * ((1+(-1+p01+p10)/sqrt(4 *c0110**2+(-1+p01+p10+2 *p11)**2))/sqrt(2 *c0110**2+p11**2+(-1+p01+p10+p11)**2+(-1+p01+p10)* sqrt(4 *c0110**2+(-1+p01+p10+2* p11)**2))+(1+(-1+p01+p10)**2/sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2* p11)**2)))/sqrt(2 *c0110**2+p00**2+p11**2+sqrt((-1+p01+p10)**2 *(4 *c0110**2+(-1+p01+p10+2 *p11)**2)))))/((p01+p10+sqrt(2 *c0110**2+p11**2+(-1+p01+p10+p11)**2+(-1+p01+p10) *sqrt(4 *c0110**2+(-1+p01+p10+2 *p11)**2))/sqrt(2)+sqrt(2* c0110**2+p00**2+p11**2+sqrt((-1+p01+p10)**2 *(4* c0110**2+(-1+p01+p10+2 *p11)**2)))/sqrt(2)) *log(2))')


    error =  sp.sqrt((error_p00*del_p00)**2+(error_p01*del_p01)**2+(error_p10*del_p10)**2 +\
                     (error_p11*del_p11)**2 + (error_c0110*del_c0110)**2)

    logNegF = sp.lambdify([c0110, p00, p01, p10, p11],logNegSimp)
    logNegU = sp.lambdify([c0110, p00, p01, p10, p11,del_c0110, del_p00,del_p01,del_p10,del_p11],error)

    c0110 = (parity_ZZ[1] + parity_YY[1])/float(4)# Define correlations as even - odd
    c0110_u = np.sqrt(parity_ZZ[2]**2 + parity_YY[2]**2)/float(4)

    p00,p01,p10,p11 = np.array(correlations_XX[0][0]),np.array(correlations_XX[0][1]),np.array(correlations_XX[0][2]),np.array(correlations_XX[0][3])
    del_p00,del_p01,del_p10,del_p11 = np.array(correlations_XX[1][0]),np.array(correlations_XX[1][1]),np.array(correlations_XX[1][2]),np.array(correlations_XX[1][3])

    ebits = []
    ebits_u = []
    for c0110_ind,p00_ind,p01_ind,p10_ind,p11_ind,c0110_u_ind,del_p00_ind,del_p01_ind,del_p10_ind,del_p11_ind in zip(c0110,p00,p01,p10,p11,c0110_u,del_p00,del_p01,del_p10,del_p11):

        ebits.append(logNegF(c0110_ind,p00_ind,p01_ind,p10_ind,p11_ind))
        ebits_u.append(logNegU(c0110_ind,p00_ind,p01_ind,p10_ind,p11_ind,c0110_u_ind,del_p00_ind,del_p01_ind,del_p10_ind,del_p11_ind))
    
    return np.array(ebits), np.array(ebits_u)

def convert_attrs_to_dict(attrManagerItems):
    data_attrs = dict([])
    for key,value in attrManagerItems:
        data_attrs[key] = value
    return data_attrs

def check_if_file_exists(filePath):
    file_exists = True
    try:
        f = h5py.File(filePath, 'r')
        f.close()
    except:
        file_exists = False

    return file_exists

def correlate_arrays(lt3,lt4):
    """
    assumes binary numpy arrays
    """
    m11 = np.sum(np.equal(    lt3[lt3 == 1],lt4[lt3 == 1]))
    m10 = np.sum(np.not_equal(lt3[lt3 == 1],lt4[lt3 == 1]))
    m01 = np.sum(np.not_equal(lt3[lt3 == 0],mt4[lt3 == 0]))
    m00 = np.sum(np.equal(    lt3[lt3 == 0],lt4[lt3 == 0]))

    return m11,m10,m01,m00

def get_1q_expectation_val(arr):
    """ 
    assumption:
    arr[0] --> bright/up state
    arr[1] --> dark/down state
    """
    return 2*(arr[0]-0.5),2*arr_u[0]


def get_2q_expectation_val(arr):
    """ 
    assumption:
    arr[0] --> brightbright/upup state
    arr[1] --> brightdark/updown state
    arr[2] --> darkbright/downup state
    arr[3] --> darkdark/downdown state
    """
    return 2*(arr[0]+arr[3]-0.5),2*np.sqrt(arr_u[0]**2+arr_u[3]**2)
