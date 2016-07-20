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
from analysis.lib.lde import sscorr; reload(sscorr)
import purify_analysis_params as analysis_params;reload(analysis_params)

#### standard parameters

save_plot_folder = r'D:\measuring\data\purification_data'


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

		self.key_list = ['/PQ_sync_number-1','/PQ_channel-1','/PQ_sync_time-1','/PQ_special-1','counted_awg_reps_w1','counted_awg_reps_w2','ssro_results','tstamp','raw_data']

		self.key_list_pq =  ['/PQ_sync_number-1','/PQ_channel-1','/PQ_sync_time-1','/PQ_special-1']
		self.key_lst_adwin_data = ['ssro_results','electron_readout_result','Phase_correction_repetitions','CR_after','CR_before','carbon_readout_result','attempts_first','attempts_second']
		self.lt4_dict = {}
		self.lt3_dict = {}

		### initialize the data dictionaries
		for key in self.key_list:
			self.lt3_dict.update({key:[]})
			self.lt4_dict.update({key:[]})

		for key in self.key_lst_adwin_data:
			self.lt3_dict.update({key:[]})
			self.lt4_dict.update({key:[]})

		self.lt3_folder = lt3_folder
		self.lt4_folder = lt4_folder

		self.ROC_lt3_tstamp = ROC_lt3_tstamp
		self.ROC_lt4_tstamp = ROC_lt4_tstamp

	def load_raw_data(self,lt3_timestamps,lt4_timestamps):
		"""
		this script takes a list of timestamps for both setups and prefilters them according to adwin filters
		creates a list of arrays associated with each time stamp for adwin_ssro,syncs, time, special, channel, marker and attaches is it to the data object for further processing.
		"""
		length = len(lt3_timestamps)
		print 'loading the data, total number of files ', length

		i = 0

		#Reinit arrays
		self.lt3_dict['counted_awg_reps_w2'] = []
		self.lt3_dict['counted_awg_reps_w1'] = []
		self.lt4_dict['counted_awg_reps_w2'] = []
		self.lt4_dict['counted_awg_reps_w1'] = []

		for t_lt3,t_lt4 in zip(lt3_timestamps,lt4_timestamps):
			
			# print self.lt3_folder
			
			a_lt3 = ppq.purifyPQAnalysis(tb.latest_data(t_lt3,folder = self.lt3_folder),hdf5_mode ='r')
			a_lt4 = ppq.purifyPQAnalysis(tb.latest_data(t_lt4,folder = self.lt4_folder),hdf5_mode ='r')

			### filter the timeharp data according to adwin events / syncs

			### need to create two sync filters here.!
			sync_filter = a_lt4.filter_pq_data_from_adwin_syncs() ### this syncfilter erases all data where from the PQ data where the adwin did NOT read out

			if len(sync_filter) == 0: # empty list --> no read outs.
				# print 'file empty, skipping these time stamps:',t_lt3,t_lt4
				# print
				continue

			### store relevant adwin results
			syncs_w1 = (np.array(a_lt3.agrp['counted_awg_reps'].value)-np.array(a_lt3.agrp['attempts_second'].value))
			syncs_w2 = np.array(a_lt3.agrp['counted_awg_reps'].value)


			if len(syncs_w1) == 0:
				print 'syncs_w1 empty, skipping these time stamps:',t_lt3,t_lt4
				print
				continue

			self.lt3_dict['counted_awg_reps_w2'].append(syncs_w2)
			self.lt3_dict['counted_awg_reps_w1'].append(syncs_w1)
			self.lt4_dict['counted_awg_reps_w2'].append(syncs_w2)
			self.lt4_dict['counted_awg_reps_w1'].append(syncs_w1)



			### need to create two sync filters here. This is done by appending them to each other and sorting the relatively short array. (should never be more than 500 events??)
			### this syncfilter erases all data from the PQ data where the adwin did NOT read out. Drastically reduces the data amount we have to handle.
			sync_filter = a_lt4.filter_pq_data_from_adwin_syncs(np.sort(np.append(syncs_w1,syncs_w2)))
			if len(sync_filter) == 0: # empty list --> no read outs.
				print 'file empty, skipping these time stamps:',t_lt3,t_lt4
				print
				continue



			self.lt3_dict['tstamp'].append(t_lt3)
			self.lt4_dict['tstamp'].append(t_lt4)
			self.lt3_dict['raw_data'].append(a_lt3)
			self.lt4_dict['raw_data'].append(a_lt4)

			for key in self.key_lst_adwin_data:
				self.lt4_dict[key].append(np.array(a_lt4.agrp[key].value))

			for key in self.key_list_pq:
				self.lt4_dict[key].append(np.array(a_lt4.pqf[key].value[sync_filter]))

			for key in self.key_lst_adwin_data:
				self.lt3_dict[key].append(np.array(a_lt3.agrp[key].value))

			#### calculate the duty cycle for that specific file.
			# print 'lde length',a_lt3.joint_grp.attrs['LDE_element_length']
			# print'first and last time',a_lt4.pqf['/PQ_time-1'].value[0],a_lt4.pqf['/PQ_time-1'][-1]
			# print 'last elapsed time in sequence vs total elapsed time',  

			time_in_LDE_sequence = a_lt3.joint_grp.attrs['LDE_element_length']*self.lt4_dict['/PQ_sync_number-1'][-1][-1]
			total_elapsed_time = (a_lt4.pqf['/PQ_time-1'].value[-1]-a_lt4.pqf['/PQ_time-1'][0])*1e-12

			print 'file no ', i+1 , ' with duty cycle of', round(100*time_in_LDE_sequence/total_elapsed_time,1), ' %'
			i += 1



	

	def check_tail_w1_w2(self,st_start = 2000e3,st_len = 50e3):
		"""
		goes through the raw_data and selects only files according to the applied sync filter: returns the measured tail in each window
		for the purification experiment there will only be one window
		"""

		i = 0
		tails_w1 = []

		for a in self.lt4_dict['raw_data']:
					
			### analysis for channel 0	&& window 1
			w1_ch0 = self.get_total_number_of_clicks_in_window(a,0,st_start,st_len)
			w1_ch1 = self.get_total_number_of_clicks_in_window(a,1,st_start,st_len)
			last_sync = a.pqf['/PQ_sync_number-1'][-1]

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


	def attach_state_filtered_syncs(self,apply_dt_filter = True, max_dt = None, verbose = True):
		"""
		checks for the signatures of psi_minus or psi_plus and returns a list of numpy arrays where each numpy array corresponds the correpsonding sync number of the event
		also has the ability to filter by the dt between the two events
		"""

		self.st_fltr_w1_ch1 	= []
		self.st_fltr_w1_ch0 	= []
		self.st_fltr_w2_ch1 	= []
		self.st_fltr_w2_ch0 	= []
		self.st_fltr_psi_plus 	= []
		self.st_fltr_psi_minus 	= []

		self.HH_sync_psi_plus 	 = []
		self.HH_sync_psi_minus	 = []

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


			#### formulate filteres according to the relevant state
			#### for psi_plus: the same detector has to click within one sync
			#### for psi_minus: different detectors have to click
			#### the filter w1 is shifted onto the filter of w2 by inserting two boolean Falses at the beginning (clicks must be consecutive accross windows)
			#### Need extra False to get past the special in between the two valid clicks
			st_fltr_w2_ch1 = np.append(st_fltr_w2_ch1,[False,False])[2:]
			st_fltr_w2_ch0 = np.append(st_fltr_w2_ch0,[False,False])[2:]

			st_fltr_psi_plus = np.logical_or(np.logical_and(st_fltr_w1_ch1,st_fltr_w2_ch1),np.logical_and(st_fltr_w1_ch0,st_fltr_w2_ch0))
			st_fltr_psi_minus = np.logical_or(np.logical_and(st_fltr_w1_ch0,st_fltr_w2_ch1),np.logical_and(st_fltr_w1_ch1,st_fltr_w2_ch0))

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

			# print np.argmax(fltr_plus)
			# print HH_s_psi_p[np.argmax(fltr_plus)]
			# print np.argmax(fltr_minus)
			# print HH_s_psi_m[np.argmax(fltr_minus)]
			# print adwin_nrs_w1

			self.RO_data_LT3_plus.append(adwin_ro_lt3[fltr_plus])
			self.RO_data_LT4_plus.append(adwin_ro_lt4[fltr_plus])
			self.RO_data_LT3_minus.append(adwin_ro_lt3[fltr_minus])	
			self.RO_data_LT4_minus.append(adwin_ro_lt4[fltr_minus])

		# print fltr_plus_lt3
		# print fltr_plus_lt4


		all_m_lt3,all_m_lt4,all_p_lt3,all_p_lt4 = np.array([]),np.array([]),np.array([]),np.array([])
		for m_lt3,m_lt4,p_lt3,p_lt4 in zip(self.RO_data_LT3_minus,self.RO_data_LT4_minus,self.RO_data_LT3_plus,self.RO_data_LT4_plus):

			# if len(m_lt3) != 0 and len(m_lt4) != 0:
			# 	if verbose:
			# 		print 'p_correlated for psi_minus', float(np.sum(np.equal(m_lt3,m_lt4)))/len(m_lt3)
			# 		print 'p_correlated for psi_plus', float(np.sum(np.equal(p_lt3,p_lt4)))/len(m_lt4)

			all_m_lt3 = np.append(all_m_lt3,m_lt3)
			all_m_lt4 = np.append(all_m_lt4,m_lt4)
			all_p_lt3 = np.append(all_p_lt3,p_lt3)
			all_p_lt4 = np.append(all_p_lt4,p_lt4)



		#### print correlation matrix for RO results
		### 


		### get overall events psi minus:
		m_correlations = [0,0,0,0]
		m_correlations[0] = np.sum(np.equal(all_m_lt3[all_m_lt3 == 1],all_m_lt4[all_m_lt3 == 1]))
		m_correlations[1] = np.sum(np.not_equal(all_m_lt3[all_m_lt3 == 1],all_m_lt4[all_m_lt3 == 1]))
		m_correlations[2] = np.sum(np.not_equal(all_m_lt3[all_m_lt3 == 0],all_m_lt4[all_m_lt3 == 0]))
		m_correlations[3] = np.sum(np.equal(all_m_lt3[all_m_lt3 == 0],all_m_lt4[all_m_lt3 == 0]))
		# print m_correlations

		# print headline_format.format("", *x)
		p_correlations = [0,0,0,0]
		p_correlations[0] = np.sum(np.equal(all_p_lt3[all_p_lt3 == 1],all_p_lt4[all_p_lt3 == 1]))
		p_correlations[1] = np.sum(np.not_equal(all_p_lt3[all_p_lt3 == 1],all_p_lt4[all_p_lt3 == 1]))
		p_correlations[2] = np.sum(np.not_equal(all_p_lt3[all_p_lt3 == 0],all_p_lt4[all_p_lt3 == 0]))
		p_correlations[3] = np.sum(np.equal(all_p_lt3[all_p_lt3 == 0],all_p_lt4[all_p_lt3 == 0]))

		if verbose:
			print 'The occurence of each event after filtering'
			print

			x = ['ms0 & ms0', 'ms0 & ms1', 'ms1 & ms0','ms1 & ms1',]
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
			F0_LT3,F1_LT3 = self.find_RO_fidelities(self.ROC_lt3_tstamp,self.lt3_dict['raw_data'][0],folder = self.lt3_folder)
			### get ssro_ROC for LT4 --> corresponds to setup A
			F0_LT4,F1_LT4 = self.find_RO_fidelities(self.ROC_lt4_tstamp,self.lt4_dict['raw_data'][0],folder = self.lt4_folder)

			### apply ROC to the results --> input arrays for this function have to be reversed!
			corrected_psi_minus,u_minus = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(m_correlations[::-1]),F0_LT4,F0_LT3,F1_LT4,F1_LT3,verbose = verbose,return_error_bars = True)
			corrected_psi_plus, u_plus  = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(p_correlations[::-1]),F0_LT4,F0_LT3,F1_LT4,F1_LT3,verbose = verbose,return_error_bars = True)

			### take care of the total number of correlations
			no_m = np.sum(m_correlations)
			no_p = np.sum(p_correlations)

			### sscorr returns the array an inverted order ms 11, ms 10, ms 01, ms 00

			m_correlations = np.array(list(reversed(corrected_psi_minus.reshape(-1)))) ### so much recasting!
			p_correlations = np.array(list(reversed(corrected_psi_plus.reshape(-1)))) ### so much recasting!
			u_minus = np.array(list(reversed(u_minus.reshape(-1))))
			u_plus = np.array(list(reversed(u_plus.reshape(-1)))) ### so much recasting!

			if return_value:
				return m_correlations,u_minus,p_correlations,u_plus,no_m,no_p
		if return_value:
			return m_correlations,[],p_correlations,[],np.sum(m_correlations),np.sum(p_correlations)

	def find_RO_fidelities(self,timestamp,raw_data,folder = ''):

		ssro_folder = tb.data_from_time(timestamp,folder = folder) 
		# print ssro_folder

		F_0,u_F0,F_1,u_F1 = ssro.get_SSRO_calibration(ssro_folder,raw_data.g.attrs['E_RO_durations'][0]) #we onlky have one RO time in E_RO_durations

		return F_0,F_1
	def sweep_filter_parameter_vs_correlations(self,parameter_name,parameter_range,apply_ROC = False):


		### initialize results lists
		no_of_minus_events,no_of_plus_events = [],[]
		minus_correlation,plus_correlation	 = [],[]
		minus_correlation_u,plus_correlation_u = [],[]
		minus_u,plus_u 							= [],[]
		for x in parameter_range:

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

			no_anti_correlations_m = float(psi_m_corrs[1]+psi_m_corrs[2])
			no_correlations_m = float(psi_m_corrs[0]+psi_m_corrs[3])
			no_anti_correlations_p = float(psi_p_corrs[1]+psi_p_corrs[2])
			no_correlations_p = float(psi_p_corrs[0]+psi_p_corrs[3])

			if no_anti_correlations_p>1:
				print psi_p_corrs
			if np.sum(psi_m_corrs) == 0:
				minus_correlation.append(0)
				minus_correlation_u.append(0)

			elif not apply_ROC:
				minus_correlation_u.append(np.sqrt((no_anti_correlations_m*(no_correlations_m**2)+no_correlations_m*(no_anti_correlations_m**2)))/((no_correlations_m + no_anti_correlations_m)**2))
				minus_correlation.append(float(no_anti_correlations_m)/np.sum(psi_m_corrs))
			else:
				if no_anti_correlations_m>no_correlations_m:
					minus_correlation.append(no_anti_correlations_m)
					minus_correlation_u.append(np.sqrt(minus_u[1]**2+minus_u[2]**2))
				else:
					minus_correlation.append(no_correlations_m)
					minus_correlation_u.append(np.sqrt(minus_u[0]**2+minus_u[3]**2))

			if np.sum(psi_p_corrs) == 0:
				plus_correlation.append(0)
				plus_correlation_u.append(0)
			elif not apply_ROC:
				plus_correlation.append(float(no_correlations_p)/np.sum(psi_p_corrs))
				plus_correlation_u.append(np.sqrt((no_anti_correlations_p*(no_correlations_p**2)+no_correlations_p*(no_anti_correlations_p**2)))/((no_correlations_p + no_anti_correlations_p)**2))
			else:
				if no_anti_correlations_p>no_correlations_p:
					plus_correlation.append(no_anti_correlations_p)
					plus_correlation_u.append(np.sqrt(plus_u[1]**2+plus_u[2]**2)) 

				else:
					plus_correlation.append(no_correlations_p)
					plus_correlation_u.append(np.sqrt(plus_u[0]**2+plus_u[3]**2)) 

			### error bars are based on poissonian statistics for correlated events vs. uncorrelated
						
		### commence plotting

		### check if we are dealing with timing or something else (like CR after...)
		if max(parameter_range) > 1000:
			x = parameter_range/1000
			x_units = ' (ns)'
		else:
			x = parameter_range
			x_units = ''


		self.create_plot(title = 'Number of events within filter', xlabel = parameter_name + x_units,ylabel = 'Occurences')
		self.plot_data(x,no_of_minus_events,label = '-')
		self.plot_data(x,no_of_plus_events,label = '+')
		plt.legend()
		self.save_and_close_plot()


		### should really have an error calculation for the minus / plus correlation
		self.create_plot(title = 'Fraction of correct correlations', xlabel = parameter_name + x_units, ylabel = 'p_right_correlation')
		self.plot_data(x,minus_correlation,y_u = minus_correlation_u,label = '-')
		self.plot_data(x,plus_correlation,y_u = plus_correlation_u,label = '+')
		plt.legend(loc = 2)
		self.save_and_close_plot()

		reload(analysis_params) ### to negate the sweep changes

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

	#######################################
	#### helper functions for plotting ####
	#######################################


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
		
		
	def save_and_close_plot(self,f = save_plot_folder, save = False, name = None):

		if name == None:
			name = 'Results'
		if save:
			plt.savefig(os.path.join(f,'name.pdf'),format='pdf')
			plt.savefig(os.path.join(f,'name.png'),format='png')
			
		plt.show()
		plt.close('all')

	def plot_data(self,x,y,**kw):
		label = kw.pop('label',None)
		y_u = kw.pop('y_u',None)
		if y_u != None:
			plt.errorbar(x,y,y_u,fmt = 'x',label = label,**kw)
		else: plt.plot(x,y,'x',label = label)

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

		print len(lt3_t_list),len(lt4_t_list)
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
