"""
Provides analysis functions for Barret Kok measurements

NK 2016
"""


import purify_pq as ppq
import numpy as np
import os,h5py
from analysis.lib.tools import toolbox as tb; reload(tb)
import copy as cp

#### standard parameters
w_start = 2790e-9
w_length = 50e-9
w_separation = 500e-9



class BK_analysis(object):
	"""
	general class that stores prefiltered data and serves as analysis suite for non-local correlation measurements
	"""

	def __init__(self,name,lt3_folder,lt4_folder,**kw):
		"""
		data is in general stored in dictionaries associated with one of the two setups. each dictionary entry contains a list of np.arrays (each entry corresponding to one given timestamp)
		"""
		self.name = name

		self.key_list = ['/PQ_sync_number-1','/PQ_channel-1','/PQ_sync_time-1','/PQ_special-1','counted_awg_reps','ssro_results']

		self.key_list_pq =  ['/PQ_sync_number-1','/PQ_channel-1','/PQ_sync_time-1','/PQ_special-1']

		self.lt4_dict = {}
		self.lt3_dict = {}

		### initialize the data dictionaries
		for key in self.key_list:
			self.lt3_dict.update({key:[]})
			self.lt4_dict.update({key:[]})

		self.lt3_folder = lt3_folder
		self.lt4_folder = lt4_folder

	def load_timestamps_lt3_lt4(self):
		"""
		loads the associated timestamps for BK.
		returns timestamps for LT3 and LT4.
		This should work via a static parameter dictionary in 'scripts', right now done by hand
		"""
		pass

	def load_raw_data(self,lt3_timestamps,lt4_timestamps):
		"""
		this script takes a list of timestamps for both setups and prefilters them according to adwin filters
		creates a list of arrays associated with each time stamp for adwin_ssro,syncs, time, special, channel, marker and attaches is it to the data object for further processing.
		"""
		length = len(lt3_timestamps)
		print 'loading the data, total number of files ', length

		i = 0

		for t_lt3,t_lt4 in zip(lt3_timestamps,lt4_timestamps):
			
			# print self.lt3_folder
			a_lt3 = ppq.purify_pq(tb.latest_data(t_lt3,folder = self.lt3_folder),hdf5_mode ='r')
			a_lt4 = ppq.purify_pq(tb.latest_data(t_lt4,folder = self.lt4_folder),hdf5_mode ='r')


			### store relevant adwin results
			self.lt3_dict['ssro_results'].append(np.array(a_lt3.agrp['ssro_results'].value))
			self.lt3_dict['counted_awg_reps'].append(np.array(a_lt3.agrp['counted_awg_reps'].value))
			self.lt4_dict['ssro_results'].append(np.array(a_lt4.agrp['ssro_results'].value))
			self.lt4_dict['counted_awg_reps'].append(np.array(a_lt4.agrp['counted_awg_reps'].value))


			### filter the timeharp data according to adwin events / syncs

			sync_filter = a_lt4.filter_pq_data_from_adwin_syncs() ### this syncfilter erases all data where from the PQ data where the adwin did NOT read out
			for key in self.key_list_pq:
				self.lt4_dict[key].append(np.array(a_lt4.pqf[key].value[sync_filter]))

			print i+1
			i += 1

		# print self.lt4_dict['/PQ_sync_number-1'][0]

	def apply_temporal_filters_to_raw_data(self,st_start = w_start,st_len = w_length,st_len_w2 = w_length,p_sep = w_separation):
		self.st_fltr_w1 = []
		self.st_fltr_w2 = []
		no_w1 = 0
		no_w2 = 0

		for st_filtered,sp_filtered in zip(self.lt4_dict['/PQ_sync_time-1'],self.lt4_dict['/PQ_special-1']):
			st_fltr_w1 = (st_filtered > st_start)  & (st_filtered < (st_start  + st_len)) & (sp_filtered == 0)
			st_fltr_w2 = (st_filtered > st_start + p_sep)  & (st_filtered < (st_start + p_sep + st_len_w2)) & (sp_filtered == 0)
			self.st_fltr_w1.append(st_fltr_w1)
			self.st_fltr_w2.append(st_fltr_w2)

			no_w1 += np.sum(st_fltr_w1)
			no_w2 += np.sum(st_fltr_w2)

		print 'number of filtered detection events in each window w1 / w2: ', no_w1, ' / ', no_w2


	def attach_state_filtered_syncs(self):
		"""
		checks for the signatures of psi_minus or psi_plus and returns a list of numpy arrays where each numpy array corresponds the correpsonding sync number of the event
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

		for fltr_w1,fltr_w2,channels,HH_sync in zip(self.st_fltr_w1,self.st_fltr_w2,self.lt4_dict['/PQ_channel-1'],self.lt4_dict['/PQ_sync_number-1']):

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
			#### the filter w1 is shifted onto the filter of w2 by inserting a boolean False at the beginning (clicks must be consecutive accross windows)

			st_fltr_w1_ch1 = np.insert(st_fltr_w1_ch1,0,False); st_fltr_w1_ch0 = np.insert(st_fltr_w1_ch0,0,False)

			st_fltr_psi_plus = np.logical_or(np.logical_and(st_fltr_w1_ch1[:-1],st_fltr_w2_ch1),np.logical_and(st_fltr_w1_ch0[:-1],st_fltr_w2_ch0))
			st_fltr_psi_minus = np.logical_or(np.logical_and(st_fltr_w1_ch0[:-1],st_fltr_w2_ch1),np.logical_and(st_fltr_w1_ch1[:-1],st_fltr_w2_ch0))

			self.st_fltr_psi_plus.append(st_fltr_psi_plus)
			self.st_fltr_psi_minus.append(st_fltr_psi_minus)
			print 'Run ', i+1
			print 'total events w1 / w2', np.sum(st_fltr_w1_ch0)+np.sum(st_fltr_w1_ch1),'/',np.sum(st_fltr_w2_ch0)+np.sum(st_fltr_w2_ch1)
			print 'psi_minus events', np.sum(st_fltr_psi_minus[1:-1])
			print 'psi_plus events', np.sum(st_fltr_psi_plus[1:-1])
			print

			### the required HH syncs are later used to get the corresponding adwin ROs
			self.HH_sync_psi_plus.append(HH_sync[st_fltr_psi_plus])
			self.HH_sync_psi_minus.append(HH_sync[st_fltr_psi_minus])

			i += 1

	def correlate_RO_results(self,apply_ROC = False):

		self.RO_data_LT3_plus = []
		self.RO_data_LT3_minus = []
		self.RO_data_LT4_plus = []
		self.RO_data_LT4_minus = []

		for HH_s_psi_p,HH_s_psi_m,adwin_syncs_lt3,adwin_syncs_lt4,adwin_ro_lt3,adwin_ro_lt4 in zip(self.HH_sync_psi_plus,self.HH_sync_psi_minus,self.lt3_dict['counted_awg_reps'],self.lt4_dict['counted_awg_reps'],self.lt3_dict['ssro_results'],self.lt4_dict['ssro_results']):


			fltr_plus_lt3  = self.filter_adwin_data_from_pq_syncs(HH_s_psi_p,adwin_syncs_lt3)
			fltr_minus_lt3 = self.filter_adwin_data_from_pq_syncs(HH_s_psi_m,adwin_syncs_lt3)
			fltr_plus_lt4  = self.filter_adwin_data_from_pq_syncs(HH_s_psi_p,adwin_syncs_lt4)
			fltr_minus_lt4 = self.filter_adwin_data_from_pq_syncs(HH_s_psi_m,adwin_syncs_lt4)

			self.RO_data_LT3_plus.append(adwin_ro_lt3[fltr_plus_lt3])
			self.RO_data_LT3_minus.append(adwin_ro_lt3[fltr_minus_lt3])
			self.RO_data_LT4_plus.append(adwin_ro_lt4[fltr_plus_lt4])
			self.RO_data_LT4_minus.append(adwin_ro_lt4[fltr_minus_lt4])

		# print fltr_plus_lt3
		# print fltr_plus_lt4

		if not apply_ROC:
			all_m_lt3,all_m_lt4,all_p_lt3,all_p_lt4 = np.array([]),np.array([]),np.array([]),np.array([])
			for m_lt3,m_lt4,p_lt3,p_lt4 in zip(self.RO_data_LT3_minus,self.RO_data_LT4_minus,self.RO_data_LT3_plus,self.RO_data_LT4_plus):
				print 'p_correlated for psi_minus', float(np.sum(np.equal(m_lt3,m_lt4)))/len(m_lt3)
				print 'p_correlated for psi_plus', float(np.sum(np.equal(p_lt3,p_lt4)))/len(m_lt4)

				all_m_lt3 = np.append(all_m_lt3,m_lt3)
				all_m_lt4 = np.append(all_m_lt4,m_lt4)
				all_p_lt3 = np.append(all_p_lt3,p_lt3)
				all_p_lt4 = np.append(all_p_lt4,p_lt4)



			#### print correlation matrix for RO results
			### 

			print 'The occurence of each event after filtering'
			print

			x = ['ms0 & ms0', 'ms0 & ms1', 'ms1 & ms0','ms1 & ms1',]
			row_format ="{:>12}" * (len(x) + 1)
			headline_format = "{:>15}"+"{:>12}" * len(x)
			print headline_format.format("", *x)

			### get overall events psi minus:
			m_correlations = [0,0,0,0]
			m_correlations[0] = np.sum(np.equal(all_m_lt3[all_m_lt3 == 1],all_m_lt4[all_m_lt3 == 1]))
			m_correlations[1] = np.sum(np.not_equal(all_m_lt3[all_m_lt3 == 1],all_m_lt4[all_m_lt3 == 1]))
			m_correlations[2] = np.sum(np.not_equal(all_m_lt3[all_m_lt3 == 0],all_m_lt4[all_m_lt3 == 0]))
			m_correlations[3] = np.sum(np.equal(all_m_lt3[all_m_lt3 == 0],all_m_lt4[all_m_lt3 == 0]))
			# print m_correlations

			

			for state, row in zip(['psi_minus'], [m_correlations]):
				print "-"*(12*4+15)
				print row_format.format(state+' |', *row)

			print

			# print headline_format.format("", *x)
			p_correlations = [0,0,0,0]
			p_correlations[0] = np.sum(np.equal(all_p_lt3[all_p_lt3 == 1],all_p_lt4[all_p_lt3 == 1]))
			p_correlations[1] = np.sum(np.not_equal(all_p_lt3[all_p_lt3 == 1],all_p_lt4[all_p_lt3 == 1]))
			p_correlations[2] = np.sum(np.not_equal(all_p_lt3[all_p_lt3 == 0],all_p_lt4[all_p_lt3 == 0]))
			p_correlations[3] = np.sum(np.equal(all_p_lt3[all_p_lt3 == 0],all_p_lt4[all_p_lt3 == 0]))

			for state, row in zip(['psi_plus'], [p_correlations]):
				print "-"*(12*4+15)
				print row_format.format(state+' |', *row)

		else: 
			print 'You need to write a routine that takes the read-out correction into account'



	def filter_adwin_data_from_pq_syncs(self,filtered_sn,counted_awg_reps):
		"""
		takes the filtered pq syncs as input and returns a boolean array.
		This array serves as filter for the adwin RO results

		TODO: generalize for arbitrary PQ data size
		"""

		# print 'elen', len(filtered_sn)
		insert_pos = np.searchsorted(counted_awg_reps,filtered_sn)
		#insert_pos = np.searchsorted(filtered_sn,adwin_syncs)

		return insert_pos
