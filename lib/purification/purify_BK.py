"""
Provides analysis functions for Barret Kok measurements

NK 2016
"""


import purify_pq as ppq
import numpy as np
import os,h5py
from analysis.lib.tools import toolbox as tb; reload(tb)
from analysis.lib.pq import pq_tools
import copy as cp
from matplotlib import pyplot as plt
import BK_analysis_params as BK_params;reload(BK_params)

#### standard parameters
w_start = 2790e-9
w_length = 50e-9
w_separation = 500e-9

save_plot_folder = r'D:\measuring\data\purification_data'


class BK_analysis(object):
	"""
	general class that stores prefiltered data and serves as analysis suite for non-local correlation measurements
	"""

	def __init__(self,name,lt3_folder,lt4_folder,**kw):
		"""
		data is in general stored in dictionaries associated with one of the two setups. each dictionary entry contains a list of np.arrays (each entry corresponding to one given timestamp)
		"""
		self.name = name

		self.key_list = ['/PQ_sync_number-1','/PQ_channel-1','/PQ_sync_time-1','/PQ_special-1','counted_awg_reps','ssro_results','tstamp','raw_data']

		self.key_list_pq =  ['/PQ_sync_number-1','/PQ_channel-1','/PQ_sync_time-1','/PQ_special-1']

		self.lt4_dict = {}
		self.lt3_dict = {}

		### initialize the data dictionaries
		for key in self.key_list:
			self.lt3_dict.update({key:[]})
			self.lt4_dict.update({key:[]})

		self.lt3_folder = lt3_folder
		self.lt4_folder = lt4_folder

	#### helper functions for plotting ####

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



	def find_tstamps_of_day(self,ts_list,day_string,analysis_folder = 'throw exception'):
		latest_t = str(int(day_string[:-1]) +1)+'_000000'
		newer_than = day_string+'_000000'

		while tb.latest_data('XX',older_than=latest_t,folder=analysis_folder,newer_than=newer_than,return_timestamp = True,raise_exc = False) != False:
			latest_t,f = tb.latest_data('XX',older_than=latest_t,folder=analysis_folder,newer_than=newer_than,return_timestamp = True,raise_exc=False)

			### debug statement that prints the full timestamp and the relevant identifier.
			#print latest_t[8:],latest_t

			### append found timestamp to list of timestamps
			ts_list.append(latest_t[8:]) 

		return ts_list

	def load_timestamps_lt3_lt4(self):
		"""
		loads the associated timestamps for BK.
		returns timestamps for LT3 and LT4.
		This should work via a static parameter dictionary in 'scripts', right now done by hand & find timestamps of day.
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

			### filter the timeharp data according to adwin events / syncs

			sync_filter = a_lt4.filter_pq_data_from_adwin_syncs() ### this syncfilter erases all data where from the PQ data where the adwin did NOT read out
			# print sync_filter
			if len(sync_filter) == 0: # empty list --> no read outs.
				# print 'file empty, skipping these time stamps:',t_lt3,t_lt4
				# print
				continue
			print 'tstamps', t_lt3,t_lt4
			### store relevant adwin results
			self.lt3_dict['ssro_results'].append(np.array(a_lt3.agrp['ssro_results'].value))
			self.lt3_dict['counted_awg_reps'].append(np.array(a_lt3.agrp['counted_awg_reps'].value))
			self.lt4_dict['ssro_results'].append(np.array(a_lt4.agrp['ssro_results'].value))
			self.lt4_dict['counted_awg_reps'].append(np.array(a_lt4.agrp['counted_awg_reps'].value))
			self.lt3_dict['tstamp'].append(t_lt3)
			self.lt4_dict['tstamp'].append(t_lt4)
			self.lt3_dict['raw_data'].append(a_lt3)
			self.lt4_dict['raw_data'].append(a_lt4)

			for key in self.key_list_pq:
				self.lt4_dict[key].append(np.array(a_lt4.pqf[key].value[sync_filter]))

			#### calculate the duty cycle for that specific file.
			# print 'lde length',a_lt3.joint_grp.attrs['LDE_element_length']
			# print'first and last time',a_lt4.pqf['/PQ_time-1'].value[0],a_lt4.pqf['/PQ_time-1'][-1]
			# print 'last elapsed time in sequence vs total elapsed time',  

			time_in_LDE_sequence = a_lt3.joint_grp.attrs['LDE_element_length']*self.lt4_dict['/PQ_sync_number-1'][-1][-1]
			total_elapsed_time = (a_lt4.pqf['/PQ_time-1'].value[-1]-a_lt4.pqf['/PQ_time-1'][0])*1e-12

			print 'file no ', i+1 , ' with duty cycle of', round(100*time_in_LDE_sequence/total_elapsed_time,1), ' %'
			i += 1

		# print self.lt4_dict['/PQ_sync_number-1'][0]

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

	def check_tail_w1_w2(self,st_start = 2000e3,p_sep = 500e3,st_len = 50e3):
		"""
		goes through the timestamp list and selects only files according to the applied sync filter: returns the measured tail in each window
		"""

		i = 0
		tails_w1,tails_w2 = [],[]

		for a in self.lt4_dict['raw_data']:
					
			### analysis for channel 0	&& window 1
			w1_ch0 = self.get_total_number_of_clicks_in_window(a,0,st_start,st_len)
			w1_ch1 = self.get_total_number_of_clicks_in_window(a,1,st_start,st_len)
			w2_ch0 = self.get_total_number_of_clicks_in_window(a,0,st_start+p_sep,st_len)
			w2_ch1 = self.get_total_number_of_clicks_in_window(a,1,st_start+p_sep,st_len)
			last_sync = a.pqf['/PQ_sync_number-1'][-1]

			tail_w1 = round(1e4*(w1_ch0+w1_ch1)/last_sync,2)
			tail_w2 = round(1e4*(w2_ch0+w2_ch1)/last_sync,2)

			# print 'tail in w1 / w2 (1e-4)    ', tail_w1, ' / ', tail_w2
			tails_w1.append(tail_w1);tails_w2.append(tail_w2)

		f,ax = self.create_plot(ylabel = 'Tail (10e-4)', title = 'Tail counts vs run number; w1_start ' +  str(round(st_start/1e3,0)))
		self.plot_data(range(len(tails_w1)),tails_w1,label = 'w1')
		self.plot_data(range(len(tails_w2)),tails_w2,label = 'w2')
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

	def apply_temporal_filters_to_prefiltered_data(self,st_start = None,st_len = None,st_len_w2 = None,p_sep =None,verbose = True):
		self.st_fltr_w1 = []
		self.st_fltr_w2 = []
		no_w1 = 0
		no_w2 = 0

		### one can also apply manual filters if one wants to deviate from the prescribed parameter dictionary
		if st_start == None:
			st_start = BK_params.BK_params['st_start']
		if st_len == None:
			st_len = BK_params.BK_params['st_len']
		if st_len_w2 == None:
			st_len_w2 = BK_params.BK_params['st_len_w2']
		if p_sep == None:
			p_sep = BK_params.BK_params['p_sep']


		for st_filtered,sp_filtered in zip(self.lt4_dict['/PQ_sync_time-1'],self.lt4_dict['/PQ_special-1']):
			st_fltr_w1 = (st_filtered > st_start)  & (st_filtered < (st_start  + st_len)) & (sp_filtered == 0)
			st_fltr_w2 = (st_filtered > st_start + p_sep)  & (st_filtered < (st_start + p_sep + st_len_w2)) & (sp_filtered == 0)
			self.st_fltr_w1.append(st_fltr_w1)
			self.st_fltr_w2.append(st_fltr_w2)

			no_w1 += np.sum(st_fltr_w1)
			no_w2 += np.sum(st_fltr_w2)
			if verbose:
				print 'number of filtered detection events in each window w1 / w2: ', no_w1, ' / ', no_w2
		return

	def attach_state_filtered_syncs(self,verbose = True):
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

			if verbose:
				print 'Run ', i+1
				print 'total events w1 / w2', np.sum(st_fltr_w1_ch0)+np.sum(st_fltr_w1_ch1),'/',np.sum(st_fltr_w2_ch0)+np.sum(st_fltr_w2_ch1)
				print 'psi_minus events', np.sum(st_fltr_psi_minus[1:-1])
				print 'psi_plus events', np.sum(st_fltr_psi_plus[1:-1])
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

		loop_array = zip(self.lt3_dict['tstamp'],self.lt4_dict['tstamp'],self.HH_sync_psi_plus,self.HH_sync_psi_minus,self.lt3_dict['counted_awg_reps'],self.lt4_dict['counted_awg_reps'],self.lt3_dict['ssro_results'],self.lt4_dict['ssro_results'])
		for t_lt3,t_lt4,HH_s_psi_p,HH_s_psi_m,adwin_syncs_lt3,adwin_syncs_lt4,adwin_ro_lt3,adwin_ro_lt4 in loop_array:


			fltr_plus_lt3  = self.filter_adwin_data_from_pq_syncs(HH_s_psi_p,adwin_syncs_lt3)
			fltr_minus_lt3 = self.filter_adwin_data_from_pq_syncs(HH_s_psi_m,adwin_syncs_lt3)
			fltr_plus_lt4  = self.filter_adwin_data_from_pq_syncs(HH_s_psi_p,adwin_syncs_lt4)
			fltr_minus_lt4 = self.filter_adwin_data_from_pq_syncs(HH_s_psi_m,adwin_syncs_lt4)

			# print t_lt3,t_lt4,adwin_syncs_lt3,fltr_plus_lt3,HH_s_psi_p
			# print t_lt3,adwin_ro_lt3,adwin_ro_lt4

			self.RO_data_LT3_plus.append(adwin_ro_lt3[fltr_plus_lt3])
			self.RO_data_LT4_plus.append(adwin_ro_lt4[fltr_plus_lt4])
			self.RO_data_LT3_minus.append(adwin_ro_lt3[fltr_minus_lt3])
			
			self.RO_data_LT4_minus.append(adwin_ro_lt4[fltr_minus_lt4])

		# print fltr_plus_lt3
		# print fltr_plus_lt4

		if not apply_ROC:
			all_m_lt3,all_m_lt4,all_p_lt3,all_p_lt4 = np.array([]),np.array([]),np.array([]),np.array([])
			for m_lt3,m_lt4,p_lt3,p_lt4 in zip(self.RO_data_LT3_minus,self.RO_data_LT4_minus,self.RO_data_LT3_plus,self.RO_data_LT4_plus):

				if len(m_lt3) != 0 and len(m_lt4) != 0:
					if verbose:
						print 'p_correlated for psi_minus', float(np.sum(np.equal(m_lt3,m_lt4)))/len(m_lt3)
						print 'p_correlated for psi_plus', float(np.sum(np.equal(p_lt3,p_lt4)))/len(m_lt4)

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

		else: 
			print 'You need to write a routine that takes the read-out correction into account'

		if return_value:
			return m_correlations,p_correlations

	def sweep_filter_parameter_vs_correlations(self,parameter_name,parameter_range):


		### initialize results lists
		no_of_minus_events,no_of_plus_events = [],[]
		minus_correlation,plus_correlation	 = [],[]
		minus_correlation_u,plus_correlation_u = [],[]

		for x in parameter_range:

			BK_params.BK_params[parameter_name] = x ### commence sweep

			self.apply_temporal_filters_to_prefiltered_data(verbose = False)
			self.attach_state_filtered_syncs(verbose = False)
			psi_m_corrs, psi_p_corrs = self.correlate_RO_results(verbose=False,return_value = True)

			no_of_minus_events.append(np.sum(psi_m_corrs))
			no_of_plus_events.append(np.sum(psi_p_corrs))

			no_anti_correlations_m = float(psi_m_corrs[1]+psi_m_corrs[2])
			no_correlations_m = float(psi_m_corrs[0]+psi_m_corrs[3])
			no_anti_correlations_p = float(psi_p_corrs[1]+psi_p_corrs[2])
			no_correlations_p = float(psi_p_corrs[0]+psi_p_corrs[3])

			minus_correlation.append(float(no_anti_correlations_m)/np.sum(psi_m_corrs))
			plus_correlation.append(float(no_correlations_p)/np.sum(psi_p_corrs))

			### error bars are based on poissonian statistics for correlated events vs. uncorrelated
			minus_correlation_u.append(np.sqrt((no_anti_correlations_m*(no_correlations_m**2)+no_correlations_m*(no_anti_correlations_m**2)))/((no_correlations_m + no_anti_correlations_m)**2))
			plus_correlation_u.append(np.sqrt((no_anti_correlations_p*(no_correlations_p**2)+no_correlations_p*(no_anti_correlations_p**2)))/((no_correlations_p + no_anti_correlations_p)**2))

		### commence plotting
		self.create_plot(title = 'Number of events within filter', xlabel = parameter_name + ' (ns)',ylabel = 'Occurences')
		self.plot_data(parameter_range/1e3,no_of_minus_events,label = '-')
		self.plot_data(parameter_range/1e3,no_of_plus_events,label = '+')
		plt.legend()
		self.save_and_close_plot()


		### should really have an error calculation for the minus / plus correlation
		self.create_plot(title = 'Fraction of correct correlations', xlabel = parameter_name + ' (ns)', ylabel = 'p_right_correlation')
		self.plot_data(parameter_range/1e3,minus_correlation,y_u = minus_correlation_u,label = '-')
		self.plot_data(parameter_range/1e3,plus_correlation,y_u = plus_correlation_u,label = '+')
		plt.legend(loc = 2)
		self.save_and_close_plot()

		reload(BK_params) ### to negate the sweep changes

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
