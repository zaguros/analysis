"""
Provides analysis functions for Barret Kok measurements

NK 2016
"""


import purify_pq as ppq
import numpy as np
import os,h5py
from analysis.lib.tools import toolbox as tb
from analysis.lib.pq import pq_tools,pq_plots
import copy as cp
from matplotlib import pyplot as plt
from analysis.lib.lde import sscorr; reload(sscorr)
from analysis.lib.m2.ssro import ssro
import BK_analysis_params as BK_params;reload(BK_params)

#### standard parameters

save_plot_folder = r'D:\measuring\data\purification_plots'


class BK_analysis(object):
	"""
	general class that stores prefiltered data and serves as analysis suite for non-local correlation measurements
	"""

	def __init__(self,name,lt3_folder,lt4_folder,ROC_lt3_tstamp,ROC_lt4_tstamp,**kw):
		"""
		data is in general stored in dictionaries associated with one of the two setups. each dictionary entry contains a list of np.arrays (each entry corresponding to one given timestamp)
		"""
		self.name = name

		self.key_list = ['/PQ_sync_number-1','/PQ_channel-1','/PQ_sync_time-1','/PQ_special-1','counted_awg_reps','ssro_results','tstamp','raw_data','CR_after']

		self.key_list_pq =  ['/PQ_sync_number-1','/PQ_channel-1','/PQ_sync_time-1','/PQ_special-1']

		self.lt4_dict = {}
		self.lt3_dict = {}

		### initialize the data dictionaries
		for key in self.key_list:
			self.lt3_dict.update({key:[]})
			self.lt4_dict.update({key:[]})

		self.lt3_folder = lt3_folder
		self.lt4_folder = lt4_folder

		self.ROC_lt3_tstamp = ROC_lt3_tstamp
		self.ROC_lt4_tstamp = ROC_lt4_tstamp

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

	def tstamps_for_both_setups(self,day_string,newest_tstamp = '235959'):
		"""
		takes a date as input and scans lt3 and lt4 for appropriate timestamps
		will throw an error if both setups have run the experiment an unequal amount of times!
		--> then you have to clean up the data folder of the setup with a longer list of timestamps
		input: day_string, e.g. '20160607'
		output: lt3_t_list,lt4_t_list
		"""

		lt3_t_list = self.find_tstamps_of_day([],day_string,analysis_folder = self.lt3_folder ,newest_tstamp = newest_tstamp)
		lt4_t_list = self.find_tstamps_of_day([],day_string,analysis_folder = self.lt4_folder, newest_tstamp = newest_tstamp)

		return self.verify_tstamp_lists(lt3_t_list,lt4_t_list,day_string)


	def find_tstamps_of_day(self,ts_list,day_string,analysis_folder = 'throw exception',newest_tstamp = '235959'):
		latest_t = day_string + newest_tstamp # where in the day do you want to begin? 235959 mean: take the whole day
		newer_than = day_string+'_000000'
		while tb.latest_data('XX',older_than=latest_t,folder=analysis_folder,newer_than=newer_than,return_timestamp = True,raise_exc = False) != False:
			latest_t,f = tb.latest_data('XX',older_than=latest_t,folder=analysis_folder,newer_than=newer_than,return_timestamp = True,raise_exc=False)
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


	def load_raw_data(self,lt3_timestamps,lt4_timestamps):
		"""
		this script takes a list of timestamps for both setups and prefilters them according to adwin filters
		creates a list of arrays associated with each time stamp for adwin_ssro,syncs, time, special, channel, marker and attaches is it to the data object for further processing.
		"""
		length = len(lt3_timestamps)
		print 'loading the data, total number of files ', length

		i = 0

		for t_lt3,t_lt4 in zip(lt3_timestamps,lt4_timestamps):
			
			# print 'tstamps', t_lt3,t_lt4
			# print tb.latest_data(t_lt4,folder = self.lt4_folder)
			
			a_lt3 = ppq.purifyPQAnalysis(tb.latest_data(t_lt3,folder = self.lt3_folder),hdf5_mode ='r')
			a_lt4 = ppq.purifyPQAnalysis(tb.latest_data(t_lt4,folder = self.lt4_folder),hdf5_mode ='r')
			# print a_lt3.agrp
			### filter the timeharp data according to adwin events / syncs

			sync_filter = a_lt4.filter_pq_data_from_adwin_syncs() ### this syncfilter erases all data where from the PQ data where the adwin did NOT read out
			# print sync_filter
			if len(sync_filter) == 0: # empty list --> no read outs.
				# print 'file empty, skipping these time stamps:',t_lt3,t_lt4
				# print
				continue


			### store relevant adwin results
			self.lt3_dict['ssro_results'].append(np.array(a_lt3.agrp['ssro_results'].value))
			self.lt3_dict['counted_awg_reps'].append(np.array(a_lt3.agrp['counted_awg_reps'].value))
			self.lt3_dict['CR_after'].append(np.array(a_lt3.agrp['CR_after'].value))
			self.lt4_dict['ssro_results'].append(np.array(a_lt4.agrp['ssro_results'].value))
			self.lt4_dict['counted_awg_reps'].append(np.array(a_lt4.agrp['counted_awg_reps'].value))
			self.lt4_dict['CR_after'].append(np.array(a_lt4.agrp['CR_after'].value))
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
			# print self.lt4_dict['/PQ_sync_number-1'][-1][0]
			time_in_LDE_sequence = a_lt3.joint_grp.attrs['LDE_element_length']*a_lt4.pqf['/PQ_sync_number-1'][-1]
			total_elapsed_time = (a_lt4.pqf['/PQ_time-1'].value[-1]-a_lt4.pqf['/PQ_time-1'][0])*1e-12
			no_of_syncs = a_lt4.pqf['/PQ_sync_number-1'][-1]
			print i+1 , ' dc: ', round(100*time_in_LDE_sequence/total_elapsed_time,1), ' %   syncs ', no_of_syncs
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

	
	def make_single_hist(self,a,channel,st_start,st_len,bin_size):

		is_ph = pq_tools.get_photons(a.pqf)[channel]
		bins = np.arange(st_start-.5,st_start+st_len,bin_size)
		y,x=np.histogram(a.pqf['/PQ_sync_time-1'].value[np.where(is_ph)], bins=bins)
		x=x[:-1]

		return x,y

	def plot_total_histogram(self,channel,st_start,st_len,bin_size,log_plot = False):

		i = 0
		for a in self.lt4_dict['raw_data']:

			if i == 0:
				x,hist = self.make_single_hist(a,channel,st_start,st_len,bin_size)
				i=1
			else:
				hist = hist + self.make_single_hist(a,channel,st_start,st_len,bin_size)[1]

		fig = plt.figure(figsize=(12,6))
		ax = fig.add_subplot(111)

		ax.plot(x/1000., hist, drawstyle='steps-post')
		### for regular plotting
		ax.set_ylim(bottom=0,top=max(hist)+100)

		#log plots look different
		if log_plot:
		    ax.set_yscale('log')
		    ax.set_ylim(bottom=0,top=max(hist)+0.5*max(hist))
		ax.set_xlabel('time (ns)')
		ax.set_ylabel('events')
		
		ax.set_xlim(min(x/1000.), max(x/1000.))
		ax.grid(True)

		plt.show()
		plt.close('all')

	def apply_temporal_filters_to_prefiltered_data(self,st_start 	= None,
														st_len 		= None,
														ch1_offset 	= None,
														st_len_w2 	= None,
														p_sep 		= None,
														dt_max 		= None,
														verbose 	= True):
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

		if dt_max == None:
			dt_max = BK_params.BK_params['dt_max']

		if ch1_offset == None:
			ch1_offset = BK_params.BK_params['ch1_offset']


		for st_filtered,sp_filtered,ch_filtered in zip(self.lt4_dict['/PQ_sync_time-1'],self.lt4_dict['/PQ_special-1'],self.lt4_dict['/PQ_channel-1']):
			### compare dt and include p_sep here
			### this copying should be avoided but is the only way that we do not have memory effects when alerting st_filtered with the ch1 offset.
			st_filtered = cp.deepcopy(st_filtered)

			### shift channel one by a certain time margin.
			ch1s = np.argwhere((ch_filtered == 1) & (sp_filtered == 0))
			st_filtered[ch1s] = st_filtered[ch1s] + ch1_offset

			dt  = np.abs(st_filtered - np.insert(st_filtered,0,0)[:-1] -p_sep)# - p_sep

			st_fltr_w1 = (st_filtered > st_start)  & (st_filtered < (st_start  + st_len)) & (sp_filtered == 0) 
			st_fltr_w2 = (st_filtered > st_start + p_sep)  & (st_filtered < (st_start + p_sep + st_len_w2)) & (sp_filtered == 0) & (dt < dt_max)
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
			#### the filter w1 is shifted onto the filter of w2 by inserting a boolean False at the beginning (clicks must be consecutive across windows)

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

	def filter_CR_after(self,threshold = None,verbose = False):
		### checks if there was low cr_after after
		### returns a filter for the specific sync numbers.
		if threshold == None:
			threshold = BK_params.BK_params['CR_after_threshold']
		### acquire filter based on CR after
		loop_array = zip(range(len(self.lt3_dict['CR_after'])),self.lt3_dict['CR_after'],self.lt4_dict['CR_after'],self.lt4_dict['counted_awg_reps'])
			
		for i,CR_lt3,CR_lt4,adwin_syncs in loop_array:

			CR_fltr = (CR_lt3 > threshold) & (CR_lt4 > threshold)
			adwin_syncs_filtered = adwin_syncs[CR_fltr]

			pq_filtered = self.lt3_dict['raw_data'][i].filter_pq_data_from_adwin_syncs(adwin_syncs = adwin_syncs_filtered,pq_syncs = self.lt4_dict['/PQ_sync_number-1'][i])
			# print 'PQ_after',pq_filtered, self.lt3_dict['tstamp'][i],self.lt4_dict['tstamp'][i]


			### neglect the data entirely
			if len(pq_filtered) == 0:
				if verbose:
					print 'empty array!',self.lt3_dict['tstamp'][i],self.lt4_dict['tstamp'][i] 
				for key in self.key_list_pq:
					self.lt4_dict[key][i] = pq_filtered
				

			### refilter the PQ data in each array
			else:
				for key in self.key_list_pq:
					self.lt4_dict[key][i] = cp.deepcopy(np.array(self.lt4_dict[key][i][pq_filtered]))




	def correlate_RO_results(self,apply_ROC = False, verbose = True,return_value = False):

		self.RO_data_LT3_plus = []
		self.RO_data_LT3_minus = []
		self.RO_data_LT4_plus = []
		self.RO_data_LT4_minus = []

		loop_array = zip(self.lt3_dict['tstamp'],self.lt4_dict['tstamp'],self.HH_sync_psi_plus,self.HH_sync_psi_minus, \
						 self.lt3_dict['counted_awg_reps'],self.lt4_dict['counted_awg_reps'],self.lt3_dict['ssro_results'],
						 self.lt4_dict['ssro_results'],self.lt3_dict['raw_data'],self.lt4_dict['raw_data'])

		for t_lt3,t_lt4,HH_s_psi_p,HH_s_psi_m,adwin_syncs_lt3,adwin_syncs_lt4,adwin_ro_lt3,adwin_ro_lt4,a_lt3,a_lt4 in loop_array:


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

		all_m_lt3,all_m_lt4,all_p_lt3,all_p_lt4 = np.array([]),np.array([]),np.array([]),np.array([])
		for m_lt3,m_lt4,p_lt3,p_lt4 in zip(self.RO_data_LT3_minus,self.RO_data_LT4_minus,self.RO_data_LT3_plus,self.RO_data_LT4_plus):

			# if verbose:
			# 	if len(m_lt3) != 0:
			# 		print 'p_correlated for psi_minus', float(np.sum(np.equal(m_lt3,m_lt4)))/len(m_lt3)
			# 	if len(p_lt4) != 0:
			# 		print 'p_correlated for psi_plus', float(np.sum(np.equal(p_lt3,p_lt4)))/len(p_lt4)

			all_m_lt3 = np.append(all_m_lt3,m_lt3)
			all_m_lt4 = np.append(all_m_lt4,m_lt4)
			all_p_lt3 = np.append(all_p_lt3,p_lt3)
			all_p_lt4 = np.append(all_p_lt4,p_lt4)



		#### print correlation matrix for RO results
		### 


		### get overall events psi minus:
		m_correlations = [0,0,0,0]

		### ROC correction should happen before this step to account for different RO durations (which we should never have in the first place.)

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
			# def ssro_correct_twoqubit_state_photon_numbers(correlations, F0a, F0b, F1a, F1b,
			#      return_error_bars=False, dF0a=0.01, dF0b=0.01, dF1a=0.01, dF1b=0.01, 
			#      verbose = True):
			### get ssro_ROC for LT3 --> corresponds to setup B
			F0_LT3,F1_LT3 = self.find_RO_fidelities(self.ROC_lt3_tstamp,a_lt3,folder = self.lt3_folder)
			### get ssro_ROC for LT4 --> corresponds to setup A
			F0_LT4,F1_LT4 = self.find_RO_fidelities(self.ROC_lt4_tstamp,a_lt4,folder = self.lt4_folder)

			### apply ROC to the results
			# m_correlations = m_correlations.tolist()
			corrected_psi_minus = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(m_correlations[::-1]),F0_LT4,F0_LT3,F1_LT4,F1_LT3)
			corrected_psi_plus  = sscorr.ssro_correct_twoqubit_state_photon_numbers(np.array(p_correlations[::-1]),F0_LT4,F0_LT3,F1_LT4,F1_LT3)

		if return_value:
			return m_correlations,p_correlations

	def find_RO_fidelities(self,timestamp,raw_data,folder = ''):

		ssro_folder = tb.data_from_time(timestamp,folder = folder) 

		F_0,u_F0,F_1,u_F1 = ssro.get_SSRO_calibration(ssro_folder,raw_data.g.attrs['E_RO_durations'][0]) #we onlky have one RO time in E_RO_durations

		return F_0,F_1

	def sweep_filter_parameter_vs_correlations(self,parameter_name,parameter_range):


		### initialize results lists
		no_of_minus_events,no_of_plus_events = [],[]
		minus_correlation,plus_correlation	 = [],[]
		minus_correlation_u,plus_correlation_u = [],[]

		for x in parameter_range:

			BK_params.BK_params[parameter_name] = x ### commence sweep

			# self.filter_CR_after(threshold = None,verbose = False)

			self.apply_temporal_filters_to_prefiltered_data(verbose = False)
			self.attach_state_filtered_syncs(verbose = False)
			psi_m_corrs, psi_p_corrs = self.correlate_RO_results(verbose=False,return_value = True)

			no_of_minus_events.append(np.sum(psi_m_corrs))
			no_of_plus_events.append(np.sum(psi_p_corrs))

			no_anti_correlations_m = float(psi_m_corrs[1]+psi_m_corrs[2])
			no_correlations_m = float(psi_m_corrs[0]+psi_m_corrs[3])
			no_anti_correlations_p = float(psi_p_corrs[1]+psi_p_corrs[2])
			no_correlations_p = float(psi_p_corrs[0]+psi_p_corrs[3])


			if np.sum(psi_m_corrs) == 0:
				minus_correlation.append(0)
				minus_correlation_u.append(0)
			else:
				minus_correlation_u.append(np.sqrt((no_anti_correlations_m*(no_correlations_m**2)+no_correlations_m*(no_anti_correlations_m**2)))/((no_correlations_m + no_anti_correlations_m)**2))
				minus_correlation.append(float(no_anti_correlations_m)/np.sum(psi_m_corrs))
			

			if np.sum(psi_p_corrs) == 0:
				plus_correlation.append(0)
				plus_correlation_u.append(0)
			else:
				plus_correlation.append(float(no_correlations_p)/np.sum(psi_p_corrs))
				plus_correlation_u.append(np.sqrt((no_anti_correlations_p*(no_correlations_p**2)+no_correlations_p*(no_anti_correlations_p**2)))/((no_correlations_p + no_anti_correlations_p)**2))

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
