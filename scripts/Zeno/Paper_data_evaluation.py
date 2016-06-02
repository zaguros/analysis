"""
This file analyses the extracted data for the Zeno paper.
Functions are executeable via an ipython notebook. 
It is assumed that the data is stored in pickle files in the folder '\Path_of_the_notebook\ZenData' . 
Imports the raw data which was extracted by the file TwoQ_Zeno_Analysis_v2.py (see ZenData)
"""

import numpy as np
import os,re
import h5py
import pickle
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import ssro
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import matplotlib as mpl
import copy as cp
import Zeno_fitting_tools as Zfits
import analysis.lib.Zeno.TwoQ_Zeno_Analysis_v2 as ZAnalysis
reload(Zfits)
reload(ZAnalysis)
reload(common)
reload(toolbox)



######################
######################
##					##
##		Globals		##
##					##
######################
######################

datafolder = 'ZenData\ '[:-1] ### folder of all the pickle/data files.

timetrace_keylist = ['evotime','fid','fid_u'] ### used in most pickle dictionaries.
datapoint = '.'
### Plot styling goes here. 

figscaling = 0.66*0.72
mpl.rcParams['font.sans-serif'] = 'Lucida Grande'
mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['mathtext.rm'] = 'sans'
# mpl.rcParams['mathtext.default'] = 'rm'
# mpl.rcParams['lines.solid_joinstyle'] = 
# mpl.rcParams['errorbar.capsize'] = 3
# font = {'family' : 'sans-serif',
#     'weight' : 'normal',
#     'size'   : 30}
GlobalErrorCapsize = 0
fit_lw = 2 *figscaling
bar_width = 0.5

color_list = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','b','g','r','m','y','black','brown','c','orange','0.75','0.0','0.5']
color_list = ['#e41a1c','#377eb8','#FAA43A','#4D4D4D','#4daf4a','','#4daf4a']


######################
######################
##					##
##		Helper		##
##		Functions	##
##					##
######################
######################

def scale_figure(figscaling):
	## writes a bunch of parameters to matplotlib config.
	## input: effective figure scaling factor
	mpl.rcParams['axes.linewidth'] = 1.8*figscaling
	mpl.rcParams['xtick.major.width'] = 1.8*figscaling
	mpl.rcParams['ytick.major.width'] = 1.8*figscaling
	mpl.rcParams['font.size'] = (22-8*(figscaling-0.66*0.72)/(1-0.66*0.72))*figscaling
	mpl.rcParams['axes.titlesize'] = mpl.rcParams['font.size']
	mpl.rcParams['legend.fontsize'] = mpl.rcParams['font.size']
	mpl.rcParams['legend.labelspacing'] = 0.5*figscaling
	mpl.rcParams['legend.columnspacing'] = 1.5*figscaling
	mpl.rcParams['legend.handletextpad'] = 0.3*figscaling
	mpl.rcParams['legend.handlelength'] = 1.*figscaling
	mpl.rcParams['legend.borderpad'] = 0.2+0.2*(figscaling-0.66*0.72)/(1-0.66*0.72)
	mpl.rcParams['lines.markersize'] = 13.5*figscaling#7*figscaling 
	mpl.rcParams['figure.figsize'] = (6*figscaling,4*figscaling)
	mpl.rcParams['lines.markeredgewidth'] = 0.3/figscaling

scale_figure(figscaling)

def open_data(filename):
	with open( datafolder+filename, "rb" ) as f:
		datadict = pickle.load(f)
	return datadict

def normalize_to_theory():
	### not implemented yet
	pass

def find_nearest_idx(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx


class data_object(object):
	"""
	This class handles all processing of data.
	It is based on an initial pickle dictionary. 
	From there further parameters/results are derived and added as attributes.
	"""

	def __init__(self,name,data_dictionary, **kw):

		self.name = name

		self.pickle_dict = data_dictionary

		self.physical_model_fit = kw.pop('physical_model_fit', True) ### assumed to be True.
		self.pheno_physical_mode_fit = False

		self.process_fidelities = kw.pop('proc_fid', False) ### do we handle state fidelities or process fidelities
		self.process_to_avg_fidelity = kw.pop('avg_fid', False) ### converts process fidelity to avg. fidelity
		self.contrast			= kw.pop('contrast', False) ### no state fidelity, no process fidelity, but we handle expectation values!


		self.save_folder = r'D:/measuring/data/Zeno_results' ### standard folder for saving plots.

		self.fancy = True 
		self.plot_title = ''
		self.do_not_plot_list = []


		self.do_fit = False ### has to be set explicitly!

		self.fit_offset = 0.5
		self.use_fixed = True
		self.use_averaged = False

		self.state_key_endings = ['X','Y','Z','0','1']
		self.state_key_start = -1





	def acquire_fit_val(self,key,xval):
		if self.use_averaged:
			fitfunc = self.pickle_dict[key]['state_dict_averaged']['X']['fit']['fitfunc']
		else:
			fitfunc = self.pickle_dict[key]['fit']['fitfunc']


		if type(xval) == list:
			ret_list = []
			for x in xval:
				ret_list.append(fitfunc(x))

			return ret_list

		else:
			return fitfunc(xval)


	def return_fidelity_from_theory(self,key,time):
		""" takes the number of measurements, 
		the type of decay and returns the expected fidelity """



	def perform_physical_model_fit(self,key, plotting = False, ax = None):

		if type(key) != list:
			key = [key]

		if len(key) == 1:
			data = self.pickle_dict[key[0]]
		elif self.use_averaged:
			data = self.pickle_dict[key[0]]['state_dict_averaged'][key[1]]
		else:
			data = self.pickle_dict[key[0]][key[1]]



		evotime = data[timetrace_keylist[0]]
		fid = data[timetrace_keylist[1]]
		fid_u = data[timetrace_keylist[2]]


		### find the position in the data dictionary
		i = self.find_index(key[0])



		### fit everything to a physical model!
		print 'this is the keylist ', key

		if key[0] == '0' or key[0] == 'single_qubit':
			t = 8.25
			p0, fitfunc, fitfunc_str = common.fit_gauss(self.fit_offset, 0.43, 0., t)
			fixed = [0,1]
			if self.process_fidelities:
				fixed = [1]

			fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)

			if plotting:
				plot.plot_fit1d(fit_result, np.linspace(0,200,1001), ax=ax,color = color_list[i], plot_data=False,add_txt = False, lw = fit_lw)

			self.results.append(key[0] + ' : T = ' + str(round(fit_result['params'][-1],2)) + ' +- ' + str(round(fit_result['error'][-1],2)))

			if key[0] == '0':
				self.t = fit_result['params_dict']['sigma'] ### use t for the following fits with multiple measurements.
				self.A = fit_result['params_dict']['A']
				if self.process_fidelities:
					self.fit_offset = fit_result['params_dict']['a']
			# print fit_result


		else:
			if '0' not in self.pickle_dict.keys(): ### in exceptional cases there might be no 0 measurements reference
				self.t = 8.25

			if not self.process_fidelities and not self.contrast: ### logical decision for state fidelities.
				fit_result, result_string = ZAnalysis.fit_State_decay(key[0],ax,self.A*2,evotime,fid, t = self.t)

			if self.process_fidelities:
					fit_result, result_string = ZAnalysis.fit_process_decay(key[0],ax,self.A,evotime,fid,'2',offset0 = self.fit_offset, t = self.t)

			if plotting:
				plot.plot_fit1d(fit_result, np.linspace(0.0,200,1001), ax=ax, plot_data=False,color = color_list[i],add_txt = False, lw = fit_lw)
			self.results.append(key[0]+' : p = ' + result_string)

		if len(key) == 1:
			self.pickle_dict[key[0]]['fit'] = fit_result
		elif self.use_averaged:
			self.pickle_dict[key[0]]['state_dict_averaged']['fit'] = fit_result
		else:
			self.pickle_dict[key[0]][key[1]] = fit_result

		

	
	def perform_phenom_fit(self,key,plotting = False, ax = None):

		"""performs a fit of time trace data to a gaussian function"""

		if len(key) == 1:
			data = self.pickle_dict[key[0]]
		else:
			data = self.pickle_dict[key[0]][key[1]]
		evotime = data[timetrace_keylist[0]]
		fid = data[timetrace_keylist[1]]
		fid_u = data[timetrace_keylist[2]]

		ii = self.find_index(key[0])

		p0, fitfunc, fitfunc_str = common.fit_gauss(self.fit_offset, 0.43, 0., 8.25)

		fixed =[1]

		if self.use_fixed:
			if int(key[0]) % 2 ==0: ### key[0] is always the number of measurements.
				fixed = [0,1]	
			else: 
				fixed = [1]

		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)

		if plotting:
			plot.plot_fit1d(fit_result, np.linspace(0,200,201), ax=ax, plot_data=False,color = color_list[ii],add_txt = False, lw = fit_lw)

		p1 = abs(round(fit_result['params'][3-len(fixed)],3))
		p1_u = round(fit_result['error'][3-len(fixed)],3)

		self.results.append(key[0] + ' : T = ' + str(round(p1,2)) + ' +- ' + str(round(p1_u,2)))

		# print 'uncertainty ', p1_u;
		if len(key) == 1:
			self.pickle_dict[key[0]]['fit'] = fit_result
		else:
			self.pickle_dict[key[0]][key[1]]['fit'] = fit_result



	def perform_semi_pheno_model_fit(self,key, plotting = False, ax = None, xlim = None):

		if len(key) == 1:
			data = self.pickle_dict[key[0]]
		elif self.use_averaged:
			data = self.pickle_dict[key[0]]['state_dict_averaged'][key[1]]
		else:
			data = self.pickle_dict[key[0]][key[1]]

		evotime = data[timetrace_keylist[0]]
		fid = data[timetrace_keylist[1]]
		fid_u = data[timetrace_keylist[2]]


		if xlim != None:
			last_plot_point = xlim[-1]
		else:
			last_plot_point = evotime[-1]
		### find the position in the data dictionary
		i = self.find_index(key[0])

		# print 'this is the ', key

		if key[0] == '0' or key[0] == 'single_qubit' or key[0] == '00':
			t = 8.25
			p0, fitfunc, fitfunc_str = common.fit_gauss(self.fit_offset, 0.43, 0., t)
			fixed = [0,1]

			if self.process_fidelities:
				fixed = [1]

			if self.contrast:
				p0, fitfunc, fitfunc_str = common.fit_gauss(0, 0.43, 0., t)
				
			
			fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

			if plotting:
				plot.plot_fit1d(fit_result, np.linspace(0.0,200,1001), ax=ax,color = color_list[i], plot_data=False,add_txt = False, lw = fit_lw)
			if key[0] == 'single_qubit':
				label_begin = '1Q'
			else: label_begin = key[0]

			self.results.append(label_begin + ' : T = ' + str(round(abs(fit_result['params'][-1]),2)) + ' +- ' + str(round(abs(fit_result['error'][-1]),2)))

			if key[0] == '0' or key[0] == '00':
				self.t = fit_result['params_dict']['sigma'] ### use t for the following fits with multiple measurements.
				self.A = fit_result['params_dict']['A']
				if self.process_fidelities:
					self.fit_offset = fit_result['params_dict']['a']


		else:
			if '0' not in self.pickle_dict.keys(): ### in exceptional cases there might be no 0 measurements reference
				self.t = 8.25

			if not self.process_fidelities and not self.contrast: ### logical decision for state fidelities.
				fit_result, result_string = ZAnalysis.fit_pheno_State_decay(key[0],ax,self.A*2,evotime,fid, t = self.t)

			elif self.process_fidelities:
				fit_result, result_string = ZAnalysis.fit_pheno_process_decay(key[0],ax,self.A,evotime,fid,'2',offset0 = self.fit_offset, t = self.t)

			elif self.contrast:
				fit_result, result_string = ZAnalysis.fit_pheno_State_decay(key[0],ax,self.A,evotime,fid,t = self.t,contrast = self.contrast)

			if plotting:
				plot.plot_fit1d(fit_result, np.linspace(0.0,200,1001), ax=ax, plot_data=False,color = color_list[i],add_txt = False, lw = fit_lw)
			
			# print fit_result['params_dict']
			result_string = str(round(abs(fit_result['params_dict']['t']),2)) + ' +- ' + str(round(fit_result['error_dict']['t'],2))
			self.results.append(key[0]+' : T = ' + result_string)

		# print 'uncertainty ', p1_u;
		if len(key) == 1:
			self.pickle_dict[key[0]]['fit'] = fit_result
		elif self.use_averaged:
			self.pickle_dict[key[0]]['state_dict_averaged'][key[1]]['fit'] = fit_result
		else:
			self.pickle_dict[key[0]][key[1]]['fit'] = fit_result



	def get_fit_values(self):

		""" fits the data. appends results to the data dictionaries"""

		self.results = []

		for key in self.pickle_dict.keys():
			if self.physical_model_fit:
				self.perform_physical_model_fit(key)

			else: self.perform_phenom_fit(key)



	def plot_timetrace(self,save_plot 			= False, 
							add_horizontal_line = None,
							add_vertical_line 	= None,
							ylim 				= None, 
							xlim 				= None,
							legend 				= True,
							plot_states			= False,**kw):
		
		xticks = kw.pop('xticks',None)
		yticks = kw.pop('yticks',None)
		ylabel = kw.pop('ylabel',None)

		save_name = 'timetrace'


		fig=plt.figure() # parameters for plots at the beginning of the code
		ax=plt.subplot()

		i = 0

		self.results = []

		keylist = sorted(self.pickle_dict.keys())

		### the single qubit data needs to move to the background and is therefore plotted first.
		# if 'single_qubit' in keylist:
		# 	keylist[-1]=keylist[0]
		# 	keylist[0] = 'single_qubit'
		print keylist
		for key in keylist:

			# print key
			# print data.keys()
			data = self.pickle_dict[key]


			## define legend entry for the data set
			legendlabel = key
			while legendlabel[0] == '0' and len(legendlabel)!=1:
				legendlabel = legendlabel[-1:]


			### can be uncommented to add a N = .. to the legend label.

			legendlabel = 'N = ' +legendlabel
			if 'single_qubit' in legendlabel:
					legendlabel = 'Spin 1'#r'C$_1$'
					# guarantess that singlequbit data is plotted in the back
					zorder = 3

			else:
				zorder = 5


			# converts process fidelities to avg fidelities	
			try:
				if self.process_to_avg_fidelity:
					data[timetrace_keylist[1]] = (2*data[timetrace_keylist[1]]+1)/3.
					data[timetrace_keylist[2]] = 2*data[timetrace_keylist[2]]/3.
				evotime = data[timetrace_keylist[0]]
				fid = data[timetrace_keylist[1]]
				fid_u = data[timetrace_keylist[2]]

				
			except KeyError: # if plot states is false then the code above is not used.
				'didnt do conversion (Line 392). No problem since plot_states = ', plot_states

			# reset do_plot and check for key in the do not plot list
			do_plot = True
			if key in self.do_not_plot_list:
					do_plot = False

			

			if plot_states:

				if self.use_averaged:
					data = data['state_dict_averaged']
					# print data['X']['fid']
					# print data['X']['fid_u']
				for key2 , state_data in data.iteritems():
					if key2[self.state_key_start:] in self.state_key_endings:
						# print key2
						# print state_data.keys()

						### get state fidelities!
						key0 = key
						key = key2 ## dirty overwrite for the remaining procedure
						# print state_data.keys()
						evotime = np.array(state_data[timetrace_keylist[0]])
						fid = np.array(state_data[timetrace_keylist[1]])
						fid_u = np.array(state_data[timetrace_keylist[2]])
						# print evotime
						# print fid
					else:
						continue ### next loop iteration		


					

					if self.do_fit:

						if self.physical_model_fit:
							
							self.perform_physical_model_fit([key0,key],plotting = do_plot, ax=ax)

						elif self.pheno_physical_mode_fit:
							print 'key and key 2', key0, key2
							

							self.perform_semi_pheno_model_fit([key0,key],plotting = do_plot,ax=ax, xlim = xlim)

						else:

							self.perform_phenom_fit([key0,key],plotting = do_plot, ax=ax)
						if do_plot:
							plt.errorbar(np.sort(evotime),fid[np.argsort(evotime)],fid_u[np.argsort(evotime)],fmt=datapoint,label=legendlabel,color=color_list[i],capsize = GlobalErrorCapsize,zorder=zorder)

					else:
						if do_plot:
							plt.errorbar(np.sort(evotime),fid[np.argsort(evotime)],fid_u[np.argsort(evotime)],marker=datapoint,color=color_list[i],label = legendlabel,capsize = GlobalErrorCapsize,zorder=zorder)

					i +=1	


			else:				

				if self.do_fit:


					if self.physical_model_fit:
						self.perform_physical_model_fit(key,plotting = do_plot, ax=ax)


					elif self.pheno_physical_mode_fit:

						self.perform_semi_pheno_model_fit([key],plotting = do_plot,ax=ax)

					else:
						self.perform_phenom_fit([key],plotting = do_plot, ax=ax)
				
					# print evotime, fid, fid_u
					if do_plot:
						plt.errorbar(np.sort(evotime),fid[np.argsort(evotime)],fid_u[np.argsort(evotime)],fmt=datapoint,label=legendlabel,color=color_list[i],capsize = GlobalErrorCapsize,zorder=zorder)

				else:
					if do_plot:
						plt.errorbar(np.sort(evotime),fid[np.argsort(evotime)],fid_u[np.argsort(evotime)],marker=datapoint,color=color_list[i],label = legendlabel,capsize = GlobalErrorCapsize,zorder=zorder)
				

			#increment for colors etc.

			if not plot_states:
				i +=1			
				

		if add_horizontal_line != None:
			plt.plot(np.linspace(-1,evotime[-1],101),[add_horizontal_line]*101,'--',lw=fit_lw,color = 'gray')

		if add_vertical_line != None:
			pass

		if ylim !=None:
			plt.ylim(ylim)

		if xlim != None:
			plt.xlim(xlim)

		### plot
		plt.xlabel(r'$\tau$ (ms)')
		
		if ylabel == None:
			if self.process_fidelities:
				ylabel = 'Process fidelity'
				save_name +='_proc'

			if not self.contrast and not self.process_fidelities or plot_states:
				ylabel = 'State fidelity'
				save_name +='_states'


			if self.contrast:
				ylabel = 'Expectation value'
				save_name += '_contrast'

		plt.ylabel(ylabel)
		
		if legend:
			# get handles
			handles, labels = ax.get_legend_handles_labels()
			# remove the errorbars
			handles = [h[0] for h in handles]
			
			
			#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
			ncol = kw.pop('legend_ncol', 2)
			# handles them in the legend
			ax.legend(handles, labels, loc='upper right',numpoints=1,frameon = False,ncol = ncol)


		if xticks != None:
			plt.xticks(xticks)

		if yticks != None:
			plt.yticks(yticks)


		plt.title(self.plot_title)

		if save_plot:
			self.save_plot(save_name)

		plt.show()
		plt.close('all')



	def plot_msmt_slices(self,x_values,save_plot=False,state = 'X',**kw):
		""" 
		takes a list of times 
		plots this list and compares theory with experiment.
		BEWARE!!!
		takes data points closest to the specified x_values!!
		the timing of that data is then given in the legend.
		"""

		yticks = kw.pop('yticks',None)
		ylim = kw.pop('ylim',None)
		legend = kw.pop('legend',False)



		fig = plt.figure()
		ax = plt.subplot()
		x_axis_ticks = []

		for key in sorted(self.pickle_dict.keys()):#self.pickle_dict.iteritems():

			if self.use_averaged:
				data = self.pickle_dict[key]['state_dict_averaged'][state]
			else:
				data = self.pickle_dict[key]
			evotime = data['evotime']

			final_x_vals = []
			final_y_vals = []
			final_u_vals = []
			for x in x_values:
				ind = find_nearest_idx(evotime,x)
				final_x_vals.append(evotime[ind])
				final_y_vals.append(data['fid'][ind])
				final_u_vals.append(data['fid_u'][ind])

			i = self.find_index(key)
			# plt.bar(i+0.*bar_width,self.acquire_fit_val(key,final_x_vals),bar_width,color = '#CD4A4A',linewidth = 0.5, edgecolor = 'black',hatch='//',alpha=0.5)
			plt.bar(i-0.5*bar_width,(np.array(final_y_vals)-0.5),bar_width,yerr=final_u_vals,bottom = 0.5,color='gray',linewidth = 1, 
																								edgecolor = 'black',
																								label = str(final_x_vals[0]), 
																								error_kw=dict(lw = 1,capthick=1,
																									ecolor = 'black'))
			
			x_axis_ticks.append(int(key))

		# plt.legend()
		if not self.contrast:
			plt.ylim(0.45)
			plt.plot(np.linspace(-1,len(self.pickle_dict.keys()),101),[0.5]*101,'-',color='black',lw=1)
		plt.xlim(-bar_width)
		if legend:
			plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

		plt.xlabel('N')
		if not self.contrast and not self.process_fidelities:
			ylabel = 'State fidelity'
		if self.process_fidelities:
			ylabel = 'Process fidelity'

		if self.contrast:
			ylabel = 'Contrast'

		if yticks != None:
			plt.yticks(yticks)
		if ylim != None:
			plt.ylim(ylim)

		plt.tick_params(axis='x',which='both',top='off')
		plt.ylabel(ylabel)
		plt.xticks(range(len(self.pickle_dict.keys())),x_axis_ticks)
		if save_plot:
			self.save_plot('msmt_sclices')


		plt.show()
		plt.close('all')

	def save_plot(self,plot_type):
		"""
		plot_type: strings which is used for the file name.
		"""

		if self.fancy:
			plt.savefig(os.path.join(self.save_folder,self.name+' '+plot_type+'.pdf'),format='pdf',bbox_inches = 'tight',pad_inches=0.2)
		else:
			plt.savefig(os.path.join(self.save_folder,self.name+' '+plot_type+'.png'),format='png',bbox_inches = 'tight',pad_inches=0.3)


	def find_index(self,key):

		### find the position in the data dictionary
		i = 0
		for keys in sorted(self.pickle_dict.keys()):
			if key == keys:
				break
			i += 1
		return i
	

	def set_do_fit(self,do_fit):

		self.do_fit = do_fit


	def average_state_dicts(self):
		"""
		averages over orthogonal states in all state dicts and creates the dictionary state_dict_averaged.
		which includes the keys 'X' 'Y' 'Z'. Each of these keys can then be plotted with the timetrace routine

		Obviously only works for data involving process fidelities.
		"""
		# timetrace_keylist
		for key, data in self.pickle_dict.iteritems():
			if key in ['0','1','2','3','4','5','6','00','02','04','08','10','12','16']:
				state_dict_averaged = {}
				state_dict = data['state_dict']
				
				i = 0
				j = 0
				k = 0

				for key2 in sorted(state_dict.keys()):
					
					if i == 0 and 'X' in key2:
						state_dict_averaged['X'] = {}
						state_dict_averaged[key2[-1]] = cp.deepcopy(state_dict[key2])
						i = 1

					elif j == 0 and 'Y' in key2:
						state_dict_averaged['Y'] = {}
						state_dict_averaged[key2[-1]] = cp.deepcopy(state_dict[key2])
						j = 1

					elif k == 0 and 'Z' in key2:
						state_dict_averaged['Z'] = {}
						state_dict_averaged[key2[-1]] = cp.deepcopy(state_dict[key2])
						k = 1

					else:
						### average the two states.
						state_dict_averaged[key2[-1]][timetrace_keylist[1]] = [(state_dict[key2][timetrace_keylist[1]][kk]+fid)/2 for kk,fid in enumerate(state_dict_averaged[key2[-1]][timetrace_keylist[1]])]
						state_dict_averaged[key2[-1]][timetrace_keylist[2]] = [(np.sqrt(state_dict[key2][timetrace_keylist[2]][kk]**2+fid_u**2)/2) for kk,fid_u in enumerate(state_dict_averaged[key2[-1]][timetrace_keylist[2]])]

					# i+=1

			self.pickle_dict[key]['state_dict_averaged'] = {}
			self.pickle_dict[key]['state_dict_averaged'] = state_dict_averaged




######################
######################
##					##
##		Fig 1A		##
##					##
######################
######################

def X_preservation():
	
	filename = 'X_preservation.p'
	pickle_dict = open_data(filename)


	#### do some after processing
	data_dict = {}

	for key in pickle_dict.keys():
		data_dict[key] = {}
		data_dict[key]['state_dict'] = {}
		data_dict[key]['state_dict']['X'] = pickle_dict[key]


	### uncomment for single x 
	THE_DATA = data_object("Fig1A",pickle_dict)
	THE_DATA.set_do_fit(True)
	THE_DATA.physical_model_fit = False
	THE_DATA.pheno_physical_mode_fit = True

	# THE_DATA.plot_timetrace(save_plot=True)
	

	# # THE_DATA.plot_msmt_slices([8],save_plot=True)
	# THE_DATA.plot_msmt_slices([45],save_plot=True)

	#######uncomment to average states.
	# filename = 'mX_preservation.p'
	# pickle_dict = open_data(filename)
	# # print pickle_dict

	# for key in pickle_dict.keys():
	# 	data_dict[key]['state_dict']['mX'] = pickle_dict[key]

	# THE_DATA = data_object("Fig1A",data_dict)
	# THE_DATA.average_state_dicts()
	# THE_DATA.use_averaged = True
	###############

	

	# print data_dict['10']['state_dict_averaged']
	THE_DATA.set_do_fit(True)
	THE_DATA.physical_model_fit = False
	THE_DATA.pheno_physical_mode_fit = True
	THE_DATA.do_not_plot_list = ['10','12']

	THE_DATA.plot_timetrace(save_plot=True,plot_states = False,xlim=[0,100],xticks=[0,40,80],yticks=[0.4,0.6,0.8,1])

	THE_DATA.plot_msmt_slices([40],legend = False,save_plot = True,yticks=[0.4,0.6,0.8,1],ylim=[0.4,1])#0.45 and 0.8

	# THE_DATA.plot_msmt_slices([8],save_plot=True)
	# THE_DATA.plot_msmt_slices([45],save_plot=True)


	# THE_DATA.physical_model_fit = False
	# THE_DATA.plot_state_timetrace(save_plot=False)
	

	### extract coherence times.
	x_axis = []
	results = []
	results_u = []


	### uncomment if you want to use average over X and -X	
	# for k in sorted(pickle_dict.keys()):

	# 	print k

	# 	try:
	# 		res = abs(THE_DATA.pickle_dict[k]['state_dict_averaged']['X']['fit']['params_dict']['t'])
	# 		res_u = THE_DATA.pickle_dict[k]['state_dict_averaged']['X']['fit']['error_dict']['t']

	# 	except:
	# 		res = abs(THE_DATA.pickle_dict[k]['state_dict_averaged']['X']['fit']['params_dict']['sigma'])
	# 		res_u = THE_DATA.pickle_dict[k]['state_dict_averaged']['X']['fit']['error_dict']['sigma']
	# 	x_axis.append(int(k))
	# 	results.append(res)
	# 	results_u.append(res_u)



	### uncomment if you want to use only one state (X or -X)		
	for k in sorted(pickle_dict.keys()):

		print k

		try:
			res = abs(THE_DATA.pickle_dict[k]['fit']['params_dict']['t'])
			res_u = THE_DATA.pickle_dict[k]['fit']['error_dict']['t']

		except:
			res = abs(THE_DATA.pickle_dict[k]['fit']['params_dict']['sigma'])
			res_u = THE_DATA.pickle_dict[k]['fit']['error_dict']['sigma']
		x_axis.append(int(k))
		results.append(res)
		results_u.append(res_u)

	results_norm, results_u_norm = ZAnalysis.normalize_and_plot_coherence(x_axis,results,results_u)

	coherence_dict = {	'msmts' 	: x_axis,
						'results'	: results_norm,
						'results_u'	: results_u_norm}

	ZAnalysis.save_data(coherence_dict,"Zeno_1Q_scaling_"+"class_bit"+".p")




######################
######################
##					##
##		Fig 1b		##
##					##
######################
######################

def Q1_proc_fid():


	filename = 'FIG1B_1Q_proc_fid.p'
	pickle_dict = open_data(filename)

	print pickle_dict['1'].keys()

	THE_DATA = data_object("Fig1B",pickle_dict,proc_fid=True)



	THE_DATA.plot_timetrace(save_plot=True,add_horizontal_line=0.5,xlim=[0,50],ylim=[0.4,0.55],legend=False,xticks=[0,25,50],yticks=[0.4,0.5])
	

	THE_DATA.plot_timetrace(save_plot=True,plot_states = True,xticks=[0,25,50],yticks=[0.4,0.6,0.8,1])

	THE_DATA.average_state_dicts()


######################
######################
##					##
##		Fig 2A		##
##					##
######################
######################

def Q2_vs_Q1():

	### even msmts


	filename = 'FIG2A_evenmsmts.p'
	pickle_dict = open_data(filename)

	plot_avg_fid = True
	if plot_avg_fid:
		yticks = [0.6,0.8,1]
		ylim = [0.55,1]
		ylabel = 'Average state fidelity'
	else:
		yticks = [0.4,0.7,1]
		ylim = None
		ylabel = 'Process fidelity'

	data_even = data_object("Fig2A_even",pickle_dict, proc_fid = True, avg_fid = plot_avg_fid,legend_ncol=1)

	### prepare fitting

	# data_even.use_fixed = False ### use if yopu want to fit gaussian functions
	data_even.set_do_fit(True)
	data_even.physical_model_fit = False
	data_even.pheno_physical_mode_fit = True

	### plot & fit

	data_even.plot_timetrace(save_plot  = True, 
								legend = True,
								xlim=[0,60],
								ylabel=ylabel,
								ylim=ylim,
								xticks=[0,30,60],
								yticks=yticks,
								legend_ncol=1,
								add_horizontal_line = 2./3.)

	data_even.average_state_dicts()

	#####################
	#####################
	#####################

	### odd msmts.

	filename = 'FIG2A_oddmsmts.p'
	pickle_dict = open_data(filename)


	data_odd = data_object("Fig2A_odd",pickle_dict, proc_fid = True, avg_fid = plot_avg_fid)

	# data_even.use_fixed = False ### use if yopu want to fit to a gaussian
	data_odd.set_do_fit(True)
	data_odd.physical_model_fit = False
	data_odd.pheno_physical_mode_fit = True

	data_odd.plot_timetrace(save_plot  = True, legend = True,xlim=[0,60],ylabel=ylabel,ylim=ylim,xticks=[0,30,60],yticks=yticks,legend_ncol=1)


######################
######################
##					##
##		Fig 2B		##
##					##
######################
######################




def TwoQFidelity():

	### this data is pre-evaluated via TwoQ_Zeno_analysis.py (see the respective ipython notebook for 2 qubits)
	
	### choose FIG2B_2Qfid_preservation.p for all 6 states (including X,X and -X,-X)
	filename = 'FIG2B_2Qfid_preservation_entangled.p'

	pickle_dict = open_data(filename)

	THE_DATA = data_object("Fig2B",pickle_dict, proc_fid = True) ### proc_fid only determines the analysis type in this case.

	THE_DATA.physical_model_fit = False
	THE_DATA.pheno_physical_mode_fit = True
	THE_DATA.set_do_fit(True)


	THE_DATA.plot_title = r''
	THE_DATA.plot_timetrace(save_plot 	= True, 
							legend 		= True,
							xlim		= [0,60],
							add_horizontal_line = 0.5,
							ylim		= [0.3,0.9],
							xticks 		= [0,30,60],
							yticks 		= [0.3,0.5,0.7,0.9],
							ylabel 		= 'Entangled state fidelity',
							legend_ncol = 1)



######################
######################
##					##
##		Fig 3		##
##					##
######################
######################


# color_list = ['b','r','0.5','m','y','black','brown','c','orange','0.75','0.0','0.5']
# color_list = ['#e41a1c','#377eb8','#4D4D4D','#FAA43A','#4daf4a','','#4daf4a']


def Protecting_XXX(states = True):
	"""
	makes four different plots: three for logical state fidelities, one for the decay of XXX averaged over the three states. 
	"""

	filename = "FIG3AB_3QZeno.p"
	pickle_dict = open_data(filename)

	THE_DATA = data_object("Fig3AB",pickle_dict)
	THE_DATA.state_key_endings = ['00p11']
	THE_DATA.state_key_start = 0
	# print pickle_dict


	#### Get fitting parameters ready
	THE_DATA.set_do_fit(True)
	THE_DATA.physical_model_fit = False
	THE_DATA.pheno_physical_mode_fit = True
	THE_DATA.fit_offset = 0.25
	THE_DATA.use_fixed = False
	THE_DATA.process_fidelities = True


	if states:
		#old config for horizontal alignment:
		### Do the fits. xlim=[0,60],ylim=[0.2,0.9],xticks=[0,25,50],yticks=[0.2,0.5,0.8]

		THE_DATA.plot_title = r'$|0,0\rangle_L$'
		THE_DATA.state_key_endings = ['00']
		THE_DATA.name = 'Fig3AB_00'
		THE_DATA.plot_timetrace(save_plot = True, plot_states = True,xlim=[0,60],ylim=[0.15,0.9],xticks=[0,25,50],yticks=[0.2,0.5,0.8],legend_ncol = 1)

		THE_DATA.name = 'Fig3AB_X0'
		THE_DATA.plot_title = r'$|X,0\rangle_L$'
		THE_DATA.state_key_endings = ['00p10']
		THE_DATA.plot_timetrace(save_plot = True, plot_states = True,xlim=[0,60],ylim=[0.3,0.9],xticks=[0,25,50],yticks=[0.4,0.6,0.8],legend_ncol = 1)

		THE_DATA.name = 'Fig3AB_00p11'
		THE_DATA.plot_title = r'$|\Phi^+ \rangle_L$'
		THE_DATA.state_key_endings = ['00p11']
		THE_DATA.plot_timetrace(save_plot = True, plot_states = True,xlim=[0,60],ylim=[0.3,0.9],xticks=[0,25,50],yticks=[0.4,0.6,0.8],legend_ncol = 1)


	else: 
		THE_DATA.name = 'Fig3AB_xxx'
		THE_DATA.process_fidelities = False
		THE_DATA.contrast = True ### sets a different label for plotting --> 'Expectation value'.
		THE_DATA.physical_model_fit = False
		THE_DATA.pheno_physical_mode_fit = True
		THE_DATA.plot_title = ''# r'$\langle XXX \rangle $ averaged over three states'
		THE_DATA.state_key_endings = ['xxx']
		THE_DATA.plot_timetrace(save_plot = True, plot_states = True,xlim=[0,60],xticks=[0,25,50],yticks=[0.0,0.4,0.8],legend_ncol = 1)




######################
######################
##					##
##		Fig 4		##
##					##
######################
######################

def coherence_scalingv2(theory = False,group_data = False,weight_XX_and_YY = False):
	### only used if theory == True
	scaling_parameters = [1,0,2.142,0,2.8374,0,3.3893,0,3.8631,0,4.28478,0,4.66859,0,5.02315,0,5.3543]

	scale_figure(0.5577)

	data2Q = []
	data2Q.append(open_data("Zeno_2Q_scaling_Z.p"))
	data2Q.append(pickle.load( open( "ZenData\Zeno_2Q_scaling_mZ.p", "rb" ) ))
	data2Q.append(pickle.load( open( "ZenData\Zeno_2Q_scaling_Y.p", "rb" ) ))
	data2Q.append(pickle.load( open( "ZenData\Zeno_2Q_scaling_mY.p", "rb" ) ))
	data2Q.append(pickle.load( open( "ZenData\Zeno_2Q_scaling_XX.p", "rb" ) ))
	data2Q.append(pickle.load( open( "ZenData\Zeno_2Q_scaling_YY.p", "rb" ) ))

	# data1Q = pickle.load( open( "ZenData\Zeno_1Q_scaling_mX.p", "rb" ) )
	data1Q = pickle.load( open( "ZenData\Zeno_1Q_scaling_class_bit.p", "rb" ) )

	data3Q = pickle.load( open( "ZenData\Zeno_3Q_scaling.p", "rb" ) )

	dict_keys = ['msmts','results','results_u']

	# print 'this is the 2Q data'
	# print 'YY:', data2Q[-1] 
	# print 'XX:', data2Q[-2] 


	if weight_XX_and_YY:
		#####################
		### we calculate a weighted average for the 2Q data of XX and YY! (the error bar of XX is much lower)
		### See also RJ Barlow Statistics p. 83 for a choice of the weights according to Maximum likelihood estimation.
		#########################

		YY = data2Q[-1] 
		XX = data2Q[-2]

		YY['weights'] = 1./np.array(YY['results_u'])**2
		XX['weights'] = 1./np.array(XX['results_u'])**2
		weight_normalization = YY['weights'] + XX['weights']
		weighted_averages = (YY['weights']*np.array(YY['results'])+XX['weights']*np.array(XX['results']))/weight_normalization
		weighted_averages_u = np.sqrt((XX['weights']*XX['results_u'])**2+(YY['weights']*YY['results_u'])**2)/weight_normalization

		print 'Results for weighing XX and YY of the 2Q data'
		print '------------------------------------------------'
		print 'XX     ', np.round(XX['results'],3)
		print 'error  ', np.round(XX['results_u'],3)
		print 'weights',np.round(XX['weights'],1)
		print 'YY     ', np.round(YY['results'],3)
		print 'error  ', np.round(YY['results_u'],3)
		print 'weights',np.round(YY['weights'],1)
		print '-------------------------------------------------'
		print 'weighted average', np.round(weighted_averages,3)
		print 'uncertainty     ', np.round(weighted_averages_u,3)

		#### overwrite XX and YY with the weighted average in data2Q for further analysis
		data2Q = data2Q[:-1]
		data2Q[-1] = {'msmts': [0,2,4,6],'results':weighted_averages,'results_u':weighted_averages_u}


	### combine all single qubit expectation values of the two qubit measurements.
	### i.e. <XI>, <YZ> etc.

	Q2_data = {'msmts': [], 'results': [], 'results_u':[]}

	Q2_data['results'] = np.array(data2Q[0]['results'])/4.
	Q2_data['results_u'] = np.array(data2Q[0]['results_u'])**2/16

	for i in range(3):
		Q2_data['results'] = Q2_data['results'] + np.array(data2Q[i+1]['results'])/4.

		Q2_data['results_u'] = np.array(data2Q[i+1]['results_u'])**2/16

	Q2_data['results_u'] = np.sqrt(Q2_data['results_u'])
	Q2_data['msmts'] = data2Q[0]['msmts']


	# print Q2_data['results']


	### combine the decay of all logical two qubit states for the 3 qubit experiments.

	Q3_data = {'msmts': [], 'results': [], 'results_u':[]}

	Q3_data['results'] = np.array(data3Q['results'][0])/3.
	Q3_data['results_u'] = np.array(data3Q['results_u'][0])**2/9

	for i in range(2):
		Q3_data['results'] = Q3_data['results'] + np.array(data3Q['results'][i+1])/3.

		Q3_data['results_u'] = np.array(data3Q['results_u'][i+1])**2/9

	Q3_data['results_u'] = np.sqrt(Q3_data['results_u'])
	Q3_data['msmts'] = data3Q['msmts'][0]

	
	# print data1Q['results']


	if theory:
		for ii,msmts in enumerate(data1Q['msmts']):
			data1Q['results'][ii] =  scaling_parameters[msmts]*data1Q['results'][ii]
			data1Q['results_u'][ii] =  scaling_parameters[msmts]*data1Q['results_u'][ii]

		for ii,msmts in enumerate(Q2_data['msmts']):
			Q2_data['results'][ii] =  scaling_parameters[msmts]*Q2_data['results'][ii]
			Q2_data['results_u'][ii] =  scaling_parameters[msmts]*Q2_data['results_u'][ii]
			
			data2Q[-1]['results'][ii] = scaling_parameters[msmts]*data2Q[-1]['results'][ii]
			data2Q[-1]['results_u'][ii] = scaling_parameters[msmts]*data2Q[-1]['results_u'][ii]

			if not weight_XX_and_YY:
				data2Q[-2]['results'][ii] = scaling_parameters[msmts]*data2Q[-2]['results'][ii]
				data2Q[-2]['results_u'][ii] = scaling_parameters[msmts]*data2Q[-2]['results_u'][ii]

		for ii,msmts in enumerate(Q3_data['msmts']):
			Q3_data['results'][ii] =  scaling_parameters[msmts]*Q3_data['results'][ii]
			Q3_data['results_u'][ii] =  scaling_parameters[msmts]*Q3_data['results_u'][ii]
			data3Q['results'][3][ii] = scaling_parameters[msmts]*data3Q['results'][3][ii]
			data3Q['results_u'][3][ii] = scaling_parameters[msmts]*data3Q['results_u'][3][ii]

	# print Q2_data
	# print Q3_data
	# print data2Q
	combined_data = [data1Q]
	combined_data.extend([Q2_data])
	
	combined_data.append(data2Q[-1]) ### YY scaling (or weighted average of XX and YY) for 2 qubit experiments

	if not weight_XX_and_YY:
		combined_data.append(data2Q[-2]) ### XX scaling for 2 qubit experiments

	combined_data.extend([Q3_data])
	# print Q3_data
	### add the XXX scaling to the data.
	combined_data.append({'msmts':data3Q['msmts'][3], 'results':data3Q['results'][3],'results_u':data3Q['results_u'][3]})

	# print combined_data
	# TwoMsmts = []

	labels = [r'$\langle X \rangle$',r'$\langle IX \rangle$ etc.',r'$\langle XX \rangle$',r'$\langle YY \rangle$',r'$\langle XIX \rangle$ etc.',r'$\langle XXX \rangle$']
	
	if group_data:
		labels = [r'1 $\mathregular{^{13}C}$',r'2 $\mathregular{^{13}C}$',r'3 $\mathregular{^{13}C}$']


		### group all data where two qubit was subject to dephasing
		ds = [data2Q[-2],data2Q[-1],Q3_data]

		if weight_XX_and_YY:
			ds = [data2Q[-1],Q3_data]


		Q2_all_data = {}
		for k in data2Q[-2].iterkeys():
			# create list of lists for a certain key
		    Q2_all_data[k] = tuple(d[k] for d in ds)
		    #flatten that list to be 1D (this is needed for plotting.)
		    Q2_all_data[k] = [item for sublist in Q2_all_data[k] for item in sublist]

		### do the same for 1 qubit/carbon
		ds = [data1Q,Q2_data]
		Q1_all_data = {}
		for k in data2Q[-2].iterkeys():
			# create list of lists for a certain key
		    Q1_all_data[k] = tuple(d[k] for d in ds)
		    #flatten that list to be 1D (this is needed for plotting.)
		    Q1_all_data[k] = [item for sublist in Q1_all_data[k] for item in sublist]

		###redefine combined_data
		combined_data = [Q1_all_data,Q2_all_data]
		### add the XXX scaling to the data.
		combined_data.append({'msmts':data3Q['msmts'][3], 'results':data3Q['results'][3],'results_u':data3Q['results_u'][3]})


	for ii,data in enumerate(combined_data):

		### extract the normalized values for two measurements.
		# if 2 in data[dict_keys[0]] :
		# 	index = data[dict_keys[0]].index(2)
		# 	TwoMsmts.append(data[dict_keys[1]][index])


		plt.errorbar(data[dict_keys[0]],data[dict_keys[1]],data[dict_keys[2]],fmt=datapoint,label =labels[ii],color = color_list[ii],capsize = 2)

	### let the theory go through the average of two measurements
	### theory predicts an enhancement of 2.14 for 2 meaurements	
	### we measure more or less 3.2
	### TODO: implement average. comment: this does not give any physical intuition.
	# print TwoMsmts
	# offset = np.mean(np.array(TwoMsmts))/2.14
	# print 'OFFSET', offset

	def f(N): return 1+0.77186*N**0.626286#np.exp(0.403+0.4628*np.log(N))
	
	f = np.vectorize(f)

	if theory:
		plt.plot(np.linspace(0,19,1001),f(np.linspace(0,19,1001)),color = 'black')

	plt.xlabel(r'$N$')
	
	# ax.set_yscale('log')
	# ax.set_xscale('log')
	if theory:
		plt.ylabel('Normalized decay time')
		plt.xlim(-0.2,18)
		plt.ylim(0.9,6.3)
	else:
		plt.ylabel(r'$T/T(N=0)$')
		plt.xlim(-0.2,18)
		plt.ylim(0.8,1.4)

	plt.xticks([0,4,8,12,16])
	# plt.legend(loc=4,frameon=False,numpoints = 1) #,



	save_folder = r'D:/measuring/data/Zeno_results'
	plt.savefig(os.path.join(save_folder,'Fig4_coherence_scaling'+'.pdf'),format='pdf',bbox_inches = 'tight',pad_inches=0.2)
	plt.savefig(os.path.join(save_folder,'Fig4_coherence_scaling'+'.png'),format='png',bbox_inches = 'tight',pad_inches=0.3)


	plt.show()
	plt.close('all')

######################
######################
##					##
##	  Supplement	##
##					##
######################
######################

### DESR after MBI: Nitrogen population
### See the executed file for details. This script sets the plot styling only.
def getNpop():
	execfile(r'D:/measuring/analysis/scripts/Zeno/MBI_population.py')

def oddmsmts():
	filename = 'FIG2A_oddmsmts.p'
	pickle_dict = open_data(filename)

	plot_avg_fid = True
	if plot_avg_fid:
		ylabel = 'Average state fidelity'
		ylim = [0.55,0.95]
		yticks = [0.6,0.7,0.8,0.9]
	else:
		ylabel = 'Logical process fidelity'
		ylim = [0.3,0.9]
		yticks = [0.3,0.4,0.5,0.6,0.7,0.8,0.9]

	data_odd = data_object("SuppM_Fig2A_odd",pickle_dict, proc_fid = True, avg_fid = plot_avg_fid)

	# data_even.use_fixed = False ### use if yopu want to fit to a gaussian
	data_odd.set_do_fit(True)
	data_odd.physical_model_fit = False
	data_odd.pheno_physical_mode_fit = True

	data_odd.plot_timetrace(save_plot  = True, legend = True,xlim=[0,90],ylabel=ylabel,ylim=ylim,xticks=[0,30,60,90],yticks=yticks)

def Protecting_XX():
	filename = 'FIG2B_XX_preservation.p'

	pickle_dict = open_data(filename)

	THE_DATA = data_object("SuppM_XX_decay",pickle_dict)

	THE_DATA.contrast = True
	THE_DATA.use_fixed = False

	THE_DATA.fit_offset = 0.0 ### force the gaussian to decay to 0.

	THE_DATA.physical_model_fit = False
	THE_DATA.pheno_physical_mode_fit = True
	THE_DATA.set_do_fit(True)


	THE_DATA.plot_title = r'$\langle XX \rangle $ averaged over six states'
	THE_DATA.plot_timetrace(save_plot = True, legend = True,xlim=[0,60],xticks=[0,20,40,60],yticks=[0.0,0.4,0.8])