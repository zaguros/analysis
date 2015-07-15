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
import matplotlib
import copy as cp
import analysis.scripts.Zeno.Zeno_fitting_tools as Zfits
reload(Zfits)
reload(common)
reload(toolbox)

"""
Analyze the zeno data
NK 2014
"""

#############################################

"Global variables"

#############################################


color_list = ['b','g','r','m','y','black','brown','c','orange','0.75']

font = {'family' : 'sans-serif',
    'weight' : 'normal',
    'size'   : 12}

matplotlib.rc('font', **font)
fit_lw = 1.5



results_folder = r'D:\measuring\data\Zeno_results' ### for final plots

timetrace_keylist = ['evotime','fid','fid_u'] ### XXXXXXX don't change!!!
### needed for saving the data and evaluating it with Paper_data_evalution.py


CarbonRO_correction = True ### or False. Your choice.

############################################


########################

""" Helper functions """

########################
 
def save_data(pickle_dict,filename,folder="ZenData\ "):
	"""
	takes a dictionary and a filename as input.
	Filename can also be a path.
	"""
	### dump to file.
	if folder != None:		
		fileOut = open(folder[:-1]+filename,"wb")
		print folder[:-1]+filename
	else:
		fileOut = open(filename,"wb")
	pickle.dump(pickle_dict,fileOut)
	fileOut.close

def calculate_max_fidelity(fid,fid_u,fid_ort,fid_u_ort):
	"""
	calculate the maximum achievable fidelity by quadratically adding two orthogonal vectors.
	input: 	2 lists of fidelities (fid and fid_ort)
			2 lists of uncertainties for these fidelity values (fid_u and fid_u_ort)

	output:
		fid_bloch maximum obtainable fidelity (list)
		fid_u_bloch corresponding uncertainty (list)
	"""

	### add both data sets together such that the Bloch-vector length in the XY plane is obtained. (geometrically)
	### it is assumed that both measurements form a consistent data set and have the same evolution times.
	fid_bloch=[]
	fid_u_bloch=[]
	for ii in range(len(fid)):
		fid_bloch.append(np.sqrt((fid[ii]-0.5)**2+(fid_ort[ii]-0.5)**2)+0.5)

		### for a derivation of the error propagation see the Zeno onenote file. 2015-04-12 NK

		err1 = (0.25+(fid[ii]-1)*fid[ii])*fid_u[ii]**2
		err2 = (0.25+(fid_ort[ii]-1)*fid_ort[ii])*fid_u_ort[ii]**2
		denominator = 0.5 + fid[ii]*(fid[ii]-1) + fid_ort[ii]*(fid_ort[ii]-1)
		err = np.sqrt((err1+err2)/denominator)
		fid_u_bloch.append(err)

	return fid_bloch, fid_u_bloch

def calc_error_propagation_fraction(N1,normalisation,u1,uNormalisation):
	"""
	calculates the error propagation onto a normalized value:
		The function under consideration is f = N1/normalisation.
		Where both values are subject to a statistical uncertainty.
	"""

	return np.sqrt((u1*normalisation)**2+(uNormalisation*N1)**2)/normalisation**2

def normalize_and_plot_coherence(x_axis,results,results_u):
	"""
	input: 3 lists --> x_axis, results and the uncertainty of the values in results
	x_axis, gives the x_axis of the plot to follow

	the function takes the minimum value of results and normalizes 
	"""
	print results
	print results_u

	normalisation = results[0]
	norm_index = results.index(results[0])
	for ii in range(len(results_u)):
				results_u[ii] = calc_error_propagation_fraction(results[ii],normalisation,results_u[ii],results_u[norm_index])
	normalised_results = [x/normalisation for x in results]

	print x_axis
	print results_u
	print normalised_results

	fig=plt.figure()
	ax=plt.subplot()

	plt.errorbar(x_axis,normalised_results,results_u,fmt='o')
	plt.xlim(x_axis[0]-0.5,x_axis[-1]+0.5)
	plt.xlabel('Number of Measurements')
	plt.ylabel('Normalized coherence time')
	plt.show()
	plt.close('all')

	return normalised_results, results_u


def calc_average(list1,list2,list_u1,list_u2):
	"""
	takes two np.arrays of measurement results and sums them elementwise.
	quadratically adds the uncertainties for these results.
	Function is predominantly used to average states like Z and -Z etc.

	returns the sum and the uncertainties as numpy array objects.
	"""

	result_list = []
	result_u_list = []

	for ii in range(len(list1)):
		result_list.append((list1[ii]+list2[ii])/2.)
		result_u_list.append((list_u1[ii]**2+list_u2[ii]**2)/4)
	result_u_list = np.sqrt(np.array(result_u_list))


	return np.array(result_list),result_u_list

def calc_state_average(fid_arr,fid_u_arr):
	"""
	takes 2D arrays of shape (6,X) where X is unspecified.
	averages the list entries then pairwise and returns 2D arrays of shape (3,X)
	"""

	result_list = []
	result_u_list= []
	for i in range(3):
		int_list, int_u_list = calc_average(fid_arr[2*i],fid_arr[2*i+1],fid_u_arr[2*i],fid_u_arr[2*i+1])
		result_list.append(int_list)
		result_u_list.append(int_u_list)

	return result_list, result_u_list

def make_pickle_dict(x_axis,y_axis,y_u,key_list):
	return {key_list[0]: x_axis, key_list[1]:y_axis, key_list[2]: y_u}

def fit_State_decay(msmts,ax,A0,evotime,fid, t = 8.25):
	"""
	takes a zeno data set for a specific state and number of measurements and returns the fitted result.
	Inputs:

	msmts 			string which signifies the number of measurements
	A0 				the amplitude for the 0 measurement case.
	evotime 		np.array with the evolution times of the data points
	fid				np.array with the measured state fidelities.

	Output:

	result_string 	a string which is used for labelling the fits
	fit_result		the fitted function for plotting.
	"""

	t1 = 9./np.sqrt(2)
	# t2 = 8.25
	t2 = t
	p = 0.25
	repump = 0.95

	# print 'A0 ', A0
	# print 't2 ', t2


	if msmts == '0':
		p0, fitfunc, fitfunc_str = common.fit_gauss(0.5, 0.43, 0., t)

	elif msmts == '1':
		### we usually do not take repumping errors into account for this fit.
		### therefore t1 is disregarded.
		p0, fitfunc,fitfunc_str = Zfits.fit_1msmt_state_fid(A0,t1,t2, p,repump,False)

		### manual option to show the intial guess in the plot.
		if False:
			ax.plot(np.linspace(0.,120.0,201), fitfunc(np.linspace(0.,120.0,201)), ':', lw=fit_lw)

	elif msmts == '2':
		p0, fitfunc,fitfunc_str = Zfits.fit_2msmt_state_fid(A0,t, p)

	elif msmts == '3':
		p0, fitfunc,fitfunc_str = Zfits.fit_3msmt_state_fid(A0,t, p)

	elif msmts == '4':
		p0, fitfunc,fitfunc_str = Zfits.fit_4msmt_state_fid(A0,t, p)

	elif msmts == '5':
		p0, fitfunc,fitfunc_str = Zfits.fit_5msmt_state_fid(A0,t, p)

	elif msmts == '6':
		p0, fitfunc,fitfunc_str = Zfits.fit_6msmt_state_fid(A0,t, p)

	elif msmts == '8':
		p0, fitfunc,fitfunc_str = Zfits.fit_8msmt_state_fid(A0,t, p)

	elif msmts == '10':
		p0, fitfunc,fitfunc_str = Zfits.fit_10msmt_state_fid(A0,t, p)

	elif msmts == '12':
		p0, fitfunc,fitfunc_str = Zfits.fit_12msmt_state_fid(A0,t, p)

	elif msmts == '16':
		p0, fitfunc,fitfunc_str = Zfits.fit_16msmt_state_fid(A0,t, p)
		if False:
			ax.plot(np.linspace(0.,120.0,201), fitfunc(np.linspace(0.,120.0,201)), ':', lw=fit_lw)

	### msmts = 0 is an exception
	fixed = [0,1]
	if msmts =='0':
		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
	else:
		fixed = [0,1] ### fixed parameter: [0,1] --> fix time and amplitude, p is the only free parameter. 
										###[0] --> fix the amplitude for 0 measurements only.

		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)

	# print fit_result['params']

	p1 = str(round(fit_result['params'][0]*100,3))
	p1_u = str(round(fit_result['error'][0]*100,3))
	# p1 = 'test'
	# p1_u = 'ing'
	result_string = p1 + ' +- ' + p1_u 

	return fit_result, result_string


def fit_process_decay(msmts,ax,A0,evotime,fid,decoded_bit,offset0 = 0.415, t = 8.2):
	"""
	takes a zeno data set for a specific number of measurements and returns the a fit to the data.
	Inputs:

	msmts 			string which signifies the number of measurements
	A0 				the amplitude for the 0 measurement case.
	evotime 		np.array with the evolution times of the data points
	fid				np.array with the measured state fidelities.

	Output:

	result_string 	a string which is used for labelling the fits
	fit_result		the fitted function for plotting.
	"""
	### decoding to carbon 2
	if decoded_bit == '2':
		# t = 8.34
		p = 0.08
		# offset0 = 0.40

	### decoding to carbon 1
	elif decoded_bit == '1':
		# t = 5.81
		p = 0.08
		# offset0 = 0.430

	if msmts == '0':
		p0, fitfunc, fitfunc_str = common.fit_gauss(0.4, 0.43, 0., t)


	elif msmts == '1':
		p0, fitfunc,fitfunc_str = Zfits.fit_1msmt_proc_fid(A0,offset0,t, p)

		### manual option to show the intial guess in the plot.
		if False:
			ax.plot(np.linspace(0.,120.0,201), fitfunc(np.linspace(0.,120.0,201)), ':', lw=fit_lw)

	elif msmts == '2':
		p0, fitfunc,fitfunc_str = Zfits.fit_2msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '3':
		p0, fitfunc,fitfunc_str = Zfits.fit_3msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '4':
		p0, fitfunc,fitfunc_str = Zfits.fit_4msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '5':
		p0, fitfunc,fitfunc_str = Zfits.fit_5msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '6':
		p0, fitfunc,fitfunc_str = Zfits.fit_6msmt_proc_fid(A0,offset0,t, p)

	### msmts = 0 is an exception
	fixed = [1]
	if msmts =='0':
		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
	else:
		fixed = [0,1,2] ### fixed parameter: [0,1,2] --> fix decay time, offset and amplitude, p is the only free parameter. 
										###[0,1] --> fix the amplitude and the offset for 0 measurements only.

		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)

	p1 = str(round(fit_result['params'][-1]*100,1))
	p1_u = str(round(fit_result['error'][-1]*100,1))

	result_string = p1 + ' +- ' + p1_u 

	return fit_result, result_string


def fit_pheno_State_decay(msmts,ax,A0,evotime,fid, t = 8.25,contrast = False):
	"""
	takes a zeno data set for a specific state and number of measurements and returns the fitted result.
	Inputs:

	msmts 			string which signifies the number of measurements
	A0 				the amplitude for the 0 measurement case.
	evotime 		np.array with the evolution times of the data points
	fid				np.array with the measured state fidelities.

	Output:

	result_string 	a string which is used for labelling the fits
	fit_result		the fitted function for plotting.
	"""

	t1 = 9./np.sqrt(2)
	# t2 = 8.25
	t2 = t
	p = 0.25
	repump = 0.95

	# print 'A0 ', A0
	# print 't2 ', t2


	if msmts == '0'or msmts == '00':
		if contrast:
			p0, fitfunc, fitfunc_str = common.fit_gauss(0.0, 0.43, 0., t)
		else:
			p0, fitfunc, fitfunc_str = common.fit_gauss(0.5, 0.43, 0., t)

	elif msmts == '1'or msmts == '01':
		### we usually do not take repumping errors into account for this fit.
		### therefore t1 is disregarded.
		p0, fitfunc,fitfunc_str = Zfits.fit_1msmt_state_fid(A0,t1,t2, p,repump,False,contrast = contrast)

		### manual option to show the intial guess in the plot.
		if False:
			ax.plot(np.linspace(0.,120.0,201), fitfunc(np.linspace(0.,120.0,201)), ':', lw=fit_lw)

	elif msmts == '2'or msmts == '02':
		p0, fitfunc,fitfunc_str = Zfits.fit_2msmt_state_fid(A0,t, p,contrast = contrast)

	elif msmts == '3'or msmts == '03':
		p0, fitfunc,fitfunc_str = Zfits.fit_3msmt_state_fid(A0,t, p,contrast = contrast)

	elif msmts == '4'or msmts == '04':
		p0, fitfunc,fitfunc_str = Zfits.fit_4msmt_state_fid(A0,t, p,contrast = contrast)

	elif msmts == '5'or msmts == '05':
		p0, fitfunc,fitfunc_str = Zfits.fit_5msmt_state_fid(A0,t, p,contrast = contrast)

	elif msmts == '6'or msmts == '06':
		p0, fitfunc,fitfunc_str = Zfits.fit_6msmt_state_fid(A0,t, p,contrast = contrast)

	elif msmts == '8'or msmts == '08':
		p0, fitfunc,fitfunc_str = Zfits.fit_8msmt_state_fid(A0,t, p,contrast = contrast)

	elif msmts == '10':
		p0, fitfunc,fitfunc_str = Zfits.fit_10msmt_state_fid(A0,t, p,contrast = contrast)

		if False:
			ax.plot(np.linspace(0.,120.0,201), fitfunc(np.linspace(0.,120.0,201)), ':', lw=fit_lw)

	elif msmts == '12':
		p0, fitfunc,fitfunc_str = Zfits.fit_12msmt_state_fid(A0,t, p,contrast = contrast)

	elif msmts == '16':
		p0, fitfunc,fitfunc_str = Zfits.fit_16msmt_state_fid(A0,t, p,contrast = contrast)

		

	### msmts = 0 is an exception
	fixed = [0,1]
	if msmts =='0':
		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
	else:
		fixed = [0] ### fixed parameter: [0,1] --> fix time and amplitude, p is the only free parameter. 
										###[0] --> fix the amplitude for 0 measurements only.

		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)

	# print fit_result['params']
	print 'params_dict ', fit_result['params_dict']
	print 'error_dict ', fit_result['error_dict']

	p1 = str(round(fit_result['params'][0]*100,3))
	p1_u = str(round(fit_result['error'][0]*100,3))
	# p1 = 'test'
	# p1_u = 'ing'
	result_string = p1 + ' +- ' + p1_u 

	return fit_result, result_string


def fit_pheno_process_decay(msmts,ax,A0,evotime,fid,decoded_bit,offset0 = 0.415, t = 8.2):
	"""
	takes a zeno data set for a specific number of measurements and returns the a fit to the data.
	Inputs:

	msmts 			string which signifies the number of measurements
	A0 				the amplitude for the 0 measurement case.
	evotime 		np.array with the evolution times of the data points
	fid				np.array with the measured state fidelities.

	Output:

	result_string 	a string which is used for labelling the fits
	fit_result		the fitted function for plotting.
	"""
	### decoding to carbon 2
	if decoded_bit == '2':
		# t = 8.34
		p = 0.08
		# offset0 = 0.40

	### decoding to carbon 1
	elif decoded_bit == '1':
		# t = 5.81
		p = 0.08
		# offset0 = 0.430

	if msmts == '0':
		p0, fitfunc, fitfunc_str = common.fit_gauss(0.4, 0.43, 0., t)


	elif msmts == '1':
		p0, fitfunc,fitfunc_str = Zfits.fit_1msmt_proc_fid(A0,offset0,t, p)

		### manual option to show the intial guess in the plot.
		if False:
			ax.plot(np.linspace(0.,120.0,201), fitfunc(np.linspace(0.,120.0,201)), ':', lw=fit_lw)

	elif msmts == '2':
		p0, fitfunc,fitfunc_str = Zfits.fit_2msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '3':
		p0, fitfunc,fitfunc_str = Zfits.fit_3msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '4':
		p0, fitfunc,fitfunc_str = Zfits.fit_4msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '5':
		p0, fitfunc,fitfunc_str = Zfits.fit_5msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '6':
		p0, fitfunc,fitfunc_str = Zfits.fit_6msmt_proc_fid(A0,offset0,t, p)

	### msmts = 0 is an exception
	fixed = [1]
	if msmts =='0':
		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
	else:
		fixed = [0,1] ### fixed parameter: [0,1,2] --> fix decay time, offset and amplitude, p is the only free parameter. 
										###[0,1] --> fix the amplitude and the offset for 0 measurements only.

		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)


	print 'params_dict ', fit_result['params_dict']
	print 'error_dict ', fit_result['error_dict']

	p1 = str(round(fit_result['params'][-1]*100,1))
	p1_u = str(round(fit_result['error'][-1]*100,1))

	result_string = p1 + ' +- ' + p1_u 

	return fit_result, result_string



########################

""" Carbon RO correction """

########################

def OneCarbonErrorpropagation(c_val,c_u,n_pop,n_pop_u):
	ret_val = np.sqrt(c_val/n_pop) ## value
	ret_u = 0.5*np.sqrt(c_u**2/(c_val*n_pop)+(n_pop_u**2)*c_val/n_pop**3) ## uncertainty

	return ret_val , ret_u

def TwoCarbonErrorpropagation(c12,c12_u,c1,c1_u,c2,c2_u,n_pop,n_pop_u):
	ret_val = np.sqrt(c12/(n_pop*c1*c2))	## value


	u_denominator = c1**2 * c2**2 * n_pop **2
	
	u_nominator = c1**2 * c2**2 * n_pop**2 * c12_u**2 + \
				c12**2 * (c1**2 * n_pop**2 * c2_u**2 + \
					c2**2 * (n_pop*c1_u**2 + c1**2 * n_pop_u**2))
	
	ret_u = np.sqrt(u_nominator)/u_denominator ## uncertainty

	return ret_val , ret_u


def ThreeCarbonError(c123,c123_u,c1,c1_u,c2,c2_u,c3,c3_u,n_pop,n_pop_u):
	ret_val = np.sqrt(c123/(n_pop*c1*c2*c3))

	u_denominator = c1**2 * c2**2 * c3**2 * n_pop**2

	u_nominator = u_denominator * c123_u**2 + c123**2 * (c1**2 * c3**2 * n_pop**2 * c2_u**2 + \
				c2**2 *(n_pop**2 * (c3**2 * c1_u**2 + c1**2 * c3_u**2) + \
				c1**2 * c3**2 * n_pop_u**2))

	ret_u = np.sqrt(u_nominator)/u_denominator

	return ret_val, ret_u



def calc_ro_correction(dict_val,dict_err,carbon_list = [1,2,5]):

	SolveDict = {}
	SolveDict['values'] = {}
	SolveDict['error'] = {}

	n_pop = dict_val['nitrogen']
	n_pop_u = dict_err['nitrogen']

	### go through all single carbon cases.
	for c in carbon_list:
		SolveDict['values'][str(c)],SolveDict['error'][str(c)] = OneCarbonErrorpropagation(dict_val[str(c)],dict_err[str(c)],n_pop,n_pop_u)

	### with these values compute the RO correction for two carbon values.

	for c in carbon_list:
		carbon_cross = cp.deepcopy(carbon_list)
		carbon_cross.remove(c)

		c1 = carbon_cross[0]
		c2 = carbon_cross[1]
		key_string = str(c1)+str(c2)

		SolveDict['values'][key_string], SolveDict['error'][key_string] = TwoCarbonErrorpropagation(
																						dict_val[key_string],
																						dict_err[key_string],
																						SolveDict['values'][str(c1)],
																						SolveDict['error'][str(c1)],
																						SolveDict['values'][str(c2)],
																						SolveDict['error'][str(c2)],
																						n_pop,
																						n_pop_u)

	### calculate the three carbon case

	key_string = ''.join(str(c) for c in carbon_list)
	print key_string
	SolveDict['values'][key_string], SolveDict['error'][key_string] = ThreeCarbonError(
																						dict_val[key_string],
																						dict_err[key_string],
																						SolveDict['values'][key_string[0]],
																						SolveDict['error'][key_string[0]],
																						SolveDict['values'][key_string[1]],
																						SolveDict['error'][key_string[1]],
																						SolveDict['values'][key_string[2]],
																						SolveDict['error'][key_string[2]],
																						n_pop,
																						n_pop_u)
	
	#### print the RO correction values
	return SolveDict

#### these values have been obtained via tomography/DESR measurements. see onenote 2015-06-09 NK.


### the derived notebook is used in the function Zeno_get_2Q_values

Dict = {'values' : {} , 'error' : {}}

Dict['values']['1']	=	0.859854829704
Dict['values']['2']	=	0.846657530075
Dict['values']['5']	=	0.88986599665
Dict['values']['12']	=	0.734843046401
Dict['values']['15']	=	0.784368785197
Dict['values']['25']	=	0.747335160652
Dict['values']['125']	=	0.666209837064
Dict['values']['nitrogen'] = 0.963484593103#0.956397901621
Dict['error']['1']	=	0.00621794744892
Dict['error']['2']	=	0.00775792880472
Dict['error']['5']	=	0.00735657065159
Dict['error']['12']	=	0.00565949773662
Dict['error']['15']	=	0.00217798138803
Dict['error']['25']	=	0.00219162359014
Dict['error']['125']	=	0.0089820755304
Dict['error']['nitrogen'] = 0.0149696842867#0.018270745207


CarbonRO_Dict = calc_ro_correction(Dict['values'],Dict['error'])

print CarbonRO_Dict


########################

""" Data acquisition """

########################

def get_Zeno_data(electron_RO=['positive'], 
				msmts='0',
				state='Z',ROBasis='IX',previous_evo=None, older_than=None,
				ssro_timestamp=None,
				single_qubit=False,
				subtract_drive_time = False,
				return_drive_time = False):
	"""
	this function finds timestamps according to an input which specifies the folder name
	Input: 	electron_RO 		is a list with e.g. positive or negative. 
								If the list is longer than 1 element then the 
								contrast values are returned.
			msmts & state 		both strings are used to find the correct folder
			older_than 			Timestamp which gives the starting point for the search 
			previous_evo		Gives the last free evolution time of the previously extracted data

			Flags for the search string: msmts

	Output: timestamp			specifies the latest evaluated folder
			loop_bit			Boolean which signals the end of data acquisition and if the output data should be evaluated.
			x_labels 			Tomographic bases. E.g. XI or ZZ
			y 					List, Read-out results (AVERAGED CONTRAST)
			y_err 				List, Read-out uncertainty
			evo_time 			List, Evolution time associated with y
			folder 				The folder of the latest timestamp
			evo_was_zero		Boolean which signifies if the first evolution time was Zero. This will then stop further data acquisition.
	"""

	if not single_qubit:
		search_string=electron_RO[0]+'_logicState_'+state+'_'+msmts+'msmt__ROBasis_'+ROBasis
	
	else:# adjust the measurement string for single qubit msmts. (no entanglement)
		search_string=electron_RO[0]+'_logicState_'+state+'_'+msmts+'msmt_singleQubit_ROBasis_'+ROBasis
	loop_bit = False
	evo_was_zero = False

	timestamp=None

	if previous_evo == None:
		previous_evo = 2000 #very high value such that all evolution times are taken into account.

	#if the desired data set exists, then read the measured values.
	if toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=older_than,
									newer_than=None,
									raise_exc=False) != False and previous_evo!=0:

		timestamp,folder=toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=older_than,
									newer_than=None,
									raise_exc=False)

		evotime,y,y_err= Zeno_get_2Q_values(timestamp,ssro_calib_timestamp=ssro_timestamp)

		x_labels = folder[folder.find('ROBasis_')+8:]
	

	else: #if the searched data does not exist --> dummy results. evotime has to be larger than previous_evo

		x_labels,y,y_err,evotime= 'II',0,0,[2001]


	if evotime[0] < previous_evo or previous_evo==None:
		loop_bit = True
		
	else:
		print 'wrong time assignment'
		x_labels,y,y_err,evotime= 'II',0,0,[2001]


	#if positive and negative RO are considered then adjust the search string and search for the same evo time.
	
	if len(electron_RO)>1 and (evotime[0] < previous_evo or previous_evo==None) and loop_bit:
		if not single_qubit:
			search_string=electron_RO[1]+'_logicState_'+state+'_'+msmts+'msmt__ROBasis_'+str(ROBasis)
		else:
			search_string=electron_RO[1]+'_logicState_'+state+'_'+msmts+'msmt_singleQubit_ROBasis_'+str(ROBasis) # adjust the measurement string for single qubit values.
		timestamp2,folder2=toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=older_than,
									newer_than=timestamp, #choose the right search direction.
									raise_exc=False)
		if evotime[0] < previous_evo or previous_evo==None:
			loop_bit = True

			evotime2,y2,y_err2= Zeno_get_2Q_values(timestamp2,ssro_calib_timestamp=ssro_timestamp)
		else:
			x_labels,y2,y_err2,evotime2= 'II',0,0,[2001]

		if electron_RO[0]== 'positive':	
			for i in range(len(y)):
				y[i]=(y[i]-y2[i])/2
				y_err[i]=((y_err[i]**2+y_err2[i]**2)**0.5)/2

		if electron_RO[0] == 'negative':
			for i in range(len(y)):
				y[i]=(-y[i]+y2[i])/2
				y_err[i]=((y_err[i]**2+y_err2[i]**2)**0.5)/2

	#if the evotime is set to zero then the data point needs to be shifted by the length of the applied parity measurements.
	#this duration is stored in the dictionary of the hdf group and is extracted here.
	if evotime[0]==0:
		# open the file, the folder/folder combination gives the complete filepath. 'r' for read only.
		Datafile=h5py.File(folder+folder[26:] + '.hdf5','r') 
		# you need to open the first datagroup. this can be accessed like a dictionary. 
		# that means file.keys[0] returns the first keyword of the dictionary
		first_dict_layer=Datafile[Datafile.keys()[0]]
		#the contents of msmt_params are dumped into the attributes object of this layer.
		#Attributes are again organized in dictionaries. --> Extract the desired value with a keyword
		if u'parity_duration' in first_dict_layer.attrs.keys():
			evotime[0]=first_dict_layer.attrs['parity_duration']*1e3
		Datafile.close()
		evo_was_zero=True

	### subtract the average drive time from the evolution time
	### 20150626 addendum, we only take the driving time of the qubits involved in the expectation value.
	subtract = 0
	if subtract_drive_time:

		Datafile=h5py.File(folder+folder[26:] + '.hdf5','r')

		first_dict_layer =  Datafile[Datafile.keys()[0]]
		c_list = first_dict_layer.attrs['carbon_list']
		parduration = first_dict_layer.attrs['parity_duration']
		no_of_carbons = float(len(c_list))
		repumper = first_dict_layer.attrs['Repump_duration']
		no_of_msmts = first_dict_layer.attrs['Nr_Zeno_parity_msmts']



		c_RO_list = [] ### will be filled with the carbons we perform tomography (X and Y valyes count) on.

		for index, t in enumerate(first_dict_layer.attrs['Tomography Bases']):
			
			if t in 'XY': c_RO_list.append(c_list[index]) ### this carbon was analysed.


		### we only did Z measurements (or <Identity> which would not make any sense)
		if c_RO_list == []:
			for index, t in enumerate(first_dict_layer.attrs['Tomography Bases']):
				if t in 'Z': c_RO_list.append(c_list[index])



		# summed_drive_time = 0
		
		# for c in c_RO_list:
		# 	### add gate times for each carbon
		# 	tau = first_dict_layer.attrs['C' + str(c) + '_Ren_tau'][0]
		# 	C_N = first_dict_layer.attrs['C' + str(c) + '_Ren_N'][0]

		#	summed_drive_time += 2*tau*C_N*no_of_msmts

		if no_of_msmts != 0: ### in this case parity_duration should also be !=0. otherwise division by 0.
			subtract = (parduration-repumper*no_of_msmts)*1e3
			subtract = subtract/no_of_carbons ## subtract AVERAGE driving time of the carbons.
			evotime = [x-subtract for x in evotime]
		else:
			subtract = 0
		Datafile.close()
		# print 'summed drive time ',summed_drive_time/len(c_RO_list)
		# print 'subtract ',subtract

	#determine the older timestamp (for the two eRO possiblities) and return that one. Check also whether or not to subtract the driving time.
	if loop_bit:
		if len(electron_RO)==1:
			if electron_RO[0]=='negative':
				if return_drive_time:
					return timestamp,loop_bit,x_labels,-1*y,y_err,evotime,folder,evo_was_zero, subtract
				else:
					return timestamp,loop_bit,x_labels,-1*y,y_err,evotime,folder,evo_was_zero
			else:
				if return_drive_time:
					return timestamp,loop_bit,x_labels,y,y_err,evotime,folder,evo_was_zero, subtract
				else:
					return timestamp,loop_bit,x_labels,y,y_err,evotime,folder,evo_was_zero

		elif toolbox.is_older(timestamp,timestamp2):
			if return_drive_time:
				return timestamp,loop_bit,x_labels,y,y_err,evotime,folder,evo_was_zero,subtract
			else:
				return timestamp,loop_bit,x_labels,y,y_err,evotime,folder,evo_was_zero

		else: 
			if return_drive_time:
				return timestamp2,loop_bit,x_labels,y,y_err,evotime,folder2,evo_was_zero,subtract
			else:
				return timestamp2,loop_bit,x_labels,y,y_err,evotime,folder2,evo_was_zero
	else:
		if return_drive_time:
			return older_than,loop_bit,x_labels,y,y_err,evotime,toolbox.data_from_time(older_than),evo_was_zero,subtract
		else:
			return older_than,loop_bit,x_labels,y,y_err,evotime,toolbox.data_from_time(older_than),evo_was_zero


def analyze_tests(older_than=None,newer_than=None,ssro_timestamp=None,N=1,q3test = False):
	"""
	Analyzes the test states. In a time window which is specified by older_than and newer_than in timestamp format.
	Analyzes a maximum of N states.
	"""

	N=2*N #we always perform a contrast measurement: positiv and negative --> Two folders.

	if q3test:
		search_string = '_1msmts_TESTSTATE_XXX'
		y_label = 'XXX contrast'
	else:
		search_string='_1msmts_TESTSTATE_ZZ'
		y_label = 'ZZ contrast'

	x_arr=[]
	y_arr=[]
	y_u_arr=[]
	old_Arr=[] #collection of relevant timestamps
	fig=plt.figure()
	ax=plt.subplot		

	for i in range(N):

		#if the desired data set exists, then read the measured values.
		if toolbox.latest_data(contains=search_string,
										return_timestamp =True,
										older_than=older_than,
										newer_than=newer_than,
										raise_exc=False) != False:

			older_than,folder=toolbox.latest_data(contains=search_string,
										return_timestamp =True,
										older_than=older_than,
										newer_than=newer_than,
										raise_exc=False)

			evotime,y,y_err= Zeno_get_2Q_values(older_than,ssro_calib_timestamp=ssro_timestamp)
			x_arr.append(i)
			y_arr.extend(y)
			y_u_arr.extend(y_err)
			if i%2==0: #append only one timestamp. we always plot the contrast of positive and negative.
				old_Arr.append(older_than)

	#condense the results for positive and negative read_out
	for i in range(len(x_arr)/2+1):
		if i==0:
			pass
		else:
			y_arr[i-1]=y_arr[2*i-1]/2.-y_arr[2*i-2]/2.
			y_u_arr[i-1]=np.sqrt(y_u_arr[2*i-1]**2+y_u_arr[2*i-2]**2)/2.

	y_arr=y_arr[:len(x_arr)/2]
	y_u_arr=y_u_arr[:len(x_arr)/2]
	x_arr=x_arr[:len(x_arr)/2]
	for i in range(len(y_arr)):
		plt.errorbar(x_arr[i],y_arr[i],y_u_arr[i],marker='o',label=old_Arr[i])
	plt.xlabel('timestamp')
	plt.ylabel(y_label)
	plt.title('test states')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.xlim(-0.5,N/2-0.5)
	print folder
	plt.savefig(os.path.join(folder,'Zeno_test_states.pdf'),format='pdf')
	plt.savefig(os.path.join(folder,'Zeno_test_states.png'),format='png')
	plt.show()
	plt.close('all')



def Zeno_get_2Q_values(timestamp=None, folder=None,folder_name='Zeno',
						measurement_name = ['adwindata'], 
						ssro_calib_timestamp =None):
	"""
	Returns the relevant CONTRAST values for a given timestamp.
	"""

	if timestamp == None and folder==None:
	    timestamp, folder   = toolbox.latest_data(folder_name,return_timestamp =True)
	elif timestamp ==None and folder!=None: 
		pass
	else:
	    folder = toolbox.data_from_time(timestamp)

	if folder != None and timestamp == None:
		d,t = toolbox.get_date_time_string_from_folder(folder)
		timestamp = toolbox.timestamp_from_datetime(t)
	### SSRO calibration
	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSRO',older_than=timestamp)


	else:
	    ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
	    ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'



	a = mbi.MBIAnalysis(folder)
	a.get_sweep_pts()
	a.get_readout_results(name='adwindata')
	try: 
		a.get_electron_ROC(ssro_calib_folder)
	except IOError:
		ssro.ssrocalib(ssro_calib_folder)
		a.get_electron_ROC(ssro_calib_folder)

	evo_time = a.sweep_pts.reshape(-1)
	y = ((a.p0.reshape(-1))-0.5)*2
	y_err = 2*a.u_p0.reshape(-1)

	if CarbonRO_correction:
		c_list = [] ### we will fill this up with the carbons which were read out

		Datafile=h5py.File(folder+folder[26:] + '.hdf5','r')

		first_dict_layer =  Datafile[Datafile.keys()[0]]

		carbons = first_dict_layer.attrs['carbon_list']
		tomos = first_dict_layer.attrs['Tomography Bases']

		#### sort carbons ascendingly. (makes it easier for the dictionary keys.)
		tomos = tomos[np.argsort(carbons)]
		carbons  =np.sort(carbons)
		
		for kk,t in enumerate(tomos):
			if t != 'I':
				c_list.append(carbons[kk])

		correction_key = ''
		for c in c_list:
			correction_key += str(c)


		### extract RO correection
		RO_corr =  CarbonRO_Dict['values'][correction_key]
		RO_corr_u = CarbonRO_Dict['error'][correction_key]
		
		### apply RO correction
		# print y_err
		y = y/RO_corr
		y_err = np.sqrt(RO_corr**2*y_err**2+RO_corr_u**2*y**2)/RO_corr**2

		# print y_err

		Datafile.close()



	return evo_time,y,y_err


########################

""" Basic plotting & analysis """

########################

def Zeno_state_fidelity(older_than_tstamp=None,msmts='0',eRO_list=['positive','negative'],
								 	state='X',
								 	plot_results=True,decoded_bit='2',
								 	ssro_timestamp=None,single_qubit=False, single_qubit_ort = False,
								 	subtract_drive_time = False):
	"""
	Plots or returns the state fidelity for a decoded qubit as a function of time (one parity expectation value)
	"""

	loop_bit = True
	evo_time=None

	y_arr=[]
	y_err_arr=[]
	evo_time_arr=[]
	ii=0

	Tomo1Dict={'X':'ZZ','mX':'ZZ',
	    'Y':'YZ',
	    'mY':'YZ',
	    'Z':'XI',
	    'mZ':'XI'}

	Tomo2Dict={'X':'ZZ','mX':'ZZ',
	    'Y':'ZY',
	    'mY':'ZY',
	    'Z':'IX',
	    'mZ':'IX'}


	### for preserving the XX expectation value for several measurements.
	Tomo3Dict={'X':'XX','mX':'XX',
	    'Y':'XX',
	    'mY':'XX',
	    'Z':'XX',
	    'mZ':'XX'}

	Tomo4Dict={'X':'YY','mX':'YY',
	    'Y':'XX',
	    'mY':'XX',
	    'Z':'XX',
	    'mZ':'XX'}



	#### this dictionary is used if we want to compute the fidelity of the two qubit state!
	Tomo2qubitDict={'X':'YY','mX':'YY',
	    'Y':'YZ',
	    'mY':'YZ',
	    'Z':'XI',
	    'mZ':'XI'}

	if single_qubit:
		### choose orthogonal bases to see the effect of a detuning
		if single_qubit_ort:
			Tomo1Dict={'X':'ZI','mX':'ZI',
			    'Y':'YI',
			    'mY':'YI',
			    'Z':'XI',
			    'mZ':'XI'}
			Tomo2Dict={'X':'IZ','mX':'IZ',
			    'Y':'IX',
			    'mY':'IX',
			    'Z':'IY',
			    'mZ':'IY'}
		### the regular tomography dictionaries for the expected bases.
		else:
			Tomo1Dict={'X':'ZI','mX':'ZI',
			    'Y':'YI',
			    'mY':'YI',
			    'Z':'XI',
			    'mZ':'XI'}
			Tomo2Dict={'X':'IZ','mX':'IZ',
			    'Y':'IY',
			    'mY':'IY',
			    'Z':'IX',
			    'mZ':'IX'}

	RODict={'1':Tomo1Dict,'2':Tomo2Dict,'3':Tomo3Dict,'4':Tomo4Dict, '2qubit': Tomo2qubitDict}
	evo_time=[2005]
	while loop_bit:
		older_than_tstamp,loop_bit,x_labels,y,y_err,evo_time,folder,evo_was_zero=get_Zeno_data(electron_RO=eRO_list,
																					state=state,
																					older_than=older_than_tstamp,
																					previous_evo=evo_time[0],
																					msmts=msmts,
																					ssro_timestamp=ssro_timestamp,ROBasis=RODict[str(decoded_bit)][state],
																					single_qubit=single_qubit, subtract_drive_time = subtract_drive_time)
		#loop_bit is true as long as new data was found.
		if loop_bit:
			y_arr=np.concatenate((y_arr,y))
			y_err_arr=np.concatenate((y_err_arr,y_err))
			evo_time_arr=np.concatenate((evo_time_arr,evo_time))
			# if an evolution time of zero was reached --> stop the evaluation.
			if evo_was_zero:
				loop_bit=False



	#select the correct expectation value and the right sign for the contrast.
	sign=1

	if 'Y' in state:
		sign=-1
	elif state=='Z':
		sign=1
	elif state=='X':
		sign=1
	
	if 'm' in state:
		sign=-1*sign


	#### if we look at the XX expectation value then the decoded bit is '3'
	#### all input states have an expected positive contrast.
	if decoded_bit == '3':
		sign = 1

	#### if we look at the YY expectation value then the decoded bit is '3'
	#### all input states have an expected negative contrast.
	if decoded_bit == '4':
		sign = -1*sign


	### in the two qubit case one of the expectation values needs to be turned around.

	if decoded_bit == '2qubit' and 'X' in state:
		sign = -1*sign


	fid_arr=(sign*np.array(y_arr)+1)/2. ### convert expectation value to fidelity

	fid_u_arr=np.array(y_err_arr)/2. ### uncertainty is therefore cut into two as well.

	if len(eRO_list)==1:
		RO_String=eRO_list[0]
	else: RO_String = 'contrast'

	fid_arr=fid_arr[np.argsort(evo_time_arr)]
	fid_u_arr=fid_u_arr[np.argsort(evo_time_arr)]
	evo_time_arr=np.sort(evo_time_arr)


	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()

		plt.errorbar(evo_time_arr,fid_arr,fid_u_arr,color='blue',marker='o')
		plt.xlabel('Evolution time (ms)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('logicState_'+state+'_stop_'+str(older_than_tstamp)+'_'+RO_String)

		plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'_'+RO_String+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'_'+RO_String+'.png'),format='png')
		plt.show()
		plt.close('all')
		print folder
	else:
		return evo_time_arr,fid_arr,fid_u_arr,older_than_tstamp,folder,RODict[str(decoded_bit)][state]

def Zeno_proc_fidelity(msmts='0',older_than_tstamp=None,**kw):
	"""
	Plots the process fidelity for a decoded qubit as a function of time

	If plot_results = False:
		Returns the AVERAGE state fidelity!
	"""	

	ssro_timestamp 			= 	kw.pop('ssro_timestamp', None)
	plot_results 			= 	kw.pop('plot_results', True)
	single_qubit 			= 	kw.pop('single_qubit',False)

	decoded_bit 			= 	kw.pop('decoded_bit','2')
	eRO_list				=	kw.pop('eRO_list',['positive','negative'])

	subtract_drive_time 	= 	kw.pop('subtract_drive_time', False)
	only_entangled 			= 	kw.pop('only_entangled',False)


	state_list=['X','mX','Y','mY','Z','mZ']

	if only_entangled: #### consider only entangled state. 
	###Only makes sense if used as a subroutine where it returns 0.5+ <XXX>/2. where XXX is a chosen observable.
		state_list = ['X','mX','Y','mY']

	if decoded_bit == '4':
		state_list = ['X','mX']

	state_dict = {}

	for state in state_list:
		state_dict[state] = {}

	fid_arr=[]
	fid_u_arr=[]

	#get individual state fidelities
	for state in state_list:
		evo_time,fid,fid_u,tstamp,folder,ROtype=Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
										eRO_list=eRO_list,
										state=state, msmts=msmts,
										plot_results=False,decoded_bit=decoded_bit,
										ssro_timestamp=ssro_timestamp,single_qubit=single_qubit, subtract_drive_time = subtract_drive_time)
		if single_qubit and ('Y' in state or 'Z' in state):
			### get the orthogonal measurements in.
			evo_time_ort,fid_ort,fid_u_ort,tstamp,folder,ROtype=Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
											eRO_list=eRO_list,
											state=state, msmts=msmts,
											plot_results=False,decoded_bit=decoded_bit,
											ssro_timestamp=ssro_timestamp,single_qubit=single_qubit , single_qubit_ort=True)

			### get maximum fidelity.
			fid, fid_u = calculate_max_fidelity(fid, fid_u, fid_ort, fid_u_ort)

		state_dict[state][timetrace_keylist[0]],state_dict[state][timetrace_keylist[1]],state_dict[state][timetrace_keylist[2]] = evo_time, fid, fid_u

		### append the state fidelities to the fidelity arrays.
		fid_arr.append(fid);fid_u_arr.append(fid_u)

	#calculate average state fidelity
	avg_fid=np.zeros(len(fid_arr[0]))
	avg_fid_u=np.zeros(len(fid_u_arr[0]))
	for i in range(len(fid_arr[0])):
		for ii in range(len(fid_arr)):
			avg_fid[i]= avg_fid[i] + fid_arr[ii][i]/len(fid_arr)
			avg_fid_u[i]= avg_fid_u[i] + fid_u_arr[ii][i]**2/len(state_list)**2
		avg_fid_u[i]=avg_fid_u[i]**0.5


	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot		

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		plt.errorbar(evo_time,(3.*avg_fid-1.)/2.,1.5*avg_fid_u,color='blue',marker='o') ### for a derivation see Equation 15 in Gilchrist et al. PRA 71, 062310 (2005) & the cited references
		plt.xlabel('Evolution time (ms)')
		plt.ylabel('Process fidelity')
		plt.title('Process_fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)

		print folder
		plt.savefig(os.path.join(folder,'Zeno1QProcFid_decBit'+str(decoded_bit)+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QProcFid_decBit'+str(decoded_bit)+'.png'),format='png')
		plt.show()
		plt.close('all')

		fig=plt.figure()
		ax=plt.subplot()

		for i,state in enumerate(state_list):
			plt.errorbar(evo_time,fid_arr[i],fid_u_arr[i],marker='o', label=state)

		plt.xlabel('Evolution time (ms)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('State_fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)
		plt.legend()

		plt.savefig(os.path.join(folder,'Zeno1QStateFidelities_decBit'+str(decoded_bit)+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QStateFidelities_decBit'+str(decoded_bit)+'.png'),format='png')
		plt.show()
		plt.close('all')

	else:
		return evo_time,avg_fid,avg_fid_u,tstamp,folder,state_dict

def Zeno_proc_list(older_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive'],decoded_bit='2',ssro_timestamp=None,
						single_qubit=False,plot_results=True,
						fitting=False,
						subtract_drive_time = False,plot_title = '',results_save = False,
						save_name = 'pickle_dict'):

	"""
	Plots the process fidelities as a function of time for a list of msmts

	If plot_results == False --> Returns the arrays for evolution times, PROCESS fidelity and uncertainty.

	"""
	fid_arr,fid=[],[]
	fid_u_arr,fid_u=[],[]
	evotime_arr,evotime=[],[]

	pickle_dict = {}

	if len(msmt_list)==0:
		print 'nothing to do here'

	else:
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder,state_dict = Zeno_proc_fidelity(older_than_tstamp=older_than_tstamp,
									eRO_list=eRO_list, msmts=msmt_list[i],ssro_timestamp=ssro_timestamp,
									decoded_bit=decoded_bit, subtract_drive_time = subtract_drive_time,
									plot_results=False)
			evotime_arr.append(np.sort(evotime))
			fid_arr.append((3*fid[np.argsort(evotime)]-1)/2)
			fid_u_arr.append(1.5*fid_u[np.argsort(evotime)])


			### append results to pickle_dict
			pickle_dict[msmt_list[i]] = make_pickle_dict(evotime,(3*fid[np.argsort(evotime)]-1)/2.,1.5*fid_u[np.argsort(evotime)],timetrace_keylist)
			pickle_dict[msmt_list[i]]['state_dict'] = {}
			pickle_dict[msmt_list[i]]['state_dict'] = state_dict

			print folder
			print tstamp
		
		### rescale evolution times to ms if it is given in seconds.

		if evotime_arr[-1][-2] < 1.0:
			for kk, timings in enumerate(evotime_arr):
				new_evo = []
				for jj in timings:
					if jj > 1.0:
						new_evo.append(jj)
					else:
						new_evo.append(jj*1e3)

				evotime_arr[kk] = np.array(new_evo)

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		if plot_results:
			fig=plt.figure()
			ax=plt.subplot()
			
			if fitting:


				result = ['0']*len(msmt_list) ### prepare the result strings.


				#### prepare fit parameters. They depend on the bit you decode to.

				if decoded_bit == '2':
					amp0 = 0.402
					offset0 = 0.397

					t = 8.34
					p = 0.09

				elif decoded_bit == '1':
					amp0 = 0.430
					offset0 = 0.351

					t = 7.34
					p = 0.09
				
				for ii,msmts in enumerate(msmt_list):
					fit_result, result[ii] = fit_process_decay(msmts,ax,amp0,evotime_arr[ii],fid_arr[ii],decoded_bit)

					plot.plot_fit1d(fit_result, np.linspace(0.0,110.0,1001), ax=ax, plot_data=False,color = color_list[ii],add_txt = False, lw =fit_lw)

					if msmts == '0':
						result[ii] = ' p = 0'

				for i in range(len(msmt_list)):
					plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],fmt='o',markersize=4,label=str(msmt_list[i]) + ' : p = ' + result[i],color=color_list[i])

			else: ### no fitting involved.
				for i in range(len(msmt_list)):
					plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],fmt='o',markersize=4,label=str(msmt_list[i])+ ' msmts',color=color_list[i])

			if single_qubit: # adds the latest single qubit measurement to the data
				evotime_single,fid_single,fid_u_single,tstamp_1q,folder,state_dict = Zeno_proc_fidelity(older_than_tstamp=older_than_tstamp,
										eRO_list=eRO_list, msmts='0',ssro_timestamp=ssro_timestamp,decoded_bit=decoded_bit,
										plot_results=False,single_qubit=True)
				single_proc = (3*fid_single[np.argsort(evotime_single)]-1)/2
				single_proc_u = 1.5*fid_u_single[np.argsort(evotime_single)]
				evotime_single = np.array([t*1e3 for t in np.sort(evotime_single)]) ## rescale to ms

				pickle_dict['single_qubit'] = make_pickle_dict(evotime_single,single_proc,single_proc_u,timetrace_keylist)



				### plot data
				plt.errorbar(evotime_single,single_proc,single_proc_u,fmt='o',markersize=4,label='1 qubit',color=color_list[len(msmt_list)])
				
				### fit data with a gaussian.
				fit_result, result[ii] = fit_process_decay('0',ax,0.5,evotime_single,single_proc,decoded_bit)

				plot.plot_fit1d(fit_result, np.linspace(0.0,110.0,1001), ax=ax, plot_data=False,color =color_list[len(msmt_list)],add_txt = False, lw =fit_lw)

			plt.xlabel('Evolution time (ms)')
			plt.ylabel('Process fidelity')
			if plot_title == '':
				plot_title = 'Process fidelity'+'_stop_'+str(tstamp)+'_'+RO_String
			plt.title(plot_title)

			# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
			plt.legend()

			if results_save:
				# print pickle_dict
				print pickle_dict.keys()
				save_data(pickle_dict,save_name+".p")
				folder = r'D:\measuring\data\Zeno_results'
			print 'Plots are saved in:'
			print folder
			plt.savefig(os.path.join(folder,'ZenoProc_decBit'+str(decoded_bit)+RO_String+'_combined.pdf'),format='pdf')
			plt.savefig(os.path.join(folder,'ZenoProc_decBit'+str(decoded_bit)+RO_String+'_combined.png'),format='png')
			plt.show()
			plt.close('all')
		else: 
			return evotime_arr,fid_arr,fid_u_arr,folder


def Zeno_state_list(older_than_tstamp=None,
						msmt_list=['0'],
						eRO_list=['positive'],
						state='Z',
						ssro_timestamp=None,
						decoded_bit='2',
						single_qubit=False,
						fitting = False,
						subtract_drive_time = False,
						physical_model_fit = False,
						pheno_phys_model_fit = False):
	
	fid=[]
	fid_u=[]
	evotime=[]
	if len(eRO_list)==0:
		print 'nothing to do here'

	else:
		evotime_arr = []
		fid_arr = []
		fid_u_arr = []

		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder,ROtype = Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
															msmts=msmt_list[i],
															eRO_list=eRO_list,decoded_bit=decoded_bit,
															subtract_drive_time = subtract_drive_time,
									plot_results=False,state=state,ssro_timestamp=ssro_timestamp)

			if sum(evotime)/len(evotime) < 1: ### scaling in seconds detected. Rescale to ms!
				for ii in range(len(evotime)-1):
					evotime[ii] = evotime[ii]*1e3

			evotime_arr.append(evotime[np.argsort(evotime)]); fid_arr.append(fid[np.argsort(evotime)]); fid_u_arr.append(fid_u[np.argsort(evotime)])
		
		if fitting:

			results = []
			results_u = []

			if physical_model_fit:
			### plot the fits to the data
			
				for kk,m in enumerate(msmt_list):
					if m == '0':
						if not 'X' in state:
							fit_result, result_string = fit_State_decay(m,ax,0.43,evotime_arr[kk],fid_arr[kk])
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							results.append('p = 0')
							
						else:
							### ZZ correlated state.
							p0, fitfunc,fitfunc_str = common.fit_poly([0])
							fit_result = fit.fit1d(evotime_arr[kk],fid_arr[kk], None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[])
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							results.append('p = 0')
							ampZ = fit_result['params'][0]*100
					else:
						if 'Z' in state:
							amp0 = (0.401+0.42) ### fitted amplitudes for mZ and Z
							fit_result, result_string = fit_State_decay(m,ax,amp0,evotime_arr[kk],fid_arr[kk])
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							results.append(result_string)
						elif 'Y' in state:
							amp0 = (0.386+0.340) ### fitted amplitudes for Y and mY
							fit_result, result_string = fit_State_decay(m,ax,amp0,evotime_arr[kk],fid_arr[kk])
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							results.append(result_string)

						elif 'X' in state:
							### fit an offset to the X state data and find the error probability due to our gates.
							if '0' not in msmt_list:
								ampZ = 100.

							p0, fitfunc, fitfunc_str = common.fit_line(ampZ/100.,0.0)
							fit_result = fit.fit1d(evotime_arr[kk],fid_arr[kk], None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=[1])
							plot.plot_fit1d(fit_result,np.linspace(0.0,110.0,len(fit_result['x'])), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)

							### calculate the error probability from the assumption that ampZ*(1-p)**m = ampNew
							### solving for p, the error probability, we get:

							p1 = str(round(100*(1-((fit_result['params'][0]*100-50)/(ampZ-50))**(1./int(m))),1))

							### calculate the error according to gaussian error propagation. See mathematica file in project folder.
							numerator = 2*(fit_result['error'][0]*np.abs(((fit_result['params'][0]*100-50.)/(ampZ-50.))**(1./int(m))))
							denominator = int(m)*np.abs(fit_result['params'][0]*2-1.)

							p1_u = str(round(numerator/denominator*100,1))

							result_string = p1 + ' +- ' + p1_u 
							results.append(result_string)


				for i in range(len(msmt_list)):
					plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],fmt='o',label=str(msmt_list[i])+ ' : ' + str(results[i]),color = color_list[i])

			elif pheno_phys_model_fit:
			### plot the fits to the data
				resultsSTR = []

				for kk,m in enumerate(msmt_list):
					if m == '0':
						if not 'X' in state:
							fit_result, result_string = fit_pheno_State_decay(m,ax,0.43,evotime_arr[kk],fid_arr[kk])
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							resultsSTR.append('p = 0')

							fit_result0 = fit_result
							# print fit_result0['params_dict']
							results.append(fit_result['params_dict']['sigma'])
							results_u.append(fit_result['error_dict']['sigma'])
							
						else:
							### ZZ correlated state.
							p0, fitfunc,fitfunc_str = common.fit_poly([0])
							fit_result = fit.fit1d(evotime_arr[kk],fid_arr[kk], None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[])
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							# results.append('p = 0')
							ampZ = fit_result['params'][0]*100
					else:

						if 'Z' in state:

							amp0 = fit_result0['params_dict']['A']*2 ### fitted amplitudes of the expectation value for mZ and Z
							fit_result, result_string = fit_pheno_State_decay(m,ax,amp0,evotime_arr[kk],fid_arr[kk])
							result_string =str(round(fit_result['params_dict']['t'],2)) + ' +- ' + str(round(fit_result['error_dict']['t'],2))
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							resultsSTR.append(result_string)

						elif 'Y' in state:

							amp0 = fit_result0['params_dict']['A']*2 ### fitted amplitudes of the expectation value for Y and mY
							fit_result, result_string = fit_pheno_State_decay(m,ax,amp0,evotime_arr[kk],fid_arr[kk])
							result_string =str(round(fit_result['params_dict']['t'],2)) + ' +- ' + str(round(fit_result['error_dict']['t'],2))
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							resultsSTR.append(result_string)
							print result_string


						else:
							print "this routine should not be used to analyze ZZ correlated states"

						results.append(fit_result['params_dict']['t'])
						results_u.append(fit_result['error_dict']['t'])

				for i in range(len(msmt_list)):
					plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],fmt='o',label=str(msmt_list[i])+ ' : ' + str(resultsSTR[i]),color = color_list[i])

			else:
				results = []
				results_u = []
				for kk,m in enumerate(msmt_list):
					### fit gaussian decays and plot a scaling law.
					p0, fitfunc, fitfunc_str = common.fit_gauss(0.5, 0.43, 0., 8.25)
					if kk == 0:
						print 'fit function: ', fitfunc_str

					fixed =[1]
					if int(msmt_list[kk]) % 2 ==0:
						print 'was here'
						fixed = [0,1]	
					else: 
						fixed = [1]

					fit_result = fit.fit1d(evotime_arr[kk],fid_arr[kk], None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
					plot.plot_fit1d(fit_result, np.linspace(0,evotime_arr[kk][-1],201), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)

					p1 = abs(round(fit_result['params'][3-len(fixed)],3))
					p1_u = round(fit_result['error'][3-len(fixed)],3)

					### fill up helper lists
					# fit_msmt_arr.append(p1)
					# fit_u_msmt_arr.append(p1_u)
					results.append(p1)
					results_u.append(p1_u)
					print 'uncertainty ', p1_u

					
					
					# print fit_msmt_arr
				for i in range(len(msmt_list)):
					plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],fmt='o',label=str(msmt_list[i])+' msmts, T = '+str(results[i])+' +- '+ str(results_u[i]),color = color_list[i])
					
					

			### plot the data itself and incorporate the fit results in the labels.
			
		else:
			### no fitting. plot the measured results
			for i in range(len(msmt_list)):
				plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],marker='o',label=str(msmt_list[i])+' msmts')

		#### treat the case of one qubit decaying.

		if single_qubit:
			evotime,fid,fid_u,tstamp,folder,ROtype = Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
															msmts='0',
															eRO_list=eRO_list,decoded_bit=decoded_bit,
									plot_results=False,state=state,ssro_timestamp=ssro_timestamp,single_qubit=True)
			### sort arrays.
			fid = fid[np.argsort(evotime)]
			fid_u = fid_u[np.argsort(evotime)]
			### plot
			plt.errorbar(np.array([t*1e3 for t in np.sort(evotime)]),fid,fid_u,marker='o',label='1 qubit')
			if 'Y' in state or 'Z' in state:
				### acquire the data for the orthogonal state.
				evotime_ort,fid_ort,fid_u_ort,tstamp,folder,ROtype = Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
																msmts='0',
																eRO_list=eRO_list,decoded_bit=decoded_bit,
										plot_results=False,state=state,ssro_timestamp=ssro_timestamp,single_qubit=True,single_qubit_ort = True)
				### sort arrays:
				fid_ort = fid_ort[np.argsort(evotime_ort)]
				fid_u_ort[np.argsort(evotime_ort)]
				plt.errorbar(np.array([t*1e3 for t in np.sort(evotime_ort)]),fid_ort,fid_u_ort,marker='o',label='1 qubit orthogonal')

				fid_bloch, fid_u_bloch = calculate_max_fidelity(fid,fid_u,fid_ort,fid_u_ort)

				### plot the maximum achievable fidelity (by quadratically adding the tomography results in X & Y)	
				plt.errorbar(np.array([t*1e3 for t in np.sort(evotime_ort)]),fid_bloch,fid_u_bloch,marker='o',label='1 qubit maximum fidelity')
		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		plt.xlabel('Evolution time (ms)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('logicState_'+state+'_stop_'+str(tstamp)+'_'+RO_String)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

		print 'Plots are saved in:'
		print folder
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.png'),format='png')
		plt.show()
		plt.close('all')

		if fitting and not physical_model_fit:
			x_axis = [int(x) for x in msmt_list]
			results_norm, results_u_norm = normalize_and_plot_coherence(x_axis,results,results_u)

			print results_norm
			print results_u_norm


			#### construct a pickle container out of the data
			pickle_dict = {	'msmts' 	: [int(x) for x in msmt_list],
							'results'	: results_norm,
							'results_u'	: results_u_norm}
			### save pickle file.
			save_data(pickle_dict,"Zeno_2Q_scaling_"+state+".p")

def Zen_Compare_Runs(msmt='0',older_timestamp_list=[None],eRO_list=['positive'],
						decoded_bit='2',ssro_timestamp_list=[None],
						single_qubit=False,plot_results=True):
	"""
	Takes a list of 'older_than' timestamps and searches for a run with the specified number of measurements.
	Raw data is then plotted.

	If plot_results = False: 	Returns the evolution time as a 1D list and a list of lists for the process fidelities and their uncertainties. 
								Additionally returns the latest folder.

	"""

	fid_arr,fid=[],[]
	fid_u_arr,fid_u=[],[]
	evotime_arr,evotime=[],[]


	for i in range(len(older_timestamp_list)):
		evotime,fid,fid_u,folder = Zeno_proc_list(older_than_tstamp=older_timestamp_list[i],
								eRO_list=eRO_list, msmt_list=[msmt],ssro_timestamp=ssro_timestamp_list[i],decoded_bit=decoded_bit,
								plot_results=False)
		evotime_arr.append(evotime[0])
		fid_arr.append(fid[0])
		fid_u_arr.append(fid_u[0])
	
	if len(eRO_list)==1:
		RO_String=eRO_list[0]
	else: RO_String = 'contrast'

	if plot_results:
		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(older_timestamp_list)):
			plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],marker='o',markersize=4,label=str(older_timestamp_list[i]))
			
		plt.xlabel('Evolution time (ms)')
		plt.ylabel('Process fidelity')
		plt.title('Process fidelity comparison'+'_'+RO_String+'_msmts_'+msmt)
		plt.legend()

		print 'Plots are saved in:'
		print folder
		plt.savefig(os.path.join(folder,'ZenoCompareRuns'+RO_String+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'ZenoCompareRuns'+RO_String+'.png'),format='png')
		plt.show()
		plt.close('all')

	### if no plot then return the results.
	else: return evotime_arr[0], fid_arr, fid_u_arr,folder


def Zeno_XX_list(older_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive','negative'],decoded_bit='3',ssro_timestamp=None,
						single_qubit=False,plot_results=True,fitting=True,pheno_phys_model_fit = False,
						subtract_drive_time = True,plot_title = '',results_save = False, save_name = 'pickle_dict'):

	"""
	Plots the XX Fidelity (averaged over all input states) for the chosen expectation values.

	decoded_bit == '3' selects XX expectation values

	Gets the values from the process fidelity routine. (This routine returns average fidelities.)

	"""
	fid_arr,fid=[],[]
	fid_u_arr,fid_u=[],[]
	evotime_arr,evotime=[],[]

	pickle_dict = {}

	if len(msmt_list)==0:
		print 'nothing to do here'

	else:
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder,state_dict = Zeno_proc_fidelity(older_than_tstamp=older_than_tstamp,
									eRO_list=eRO_list, msmts=msmt_list[i],ssro_timestamp=ssro_timestamp,
									decoded_bit=decoded_bit,
									plot_results=False, subtract_drive_time = subtract_drive_time)
			evotime_arr.append(np.sort(evotime))
			fid_arr.append(fid[np.argsort(evotime)])
			fid_u_arr.append(fid_u[np.argsort(evotime)])

			### make pickle dict with expectation values.
			pickle_dict[msmt_list[i]] = make_pickle_dict(evotime_arr[i],2*(np.array(fid_arr[i])-0.5),2*np.array(fid_u_arr[i]),timetrace_keylist)
			pickle_dict[msmt_list[i]]['state_dict'] = {}
			pickle_dict[msmt_list[i]]['state_dict'] = state_dict
		### rescale evolution times to ms if it is given in seconds.

		if evotime_arr[-1][-2] < 1.0: ### pick an entry in the middle of the evolution time list (do not accidentally pick 0..)
			for kk, timings in enumerate(evotime_arr):
				new_evo = []
				for jj in timings:
					if jj > 1.0:
						new_evo.append(jj)
					else:
						new_evo.append(jj*1e3)

				evotime_arr[kk] = np.array(new_evo)

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		expectation_arr = (np.array(fid_arr)-0.5)*2
		expectation_u_arr = np.array(fid_u_arr)*2

		if plot_results:
			fig=plt.figure()
			ax=plt.subplot()
			

			if fitting:
				results = []
				results_u = []
				if pheno_phys_model_fit:
				### plot the fits to the data
					resultsSTR = []

					for kk,m in enumerate(msmt_list):
						if m == '0':

							fit_result, result_string = fit_pheno_State_decay(m,ax,0.43,evotime_arr[kk],fid_arr[kk])
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							resultsSTR.append('p = 0')
							fit_result0 = fit_result
							# print fit_result0['params_dict']
							results.append(fit_result['params_dict']['sigma'])
							results_u.append(fit_result['error_dict']['sigma'])	

						else:

							amp0 = fit_result0['params_dict']['A']*2 ### fitted amplitudes of the expectation value for mZ and Z
							fit_result, result_string = fit_pheno_State_decay(m,ax,amp0,evotime_arr[kk],fid_arr[kk])
							result_string =str(round(fit_result['params_dict']['t'],2)) + ' +- ' + str(round(fit_result['error_dict']['t'],2))
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							resultsSTR.append(result_string)


							results.append(fit_result['params_dict']['t'])
							results_u.append(fit_result['error_dict']['t'])	

					for i in range(len(msmt_list)):
						plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],fmt='o',label=str(msmt_list[i])+ ' : ' + str(resultsSTR[i]),color = color_list[i])

				else:
					for kk,m in enumerate(msmt_list):
						### fit gaussian decays and plot a scaling law.
						p0, fitfunc, fitfunc_str = common.fit_gauss(0.0, 0.86, 0., 8.25)
						if kk == 0:
							print 'fit function: ', fitfunc_str

						fixed =[1]
						if int(msmt_list[kk]) % 2 ==0:
							print 'was here'
							fixed = [0,1]	
						else: 
							fixed = [1]

						fit_result = fit.fit1d(evotime_arr[kk],expectation_arr[kk], None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
						plot.plot_fit1d(fit_result, np.linspace(0,evotime_arr[kk][-1],201), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)

						p1 = abs(round(fit_result['params'][3-len(fixed)],3))
						p1_u = round(fit_result['error'][3-len(fixed)],3)

						### fill up helper lists
						# fit_msmt_arr.append(p1)
						# fit_u_msmt_arr.append(p1_u)
						results.append(p1)
						results_u.append(p1_u)
						print 'uncertainty ', p1_u

						
						
						# print fit_msmt_arr
					for i in range(len(msmt_list)):
						plt.errorbar(evotime_arr[i],expectation_arr[i],expectation_u_arr[i],fmt='o',label=str(msmt_list[i])+' msmts, T = '+str(results[i])+' +- '+ str(results_u[i]),color = color_list[i])

			else:
				for i in range(len(msmt_list)):
					plt.errorbar(evotime_arr[i],expectation_arr[i],expectaiton_u_arr[i],fmt='o',markersize=4,label=str(msmt_list[i])+ ' msmts',color = color_list[i])

			# if single_qubit: # adds the latest single qubit measurement to the data
			# 	evotime_single,fid_single,fid_u_single,tstamp,folder = Zeno_proc_fidelity(older_than_tstamp=older_than_tstamp,
			# 							eRO_list=eRO_list, msmts='0',ssro_timestamp=ssro_timestamp,decoded_bit=decoded_bit,
			# 							plot_results=False,single_qubit=True)
			# 	plt.errorbar(np.array([t*1e3 for t in np.sort(evotime_single)]),(3*fid_single[np.argsort(evotime_single)]-1)/2,1.5*fid_u_single[np.argsort(evotime_single)],marker='o',markersize=4,label='1 qubit')

			plt.xlabel('Evolution time (ms)')
			plt.ylabel('<XX>')
			if plot_title == '':
				plot_title = 'Average fidelity'+'_stop_'+str(tstamp)+'_'+RO_String
			plt.title(plot_title)
			plt.legend()

			print 'Plots are saved in:'
			print folder
			plt.savefig(os.path.join(folder,'ZenoProc_decBit'+str(decoded_bit)+RO_String+'_combined.pdf'),format='pdf')
			plt.savefig(os.path.join(folder,'ZenoProc_decBit'+str(decoded_bit)+RO_String+'_combined.png'),format='png')
			plt.show()
			plt.close('all')


			if results_save:
				save_data(pickle_dict,save_name+".p")


			if fitting:
				x_axis = [int(x) for x in msmt_list]
				results_norm, results_u_norm = normalize_and_plot_coherence(x_axis,results,results_u)


				#### construct a pickle container out of the data
				pickle_dict = {	'msmts' 	: [int(x) for x in msmt_list],
								'results'	: results_norm,
								'results_u'	: results_u_norm}
				### save pickle file.
				save_data(pickle_dict,"Zeno_2Q_scaling_XX.p")
		else: 
			return evotime_arr,fid_arr,fid_u_arr,folder


def Zeno_YY_list(older_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive','negative'],decoded_bit='4',ssro_timestamp=None,
						plot_results=True,fitting=True,pheno_phys_model_fit = False,
						subtract_drive_time = False,plot_title = ''):

	"""
	Plots the XX Fidelity (averaged over all input states) for the chosen expectation values.

	decoded_bit == '3' selects XX expectation values

	Gets the values from the process fidelity routine. (This routine returns average fidelities.)

	"""
	fid_arr,fid=[],[]
	fid_u_arr,fid_u=[],[]
	evotime_arr,evotime=[],[]
	if len(msmt_list)==0:
		print 'nothing to do here'

	else:
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder,state_dict = Zeno_proc_fidelity(older_than_tstamp=older_than_tstamp,
									eRO_list=eRO_list, msmts=msmt_list[i],ssro_timestamp=ssro_timestamp,
									decoded_bit=decoded_bit,
									plot_results=False, subtract_drive_time = subtract_drive_time)
			evotime_arr.append(np.sort(evotime))
			fid_arr.append(fid[np.argsort(evotime)])
			fid_u_arr.append(fid_u[np.argsort(evotime)])

		
		### rescale evolution times to ms if it is given in seconds.

		if evotime_arr[-1][-2] < 1.0: ### pick an entry in the middle of the evolution time list (do not accidentally pick 0..)
			for kk, timings in enumerate(evotime_arr):
				new_evo = []
				for jj in timings:
					if jj > 1.0:
						new_evo.append(jj)
					else:
						new_evo.append(jj*1e3)

				evotime_arr[kk] = np.array(new_evo)

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		expectation_arr = (np.array(fid_arr)-0.5)*2
		expectation_u_arr = np.array(fid_u_arr)*2

		if plot_results:
			fig=plt.figure()
			ax=plt.subplot()
			

			if fitting:
				results = []
				results_u = []

				if pheno_phys_model_fit:
				### plot the fits to the data
					resultsSTR = []

					for kk,m in enumerate(msmt_list):
						if m == '0':

							fit_result, result_string = fit_pheno_State_decay(m,ax,0.43,evotime_arr[kk],fid_arr[kk],contrast = False)
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							resultsSTR.append('p = 0')
							fit_result0 = fit_result
							# print fit_result0['params_dict']
							results.append(fit_result['params_dict']['sigma'])
							results_u.append(fit_result['error_dict']['sigma'])	

						else:

							amp0 = fit_result0['params_dict']['A']*2 ### fitted amplitudes of the expectation value for mZ and Z
							fit_result, result_string = fit_pheno_State_decay(m,ax,amp0,evotime_arr[kk],fid_arr[kk],contrast = False)
							result_string =str(round(fit_result['params_dict']['t'],2)) + ' +- ' + str(round(fit_result['error_dict']['t'],2))
							plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)
							resultsSTR.append(result_string)


							results.append(fit_result['params_dict']['t'])
							results_u.append(fit_result['error_dict']['t'])	
							
					for i in range(len(msmt_list)):
						plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],fmt='o',label=str(msmt_list[i])+ ' : ' + str(resultsSTR[i]),color = color_list[i])
				else:
					for kk,m in enumerate(msmt_list):
						### fit gaussian decays and plot a scaling law.
						p0, fitfunc, fitfunc_str = common.fit_gauss(0.0, 0.86, 0., 8.25)
						if kk == 0:
							print 'fit function: ', fitfunc_str

						fixed =[1]
						if int(msmt_list[kk]) % 2 ==0:
							print 'was here'
							fixed = [0,1]	
						else: 
							fixed = [1]

						fit_result = fit.fit1d(evotime_arr[kk],expectation_arr[kk], None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
						plot.plot_fit1d(fit_result, np.linspace(0,evotime_arr[kk][-1],201), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = fit_lw)

						p1 = abs(round(fit_result['params'][3-len(fixed)],3))
						p1_u = round(fit_result['error'][3-len(fixed)],3)

						### fill up helper lists
						# fit_msmt_arr.append(p1)
						# fit_u_msmt_arr.append(p1_u)
						results.append(p1)
						results_u.append(p1_u)
						print 'uncertainty ', p1_u
					
					# print fit_msmt_arr
					for i in range(len(msmt_list)):
						plt.errorbar(evotime_arr[i],expectation_arr[i],expectation_u_arr[i],fmt='o',label=str(msmt_list[i])+' msmts, T = '+str(results[i])+' +- '+ str(results_u[i]),color = color_list[i])

			else:
				for i in range(len(msmt_list)):
					plt.errorbar(evotime_arr[i],expectation_arr[i],expectaiton_u_arr[i],fmt='o',markersize=4,label=str(msmt_list[i])+ ' msmts',color = color_list[i])

			# if single_qubit: # adds the latest single qubit measurement to the data
			# 	evotime_single,fid_single,fid_u_single,tstamp,folder = Zeno_proc_fidelity(older_than_tstamp=older_than_tstamp,
			# 							eRO_list=eRO_list, msmts='0',ssro_timestamp=ssro_timestamp,decoded_bit=decoded_bit,
			# 							plot_results=False,single_qubit=True)
			# 	plt.errorbar(np.array([t*1e3 for t in np.sort(evotime_single)]),(3*fid_single[np.argsort(evotime_single)]-1)/2,1.5*fid_u_single[np.argsort(evotime_single)],marker='o',markersize=4,label='1 qubit')

			plt.xlabel('Evolution time (ms)')
			plt.ylabel('<YY>')
			if plot_title == '':
				plot_title = 'Average fidelity'+'_stop_'+str(tstamp)+'_'+RO_String
			plt.title(plot_title)
			plt.legend()

			print 'Plots are saved in:'
			print folder
			plt.savefig(os.path.join(folder,'YY_expectation.pdf'),format='pdf')
			plt.savefig(os.path.join(folder,'YY_expectation.png'),format='png')
			plt.show()
			plt.close('all')

			if fitting:
				x_axis = [int(x) for x in msmt_list]
				results_norm, results_u_norm = normalize_and_plot_coherence(x_axis,results,results_u)


				#### construct a pickle container out of the data
				pickle_dict = {	'msmts' 	: [int(x) for x in msmt_list],
								'results'	: results_norm,
								'results_u'	: results_u_norm}
				### save pickle file.
				save_data(pickle_dict,"Zeno_2Q_scaling_YY.p")
		else: 
			return evotime_arr,fid_arr,fid_u_arr,folder

def Zeno_2qubit_list(older_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive','negative'],ssro_timestamp=None,plot_results=True,
						subtract_drive_time = False,only_entangled = True):

	"""
	Plots the 2-qubit state fidelity (averaged over all input states)

	Gets the values from the process fidelity routine. (This routine returns average fidelities.)

	"""
	fid_arr,fid=[],[]
	fid_u_arr,fid_u=[],[]
	evotime_arr,evotime=[],[]


	#RODict={'1':Tomo1Dict,'2':Tomo2Dict,'3':Tomo3Dict, '2qubit': Tomo2qubitDict} --> dictionaries used in the state fidelity routine.

	decode_list = ['2','3','2qubit'] ### this list is needed to choose the corresponding read-out values from the state fidelity routine

	### Zeno_proc_fidelity returns average fidelities which get calculated from expectation values via F = 0.5+<observable>/2
	### need to be calculated backwards and obtain the fidelity from all 3 expectation values.
	### <obervable> = 2*(F-0.5)

	for ii in range(len(msmt_list)):
		for mm, decode in enumerate(decode_list):
			evotime,fid_int,fid_u_int,tstamp,folder,state_dict = Zeno_proc_fidelity(older_than_tstamp=older_than_tstamp,
									eRO_list=eRO_list, msmts=msmt_list[ii],ssro_timestamp=ssro_timestamp,decoded_bit=decode,
									plot_results=False, subtract_drive_time = subtract_drive_time,only_entangled = only_entangled)
			if mm == 0:
				fid = 2*(fid_int-0.5)
				print fid
				fid_u = (2*fid_u_int)**2/16 ### errors are summed quadratically and divided by 4**2 = 16

			else:
				fid = np.add(fid,2*(fid_int-0.5)) ### add absolute value of the expectation values together
				fid_u = np.add(fid_u,(2*fid_u_int)**2/16) ### add statistical errors together.


		#######################
		# calculate the 2qubit state fidelity from the accumulated average contrast values:
		#####################

		fid = 0.25*fid+0.25

		evotime_arr.append(np.sort(evotime))
		fid_arr.append(fid[np.argsort(evotime)])
		fid_u_arr.append(np.sqrt(fid_u[np.argsort(evotime)])) ### need to take sqrt for the errors.
	
	### rescale evolution times to ms if it is given in seconds.

	for kk, timings in enumerate(evotime_arr):
		new_evo = []
		for jj in timings:
			if jj > 1.0:
				new_evo.append(jj)
			else:
				new_evo.append(jj*1e3)

		evotime_arr[kk] = np.array(new_evo)

	if len(eRO_list)==1:
		RO_String=eRO_list[0]
	else: RO_String = 'contrast'

	if plot_results:
		fig=plt.figure()
		ax=plt.subplot()
		
		########################
		#### TODO fitting!!!####
		########################

		result = ['0']*len(msmt_list) ### prepare the result strings.


		#### prepare fit parameters. They depend on the bit you decode to.


		### first guesses
		amp0 = 0.502
		offset0 = 0.35
		t = 8.34

		
		for ii,msmts in enumerate(msmt_list):
			fit_result, result[ii] = fit_pheno_process_decay(msmts,ax,amp0,evotime_arr[ii],fid_arr[ii],'2', offset0 = offset0, t=t)

			plot.plot_fit1d(fit_result, np.linspace(0.0,110.0,1001), ax=ax, plot_data=False,color = color_list[ii],add_txt = False, lw =fit_lw)

			if msmts == '0':
				result[ii] = ' p = 0'
				### update parameters:
				amp0 = fit_result['params_dict']['A']
				offset0 = fit_result['params_dict']['a']
				t = fit_result['params_dict']['sigma']

		for i in range(len(msmt_list)):
			plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],fmt='o',markersize=4,label=str(msmt_list[i]) + ' : p = ' + result[i],color=color_list[i])

		plt.xlabel('Evolution time (ms)')
		plt.ylabel('Average state fidelity')
		plt.title('Average state fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)
		plt.legend()#bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.xlim([-2,60])
		plt.ylim([0.28,0.86])

		print 'Plots are saved in:'
		print folder
		plt.savefig(os.path.join(folder,'Zeno_2qubit_fidelity_'+RO_String+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno_2qubit_fidelity_'+RO_String+'.png'),format='png')
		plt.show()
		plt.close('all')
	else: 
		return evotime_arr,fid_arr,fid_u_arr,folder

def	Zeno_SingleQubit(older_than_tstamp=None,msmts='0',eRO_list=['positive'],
									 	state='X',
									 	plot_results=True,
									 	ssro_timestamp=None,subtract_drive_time = False,optimum_strategy = False):
		"""
		Plots or returns the state fidelity for a qubit as a function of time (one parity expectation value)
		"""

		RO_dict = {
		'X' : 'X',
		'mX' : 'X',
		'Y'  : 'Y',
		'mY'  : 'Y',
		'Z'  : 'Z',
		'mZ'  : 'Z',
		}

		loop_bit = True
		evo_time=None

		y_arr=[]
		y_err_arr=[]
		evo_time_arr=[]
		ii=0

		evotime_dict = {
		'0' : 0.0,
		'1' : 1.0,
		'2' : 6.0,
		'4' : 14.0,
		'6' : 16.0,
		'8' : 16.0,
		'10' : 18.0,
		'12' : 18.0,
		'16' : 18.0
		}

		evo_time=[2005]
		while loop_bit:
			older_than_tstamp,loop_bit,x_labels,y,y_err,evo_time,folder,evo_was_zero,subtracted_time=get_Zeno_data(electron_RO=eRO_list,
																						state=state,
																						older_than=older_than_tstamp,
																						previous_evo=evo_time[0],
																						msmts=msmts,
																						ssro_timestamp=ssro_timestamp,
																						ROBasis=RO_dict[state], 
																						subtract_drive_time = subtract_drive_time, return_drive_time = subtract_drive_time)
			#loop_bit is true as long as new data was found.
			if loop_bit:
				y_arr=np.concatenate((y_arr,y))
				y_err_arr=np.concatenate((y_err_arr,y_err))
				evo_time_arr=np.concatenate((evo_time_arr,evo_time))
				# if an evolution time of zero was reached --> stop the evaluation.
				if evo_was_zero or evotime_dict[msmts] in (evo_time+subtracted_time):
					loop_bit=False



		#select the correct expectation value and the right sign for the contrast.
		sign=1

		if 'Y' in state:
			sign=-1
		elif state=='Z':
			sign=1
		elif state=='X':
			sign=1
		
		if 'm' in state:
			sign=-1*sign

		if optimum_strategy and 'Y' in state:
			sign = -1*sign

		fid_arr=(sign*np.array(y_arr)+1)/2.
		fid_u_arr=np.array(y_err_arr)/2.

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		fid_arr=fid_arr[np.argsort(evo_time_arr)]
		fid_u_arr=fid_u_arr[np.argsort(evo_time_arr)]
		evo_time_arr=np.sort(evo_time_arr)


		if plot_results==True:
			fig=plt.figure()
			ax=plt.subplot()

			plt.errorbar(evo_time_arr,fid_arr,fid_u_arr,color='blue',marker='o')
			plt.xlabel('Evolution time (ms)')
			plt.ylabel('logical qubit state fidelity')
			plt.title('logicState_'+state+'_stop_'+str(older_than_tstamp)+'_'+RO_String)

			plt.savefig(os.path.join(folder,'Zeno1QDecay'+'_'+RO_String+'.pdf'),format='pdf')
			plt.savefig(os.path.join(folder,'Zeno1QDecay_'+RO_String+'.png'),format='png')
			plt.show()
			plt.close('all')
			print folder
		else:
			return evo_time_arr,fid_arr,fid_u_arr,older_than_tstamp,folder

def Zeno_1Q_msmt_list(older_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive','negative'],state='Z',
						ssro_timestamp=None,
						fitting = True, pheno_phys_model_fit = False,
						subtract_drive_time = False, physical_model_fit = False, plot_title = '',results_save = False):
	fid=[]
	fid_u=[]
	evotime=[]

	results = [] ## list that will collect the fit results.
	results_u = []
	resultsSTR = []


	pickle_dict = {} ## empty dictionary for data saving.

	if len(eRO_list)==0:
		print 'nothing to do here'

	else:
		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder = Zeno_SingleQubit(older_than_tstamp=older_than_tstamp,
															msmts=msmt_list[i],
															eRO_list=eRO_list,
									plot_results=False,state=state,ssro_timestamp=ssro_timestamp,subtract_drive_time = subtract_drive_time)
			if fitting:
				if physical_model_fit:
					### fit everything to a physical model!
					if msmt_list[i] == '0':
						t = 8.25
						p0, fitfunc, fitfunc_str = common.fit_gauss(0.5, 0.43, 0., t)
						fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,1])
						plot.plot_fit1d(fit_result, np.linspace(evotime[0],evotime[-1],1001), ax=ax,color = color_list[i], plot_data=False,add_txt = False, lw = fit_lw)
						results.append('0 : p = 0')

					else:
						fit_result, result_string = fit_State_decay(msmt_list[i],ax,0.415*2,evotime,fid)
						plot.plot_fit1d(fit_result, np.linspace(0.0,evotime[-1],1001), ax=ax, plot_data=False,color = color_list[i],add_txt = False, lw = fit_lw)
						results.append(result_string)

					plt.errorbar(np.sort(evotime),fid[np.argsort(evotime)],fid_u[np.argsort(evotime)],fmt='o',label=str(msmt_list[i])+' : p = ' + str(results[i]),color=color_list[i])
					
				if pheno_phys_model_fit:
				### plot the fits to the data
					

					if msmt_list[i] == '0':

						fit_result, result_string = fit_pheno_State_decay(msmt_list[i],ax,0.43,evotime,fid,contrast = False)
						plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[i],add_txt = False, lw = fit_lw)
						resultsSTR.append('p = 0')
						fit_result0 = fit_result
						# print fit_result0['params_dict']
						results.append(fit_result['params_dict']['sigma'])
						results_u.append(fit_result['error_dict']['sigma'])	

					else:

						amp0 = fit_result0['params_dict']['A']*2 ### fitted amplitudes of the expectation value for mZ and Z
						fit_result, result_string = fit_pheno_State_decay(msmt_list[i],ax,amp0,evotime,fid,contrast = False)
						result_string =str(round(abs(fit_result['params_dict']['t']),2)) + ' +- ' + str(round(fit_result['error_dict']['t'],2))
						plot.plot_fit1d(fit_result, np.linspace(0.0,evotime[-1],1001), ax=ax, plot_data=False,color = color_list[i],add_txt = False, lw = fit_lw)
						resultsSTR.append(result_string)


						results.append(abs(fit_result['params_dict']['t']))
						results_u.append(abs(fit_result['error_dict']['t']))


					plt.errorbar(evotime,fid,fid_u,fmt='o',label=str(msmt_list[i])+ ' : ' + str(resultsSTR[i]),color = color_list[i])

				else:
					### fit gaussian decays and plot a scaling law.
					p0, fitfunc, fitfunc_str = common.fit_gauss(0.5, 0.415, 0., 8.25)
					if i == 0:
						print 'fit function: ', fitfunc_str

					# fixed =[1]
					if int(msmt_list[i]) % 2 ==0:
						fixed = [0,1]
					fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
					plot.plot_fit1d(fit_result, np.linspace(0,evotime[-1],201), ax=ax, plot_data=False,color = color_list[i],add_txt = False, lw = fit_lw)

					p1 = abs(round(fit_result['params'][3-len(fixed)],3))
					p1_u = round(fit_result['error'][3-len(fixed)],3)

					### fill up helper lists
					# fit_msmt_arr.append(p1)
					# fit_u_msmt_arr.append(p1_u)
					results.append(p1)
					results_u.append(p1_u)
					print 'uncertainty ', p1_u
					# print fit_msmt_arr
					plt.errorbar(np.sort(evotime),fid[np.argsort(evotime)],fid_u[np.argsort(evotime)],fmt='o',label=str(msmt_list[i])+' msmts, T = '+str(results[i])+' +- '+ str(results_u[i]),color=color_list[i])
					

				

			else:
				plt.errorbar(np.sort(evotime),fid[np.argsort(evotime)],fid_u[np.argsort(evotime)],marker='o',label=str(msmt_list[i]),color=color_list[i])
			
			### fill up data dictionary
			if len(msmt_list[i]) == 1:
				pickle_dict['0'+str(msmt_list[i])] = make_pickle_dict(evotime,fid,fid_u,timetrace_keylist)
			else:
				pickle_dict[str(msmt_list[i])] = make_pickle_dict(evotime,fid,fid_u,timetrace_keylist)

		### fitting routine.
		
		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		plt.xlabel('Evolution time (ms)')
		plt.ylabel('F(|x>)')
		if plot_title == '':
			plt.title('State_'+state+'_stop_'+str(tstamp)+'_'+RO_String)
		else: plt.title(plot_title)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		# plt.legend()

		print 'Plots are saved in:'
		print folder

		filename = 'Zeno1QAvgDecays_'+RO_String+'_combined'

		if results_save:
			folder = r'D:\measuring\data\Zeno_results'
			save_data(pickle_dict,"X_preservation.p")

		plt.savefig(os.path.join(folder,filename + '.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,filename + '.png'),format='png')
		plt.show()
		plt.close('all')



		if fitting and not physical_model_fit:
			x_axis = [int(x) for x in msmt_list]
			results_norm, results_u_norm = normalize_and_plot_coherence(x_axis,results,results_u)

			#### construct a pickle container out of the data
			### we start off with a dictionary
			pickle_dict = {	'msmts' 	: [int(x) for x in msmt_list],
							'results'	: results_norm,
							'results_u'	: results_u_norm}
			save_data(pickle_dict,"Zeno_1Q_scaling_"+state+".p")


def Zeno_1Q_proc_fidelity(msmts='1',
								eRO_list=['positive','negative'],
								older_than_tstamp=None,
								plot_results=True,
								ssro_timestamp=None, subtract_drive_time = False, plot_title = '',results_save = False,optimum_strategy = False):
	
	"""
	Plots the process fidelity for a single qubit as a function of time

	If plot_Results = False:
		Returns the AVERAGE state fidelity!
	"""	
	state_list=['X','mX','Y','mY','Z','mZ']


	fid_arr=[]
	fid_u_arr=[]

	pickle_dict = {} ## empty dictionary for data saving.
	state_dict = {}
	for state in state_list:
		state_dict[state] = {}
	#get individual state fidelities
	for state in state_list:
		evo_time,fid,fid_u,tstamp,folder=Zeno_SingleQubit(older_than_tstamp=older_than_tstamp,
										eRO_list=eRO_list,
										state=state, msmts=msmts,
										plot_results=False,
										ssro_timestamp=ssro_timestamp, subtract_drive_time = subtract_drive_time, optimum_strategy = optimum_strategy)

		### append the state fidelities to the fidelity arrays.
		fid_arr.append(fid);fid_u_arr.append(fid_u)

		### append to state dict
		state_dict[state][timetrace_keylist[0]],state_dict[state][timetrace_keylist[1]],state_dict[state][timetrace_keylist[2]] = evo_time, fid, fid_u

	#calculate average state fidelity
	avg_fid=np.zeros(len(fid_arr[0]))
	avg_fid_u=np.zeros(len(fid_u_arr[0]))
	for i in range(len(fid_arr[0])): 	### sum over free evolution times
		for ii in range(len(fid_arr)): 	### sum over measured states.
			avg_fid[i]= avg_fid[i] + fid_arr[ii][i]/len(fid_arr)
			avg_fid_u[i]= avg_fid_u[i] + fid_u_arr[ii][i]**2/len(state_list)**2
		avg_fid_u[i]=avg_fid_u[i]**0.5


	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot		

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		plt.errorbar(evo_time,(3*avg_fid-1)/2.,1.5*avg_fid_u,color='blue',marker='o')

		pickle_dict[msmts] = make_pickle_dict(evo_time,(3*avg_fid-1)/2.,1.5*avg_fid_u,timetrace_keylist)


		plt.xlabel('Evolution time (ms)')
		plt.ylabel('Process fidelity')
		if plot_title == '':
			plot_title = 'Process_fidelity'+'_stop_'+str(tstamp)+'_'+RO_String

		plt.title(plot_title)

		if results_save:
			folder = r'D:\measuring\data\Zeno_results'
			

		plt.ylim(0.35,0.55)
		plt.xlim(-0.2,evo_time[-1]+2)

		plt.plot(np.linspace(0,evo_time[-1]+5,101),[0.5]*101,'--',lw=2,color='grey')


		print folder
		plt.savefig(os.path.join(folder,'Zeno1QProcFid_'+msmts+'msmts.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QProcFid_'+msmts+'msmts.png'),format='png')
		plt.show()
		plt.close('all')

		fig=plt.figure()
		ax=plt.subplot()

		### average the logical states with each other. I.e. mZ with Z etc.

		fid_state_avg = []
		fid_u_state_avg = []

		fid_state_avg,fid_u_state_avg = calc_state_average(fid_arr,fid_u_arr)

		### plot the states.

		for i in range(3):
			plt.errorbar(evo_time,fid_state_avg[i],fid_u_state_avg[i],fmt='o', label=['X','Y','Z'][i],color=color_list[i])


		for ii, state in enumerate(['X','Y','Z']):
			pickle_dict[msmts][state] = make_pickle_dict(evo_time,fid_state_avg[ii],fid_u_state_avg[ii],timetrace_keylist)
		
		if results_save:
			pickle_dict[msmts]['state_dict'] = {}
			pickle_dict[msmts]['state_dict'] = state_dict
			save_data(pickle_dict,"FIG1B_1Q_proc_fid.p")

		plt.xlabel('Evolution time (ms)')
		plt.ylabel('State fidelity')
		plt.xlim(-1,evo_time[-1]+2)
		plt.title('Individual state fidelities')
		plt.legend()



		plt.savefig(os.path.join(folder,'Zeno1QStateFidelities_'+msmts+'msmts.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QStateFidelities_'+msmts+'msmts.png'),format='png')
		plt.show()
		plt.close('all')

	else:
		return evo_time,avg_fid,avg_fid_u,tstamp,folder

def Zeno_3q_state_expectation(older_than_tstamp=None,msmts='0',eRO_list=['positive','negative'],
								 	state='00',RO_basis = 'XXX',
								 	plot_results=True,
								 	ssro_timestamp=None,
								 	subtract_drive_time = False):
	"""
	Plots or returns the 'absolute' value of the expectation value of an input state as a function of time.
	absolute value is not exactly correct here, we multiply every expectation value with Y with -1...
	"""

	loop_bit = True
	evo_time=None

	y_arr=[]
	y_err_arr=[]
	evo_time_arr=[]
	ii=0

	evo_time=[2005]
	drive_time = 0.
	while loop_bit:
		older_than_tstamp,loop_bit,x_labels,y,y_err,evo_time,folder,evo_was_zero,drive_time=get_Zeno_data(electron_RO=eRO_list,
																					state=state,
																					older_than=older_than_tstamp,
																					previous_evo=evo_time[0]+drive_time,
																					msmts=msmts,
																					ssro_timestamp=ssro_timestamp,ROBasis=RO_basis,
																					single_qubit=False,
																					subtract_drive_time = subtract_drive_time,
																					return_drive_time = True)
		#loop_bit is true as long as new data was found.
		if loop_bit:
			y_arr=np.concatenate((y_arr,y))
			y_err_arr=np.concatenate((y_err_arr,y_err))
			evo_time_arr=np.concatenate((evo_time_arr,evo_time))
			# if an evolution time of zero was reached --> stop the evaluation.
			if evo_was_zero:
				loop_bit=False


	#select the correct expectation value and the right sign for the contrast.
	sign=1

	if 'Y' in RO_basis:
		sign=-1

	exp_arr=sign*np.array(y_arr) 

	exp_u_arr=np.array(y_err_arr) 

	if len(eRO_list)==1:
		RO_String=eRO_list[0]
	else: RO_String = 'contrast'

	exp_arr=exp_arr[np.argsort(evo_time_arr)]
	exp_u_arr=exp_u_arr[np.argsort(evo_time_arr)]
	evo_time_arr=np.sort(evo_time_arr)

	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()

		plt.errorbar(evo_time_arr,exp_arr,exp_u_arr,color='blue',marker='o')
		plt.xlabel('Evolution time (ms)')
		plt.ylabel('Expectation value')
		plt.title('logicState_'+state+'_stop_'+str(older_than_tstamp)+'_basis_' + RO_basis)

		plt.savefig(os.path.join(folder,'Zeno3Q_msmt'+str(msmts)+'_basis_'+RO_basis + '.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno3Q_msmt'+str(msmts)+'_basis_'+RO_basis + '.png'),format='png')
		plt.show()
		plt.close('all')
		print folder
	else:
		return evo_time_arr,exp_arr,exp_u_arr,older_than_tstamp,folder

def Zeno_3q_state(older_than_tstamp=None,msmts='0',eRO_list=['positive','negative'],
								 	state='00',
								 	plot_results=True,
								 	ssro_timestamp=None,
								 	subtract_drive_time = False):
	"""
	Plots or returns the 'absolute' value of the expectation value of an input state as a function of time
	Plots all the 3 qubit expectation values for a specific state.
	"""

	State_dict = {
		'00'	:	['XIX','IXX','XXI','XXX'],
		'00p10' :	['IZZ','IXX','IYY','XXX'],
		'00p11' :	['XXI','YYI','ZZI','XXX']
	}

	loop_bit = True
	evo_time=None

	y_arr=[0]*len(State_dict[state])
	y_u_arr=[0]*len(State_dict[state])
	evo_time_arr=[0]*(len(State_dict[state]))

	evo_time=[2005]

	for kk,RO_basis in enumerate(State_dict[state]):
		evo_time,exp,exp_u,tstamp,folder = Zeno_3q_state_expectation(eRO_list = eRO_list,
																					state = state,
																					older_than_tstamp=older_than_tstamp,
																					msmts=msmts,
																					ssro_timestamp=ssro_timestamp,RO_basis=RO_basis,
																					plot_results=False,
																					subtract_drive_time = subtract_drive_time)
		y_arr[kk] = exp[np.argsort(evo_time)]
		y_u_arr[kk] = exp_u[np.argsort(evo_time)]
		evo_time_arr[kk] = np.sort(evo_time)




	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()

		for kk, label in enumerate(State_dict[state]):
			plt.errorbar(evo_time_arr[kk],y_arr[kk],y_u_arr[kk],marker='o',label = label)

		plt.xlabel('Evolution time (ms)')
		plt.ylabel('Measured Expectation value')
		plt.title('logicState_'+state+'_stop_'+str(tstamp))
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

		plt.savefig(os.path.join(folder,'Zeno3Q_msmt'+str(msmts) + '.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno3Q_msmt'+str(msmts) + '.png'),format='png')
		plt.show()
		plt.close('all')
		print folder

	else:
		return evo_time_arr,y_arr,y_u_arr,tstamp,folder

def Zeno_3q_state_list(older_than_tstamp=None,msmts_list=['0'],eRO_list=['positive','negative'],
								 	state='00',
								 	plot_results=True,
								 	ssro_timestamp=None,
								 	subtract_drive_time = False):
	"""
	Plots or returns the 'absolute' value of the expectation value of an input state as a function of time
	Takes a list for the number of zeno measurements and plots all data for that measurement for comparison.
	"""

	State_dict = {
		'00'	:	['XIX','IXX','XXI','XXX'],
		'00p10' :	['IZZ','IXX','IYY','XXX'],
		'00p11' :	['XXI','YYI','ZZI','XXX']
	}


	# State_dict = {
	# 	'00'	:	['XXX'],
	# 	'00p10' :	['IXX','XXX'],
	# 	'00p11' :	['XXX']
	# }


	loop_bit = True
	evo_time=None

	y_arr=[0]*len(msmts_list)
	y_u_arr=[0]*len(msmts_list)
	evo_time_arr=[0]*len(msmts_list)

	evo_time=[2005]

	for kk,msmts in enumerate(msmts_list):
		evo_time,exp,exp_u,tstamp,folder = Zeno_3q_state(eRO_list = eRO_list, state = state,
																older_than_tstamp=older_than_tstamp,
																msmts=msmts,
																ssro_timestamp=ssro_timestamp,
																plot_results=False,
																subtract_drive_time = subtract_drive_time)
		y_arr[kk] = exp
		y_u_arr[kk] = exp_u
		evo_time_arr[kk] = evo_time




	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()
		for ii, msmt in enumerate(msmts_list):
			for kk, label in enumerate(State_dict[state]):
				if label == 'XXX':
					plt.errorbar(evo_time_arr[ii][kk],y_arr[ii][kk],y_u_arr[ii][kk],fmt='o',label = msmt + ' msmt : '+ label)

		plt.xlabel('Evolution time (ms)')
		plt.ylabel('Measured Expectation value')
		plt.title('logicState_'+state+'_stop_'+str(tstamp))
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

		plt.savefig(os.path.join(folder,'Zeno3Q_msmt_list' + '.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno3Q_msmt_list' + '.png'),format='png')
		plt.show()
		plt.close('all')
		print folder

	else:
		return evo_time_arr,y_arr,y_u_arr,tstamp,folder


def Zeno_3q_logicalFidelity_list(older_than_tstamp=None,msmts_list=['0'],eRO_list=['positive','negative'],
								 	state_list=['00'],
								 	plot_results=True,
								 	ssro_timestamp=None,
								 	subtract_drive_time = False,fitting = False,pheno_phys_model_fit = False,
								 	plot_title = '',results_save = False, save_name = 'pickle_dict'):
	"""
	Plots or returns the logical 2 qubit state fidelity.
	Takes a list for the number of zeno measurements and plots all data for that measurement for comparison.
	"""

	State_dict = {
		'00'	:	['XIX','IXX','XXI','XXX'],
		'00p10' :	['IZZ','IXX','IYY','XXX'],
		'00p11' :	['XXI','YYI','ZZI','XXX']
	}

	State_labeling_dict = {'00' : r'|$0,0 \rangle$', '00p10' : r'|$X,0 \rangle$', '00p11' : r'|$\Phi^+ \rangle$'}

	### initialize arrays for final results.
	fid_arr = [[0]*len(msmts_list)]*len(state_list)
	fid_u_arr = [[0]*len(msmts_list)]*len(state_list)
	fid_evo_arr = [[0]*len(msmts_list)]*len(state_list)

	xxx_arr = [0]*len(msmts_list)
	xxx_u_arr = [0]*len(msmts_list)

	pickle_dict_fid = {}



	for jj,state in enumerate(state_list):
	
		### initialize data loop

		loop_bit = True
		evo_time=None

		y_arr=[0]*len(msmts_list)
		y_u_arr=[0]*len(msmts_list)
		evo_time_arr=[0]*len(msmts_list)

		evo_time=[2005] ### dummy evo time for get_data

		for kk,msmts in enumerate(msmts_list):
			evo_time,exp,exp_u,tstamp,folder = Zeno_3q_state(eRO_list = eRO_list, state = state,
																	older_than_tstamp=older_than_tstamp,
																	msmts=msmts,
																	ssro_timestamp=ssro_timestamp,
																	plot_results=False,
																	subtract_drive_time = subtract_drive_time)
			y_arr[kk] = exp
			y_u_arr[kk] = exp_u
			evo_time_arr[kk] = evo_time[kk]
			fid_evo_arr[jj][kk] = evo_time[kk]

		fid_msmt_arr = [0]*len(msmts_list)
		fid_u_msmt_arr = [0]*len(msmts_list)

		for mm in range(len(msmts_list)):

			fid_msmt_arr[mm] = [0.25 + 0.25*(y_arr[mm][0][nn]+y_arr[mm][1][nn]+y_arr[mm][2][nn]) for nn in range(len(y_arr[mm][0]))]
			fid_u_msmt_arr[mm] = [0.25*np.sqrt(y_u_arr[mm][0][nn]**2+y_u_arr[mm][1][nn]**2+y_u_arr[mm][2][nn]**2) for nn in range(len(y_arr[mm][0]))]

			if jj==0: 
				xxx_arr[mm] = [x/len(state_list) for x in y_arr[mm][3]]
				xxx_u_arr[mm] = [x**2 for x in y_u_arr[mm][3]]
				# print "i passed here!! if"
			else:
				xxx_arr[mm] = [xxx_arr[mm][ll] + x/len(state_list) for ll,x in enumerate(y_arr[mm][3])]
				xxx_u_arr[mm] = [xxx_u_arr[mm][ll] + x**2 for ll,x in enumerate(y_u_arr[mm][3])]
				# print "i passed here!! (else)"
		### values of the list pointer need to be deepcopied into the other array. 

		fid_arr[jj] = cp.deepcopy(fid_msmt_arr)
		fid_u_arr[jj] = cp.deepcopy(fid_u_msmt_arr)


		

	for jj in range(len(xxx_u_arr)):
		xxx_u_arr[jj] = [np.sqrt(x)/len(state_list) for x in xxx_u_arr[jj]] ### calculate the error on the average xxx value

	### only takes the fidelity of the last state into account
	for key in range(len(msmts_list)):
		pickle_dict_fid[msmts_list[key]] = {} ### dict in dict.

	print pickle_dict_fid
	for kk, state in enumerate(state_list):
		for jj, msmt in enumerate(msmts_list): 
			pickle_dict_fid[msmt][state] = make_pickle_dict(fid_evo_arr[kk][jj],fid_arr[kk][jj],fid_u_arr[kk][jj],timetrace_keylist)
			# print pickle_dict_fid
			pickle_dict_fid[msmt]['xxx'] = make_pickle_dict(fid_evo_arr[-1][jj],xxx_arr[jj],xxx_u_arr[jj],timetrace_keylist) ###overwrites a few times...
	
	if results_save:
		save_data(pickle_dict_fid,save_name+".p")



	if fitting:
		amp_dict = {}
		off_dict = {}
	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()

		fit_arr = [[0]*len(msmts_list)]*(len(state_list)+1)
		fit_u_arr = [[0]*len(msmts_list)]*(len(state_list)+1)

		for kk, state in enumerate(state_list):

			### initialize helper lists
			fit_msmt_arr = []
			fit_u_msmt_arr = []

			for ii, msmt in enumerate(msmts_list):

				if fitting:

					if pheno_phys_model_fit:
					### plot the fits to the data
						
						if msmt == '0':

							fit_result, result_string = fit_pheno_process_decay(msmt,ax,0.43,fid_evo_arr[kk][ii],fid_arr[kk][ii],decoded_bit='2')
							plot.plot_fit1d(fit_result, np.linspace(0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk*len(state_list)+ii],add_txt = False, lw = fit_lw)

							amp_dict[state] = fit_result['params_dict']['A']
							off_dict[state] = fit_result['params_dict']['a']

							p1 = fit_result['params_dict']['sigma']
							p1_u = fit_result['error_dict']['sigma']
							fit_msmt_arr.append(p1)
							fit_u_msmt_arr.append(p1_u)

						else:

							fit_result, result_string = fit_pheno_process_decay(msmt,ax,amp_dict[state],fid_evo_arr[kk][ii],fid_arr[kk][ii],decoded_bit='2',offset0=off_dict[state])
							result_string =str(round(fit_result['params_dict']['t'],2)) + ' +- ' + str(round(fit_result['error_dict']['t'],2))
							plot.plot_fit1d(fit_result, np.linspace(0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk*len(state_list)+ii],add_txt = False, lw = fit_lw)
							

							p1 = fit_result['params_dict']['t']
							p1_u = fit_result['error_dict']['t']
							fit_msmt_arr.append(p1)
							fit_u_msmt_arr.append(p1_u)	

						if len(state_list) != 1:
							label_string = State_labeling_dict[state]  + ' : '+ msmt + ' msmts. t = '+str(round(p1,2))+'+-'+str(round(p1_u,2))
						else: label_string = msmt + ' msmts. T = '+str(round(p1,2))+'+-'+str(round(p1_u,2))

						print ii
						print kk
						print len(state_list)

						plt.errorbar(fid_evo_arr[kk][ii],fid_arr[kk][ii],fid_u_arr[kk][ii],
														fmt 	= 'o',
														label 	= label_string,
														color 	= color_list[kk*len(state_list)+ii])

					else:
						p0, fitfunc, fitfunc_str = common.fit_gauss(0.25, 0.43, 0., 8)

						# if int(msmt) % 2 ==0:
						# 	fixed = [0,1]	
						# else: 
						fixed = [1]

						fit_result = fit.fit1d(fid_evo_arr[kk][ii],fid_arr[kk][ii], None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed= fixed)
						plot.plot_fit1d(fit_result, np.linspace(0,fid_evo_arr[kk][ii][-1],201), ax=ax, plot_data=False,color = color_list[kk*len(state_list)+ii],add_txt = False, lw =fit_lw)

						p1 = round(fit_result['params'][3-len(fixed)],3)
						p1_u = round(fit_result['error'][3-len(fixed)],3)

						### fill up helper lists
						fit_msmt_arr.append(p1)
						fit_u_msmt_arr.append(p1_u)

						print 'uncertainty ', p1_u
						print fit_msmt_arr
						if len(state_list) != 1:
							label_string = State_labeling_dict[state]  + ' : '+ msmt + ' msmts. T = '+str(round(p1,1))+'+-'+str(round(p1_u))
						else: label_string = msmt + ' msmts. T = '+str(round(p1,1))+'+-'+str(round(p1_u))

						plt.errorbar(fid_evo_arr[kk][ii],fid_arr[kk][ii],fid_u_arr[kk][ii],
														fmt 	= 'o',
														label 	= label_string,
														color 	= color_list[kk*len(state_list)+ii])
				else:
					plt.errorbar(fid_evo_arr[kk][ii],fid_arr[kk][ii],fid_u_arr[kk][ii],fmt='o',label = msmt + ' msmts_state: '+ state)

			### dump helper list into result array.

			fit_arr[kk] = cp.deepcopy(fit_msmt_arr)
			fit_u_arr[kk] = cp.deepcopy(fit_u_msmt_arr)
			
		#############
		# PLOT Logical 2Q STATE FIDELITIES!
		############

		plt.xlabel('Evolution time (ms)')
		plt.ylabel('State fidelity')
		if plot_title == '':
			plot_title = 'Logical 2-qubit state fidelity_stop_'+str(tstamp)
		if len(state_list) == 1:
			plot_title = plot_title + ' with ' + State_labeling_dict[state]
		plt.title(plot_title)
		if len(state_list)*len(msmts_list) < 4:
			plt.legend()
		else:
			plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

		plt.savefig(os.path.join(folder,'Zeno3Q_logical2QFidelity_msmt_list' + '.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno3Q_logical2QFidelity_msmt_list' + '.png'),format='png')
		plt.show()
		plt.close('all')


		#######

		# PLOT averaged XXX decays

		######

		fig=plt.figure()
		ax=plt.subplot()

		fit_msmt_arr = []
		fit_u_msmt_arr = []

		for ii in range(len(msmts_list)):


			if pheno_phys_model_fit:

			### plot the fits to the data
				
				if msmts_list[ii] == '0':

					# print xxx_arr[ii]
					fit_result, result_string = fit_pheno_State_decay(msmts_list[ii],ax,0.43,fid_evo_arr[0][ii],xxx_arr[ii],contrast = True,t=4.)
					plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[ii],add_txt = False, lw = fit_lw)

					fit_result0 = fit_result ### fitted amplitudes of the expectation value for mZ and Z
					# print fit_result0['params_dict']
					p1_xxx = fit_result['params_dict']['sigma']
					p1_u_xxx = fit_result['error_dict']['sigma']
					fit_msmt_arr.append(p1_xxx)
					fit_u_msmt_arr.append(p1_u_xxx)

				else:

					amp0 = abs(fit_result0['params_dict']['A'])
					fit_result, result_string = fit_pheno_State_decay(msmts_list[ii],ax,amp0,fid_evo_arr[0][ii],xxx_arr[ii],contrast = True)
					result_string =str(round(fit_result['params_dict']['t'],2)) + ' +- ' + str(round(fit_result['error_dict']['t'],2))
					plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[ii],add_txt = False, lw = fit_lw)
					

					p1_xxx = fit_result['params_dict']['t']
					p1_u_xxx = fit_result['error_dict']['t']
					fit_msmt_arr.append(p1_xxx)
					fit_u_msmt_arr.append(p1_u_xxx)	

				plt.errorbar(fid_evo_arr[0][ii] , xxx_arr[ii] , xxx_u_arr[ii],
												fmt 	= 'o',
												label 	= str(msmts_list[ii])+' msmts, T = '+str(round(p1_xxx,2))+' +- '+ str(round(p1_u_xxx,2)) ,
												color 	= color_list[ii])

			else:
				p0, fitfunc, fitfunc_str = common.fit_gauss(0.0, 0.43, 0., 8)

				if int(msmts_list[ii]) % 2 ==0:
					fixed = [0,1]	
				else: 
					fixed = [1]

				fit_result = fit.fit1d(fid_evo_arr[0][ii],xxx_arr[ii], None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
				plot.plot_fit1d(fit_result, np.linspace(0,fid_evo_arr[0][ii][-1],201), ax=ax, plot_data=False,color = color_list[ii],add_txt = False, lw = fit_lw)

				p1_xxx = round(fit_result['params'][3-len(fixed)],3)
				p1_u_xxx = round(fit_result['error'][3-len(fixed)],3)

				# ### fill up helper lists
				fit_msmt_arr.append(p1_xxx)
				fit_u_msmt_arr.append(p1_u_xxx)
				

				plt.errorbar(fid_evo_arr[0][ii] , xxx_arr[ii] , xxx_u_arr[ii],	
															fmt 	= 'o',
															label 	= str(msmts_list[ii])+' msmts, T = '+str(p1_xxx)+' +- '+ str(p1_u_xxx) , 
															color 	= color_list[ii])

		fit_arr[len(state_list)] = cp.deepcopy(fit_msmt_arr) ### remember that python starts counting at 0
		fit_u_arr[len(state_list)] = cp.deepcopy(fit_u_msmt_arr) ### same as above.

		plt.title('Preservation of <XXX>')
		plt.legend()
		plt.xlabel('Evolution time (ms)')
		plt.ylabel('<XXX>')
		plt.show()
		plt.close('all')


		###############

		# PLOT normalized coherence times

		##############
		fig=plt.figure()
		ax=plt.subplot()

		### plot the fitted coherence times
		x_axis = [int(x) for x in msmts_list]


		#### normalize the fitted coherence times

		### have to recalculate the uncertainties via error propagation first! Then normalise fitted coherence times.
		
		for kk in range(len(state_list)+1):
			normalisation = cp.deepcopy(fit_arr[kk][0])
			for ii in range(len(msmts_list)):
				fit_u_arr[kk][ii] = calc_error_propagation_fraction(fit_arr[kk][ii],normalisation,fit_u_arr[kk][ii],fit_u_arr[kk][0])

			fit_arr[kk] = [x/normalisation for x in cp.deepcopy(fit_arr[kk])]

		#### construct a pickle container out of the data
		### we start off with a dictionary
		pickle_dict = {	'msmts' 	: [x_axis]*(len(state_list)+1),
						'results'	: fit_arr,
						'results_u'	: fit_u_arr}
		### save fit values.
		if '1' not in msmts_list:
			save_data(pickle_dict,"Zeno_3Q_scaling.p")

		for kk, state in enumerate(state_list+['XXX']):
			plt.errorbar(x_axis,fit_arr[kk],fit_u_arr[kk],fmt = 'o',label = state)

		plt.xlabel('Number of Measurements')
		plt.ylabel('Normalized coherence time')
		plt.ylim(0,1.8)
		plt.xlim(-0.5,4.5)

		## put legend in place
		if len(state_list)*len(msmts_list) < 4:
			plt.legend()
		else:
			plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.show()
		plt.close('all')


		print 'Files are saved in: '
		print folder

	else:
		return evo_time_arr,y_arr,y_u_arr,tstamp,folder

def coherence_scaling():
	data2Q = []
	data2Q.append(pickle.load( open( "ZenData\Zeno_2Q_scaling_Z.p", "rb" ) ))
	data2Q.append(pickle.load( open( "ZenData\Zeno_2Q_scaling_mZ.p", "rb" ) ))
	data2Q.append(pickle.load( open( "ZenData\Zeno_2Q_scaling_Y.p", "rb" ) ))
	data2Q.append(pickle.load( open( "ZenData\Zeno_2Q_scaling_mY.p", "rb" ) ))
	data2Q.append(pickle.load( open( "ZenData\Zeno_2Q_scaling_XX.p", "rb" ) ))
	data2Q.append(pickle.load( open( "ZenData\Zeno_2Q_scaling_YY.p", "rb" ) ))

	data1Q = pickle.load( open( "ZenData\Zeno_1Q_scaling_X.p", "rb" ) )

	data3Q = pickle.load( open( "ZenData\Zeno_3Q_scaling.p", "rb" ) )

	dict_keys = ['msmts','results','results_u']


	

	combined_data = [data1Q]
	combined_data.extend(data2Q)


	for ii in range(len(data3Q['msmts'])):
		combined_data.append({'msmts':data3Q['msmts'][ii], 'results':data3Q['results'][ii],'results_u':data3Q['results_u'][ii]})
	


	# fig=plt.figure()
	# ax=plt.subplot()

	# labels1Q = ['1Q state |X>']
	# labels2Q = ['2Q Z_Logic','2Q mZ_Logic','2Q Y_Logic','2Q mY_Logic','2Q average XX']
	# labels3Q = ['3Q 00_Logic','3Q 00+10_Logic','3Q 00+11_Logic','3Q average XXX']
	# labels = labels1Q + labels2Q + labels3Q


	# for ii,data in enumerate(combined_data):
	# 	print ii
	# 	plt.errorbar(data[dict_keys[0]],data[dict_keys[1]],data[dict_keys[2]],fmt='o',label =labels[ii],color = color_list[ii])

	# def f(N): return (1.*N)**0.44
	
	# f = np.vectorize(f)

	# print f(0.)

	# plt.plot(np.linspace(0,15,1001),f(np.linspace(0,15,1001)), label = "Theory")
	# plt.xlabel('Number of Measurements')
	# plt.ylabel('Normalized coherence time')
	# # ax.set_yscale('log')
	# # ax.set_xscale('log')
	# plt.xlim(1,9)
	# plt.ylim(0,8)
	# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	# plt.show()
	# plt.close('all')


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


	### combine the decay of all logical two qubit states for the 3 qubit experiments.

	Q3_data = {'msmts': [], 'results': [], 'results_u':[]}

	Q3_data['results'] = np.array(data3Q['results'][0])/3.
	Q3_data['results_u'] = np.array(data3Q['results_u'][0])**2/9

	for i in range(2):
		Q3_data['results'] = Q3_data['results'] + np.array(data3Q['results'][i+1])/3.

		Q3_data['results_u'] = np.array(data3Q['results_u'][i+1])**2/9

	Q3_data['results_u'] = np.sqrt(Q3_data['results_u'])
	Q3_data['msmts'] = data3Q['msmts'][0]


	# print Q2_data
	# print Q3_data
	# print data2Q
	combined_data = [data1Q]
	combined_data.extend([Q2_data])
	combined_data.append(data2Q[-2]) ### XX scaling for 2 qubit experiments
	combined_data.append(data2Q[-1]) ### YY scaling for 2 qubit experiments
	combined_data.extend([Q3_data])
	combined_data.append({'msmts':data3Q['msmts'][3], 'results':data3Q['results'][3],'results_u':data3Q['results_u'][3]})

	# print combined_data


	labels = [r'$\langle X \rangle$',r'$\langle IX \rangle$ etc.',r'$\langle XX \rangle$',r'$\langle YY \rangle$',r'$\langle XIX \rangle$ etc.',r'$\langle XXX \rangle$']
	for ii,data in enumerate(combined_data):
		# print ii
		plt.errorbar(data[dict_keys[0]],data[dict_keys[1]],data[dict_keys[2]],fmt='o',label =labels[ii],color = color_list[ii])

	# def f(N): return (6.*N)**0.44
	
	# f = np.vectorize(f)

	# print f(0.)

	# plt.plot(np.linspace(0,15,1001),f(np.linspace(0,15,1001)), label = "Theory")
	plt.xlabel('Number of Measurements')
	plt.ylabel('Normalized coherence time')
	# ax.set_yscale('log')
	# ax.set_xscale('log')
	plt.xlim(1,9)
	plt.ylim(2,8)
	plt.legend(loc=4,prop={'size':10}) #bbox_to_anchor=(1.05, 1),
	plt.show()
	plt.close('all')





# def ShowResults():
# 	### this is a static function (no input parameters)
# 	### combines two sets (see older_timestamp_lists) of zeno measurements and plots them in a graph.
# 	### extendible for more than two measurement sets given that the evolution times stay the same.
# 	### saves a pickle file which contains the measured values as a dictionary.
# 	### the specification of 'None' for SSRO timestamps slets the analyi
# 	evo0,fid0,fid_u0,folder = Zen_Compare_Runs(msmt='0',older_timestamp_list=['20150404_175544','20150403_230000'],eRO_list=['positive','negative'],
# 	  						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)
# 	evo1,fid1,fid_u1,folder = Zen_Compare_Runs(msmt='1',older_timestamp_list=['20150404_175544','20150330_161704'],eRO_list=['positive','negative'],
# 	 						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)
# 	evo2,fid2,fid_u2,folder = Zen_Compare_Runs(msmt='2',older_timestamp_list=['20150404_175544','20150330_161704'],eRO_list=['positive','negative'],
# 	  						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)
# 	evo3,fid3,fid_u3,folder = Zen_Compare_Runs(msmt='3',older_timestamp_list=['20150404_175544','20150330_233000'],eRO_list=['positive','negative'],
# 	 						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)
# 	evo4,fid4,fid_u4,folder = Zen_Compare_Runs(msmt='4',older_timestamp_list=['20150404_175544','20150330_233000'],eRO_list=['positive','negative'],
# 	 						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)

# 	evo5,fid5,fid_u5,folder = Zen_Compare_Runs(msmt='5',older_timestamp_list=['20150414_090000','20150414_023527'],eRO_list=['positive','negative'],
# 	 						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)
# 	evo6,fid6,fid_u6,folder = Zen_Compare_Runs(msmt='6',older_timestamp_list=['20150414_090000','20150414_023527'],eRO_list=['positive','negative'],
# 	 						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)

# 	######################################
# 	### get the single qubit results. Beware, the function Zeno_proc_fidelity returns average state fidelities.

# 	evo_single, fid_single_avg, fid_u_single_avg, tstamp, folder = Zeno_proc_fidelity(older_than_tstamp='20150405_043659',
# 																		msmts='0',eRO_list=['positive','negative'],
# 																		ssro_timestamp=None,decoded_bit='2', plot_results = False, 
# 																		single_qubit = True)

# 	fid_single_proc = (3*fid_single_avg-1)/2.
# 	fid_u_single_proc = 1.5*fid_u_single_avg

# 	evo_single = evo_single*1000 ### rescale from seconds to ms.

# 	###########################


# 	fid_arr = [fid0,fid1,fid2,fid3,fid4,fid5,fid6]
# 	fid_u_arr = [fid_u0,fid_u1,fid_u2,fid_u3,fid_u4,fid_u5,fid_u6]
# 	evo_arr = [evo0,evo1,evo2,evo3,evo4,evo5,evo6]

# 	### This loop fuses the measurement runs
# 	for i in range(len(fid_arr)):
# 		f_new = []
# 		f_u_new = []
# 		for j,fids in enumerate(fid_arr[i]):

# 			if j == 0:
# 				f_new = [fid/len(fid_arr[i]) for kk,fid in enumerate(fids)]
# 				f_u_new = [fid_u**2 for kk,fid_u in enumerate(fid_u_arr[i][j])]
# 			else:
# 				### add up fidelities
# 				f_new = [f_new[kk]+fid/len(fid_arr[i]) for kk,fid in enumerate(fids)]
# 				### start constructing the error bars. It is now a sum of squared errors.
# 				f_u_new = [f_u_new[kk]+fid_u**2 for kk,fid_u in enumerate(fid_u_arr[i][j])]
# 		### Finally, take the square root for the error entries and devide by the number of traces taken into account
# 		f_u_new = [np.sqrt(f_u_new[kk])/len(fid_arr[i]) for kk in range(len(f_u_new))]

# 		### overwrite the entries in fid_arr and fid_u_arr in order to make them ready for plotting.
# 		fid_arr[i] = f_new
# 		fid_u_arr[i] = f_u_new

# 	### start to plot the results
# 	fig=plt.figure()
# 	ax=plt.subplot()


# 	for ii in range(len(fid_arr)):

# 		plt.errorbar(evo_arr[ii],fid_arr[ii],fid_u_arr[ii],marker='o',markersize=4,label=str(ii))
	
# 	plt.xlabel('Evolution time (ms)')
# 	plt.ylabel('Process fidelity')
# 	plt.title('Process fidelity averaged over ' + str(2) + ' Zeno sets')
# 	plt.legend()

# 	print 'Plots are saved in:'
# 	print folder
# 	plt.savefig(os.path.join(folder,'ZenoRuns_averaged'+'.pdf'),format='pdf')
# 	plt.savefig(os.path.join(folder,'ZenoRuns_averaged'+'.png'),format='png')
# 	plt.show()
# 	plt.close('all')

# 	#### construct a pickle container out of the data
# 	### we start off with a dictionary
# 	pickle_dict = {'0':[evo_arr[0],fid_arr[0],fid_u_arr[0]],
# 					'1':[evo_arr[1],fid_arr[1],fid_u_arr[1]],
# 					'2':[evo_arr[2],fid_arr[2],fid_u_arr[2]],
# 					'3':[evo_arr[3],fid_arr[3],fid_u_arr[3]],
# 					'4':[evo_arr[4],fid_arr[4],fid_u_arr[4]],
# 					'5':[evo_arr[5],fid_arr[5],fid_u_arr[5]],
# 					'6':[evo_arr[6],fid_arr[6],fid_u_arr[6]],
# 					'single':[evo_single,fid_single_proc,fid_u_single_proc]}
	
# 	### dump to file.			
# 	fileOut = open("Zeno_data.p","wb")
# 	pickle.dump(pickle_dict,fileOut)
# 	fileOut.close
