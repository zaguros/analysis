"""
This script analyses the Z initializations of the specific carbons annd prints the most recent values.
"""

import numpy as np
import os,re
import h5py
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import ssro
from analysis.lib.tools import plot

from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import matplotlib
import copy as cp

reload(common)
reload(toolbox)


#### this block reads tomography files and gets <Z>, <ZZ> and <ZZZ> expectation values.

ssro_calib_timestamp = '20150610_174113'
ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'


if True:
	carbons  =[1,2,5]
	RO_list=  ['negative','positive']
	search_strings = []
	for c in carbons:
		search_strings.append('CarbonRO_'+str(c))

		carbon_cross = cp.deepcopy(carbons)
		carbon_cross.remove(c)
		for cc in carbon_cross:
			cc_string = str(carbon_cross[0])+str(carbon_cross[1])
		search_strings.append('CarbonRO_'+cc_string)

	search_strings.append('CarbonRO_'+'125')

	print search_strings

	newer_than = '20150610_180000'
	older_than = '20150611_230000'
	RO = ''
	for s in search_strings:
		folder_list = []
		timestamp = older_than ## reset.
		while toolbox.latest_data(	contains = s+'_'+RO,
									return_timestamp = True,
									newer_than  =newer_than,
									older_than = timestamp,
									raise_exc = False) != False:
			for RO in RO_list:
				timestamp,folder = toolbox.latest_data(contains = s+'_'+RO,
					return_timestamp = True,
					newer_than  =newer_than,
					older_than = timestamp)
				folder_list.append(folder)

		# print folder_list
		res_list = []
		res_u_list = []

		for folder in folder_list:
			a = mbi.MBIAnalysis(folder)
			a.get_sweep_pts()
			a.get_readout_results(name='adwindata')
			a.get_electron_ROC(ssro_calib_folder)

			res_list.append(2*abs(a.p0[0]-0.5))
			res_u_list.append(abs(2*a.u_p0[0])**2)


		res = np.sum(np.array(res_list))/len(res_list)
		res_u = np.sqrt(np.sum(np.array(res_u_list)))/len(res_list)

		print s
		print 'number of measurements', len(res_list)
		print 'average measured expectation value: ', res
		print 'uncertainty: ', res_u


############################################################################################

# Second block of this script: calculate the correction values for individual read-outs.

############################################################################################

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

	# n_pop_u = 0.0000002

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

#### the correction values from the first block of this script


CarbonROC_Dict = {'values' : {} , 'error' : {}}

CarbonROC_Dict['values']['1']	=	0.859854829704
CarbonROC_Dict['values']['2']	=	0.846657530075
CarbonROC_Dict['values']['5']	=	0.88986599665
CarbonROC_Dict['values']['12']	=	0.734843046401
CarbonROC_Dict['values']['15']	=	0.784368785197
CarbonROC_Dict['values']['25']	=	0.747335160652
CarbonROC_Dict['values']['125']	=	0.666209837064
CarbonROC_Dict['values']['nitrogen'] = 0.963484593103#0.963891440771
CarbonROC_Dict['error']['1']	=	0.00621794744892
CarbonROC_Dict['error']['2']	=	0.00775792880472
CarbonROC_Dict['error']['5']	=	0.00735657065159
CarbonROC_Dict['error']['12']	=	0.00565949773662
CarbonROC_Dict['error']['15']	=	0.00217798138803
CarbonROC_Dict['error']['25']	=	0.00219162359014
CarbonROC_Dict['error']['125']	=	0.0089820755304
CarbonROC_Dict['error']['nitrogen'] = 0.0149696842867#0.0148976423135


SolveDict = calc_ro_correction(CarbonROC_Dict['values'],CarbonROC_Dict['error'])

val_dict,err_dict = SolveDict['values'],SolveDict['error']

for key in val_dict.keys():
	print 'Correction for the RO combination ' + key + ' ' +str(round(val_dict[key],5)) +' +- ' + str(round(err_dict[key],5))