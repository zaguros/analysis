"""
Contains analysis parameters for purification measurements. Such as temporal filters
"""
import numpy as np


SPCorr_settings = {
	'st_start'		: 2826e3,
	'st_len'		: 26e3,
	'photon_channel'		: 2, ### 2 for both channels. 0 or 1 one corresponds to one HH channel.
	'ch1_offset'	: 17.5e3 
}

"""
passes adwin_filter_params in the form {'setup_key (e.g. lt3)' : {'param (e.g. CR_after)' : [enabled, min, max]}}
"""

SPSP_fltr_adwin_settings = {}

SPSP_fltr_adwin_settings['fltr_dict_lt3'] = {
				'CR_before'  : [0,2,2000],
				'CR_after'  : [0,10,2000],
			}

SPSP_fltr_adwin_settings['fltr_dict_lt4'] = 	{
				'CR_before'  : [0,2,2000],
				'CR_after'  : [0,10,2000],
				'elapsed_since_phase_stab' : [0,0e3,50e3],
				'repetition_number'		   : [0,00,1000],
				'pst_msmt_phase'		   : [0,80,100],
				'DD_repetitions' 		   : [0,10,400]		
			}



data_settings = {

	'base_folder_lt3' : r'D:\measuring\data\Single_click_expm\LT3_data\EntOnDemand20170705',
	'base_folder_lt4' : r'D:\measuring\data\Single_click_expm\LT4_data\EntOnDemand20170705',
	
	'filenames_for_expms' : {'SweepPhase0p05' : r'EntangleXsweepY',
							 'SweepTheta'		: r'EntangleSweepEverything'},


}