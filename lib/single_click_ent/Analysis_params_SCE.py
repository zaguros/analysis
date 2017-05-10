"""
Contains analysis parameters for purification measurements. Such as temporal filters
"""
import numpy as np


SPCorr_settings = {
	'st_start'		: 1866e3,
	'st_len'		: 30e3,
	'photon_channel'		: 0, ### 2 for both channels. 0 or 1 one corresponds to one HH channel.
	'ch1_offset'	: 35e3 
}

SPSP_fltr_adwin_settings = {
	'list_of_adwin_params' : ['CR_after','CR_before','elapsed_since_phase_stab']

	}

SPSP_fltr_adwin_settings['fltr_dict_lt3'] = {
				'CR_before'  : [2,2000],
				'CR_after'  : [20,2000],
				'elapsed_since_phase_stab' : [0e3,1000e3]
			}

SPSP_fltr_adwin_settings['fltr_dict_lt4'] = 	{
				'CR_before'  : [2,2000],
				'CR_after'  : [20,2000],
				'elapsed_since_phase_stab' : [0e3,1000e3]
			}
