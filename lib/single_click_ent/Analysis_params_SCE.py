"""
Contains analysis parameters for purification measurements. Such as temporal filters
"""
import numpy as np


SPCorr_settings = {
	'st_start'		: 1865e3,
	'st_len'		: 30e3,
	'photon_channel'		: 1, ### 2 for both channels. 0 or 1 one corresponds to one HH channel.
	'ch1_offset'	: 34.95e3 
}

"""
passes adwin_filter_params in the form {'setup_key (e.g. lt3)' : {'param (e.g. CR_after)' : [enabled, min, max]}}
"""

SPSP_fltr_adwin_settings = {}

SPSP_fltr_adwin_settings['fltr_dict_lt3'] = {
				'CR_before'  : [1,2,2000],
				'CR_after'  : [1,10,2000],
			}

SPSP_fltr_adwin_settings['fltr_dict_lt4'] = 	{
				'CR_before'  : [1,2,2000],
				'CR_after'  : [1,10,2000],
				'elapsed_since_phase_stab' : [0,0e3,200e3],
				'repetition_number'		   : [0,0,2000],
				'pst_msmt_phase'		   : [0,85,95]
			}
