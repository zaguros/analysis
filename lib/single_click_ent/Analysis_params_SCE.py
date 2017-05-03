"""
Contains analysis parameters for purification measurements. Such as temporal filters
"""
import numpy as np


SPCorr_settings = {
	'st_start'		: 1864e3,
	'st_len'		: 30e3,
	'photon_channel'		: 0, ### 2 for both channels. 0 or 1 one corresponds to one HH channel.
	'ch1_offset'	: 35e3 
}

SPSP_fltr_adwin_settings = {
	'list_of_adwin_params' : ['CR_after'],
	'fltr_minima_dict_lt3'  : 
			{
				'CR_after'  : 3,
			}

	'fltr_maxima_dict_lt3'  : 			
			{
				'CR_after'  : 2000,
			}

	'fltr_minima_dict_lt4'  : 
			{
				'CR_after'  : 3,
			}

	'fltr_maxima_dict_lt4'  : 			
			{
				'CR_after'  : 2000,
			}

}