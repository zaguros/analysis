"""
Contains analysis parameters for purification measurements. Such as temporal filters
"""
import numpy as np


SPCorr_settings = {
	'st_start'		: 1843e3,
	'st_len'		: 60e3,
	'photon_channel'		: 2, ### 2 for both channels. 0 or 1 one corresponds to one HH channel.
	'ch1_offset'	: 35e3 ### not in use right now
}