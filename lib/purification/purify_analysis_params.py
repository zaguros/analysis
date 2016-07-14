"""
Contains analysis parameters for purification measurements. Such as temporal filters
"""

temporal_filter = {
	
	'st_start' 		: 2752e3,
	'st_len'		: 20e3,
	'st_len_w2' 	: 20e3,
}

decoherence_filter = {
	'max_reps_w1'	: 250,
	'min_reps_w1'	: 1,
	'max_reps_w2'	: 250,
	'min_reps_w2'	: 1,
}
