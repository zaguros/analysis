"""
Contains analysis parameters for purification measurements. Such as temporal filters
"""

temporal_filter = {
	
	'st_start' 		: 2465e3,
	'st_len'		: 50e3,
	'st_len_w2' 	: 50e3,
}

decoherence_filter = {
	'max_reps_w1'	: 250,
	'min_reps_w1'	: 1,
	'max_reps_w2'	: 250,
	'min_reps_w2'	: 1,
}
