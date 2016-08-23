"""
Contains analysis parameters for purification measurements. Such as temporal filters
"""

filter_settings = {
	
	'st_start' 		: 2772.5e3,#4.5e3,
	'st_len'		: 40e3,
	'st_len_w2' 	: 40e3,
	'min_cr_lt3_before'	: 1,
	'min_cr_lt4_before'	: 1,
	'min_cr_lt3_after'	: 1,
	'min_cr_lt4_after'	: 1,
	'max_reps_w1'	: 1000,
	'min_reps_w1'	: 1,
	'max_reps_w2'	: 50,
	'min_reps_w2'	: 1,
	'max_dt'		: 40e3,
}
