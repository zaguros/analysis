import datetime

bs_process_params = {

    'st_start_ch0'            : 5420000,
    'st_len'                  : 250000 ,
    'st_start_ch1'            : 5420000,
    'pulse_sep'               : 250000 ,

    'ent_marker_channel_bs' : 1,

    'st_pulse_start_ch0'      : 5443000,
    'st_pulse_start_ch1'      : 5444000, #for ch 1
    'st_pulse_len'            : 5000,
    'pulse_max_sn_diff'       : 6000000, #6 million syncs ~ 120 secs

    'hist_binsize_ps'		  : 100, #ps 
    
    }
after_pulse_process_params = {
    
    'first_st_start'         : bs_process_params['st_start_ch0'] + bs_process_params['pulse_sep'],
    'first_st_len'           : 50000,
    'after_pulse_st_start'   : bs_process_params['st_start_ch0'] + 2*bs_process_params['pulse_sep'],
    'after_pulse_st_len'     : 50000,
    
    }

lt_process_params = {

    'ent_marker_channel_lt3'         : 4,
    'ent_marker_channel_lt4'         : 4,
    'ent_marker_lt_timebin_limit'    : 10000,
    'sn_diff_marker_ent_early'       : -1,
    'sn_diff_marker_ent_late'        : 0,
    'invalid_marker_channel_lt'      : 8,
    'invalid_marker_max_sn_diff'     : 251,
    
    'ro_channel'                    : 0,
    'ro_start'                      : 10620,
    'ro_length'                     : 3700,
    
    'rnd_channel'                   : 1,
    'rnd_start'                     : 10000,
    'rnd_length'                    : 1000,
    'rnd_0_channel'                 : 1,
    'rnd_1_channel'                 : 2,
    
    'psb_tail_start_lt3'            : 7480,
    'psb_tail_start_lt4'            : 5350,
    'psb_tail_len'                  : 200,
    'pulse_sep'              	    : bs_process_params['pulse_sep']/1000,
    'psb_prepulse_start_lt3'        : 7570,
    'psb_prepulse_start_lt4'        : 5366,
    'psb_prepulse_len'              : 20,


}

analysis_params = {

    'st_start_ch0'             : 5425300,
    'st_len'                   : 55000 ,
    'st_start_ch1'             : 5426000-900,
    'st_len_w2_00'             : 4000,
    'st_len_w2_11'             : 2500,
    'pulse_sep'                : bs_process_params['pulse_sep'],
    'st_pulse_start'           : bs_process_params['st_pulse_start_ch0'],
    'st_pulse_len'             : bs_process_params['st_pulse_len'],
    'hist_binsize_ps'		   : 100, #ps 
    'ro_start'                 : 10620,
    'ro_length'                : 3700,
    'invalid_marker_max_sn_diff': 251,

    'F0A' : 0.9550,
    'F1A' : 0.9786,
    'F0B' : 0.9376, 
    'F1B' : 0.9952,
}
    #XX-high_str F0A=0.9320, F1A=0.9966, F0B=0.915, F1B=0.9976):
    #XX_lock F0A=0.9380, F1A=0.9908, F0B=0.9096, F1B=0.9976):
    #ZZ F0A=0.9550, F1A=0.9786, F0B=0.9376, F1B=0.9952):
    #XXlotr F0A=0.96, F1A=0.99, F0B=0.943, F1B=0.995):
