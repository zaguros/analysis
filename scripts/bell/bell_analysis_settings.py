import datetime

bs_process_params = {

    'st_start_ch0'            : 5430000,
    'st_len'                  : 250000 ,
    'st_start_ch1'            : 5430000,
    'pulse_sep'               : 250000 ,  #XXXX

    'ent_marker_channel_bs' : 1,

    'st_pulse_start'          : 5438000,
    'st_pulse_len'            : 4000,
    'pulse_max_sn_diff'       : 3000000, #3 million syncs ~ 60 secs

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
    'psb_prepulse_start_lt3'        : 7470,
    'psb_prepulse_start_lt4'        : 5366,
    'psb_prepulse_len'              : 20,


}

analysis_params = {

    'st_start_ch0'             : 5444300,
    'st_len'                   : 55000 ,
    'st_start_ch1'             : 5444300+1000,
    'st_len_w2_00'             : 1000,
    'st_len_w2_11'             : 15000,
    'pulse_sep'                : bs_process_params['pulse_sep'],
    'st_pulse_start'           : bs_process_params['st_pulse_start'],
    'st_pulse_len'             : bs_process_params['st_pulse_len'],
    'hist_binsize_ps'		  : 100, #ps 

    'F0A' : 0.950,
    'F1A' : 0.980,
    'F0B' : 0.935, 
    'F1B' : 0.995,
}
    #XX-high_str F0A=0.9320, F1A=0.9966, F0B=0.915, F1B=0.9976):
    #XX_lock F0A=0.9380, F1A=0.9908, F0B=0.9096, F1B=0.9976):
    #ZZ F0A=0.9550, F1A=0.9786, F0B=0.9376, F1B=0.9952):
    #XXlotr F0A=0.96, F1A=0.99, F0B=0.943, F1B=0.995):
#Beam splitter light on on webcam or wavemeter not working (night of 19/05/2015)
bad_time_ranges = [[datetime.datetime(2015,3, 4,10, 0,33),      datetime.datetime(2015,3, 4,10,14, 4)], 
                   [datetime.datetime(2015,3, 9,15,50,49),      datetime.datetime(2015,3, 9,16, 4, 3)],
                   [datetime.datetime(2015,3,10,11,13, 6),      datetime.datetime(2015,3,10,11,19,22)],
                   [datetime.datetime(2015,3,11, 8,13,31),      datetime.datetime(2015,3,11, 8,18,57)],
                   [datetime.datetime(2015,3,11,17,32,42),      datetime.datetime(2015,3,11,17,56, 4)],
                   [datetime.datetime(2015,3,12,12,55,24),      datetime.datetime(2015,3,12,12,57,22)],
                   [datetime.datetime(2015,3,19,12,31,54),      datetime.datetime(2015,3,19,12,43,26)],
                   [datetime.datetime(2015,3,25, 5,51, 0),      datetime.datetime(2015,3,25, 7, 0, 0)],
                   [datetime.datetime(2015,3,26, 7,48,54),      datetime.datetime(2015,3,26, 7,49,13)],
                   [datetime.datetime(2015,3,26,14,56,58),      datetime.datetime(2015,3,26,15, 8,43)],
                   [datetime.datetime(2015,4,23,14,31, 8),      datetime.datetime(2015,4,23,14,36,22)],
                   [datetime.datetime(2015,4,24,11,14,54),      datetime.datetime(2015,4,24,11,23,31)],
                   [datetime.datetime(2015,4,28, 1,30, 0),      datetime.datetime(2015,4,28,20, 0, 0)], #ch1-ch0 offset 4th
                   [datetime.datetime(2015,4,29, 6,40, 0),      datetime.datetime(2015,4,29, 7,40, 0)], #ch1-ch0 offset 4th
                   [datetime.datetime(2015,4,29,16,13, 0),      datetime.datetime(2015,4,29,23,30, 0)], #ch1-ch0 offset 4th
                   [datetime.datetime(2015,5,15, 8,50,00),      datetime.datetime(2015,5,15,10, 0, 0)], #ch1-ch0 offset ZZ
                   [datetime.datetime(2015,5,19,23,22,13),      datetime.datetime(2015,5,20,10,25,10)],   #lasermeister crash ZZ
                   [datetime.datetime(2015,6,12,21,11,00),      datetime.datetime(2015,6,13,04,25,10)], #ch1-ch0 offset XX
                   [datetime.datetime(2015,6,13,15,31,44),      datetime.datetime(2015,6,13,18,22,23)], #Bell_optimizer not running on LT3 XX
                   [datetime.datetime(2015,6,13,18,22,23),      datetime.datetime(2015,6,13,20, 8,23)], # Bell_optimizer not ruuning on LT3 XX
                   [datetime.datetime(2015,6,14, 4,34,07),      datetime.datetime(2015,6,14, 9, 0, 0)], #trouble with yellow laser or WM from LT4 XX
                   ]