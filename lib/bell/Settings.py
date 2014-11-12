import numpy as np

# name of the group (top-level) in the data that contains all analysis results
ANALYSISGRP = 'analysis'

# Color cycle for plotting
# same colors as in matplotlibrc -- handy for barplots, since they don't use the rc colorcycle
COLORS = 'RoyalBlue', 'Crimson', 'LimeGreen', 'DarkOrange', 'BlueViolet'

## setting the photon windows
first_win_min = 0 # The minimum of the first window
first_win_min_ch0 = 0 # The minimum of the first window
first_win_max = 10 # The maximum of the first window, also the minimum of the second window
second_win_min = 20
second_win_max = 30 # The maximum of the second window

## Sets if the selection of the entanglement events should be re-executed.
Force_eval = True

## Sets print commands
VERBOSE = True

## Set filters
WINDOW_LENGTH = 100.
ch0_start = 5378.
ch1_start = 5378.
ch0_stop = ch0_start + WINDOW_LENGTH
ch1_stop = ch1_start + WINDOW_LENGTH
dif_win1_win2 = 600

CNTR_PEAK_RANGE = 0

dif_ch0_ch1 = 10

"""
Returns all entanglement events. 
Colums are:
Sync Number | Sync Time Photon 1 | Sync Time photon 2 | Photon 1 Channel | Photon 2 Channel | Attempts
| Amount of Photons LT1 | Amount of Photons LT 3 | CR Check Before LT1 | CR Check After LT1 | CR Check Before LT3 
| CR Check After LT3 | psiminus | absolute time
"""

Sync_num_BS = 0
Sync_time_ph_win_1_BS = 1
Sync_time_ph_win_2_BS = 2
Chan_ph_win_1_BS = 3
Chan_ph_win_2_BS = 4
Psiminus_event = 5
Abs_Time = 6


Sync_num_LT3 = 7
Num_phot_LT3 = 8
RND_num_ind_LT3 = 9
RND_num_LT3 = 10
Sync_Time_RND_num_LT3 = 11
Sync_Time_photon_1_LT3 = 12
Sync_Time_photon_2_LT3 = 13
Sync_Time_photon_3_LT3 = 14
Sync_Time_photon_4_LT3 = 15
Sync_Time_photon_5_LT3 = 16
Sync_Time_photon_6_LT3 = 17
Sync_Time_photon_7_LT3 = 18
Sync_Time_photon_8_LT3 = 19
Sync_Time_photon_9_LT3 = 20
Sync_Time_photon_10_LT3 = 21
Sync_Time_photon_11_LT3 = 22
Sync_Time_photon_12_LT3 = 23
Sync_Time_photon_13_LT3 = 24
Sync_Time_photon_14_LT3 = 25
Sync_Time_photon_15_LT3 = 26
Sync_Time_photon_16_LT3 = 27
Sync_Time_photon_17_LT3 = 28
Sync_Time_photon_18_LT3 = 29
Sync_Time_photon_19_LT3 = 30
Sync_Time_photon_20_LT3 = 31
Sync_Time_photon_21_LT3 = 32
Sync_Time_photon_22_LT3 = 33
Sync_Time_photon_23_LT3 = 34
Sync_Time_photon_24_LT3 = 35

Sync_num_LT4 = 36
Num_phot_LT4 = 37
RND_num_ind_LT4 = 38
RND_num_LT4 = 39
Sync_Time_RND_num_LT4 = 40
Sync_Time_photon_1_LT4 = 41
Sync_Time_photon_2_LT4 = 42
Sync_Time_photon_3_LT4 = 43
Sync_Time_photon_4_LT4 = 44
Sync_Time_photon_5_LT4 = 45
Sync_Time_photon_6_LT4 = 46
Sync_Time_photon_7_LT4 = 47
Sync_Time_photon_8_LT4 = 48
Sync_Time_photon_9_LT4 = 49
Sync_Time_photon_10_LT4 = 50
Sync_Time_photon_11_LT4 = 51
Sync_Time_photon_12_LT4 = 52
Sync_Time_photon_13_LT4 = 53
Sync_Time_photon_14_LT4 = 54
Sync_Time_photon_15_LT4 = 55
Sync_Time_photon_16_LT4 = 56
Sync_Time_photon_17_LT4 = 57
Sync_Time_photon_18_LT4 = 58
Sync_Time_photon_19_LT4 = 59
Sync_Time_photon_20_LT4 = 60
Sync_Time_photon_21_LT4 = 61
Sync_Time_photon_22_LT4 = 62
Sync_Time_photon_23_LT4 = 63
Sync_Time_photon_24_LT4 = 64