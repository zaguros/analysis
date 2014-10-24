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

Sync_num = 0
Sync_time_ph_win_1 = 1
Sync_time_ph_win_2 = 2
Chan_ph_win_1 = 3
Chan_ph_win_2 = 4
Attempts = 5
Noof_ph_LT1 = 6
Noof_ph_LT2 = 7
CR_check_bef_LT1 = 8
CR_check_aft_LT1 = 9
CR_check_bef_LT3 = 10
CR_check_aft_LT3 = 11
Psiminus_event = 12
Abs_Time = 13

