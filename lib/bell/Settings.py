import numpy as np

# name of the group (top-level) in the data that contains all analysis results
ANALYSISGRP = 'analysis'

# Color cycle for plotting
# same colors as in matplotlibrc -- handy for barplots, since they don't use the rc colorcycle
COLORS = 'RoyalBlue', 'Crimson', 'LimeGreen', 'DarkOrange', 'BlueViolet'

## setting the photon windows
first_win_min = 0 # The minimum of the first window
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

dif_ch0_ch_1 = 10