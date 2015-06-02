import numpy as np
import h5py, os
from analysis.lib.pq import pq_tools
from analysis.lib.tools import toolbox as tb
import Settings, Filter
import matplotlib.pyplot as plt

def _aggregated_coincidences(Base_Folder):
    in_coincidences = np.empty((0,4))
    in_coincidences = np.vstack((in_coincidences, pq_tools.get_coincidences_from_folder(Base_Folder)))
    return in_coincidences

def filter_coincidences(coincidences, ch0_start, ch1_start, WINDOW_LENGTH, dif_win1_win2, noof_pulses,
                         column_st_0, column_st_1):

    tail_length = WINDOW_LENGTH
    pulse_sep = dif_win1_win2

    
    f_st0 = pq_tools.filter_synctimes(coincidences[:,column_st_0], ch0_start, ch0_start + tail_length, noof_pulses, pulse_sep, pq_file = False)
    f_st1 = pq_tools.filter_synctimes(coincidences[:,column_st_1], ch1_start, ch1_start +  tail_length, noof_pulses, pulse_sep, pq_file = False)
        
    return f_st0 & f_st1 # & f_dt

def TPQI_analysis(Base_Folder, ch0_start, ch1_start, WINDOW_LENGTH, dif_win1_win2, noof_pulses, Verbose = True):
    # Gets coincident photons from Hydraharp data
    coincidences = _aggregated_coincidences(Base_Folder)

    dt_index = 0
    column_st_0 = 1
    column_st_1 = 2
    column_sync_num_ch0 = 3

    # Defines the difference in arrival time between the coincident photons
    dts = coincidences[:,dt_index] + (ch1_start - ch0_start)
    dts = dts * 10**(-3)

    if Verbose:
        print
        print 'Found {} coincident photons in all runs.'.format(int(len(coincidences)))
        print '===================================='
        print

    #Filters the coincident photons by selecting only the photons emitted by the NV center
    is_sync_time_filter = filter_coincidences(coincidences, ch0_start, ch1_start, WINDOW_LENGTH,
                                                 dif_win1_win2, noof_pulses, column_st_0, column_st_1)
    filtered_dts = dts[is_sync_time_filter]

    if Verbose:
        print
        print 'Found {} coincident photons after filtering.'.format(int(sum(is_sync_time_filter)))
        print '===================================='
        print
        
    return dts, filtered_dts