import numpy as np
import h5py, os
from analysis.lib.pq import pq_tools
from analysis.lib.tools import toolbox as tb
import Settings, files, Filter
import matplotlib.pyplot as plt

def _aggregated_coincidences(Base_Folder):
    in_coincidences = np.empty((0,4))
    in_coincidences = np.vstack((in_coincidences, pq_tools.get_coincidences_from_folder(Base_Folder)))
    return in_coincidences

def filter_coincidences(coincidences, Filter_settings):
    ch0_start = Filter_settings[0]
    ch1_start = Filter_settings[1]
    tail_length = Filter_settings[2]
    pulse_sep = Filter_settings[3]
    noof_pulses = Filter_settings[4]
    
    f_st0 = pq_tools.filter_synctimes(coincidences[:,1], ch0_start, ch0_start + tail_length, noof_pulses, pulse_sep, pq_file = False)
    f_st1 = pq_tools.filter_synctimes(coincidences[:,2], ch1_start, ch1_start +  tail_length, noof_pulses, pulse_sep, pq_file = False)
        
    return f_st0 & f_st1 # & f_dt

def TPQI_analysis(Base_Folder, Filter_settings, Verbose = True):
    # Gets coincident photons from Hydraharp data
    coincidences = _aggregated_coincidences(Base_Folder)


    # Defines the difference in arrival time between the coincident photons
    dts = coincidences[:,0] + abs(Filter_settings[1]-Filter_settings[0])
    dts = dts * 10**(-3)

    if Verbose:
        print
        print 'Found {} coincident photons in all runs.'.format(int(len(coincidences)))
        print '===================================='
        print

    #Filters the coincident photons by selecting only the photons emitted by the NV center
    is_sync_time_filter = filter_coincidences(coincidences, Filter_settings)
    filtered_dts = dts[is_sync_time_filter]

    if Verbose:
        print
        print 'Found {} coincident photons after filtering.'.format(int(sum(is_sync_time_filter)))
        print '===================================='
        print
        
    return dts, filtered_dts