import numpy as np
import h5py
import time
import os
import gc as gc
import Analysis
from datetime import datetime
from analysis.lib.tools import toolbox as tb
from analysis.lib.pq import pq_tools, pq_plots


def get_Bell_events(fp_BS,fp_LT3,fp_LT4, BS_marker_chan, first_win_min, 
        first_win_max, second_win_min, second_win_max, force_eval=False, VERBOSE = True):
    
    """
    Returns either the entanglement events already saved; the corresponding attributes
    and save = False, or returns the newly calculated entanglement events;
    the corresponding attributes (column names) and the save = True. Put in the file path 
    of the BS then of LT3 then of LT4 then the marker channel and then the repitition if
    it's looped. Also put in the windows to determine if an photon is the first photon to arrive.
    (first_win_min,first_win_max,second_win_min, second_win_max)
    """
    folder, name =  os.path.split(tb.get_msmt_header(fp_BS))
    if VERBOSE:
        print name
 
    if tb.has_analysis_data(fp_LT3, 'Total_SSRO_events'):
        Total_SSRO_events_LT3, _att_LT3 = tb.get_analysis_data(fp_LT3, 'Total_SSRO_events')

    if tb.has_analysis_data(fp_LT4, 'Total_SSRO_events'):
        Total_SSRO_events_LT4, _att_LT4 = tb.get_analysis_data(fp_LT4, 'Total_SSRO_events')

    if tb.has_analysis_data(fp_BS, 'Entanglement_events') and not force_eval:
        entanglement_events, _a = tb.get_analysis_data(fp_BS, 'Entanglement_events')

    if force_eval or not tb.has_analysis_data(fp_BS, 'Entanglement_events'):
        # Opens beamsplitter data 
        f = h5py.File(fp_BS, 'r')
        sync_times = f['/PQ_sync_time-1'].value 
        sync_numbers = f['/PQ_sync_number-1'].value
        channel = f['/PQ_channel-1'].value
        abs_times = f['/PQ_time-1'].value 
        f.close()

        sync_num_with_markers = sync_numbers[pq_tools.filter_marker(fp_BS,BS_marker_chan, VERBOSE = False)]
        unique_sync_num_with_markers = np.unique(sync_num_with_markers)

        if (len(Total_SSRO_events_LT3[:,2]) == len(unique_sync_num_with_markers)) & \
                                        (len(Total_SSRO_events_LT4[:,2]) == len(unique_sync_num_with_markers)):
            print 
            print 
            print "The number of markers matches and is:", len(unique_sync_num_with_markers)
            print "======================================================================="
            print      
        elif '215201_Bell_BS_full_BellLFBT_day2_Run8' in fp_BS:
            sync_times_with_marker = sync_times[pq_tools.filter_marker(fp_BS,BS_marker_chan, VERBOSE = VERBOSE)]
            last_sync_time = sync_times_with_marker[len(sync_times_with_marker)-1]
            print 
            print 
            print "Filepath:", fp_BS
            print "The number of markers does not match!!!!"
            print "The number of BS markers is:", len(unique_sync_num_with_markers)
            print "We have lookend into this in the data and there is a corrupt marker in the BS data which has a sync time of", last_sync_time
            print "The number of LT3 markers is:", len(Total_SSRO_events_LT3[:,2])
            print "The number of LT4 markers is:", len(Total_SSRO_events_LT4[:,2])
            print "======================================================================="
            print
            print
        else:
            print 
            print 
            print "Filepath:", fp_BS
            print "The number of markers does not match!!!!"
            print "The number of BS markers is:", len(unique_sync_num_with_markers)
            print "The number of LT3 markers is:", len(Total_SSRO_events_LT3[:,2])
            print "The number of LT4 markers is:", len(Total_SSRO_events_LT4[:,2])
            print "======================================================================="
            print
            print
            #raise

        # Gets filters for photons with markers in the first and second window
        # from the Filter file
        is_photon_1st_window_with_markers, is_photon_2nd_window_with_markers =\
                                        pq_tools.get_photons_with_markers(fp_BS, BS_marker_chan,
                                            first_win_min, first_win_max, second_win_min, second_win_max, VERBOSE = VERBOSE)

        # Gets filters for photons with markers in the first and second window
        # from the Filter file
        is_photon_1st_window_with_markers, is_photon_2nd_window_with_markers =\
                                        pq_tools.get_photons_with_markers(fp_BS, BS_marker_chan,
                                            first_win_min, first_win_max, second_win_min, second_win_max)

        # Retrieves sync numbers and sync times for photons both in the first
        # and 2nd window
        Sync_num_1st_window_with_markers = sync_numbers[is_photon_1st_window_with_markers]
        channel_1st_window_with_markers = channel[is_photon_1st_window_with_markers]
        Sync_times_1st_window_with_markers = sync_times[is_photon_1st_window_with_markers]

        Sync_num_2nd_window_with_markers = sync_numbers[is_photon_2nd_window_with_markers]
        channel_2nd_window_with_markers =  channel[is_photon_2nd_window_with_markers]
        Sync_times_2nd_window_with_markers = sync_times[is_photon_2nd_window_with_markers]
       
       
        # Defines a filter for all events with markers
        is_all_markers = is_photon_1st_window_with_markers | is_photon_2nd_window_with_markers
        
        # Gets the absolute times for all events with makers
        PLU_mrkr_abs_times = abs_times[is_all_markers]
        
        #Initializes the final array of entanglement events
        entanglement_events = np.empty((0,7), np.uint64)

        columns = "Sync_Number, Sync_Time_photon_1, Sync_Time_photon_2, Channel_photon_1,\
            Channel_photon_2, psiminus, abs_time"

        _a = {'Columns': columns}

        # Get all real entanglement events, loops over sync numbers
        for i,s in enumerate(unique_sync_num_with_markers):
            
            # The attempt is defined as the sync number modulo 250 
            #(250 = the repitition rate)
            
            # Return filters for specific sync number s
            is_ph_1st_win_sync_num_s = Sync_num_1st_window_with_markers == s
            is_ph_2nd_win_sync_num_s = Sync_num_2nd_window_with_markers == s
            
            # Test if there is one photon in both windows
            if len(Sync_num_1st_window_with_markers[is_ph_1st_win_sync_num_s]) == 1 \
                and len(Sync_num_2nd_window_with_markers[is_ph_2nd_win_sync_num_s]) == 1:

                # Saves sync times an channels of both photons
                stimes = np.array([ Sync_times_1st_window_with_markers[is_ph_1st_win_sync_num_s],\
                    Sync_times_2nd_window_with_markers[is_ph_2nd_win_sync_num_s]]).reshape(-1)
                channel_1 = channel_1st_window_with_markers[is_ph_1st_win_sync_num_s]
                channel_2 = channel_2nd_window_with_markers[is_ph_2nd_win_sync_num_s]
                chans = np.array([channel_1,channel_2])            

                # Determines if event is psiplus or psiminus
                if channel_1 == channel_2:
                    psiminus = 0
                else:
                    psiminus = 1
                
            # Test if there are two photons in the first window
            elif len(Sync_num_1st_window_with_markers[is_ph_1st_win_sync_num_s]) == 2 and \
                        len(Sync_num_2nd_window_with_markers[is_ph_2nd_win_sync_num_s]) == 0:
                
                # Saves sync times an channels of both photons
                stimes = Sync_times_1st_window_with_markers[is_ph_1st_win_sync_num_s]
                chans = Channel_1st_window_with_markers[is_ph_1st_win_sync_num_s]
                
                # Set psiminus to two meaning that there is no entanglement since both photons
                # are in first window
                psiminus  = 2 
         
            # Test if there are two photons in the second window   
            elif len(Sync_num_1st_window_with_markers[is_ph_1st_win_sync_num_s]) == 0 and \
                        len(Sync_num_2nd_window_with_markers[is_ph_2nd_win_sync_num_s]) == 2:
                
                # Saves sync times an channels of both photons
                stimes = Sync_times_2nd_window_with_markers[is_ph_2nd_win_sync_num_s]
                chans = Channel_2nd_window_with_markers[is_ph_2nd_win_sync_num_s]

                # Set psiminus to two meaning that there is no entanglement since both photons
                # are in second window
                psiminus  = 3 

            # Disregards events with more than two photons
            else:
                    continue

                
            """
            Returns all entanglement events. 
            Colums are:
            Sync Nymber BS | Sync Time Photon 1 BS | Sync Time photon 2 BS | Photon 1 Channel BS | 
            Photon 2 Channel BS | psiminus | absolute time BS
            """
           
                  
            _event = np.array([s, 
                            stimes[0],
                            stimes[1], 
                            chans[0], 
                            chans[1], 
                            psiminus, 
                            PLU_mrkr_abs_times[i]], dtype = np.uint64)
                            
            entanglement_events = np.vstack((entanglement_events, _event))

        if tb.has_analysis_data(fp_BS,'Entanglement_events'):
            tb.clear_analysis_data(fp_BS)
            tb.set_analysis_data(fp_BS, 'Entanglement_events', entanglement_events, _a)
        else:
            tb.set_analysis_data(fp_BS, 'Entanglement_events', entanglement_events, _a)


    if VERBOSE:
        print
        print 'Found {} valid entanglement events.'.format(int(len(entanglement_events)))
        print '===================================='
        print
        if '215201_Bell_BS_full_BellLFBT_day2_Run8' in fp_BS:
            print
            print 'Found {} valid entanglement events.'.format(int(len(entanglement_events)))
            print 'This does not correspond with the number of markers because there is one corrupt marker'
            print '===================================='
            print


    BS_LT3_data = np.hstack((entanglement_events, Total_SSRO_events_LT3))
    All_combined_data = np.hstack((BS_LT3_data,Total_SSRO_events_LT4))

    Combined_attributes = _a['Columns'] + ', ' + _att_LT3['Columns'] + ', ' + _att_LT4['Columns']
    _combined_attributes = {'Columns': Combined_attributes}

    return All_combined_data, _combined_attributes

#######################  SSRO events #################################

def get_total_SSRO_events(pqf, RO_start, RO_length, marker_chan, chan_rnd_0, chan_rnd_1, sync_time_lim, VERBOSE = True):
    """
    Returns SSRO data for all marked events. 
    Colums are:
    Sync Nymber | number of photons | RND number indicator | RND number | Sync Time RND number | Sync Times photon 1-24 |
    """

    _a = get_attributes(pqf)

    # Gets the number of blocks in the data
    num_blocks = tb.get_num_blocks(pqf)
   
    if VERBOSE or (num_blocks > 1):
        print 'The total number of blocks is:', num_blocks
 
    # Initializes arrays to save the PQ-data
    PQ_sync_number = np.empty((0,), dtype = np.uint32) 
    PQ_special = np.empty((0,), dtype = np.uint32)         
    PQ_sync_time = np.empty((0,), dtype = np.uint64)
    PQ_time = np.empty((0,), dtype = np.uint64)      
    PQ_channel = np.empty((0,), dtype = np.uint32)

    # Initializes an array to save the SSRO data
    total_SSRO_events = np.empty((0,30), dtype = np.uint64)

    # Loops over every block
    for i in range(num_blocks):
        # Get the SSRO events and PQ data for these sync numbers
        gc.collect()
        _events, _PQ_sync_number, _PQ_special, _PQ_sync_time, _PQ_time, _PQ_channel = \
                get_SSRO_events(pqf, marker_chan, RO_start, RO_length, chan_rnd_0, chan_rnd_1, sync_time_lim = sync_time_lim, index = i+1, VERBOSE = VERBOSE)

        # Concatenates all PQ data
        PQ_sync_number = np.hstack((PQ_sync_number,_PQ_sync_number)) 
        PQ_special = np.hstack((PQ_special, _PQ_special))         
        PQ_sync_time = np.hstack((PQ_sync_time, _PQ_sync_time)) 
        PQ_time = np.hstack((PQ_time, _PQ_time))      
        PQ_channel = np.hstack((PQ_channel, _PQ_channel)) 
                
        # Stacks all SSRO data    
        total_SSRO_events = np.vstack((total_SSRO_events, _events))
        
        if VERBOSE:
            print
            print 'Found {} valid marked SSRO events in block'.format(int(len(_events))), i+1
            print '===================================='
            print

    print
    print 'Found {} valid marked SSRO events in all blocks'.format(int(len(total_SSRO_events)))
    print '===================================='
    print       

    return total_SSRO_events, _a, PQ_sync_number, PQ_special, PQ_sync_time, PQ_time, PQ_channel


def get_total_SSRO_events_quick(pqf, RO_start, RO_length, marker_chan, sync_time_lim, VERBOSE = True):
    """
    Returns quick SSRO data for all marked events. 
    Sync Nymber | number of photons | Sync time first photon
    """
    
    columns = "Sync_Number,Number of photons,Sync_Time_photon_1"
    _a = {'Columns': columns}

    # Gets the number of blocks for the data
    num_blocks = tb.get_num_blocks(pqf)

    if VERBOSE:
        print 'The total number of blocks is:', num_blocks

    # Initializes the array to save all SSRO events
    total_SSRO_events = np.empty((0,3))

    # Loop over all blocks
    for i in range(num_blocks):
        _events = get_SSRO_events_quick(pqf, marker_chan,  RO_start, RO_length, sync_time_lim = sync_time_lim, index = i+1)

        # Stacks all events for several blocks
        total_SSRO_events = np.vstack((total_SSRO_events, _events))

        if VERBOSE:
            print
            print 'Found {} valid marked SSRO events in block'.format(int(len(_events))), i+1
            print '===================================='
            print


    if VERBOSE:
        print
        print 'Found {} valid marked SSRO events in all blocks'.format(int(len(total_SSRO_events)))
        print '===================================='
        print       

    return total_SSRO_events, _a


def get_SSRO_events(pqf, marker_chan ,RO_start, RO_length, chan_rnd_0, chan_rnd_1, sync_time_lim = 0, index = 1, VERBOSE = True):
    """
    Returns an array with sync numbers in the first column, the number of photons in the readout window
    in the second column, a random number generation check in the third colum (1 if one is generated 0 if not),
    the random number itself in the fourth column, the sync time of the random number in the fifth column, 
    and the sync times of the 1st to the 24th photon.
    """

    # Define all block names
    sync_time_name = '/PQ_sync_time-' + str(index)
    sync_num_name = '/PQ_sync_number-' + str(index)
    spec_name = '/PQ_special-' + str(index)
    chan_name = '/PQ_channel-' + str(index)
    time_name = '/PQ_time-' + str(index)
  

    # Open files to determine if there are markers

    if type(pqf) == h5py._hl.files.File: 
        special_RO =pqf[spec_name].value
        channel_RO = pqf[chan_name].value

    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        special_RO = f[spec_name].value
        channel_RO = f[chan_name].value

        f.close()
    else:
        print "Neither filepath nor file enetered in function please check:", pqf

    # Initializes an array to save all SSRO data
    SSRO_events = np.empty((0,30), dtype = np.uint64)

    # Initializes arrays to save all PQ data
    PQ_sync_number = np.empty((0,), dtype = np.uint32) 
    PQ_special = np.empty((0,), dtype = np.uint32)         
    PQ_sync_time = np.empty((0,), dtype = np.uint64)
    PQ_time = np.empty((0,), dtype = np.uint64) 
    PQ_channel = np.empty((0,), dtype = np.uint32)
   
    is_mrkr, num_mrkr = pq_tools.get_markers_with_num_markers(pqf, marker_chan, index = index)
    
    if VERBOSE:
        print "The number of markers is:", num_mrkr

    # If there are no markers, no further analysis is necessary and empty arrays are returned. 
    if num_mrkr > 0: 
        # Get a list with the unique sync numbers
        unique_sync_num_with_markers = get_unique_sync_with_mrkr_quick(pqf, is_mrkr,sync_num_name, sync_time_name, sync_time_lim, num_mrkr)

        PQ_save_sync_numbers = np.empty((num_mrkr * 3), dtype = np.uint32)

        for i,j in enumerate(unique_sync_num_with_markers):
            PQ_save_sync_numbers[i*3] = j - 1
            PQ_save_sync_numbers[(i*3) + 1] = j
            PQ_save_sync_numbers[(i*3) + 2] = j + 1

        num_sync_num = pq_tools.get_length_block(pqf, index = index)

        indices_PQ_save_sync_num = get_indices_PQ_save_sync_num(pqf, PQ_save_sync_numbers, sync_num_name, num_sync_num)
        
        if VERBOSE:
            print "The number of unique sync numbers that have a marker is:", len(unique_sync_num_with_markers)
            print "The number of sync numbers will thus be save is:", len(PQ_save_sync_numbers)
            print "The number of events with these sync numbers is:", len(indices_PQ_save_sync_num)

        PQ_sync_number, PQ_special, PQ_sync_time, PQ_time, PQ_channel = \
                                                get_PQ_save_data(pqf,indices_PQ_save_sync_num, index = index)
        
        if type(pqf) == h5py._hl.files.File: 
            sync_num_RO_data = pqf[sync_num_name]
            first_sync_num_block = sync_num_RO_data[0]
            last_sync_num_block = sync_num_RO_data[-1]
    
        elif type(pqf) == str:
            f = h5py.File(pqf, 'r')
            sync_num_RO_data = f[sync_num_name]
            first_sync_num_block = sync_num_RO_data[0]
            last_sync_num_block = sync_num_RO_data[-1]

            f.close()
        else:
            print "Neither filepath nor file enetered in function please check:", pqf
            raise 
  
        # Check if a marked sync number is at the beginning or end of block to determine if some additional
        # analysis should be performed. 
        first_sync_num_with_marker = unique_sync_num_with_markers[0]
        first_sync_num_with_marker_2 = unique_sync_num_with_markers[0] -1
        last_sync_num_with_marker = unique_sync_num_with_markers[len(unique_sync_num_with_markers)-1]
        last_sync_num_with_marker_2 = unique_sync_num_with_markers[len(unique_sync_num_with_markers)-1] + 1

        if first_sync_num_with_marker == first_sync_num_block:
            print
            print "Something goes wrong with the first number, think of a way to fix this!"
            print
            #raise

        if first_sync_num_with_marker_2 == first_sync_num_block:
            print
            print "Something goes wrong with the first number, think of a way to fix this!"
            print
            #raise
            

        if last_sync_num_with_marker == last_sync_num_block:
            print
            print "Something goes wrong with the last number, think of a way to fix this!"
            print 
            #raise
            
        if last_sync_num_with_marker_2 == last_sync_num_block:
            print
            print "Something goes wrong with the last number, think of a way to fix this!"
            print 
            #raise

        SSRO_events = Analysis.Analyze_SSRO_data(PQ_sync_number, PQ_special, PQ_sync_time, \
                                        PQ_time, PQ_channel, RO_start, RO_length, marker_chan, \
                                        unique_sync_num_with_markers, chan_rnd_0, chan_rnd_1, sync_time_lim = sync_time_lim)

    return SSRO_events, PQ_sync_number, PQ_special, PQ_sync_time, PQ_time, PQ_channel

def get_SSRO_events_quick(pqf, marker_chan, RO_start, RO_length, sync_time_lim = 0, index = 1):
    """
    Returns an array with sync numbers in the first row, the number of photons in the readout window
    in the second column and the sync time of the first photon(the lowest sync time) in the third column.
    """

    sync_time_name = '/PQ_sync_time-' + str(index)
    sync_num_name = '/PQ_sync_number-' + str(index)
    spec_name = '/PQ_special-' + str(index)
    chan_name = '/PQ_channel-' + str(index)
    time_name = '/PQ_time-' + str(index)

    if type(pqf) == h5py._hl.files.File: 
        special_RO =pqf[spec_name].value
        channel_RO = pqf[chan_name].value

    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        special_RO = f[spec_name].value
        channel_RO = f[chan_name].value
        f.close()
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise 

    # checks if there are markers
    is_special = special_RO == 1
    is_channel = channel_RO == marker_chan
    is_mrkr = is_special & is_channel

    del is_special
    del is_channel
    
    num_mrkr = len(channel_RO[is_mrkr])

    # Get the indices of the marker events
    indices_mrkr = np.where(is_mrkr == True)[0]
    del is_mrkr
    marker_sync_numbers = np.empty((num_mrkr), dtype = np.uint32)

    # Get the list of sync numbers corresponding to markerrs
    for i,j in enumerate(indices_mrkr): 

        if type(pqf) == h5py._hl.files.File: 
            sync_num_RO_data = pqf[sync_num_name]
            sync_time_RO_data = pqf[sync_time_name]

            sync_num_RO_temp = sync_num_RO_data[j]
            sync_time_RO_temp = sync_time_RO_data[j]

        elif type(pqf) == str:
            f = h5py.File(pqf, 'r')
            sync_num_RO_data = f[sync_num_name]
            sync_time_RO_data = f[sync_time_name]

            sync_num_RO_temp = sync_num_RO_data[j]
            sync_time_RO_temp = sync_time_RO_data[j]

            f.close()

        # Substracts 1 of the sync number if the sync time is lower than the limit
        if sync_time_lim > 0:
            is_small_sync_time = sync_time_RO_temp <= sync_time_lim
            is_large_sync_time = sync_time_RO_temp > sync_time_lim

            if is_small_sync_time:
                marker_sync_numbers[i] = sync_num_RO_temp - 1
            elif is_large_sync_time:
                marker_sync_numbers[i] = sync_num_RO_temp
        else:
            marker_sync_numbers[i] = sync_num_RO_temp

    # Get a list with the unique sync numbers
    unique_sync_num_with_markers = np.unique(marker_sync_numbers)

    # Initializes a list with the indices of the sync numbers which have a marker
    indices_sync_num = []

    # Loops over the subblock of sync numbers to prevent memory errors
    if type(pqf) == h5py._hl.files.File: 
        sync_num_RO = pqf[sync_num_name].value

    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        sync_num_RO = f[sync_num_name].value
        f.close()
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise 

    is_sync_num_with_mrkr = np.array(([False] * len(sync_num_RO)))

    for i,j in enumerate(unique_sync_num_with_markers):
        is_sync_num_with_mrkr = is_sync_num_with_mrkr | (sync_num_RO == j)
        
    indices_sync_num.extend(np.where(is_sync_num_with_mrkr == True)[0])

    del sync_num_RO

    # Initializes arrays to save all PQ data
    PQ_sync_number = np.empty((0,), dtype = np.uint32) 
    PQ_special = np.empty((0,), dtype = np.uint32)         
    PQ_sync_time = np.empty((0,), dtype = np.uint64)
    PQ_time = np.empty((0,), dtype = np.uint64) 
    PQ_channel = np.empty((0,), dtype = np.uint32)


    # Get PQ data for all marked events in this blocks
    for i,j in enumerate(indices_sync_num):

        if type(pqf) == h5py._hl.files.File: 
            sync_num_RO_data = pqf[sync_num_name]
            sync_time_RO_data = pqf[sync_time_name]
            special_RO_data = pqf[spec_name]

            sync_num_RO = sync_num_RO_data[j]
            sync_time_RO = sync_time_RO_data[j]
            special_RO = special_RO_data[j]

            first_sync_num_block = sync_num_RO_data[0]
            last_sync_num_block = sync_num_RO_data[len(sync_num_RO_data)-1]

    
        elif type(pqf) == str:
            f = h5py.File(pqf, 'r')
            sync_num_RO_data = f[sync_num_name]
            sync_time_RO_data = f[sync_time_name]
            special_RO_data = f[spec_name]

            sync_num_RO = sync_num_RO_data[j]
            sync_time_RO = sync_time_RO_data[j]
            special_RO = special_RO_data[j]

            first_sync_num_block = sync_num_RO_data[0]
            last_sync_num_block = sync_num_RO_data[len(sync_num_RO_data)-1]

            f.close()
        else:
            print "Neither filepath nor file enetered in function please check:", pqf
            raise 

        # Check if a marked sync number is at the beginning or end of block to determine if some additional
        # analysis should be performed. 
        first_sync_num_with_marker = unique_sync_num_with_markers[0]
        first_sync_num_with_marker_2 = unique_sync_num_with_markers[0] -1
        last_sync_num_with_marker = unique_sync_num_with_markers[-1]
        last_sync_num_with_marker_2 = unique_sync_num_with_markers[-1] + 1

        if first_sync_num_with_marker == first_sync_num_block:
            print
            print "Something goes wrong with the first number, think of a way to fix this!"
            print
            raise

        if first_sync_num_with_marker_2 == first_sync_num_block:
            print
            print "Something goes wrong with the first number, think of a way to fix this!"
            print
            raise

        if last_sync_num_with_marker == last_sync_num_block:
            print
            print "Something goes wrong with the last number, think of a way to fix this!"
            print
            raise

        if last_sync_num_with_marker_2 == last_sync_num_block:
            print
            print "Something goes wrong with the last number, think of a way to fix this!"
            print
            raise


        # Concatenates the PQ data for several sync numbers
        PQ_sync_number = np.hstack((PQ_sync_number, sync_num_RO)) 
        PQ_special = np.hstack((PQ_special, special_RO))         
        PQ_sync_time = np.hstack((PQ_sync_time, sync_time_RO)) 

    Quick_SSRO_events = np.empty((0,3))
    is_ph_RO = PQ_special == 0   
    is_in_window = (RO_start  <= PQ_sync_time) & (PQ_sync_time < (RO_start + RO_length))
    is_ph_RO_in_ro_window = is_in_window & is_ph_RO

    for i,s in enumerate(unique_sync_num_with_markers):

        is_sync_num_s = PQ_sync_number == s
        is_photons_RO = is_sync_num_s & is_ph_RO_in_ro_window

        sync_time_RO_photons = PQ_sync_time[is_photons_RO]

        num_phot = len(sync_time_RO_photons)
        if len(sync_time_RO_photons) > 0:
            arr_time_first_phot = min(sync_time_RO_photons)
        else:
            arr_time_first_phot = 0

        _event = np.array([s, num_phot, arr_time_first_phot])
        Quick_SSRO_events = np.vstack((Quick_SSRO_events, _event))

    return Quick_SSRO_events


############################ Re analyze SSRO data ################################################

def re_analyze_SSRO_data(pqf, marker_chan ,RO_start, RO_length, chan_rnd_0, chan_rnd_1, sync_time_lim = 0,  VERBOSE = True ):
    # Define all block names
    sync_time_name = '/PQ_sync_time'
    sync_num_name = '/PQ_sync_number'
    spec_name = '/PQ_special'
    chan_name = '/PQ_channel'
    time_name = '/PQ_time'


    if type(pqf) == h5py._hl.files.File: 
        PQ_sync_number = pqf[sync_num_name].value
        PQ_sync_time = pqf[sync_time_name].value
        PQ_time = pqf[time_name].value
        PQ_channel = pqf[chan_name].value
        PQ_special = pqf[spec_name].value


    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        PQ_sync_number = f[sync_num_name].value
        PQ_sync_time = f[sync_time_name].value
        PQ_time = f[time_name].value
        PQ_channel = f[chan_name].value
        PQ_special = f[spec_name].value
        f.close()
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise 
    # Initializes an array to save all SSRO data
    SSRO_events = np.empty((0,30), dtype = np.uint64)

    _a = get_attributes(pqf)

    unique_sync_numbers = np.unique(PQ_sync_number)

    unique_sync_num_with_markers = unique_sync_numbers[1::3]



    # Create a filter which is True for all photons
    is_ph_RO = (PQ_special == 0 ) & (PQ_channel == 0)

    # Filters if sync times are in the defined read out window
    is_in_window = (RO_start  <= PQ_sync_time) & (PQ_sync_time < (RO_start + RO_length))
    
    # Defines a filter for photons in the readout window
    is_ph_RO_in_ro_window = is_in_window & is_ph_RO

    is_special = PQ_special == 1
    is_channel_rnd_0 = PQ_channel == chan_rnd_0
    is_channel_rnd_1 = PQ_channel == chan_rnd_1

    is_rnd = (PQ_special == 0 ) & (PQ_channel == 1)
    is_rnd_0 = is_special & is_channel_rnd_0
    is_rnd_1 = is_special & is_channel_rnd_1

    is_PQ_mrkr =  (PQ_special == 1)  & (PQ_channel == marker_chan)
    print np.sum(is_PQ_mrkr)
    print len(is_PQ_mrkr)

    # Loop over all sync numbers with markers
    for i,s in enumerate(unique_sync_num_with_markers):

        # Create a filter which filters on a specific sync number
        is_sync_num_s = PQ_sync_number == s
        is_sync_num_s_plus_1 = PQ_sync_number == (s + 1)
        print s + 1
        print PQ_sync_number[is_sync_num_s_plus_1]
        print len(is_sync_num_s_plus_1)
        print is_sync_num_s_plus_1 
        print is_PQ_mrkr
        print np.sum(is_sync_num_s_plus_1)

        # Creates a filter which filterse if there is a photon in the readout window for a specific sync number
        is_photons_RO = is_sync_num_s & is_ph_RO_in_ro_window

        # Gets the sync times for photons in the readout window corresponding to a certain sync number
        sync_time_RO_photons = PQ_sync_time[is_photons_RO]

        # Defines absolute time of the marker
        sync_time_mrkr_s = PQ_sync_time[(is_sync_num_s & is_PQ_mrkr)]
        print sync_time_mrkr_s
        abs_time_mrkr_s = PQ_time[(is_sync_num_s & is_PQ_mrkr)]
        sync_time_mrkr_s_plus_1 = PQ_sync_time[(is_sync_num_s_plus_1 & is_PQ_mrkr)]
        print sync_time_mrkr_s_plus_1
        abs_time_mrkr_s_plus_1 = PQ_time[(is_sync_num_s_plus_1 & is_PQ_mrkr)]

        if (len(sync_time_mrkr_s) > 0) & (sync_time_mrkr_s > sync_time_lim):
            abs_time_mrkr = abs_time_mrkr_s
        elif (sync_time_mrkr_s_plus_1 <= sync_time_lim):
            abs_time_mrkr = abs_time_mrkr_s_plus_1 
        else:
            print "Something strange is hapening with the markers or I programmed it bad"
            raise

        # Makes two boolean list for random channel 1 and 2 (named 0 & 1) and filters them on the sync number
        # One of these list should be False completely and the other one should be True once, indicating there
        # is one random number for each markers
        rnd_0 = is_rnd_0[is_sync_num_s]
        rnd_1 = is_rnd_1[is_sync_num_s]
        rnd_generated = is_rnd[is_sync_num_s]

        # Filters for the sync time are created, there should be only one True in both filters which gives
        # the sync time of the random number
        is_sync_time_rnd = is_rnd & is_sync_num_s

        # Checks if it is a random number in marker channel 1
        if sum(rnd_generated) > 0:
            if (sum(rnd_0) == 1):
                # States that ther is a random number, that it is zero and gets the sync time of this number
                rnd_gen_check = 1
                rnd_num = 0
                sync_time_rnd_num = PQ_sync_time[is_sync_time_rnd][-1]
            # Checks if it is a random number in marker channel 2
            elif (sum(rnd_1) == 1):
                # States that ther is a random number, that it is zero and gets the sync time of this number
                rnd_gen_check = 1
                rnd_num = 1
                sync_time_rnd_num = PQ_sync_time[is_sync_time_rnd][-1]
            else:
                print "Somehow there is now marker generated, look at the RND number generator"
                print"========================================================================"
                print 
                raise
        # States that no random number 
        else:
            print "There are events for which no random number is generated"
            rnd_gen_check = 0
            rnd_num = 2
            sync_time_rnd_num = 0

        # Define the number of readout photons
        num_phot = len(sync_time_RO_photons)

        # Makes the array for the arrival times of the photons for different numbers of photons
        if (len(sync_time_RO_photons) > 0) & (len(sync_time_RO_photons) == 24):
            arr_times = sync_time_RO_photons
        elif (len(sync_time_RO_photons) > 0) & (len(sync_time_RO_photons) < 24):
            zero_addition = np.zeros((24-len(sync_time_RO_photons),), dtype = np.uint64)
            arr_times = np.concatenate((sync_time_RO_photons,zero_addition))
        else:
            arr_times = np.zeros((24,), dtype = np.uint64)

        # Stacks all SSRO events for different sync numbers
        _event = np.concatenate((np.array([s, abs_time_mrkr, num_phot, rnd_gen_check, rnd_num, sync_time_rnd_num], dtype = np.uint64) , arr_times))
        SSRO_events = np.vstack((SSRO_events, _event))

    tb.clear_analysis_data(pqf, 'Total_SSRO_events')
    tb.set_analysis_data(pqf, 'Total_SSRO_events', SSRO_events, _a) 

    return SSRO_events

###################### Subfunctions ############################################################################
def get_attributes(pqf):
    """
    Get the attributes for the SSRO events from LT3 and LT4
    """
    if 'lt3' in pqf:
        columns = "Sync_Number_LT3,Absolute time marker sync num LT3,Number of photons_LT3,\
        Random Number Indicator LT3,Random Number LT3,\
        Sync_Time Random Number LT3,Sync_Time_photon_1_LT3,Sync_Time_photon_2_LT3,Sync_Time_photon_3_LT3,\
        Sync_Time_photon_4_LT3,Sync_Time_photon_5_LT3,Sync_Time_photon_6_LT3,Sync_Time_photon_7_LT3,\
        Sync_Time_photon_8_LT3,Sync_Time_photon_9_LT3,Sync_Time_photon_10_LT3,Sync_Time_photon_11_LT3,\
        Sync_Time_photon_12_LT3,Sync_Time_photon_13_LT3,Sync_Time_photon_14_LT3,Sync_Time_photon_15_LT3,\
        Sync_Time_photon_16_LT3,Sync_Time_photon_17_LT3,Sync_Time_photon_18_LT3,Sync_Time_photon_19_LT3,\
        Sync_Time_photon_20_LT3,Sync_Time_photon_21_LT3,Sync_Time_photon_22_LT3,Sync_Time_photon_23_LT3,\
        Sync_Time_photon_24_LT3"
    elif 'lt4' in pqf:
        columns = "Sync_Number_LT4,Absolute time marker sync num LT4,Number of photons_LT4,\
        Random Number Indicator LT4,Random Number LT4,\
        Sync_Time Random Number LT4,Sync_Time_photon_1_LT4,Sync_Time_photon_2_LT4,Sync_Time_photon_3_LT4,\
        Sync_Time_photon_4_LT4,Sync_Time_photon_5_LT4,Sync_Time_photon_6_LT4,Sync_Time_photon_7_LT4,\
        Sync_Time_photon_8_LT4,Sync_Time_photon_9_LT4,Sync_Time_photon_10_LT4,Sync_Time_photon_11_LT4,\
        Sync_Time_photon_12_LT4,Sync_Time_photon_13_LT4,Sync_Time_photon_14_LT4,Sync_Time_photon_15_LT4,\
        Sync_Time_photon_16_LT4,Sync_Time_photon_17_LT4,Sync_Time_photon_18_LT4,Sync_Time_photon_19_LT4,\
        Sync_Time_photon_20_LT4,Sync_Time_photon_21_LT4,Sync_Time_photon_22_LT4,Sync_Time_photon_23_LT4,\
        Sync_Time_photon_24_LT4"
    else: 
        columns = ""

    _a = {'Columns': columns}

    return _a

def get_unique_sync_with_mrkr_quick(pqf, is_mrkr,sync_num_name, sync_time_name, sync_time_lim, num_mrkr):
    # Get the indices of the marker events
    indices_mrkr = np.where(is_mrkr == True)[0]
    marker_sync_numbers = np.empty((num_mrkr), dtype = np.uint32)

    # Get the list of sync numbers corresponding to markerrs
    for i,j in enumerate(indices_mrkr): 

        if type(pqf) == h5py._hl.files.File: 
            sync_num_RO_data = pqf[sync_num_name]
            sync_time_RO_data = pqf[sync_time_name]

            sync_num_RO = sync_num_RO_data[j]
            sync_time_RO = sync_time_RO_data[j]

        elif type(pqf) == str:
            f = h5py.File(pqf, 'r')
            sync_num_RO_data = f[sync_num_name]
            sync_time_RO_data = f[sync_time_name]

            sync_num_RO = sync_num_RO_data[j]
            sync_time_RO = sync_time_RO_data[j]

            f.close()

        # Substracts 1 of the sync number if the sync time is lower than the limit
        if sync_time_lim > 0:
            is_small_sync_time = sync_time_RO <= sync_time_lim
            is_large_sync_time = sync_time_RO > sync_time_lim

            if is_small_sync_time:
                marker_sync_numbers[i] = sync_num_RO - 1
            elif is_large_sync_time:
                marker_sync_numbers[i] = sync_num_RO
        else:
            marker_sync_numbers[i] = sync_num_RO

    return marker_sync_numbers

def get_indices_PQ_save_sync_num(pqf, PQ_save_sync_numbers, sync_num_name, num_sync_num):
    """
    Returns the indices of all data points corresponding to the PQ_save_sync_numbers
    """
    # Spits the length of the block up in 10 parts
    len_RO_num = num_sync_num/10
    if len_RO_num > 20000000:
        print "Something might go wrong with the memory"
        print "The length of the array opened is", len_RO_num

    # Initializes a list with the indices of the sync numbers which have a marker
    indices_PQ_save_sync_num = []

    # Loops over the subblock of sync numbers to prevent memory errors
    for i in range(10):
        start = i * len_RO_num
        if i != 9:
            stop = (i+1) * len_RO_num
        else:
            stop = num_sync_num

        if type(pqf) == h5py._hl.files.File: 
            sync_num_RO = pqf[sync_num_name]
            sync_num_RO_block = sync_num_RO[start:stop]

        elif type(pqf) == str:
            f = h5py.File(pqf, 'r')
            sync_num_RO = f[sync_num_name]
            sync_num_RO_block = sync_num_RO[start:stop]
            f.close()
        else:
            print "Neither filepath nor file enetered in function please check:", pqf
            raise 

        for i,j in enumerate(PQ_save_sync_numbers):
            indices_PQ_save_sync_num_j = np.where(sync_num_RO_block == j)[0] + start
            indices_PQ_save_sync_num.extend(indices_PQ_save_sync_num_j)
    return indices_PQ_save_sync_num

def get_PQ_save_data(pqf,indices_PQ_save_sync_num, index = 1):
    """
    Returns the PQ data for all indices corresponding to the sync numbers we want to save.
    """
    # Define all block names
    sync_time_name = '/PQ_sync_time-' + str(index)
    sync_num_name = '/PQ_sync_number-' + str(index)
    spec_name = '/PQ_special-' + str(index)
    chan_name = '/PQ_channel-' + str(index)
    time_name = '/PQ_time-' + str(index)

    # Initializes arrays to save all PQ data
    PQ_sync_number = np.empty((0,), dtype = np.uint32) 
    PQ_special = np.empty((0,), dtype = np.uint32)         
    PQ_sync_time = np.empty((0,), dtype = np.uint64)
    PQ_time = np.empty((0,), dtype = np.uint64) 
    PQ_channel = np.empty((0,), dtype = np.uint32)

    # Get PQ data for all marked events in this blocks
    for i,j in enumerate(indices_PQ_save_sync_num):
        if type(pqf) == h5py._hl.files.File: 
            sync_num_RO_data = pqf[sync_num_name]
            sync_time_RO_data = pqf[sync_time_name]
            time_RO_data = pqf[time_name]
            channel_RO_data = pqf[chan_name]
            special_RO_data = pqf[spec_name]

            sync_num_RO = sync_num_RO_data[j]
            sync_time_RO = sync_time_RO_data[j]
            time_RO = time_RO_data[j]
            channel_RO = channel_RO_data[j]
            special_RO = special_RO_data[j]
    
        elif type(pqf) == str:
            f = h5py.File(pqf, 'r')
            sync_num_RO_data = f[sync_num_name]
            sync_time_RO_data = f[sync_time_name]
            time_RO_data = f[time_name]
            channel_RO_data = f[chan_name]
            special_RO_data = f[spec_name]

            sync_num_RO = sync_num_RO_data[j]
            sync_time_RO = sync_time_RO_data[j]
            time_RO = time_RO_data[j]
            channel_RO = channel_RO_data[j]
            special_RO = special_RO_data[j]

            f.close()
        else:
            print "Neither filepath nor file enetered in function please check:", pqf
            raise 

        # Concatenates the PQ data for several sync numbers
        PQ_sync_number = np.hstack((PQ_sync_number, sync_num_RO)) 
        PQ_special = np.hstack((PQ_special, special_RO))         
        PQ_sync_time = np.hstack((PQ_sync_time, sync_time_RO)) 
        PQ_time = np.hstack((PQ_time, time_RO))   
        PQ_channel = np.hstack((PQ_channel, channel_RO))

    return PQ_sync_number, PQ_special, PQ_sync_time, PQ_time, PQ_channel