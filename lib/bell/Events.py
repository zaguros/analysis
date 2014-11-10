import numpy as np
import h5py
import time
from datetime import datetime
from analysis.lib.tools import toolbox as tb
from analysis.lib.pq import pq_tools, pq_plots


def get_entanglement_events(fp_BS,fp_LT3,fp_LT4,chan, i, first_win_min, 
        first_win_max, second_win_min, second_win_max, force_eval=False, VERBOSE = True):
    
    """
    Returns either the entanglement events already saved; the corresponding attributes
    and save = False, or returns the newly calculated entanglement events;
    the corresponding attributes (column names) and the save = True. Put in the file path 
    of the BS then of LT3 then of LT4 then the marker channel and then the repitition if
    it's looped. Also put in the windows to determine if an photon is the first photon to arrive.
    (first_win_min,first_win_max,second_win_min, second_win_max)
    """

    
    if tb.has_analysis_data(fp_BS, 'Entanglement_events') and not force_eval:
        tev, _a = tb.get_analysis_data(fp_BS, 'Entanglement_events')
        unique_sync_num_with_markers = tev[:,0]
        if VERBOSE:
            print
            string = 'The number of events with PLU markers in run ' + str(i+1) + ' is:'
            print string, len(unique_sync_num_with_markers)
            print
            print 'Found {} saved valid entanglement events.'.format(int(len(tev)))
            print '===================================='
            print
        return tev, _a 

    # Opens beamsplitter data 
    f = h5py.File(fp_BS, 'r')
    sync_times = f['/PQ_sync_time-1'].value * 1e-3 # we prefer ns over ps
    sync_numbers = f['/PQ_sync_number-1'].value
    Channel = f['/PQ_channel-1'].value
    abs_times = f['/PQ_time-1'].value
    f.close()
        
    # Opens LT3 data  
    g = h5py.File(fp_LT3,'r')
    for k in g.keys():
        if type(g[k])==h5py.Group:
            ad3_reps = g[('/'+ str(k) + '/ssro/entanglement_events')].value
            ad3_ssro = g[('/'+ str(k) + '/ssro/RO_data')].value
            ad3_CR_before = g[('/'+ str(k) + '/ssro/CR_before')].value
            ad3_CR_after = g[('/'+ str(k) + '/ssro/CR_after')].value
    g.close()
    
    # Opens LT4data    
    h = h5py.File(fp_LT4,'r')
    for k in h.keys():
        if type(h[k])==h5py.Group:
            ad4_reps = h[('/'+ str(k) + '/ssro/entanglement_events')].value
            ad4_ssro = h[('/'+ str(k) + '/ssro/RO_data')].value
            ad4_CR_before = h[('/'+ str(k) + '/ssro/CR_before')].value
            ad4_CR_after = h[('/'+ str(k) + '/ssro/CR_after')].value
    h.close()    
    

    sync_num_with_markers = sync_numbers[pq_tools.filter_marker(fp_BS,chan)]
    print sync_num_with_markers
    unique_sync_num_with_markers = np.unique(sync_num_with_markers)
    print unique_sync_num_with_markers

    if VERBOSE:

        print 
        string = 'The number of events with PLU markers in run ' + str(i+1) + ' is:'
        print string, len(unique_sync_num_with_markers)
        print
        print 'Adwin LT3'
        print '---------'
        print 'Number of events:', ad3_reps,
        if ad3_reps != len(unique_sync_num_with_markers):
            print 'number of Adwin LT3 events does not match the PLU marker \
                    events - data set seems faulty :('
        else:
            print 'OK :)'
            
        print 
        print 'Adwin LT4'
        print '---------'
        print 'Number of events:', ad4_reps
        if ad4_reps != len(unique_sync_num_with_markers):
            print 'number of Adwin LT3 events does not match the PLU marker \
                    events - data set seems faulty :('
        else:
            print 'OK :)'
    
    # Gets filters for photons with markers in the first and second window
    # from the Filter file
    is_photon_1st_window_with_markers, is_photon_2nd_window_with_markers =\
                                    pq_tools.get_photons_with_markers(fp_BS,chan,
                                        first_win_min, first_win_max, second_win_min, second_win_max)

    # Retrieves sync numbers and sync times for photons both in the first
    # and 2nd window
    Sync_num_1st_window_with_markers = sync_numbers[is_photon_1st_window_with_markers]
    Channel_1st_window_with_markers = Channel[is_photon_1st_window_with_markers]#XXXXXXXXXXXXXXX
    Sync_times_1st_window_with_markers = sync_times[is_photon_1st_window_with_markers]
    Sync_num_2nd_window_with_markers = sync_numbers[is_photon_2nd_window_with_markers]
    Channel_2nd_window_with_markers =  Channel[is_photon_2nd_window_with_markers]
    Sync_times_2nd_window_with_markers = sync_times[is_photon_2nd_window_with_markers]
   
   
    # Defines a filter for all events with markers
    is_all_markers = is_photon_1st_window_with_markers | is_photon_2nd_window_with_markers
    
    # Gets the absolute times for all events with makers
    PLU_mrkr_abs_times = abs_times[is_all_markers]
    
    #Initializes the final array of entanglement events
    entanglement_events = np.empty((0,14))

    # Get all real entanglement events, loops over sync numbers
    for i,s in enumerate(unique_sync_num_with_markers):
        
        # The attempt is defined as the sync number modulo 250 
        #(250 = the repitition rate)
        attempt = np.mod(s,250)
        
        # Return filters for specific sync number s
        is_ph_1st_win_sync_num_s = Sync_num_1st_window_with_markers == s
        is_ph_2nd_win_sync_num_s = Sync_num_2nd_window_with_markers == s
        
        # Test if there is one photon in both windows
        if len(Sync_num_1st_window_with_markers[is_ph_1st_win_sync_num_s]) == 1 \
            and len(Sync_num_2nd_window_with_markers[is_ph_2nd_win_sync_num_s]) == 1:

            # Saves sync times an channels of both photons
            stimes = np.array([ Sync_times_1st_window_with_markers[is_ph_1st_win_sync_num_s],\
                Sync_times_2nd_window_with_markers[is_ph_2nd_win_sync_num_s]]).reshape(-1)
            channel_1 = Channel_1st_window_with_markers[is_ph_1st_win_sync_num_s]
            channel_2 = Channel_2nd_window_with_markers[is_ph_2nd_win_sync_num_s]
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
        Sync Nymber | Sync Time Photon 1 | Sync Time photon 2 | Photon 1 Channel | 
        Photon 2 Channel | Attempts | Amout of Photons LT3 | Amount of Photons LT 3 | 
        CR Check Before LT3 | CR Check After LT3 | CR Check Before LT4 | 
        CR Check After LT4 | psiminus | absolute time
        """
        
        columns = "Sync_Number, Sync_Time_photon_1, Sync_Time_photon_2, Channel_photon_1,\
        Channel_photon_2, Attempts, Amount_of_ph_LT3, CR_check_before_LT3,CR_check_after_LT3,\
        Amount_of_ph_LT4, CR_check_before_LT4, CR_check_after_LT4, psiminus, abs_time"

        _a = {'Columns': columns}
                
        _event = np.array([s, 
                        stimes[0],
                        stimes[1], 
                        chans[0], 
                        chans[1], 
                        attempt, 
                        ad3_ssro[i],
                        ad4_ssro[i],
                        ad3_CR_before[i],
                        ad3_CR_after[i],
                        ad4_CR_before[i],
                        ad4_CR_after[i], 
                        psiminus, 
                        PLU_mrkr_abs_times[i]])
                        
        entanglement_events = np.vstack((entanglement_events, _event))

    if VERBOSE:
        print
        print 'Found {} valid entanglement events.'.format(int(len(entanglement_events)))
        print '===================================='
        print
    
    return entanglement_events, _a

#######################  SSRO events #################################

def get_total_SSRO_events(pqf, RO_start, marker_chan, chan_rnd_0, chan_rnd_1, sync_time_lim, VERBOSE = True):
    """
    Returns SSRO data for all marked events. 
    Colums are:
    Sync Nymber | number of photons | RND number indicator | RND number | Sync Time RND number | Sync Times photon 1-24 |
    """
    
    columns = "Sync_Number, Number of photons, Random Number Indicator, Random Number, Sync_Time Random Number, Sync_Time_photon_1, \
    Sync_Time_photon_2, Sync_Time_photon_3, Sync_Time_photon_4, Sync_Time_photon_5, Sync_Time_photon_6, \
    Sync_Time_photon_7, Sync_Time_photon_8, Sync_Time_photon_9, Sync_Time_photon_10, Sync_Time_photon_11, \
    Sync_Time_photon_12, Sync_Time_photon_13, Sync_Time_photon_14, Sync_Time_photon_15, Sync_Time_photon_16,\
    Sync_Time_photon_17, Sync_Time_photon_18, Sync_Time_photon_19, Sync_Time_photon_20, Sync_Time_photon_21, \
    Sync_Time_photon_22, Sync_Time_photon_23, Sync_Time_photon_24"
    _a = {'Columns': columns}

    # Gets the number of blocks in the data
    num_blocks = tb.get_num_blocks(pqf)

    if VERBOSE:
        print 'The total number of blocks is:', num_blocks

    # Initializes arrays to save the PQ-data
    PQ_sync_number = np.empty((0,), dtype = int) 
    PQ_special = np.empty((0,), dtype = int)         
    PQ_sync_time = np.empty((0,), dtype = long)
    PQ_time = np.empty((0,), dtype = long)      
    PQ_channel = np.empty((0,), dtype = int)

    # Initializes an array to save the SSRO data
    total_SSRO_events = np.empty((0,29))


    # Loops over every block
    for i in range(num_blocks):

        print "Main Loop", datetime.now()
        # Get a list of sync numbers which correspond to events for which a marker has arrived. 
        unique_sync_num_with_markers = \
            pq_tools.get_un_sync_num_with_markers(pqf, marker_chan, sync_time_lim = sync_time_lim, index = i+1, VERBOSE = VERBOSE)
        
        # Get the SSRO events and PQ data for these sync numbers
        _events, _PQ_sync_number, _PQ_special, _PQ_sync_time, _PQ_time, _PQ_channel = \
                get_SSRO_events(pqf, unique_sync_num_with_markers, RO_start, chan_rnd_0, chan_rnd_1, index = i+1)
                    
        
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

    if VERBOSE:
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
    
    columns = "Sync_Number, Number of photons,  Sync_Time_photon_1"
    _a = {'Columns': columns}

    # Gets the number of blocks for the data
    num_blocks = tb.get_num_blocks(pqf)

    if VERBOSE:
        print 'The total number of blocks is:', num_blocks

    # Initializes the array to save all SSRO events
    total_SSRO_events = np.empty((0,3))

    # Loop over all blocks
    for i in range(num_blocks):

        # Returns a list with the sync numbers for which a marker is observed
        unique_sync_num_with_markers = \
            pq_tools.get_un_sync_num_with_markers(pqf, marker_chan, sync_time_lim = sync_time_lim, index = i+1, VERBOSE = VERBOSE)
        
        # Gets all events for a block
        _events = get_SSRO_events_quick(pqf, unique_sync_num_with_markers, RO_start, RO_length, index = i+1)

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


def get_SSRO_events(pqf, unique_sync_num_with_markers,RO_start, chan_rnd_0, chan_rnd_1, index = 1):
    """
    Returns an array with sync numbers in the first column, the number of photons in the readout window
    in the second column, a random number generation check in the third colum (1 if one is generated 0 if not),
    the random number itself in the fourth column, the sync time of the random number in the fifth column, 
    and the sync times of the 1st to the 24th photon.
    """

    print "Before opening the files", datetime.now()

    # Define all block names
    sync_time_name = '/PQ_sync_time-' + str(index)
    sync_num_name = '/PQ_sync_number-' + str(index)
    spec_name = '/PQ_special-' + str(index)
    chan_name = '/PQ_channel-' + str(index)
    time_name = '/PQ_time-' + str(index)

    if type(pqf) == h5py._hl.files.File: 
        sync_numbers = pqf[sync_num_name].value
        special_RO =pqf[spec_name].value
        sync_time_RO =pqf[sync_time_name].value
        time_RO = pqf[time_name].value
        channel_RO = pqf[chan_name].value

        # Get name of the group to find read out length
        group = tb.get_msmt_name(pqf)
        total_string_name = '/' + group + '/joint_params'
        RO_length =pqf[total_string_name].attrs['LDE_RO_duration']  * 1e9

    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        sync_num_RO = f[sync_num_name].value
        special_RO = f[spec_name].value
        sync_time_RO = f[sync_time_name].value
        time_RO = f[time_name].value
        channel_RO = f[chan_name].value

        # Get name of the group to find read out length
        group = tb.get_msmt_name(pqf)
        total_string_name = '/' + group + '/joint_params'
        RO_length = f[total_string_name].attrs['LDE_RO_duration']  * 1e9
        f.close()
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise 

    print "Files opened", datetime.now()

    # Initializes an array to save all SSRO data
    SSRO_events = np.empty((0,29))

    # Initializes arrays to save all PQ data
    PQ_sync_number = np.empty((0,), dtype = int) 
    PQ_special = np.empty((0,), dtype = int)         
    PQ_sync_time = np.empty((0,), dtype = long)
    PQ_time = np.empty((0,), dtype = long)      
    PQ_channel = np.empty((0,), dtype = int)

    # Create a filter which is True for all photons
    is_ph_RO = special_RO == 0   

    # Filters if sync times are in the defined read out window
    is_in_window = (RO_start  <= sync_time_RO) & (sync_time_RO < (RO_start + RO_length))
    
    # Defines a filter for photons in the readout window
    is_ph_RO_in_ro_window = is_in_window & is_ph_RO

    # Gets filters for if an event is a random number in the first channel or in the second channel
    # Channels are inputs for the function
    is_rnd_0, is_rnd_1 = pq_tools.get_rndm_num(pqf, chan_rnd_0, chan_rnd_1, index = index)

    print "Check if something goes wrong with makers", datetime.now()

    if len(unique_sync_num_with_markers) > 0:
        first_sync_num_with_marker = unique_sync_num_with_markers[0]
        first_sync_num_with_marker_2 = unique_sync_num_with_markers[0] -1
        irst_sync_num_with_marker = unique_sync_num_with_markers[0]
        last_sync_num_with_marker = unique_sync_num_with_markers[len(unique_sync_num_with_markers)-1]
        last_sync_num_with_marker_2 = unique_sync_num_with_markers[len(unique_sync_num_with_markers)-1] + 1
        first_sync_num_block = sync_num_RO[0]
        last_sync_num_block = sync_num_RO[len(sync_num_RO)-1]

        if first_sync_num_with_marker == first_sync_num_block:
            print 
            print
            print
            print "Something goes wrong with the first number, think of a way to fix this!"
            print 
            print
            print

        if first_sync_num_with_marker_2 == first_sync_num_block:
            print 
            print
            print
            print "Something goes wrong with the first number, think of a way to fix this!"
            print 
            print
            print
            

        if last_sync_num_with_marker == last_sync_num_block:
            print 
            print
            print
            print "Something goes wrong with the last number, think of a way to fix this!"
            print 
            print
            print
            
        if last_sync_num_with_marker_2 == last_sync_num_block:
            print 
            print
            print
            print "Something goes wrong with the last number, think of a way to fix this!"
            print 
            print
            print

    print "Check finished", datetime.now()


    # Loop over all sync numbers with markers
    for i,s in enumerate(unique_sync_num_with_markers):
        print "Loop over marked events", datetime.now()

        # Create a filter which filters on a specific sync number
        is_sync_num_s = sync_num_RO == s

        # Creates a filter which filterse if there is a photon in the readout window for a specific sync number
        is_photons_RO = is_sync_num_s & is_ph_RO_in_ro_window

        # Gets the sync times for photons in the readout window corresponding to a certain sync number
        sync_time_RO_photons = sync_time_RO[is_photons_RO]

        # Makes two boolean list for random channel 1 and 2 (named 0 & 1) and filters them on the sync number
        # One of these list should be False completely and the other one should be True once, indicating there
        # is one random number for each markers
        rnd_0 = is_rnd_0[is_sync_num_s]
        rnd_1 = is_rnd_1[is_sync_num_s]
        
        # Filters for the sync time are created, there should be only one True in both filters which gives
        # the sync time of the random number
        is_sync_time_rnd_num_0 = is_rnd_0 & is_sync_num_s
        is_sync_time_rnd_num_1 = is_rnd_1 & is_sync_num_s

        # Checks if it is a random number in marker channel 1
        if (sum(rnd_0) == 1):
            # States that ther is a random number, that it is zero and gets the sync time of this number
            rnd_gen_check = 1
            rnd_num = 0
            sync_time_rnd_num = sync_time_RO[is_sync_time_rnd_num_0]
        # Checks if it is a random number in marker channel 2
        elif (sum(rnd_1) == 1):
            # States that ther is a random number, that it is zero and gets the sync time of this number
            rnd_gen_check = 1
            rnd_num = 1
            sync_time_rnd_num = sync_time_RO[is_sync_time_rnd_num_1]
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
            zero_addition = np.zeros(24-len(sync_time_RO_photons))
            arr_times = np.concatenate((sync_time_RO_photons,zero_addition))
        else:
            arr_times = np.zeros(24)

        # Get the PQ data for a specific sync number
        _PQ_sync_number = sync_num_RO[is_sync_num_s]
        _PQ_special = special_RO[is_sync_num_s]
        _PQ_sync_time = sync_time_RO[is_sync_num_s]
        _PQ_time = time_RO[is_sync_num_s]
        _PQ_channel = channel_RO[is_sync_num_s]

        # Concatenates the PQ data for several sync numbers
        PQ_sync_number = np.hstack((PQ_sync_number,_PQ_sync_number)) 
        PQ_special = np.hstack((PQ_special, _PQ_special))         
        PQ_sync_time = np.hstack((PQ_sync_time, _PQ_sync_time)) 
        PQ_time = np.hstack((PQ_time, _PQ_time))      
        PQ_channel = np.hstack((PQ_channel, _PQ_channel)) 

        # Stacks all SSRO events for different sync numbers
        _event = np.concatenate((np.array([s, num_phot, rnd_gen_check, rnd_num, sync_time_rnd_num]) , arr_times))
        SSRO_events = np.vstack((SSRO_events, _event))

    return SSRO_events, PQ_sync_number, PQ_special, PQ_sync_time, PQ_time, PQ_channel



def get_SSRO_events_quick(pqf, unique_sync_num_with_markers,RO_start, RO_length, index = 1):
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
        sync_numbers = pqf[sync_num_name].value
        special_RO =pqf[spec_name].value
        sync_time_RO =pqf[sync_time_name].value
        time_RO = pqf[time_name].value
        channel_RO = pqf[chan_name].value


        # Get name of the group to find read out length
        group = tb.get_msmt_name(pqf)
        total_string_name = '/' + group + '/joint_params'
        #RO_length =pqf[total_string_name].attrs['LDE_RO_duration']  * 1e9

    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        sync_num_RO = f[sync_num_name].value
        special_RO = f[spec_name].value
        sync_time_RO = f[sync_time_name].value
        time_RO = f[time_name].value
        channel_RO = f[chan_name].value

        # Get name of the group to find read out length
        group = tb.get_msmt_name(pqf)
        total_string_name = '/' + group + '/joint_params'
        #RO_length = f[total_string_name].attrs['LDE_RO_duration']  * 1e9
        f.close()
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise 

    Quick_SSRO_events = np.empty((0,3))
    is_ph_RO = special_RO == 0   
    is_in_window = (RO_start  <= sync_time_RO) & (sync_time_RO < (RO_start + RO_length))
    is_ph_RO_in_ro_window = is_in_window & is_ph_RO

    for i,s in enumerate(unique_sync_num_with_markers):

        is_sync_num_s = sync_num_RO == s
        is_photons_RO = is_sync_num_s & is_ph_RO_in_ro_window
        sync_time_RO_photons = sync_time_RO[is_photons_RO]

        num_phot = len(sync_time_RO_photons)
        if len(sync_time_RO_photons) > 0:
            arr_time_first_phot = min(sync_time_RO_photons)
        else:
            arr_time_first_phot = 0

        _event = np.array([s, num_phot, arr_time_first_phot])
        Quick_SSRO_events = np.vstack((Quick_SSRO_events, _event))


    return Quick_SSRO_events