import numpy as np
import h5py
from analysis.lib.tools import toolbox as tb
from analysis.lib.pq import pq_tools, pq_plots

def get_entanglement_events(fp_BS,fp_LT1,fp_LT3,chan, i, first_win_min, 
        first_win_max, second_win_min, second_win_max, force_eval=False, VERBOSE = True):
    
    """
    Returns either the entanglement events already saved; the corresponding attributes
    and save = False, or returns the newly calculated entanglement events;
    the corresponding attributes (column names) and the save = True. Put in the file path 
    of the BS then of LT1 then of LT3 then the marker channel and then the repitition if
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
        
    # Opens LT1 data  
    g = h5py.File(fp_LT1,'r')
    for k in g.keys():
        if type(g[k])==h5py.Group:
            ad1_reps = g[('/'+ str(k) + '/ssro/entanglement_events')].value
            ad1_ssro = g[('/'+ str(k) + '/ssro/RO_data')].value
            ad1_CR_before = g[('/'+ str(k) + '/ssro/CR_before')].value
            ad1_CR_after = g[('/'+ str(k) + '/ssro/CR_after')].value
    g.close()
    
    # Opens LT3data    
    h = h5py.File(fp_LT3,'r')
    for k in h.keys():
        if type(h[k])==h5py.Group:
            ad3_reps = h[('/'+ str(k) + '/ssro/entanglement_events')].value
            ad3_ssro = h[('/'+ str(k) + '/ssro/RO_data')].value
            ad3_CR_before = h[('/'+ str(k) + '/ssro/CR_before')].value
            ad3_CR_after = h[('/'+ str(k) + '/ssro/CR_after')].value
    h.close()    
    

    sync_num_with_markers = sync_numbers[pq_tools.filter_marker(fp_BS,chan)]
    unique_sync_num_with_markers = np.unique(sync_num_with_markers)

    if VERBOSE:

        print 
        string = 'The number of events with PLU markers in run ' + str(i+1) + ' is:'
        print string, len(unique_sync_num_with_markers)
        print
        print 'Adwin LT1'
        print '---------'
        print 'Number of events:', ad1_reps,
        if ad1_reps != len(unique_sync_num_with_markers):
            print 'number of Adwin LT1 events does not match the PLU marker \
                    events - data set seems faulty :('
        else:
            print 'OK :)'
            
        print 
        print 'Adwin LT3'
        print '---------'
        print 'Number of events:', ad3_reps,
        if ad3_reps != len(unique_sync_num_with_markers):
            print 'number of Adwin LT1 events does not match the PLU marker \
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
        Photon 2 Channel | Attempts | Amout of Photons LT1 | Amount of Photons LT 3 | 
        CR Check Before LT1 | CR Check After LT1 | CR Check Before LT3 | 
        CR Check After LT3 | psiminus | absolute time
        """
        
        columns = "Sync_Number, Sync_Time_photon_1, Sync_Time_photon_2, Channel_photon_1,\
Channel_photon_2, Attempts, Amount_of_ph_LT1, CR_check_before_LT1,CR_check_after_LT1,\
Amount_of_ph_LT3, CR_check_before_LT3, CR_check_after_LT3, psiminus, abs_time"

        _a = {'Columns': columns}
                
        _event = np.array([s, 
                        stimes[0],
                        stimes[1], 
                        chans[0], 
                        chans[1], 
                        attempt, 
                        ad1_ssro[i],
                        ad3_ssro[i],
                        ad1_CR_before[i],
                        ad1_CR_after[i],
                        ad3_CR_before[i],
                        ad3_CR_after[i], 
                        psiminus, 
                        PLU_mrkr_abs_times[i]])
                        
        entanglement_events = np.vstack((entanglement_events, _event))

    if VERBOSE:
        print
        print 'Found {} valid entanglement events.'.format(int(len(entanglement_events)))
        print '===================================='
        print
    
    return entanglement_events, _a