import numpy as np
import h5py
import Settings, files, Filter
  

def get_entanglement_events(fp_BS,fp_LT1,fp_LT3,chan,unique_sync_num, force_eval=False):
    
    """
    Returns either the entanglement events already saved; the corresponding attributes
    and save = False, or returns the newly calculated entanglement events;
    the corresponding attributes (column names) and the save = True    
    """

   
    if files.has_analysis_data(fp_BS, 'Entanglement_events') and not force_eval:
        tev, _a = files.get_analysis_data(fp_BS, 'Entanglement_events')
        if Settings.VERBOSE:
            print
            print 'Found {} saved valid entanglement events.'.format(int(len(tev)))
            print '===================================='
            print
        return tev, _a , False

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
    
    if Settings.VERBOSE:

        print 
        print 'Adwin LT1'
        print '---------'
        print 'Number of events:', ad1_reps,
        if ad1_reps != len(unique_sync_num):
            print 'number of Adwin LT1 events does not match the PLU marker \
                    events - data set seems faulty :('
        else:
            print 'OK :)'
            
        print 
        print 'Adwin LT3'
        print '---------'
        print 'Number of events:', ad3_reps,
        if ad3_reps != len(unique_sync_num):
            print 'number of Adwin LT1 events does not match the PLU marker \
                    events - data set seems faulty :('
        else:
            print 'OK :)'
    
    # Gets filters for photons with markers in the first and second window
    # from the Filter file
    is_photon_1st_window_with_markers, is_photon_2nd_window_with_markers =\
                                    Filter.get_photons_with_markers(fp_BS,chan)
    
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
    for i,s in enumerate(unique_sync_num):
        
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
        
        _a = "Sync_Number, Sync_Time_photon_1, Sync_Time_photon_2, Channel_photon_1,\
Channel_photon_2, Attempts, Amount_of_ph_LT1, CR_check_before_LT1,CR_check_after_LT1,\
Amount_of_ph_LT3, CR_check_before_LT3, CR_check_after_LT3, psiminus, abs_time"
                
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

    if Settings.VERBOSE:
        print
        print 'Found {} valid entanglement events.'.format(int(len(entanglement_events)))
        print '===================================='
        print
    
    return entanglement_events, _a, True      