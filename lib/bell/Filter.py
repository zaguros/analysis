import numpy as np
import h5py
import Settings, files


########### Determine which events have PLU Makers ###########

def sync_num_ent_events_HH(fp):
    f = h5py.File(fp, 'r')
    sync_numbers = f['/PQ_sync_number-1'].value
    f.close()

    is_ph_win1, is_ph_win2 = get_photons_in_sync_windows(fp)
    #print sync_numbers
    #print sum(is_ph_win1)
    sync_num_win1 = sync_numbers[is_ph_win1]
    #print sync_num_win1
    sync_num_win2 = sync_numbers[is_ph_win2]
    #print sync_num_win2
    
    _is_un_sync_num = np.in1d(sync_num_win1,sync_num_win2)
    unique_sync_ent_events = sync_num_win1[_is_un_sync_num]
    
    return unique_sync_ent_events        

def get_markers(fp, chan):
    """
    returns a filter (1d-array): whether events are markers on the given channel
    """
    f = h5py.File(fp, 'r')
    channel = f['/PQ_channel-1'].value
    special = f['/PQ_special-1'].value
    f.close()
    
    is_special = special==1
    is_channel = channel==chan
    
    return (is_special & is_channel)

def filter_on_same_sync_number(source_sync_numbers, target_sync_numbers):
    """
    returns a filter for target_sync_numbers that's true for all sync numbers
    that are also in source_sync_numbers.
    """
    return np.in1d(target_sync_numbers, source_sync_numbers)
  
def filter_marker(fp, chan):
    """
    Note: at the moment this filter includes the marker events
    on which we filter.
    """
    is_mrkr = get_markers(fp, chan)
        
    f = h5py.File(fp, 'r')
    sync_numbers = f['/PQ_sync_number-1'].value
    f.close()
    
    marker_sync_numbers = sync_numbers[is_mrkr]
    
    return filter_on_same_sync_number(marker_sync_numbers, sync_numbers)
    
def sync_num_with_markers(fp, chan, Force_eval = True):
    """
    Returns a list of the unique sync numbers with markers
    """

    if files.has_analysis_data(fp, 'Entanglement_events') and not Force_eval:
        tev, _a  = files.get_analysis_data(fp, 'Entanglement_events')
        unique_sync_num_with_markers = tev[:,0]
        return unique_sync_num_with_markers
    
    f = h5py.File(fp, 'r')
    sync_numbers = f['/PQ_sync_number-1'].value
    f.close()
  
    sync_num_with_markers = sync_numbers[filter_marker(fp,chan)]
    unique_sync_num_with_markers = np.unique(sync_num_with_markers)
    
    return unique_sync_num_with_markers

##### Returns filters for photon windows (in combination with PLU markers) #####

def get_photons_in_sync_windows(fp):
    """
    Returns two filters whether events are in the first or 
    in the second time window.
    """
    f = h5py.File(fp, 'r')
    sync_time = f['/PQ_sync_time-1'].value
    special = f['/PQ_special-1'].value
    f.close()
    

    is_photon = special == 0
    
    is_event_first_window = (sync_time > Settings.first_win_min) & \
                         (sync_time <= Settings.first_win_max)
    is_event_second_window = (sync_time > Settings.second_win_min) & \
                          (sync_time <= Settings.second_win_max)

    is_photon_first_window = is_photon & is_event_first_window
    is_photon_second_window = is_photon & is_event_second_window
    
    is_photon_check = is_photon_first_window | is_photon_second_window
        
    if sum(is_photon_check) != sum(is_photon):
        print "Not all detected photons are in the broad windows set"
        
    return is_photon_first_window, is_photon_second_window
   
         
def get_photons_with_markers(fp, chan):
    """
    Return two filters whether events are first window photons or second window
    photons with markers/
    """
    
    is_photon_first_window, is_photon_second_window = get_photons_in_sync_windows(fp)
    is_events_with_marker = filter_marker(fp,chan)
    
    is_photon_first_window_with_markers = is_photon_first_window & \
                                            is_events_with_marker
    is_photon_second_window_with_markers = is_photon_second_window &\
                                            is_events_with_marker

    return is_photon_first_window_with_markers, is_photon_second_window_with_markers
    
###### Filters  data on sync time #############################################

def get_events_in_correct_range(Total_entanglement_events, **kw):
   
    VERBOSE = kw.pop('VERBOSE',False)
        
    # Initialize values for bar plot
    psiminus_filt_up_up = 0
    psiminus_filt_up_down = 0
    psiminus_filt_down_up = 0
    psiminus_filt_down_down = 0

    psiplus_filt_up_up = 0
    psiplus_filt_up_down = 0
    psiplus_filt_down_up = 0
    psiplus_filt_down_down = 0
    
    # Makes boolean filters which determine if the SSRO correspond with the up 
    # (photons are emitted) or down (no photons are emitted) state
    is_up_LT1 = Total_entanglement_events[:,6] > 0
    is_down_LT1 = Total_entanglement_events[:,6] == 0
    is_up_LT3 = Total_entanglement_events[:,7] > 0
    is_down_LT3 = Total_entanglement_events[:,7] == 0
    
    # Makes boolean filters for up, up; up, down; down, up; down, down events
    is_upLT1_upLT3 = is_up_LT1 & is_up_LT3
    is_upLT1_downLT3 = is_up_LT1 & is_down_LT3
    is_downLT1_upLT3 = is_down_LT1 & is_up_LT3
    is_downLT1_downLT3 = is_down_LT1 & is_down_LT3

    # The final filter.
    for i in np.arange(len(Total_entanglement_events)):
        if Total_entanglement_events[i,3] == 0:
            Filt_1st_win_start = Settings.ch0_start
            Filt_1st_win_stop = Settings.ch0_stop
        else:
            Filt_1st_win_start = Settings.ch1_start
            Filt_1st_win_stop = Settings.ch1_stop
        if Total_entanglement_events[i,3] == 0:
            Filt_2nd_win_start = Settings.ch0_start + Settings.dif_win1_win2
            Filt_2nd_win_stop = Settings.ch0_stop + Settings.dif_win1_win2
        else:
            Filt_2nd_win_start = Settings.ch1_start + Settings.dif_win1_win2
            Filt_2nd_win_stop = Settings.ch1_stop + Settings.dif_win1_win2
        
        if Total_entanglement_events[i,12] == 1: # Then it is a psi minus event.
            # Now it is checked if the sync times of the first and second photon
            # lie within the window set above
            if ((Total_entanglement_events[i,1] >= Filt_1st_win_start) and \
            (Total_entanglement_events[i,1] <= Filt_1st_win_stop)) and \
            ((Total_entanglement_events[i,2] >= Filt_2nd_win_start) and \
            (Total_entanglement_events[i,2] <= Filt_2nd_win_stop)):
                # Now we check what kind of event we have
                if is_up_LT1[i] & is_up_LT3[i]:
                    psiminus_filt_up_up = psiminus_filt_up_up + 1
                elif is_up_LT1[i] & is_down_LT3[i]:
                    psiminus_filt_up_down = psiminus_filt_up_down + 1
                elif is_down_LT1[i] & is_up_LT3[i]:
                    psiminus_filt_down_up = psiminus_filt_down_up + 1
                elif is_down_LT1[i] & is_down_LT3[i]:
                    psiminus_filt_down_down = psiminus_filt_down_down + 1
        else: # These are the psi plus events. The filter is exactly the same
            if ((Total_entanglement_events[i,1] >= Filt_1st_win_start) and \
            (Total_entanglement_events[i,1] <= Filt_1st_win_stop)) and \
            ((Total_entanglement_events[i,2] >= Filt_2nd_win_start) and \
            (Total_entanglement_events[i,2] <= Filt_2nd_win_stop)):
                if is_up_LT1[i] & is_up_LT3[i]:
                    psiplus_filt_up_up = psiplus_filt_up_up + 1
                elif is_up_LT1[i] & is_down_LT3[i]:
                    psiplus_filt_up_down = psiplus_filt_up_down + 1
                elif is_down_LT1[i] & is_up_LT3[i]:
                    psiplus_filt_down_up = psiplus_filt_down_up + 1
                elif is_down_LT1[i] & is_down_LT3[i]:
                    psiplus_filt_down_down = psiplus_filt_down_down + 1
                
    psiplus_filt_bars = np.array([psiplus_filt_up_up, psiplus_filt_up_down, \
                                  psiplus_filt_down_up, psiplus_filt_down_down])
    psiplus_filt_bars_norm = psiplus_filt_bars/float(sum(psiplus_filt_bars))

    psiminus_filt_bars = np.array([psiminus_filt_up_up, psiminus_filt_up_down,\
                                psiminus_filt_down_up, psiminus_filt_down_down])
    psiminus_filt_bars_norm = psiminus_filt_bars/float(sum(psiminus_filt_bars))
    
    
    if VERBOSE:
        print 'There are {} psiplus entanglement events after \
filtering.'.format(sum(psiplus_filt_bars))
        print 'There are {} psiminus entanglement events after \
filtering.'.format(sum(psiminus_filt_bars))
        
    return psiplus_filt_bars, psiplus_filt_bars_norm, psiminus_filt_bars,\
            psiminus_filt_bars_norm
            
############# DT filters ########################################################

def DT_filter_max(Total_entanglement_events, dt):
      
    Sync_time_photon1 = Total_entanglement_events[:,1]
    Sync_time_photon2 = Total_entanglement_events[:,2]
    Delta_T = Sync_time_photon2 - Sync_time_photon1 - 600.
    
    for i in np.arange(len(Sync_time_photon1)):
        if Total_entanglement_events[i,3] == 0 and Total_entanglement_events[i,4] == 1:
            Delta_T[i] = Delta_T[i] - 1.
        if Total_entanglement_events[i,3] == 1 and Total_entanglement_events[i,4] == 0:
            Delta_T[i] = Delta_T[i] + 1.
            
    is_filtered_events = np.abs(Delta_T) <= dt

    
    return Total_entanglement_events[is_filtered_events]
    
def DT_filter_interval(Total_entanglement_events, dt_min, dt_max):
      
    Sync_time_photon1 = Total_entanglement_events[:,1]
    Sync_time_photon2 = Total_entanglement_events[:,2]
    Delta_T = Sync_time_photon2 - Sync_time_photon1 - 600.
    
    for i in np.arange(len(Sync_time_photon1)):
        if Total_entanglement_events[i,3] == 0 and Total_entanglement_events[i,4] == 1:
            Delta_T[i] = Delta_T[i] - 1.
        if Total_entanglement_events[i,3] == 1 and Total_entanglement_events[i,4] == 0:
            Delta_T[i] = Delta_T[i] + 1.
            
    is_filtered_events_dtmax = np.abs(Delta_T) <= dt_max 
    is_filtered_events_dtmin = np.abs(Delta_T) >= dt_min
    is_filtered_events = is_filtered_events_dtmax & is_filtered_events_dtmin

    
    return Total_entanglement_events[is_filtered_events]