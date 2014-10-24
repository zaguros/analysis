import numpy as np
import h5py
import Settings, files
from analysis.lib.tools import toolbox as tb
from analysis.lib.pq import pq_tools, pq_plots


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