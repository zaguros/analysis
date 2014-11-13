import numpy as np
import h5py
import Settings, files, Filter
from analysis.lib.tools import toolbox as tb
from analysis.lib.lde import sscorr
from analysis.lib.pq import pq_tools, pq_plots

## Analysis of raw_data ##

def analyze_raw_data(Total_entanglement_events, Psiminus_event, Noof_ph_LT1, Noof_ph_LT2):
    """
    Returns boolean filters which determine if an event is psiminus or psiplus 
    (see Events.get_entanglement_events for definition of both). Put in a Total entanglement
    event table, the column number where information about Psiminus_events are saved, the column
    number where the number of photons in LT1 is defined and the column number where the number
    of photons in LT2 is defined.
    """
    is_psiminus = Total_entanglement_events[:,Psiminus_event] == 1
    is_psiplus = Total_entanglement_events[:,Psiminus_event] == 0
    
    if Settings.VERBOSE:
        print 'Number of psiminus events:', sum(is_psiminus)
        print 'Number of psiplus event:', sum(is_psiplus)

    # Makes boolean filters which determine if the SSRO correspond with the up (photons are emitted) or down (no photons are emitted) state
    is_up_LT1 = Total_entanglement_events[:,Noof_ph_LT1] > 0
    is_down_LT1 = Total_entanglement_events[:,Noof_ph_LT1] == 0
    is_up_LT3 = Total_entanglement_events[:,Noof_ph_LT2] > 0
    is_down_LT3 = Total_entanglement_events[:,Noof_ph_LT2] == 0

    # Makes boolean filters for up, up; up, down; down, up; down, down events
    is_upLT1_upLT3 = is_up_LT1 & is_up_LT3
    is_upLT1_downLT3 = is_up_LT1 & is_down_LT3
    is_downLT1_upLT3 = is_down_LT1 & is_up_LT3
    is_downLT1_downLT3 = is_down_LT1 & is_down_LT3

    # Makes boolean filters for psiplus state
    is_psiplus_up_up = is_upLT1_upLT3 & is_psiplus
    is_psiplus_up_down = is_upLT1_downLT3 & is_psiplus
    is_psiplus_down_up = is_downLT1_upLT3 & is_psiplus
    is_psiplus_down_down = is_downLT1_downLT3 & is_psiplus

    # Makes boolean filters for psiminus state
    is_psiminus_up_up = is_upLT1_upLT3 & is_psiminus
    is_psiminus_up_down = is_upLT1_downLT3 & is_psiminus
    is_psiminus_down_up = is_downLT1_upLT3 & is_psiminus
    is_psiminus_down_down = is_downLT1_downLT3 & is_psiminus

    # Determines the amount of up, up; up, down; down, up and down, down events for the psiplus state
    psiplus_up_up = sum(is_psiplus_up_up)
    psiplus_up_down = sum(is_psiplus_up_down)
    psiplus_down_up = sum(is_psiplus_down_up)
    psiplus_down_down = sum(is_psiplus_down_down)

    psiplus_bars = np.array([psiplus_up_up, psiplus_up_down, psiplus_down_up, psiplus_down_down])
    psiplus_bars_norm = psiplus_bars/float(sum(psiplus_bars))


    # Determines the amount of up, up; up, down; down, up and down, down events for the psiminus state
    psiminus_up_up = sum(is_psiminus_up_up)
    psiminus_up_down = sum(is_psiminus_up_down)
    psiminus_down_up = sum(is_psiminus_down_up)
    psiminus_down_down = sum(is_psiminus_down_down)

    psiminus_bars = np.array([psiminus_up_up, psiminus_up_down, psiminus_down_up, psiminus_down_down])
    psiminus_bars_norm = psiminus_bars/float(sum(psiminus_bars))
    
    return psiplus_bars, psiplus_bars_norm, psiminus_bars, psiminus_bars_norm
    
## Analysis of filtered data ##


def Ent_fid_vs_dt_max(Total_entanglement_events, pts, Start_dt, Last_dt, F1a, F1b, F0a, F0b):
    # Defines range for dt
    x=np.linspace(Start_dt,Last_dt,pts)

    # Initializes Psiplus and Psiminus fidelity and errors
    Psiplus_Fid = np.zeros(pts)
    Psiminus_Fid = np.zeros(pts)
    err_Psiplus_Fid = np.zeros(pts)
    err_Psiminus_Fid = np.zeros(pts)

    # Initializes Corrected Psiplus and Psiminus fidelities and errors
    Psiplus_Fid_corr = np.zeros(pts)
    Psiminus_Fid_corr = np.zeros(pts)
    err_Psiplus_Fid_corr = np.zeros(pts)
    err_Psiminus_Fid_corr = np.zeros(pts)

    # Initializes numbers of Psiplus and Psiminus events
    psiplus_events= np.zeros(pts)
    psiminus_events = np.zeros(pts)
    

    for i,a in enumerate(x):

        # Filter on photons with dt = a
        Total_entanglement_events_dt_filter = DT_filter_max(Total_entanglement_events, a)
    
        # Filters on photons with sync times within the tail
        psiplus_filt_bars, psiplus_filt_bars_norm, psiminus_filt_bars, psiminus_filt_bars_norm = \
                        get_events_in_correct_range(Total_entanglement_events_dt_filter, VERBOSE = False)
   
        # Determines total number of Psiplus and Psiminus events
        psiplus_events[i] = sum(psiplus_filt_bars)
        psiminus_events[i] = sum(psiminus_filt_bars)
    
        # Sets Readout Fidelities and calculates errors
        Psiplus_Fid[i] = (psiplus_filt_bars_norm[1]+psiplus_filt_bars_norm[2])
        Psiminus_Fid[i] = (psiminus_filt_bars_norm[0]+psiminus_filt_bars_norm[3])
        err_Psiplus_Fid[i] = Psiplus_Fid[i]/np.sqrt(psiplus_events[i])
        err_Psiminus_Fid[i] = Psiminus_Fid[i]/np.sqrt(psiminus_events[i])

   
        # Caclculates Corrected Readout Fidelities
        psiplus_filt_bars_norm_temp = psiplus_filt_bars_norm[::-1]
        psiminus_filt_bars_norm_temp = psiminus_filt_bars_norm[::-1]    
    
        psiplus_filt_bars_norm_corr_temp = sscorr.ssro_correct_twoqubit_state_photon_numbers(psiplus_filt_bars_norm_temp, F0a, F0b, F1a, F1b,
            return_error_bars=False, dF0a=0.01, dF0b=0.01, dF1a=0.01, dF1b=0.01, 
            verbose = False)
        psiminus_filt_bars_norm_corr_temp = sscorr.ssro_correct_twoqubit_state_photon_numbers(psiminus_filt_bars_norm_temp, F0a, F0b, F1a, F1b,
            return_error_bars=False, dF0a=0.01, dF0b=0.01, dF1a=0.01, dF1b=0.01, 
            verbose = False)
    
        psiplus_filt_bars_norm_corr = psiplus_filt_bars_norm_corr_temp[::-1]
        psiminus_filt_bars_norm_corr = psiminus_filt_bars_norm_corr_temp[::-1]
    
        # Sets corrected Readout Fidelities and errors
        Psiplus_Fid_corr[i] = (psiplus_filt_bars_norm_corr[1]+psiplus_filt_bars_norm_corr[2])
        Psiminus_Fid_corr[i] = (psiminus_filt_bars_norm_corr[0]+psiminus_filt_bars_norm_corr[3])
        err_Psiplus_Fid_corr[i] = Psiplus_Fid_corr[i]/np.sqrt(psiplus_events[i])
        err_Psiminus_Fid_corr[i] = Psiminus_Fid_corr[i]/np.sqrt(psiminus_events[i])
        
    return x, Psiplus_Fid, Psiminus_Fid, err_Psiplus_Fid, err_Psiminus_Fid, Psiplus_Fid_corr, Psiminus_Fid_corr, err_Psiplus_Fid_corr, err_Psiminus_Fid_corr, psiplus_events, psiminus_events
    
def Ent_fid_vs_dt_interval(pts, Interval_length, dt_start, Total_entanglement_events, F1a, F1b, F0a, F0b, Verbose = True):
    # Set x-axis
    x= np.zeros(pts)

    for i in np.arange(pts):
        x[i] = 0.5 + 3*i

    # Initializes Psiplus and Psiminus fidelity and errors
    Psiplus_Fid = np.zeros(pts)
    Psiminus_Fid = np.zeros(pts)
    err_Psiplus_Fid = np.zeros(pts)
    err_Psiminus_Fid = np.zeros(pts)

    # Initializes Corrected Psiplus and Psiminus fidelities and errors
    Psiplus_Fid_corr = np.zeros(pts)
    Psiminus_Fid_corr = np.zeros(pts)
    err_Psiplus_Fid_corr = np.zeros(pts)
    err_Psiminus_Fid_corr = np.zeros(pts)

    # Initializes numbers of Psiplus and Psiminus events
    psiplus_events= np.zeros(pts)
    psiminus_events = np.zeros(pts)


    for i in np.arange(pts):
   
        Start_dt = dt_start + i * Interval_length
        Stop_dt = Start_dt + Interval_length

        # Filter on photons with dt = a
        Total_entanglement_events_dt_filter = Filter.DT_filter_interval(Total_entanglement_events, Start_dt, Stop_dt)
    
        # Filters on photons with sync times within the tail
        psiplus_filt_bars, psiplus_filt_bars_norm, psiminus_filt_bars, psiminus_filt_bars_norm = \
                        get_events_in_correct_range(Total_entanglement_events_dt_filter, VERBOSE = False)
   
        # Determines total number of Psiplus and Psiminus events
        psiplus_events[i] = sum(psiplus_filt_bars)
        psiminus_events[i] = sum(psiminus_filt_bars)
        
        if Verbose:
            print  "Interval: ", i+1 
            print "Window is set correctly"
            if psiplus_events[i] == 0:
                print " Window is set too small no psiplus events found"
                Psiplus_Fid[i] = 0
            if psiminus_events[i] == 0:
                print " Window is set too small no psiminus events found"
                Psiminus_Fid[i] = 0
    
    
        # Sets Readout Fidelities and calculates errors
        Psiplus_Fid[i] = (psiplus_filt_bars_norm[1]+psiplus_filt_bars_norm[2])
        Psiminus_Fid[i] = (psiminus_filt_bars_norm[0]+psiminus_filt_bars_norm[3])
        err_Psiplus_Fid[i] = Psiplus_Fid[i]/np.sqrt(psiplus_events[i])
        err_Psiminus_Fid[i] = Psiminus_Fid[i]/np.sqrt(psiminus_events[i])

    
        # Caclculates Corrected Readout Fidelities
        psiplus_filt_bars_norm_temp = psiplus_filt_bars_norm[::-1]
        psiminus_filt_bars_norm_temp = psiminus_filt_bars_norm[::-1]    
    
        psiplus_filt_bars_norm_corr_temp = sscorr.ssro_correct_twoqubit_state_photon_numbers(psiplus_filt_bars_norm_temp, F0a, F0b, F1a, F1b,
            return_error_bars=False, dF0a=0.01, dF0b=0.01, dF1a=0.01, dF1b=0.01, 
            verbose = False)
        psiminus_filt_bars_norm_corr_temp = sscorr.ssro_correct_twoqubit_state_photon_numbers(psiminus_filt_bars_norm_temp, F0a, F0b, F1a, F1b,
            return_error_bars=False, dF0a=0.01, dF0b=0.01, dF1a=0.01, dF1b=0.01, 
            verbose = False)
    
        psiplus_filt_bars_norm_corr = psiplus_filt_bars_norm_corr_temp[::-1]
        psiminus_filt_bars_norm_corr = psiminus_filt_bars_norm_corr_temp[::-1]
    
        # Sets corrected Readout Fidelities and errors
        Psiplus_Fid_corr[i] = (psiplus_filt_bars_norm_corr[1]+psiplus_filt_bars_norm_corr[2])
        Psiminus_Fid_corr[i] = (psiminus_filt_bars_norm_corr[0]+psiminus_filt_bars_norm_corr[3])
        err_Psiplus_Fid_corr[i] = Psiplus_Fid_corr[i]/np.sqrt(psiplus_events[i])
        err_Psiminus_Fid_corr[i] = Psiminus_Fid_corr[i]/np.sqrt(psiminus_events[i])
        
    return x, Psiplus_Fid, Psiminus_Fid, err_Psiplus_Fid, err_Psiminus_Fid, Psiplus_Fid_corr, Psiminus_Fid_corr, err_Psiplus_Fid_corr, err_Psiminus_Fid_corr, psiplus_events, psiminus_events
  
## Tail Analysis ##
    
def ZPL_tail_analysis_per_run(BS_fp_len, fp, fp_LT1, ch0_start, ch1_start, ch0_stop, ch1_stop, dif_win1_win2,Verbose = True):

    # Initalizes file in which data is saved
    Tail_analysis = np.empty([5,BS_fp_len])

    for i in np.arange(BS_fp_len):
    
        # Opens file to retrieve BS sync times, numbers and filters for photons in channel 0 or 1
        f = h5py.File(fp[i], 'r')
        sync_numbers = f['/PQ_sync_number-1'].value
        sync_times = f['/PQ_sync_time-1'].value
        is_ph0, is_ph1 = pq_tools.get_photons(f)
        f.close()

        # Opens file to retrieve number of repetitions
        g = h5py.File(fp_LT1[i],'r')
        for k in g.keys():
            if type(g[k])==h5py.Group and not 'analysis' in k:
                Completed_reps = g[('/'+ str(k) + '/ssro/completed_reps')].value
        g.close()

        # Calculates the total repetitions
        Total_reps =  Completed_reps * 250
    
        # Makes a filter for events that are photons
        is_photon = is_ph0 | is_ph1
    
        # Makes a list of sync times for photons in channel 0 and 1
        sync_times_ph0 = sync_times[is_ph0]
        sync_times_ph1 = sync_times[is_ph1]
    
        """
        Note that all filters below can only be applied on sync_times_ph0 or sync_times_ph1 arrays with the same length
        """
    
        # Filters events in channel 0 on being in the first window or the second window, windows are set above
        is_event_first_window_ph0 = (sync_times_ph0 > ch0_start * 1e3) & \
                         (sync_times_ph0 <= ch0_stop * 1e3)
        is_event_second_window_ph0 = (sync_times_ph0 > ch0_start * 1e3 + dif_win1_win2 * 1e3) & \
                          (sync_times_ph0 <= ch0_stop * 1e3 + dif_win1_win2 * 1e3)
        # Filters events in channel 1 on being in the first window or the second window, windows are set above
        is_event_first_window_ph1 = (sync_times_ph1 > ch1_start * 1e3) & \
                         (sync_times_ph1 <= ch1_stop * 1e3)
        is_event_second_window_ph1 = (sync_times_ph1 > ch1_start *1e3 + dif_win1_win2 * 1e3) & \
                          (sync_times_ph1 <= ch1_stop * 1e3 + dif_win1_win2 * 1e3)
        
        # Calculates the amount of photons in the first or second window for both channel
        photons_first_window_ph0 = sum(is_event_first_window_ph0)
        photons_second_window_ph0 = sum(is_event_second_window_ph0)
        photons_first_window_ph1 = sum(is_event_first_window_ph1)
        photons_second_window_ph1 = sum(is_event_second_window_ph1)

        # Combines the photons from both channels
        photons_first_window = photons_first_window_ph0 + photons_first_window_ph1
        photons_second_window = photons_second_window_ph0 + photons_second_window_ph1
        Total_number_of_photons = sum(is_photon)

        # Normalize the amount of photons with respect to the total amount of repetitions
        photons_first_window_per_rep = photons_first_window / float(Total_reps)
        photons_second_window_per_rep = photons_second_window / float(Total_reps)
    
        # Saves calculated data
        Tail_analysis[:,i] = [Total_number_of_photons, photons_first_window, photons_second_window,\
                          photons_first_window_per_rep, photons_second_window_per_rep]

        if Verbose:
            print "Run: ", i+1
            print "Total number of photons: ", Total_number_of_photons
            print "Number of photons in first window", photons_first_window
            print "Number of photons in second window", photons_second_window
            print "Number of photons in first window per repetetion", photons_first_window_per_rep
            print "Number of photons in second window per repetetion", photons_second_window_per_rep
            print
    
    return Tail_analysis 
    
def PSB_tail_analysis_per_run(fp, first_window_start, first_window_stop, second_window_start, second_window_stop, Verbose = True):
    
    Tail_analysis = np.empty([5,len(fp)])
    
    for i in np.arange(len(fp)):
   
        if tb.has_data(fp[i], 'PQ_sync_time-1'):

            # Opens file to retrieve number of repetitions
            g = h5py.File(fp[i],'r')
            sync_numbers = g['/PQ_sync_number-1'].value
            sync_times = g['/PQ_sync_time-1'].value
            is_ph0, is_ph1 = pq_tools.get_photons(g) #Note that there are only photons in channel 0
            for k in g.keys():
                if type(g[k])==h5py.Group:
                    Completed_reps = g[('/'+ str(k) + '/ssro/completed_reps')].value
            g.close()

            # Calculates the total repetitions
            Total_reps =  Completed_reps * 250
    
            # Makes a list of sync times for photons 
            sync_times_ph0 = sync_times[is_ph0]

            """
            Note that all filters below can only be applied on sync_times_ph0 or sync_times_ph1 arrays with the same length
            """
    
            # Filters events in channel 0 on being in the first window or the second window, windows are set above
            is_event_first_window_ph0 = (sync_times_ph0 > first_window_start) & \
                             (sync_times_ph0 <= first_window_stop)
            is_event_second_window_ph0 = (sync_times_ph0 > second_window_start) & \
                              (sync_times_ph0 <= second_window_stop)


            # Calculates the amount of photons in the first or second window for both channel
            photons_first_window = sum(is_event_first_window_ph0)
            photons_second_window = sum(is_event_second_window_ph0)
            Total_number_of_photons = sum(is_ph0)

            # Normalize the amount of photons with respect to the total amount of repetitions
            photons_first_window_per_rep = photons_first_window / float(Total_reps)
            photons_second_window_per_rep = photons_second_window / float(Total_reps)
    
            # Saves calculated data
            Tail_analysis[:,i] = [Total_number_of_photons, photons_first_window, photons_second_window,\
                          photons_first_window_per_rep, photons_second_window_per_rep]

            if Settings.VERBOSE:
                print "Run: ", i+1
                print "Total number of photons: ", Total_number_of_photons
                print "Number of photons in first window", photons_first_window
                print "Number of photons in second window", photons_second_window
                print "Number of photons in first window per repetition", photons_first_window_per_rep
                print "Number of photons in second window per repetition", photons_second_window_per_rep
                print
            
        else:
            #Tail_analysis_LT1[:,i] = [0,0,0,0,0]
            if Settings.VERBOSE:
                print "Run: ", i+1
                print "There is no data for this run in LT1"
                Tail_analysis[:,i] = np.array([0,0,0,0,0])
                print
                
    return Tail_analysis
    
## Check if  HH entanglement events correspond with PLU Markers
   
def HH_event_PLU_Marker_Check(BS_fp_len, fp_BS, fp_LT1, fp_LT3, first_win_min, first_win_max, second_win_min, second_win_max):
    """
    This function checks the number of entanglement events which is obtained from the 
    Hydraharp data and compares it to the entanglement events found by markers.
    """


    for i in np.arange(BS_fp_len):
    
        # Checks if there is enough data for the analysis
        if tb.has_data(fp_LT1[i], 'PQ_sync_time-1') and tb.has_data(fp_LT3[i], 'PQ_sync_time-1'):
            print "Run: ", str(i+1)
            
            # Opens LT1 data  
            g = h5py.File(fp_LT1[i],'r')
            for k in g.keys():
                if type(g[k])==h5py.Group:
                    ad1_reps = g[('/'+ str(k) + '/ssro/entanglement_events')].value
            g.close()
    
            # Opens LT3 data    
            h = h5py.File(fp_LT3[i],'r')
            for k in h.keys():
                if type(h[k])==h5py.Group:
                    ad3_reps = h[('/'+ str(k) + '/ssro/entanglement_events')].value
            h.close()  
            
            # Open BS data
            pqf = h5py.File(fp_BS[i], 'r')
            sync_numbers = pqf['/PQ_sync_number-1'].value

            # Defines two filters for photons arriving in channel 0 and 1
            is_ph0, is_ph1 = pq_tools.get_photons(pqf)
            
            test = second_win_max - first_win_max

            # Defines filters for if photons in channel 0 and 1 are arriving in the first or scond window  
            print first_win_min
            print first_win_max
            print second_win_min
            print second_win_max       
            is_valid_ph0_w1 = is_ph0 & pq_tools.filter_synctimes(pqf,first_win_min ,first_win_max)
            is_valid_ph0_w2 = is_ph0 & pq_tools.filter_synctimes(pqf,second_win_min ,second_win_max)
            is_valid_ph1_w1 = is_ph1 & pq_tools.filter_synctimes(pqf,first_win_min ,first_win_max)
            is_valid_ph1_w2 = is_ph1 & pq_tools.filter_synctimes(pqf,second_win_min ,second_win_max)
            pqf.close()
    
    
            # Applies filters above to get a list of sync numbers
            valid_ph0_w1_sn = sync_numbers[is_valid_ph0_w1]
            print len(valid_ph0_w1_sn)
            print len(np.unique(valid_ph0_w1_sn))
            valid_ph0_w2_sn = sync_numbers[is_valid_ph0_w2]
            valid_ph1_w1_sn = sync_numbers[is_valid_ph1_w1]
            valid_ph1_w2_sn = sync_numbers[is_valid_ph1_w2]
            

            # Defines type of events by checking if sync numbers are in two lists
            is_00_sn = np.in1d(valid_ph0_w1_sn,valid_ph0_w2_sn)
            w0w0_sn = valid_ph0_w1_sn[is_00_sn]
            is_01_sn = np.in1d(valid_ph0_w1_sn,valid_ph1_w2_sn)
            w0w1_sn = valid_ph0_w1_sn[is_01_sn]
            is_10_sn = np.in1d(valid_ph1_w1_sn,valid_ph0_w2_sn)
            w1w0_sn = valid_ph1_w1_sn[is_10_sn]
            is_11_sn = np.in1d(valid_ph1_w1_sn,valid_ph1_w2_sn)
            w1w1_sn = valid_ph1_w1_sn[is_11_sn]
    
            # Defines a list with psiplus and psiminus events
            psi_plus_sn = np.unique(np.append(w0w0_sn,w1w1_sn))
            psi_minus_sn = np.unique(np.append(w0w1_sn,w1w0_sn))
    
            if Settings.VERBOSE == True:
                print "Number of Hydraharp psi_plus events: ", len(psi_plus_sn)
                print "Number of Hydraharp psi_minus events: ", len(psi_minus_sn)
                print "Number of PLU Markers LT1: ", ad1_reps
                print "Number of PLU Markers LT3: ", ad3_reps
                print
        else:
            print "There is no analysis data for run", i+1

###### Filters  data on sync time #############################################

def get_events_in_correct_range(Total_entanglement_events, Noof_ph_LT1, Noof_ph_LT2, chan_ph_win_1, chan_ph_win_2, sync_time_ph_win_1, sync_time_ph_win_2, psiminus_event, ch0_start, ch0_stop, ch1_start, ch1_stop, dif_win1_win2, **kw):
   
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
    is_up_LT1 = Total_entanglement_events[:,Noof_ph_LT1] > 0
    is_down_LT1 = Total_entanglement_events[:,Noof_ph_LT1] == 0
    is_up_LT3 = Total_entanglement_events[:,Noof_ph_LT2] > 0
    is_down_LT3 = Total_entanglement_events[:,Noof_ph_LT2] == 0
    
    # Makes boolean filters for up, up; up, down; down, up; down, down events
    is_upLT1_upLT3 = is_up_LT1 & is_up_LT3
    is_upLT1_downLT3 = is_up_LT1 & is_down_LT3
    is_downLT1_upLT3 = is_down_LT1 & is_up_LT3
    is_downLT1_downLT3 = is_down_LT1 & is_down_LT3

    # The final filter.
    for i in np.arange(len(Total_entanglement_events)):
        if Total_entanglement_events[i,chan_ph_win_1] == 0:
            Filt_1st_win_start = ch0_start
            Filt_1st_win_stop = ch0_stop
        else:
            Filt_1st_win_start = ch1_start
            Filt_1st_win_stop = ch1_stop
        if Total_entanglement_events[i,chan_ph_win_1] == 0:
            Filt_2nd_win_start = ch0_start + dif_win1_win2
            Filt_2nd_win_stop = ch0_stop + dif_win1_win2
        else:
            Filt_2nd_win_start = ch1_start + dif_win1_win2
            Filt_2nd_win_stop = ch1_stop + dif_win1_win2
        
        if Total_entanglement_events[i, psiminus_event ] == 1: # Then it is a psi minus event.
            # Now it is checked if the sync times of the first and second photon
            # lie within the window set above
            if ((Total_entanglement_events[i,sync_time_ph_win_1] >= Filt_1st_win_start) and \
            (Total_entanglement_events[i,sync_time_ph_win_1] <= Filt_1st_win_stop)) and \
            ((Total_entanglement_events[i,sync_time_ph_win_2] >= Filt_2nd_win_start) and \
            (Total_entanglement_events[i,sync_time_ph_win_2] <= Filt_2nd_win_stop)):
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
            if ((Total_entanglement_events[i,sync_time_ph_win_1] >= Filt_1st_win_start) and \
            (Total_entanglement_events[i,sync_time_ph_win_1] <= Filt_1st_win_stop)) and \
            ((Total_entanglement_events[i,sync_time_ph_win_2] >= Filt_2nd_win_start) and \
            (Total_entanglement_events[i,sync_time_ph_win_2] <= Filt_2nd_win_stop)):
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

def DT_filter_max(Total_entanglement_events, dt, sync_time_ph_win_1, sync_time_ph_win_2, chan_ph_win_1, chan_ph_win_2):
      
    Sync_time_photon1 = Total_entanglement_events[:,sync_time_ph_win_1]
    Sync_time_photon2 = Total_entanglement_events[:,sync_time_ph_win_2]
    Delta_T = Sync_time_photon2 - Sync_time_photon1 - 600.
    
    for i in np.arange(len(Sync_time_photon1)):
        if Total_entanglement_events[i,chan_ph_win_1] == 0 and Total_entanglement_events[i,chan_ph_win_2] == 1:
            Delta_T[i] = Delta_T[i] - 1.
        if Total_entanglement_events[i,chan_ph_win_1] == 1 and Total_entanglement_events[i,chan_ph_win_2] == 0:
            Delta_T[i] = Delta_T[i] + 1.
            
    is_filtered_events = np.abs(Delta_T) <= dt

    
    return Total_entanglement_events[is_filtered_events]

def Laser_filter(Total_Bell_events, first_win_min_ch0, window_length, dif_win1_win2, dif_ch0_ch1, \
                                     sync_time_ph_win_1, sync_time_ph_win_2, channel_ph_win_1, channel_ph_win_2):

    ch0_start = first_win_min_ch0
    ch0_stop = ch0_start + window_length
    ch0_start_w2 = ch0_start + dif_win1_win2
    ch0_stop_w2 = ch0_stop + dif_win1_win2

    sync_time_ph_1 = Total_Bell_events[:,sync_time_ph_win_1]
    channel_ph_1 = Total_Bell_events[:,channel_ph_win_1]

    sync_time_ph_2 = Total_Bell_events[:,sync_time_ph_win_2]
    channel_ph_2 = Total_Bell_events[:,channel_ph_win_2]

    first_window_ch0_filt = (sync_time_ph_1 > ch0_start) & (sync_time_ph_1 <= ch0_stop)
    second_window_ch0_filt = (sync_time_ph_2 > ch0_start_w2) & (sync_time_ph_2 <= ch0_stop_w2)


    first_window_ch1_filt = (sync_time_ph_1 > (ch0_start + dif_ch0_ch1)) & (sync_time_ph_1 <= (ch0_stop + dif_ch0_ch1))
    second_window_ch1_filt = (sync_time_ph_2 > (ch0_start_w2 + dif_ch0_ch1)) & (sync_time_ph_2 <= (ch0_stop_w2 + dif_ch0_ch1))

    is_photon_first_tail_ch0 = channel_ph_1 == 0
    is_photon_first_tail_ch1 = channel_ph_1 == 1

    is_photon_second_tail_ch0 = channel_ph_2 == 0
    is_photon_second_tail_ch1 = channel_ph_2 == 1

    is_ph_first_tail = (is_photon_first_tail_ch0 & first_window_ch0_filt) | \
                                                                (is_photon_first_tail_ch1 & first_window_ch1_filt)


    is_ph_second_tail = (is_photon_second_tail_ch0 & second_window_ch0_filt) | \
                                                            (is_photon_second_tail_ch1 & second_window_ch1_filt)

    is_tail_photons = is_ph_first_tail & is_ph_second_tail


    Filtered_Bell_events = Total_Bell_events[is_tail_photons]

    return Filtered_Bell_events

############## Bell ################################################################

def CHSH_value(Total_Bell_events, psiminus, rnd_num_LT3, rnd_num_LT4, num_phot_LT3, num_phot_LT4, VERBOSE = True):

    """
    The up state is defined as ms =0 (emitting photons) the down state is defined as ms = -1 no photons
    """
    # Psiminus Initialization
    N_psi_min_a_b = 0

    N_psi_min_a_b_uu = 0
    N_psi_min_a_b_ud = 0
    N_psi_min_a_b_du = 0
    N_psi_min_a_b_dd = 0

    N_psi_min_a_pr_b = 0

    N_psi_min_a_pr_b_uu = 0
    N_psi_min_a_pr_b_ud = 0
    N_psi_min_a_pr_b_du = 0
    N_psi_min_a_pr_b_dd = 0

    N_psi_min_a_b_pr = 0

    N_psi_min_a_b_pr_uu = 0
    N_psi_min_a_b_pr_ud = 0
    N_psi_min_a_b_pr_du = 0
    N_psi_min_a_b_pr_dd = 0

    N_psi_min_a_pr_b_pr = 0

    N_psi_min_a_pr_b_pr_uu = 0
    N_psi_min_a_pr_b_pr_ud = 0
    N_psi_min_a_pr_b_pr_du = 0
    N_psi_min_a_pr_b_pr_dd = 0

    # Psi plus initialization
    N_psi_plus_a_b = 0

    N_psi_plus_a_b_uu = 0
    N_psi_plus_a_b_ud = 0
    N_psi_plus_a_b_du = 0
    N_psi_plus_a_b_dd = 0

    N_psi_plus_a_pr_b = 0

    N_psi_plus_a_pr_b_uu = 0
    N_psi_plus_a_pr_b_ud = 0
    N_psi_plus_a_pr_b_du = 0
    N_psi_plus_a_pr_b_dd = 0

    N_psi_plus_a_b_pr = 0

    N_psi_plus_a_b_pr_uu = 0
    N_psi_plus_a_b_pr_ud = 0
    N_psi_plus_a_b_pr_du = 0
    N_psi_plus_a_b_pr_dd = 0

    N_psi_plus_a_pr_b_pr = 0

    N_psi_plus_a_pr_b_pr_uu = 0
    N_psi_plus_a_pr_b_pr_ud = 0
    N_psi_plus_a_pr_b_pr_du = 0
    N_psi_plus_a_pr_b_pr_dd = 0


    for i in range(len(Total_Bell_events)):

        if Total_Bell_events[i,psiminus] ==  1:
        # Now we have a psiminus event
            if (Total_Bell_events[i,rnd_num_LT3] == 0) & (Total_Bell_events[i,rnd_num_LT4] == 0):
                # Now we have a psiminus a,b event
                N_psi_min_a_b = N_psi_min_a_b + 1

                if (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_min_a_b_uu = N_psi_min_a_b_uu + 1 

                elif (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_min_a_b_du = N_psi_min_a_b_du + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_min_a_b_ud = N_psi_min_a_b_ud + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_min_a_b_dd = N_psi_min_a_b_dd + 1 

            elif (Total_Bell_events[i,rnd_num_LT3] == 1) & (Total_Bell_events[i,rnd_num_LT4] == 0):
                # Now we have a psiminus a', b event
                N_psi_min_a_pr_b = N_psi_min_a_pr_b + 1

                if (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_min_a_pr_b_uu = N_psi_min_a_pr_b_uu + 1 

                elif (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_min_a_pr_b_du = N_psi_min_a_pr_b_du + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_min_a_pr_b_ud = N_psi_min_a_pr_b_ud + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_min_a_pr_b_dd = N_psi_min_a_pr_b_dd + 1 

            elif (Total_Bell_events[i,rnd_num_LT3] == 0) & (Total_Bell_events[i,rnd_num_LT4] == 1):
                # Now we have a psiminus a,b' event
                N_psi_min_a_b_pr = N_psi_min_a_b_pr + 1

                if (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_min_a_b_pr_uu = N_psi_min_a_b_pr_uu + 1 

                elif (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_min_a_b_pr_du = N_psi_min_a_b_pr_du + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_min_a_b_pr_ud = N_psi_min_a_b_pr_ud + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_min_a_b_pr_dd = N_psi_min_a_b_pr_dd + 1

            elif (Total_Bell_events[i,rnd_num_LT3] == 1) & (Total_Bell_events[i,rnd_num_LT4] == 1):
                # Not we have a psiminus a',b' event
                N_psi_min_a_pr_b_pr = N_psi_min_a_pr_b_pr + 1

                if (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_min_a_pr_b_pr_uu = N_psi_min_a_pr_b_pr_uu + 1 

                elif (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_min_a_pr_b_pr_du = N_psi_min_a_pr_b_pr_du + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_min_a_pr_b_pr_ud = N_psi_min_a_pr_b_pr_ud + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_min_a_pr_b_pr_dd = N_psi_min_a_pr_b_pr_dd + 1


        elif Total_Bell_events[i,psiminus] ==  0:
            # Now we have a psiplus event

            if (Total_Bell_events[i,rnd_num_LT3] == 0) & (Total_Bell_events[i,rnd_num_LT4] == 0):
                # Now we have a psiplus a,b event
                N_psi_plus_a_b = N_psi_plus_a_b + 1

                if (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_plus_a_b_uu = N_psi_plus_a_b_uu + 1 

                elif (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_plus_a_b_du = N_psi_plus_a_b_du + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_plus_a_b_ud = N_psi_plus_a_b_ud + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_plus_a_b_dd = N_psi_plus_a_b_dd + 1 

            elif (Total_Bell_events[i,rnd_num_LT3] == 1) & (Total_Bell_events[i,rnd_num_LT4] == 0):
                # Now we have a psiplus a', b event
                N_psi_plus_a_pr_b = N_psi_plus_a_pr_b + 1

                if (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_plus_a_pr_b_uu = N_psi_plus_a_pr_b_uu + 1 

                elif (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_plus_a_pr_b_du = N_psi_plus_a_pr_b_du + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_plus_a_pr_b_ud = N_psi_plus_a_pr_b_ud + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_plus_a_pr_b_dd = N_psi_plus_a_pr_b_dd + 1 

            elif (Total_Bell_events[i,rnd_num_LT3] == 0) & (Total_Bell_events[i,rnd_num_LT4] == 1):
                # Now we have a psiplus a,b' event
                N_psi_plus_a_b_pr = N_psi_plus_a_b_pr + 1

                if (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_plus_a_b_pr_uu = N_psi_plus_a_b_pr_uu + 1 

                elif (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_plus_a_b_pr_du = N_psi_plus_a_b_pr_du + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_plus_a_b_pr_ud = N_psi_plus_a_b_pr_ud + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_plus_a_b_pr_dd = N_psi_plus_a_b_pr_dd + 1

            elif (Total_Bell_events[i,rnd_num_LT3] == 1) & (Total_Bell_events[i,rnd_num_LT4] == 1):
                # Not we have a psiplus a',b' event
                N_psi_plus_a_pr_b_pr = N_psi_plus_a_pr_b_pr + 1

                if (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_plus_a_pr_b_pr_uu = N_psi_plus_a_pr_b_pr_uu + 1 

                elif (Total_Bell_events[i,num_phot_LT3] > 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_plus_a_pr_b_pr_du = N_psi_plus_a_pr_b_pr_du + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] > 0):
                    N_psi_plus_a_pr_b_pr_ud = N_psi_plus_a_pr_b_pr_ud + 1 

                elif (Total_Bell_events[i,num_phot_LT3] == 0) & (Total_Bell_events[i,num_phot_LT4] == 0):
                    N_psi_plus_a_pr_b_pr_dd = N_psi_plus_a_pr_b_pr_dd + 1

    if N_psi_min_a_b != (N_psi_min_a_b_uu + N_psi_min_a_b_dd + N_psi_min_a_b_du + N_psi_min_a_b_ud):
        print "Something weird is hapening with psi minus a,b"

    if N_psi_min_a_pr_b != (N_psi_min_a_pr_b_uu + N_psi_min_a_pr_b_dd + N_psi_min_a_pr_b_du + N_psi_min_a_pr_b_ud):
        print "Something weird is hapening with psi minus a',b"

    if N_psi_min_a_b_pr != (N_psi_min_a_b_pr_uu + N_psi_min_a_b_pr_dd + N_psi_min_a_b_pr_du + N_psi_min_a_b_pr_ud):
        print "Something weird is hapening with psi minus a,b'"

    if N_psi_min_a_pr_b_pr != (N_psi_min_a_pr_b_pr_uu + N_psi_min_a_pr_b_pr_dd + 
                                                N_psi_min_a_pr_b_pr_du + N_psi_min_a_pr_b_pr_ud):
        print "Something weird is hapening with psi minus a',b'"

    if N_psi_plus_a_b != (N_psi_plus_a_b_uu + N_psi_plus_a_b_dd + N_psi_plus_a_b_du + N_psi_plus_a_b_ud):
        print "Something weird is hapening with psi plus a,b"

    if N_psi_plus_a_pr_b != (N_psi_plus_a_pr_b_uu + N_psi_plus_a_pr_b_dd + N_psi_plus_a_pr_b_du + N_psi_plus_a_pr_b_ud):
        print "Something weird is hapening with psi plus a',b"

    if N_psi_plus_a_b_pr != (N_psi_plus_a_b_pr_uu + N_psi_plus_a_b_pr_dd + N_psi_plus_a_b_pr_du + N_psi_plus_a_b_pr_ud):
        print "Something weird is hapening with psi plus a,b'"

    if N_psi_plus_a_pr_b_pr != (N_psi_plus_a_pr_b_pr_uu + N_psi_plus_a_pr_b_pr_dd + 
                                                N_psi_plus_a_pr_b_pr_du + N_psi_plus_a_pr_b_pr_ud):
        print "Something weird is hapening with psi plus a',b'"

    if VERBOSE:
        print "Number of events used for psimin a,b:", N_psi_min_a_b
        print "Number of events used for psimin a',b':",  N_psi_min_a_pr_b_pr
        print "Number of events used for psimin a',b:", N_psi_min_a_pr_b
        print "Number of events used for psimin a,b':", N_psi_min_a_b_pr
        print
        print 
        print "Number of events used for psiplus a,b:", N_psi_plus_a_b
        print "Number of events used for psiplus a',b':",  N_psi_plus_a_pr_b_pr
        print "Number of events used for psiplus a',b:", N_psi_plus_a_pr_b
        print "Number of events used for psiplus a,b':", N_psi_plus_a_b_pr

    E_a_b_min = (N_psi_min_a_b_uu - N_psi_min_a_b_du - N_psi_min_a_b_ud + N_psi_min_a_b_dd)/float(N_psi_min_a_b)
    E_a_pr_b_min = (N_psi_min_a_pr_b_uu - N_psi_min_a_pr_b_du - 
                                                N_psi_min_a_pr_b_ud + N_psi_min_a_pr_b_dd)/float(N_psi_min_a_pr_b)
    E_a_b_pr_min = (N_psi_min_a_b_pr_uu - N_psi_min_a_b_pr_du - 
                                                N_psi_min_a_b_pr_ud + N_psi_min_a_b_pr_dd)/float(N_psi_min_a_b_pr)
    E_a_pr_b_pr_min = (N_psi_min_a_pr_b_pr_uu - N_psi_min_a_pr_b_pr_du - 
                                            N_psi_min_a_pr_b_pr_ud + N_psi_min_a_pr_b_pr_dd)/float(N_psi_min_a_pr_b_pr)

    E_a_b_plus = (N_psi_plus_a_b_uu - N_psi_plus_a_b_du - N_psi_plus_a_b_ud + N_psi_plus_a_b_dd)/float(N_psi_plus_a_b)
    E_a_pr_b_plus = (N_psi_plus_a_pr_b_uu - N_psi_plus_a_pr_b_du - 
                                                N_psi_plus_a_pr_b_ud + N_psi_plus_a_pr_b_dd)/float(N_psi_plus_a_pr_b)
    E_a_b_pr_plus = (N_psi_plus_a_b_pr_uu - N_psi_plus_a_b_pr_du - 
                                                N_psi_plus_a_b_pr_ud + N_psi_plus_a_b_pr_dd)/float(N_psi_plus_a_b_pr)
    E_a_pr_b_pr_plus = (N_psi_plus_a_pr_b_pr_uu - N_psi_plus_a_pr_b_pr_du - 
                                        N_psi_plus_a_pr_b_pr_ud + N_psi_plus_a_pr_b_pr_dd)/float(N_psi_plus_a_pr_b_pr)

    CHSH_psimin = E_a_b_min + E_a_pr_b_min - E_a_b_pr_min + E_a_pr_b_pr_min
    CHSH_psiplus = -E_a_b_plus + E_a_pr_b_plus + E_a_b_pr_plus + E_a_pr_b_pr_plus

    return CHSH_psimin, CHSH_psiplus, E_a_b_min, E_a_pr_b_min, E_a_b_pr_min, E_a_pr_b_pr_min, E_a_b_plus, E_a_pr_b_plus, E_a_b_pr_plus, E_a_pr_b_pr_plus

