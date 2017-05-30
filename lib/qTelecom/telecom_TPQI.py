import numpy as np
import h5py, os
from analysis.lib.pq import pq_tools;reload(pq_tools)
from analysis.lib.m2.ssro import pqsequence; reload(pqsequence)
from analysis.lib.tools import toolbox as tb
from analysis.lib.purification import purify_pq; reload(purify_pq)
#import Settings, Filter
import matplotlib.pyplot as plt

def get_coincidences_and_adwin_data_from_folder(folder_primary,use_invalid_data_marker = False, folder_secondary = None,syncs_per_CR_check = 1, index = 1,force_coincidence_evaluation = False,contains = '', save = True,Verbose=False):

    sync_num_name = 'PQ_sync_number-' + str(index)
    # print 'this is the save!', save
    filepaths_primary = tb.get_all_msmt_filepaths(folder_primary, pattern = contains)

    if folder_secondary != None:
        load_secondary = True
        filepaths_secondary = tb.get_all_msmt_filepaths(folder_secondary, pattern = contains) 
    else:
        load_secondary = False

    first_data = True 

    # print filepaths
    for i,(fp,fp2) in enumerate(zip(filepaths_primary,filepaths_secondary)):
        if Verbose:
            print fp
            print fp2
        if abs(int(os.path.split(fp)[1][:6]) - int(os.path.split(fp2)[1][:6])) > 100:
            print fp
            print fp2
            raise Exception('Time difference too great!')
           
        
        pqf = pq_tools.pqf_from_fp(fp)
        if load_secondary:
            fs = os.path.split(fp2)[0]
            purifyPQ_s = purify_pq.purifyPQAnalysis(fs)

        if sync_num_name in pqf.keys():

            coinc = pq_tools.get_coincidences(pqf,force_coincidence_evaluation =force_coincidence_evaluation,save = save)

            if coinc.size:

                new_vars = np.array([])

                if use_invalid_data_marker:
                    # Fixed number of LDE attempts per repetition, so straightforward to pull CR check value index
                    ind = (np.floor(coinc[:,3]/syncs_per_CR_check)).astype(int)

                    # print cr_check_ind
                    invalid_data_marker = (purifyPQ_s.agrp['invalid_data_markers'].value)[ind]
                    new_vars = np.hstack((new_vars, invalid_data_marker[:, np.newaxis])) if new_vars.size else invalid_data_marker[:, np.newaxis]

                if first_data:
            
                    co = np.hstack((coinc,new_vars)) if new_vars.size else coinc
                    first_data = False

                else:
                       
                    co = np.vstack((co, np.hstack((coinc,new_vars)))) if new_vars.size else np.vstack((co, coinc))

        pqf.close()
        if load_secondary:
            purifyPQ_s.finish()

    
    return co



def _aggregated_coincidences(Base_Folder,**kw):
    
    force_coincidence_evaluation = kw.pop('force_coincidence_evaluation',False)
    save = kw.pop('save', True)
    use_invalid_data_marker = kw.pop('use_invalid_data_marker',False)
    folder_secondary = kw.pop('folder_secondary', None)
    syncs_per_CR_check = kw.pop('syncs_per_CR_check', 1)
    Verbose = kw.pop('Verbose', False)
    contains = kw.pop('contains', '')
    older_than = kw.pop('older_than', None)
    newer_than = kw.pop('newer_than', None)

    if use_invalid_data_marker and folder_secondary == None:
        raise Exception('No secondary folder passed')

    if use_invalid_data_marker:
        in_coincidences = get_coincidences_and_adwin_data_from_folder(Base_Folder,folder_secondary = folder_secondary,syncs_per_CR_check = syncs_per_CR_check,
            use_invalid_data_marker = True,force_coincidence_evaluation =force_coincidence_evaluation,save = save,Verbose = Verbose, contains = contains)
    else:
        in_coincidences = pq_tools.get_coincidences_from_folder(Base_Folder,force_coincidence_evaluation =force_coincidence_evaluation,save = save, contains  = contains, older_than = older_than, newer_than = newer_than,**kw)

    return in_coincidences

def filter_coincidences(coincidences, ch0_start, ch1_start, WINDOW_LENGTH, dif_win1_win2, noof_pulses,
                         column_st_0, column_st_1):

    tail_length = WINDOW_LENGTH
    pulse_sep = dif_win1_win2

    
    f_st0 = pq_tools.filter_synctimes(coincidences[:,column_st_0], ch0_start, ch0_start + tail_length, noof_pulses, pulse_sep, pq_file = False)
    f_st1 = pq_tools.filter_synctimes(coincidences[:,column_st_1], ch1_start, ch1_start +  tail_length, noof_pulses, pulse_sep, pq_file = False)
        
    return f_st0 & f_st1 # & f_dt

def filter_cr_check(load_cr_check,columns_cr, coincidences, cr_check_min_p, cr_check_min_s, filter_before_or_after, load_secondary):
    
    column_cr_before_p, column_cr_after_p, column_cr_before_s, column_cr_after_s = columns_cr

    if load_cr_check == False:
        return np.ones(coincidences.shape[0], dtype=bool)
    
    filter_col = coincidences[:,(column_cr_before_p if filter_before_or_after == 'before' else column_cr_after_p)] >= cr_check_min_p

    if load_secondary:
        filter_col_s = coincidences[:,(column_cr_before_s if filter_before_or_after == 'before' else column_cr_after_s)] >= cr_check_min_s
        filter_col = filter_col & filter_col_s

    return filter_col


def filter_syncnum(coincidences, syncs_per_CR_check, min_filter_attempts = 1, delta_attempts = 50, column_sn = 3):
    """
    Return a filter for syncs based on the number of attempts since a CR check
    """

    sync_num = coincidences[:,column_sn]

    fltr = ((((sync_num-1) % syncs_per_CR_check) + 1) >= min_filter_attempts ) & ((((sync_num-1) % syncs_per_CR_check) + 1) < (min_filter_attempts + delta_attempts) )

    return fltr
 


def filter_adwin_comm_time(load_adwin_comm_time,columns_adwin_com, coincidences, adwin_com_max_p, adwin_com_max_s, load_secondary):
    
    column_adwin_com_p, column_adwin_com_s = columns_adwin_com

    if load_adwin_comm_time == False:
        return np.ones(coincidences.shape[0], dtype=bool)
    
    filter_col = np.logical_and(coincidences[:,column_adwin_com_p] <= adwin_com_max_p , coincidences[:,column_adwin_com_p] > 0)

    if load_secondary:
        filter_col_s = np.logical_and(coincidences[:,column_adwin_com_s] <= adwin_com_max_s , coincidences[:,column_adwin_com_s] > 0)
        filter_col = filter_col & filter_col_s

    return filter_col

def filter_no_of_attempts(load_TPQI_attemtps,column_no_of_sequences,coincidences,min_attempts,max_attempts):
    if load_TPQI_attemtps == False:
        return np.ones(coincidences.shape[0], dtype=bool)

    return (coincidences[:,column_no_of_sequences] >= min_attempts) & (coincidences[:,column_no_of_sequences] <= max_attempts)

def TPQI_analysis(Base_Folder_primary, ch0_start, ch1_start, WINDOW_LENGTH, dif_win1_win2, noof_pulses, 
                                    filter_attempts = False, syncs_per_CR_check = 50, delta_attempts = 50, min_filter_attempts = 1, max_dt = -1,
                                    return_sn = False , 
                                    contains = 'TPQI',
                                    folder_is_daydir = False,
                                    older_than = None,
                                    newer_than = None,
                                    force_coincidence_evaluation = False,
                                    Verbose = True):
    # Gets coincident photons from Hydraharp data
    coincidences = _aggregated_coincidences(Base_Folder_primary,contains=contains,force_coincidence_evaluation = force_coincidence_evaluation,older_than = older_than, newer_than = newer_than,folder_is_daydir = folder_is_daydir)

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

    if filter_attempts:
        filtered_attempts = filter_syncnum(coincidences, syncs_per_CR_check, min_filter_attempts = min_filter_attempts, delta_attempts = delta_attempts, column_sn = column_sync_num_ch0)
        is_sync_time_filter = is_sync_time_filter & filtered_attempts

    filtered_dts = dts[is_sync_time_filter]

    if max_dt != -1:
        dt_filter = np.abs(np.mod(filtered_dts+0.5*dif_win1_win2/1e3,dif_win1_win2/1e3)-0.5*dif_win1_win2/1e3) < max_dt/1e3
        filtered_dts =filtered_dts[dt_filter]

    if Verbose:
        print
        print 'Found {} coincident photons after filtering.'.format(int(len(filtered_dts)))
        print '===================================='
        print

    return dts, filtered_dts

def get_photons_from_folder(folder,ch0_start,ch1_start, WINDOW_LENGTH,dif_win1_win2, noof_pulses,syncs_per_CR_check, min_filter_attempts,delta_attempts,index = 1, contains = '', pq_device = ''):
    sync_num_name = pq_device + 'PQ_sync_number-' + str(index)
    # print 'this is the save!', save
    filepaths = tb.get_all_msmt_filepaths(folder) 

    if contains != '':
        new_fps = []
        for f in filepaths:
            if contains in f:
                new_fps.append(f)
        filepaths = new_fps

    first_run = True
    ph0, ph1 = 0,0

    # print filepaths
    for i,f in enumerate(filepaths):
        
        pqf = pq_tools.pqf_from_fp(f, rights = 'r+')

        if sync_num_name in pqf.keys():
            ph0t,ph1t = get_clicks(pqf) 
            
            filt_phot = filter_photons(ph0t, ch0_start, WINDOW_LENGTH,dif_win1_win2, noof_pulses, 0)
            filt_attempt = filter_syncnum(ph0t, syncs_per_CR_check, min_filter_attempts = min_filter_attempts, delta_attempts = delta_attempts, column_sn = 1)
            ph0 += np.sum(ph0t[filt_phot & filt_attempt])
            
            filt_phot = filter_photons(ph1t, ch1_start, WINDOW_LENGTH,dif_win1_win2, noof_pulses, 0)
            filt_attempt = filter_syncnum(ph1t, syncs_per_CR_check, min_filter_attempts = min_filter_attempts, delta_attempts = delta_attempts, column_sn = 1)
            ph1 += np.sum(ph1t[filt_phot & filt_attempt])

        pqf.close()
          
    return ph0,ph1


def get_clicks(pqf, index = 1, pq_device = ''):


    sync_time_name = pq_device + '/PQ_sync_time-' + str(index)
    tot_time_name =  pq_device + '/PQ_time-' + str(index)
    sync_num_name = pq_device + '/PQ_sync_number-' + str(index)

    sync_time = pqf[sync_time_name].value
    total_time = pqf[tot_time_name].value
    sync_number = pqf[sync_num_name].value

    is_ph0, is_ph1 = pq_tools.get_photons(pqf)
    # thin down a bit with loose filtering

    st0 = sync_time[is_ph0]
    sn0 = sync_number[is_ph0]
    
    st1 = sync_time[is_ph1]
    sn1 = sync_number[is_ph1]

    return np.array([st0,sn0]).T,np.array([st1,sn1]).T

def filter_photons(photons, ch_start, WINDOW_LENGTH, dif_win1_win2, noof_pulses,
                         column_st):

    tail_length = WINDOW_LENGTH
    pulse_sep = dif_win1_win2

    f_st = pq_tools.filter_synctimes(photons[:,column_st], ch_start, ch_start + tail_length, noof_pulses, pulse_sep, pq_file = False)
        
    return f_st


def TPQI_analysis_w_invalid_data_marker(Base_Folder_primary, ch0_start, ch1_start, WINDOW_LENGTH, dif_win1_win2, noof_pulses,Base_Folder_secondary = None,
                                    use_invalid_data_marker = False, filter_attempts = False, syncs_per_CR_check = 50, delta_attempts = 50, min_filter_attempts = 1,
                                    return_sn = False , 
                                    contains = 'TPQI',
                                    force_coincidence_evaluation = False,
                                    Verbose = True):
    # Gets coincident photons from Hydraharp data
    coincidences = _aggregated_coincidences(Base_Folder_primary,folder_secondary = Base_Folder_secondary,use_invalid_data_marker = use_invalid_data_marker,syncs_per_CR_check = syncs_per_CR_check,contains=contains,force_coincidence_evaluation = force_coincidence_evaluation, Verbose = Verbose)

    dt_index = 0
    column_st_0 = 1
    column_st_1 = 2
    column_sync_num_ch0 = 3
    colum_invalid = 4

    # Defines the difference in arrival time between the coincident photons
    dts = coincidences[:,dt_index] + (ch1_start - ch0_start)
    dts = dts * 10**(-3)

    if Verbose:
        print
        print 'Found {} coincident photons in all runs.'.format(int(len(coincidences)))
        print '===================================='
        print

    #Filters the coincident photons by selecting only the photons emitted by the NV center
    filtered = filter_coincidences(coincidences, ch0_start, ch1_start, WINDOW_LENGTH,
                                                 dif_win1_win2, noof_pulses, column_st_0, column_st_1)

    if filter_attempts:
        filtered_attempts = filter_syncnum(coincidences, syncs_per_CR_check, min_filter_attempts = min_filter_attempts, delta_attempts = delta_attempts, column_sn = column_sync_num_ch0)
        filtered = filtered & filtered_attempts

    if use_invalid_data_marker:
        filtered_invalid = coincidences[:,colum_invalid] == 0
        filtered = filtered & filtered_invalid

    filtered_dts = dts[filtered]

    if Verbose:
        print
        print 'Found {} coincident photons after filtering.'.format(int(sum(filtered)))
        print '===================================='
        print

    return dts, filtered_dts

# Apply the same filters as for the TPQI, but return the total tail counts
def TPQI_analysis_tail(Base_Folder_primary, ch0_start, ch1_start, WINDOW_LENGTH, dif_win1_win2, noof_pulses, 
                                    filter_attempts = False, syncs_per_CR_check = 50, delta_attempts = 50, min_filter_attempts = 1,
                                    return_sn = False , 
                                    contains = 'TPQI',
                                    force_coincidence_evaluation = False,
                                    Verbose = True):
    # Gets coincident photons from Hydraharp data
    return get_photons_from_folder(Base_Folder_primary,ch0_start,ch1_start, WINDOW_LENGTH,dif_win1_win2, noof_pulses,syncs_per_CR_check, min_filter_attempts,delta_attempts,contains=contains)

# Added code to analyse based on CR check counts
def TPQI_analysis_w_CR_check_filtering(Base_Folder_primary, ch0_start, ch1_start, WINDOW_LENGTH, dif_win1_win2, 
                                        noof_pulses, contains = 'TPQI',return_sn = False , Verbose = True, load_cr_check = False, load_TPQI_attempts = False,
                                        load_adwin_comm_time =False, Base_Folder_secondary = None, No_of_LDE_attempts = 100, filter_before_or_after = 'before', 
                                        cr_check_min_s = 0, cr_check_min_p = 0, adwin_com_max_p = 1e6, adwin_com_max_s = 1e6,
                                        min_attempts = 0, max_attempts = 99,):
    # Gets coincident photons from Hydraharp data
    coincidences = _aggregated_coincidences(Base_Folder_primary,folder_secondary = Base_Folder_secondary,contains = contains,
                                                                load_cr_check = load_cr_check, 
                                                                load_adwin_comm_time = load_adwin_comm_time,
                                                                No_of_LDE_attempts = No_of_LDE_attempts,
                                                                load_TPQI_attempts = load_TPQI_attempts)

    dt_index = 0
    column_st_0 = 1
    column_st_1 = 2
    column_sync_num_ch0 = 3
    
    column_cr_before_p = 4
    column_cr_after_p = 5
    column_cr_before_s = 6
    column_cr_after_s = 7
    column_no_of_sequences = 8
    columns_cr = [column_cr_before_p, column_cr_after_p, column_cr_before_s, column_cr_after_s]

    offset = (4 if load_cr_check else 0)
    column_adwin_com_p = 4 + offset
    column_adwin_com_s = 5 + offset
    columns_adwin_com = [column_adwin_com_p, column_adwin_com_s] 

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

    #Filters the coincident photons by setting cr check thresholds
    load_secondary = Base_Folder_secondary != None
    cr_check_filter = filter_cr_check(load_cr_check, columns_cr, coincidences, cr_check_min_p, cr_check_min_s, filter_before_or_after, load_secondary)
    TPQI_attempt_filter = filter_no_of_attempts(load_TPQI_attempts,column_no_of_sequences,coincidences,min_attempts,max_attempts)
    Adwin_comm_time_filter = filter_adwin_comm_time(load_adwin_comm_time,columns_adwin_com, coincidences, adwin_com_max_p, adwin_com_max_s, load_secondary)
    
    overall_filter = is_sync_time_filter & cr_check_filter & TPQI_attempt_filter & Adwin_comm_time_filter

    filtered_dts = dts[overall_filter]

    if Verbose:
        print
        print 'Found {} coincident photons after filtering.'.format(int(sum(overall_filter)))
        print '===================================='
        print

    cr_checks = []
    adwin_comm_times = []

    if load_cr_check:
        if not load_secondary:
            cr_checks = coincidences[:,[column_cr_before_p,column_cr_after_p]][overall_filter]
        else:
            cr_checks = coincidences[:,[column_cr_before_p,column_cr_after_p,column_cr_before_s,column_cr_after_s]][overall_filter]

    if load_adwin_comm_time:
        if not load_secondary:
            adwin_comm_times = coincidences[:,[column_adwin_com_p]][overall_filter]
        else:
            adwin_comm_times = coincidences[:,[column_adwin_com_p,column_adwin_com_s]][overall_filter]


    return dts, filtered_dts, cr_checks, adwin_comm_times