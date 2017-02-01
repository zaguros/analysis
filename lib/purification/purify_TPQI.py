import numpy as np
import h5py, os
from analysis.lib.pq import pq_tools;reload(pq_tools)
from analysis.lib.m2.ssro import pqsequence; reload(pqsequence)
from analysis.lib.tools import toolbox as tb
from analysis.lib.purification import purify_pq; reload(purify_pq)
#import Settings, Filter
import matplotlib.pyplot as plt

def get_coincidences_and_adwin_data_from_folder(folder_primary,save = True,folder_secondary = None, index = 1, 
                                                load_cr_check = False, load_TPQI_attempts = False, load_adwin_comm_time = False,
                                                No_of_LDE_attempts = 100):

    sync_num_name = 'PQ_sync_number-' + str(index)
    # print 'this is the save!', save
    filepaths_primary = tb.get_all_msmt_filepaths(folder_primary)

    if folder_secondary != None:
        load_secondary = True
        filepaths_secondary = tb.get_all_msmt_filepaths(folder_secondary) 
    else:
        load_secondary = False

    co = np.ones([1,4])

    # print filepaths
    for i,fp in enumerate(filepaths_primary):
        # print f
        fp = os.path.split(fp)[0]
        
        purifyPQ_p = purify_pq.purifyPQAnalysis(fp)
        if load_secondary:
            fs = os.path.split(filepaths_secondary[i])[0]
            purifyPQ_s = purify_pq.purifyPQAnalysis(fs)

        if sync_num_name in purifyPQ_p.pqf.keys():

            if i == 0:
                coinc = pq_tools.get_coincidences(purifyPQ_p.pqf)
            else:
                coinc = pq_tools.get_coincidences(purifyPQ_p.pqf,save = save,force_coincidence_evaluation=False)

            new_vars = np.array([])

            if load_cr_check:
                # Fixed number of LDE attempts per repetition, so straightforward to pull CR check value index
                cr_check_ind = (np.floor(coinc[:,3]/No_of_LDE_attempts)).astype(int)

                # print cr_check_ind
                cr_before = (purifyPQ_p.agrp['CR_before'].value)[cr_check_ind]
                cr_after = (purifyPQ_p.agrp['CR_after'].value)[cr_check_ind]

                new_vars = np.hstack((new_vars, cr_before[:, np.newaxis],cr_after[:, np.newaxis])) if new_vars.size else cr_before[:, np.newaxis],cr_after[:, np.newaxis]

                if load_secondary:
                    
                    cr_before = (purifyPQ_s.agrp['CR_before'].value)[cr_check_ind]
                    cr_after = (purifyPQ_s.agrp['CR_after'].value)[cr_check_ind]
                    new_vars = np.hstack((new_vars, cr_before[:, np.newaxis],cr_after[:, np.newaxis]))

            if load_adwin_comm_time:
                # Fixed number of LDE attempts per repetition, so straightforward to pull CR check value index
                ind = (np.floor(coinc[:,3]/No_of_LDE_attempts)).astype(int)

                adwin_comm_time = (purifyPQ_p.agrp['adwin_communication_time'].value)[ind]

                new_vars = np.hstack((new_vars, adwin_comm_time[:, np.newaxis])) if new_vars.size else adwin_comm_time[:, np.newaxis]

                if load_secondary:
                    
                    adwin_comm_time = (purifyPQ_s.agrp['adwin_communication_time'].value)[ind]
                    new_vars = np.hstack((new_vars,adwin_comm_time[:, np.newaxis]))

                # add_vars = np.hstack((add_vars, new_vars)) if add_vars.size else new_vars

            if load_TPQI_attempts:
                TPQI_attempts = (coinc[:,3]%No_of_LDE_attempts).astype(int)

                new_vars = np.hstack((new_vars, TPQI_attempts[:, np.newaxis])) if new_vars.size else TPQI_attempts[:,np.newaxis]


            if i == 0:
        
                co = np.hstack((coinc,new_vars)) if new_vars.size else coinc

            else:

                if co[0,3] == 1:
                    
                    co = np.hstack((coinc,new_vars)) if new_vars.size else coinc
                else:
                   
                    co = np.vstack((co, np.hstack((coinc,new_vars)))) if new_vars.size else np.vstack((co, coinc))
    
    return co


def filter_syncnum(pqf, attempts, max_attempts, first, pq_file = True, index = 1):
    """
    Return a filter for sync times in between t0 and t1
    """

    sync_num_name = '/PQ_sync_number-' + str(index)

    if pq_file:
        sync_num = pqf[sync_num_name].value
    else:
        sync_num = pqf
    if first ==1:
        fltr = sync_num % attempts <= max_attempts 
    elif first ==2:
        fltr = sync_num % attempts >= max_attempts
    return fltr
 

def get_coincidences_from_folder(folder, attempts, max_attempts, first, index = 1,save = True,contains = '', force_coincidence_evaluation = False):
    sync_num_name = 'PQ_sync_number-' + str(index)
    # print 'this is the save!', save
    filepaths = tb.get_all_msmt_filepaths(folder) 

    if contains != '':
        new_fps = []
        for f in filepaths:
            if contains in f:
                new_fps.append(f)
        filepaths = new_fps

    co = np.ones([1,4])
    # print filepaths
    for i,f in enumerate(filepaths):
        if i == 0:
            pqf = pq_tools.pqf_from_fp(f, rights = 'r+')

            makeValeriaFilter = filter_syncnum(pqf,attempts,max_attempts, first)
            # print 'the filter is', makeValeriaFilter
            if sync_num_name in pqf.keys():
                co = pq_tools.get_coincidences(pqf, fltr0 = makeValeriaFilter, fltr1 = makeValeriaFilter, force_coincidence_evaluation=force_coincidence_evaluation)  
            # for j in makeValeriaFilter:
            #     if j:
            #         if sync_num_name in pqf.keys():
            #             co = pq_tools.get_coincidences(pqf)                   
        else:
            pqf = pq_tools.pqf_from_fp(f, rights = 'r+')
            makeValeriaFilter = filter_syncnum(pqf,attempts,max_attempts,first)
            # print 'the filter is', makeValeriaFilter

            if sync_num_name in pqf.keys():
                if co[0,3] == 1:
                    co = pq_tools.get_coincidences(pqf, fltr0 = makeValeriaFilter, fltr1 = makeValeriaFilter,force_coincidence_evaluation=force_coincidence_evaluation,save = save)
                else:
                    co = np.vstack((co, pq_tools.get_coincidences(pqf, fltr0 = makeValeriaFilter, fltr1 = makeValeriaFilter,force_coincidence_evaluation=force_coincidence_evaluation,save = save)))
               

    return co



def _aggregated_coincidences(Base_Folder, force_coincidence_evaluation = False,**kw):

    load_cr_check = kw.get('load_cr_check',False)
    load_adwin_comm_time = kw.get('load_adwin_comm_time',False)
    contains = kw.pop('contains','TPQI')
    if load_cr_check or load_adwin_comm_time:
        in_coincidences = get_coincidences_and_adwin_data_from_folder(Base_Folder,save=True, **kw)

    else:
        in_coincidences = np.empty((0,4))
        in_coincidences = np.vstack((in_coincidences, pq_tools.get_coincidences_from_folder(Base_Folder,contains = contains,save=True,force_coincidence_evaluation = force_coincidence_evaluation)))

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
                                    return_sn = False , 
                                    contains = 'TPQI',
                                    force_coincidence_evaluation = False,
                                    Verbose = True):
    # Gets coincident photons from Hydraharp data
    coincidences = _aggregated_coincidences(Base_Folder_primary,contains=contains,force_coincidence_evaluation = force_coincidence_evaluation)

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

def TPQI_analysis_w_extra_filtering(Base_Folder_primary, ch0_start, ch1_start, WINDOW_LENGTH, dif_win1_win2, noof_pulses, attempts = 50, max_attempts = 1, first = 1,
                                    return_sn = False , 
                                    contains = 'TPQI', force_coincidence_evaluation = False,
                                    Verbose = True ):
    # Gets coincident photons from Hydraharp data
    coincidences = get_coincidences_from_folder(Base_Folder_primary, attempts, max_attempts, first, force_coincidence_evaluation =force_coincidence_evaluation,contains=contains)

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