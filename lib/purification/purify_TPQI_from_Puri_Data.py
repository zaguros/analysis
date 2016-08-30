import numpy as np
import h5py, os
from analysis.lib.pq import pq_tools;reload(pq_tools)
from analysis.lib.m2.ssro import pqsequence; reload(pqsequence)
from analysis.lib.tools import toolbox as tb
from analysis.lib.purification import purify_pq; reload(purify_pq)
#import Settings, Filter
import matplotlib.pyplot as plt

import purify_analysis_params as analysis_params;reload(analysis_params)

def get_coincidences(pqf, index = 1, fltr0=None, fltr1=None, pulse_offset = 0, force_coincidence_evaluation = False, save = True):

    sync_time_name = '/PQ_sync_time-' + str(index)
    tot_time_name =  '/PQ_time-' + str(index)
    sync_num_name = '/PQ_sync_number-' + str(index)

    if pulse_offset == 0:
        if pq_tools.has_analysis_data(pqf, 'coincidences') and not force_coincidence_evaluation:
            c, c_attrs = pq_tools.get_analysis_data(pqf, 'coincidences')
            return c    
    else:
        if pq_tools.has_analysis_data(pqf, 'coincidences_offset_'+str(pulse_offset)) and not force_coincidence_evaluation:
            c, c_attrs = pq_tools.get_analysis_data(pqf,'coincidences_offset_'+str(pulse_offset))
            return c   

    sync_time = pqf[sync_time_name].value
    total_time = pqf[tot_time_name].value
    sync_number = pqf[sync_num_name].value

    is_ph0, is_ph1 = pq_tools.get_photons(pqf)

    # thin down a bit with loose filtering
    if fltr0 != None:
        fltr0 = is_ph0 & fltr0
    else:
        fltr0 = is_ph0

    if fltr1 != None:
        fltr1 = is_ph1 & fltr1
    else:
        fltr1 = is_ph1

    st0 = sync_time[fltr0]
    t0  = total_time[fltr0]
    sn0 = sync_number[fltr0]
    
    st1 = sync_time[fltr1]
    t1  = total_time[fltr1]
    sn1 = sync_number[fltr1] + pulse_offset
   
    samesync0 = np.in1d(sn0, sn1)
    samesync1 = np.in1d(sn1, sn0)

    c_st0 = st0[samesync0]
    c_st1 = st1[samesync1]
    c_t0 = t0[samesync0]
    c_sn0 = sn0[samesync0]
    c_t1 = t1[samesync1]
    c_sn1 = sn1[samesync1]

    ### Code for afterpulsing (PH 16) - also need to uncomment one line in for loop
    # unique_vals, uniq_counts = np.unique(sn1, return_counts = True)
    # repeated = np.in1d(sn1,unique_vals[uniq_counts > 1])
    # c_st0 = st1[repeated]
    # c_st1 = st1[repeated]
    # c_t0 = t1[repeated]
    # c_sn0 = sn1[repeated]
    # c_t1 = t1[repeated]
    # c_sn1 = sn1[repeated]

    coincidences = np.empty((0,4))
    for i, (_sn0, _t0, _st0) in enumerate(zip(c_sn0, c_t0, c_st0)):

        _c = c_sn1==_sn0
        # _c[i] = 0 # Not the same entry (for afterpulsing)
        
        for _t1, _st1 in zip(c_t1[_c], c_st1[_c]):
            dt = int(_st0) - int(_st1)
            coincidences = np.vstack((coincidences, np.array([dt, _st0, _st1, _sn0])))

    # coincidences = np.empty((0,4))
    # for _sn0, _t0, _st0 in zip(c_sn0, c_t0, c_st0):
    #     _c = c_sn1==_sn0
        
    #     for _t1, _st1 in zip(c_t1[_c], c_st1[_c]):
    #         dt = int(_t0) - int(_t1)
    #         coincidences = np.vstack((coincidences, np.array([dt, _st0, _st1, _sn0])))
    
    if save:
        
        if pulse_offset == 0:
            pq_tools.set_analysis_data(pqf, 'coincidences', coincidences,
                          columns=('dt = ch0-ch1 (bins)', 'sync_time ch0 (bins)', 'sync_time ch1 (bins)', 'sync number'))
        else:
            pq_tools.set_analysis_data(pqf, 'coincidences_offset_'+str(pulse_offset), coincidences,
                          columns=('dt = ch0-ch1 (bins)', 'sync_time ch0 (bins)', 'sync_time ch1 (bins)', 'sync number'))    

    return coincidences


def get_coincidences_from_folder(folder,lt4_timestamps,timing_offsets,offsets_ch1, index = 1,save = True,pulse_offset=0):

    sync_num_name = 'PQ_sync_number-' + str(index)
    # print 'this is the save!', save

    co = np.ones([1,4])

    for i,t_lt4 in enumerate(lt4_timestamps):
        f = tb.data_from_time(t_lt4,folder = folder)
        f = tb.get_msmt_fp(f)
        pqf = pq_tools.pqf_from_fp(f, rights = 'r+')

        if sync_num_name in pqf.keys():
            if i == 0:
                co = get_coincidences(pqf,pulse_offset=pulse_offset) 
                co[:,(1,2)]=  co[:,(1,2)] + timing_offsets[i]
                co[:,2]=  co[:,2] + offsets_ch1[i]        
            else:
                co_temp = get_coincidences(pqf,pulse_offset=pulse_offset,save = save)
                co_temp[:,(1,2)]=  co_temp[:,(1,2)] + timing_offsets[i]
                co_temp[:,2]=  co_temp[:,2] + offsets_ch1[i]
                co = np.vstack((co, co_temp))
                    
    return co


def filter_coincidences(coincidences, ch0_start, ch1_start, WINDOW_LENGTH,
                         column_st_0, column_st_1):

    tail_length = WINDOW_LENGTH
    
    
    f_st0 = pq_tools.filter_synctimes(coincidences[:,column_st_0], ch0_start, ch0_start + tail_length, 1, 0, pq_file = False)
    f_st1 = pq_tools.filter_synctimes(coincidences[:,column_st_1], ch1_start, ch1_start +  tail_length, 1, 0, pq_file = False)
        
    return f_st0 & f_st1 # & f_dt

def TPQI_analysis(Base_Folder_primary,ts_list,timing_offsets,offsets_ch1, ch0_start = None, ch1_start = None, WINDOW_LENGTH = None, 
                                    return_sn = False , Verbose = True, ):


    if ch0_start == None:
        ch0_start = analysis_params.filter_settings['st_start']
    if ch1_start == None:
        ch1_start = analysis_params.filter_settings['st_start']
    if WINDOW_LENGTH == None:
        WINDOW_LENGTH = analysis_params.filter_settings['st_len']
    
    # Gets coincident photons from Hydraharp data
    coincidences = get_coincidences_from_folder(Base_Folder_primary,ts_list,timing_offsets,offsets_ch1)
    offset_coincidences = get_coincidences_from_folder(Base_Folder_primary,ts_list,timing_offsets,offsets_ch1,pulse_offset=1)

    dt_index = 0
    column_st_0 = 1
    column_st_1 = 2
    column_sync_num_ch0 = 3

    # Defines the difference in arrival time between the coincident photons
    dts = coincidences[:,dt_index] + (ch1_start - ch0_start)
    dts = dts * 10**(-3)

    # Defines the difference in arrival time between the coincident photons
    dts_offset = offset_coincidences[:,dt_index] + (ch1_start - ch0_start)
    dts_offset = dts_offset * 10**(-3)


    if Verbose:
        print
        print 'Found {} coincident photons in all runs.'.format(int(len(coincidences)))
        print '===================================='
        print


    #Filters the coincident photons by selecting only the photons emitted by the NV center
    is_sync_time_filter = filter_coincidences(coincidences, ch0_start, ch1_start, WINDOW_LENGTH,
                                                column_st_0, column_st_1)
    is_sync_time_filter_offset = filter_coincidences(offset_coincidences, ch0_start, ch1_start, WINDOW_LENGTH,
                                                column_st_0, column_st_1)

    filtered_dts = dts[is_sync_time_filter]
    filtered_dts_offset = dts_offset[is_sync_time_filter_offset]

    if Verbose:
        print
        print 'Found {} coincident photons after filtering.'.format(int(sum(is_sync_time_filter)))
        print '===================================='
        print

    return dts,dts_offset, filtered_dts,filtered_dts_offset

def find_tstamps_of_day(ts_list,day_string,contains='XX',analysis_folder = 'throw exception',newest_tstamp = '235959'):

    latest_t = day_string + newest_tstamp # where in the day do you want to begin? 235959 mean: take the whole day
    newer_than = day_string+'_000000'

    while tb.latest_data(contains,older_than=latest_t,folder=analysis_folder,newer_than=newer_than,return_timestamp = True,raise_exc = False) != False:

        latest_t,f = tb.latest_data(contains,older_than=latest_t,folder=analysis_folder,newer_than=newer_than,return_timestamp = True,raise_exc=False)
        
        ### debug statement that prints the full timestamp and the relevant identifier.
        # print latest_t[8:],latest_t

        ### append found timestamp to list of timestamps
        ts_list.append(latest_t) 

    return ts_list


def get_tstamps_and_offsets(analysis_folder, contains = 'Purify', verbose = False, unshifted_days = None,shifted_days = None,shifted_data_correction_time = None, shifted_data_start_offset_ch1 = None):
    all_lt4 = []
    offsets,offsets_ch1 = [],[] # Hold offset to compensate for timing change for new APD

    if shifted_data_correction_time == None:
        shifted_data_correction_time = analysis_params.data_settings['shifted_data_correction_time']

    if shifted_data_start_offset_ch1 == None:
        shifted_data_start_offset_ch1 = analysis_params.data_settings['shifted_data_start_offset_ch1']

    if unshifted_days == None:
        unshifted_days = analysis_params.data_settings['unshifted_days']

    if shifted_days == None:
        shifted_days = analysis_params.data_settings['shifted_days']

    for d in unshifted_days+shifted_days:
            if verbose:
                    print d
            tstamp_lt4 = find_tstamps_of_day([],d,contains = contains,analysis_folder = analysis_folder)
            #,newest_tstamp = '110000') ### newest timestamp allows for only taking parts of a day.
            all_lt4.extend(tstamp_lt4)
            if d in shifted_days:
                if verbose:
                    print 'shifting ',d
                offsets.extend(np.zeros(np.shape(tstamp_lt4))+ shifted_data_correction_time)
                offsets_ch1.extend(np.zeros(np.shape(tstamp_lt4))+ shifted_data_start_offset_ch1)           
            else:
                offsets.extend(np.zeros(np.shape(tstamp_lt4)))
                offsets_ch1.extend(np.zeros(np.shape(tstamp_lt4)))
            if verbose:
                print 'Found ' + str(len(tstamp_lt4)) + ' timestamps!'

    return all_lt4, offsets, offsets_ch1

