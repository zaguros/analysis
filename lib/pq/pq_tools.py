import numpy as np
from matplotlib import pyplot as plt
import h5py
from analysis.lib.tools import toolbox as tb

def get_photons(pqf, index = 1, pq_device = ''):
    """
    returns two filters (1d-arrays): whether events are ch0-photons/ch1-photons
    """

    chan_name = pq_device + '/PQ_channel-' + str(index)
    spec_name = pq_device + '/PQ_special-' + str(index)

    if type(pqf) == h5py._hl.files.File:

        channel = pqf[chan_name].value
        special = pqf[spec_name].value

    elif type(pqf) == str:

        f = h5py.File(pqf,'r')
        channel = f[chan_name].value
        special = f[spec_name].value
        f.close()

    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise    


    is_not_special = special == 0
    is_channel_0 = channel == 0
    is_channel_1 = channel == 1
    
    is_photon_0 = np.logical_and(is_not_special, is_channel_0)
    is_photon_1 = np.logical_and(is_not_special, is_channel_1)
    
    return is_photon_0, is_photon_1

def get_length_block(pqf, index = 1,  pq_device = ''):
    """
    Returns the number of elements in one block
    """
    sync_num_name =  pq_device + '/PQ_sync_number-' + str(index)
    
    # Determines how long the block is
    if type(pqf) == h5py._hl.files.File: 
        sync_num = pqf[sync_num_name]
        num_elements_block = len(sync_num)
   
    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        sync_num = f[sync_num_name]
        num_elements_block = len(sync_num)
        f.close()
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise 

    return num_elements_block



def get_markers(pqf, chan, index = 1, pq_device = ''):
    """
    returns a filter (1d-array): whether events are markers on the given channel
    """

    chan_name = pq_device + '/PQ_channel-' + str(index)
    spec_name = pq_device + '/PQ_special-' + str(index)
    
    if type(pqf) == h5py._hl.files.File:

        channel = pqf[chan_name].value
        special = pqf[spec_name].value

    elif type(pqf) == str:

        f = h5py.File(pqf,'r')
        channel = f[chan_name].value
        special = f[spec_name].value
        f.close()
    
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise

    is_special = special == 1
    is_channel = channel == chan

    return (is_special & is_channel)

def get_markers_with_num_markers(pqf, chan, index = 1, pq_device = ''):
    """
    returns a filter (1d-array): whether events are markers on the given channel
    """

    chan_name = pq_device + '/PQ_channel-' + str(index)
    spec_name = pq_device + '/PQ_special-' + str(index)
    
    if type(pqf) == h5py._hl.files.File:

        channel = pqf[chan_name].value
        special = pqf[spec_name].value

    elif type(pqf) == str:

        f = h5py.File(pqf,'r')
        channel = f[chan_name].value
        special = f[spec_name].value
        f.close()
    
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise

    is_special = special == 1
    is_channel = channel == chan
    is_marker = is_special & is_channel
    num_markers = len(channel[is_marker])

    return is_marker, num_markers


def get_rndm_num(pqf, chan_rnd_0, chan_rnd_1, num_blocks, add_seg_length_old, index = 1, pq_device = ''):
    """
    returns a filter (1d-array): whether events are markers on the given channel
    """

    add_seg_length = 0

    chan_name = pq_device + '/PQ_channel-' + str(index)
    spec_name = pq_device + '/PQ_special-' + str(index)
    sync_num_name = pq_device + '/PQ_sync_number-' + str(index)
    
    if index != num_blocks:
        if type(pqf) == h5py._hl.files.File:

            channel = pqf[chan_name].value
            special = pqf[spec_name].value
            sync_num = pqf[sync_num_name].value

        elif type(pqf) == str:

            f = h5py.File(pqf,'r')
            channel = f[chan_name].value
            special = f[spec_name].value
            sync_num = f[sync_num_name].value
            f.close()
        
        else:
            print "Neither filepath nor file enetered in function please check:", pqf
            raise

        if add_seg_length_old > 0:
            channel = channel[add_seg_length_old::]
            special = special[add_seg_length_old::]
            sync_num = sync_num[add_seg_length_old::]

        last_sync_num = sync_num[-1]

        index = index + 1
        chan_name = pq_device + '/PQ_channel-' + str(index)
        spec_name = pq_device + '/PQ_special-' + str(index)
        sync_num_name = pq_device + '/PQ_sync_number-' + str(index)

        
        if type(pqf) == h5py._hl.files.File:

            channel2 = pqf[chan_name].value
            special2 = pqf[spec_name].value
            sync_num2 = pqf[sync_num_name].value

        elif type(pqf) == str:

            f = h5py.File(pqf,'r')
            channel2 = f[chan_name].value
            special2 = f[spec_name].value
            sync_num2 = f[sync_num_name].value
            f.close()
        
        else:
            print "Neither filepath nor file enetered in function please check:", pqf
            raise

        ad_channel = []

        if sync_num2[0] == last_sync_num:
            is_last_sync_num = sync_num2 == last_sync_num
            ad_sync_num = sync_num2[is_last_sync_num]
            ad_channel = channel2[is_last_sync_num]
            ad_special = special2[is_last_sync_num]

            channel = np.hstack((channel, ad_channel))
            special = np.hstack((special, ad_special))
            sync_num = np.hstack((sync_num, ad_sync_num))

        add_seg_length  = len(ad_channel)

    else:
        if type(pqf) == h5py._hl.files.File:

            channel = pqf[chan_name].value
            special = pqf[spec_name].value
            sync_num = pqf[sync_num_name].value

        elif type(pqf) == str:

            f = h5py.File(pqf,'r')
            channel = f[chan_name].value
            special = f[spec_name].value
            sync_num = f[sync_num_name].value
            f.close()
        
        else:
            print "Neither filepath nor file enetered in function please check:", pqf
            raise

        if add_seg_length_old > 0:
            channel = channel[add_seg_length_old::]
            special = special[add_seg_length_old::]
            sync_num = sync_num[add_seg_length_old::]

    is_rnd = (special == 0 ) & (channel == 1)

    is_special = special == 1
    is_channel_rnd_0 = channel == chan_rnd_0
    is_channel_rnd_1 = channel == chan_rnd_1

    # returns a filter for the raw data which is true for every data point for which there is a random number
    sync_num_rnd = np.unique(sync_num[is_rnd])
    is_data_with_rndm_num = np.in1d(sync_num, sync_num_rnd)
    
    is_rnd_0 = (is_special & is_channel_rnd_0) & is_data_with_rndm_num
    is_rnd_1 = (is_special & is_channel_rnd_1) & is_data_with_rndm_num
    dif_should_and_get = -len(sync_num_rnd) + np.sum(is_rnd_1) + np.sum(is_rnd_0)

    print "Dif what it should be what we get", dif_should_and_get

    sync_num_final_rnd = sync_num[(is_rnd_0 | is_rnd_1)]
    print len(sync_num_rnd), len(np.unique(sync_num_rnd))
    print len(sync_num_final_rnd), np.sum(is_rnd_1) + np.sum(is_rnd_0)

    len_check = int(len(sync_num_final_rnd)/10.)
    num_unique_sync_nums = 0

    last_number = 0
    indices_double_sync_nums = []
    indices_to_be_saved = []


    for i in range(10):
        if i != 9:
            temp_check = sync_num_final_rnd[(i*len_check):((i+1)*len_check)]
            first_number = temp_check[0]
            if first_number == last_number:
                print "One less unique number"
            last_number = temp_check[len(temp_check)-1]   
            unique_sync_nums, index_sync_nums = np.unique(temp_check, return_index = True)
            index_sync_nums = np.sort(index_sync_nums + (i* len_check))
            for j in range(len(index_sync_nums)):
                if j != len(index_sync_nums) - 1:
                    test = index_sync_nums[j+1] - index_sync_nums[j]
                    if test > 1:
                        for k in range(test-1):
                            indices_to_be_saved.append(index_sync_nums[j]+k+1)             
            num_unique_sync_nums = num_unique_sync_nums + len(unique_sync_nums)
        else:
            temp_check = sync_num_final_rnd[(i*len_check)::]
            first_number = temp_check[0]
            if first_number == last_number:
                print "One less unique number"
            last_number = temp_check[len(temp_check)-1] 
            unique_sync_nums, index_sync_nums = np.unique(temp_check, return_index = True)
            index_sync_nums = np.sort(index_sync_nums + (i* len_check))
            for j in range(len(index_sync_nums)):
                if j != len(index_sync_nums) - 1:
                    test = index_sync_nums[j+1] - index_sync_nums[j]
                    if test > 1:
                        for k in range(test-1):
                            indices_to_be_saved.append(index_sync_nums[j]+k+1)
            num_unique_sync_nums = num_unique_sync_nums + len(unique_sync_nums)
            
           
    print "Length unique sync numbers", num_unique_sync_nums
    print "The extra sync numbers", len(sync_num_final_rnd) - num_unique_sync_nums

    print len(indices_to_be_saved)
    print len(np.unique(indices_to_be_saved))
    corresp_sync_num = []

    for i,j in enumerate(indices_to_be_saved):
        corresp_sync_num.append(sync_num_final_rnd[j])

    Additional_sync_num, index_un_add_sync_num = np.unique(corresp_sync_num, return_index = True)
    print len(Additional_sync_num)
    print Additional_sync_num

    return is_rnd_0, is_rnd_1, add_seg_length

def get_multiple_photon_syncs(pqf, index = 1, pq_device = ''):

    spec_name = pq_device + '/PQ_special-' + str(index)
    sync_num_name = pq_device + '/PQ_sync_number-' + str(index)

    special = pqf[spec_name].value
    sync_numbers = pqf[sync_num_name].value

    is_photon = special == 0
    photon_sync_numbers = sync_numbers[is_photon]
    
    # this works nicely on sorted arrays
    is_multiple_photon_sync = photon_sync_numbers[1:] == photon_sync_numbers[:-1]
    #multiple_photon_sync_numbers = photon_sync_numbers[is_multiple_photon_sync]

    return is_multiple_photon_sync

def get_coincidences(pqf, index = 1, fltr0=None, fltr1=None, force_coincidence_evaluation = False, save = True, pq_device = ''):
    sync_time_name = pq_device + '/PQ_sync_time-' + str(index)
    tot_time_name =  pq_device + '/PQ_time-' + str(index)
    sync_num_name = pq_device + '/PQ_sync_number-' + str(index)

    if has_analysis_data(pqf, 'coincidences') and not force_coincidence_evaluation:
        c, c_attrs = get_analysis_data(pqf, 'coincidences')
        return c    

    sync_time = pqf[sync_time_name].value
    total_time = pqf[tot_time_name].value
    sync_number = pqf[sync_num_name].value

    is_ph0, is_ph1 = get_photons(pqf)
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
    sn1 = sync_number[fltr1]
    # print len(st0),len(t0),len(sn0),len(st1),len(t1),len(sn1),

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
            dt = int(_t0) - int(_t1)
            coincidences = np.vstack((coincidences, np.array([dt, _st0, _st1, _sn0])))

    # coincidences = np.empty((0,4))
    # for _sn0, _t0, _st0 in zip(c_sn0, c_t0, c_st0):
    #     _c = c_sn1==_sn0
        
    #     for _t1, _st1 in zip(c_t1[_c], c_st1[_c]):
    #         dt = int(_t0) - int(_t1)
    #         coincidences = np.vstack((coincidences, np.array([dt, _st0, _st1, _sn0])))
    
    if save:
        set_analysis_data(pqf, 'coincidences', coincidences,
                          columns=('dt = ch0-ch1 (bins)', 'sync_time ch0 (bins)', 'sync_time ch1 (bins)', 'sync number'))
                       
    return coincidences


def get_coincidences_from_folder(folder, index = 1, save = True, contains = '', force_coincidence_evaluation = False, pq_device = ''):

    sync_num_name = pq_device + 'PQ_sync_number-' + str(index)
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
            pqf = pqf_from_fp(f, rights = 'r+')
            if sync_num_name in pqf.keys():
                co = get_coincidences(pqf)            
        else:
            pqf = pqf_from_fp(f, rights = 'r+')

            if sync_num_name in pqf.keys():
                if co[0,3] == 1:
                    co = get_coincidences(pqf,save = save, force_coincidence_evaluation = force_coincidence_evaluation)
                else:
                    co = np.vstack((co, get_coincidences(pqf,save = save, force_coincidence_evaluation = force_coincidence_evaluation)))
                    
    return co

def get_photons_in_sync_windows(pqf, first_win_min, first_win_max, second_win_min, second_win_max, index = 1, pq_device = '', VERBOSE = True):
    """
    Returns two filters whether events are in the first or 
    in the second time window.
    """

    sync_time_name = pq_device + '/PQ_sync_time-' + str(index)
    spec_name = pq_device + '/PQ_special-' + str(index)

    if type(pqf) == h5py._hl.files.File: 
        sync_time = pqf[sync_time_name].value
        special = pqf[spec_name].value

    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        sync_time = f[sync_time_name].value
        special = f[spec_name].value
        f.close()
    
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise   

    is_photon = special == 0
    
    is_event_first_window = (sync_time > first_win_min) & \
                         (sync_time <= first_win_max)
    is_event_second_window = (sync_time > second_win_min) & \
                          (sync_time <= second_win_max)

    is_photon_first_window = is_photon & is_event_first_window
    is_photon_second_window = is_photon & is_event_second_window
    
    is_photon_check = is_photon_first_window | is_photon_second_window
        
    if VERBOSE:
        if sum(is_photon_check) != sum(is_photon):
            print "Not all detected photons are in the broad windows set"
    
    return is_photon_first_window, is_photon_second_window

def get_tail_filtered_photons(pqf, first_win_min_ch0, dif_win1_win2, window_length, dif_ch0_ch1, index = 1, pq_device = '', VERBOSE = True):
    """
    Returns two filters whether events are in the first or 
    in the second tail.
    """

    chan_name = pq_device + '/PQ_channel-' + str(index)
    sync_time_name = pq_device + '/PQ_sync_time-' + str(index)
    spec_name = pq_device + '/PQ_special-' + str(index)

    if type(pqf) == h5py._hl.files.File: 
        channel = pqf[chan_name].value
        sync_time = pqf[sync_time_name].value
        special = pqf[spec_name].value

    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        channel = f[chan_name].value
        sync_time = f[sync_time_name].value
        special = f[spec_name].value
        f.close()
    
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise   

    is_photon = special == 0
    is_ch0 = channel == 0
    is_ch1 = channel == 1
    
    first_win_min = first_win_min_ch0
    first_win_max = first_win_min + window_length
    second_win_min = first_win_min + dif_win1_win2
    second_win_max = first_win_max + dif_win1_win2

    is_event_first_tail_ch0_filt = (sync_time > first_win_min) & \
                         (sync_time <= first_win_max)
    is_event_second_tail_ch0_filt = (sync_time > second_win_min) & \
                          (sync_time <= second_win_max)

    is_event_first_tail_ch1_filt = (sync_time > (first_win_min + dif_ch0_ch1)) & \
                         (sync_time <= (first_win_max + dif_ch0_ch1))
    is_event_second_tail_ch1_filt = (sync_time > (second_win_min + dif_ch0_ch1)) & \
                          (sync_time <= (second_win_max + dif_ch0_ch1))


    is_photon_ch0 = is_photon & is_ch0
    is_photon_ch1 = is_photon & is_ch1

    is_ph_first_tail = (is_photon_ch0 & is_event_first_tail_ch0_filt) | (is_photon_ch1 & is_event_first_tail_ch1_filt)
    is_ph_second_tail = (is_photon_ch0 & is_event_second_tail_ch0_filt) | (is_photon_ch1 & is_event_second_tail_ch1_filt)
    
    if VERBOSE:
        print "The total number of photons detected in the first tail is:", sum(is_ph_first_tail)
        print "The total number of photons detected in the second tail is:", sum(is_ph_second_tail)

    return is_ph_first_tail, is_ph_second_tail

def get_un_sync_num_with_markers(pqf, marker_chan, sync_time_lim = 0, index = 1, pq_device = '', VERBOSE = True):
    """
    Returns a list with the unique sync numbers with a marker on a specific marker channel.
    """
    sync_num_name = pq_device + '/PQ_sync_number-' + str(index)

    if type(pqf) == h5py._hl.files.File: 
        sync_numbers = pqf[sync_num_name].value

    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        sync_numbers = f[sync_num_name].value
        f.close()
    
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise 


    if sync_time_lim > 0:
        sync_num_with_markers = sync_numbers[filter_marker_time_lim(pqf,marker_chan, sync_time_lim, index = index, pq_device = pq_device, VERBOSE = VERBOSE)]
    else:
        sync_num_with_markers = sync_numbers[filter_marker(pqf,marker_chan, index = index, pq_device = pq_device, VERBOSE = VERBOSE)]

    unique_sync_num_with_markers = np.unique(sync_num_with_markers)

    if VERBOSE:
        print "The number of events with a sync number that has a marker is:", len(sync_num_with_markers)
        print "The number of unique sync numbers that have a marker is:", len(unique_sync_num_with_markers)

    return unique_sync_num_with_markers

def get_combined_tail_counts_per_shot(pqf, first_win_min, first_win_max, second_win_min, second_win_max, first_win_min_ch0, dif_win1_win2, window_length, dif_ch0_ch1, index = 1, pq_device = '', VERBOSE = True):
    """
    Retuns the total number of shots, total number of photons per shot, and the tailcounts per shot for the different tails which are filtered using
    the arguments of the function. Finally it returns the total tailcounts per shot
    """

    sync_num_name = pq_device + '/PQ_sync_number-' + str(index)

    if type(pqf) == h5py._hl.files.File: 
        sync_num = pqf[sync_num_name]
        Total_shots = sync_num[len(sync_num)-1]

    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        sync_num = f[sync_num_name]
        Total_shots = sync_num[len(sync_num)-1]
        f.close()
    
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise


    is_photon_first_window, is_photon_second_window = get_photons_in_sync_windows(pqf,
                                                                                    first_win_min,
                                                                                    first_win_max,
                                                                                    second_win_min,
                                                                                    second_win_max,
                                                                                    index = index,
                                                                                    pq_device = pq_device,
                                                                                    VERBOSE = VERBOSE)

    if VERBOSE:
        print "Total number of photons in the first window", np.sum(is_photon_first_window)
        print "Total number of photons in the second window", np.sum(is_photon_second_window)

    Tot_ph_per_shot = ((np.sum(is_photon_first_window)+ np.sum(is_photon_second_window))/float(Total_shots))

    is_photon_first_tail, is_photon_second_tail = get_tail_filtered_photons(pqf, 
                                                                                 first_win_min_ch0, 
                                                                                 dif_win1_win2,
                                                                                 window_length,
                                                                                 dif_ch0_ch1,
                                                                                 index = index,
                                                                                 pq_device = pq_device,
                                                                                 VERBOSE = VERBOSE)

    TC_p_shot_first_tail = ((np.sum(is_photon_first_tail))/float(Total_shots))
    TC_p_shot_second_tail = ((np.sum(is_photon_second_tail))/float(Total_shots))
    TC_p_shot = TC_p_shot_first_tail + TC_p_shot_second_tail

    return Total_shots, Tot_ph_per_shot, TC_p_shot_first_tail, TC_p_shot_second_tail, TC_p_shot
    

def get_tail_filtered_ph_sync_num(pqf, first_win_min_ch0, dif_win1_win2, window_length, dif_ch0_ch1, index = 1, pq_device = '', VERBOSE = True):
    """
    Returns the sync numbers of the photons in the first and second tail. The input necessary is the start of the first
    tail of channel 0, the difference between the two tails, the difference between the two channels and the length of the tial.
    """
    sync_num_name = pq_device + '/PQ_sync_number-' + str(index)

    if type(pqf) == h5py._hl.files.File: 
        sync_num = pqf[sync_num_name].value

    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        sync_num = f[sync_num_name].value
        f.close()
    
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise

    is_photon_first_tail, is_photon_second_tail = get_tail_filtered_photons(pqf,
                                                                       first_win_min_ch0,
                                                                       dif_win1_win2,
                                                                       window_length,
                                                                       dif_ch0_ch1,
                                                                       VERBOSE = VERBOSE)



    is_tail_photons = is_photon_first_tail | is_photon_second_tail
    sync_num_ph_first_tail = sync_num[is_photon_first_tail]
    sync_num_ph_second_tail = sync_num[is_photon_second_tail]
  

    overlapping_sync_num_12 = filter_on_same_sync_number(sync_num_ph_second_tail, sync_num_ph_first_tail)
    overlapping_sync_num_21 = filter_on_same_sync_number(sync_num_ph_first_tail, sync_num_ph_second_tail)

    non_overlapping_sync_num_12 = np.array([not x for x in overlapping_sync_num_12])
    non_overlapping_sync_num_21 = np.array([not x for x in overlapping_sync_num_21])

    if VERBOSE:
        print "The ratio of photons that occur in the first tail and also in the second tail is:", 1. - sum(non_overlapping_sync_num_21)/float(len(non_overlapping_sync_num_21))
        print "The ratio of photons that occur in the second tail and also in the first tail is:", 1. - sum(non_overlapping_sync_num_12)/float(len(non_overlapping_sync_num_12))

    unique_sync_num_first_tail = np.unique(sync_num_ph_first_tail[non_overlapping_sync_num_12])
    unique_sync_num_second_tail = np.unique(sync_num_ph_second_tail[non_overlapping_sync_num_21])

    return unique_sync_num_first_tail, unique_sync_num_second_tail, is_tail_photons



##############################################################################
### Filters
##############################################################################


def filter_synctimes(pqf, t0, t1, window_reps=1, window_period=None, pq_file = True, pq_device = '', index = 1):
    """
    Return a filter for sync times in between t0 and t1
    """

    sync_time_name = pq_device + '/PQ_sync_time-' + str(index)

    if pq_file:
        sync_time = pqf[sync_time_name].value
    else:
        sync_time = pqf

    for r in range(window_reps):
        if r == 0:
            fltr = (sync_time >= t0) & (sync_time <= t1)
        else:
            fltr = fltr | (sync_time >= (r * window_period + t0)) & (sync_time <= (r * window_period + t1))
    return fltr



def filter_on_same_sync_number(source_sync_numbers, target_sync_numbers):
    """
    returns a filter for target_sync_numbers that's true for all sync numbers that are also
    in source_sync_numbers.
    """
    return np.in1d(target_sync_numbers, source_sync_numbers)

def filter_marker(pqf, chan, index = 1, pq_device = '', VERBOSE = True):
    """
    Note: at the moment this filter includes the marker events on which we filter.
    """

    sync_num_name = pq_device + '/PQ_sync_number-' + str(index)
    
    if type(pqf) == h5py._hl.files.File: 
        sync_numbers = pqf[sync_num_name].value
    elif type(pqf) == str:        
        f = h5py.File(pqf, 'r')
        sync_numbers = f[sync_num_name].value
        f.close()
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise

    is_mrkr = get_markers(pqf, chan, index = index, pq_device = pq_device)
    marker_sync_numbers = sync_numbers[is_mrkr]

    if VERBOSE:
        print "The number of markers is:", len(marker_sync_numbers)

    return filter_on_same_sync_number(marker_sync_numbers, sync_numbers)

def filter_marker_time_lim(pqf, chan, sync_time_lim, index = 1, pq_device = '', VERBOSE = True):
    """
    Note: at the moment this filter includes the marker events on which we filter.
    """

    sync_time_name = pq_device + '/PQ_sync_time-' + str(index)
    sync_num_name = pq_device + '/PQ_sync_number-' + str(index)
    
    if type(pqf) == h5py._hl.files.File: 
        sync_numbers = pqf[sync_num_name].value
        sync_times = pqf[sync_time_name].value
    elif type(pqf) == str:        
        f = h5py.File(pqf, 'r')
        sync_numbers = f[sync_num_name].value
        sync_times = f[sync_time_name].value
        f.close()
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise

    is_small_sync_time = sync_times <= sync_time_lim
    is_large_sync_time = sync_times > sync_time_lim
    is_mrkr = get_markers(pqf, chan, index = index, pq_device = pq_device)
    if VERBOSE:
        print "The number of markers is:", len(sync_numbers[is_mrkr])


    marker_sync_num_small = sync_numbers[(is_mrkr & is_small_sync_time)] - 1
    marker_sync_num_large = sync_numbers[(is_mrkr & is_large_sync_time)]
    if len(marker_sync_num_large) > 0:
        marker_sync_numbers = np.concatenate((marker_sync_num_small,marker_sync_num_large))
    else:
        marker_sync_numbers = marker_sync_num_small

    return filter_on_same_sync_number(marker_sync_numbers, sync_numbers)

def get_photons_with_markers(pqf, chan, first_win_min, first_win_max, second_win_min, second_win_max, VERBOSE = False):
    """
    Return two filters whether events are first window photons or second window
    photons with markers.
    """
    
    is_photon_first_window, is_photon_second_window = get_photons_in_sync_windows(pqf,
                            first_win_min, first_win_max, second_win_min, second_win_max, VERBOSE = VERBOSE)
    is_events_with_marker = filter_marker(pqf,chan, VERBOSE = VERBOSE)
    
    is_photon_first_window_with_markers = is_photon_first_window & \
                                            is_events_with_marker
    is_photon_second_window_with_markers = is_photon_second_window &\
                                            is_events_with_marker

    return is_photon_first_window_with_markers, is_photon_second_window_with_markers


##############################################################################
### File management
##############################################################################


def fp_from_pqf(pqf):
    return pqf.filename

def pqf_from_fp(fp , rights = 'r'):
    pqf = h5py.File(fp, rights)
    return pqf   


def set_analysis_data(pqf, name, data, analysisgrp = 'analysis', subgroup=None, **kw):

    agrp = pqf.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))
    if name in agrp.keys():
        del agrp[name]
    agrp[name] = data
    pqf.flush()
    
    for k in kw:
        agrp[name].attrs[k] = kw[k]
    
    pqf.flush()
       

def has_analysis_data(pqf, name, analysisgrp = 'analysis', subgroup=None):
    
    agrp = pqf.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))
    if name in agrp.keys():
        return True
    else:
        return False
    
def get_analysis_data(pqf, name, analysisgrp = 'analysis', subgroup=None):
        
    agrp = pqf.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))

    if name not in agrp.keys():
        return None
    
    dat = agrp[name].value
    attrs = {}
    for (an, av) in agrp[name].attrs.items():
        attrs[an] = av

    return dat, attrs


def delete_analysis_data(pqf, name, analysisgrp = 'analysis', subgroup=None):

    agrp = pqf.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))
    if name in agrp.keys():
        del agrp[name]
    pqf.flush()


##############################################################################
### Photon histograms in time
##############################################################################


def get_photon_hist(pqf, index = 1, **kw):
    
    pq_device = kw.pop('pq_device', '')
    sync_time_name = pq_device + '/PQ_sync_time-' + str(index)

    save = kw.pop('save', False)
    fltr = kw.pop('fltr', None)
    force_eval = kw.pop('force_eval', True)
    start = kw.pop('start', 0) 
    length = kw.pop('length', 1e6) 
    hist_binsize= kw.pop('hist_binsize', 1e3)
    offset = kw.pop('offset',0) 
    offset_ch1 = kw.pop('offset_ch1',0)

    if not force_eval and has_analysis_data(pqf, 'photon_histogram'):
        h, h_attrs = get_analysis_data(pqf, 'photon_histogram')
        be, be_attrs = get_analysis_data(pqf, 'photon_histogram_binedges')
        h0 = h[:,0]
        h1 = h[:,1]
        return (h0, be), (h1, be)
    
    sync_time = pqf[sync_time_name].value
    
    ph0, ph1 = get_photons(pqf, index = index, pq_device = pq_device)

    if fltr != None:
        _fltr0 = (ph0 & fltr)
        _fltr1 = (ph1 & fltr)
    else:
        _fltr0 = ph0
        _fltr1 = ph1
    
    st0 = sync_time[_fltr0] + offset
    st1 = sync_time[_fltr1] + offset + offset_ch1

    binedges = np.arange(start,start+length, hist_binsize)

    h0, b0 = np.histogram(st0, bins=binedges)
    h1, b1 = np.histogram(st1, bins=binedges)
    
    if save:
        set_analysis_data(pqf, 'photon_histogram', np.vstack((h0,h1)).transpose(),
                          columns=('channel_0', 'channel_1'))
        set_analysis_data(pqf, 'photon_histogram_binedges_ns', b0)
        delete_analysis_data(pqf, 'photon_histogram_event_filter')
        if fltr != None:
            set_analysis_data(pqf, 'photon_histogram_event_filter', fltr)
        
    return (h0, b0), (h1, b1)


def get_photon_hists_from_folder(folder, **kw):
    '''
    return the cumulative photon histogram from all data contained in a folder
    (all sub-levels are searched).
    '''
    filepaths = tb.get_all_msmt_filepaths(folder)
    for i,f in enumerate(filepaths):
        if i == 0:
            pqf = pqf_from_fp(f, rights = 'r+')
            (h0,b0),(h1,b1) = get_photon_hist(pqf, **kw)
        else:
            pqf = pqf_from_fp(f, rights = 'r+')
            (_h0,_b0),(_h1,_b1) = get_photon_hist(pqf, **kw)
            h0 += _h0
            h1 += _h1
    return (h0, b0), (h1, b1)


##############################################################################
### Plotting photon histograms
##############################################################################

def _plot_photon_hist(ax, h, b, log=True, **kw):
    label = kw.pop('label', '')

    _h = h.astype(float)
    _h[_h<=1e-1] = 1e-1
    _h = np.append(_h, _h[-1])
           
    ax.plot(b, _h, drawstyle='steps-post', label=label)
    if log:
        ax.set_yscale('log')
    ax.set_xlabel('time (ps)')
    ax.set_ylabel('events')
    ax.set_ylim(bottom=0.1)
    ax.set_xlim(min(b), max(b))

def plot_photon_hist(pqf, **kw):    
    ret = kw.pop('ret', 'subplots')

    (h0, b0), (h1, b1) = get_photon_hist(pqf, **kw)
   
    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,8))
    _plot_photon_hist(ax0, h0, b0)
    _plot_photon_hist(ax1, h1, b1)

    ax0.set_title('photons channel 0')
    ax1.set_title('photons channepql 1')

    fp = fp_from_pqf(pqf)
    
    fig.suptitle(tb.get_msmt_header(fp) + ' -- Photon histogram')
    
    if ret == 'subplots':
        return fig, (ax0, ax1)


def plot_photon_hist_folder(folder, **kw):    
    ret = kw.pop('ret', 'subplots')

    (h0, b0), (h1, b1) = get_photon_hists_from_folder(folder, **kw)
   
    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,8))
    _plot_photon_hist(ax0, h0, b0)
    _plot_photon_hist(ax1, h1, b1)

    ax0.set_title('photons channel 0')
    ax1.set_title('photons channel 1')

    # fp = fp_from_pqf(pqf)
    
    # fig.suptitle(tb.get_msmt_header(fp) + ' -- Photon histogram')
    
    if ret == 'subplots':
        return fig, (ax0, ax1)        
