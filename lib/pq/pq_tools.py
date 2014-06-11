import numpy as np

def get_photons(pqf):
    """
    returns two filters (1d-arrays): whether events are ch0-photons/ch1-photons
    """
    channel = pqf['/PQ_channel-1'].value
    special = pqf['/PQ_special-1'].value
    
    is_not_special = special==0
    is_channel_0 = channel==0
    is_channel_1 = channel==1
    
    is_photon_0 = np.logical_and(is_not_special, is_channel_0)
    is_photon_1 = np.logical_and(is_not_special, is_channel_1)
    
    return is_photon_0, is_photon_1

def filter_on_same_sync_number(source_sync_numbers, target_sync_numbers):
    """
    returns a filter for target_sync_numbers that's true for all sync numbers that are also
    in source_sync_numbers.
    """
    return np.in1d(target_sync_numbers, source_sync_numbers)

def get_markers(pqf, chan):
    """
    returns a filter (1d-array): whether events are markers on the given channel
    """
    channel = pqf['/PQ_channel-1'].value
    special = pqf['/PQ_special-1'].value
    
    is_special = special==1
    is_channel = channel==chan

    return (is_special & is_channel)

def get_multiple_photon_syncs(pqf):
    special = pqf['/PQ_special-1'].value
    sync_numbers = pqf['/PQ_sync_number-1'].value

    is_photon = special == 0
    photon_sync_numbers = sync_numbers[is_photon]
    
    # this works nicely on sorted arrays
    is_multiple_photon_sync = photon_sync_numbers[1:] == photon_sync_numbers[:-1]
    #multiple_photon_sync_numbers = photon_sync_numbers[is_multiple_photon_sync]

    return is_multiple_photon_sync

def get_coincidences(pqf, fltr0=None, fltr1=None):

    sync_time = pqf['/PQ_sync_time-1'].value
    total_time = pqf['/PQ_time-1'].value
    sync_number = pqf['/PQ_sync_number-1'].value
    
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
    t1 = total_time[fltr1]
    sn1 = sync_number[fltr1]
    #print len(st0),len(t0),len(sn0),len(st1),len(t1),len(sn1),

    samesync0 = np.in1d(sn0, sn1)
    samesync1 = np.in1d(sn1, sn0)
    
    c_st0 = st0[samesync0]
    c_st1 = st1[samesync1]
    c_t0 = t0[samesync0]
    c_sn0 = sn0[samesync0]
    c_t1 = t1[samesync1]
    c_sn1 = sn1[samesync1]
    
    coincidences = np.empty((0,4))
    for _sn0, _t0, _st0 in zip(c_sn0, c_t0, c_st0):
        _c = c_sn1==_sn0
        
        for _t1, _st1 in zip(c_t1[_c], c_st1[_c]):
            dt = int(_t0) - int(_t1)
            coincidences = np.vstack((coincidences, np.array([dt, _st0, _st1, _sn0])))
                       
    return coincidences

def filter_synctimes(pqf, t0, t1, window_reps=1, window_period=None):
    sync_time = pqf['/PQ_sync_time-1'].value
    for r in range(window_reps):
        if r == 0:
            fltr = (sync_time >= t0) & (sync_time <= t1)
        else:
            fltr = fltr | (sync_time >= (r * window_period + t0)) & (sync_time <= (r * window_period + t1))
    return fltr