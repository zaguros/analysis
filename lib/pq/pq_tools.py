import numpy as np
from matplotlib import pyplot as plt
import h5py
from analysis.lib.tools import toolbox as tb

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

def get_markers(pqf, chan):
    """
    returns a filter (1d-array): whether events are markers on the given channel
    """
    
    if type(pqf) == h5py._hl.files.File:

        channel = pqf['/PQ_channel-1'].value
        special = pqf['/PQ_special-1'].value
    
        is_special = special==1
        is_channel = channel==chan

        return (is_special & is_channel)

    elif type(pqf) == str:

        f = h5py.File(pqf,'r')
        channel = f['/PQ_channel-1'].value
        special = f['/PQ_special-1'].value
        f.close()
    
        is_special = special==1
        is_channel = channel==chan

        return (is_special & is_channel)
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise

def get_multiple_photon_syncs(pqf):
    special = pqf['/PQ_special-1'].value
    sync_numbers = pqf['/PQ_sync_number-1'].value

    is_photon = special == 0
    photon_sync_numbers = sync_numbers[is_photon]
    
    # this works nicely on sorted arrays
    is_multiple_photon_sync = photon_sync_numbers[1:] == photon_sync_numbers[:-1]
    #multiple_photon_sync_numbers = photon_sync_numbers[is_multiple_photon_sync]

    return is_multiple_photon_sync

def get_coincidences(pqf, fltr0=None, fltr1=None, force_coincidence_evaluation = False, save = True):

    if has_analysis_data(pqf, 'coincidences') and not force_coincidence_evaluation:
        c, c_attrs = get_analysis_data(pqf, 'coincidences')
        return c    

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

    if save:
        set_analysis_data(pqf, 'coincidences', coincidences,
                          columns=('dt = ch0-ch1 (bins)', 'sync_time ch0 (bins)', 'sync_time ch1 (bins)', 'sync number'))
                       
    return coincidences

def get_coincidences_from_folder(folder):

    filepaths = tb.get_all_msmt_filepaths(folder) 
    co = np.ones([1,4])

    for i,f in enumerate(filepaths):
        if i == 0:
            pqf = pqf_from_fp(f, rights = 'r+')
            if 'PQ_sync_number-1' in pqf.keys():
                co = get_coincidences(pqf)           
        else:
            pqf = pqf_from_fp(f, rights = 'r+')
            if 'PQ_sync_number-1' in pqf.keys():
                if co[0,3] == 1:
                    co = get_coincidences(pqf)
                else:
                    co = np.vstack((co, get_coincidences(pqf)))
                    
    return co

def get_photons_in_sync_windows(pqf, first_win_min, first_win_max, second_win_min, second_win_max):
    """
    Returns two filters whether events are in the first or 
    in the second time window.
    """
    if type(pqf) == h5py._hl.files.File: 
        sync_time = pqf['/PQ_sync_time-1'].value
        special = pqf['/PQ_special-1'].value

    elif type(pqf) == str:
        f = h5py.File(pqf, 'r')
        sync_time = f['/PQ_sync_time-1'].value
        special = f['/PQ_special-1'].value
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
        
    if sum(is_photon_check) != sum(is_photon):
        print "Not all detected photons are in the broad windows set"
    
    return is_photon_first_window, is_photon_second_window




##############################################################################
### Filters
##############################################################################


def filter_synctimes(pqf, t0, t1, window_reps=1, window_period=None, pq_file = True):
    if pq_file:
        sync_time = pqf['/PQ_sync_time-1'].value
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

def filter_marker(pqf, chan):
    """
    Note: at the moment this filter includes the marker events on which we filter.
    """
    
    if type(pqf) == h5py._hl.files.File: 
        is_mrkr = get_markers(pqf, chan)

        sync_numbers = pqf['/PQ_sync_number-1'].value

        marker_sync_numbers = sync_numbers[is_mrkr]
    
        return filter_on_same_sync_number(marker_sync_numbers, sync_numbers)  

    elif type(pqf) == str:
        is_mrkr = get_markers(pqf, chan)
        
        f = h5py.File(pqf, 'r')
        sync_numbers = f['/PQ_sync_number-1'].value
        f.close()
    
        marker_sync_numbers = sync_numbers[is_mrkr]
    
        return filter_on_same_sync_number(marker_sync_numbers, sync_numbers)

    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise

def get_photons_with_markers(pqf, chan, first_win_min, first_win_max, second_win_min, second_win_max):
    """
    Return two filters whether events are first window photons or second window
    photons with markers.
    """
    
    is_photon_first_window, is_photon_second_window = get_photons_in_sync_windows(pqf,
                            first_win_min, first_win_max, second_win_min, second_win_max)
    is_events_with_marker = filter_marker(pqf,chan)
    
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


def get_photon_hist(pqf, **kw):
    save = kw.pop('save', False)
    fltr = kw.pop('fltr', None)
    force_eval = kw.pop('force_eval', True)
    start = kw.pop('start', 0) 
    length = kw.pop('length', 1e6) 
    hist_binsize= kw.pop('hist_binsize', 1e3) 
    
    if not force_eval and has_analysis_data(pqf, 'photon_histogram'):
        h, h_attrs = get_analysis_data(pqf, 'photon_histogram')
        be, be_attrs = get_analysis_data(pqf, 'photon_histogram_binedges')
        h0 = h[:,0]
        h1 = h[:,1]
        return (h0, be), (h1, be)
    
    sync_time = pqf['/PQ_sync_time-1'].value
    
    ph0, ph1 = get_photons(pqf)
    if fltr != None:
        _fltr0 = (ph0 & fltr)
        _fltr1 = (ph1 & fltr)
    else:
        _fltr0 = ph0
        _fltr1 = ph1
    
    st0 = sync_time[_fltr0]
    st1 = sync_time[_fltr1]

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
            (h0,b0),(h1,b1) = get_photon_hist(f, **kw)
        else:
            (_h0,_b0),(_h1,_b1) = get_photon_hist(f, **kw)
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