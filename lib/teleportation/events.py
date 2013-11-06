import numpy as np
import h5py

from analysis.lib.fitting import fit
from analysis.lib import tools

import files, stats, settings

### identification of events
def get_photons(fp):
    """
    returns two filters (1d-arrays): whether events are ch0-photons/ch1-photons
    """
    f = h5py.File(fp, 'r')
    channel = f['/HH_channel-1'].value
    special = f['/HH_special-1'].value
    f.close()

    is_not_special = special==0
    is_channel_0 = channel==0
    is_channel_1 = channel==1

    is_photon_0 = np.logical_and(is_not_special, is_channel_0)
    is_photon_1 = np.logical_and(is_not_special, is_channel_1)

    return is_photon_0, is_photon_1

def get_markers(fp, chan):
    """
    returns a filter (1-array): whether events are markers on the given channel
    """
    f = h5py.File(fp, 'r')
    channel = f['/HH_channel-1'].value
    special = f['/HH_special-1'].value
    f.close()
    
    is_special = special==1
    is_channel = channel==chan
    
    return (is_special & is_channel)

### photon histograms
def get_photon_hist(fp, **kw):
    save = kw.pop('save', False)
    fltr = kw.pop('fltr', None)
    force_eval = kw.pop('force_eval', True)
    binedges = kw.pop('binedges', settings.PHOTONHIST_BINEDGES)
    
    if not force_eval and files.has_analysis_data(fp, 'photon_histogram'):
        h, h_attrs = files.get_analysis_data(fp, 'photon_histogram')
        be, be_attrs = files.get_analysis_data(fp, 'photon_histogram_binedges_ns')
        h0 = h[:,0]
        h1 = h[:,1]
        return (h0, be), (h1, be)
    
    f = h5py.File(fp, 'r')
    sync_time_ns = f['/HH_sync_time-1'].value * 1e-3
    f.close()
    
    ph0, ph1 = get_photons(fp)
    if fltr != None:
        _fltr0 = (ph0 & fltr)
        _fltr1 = (ph1 & fltr)
    else:
        _fltr0 = ph0
        _fltr1 = ph1
    
    st0 = sync_time_ns[_fltr0]
    st1 = sync_time_ns[_fltr1]
    
    h0, b0 = np.histogram(st0, bins=binedges)
    h1, b1 = np.histogram(st1, bins=binedges)
    
    if save:
        files.set_analysis_data(fp, 'photon_histogram', vstack((h0,h1)).transpose(),
                          columns=('channel_0', 'channel_1'))
        files.set_analysis_data(fp, 'photon_histogram_binedges_ns', b0)
        files.delete_analysis_data(fp, 'photon_histogram_event_filter')
        if fltr != None:
            files.set_analysis(fp, 'photon_histogram_event_filter', fltr)
        
    return (h0, b0), (h1, b1)


def get_photon_hists_from_folder(folder, **kw):
    '''
    return the cumulative photon histogram from all data contained in a folder
    (all sub-levels are searched).
    '''
    filepaths = files.get_all_msmt_filepaths(folder)
    for i,f in enumerate(filepaths):
        if i == 0:
            (h0,b0),(h1,b1) = get_photon_hist(f, **kw)
        else:
            (_h0,_b0),(_h1,_b1) = get_photon_hist(f, **kw)
            h0 += _h0
            h1 += _h1
    return (h0, b0), (h1, b1)

##############################################################################
### filtering
##############################################################################

def filter_synctime(sync_times, t0, t1):
    """
    return a filter that's true for all sync times that are in between t0 and t1
    (t0 is the lower bound, t1 is the upper bound)
    """
    fltr = ((sync_times >= t0) & (sync_times <= t1))
    return fltr

def filter_on_same_sync_number(source_sync_numbers, target_sync_numbers):
    """
    returns a filter for target_sync_numbers that's true for all sync numbers that are also
    in source_sync_numbers.
    """
    return np.in1d(target_sync_numbers, source_sync_numbers)

def filter_marker(fp, chan):
    """
    Note: at the moment this filter includes the marker events on which we filter.
    """
    is_mrkr = events.get_markers(fp, chan)
        
    f = h5py.File(fp, 'r')
    sync_numbers = f['/HH_sync_number-1'].value
    f.close()
    
    marker_sync_numbers = sync_numbers[is_mrkr]
    
    return filter_on_same_sync_number(marker_sync_numbers, sync_numbers)


##############################################################################
### Tail analysis
##############################################################################

def _fit_tail(h, b):
    binwidth = b[1] - b[0]
    _b = b-b[0]
    
    A = fit.Parameter(h[0], 'A')
    T = fit.Parameter(11, 'T')
    o = fit.Parameter(0, 'o')
    p0 = [A, T, o]
    fitfunc_str = 'A exp(-t/T) + o'
    
    def ff(x):
        return o() + A() * np.exp(-x/T())
    fitres = fit.fit1d(_b, h, None, p0=p0, fitfunc=ff, fitfunc_str=fitfunc_str,
                    do_print=True, ret=True)
    return fitres

def fit_tail(fp, **kw):
    (h0, b0), (h1, b1) = get_photon_hist(fp, **kw)
        
    
    _x0 = (b0[:-1]+b0[1:])/2.
    _x1 = (b1[:-1]+b1[1:])/2.
    binwidth = _x0[1] - _x0[0]
    tail0 = (_x0 >= settings.CH0_START) & (_x0 <= settings.CH0_START + settings.FIT_TAIL_LENGTH)
    tail1 = (_x1 >= settings.CH1_START) & (_x1 <= settings.CH1_START + settings.FIT_TAIL_LENGTH)
    x0, y0 = _x0[tail0], h0[tail0]
    x1, y1 = _x1[tail1], h1[tail1]
    
    fitres0 = _fit_tail(y0, x0)
    fitres1 = _fit_tail(y1, x1)
    
    starts = stats.get_sequence_starts(fp) * settings.SEQREPS
    tcpsh0 = fitres0['params_dict']['A'] * fitres0['params_dict']['T'] / binwidth / starts
    tcpsh1 = fitres1['params_dict']['A'] * fitres1['params_dict']['T'] / binwidth / starts
    
    _x0 = np.linspace(x0[0], x0[-1], 201)
    ret0 = { 'params' : fitres0['params_dict'], 
            'errors' : fitres0['error_dict'], 
            'data_x' : x0,
            'data_y' : y0, 
            'fit_x' : _x0,
            'fit_y' : fitres0['fitfunc'](_x0-_x0[0]),
            'info' : fitres0['fitfunc_str'] + '\n' + fit.str_fit_params(fitres0),
            'tcpsh' : tcpsh0,
            }
    
    _x1 = np.linspace(x1[0], x1[-1], 201)
    ret1 = { 'params' : fitres1['params_dict'], 
            'errors' : fitres1['error_dict'], 
            'data_x' : x1,
            'data_y' : y1,
            'fit_x' : _x1,
            'fit_y' : fitres1['fitfunc'](_x1-_x1[0]),
            'info' : fitres1['fitfunc_str'] + '\n' + fit.str_fit_params(fitres1),
            'tcpsh' : tcpsh1,
            }
    
    keys0 = []
    params0 = np.empty((0,2))
    for k in ret0['params']:
        keys0.append(k)
        params0 = np.vstack((params0, np.array([ret0['params'][k], ret0['errors'][k]])))
    
    files.set_analysis_data(fp, 'ch0_fitparams', params0, subgroup='tail',
                      info=ret0['info'], columns=('fitvalue', 'error'),
                      rows=keys0)
    files.set_analysis_data(fp, 'ch0_fit_data', np.vstack((ret0['data_x'], ret0['data_y'])).transpose(),
                      subgroup='tail', columns=('time bin center (ns)', 'counts per bin'))
    files.set_analysis_data(fp, 'ch0_fit_curve', np.vstack((ret0['fit_x'], ret0['fit_y'])).transpose(),
                      subgroup='tail', columns=('time (ns)', 'counts per bin'))
    files.set_analysis_data(fp, 'ch0_tcpsh', ret0['tcpsh'], subgroup='tail',
                      starts = starts)
    
    keys1 = []
    params1 = np.empty((0,2))
    for k in ret1['params']:
        keys1.append(k)
        params1 = np.vstack((params1, np.array([ret1['params'][k], ret1['errors'][k]])))
    
    files.set_analysis_data(fp, 'ch1_fitparams', params1, subgroup='tail',
                      info=ret1['info'], columns=('fitvalue', 'error'),
                      rows=keys1)
    files.set_analysis_data(fp, 'ch1_fit_data', np.vstack((ret1['data_x'], ret1['data_y'])).transpose(),
                      subgroup='tail', columns=('time bin center (ns)', 'counts per bin'))
    files.set_analysis_data(fp, 'ch1_fit_curve', np.vstack((ret1['fit_x'], ret1['fit_y'])).transpose(),
                      subgroup='tail', columns=('time (ns)', 'counts per bin'))
    files.set_analysis_data(fp, 'ch1_tcpsh', ret1['tcpsh'], subgroup='tail',
                      starts = starts)
    
    return ret0, ret1