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
    returns a filter (1d-array): whether events are markers on the given channel
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
            files.set_analysis_data(fp, 'photon_histogram_event_filter', fltr)
        
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
    is_mrkr = get_markers(fp, chan)
        
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

##############################################################################
### Teleportation events
##############################################################################

def get_teleportation_events(folder, force_eval=False, verbose=True):
    fp = files.get_msmt_fp(folder)
    fname = files.get_msmt_name(fp)

    if files.has_analysis_data(fp, 'LDE_ph0_sync_times') and \
        files.has_analysis_data(fp, 'LDE_ph1_sync_times') and \
        files.has_analysis_data(fp, 'LDE_ph0_sync_numbers') and \
        files.has_analysis_data(fp, 'LDE_ph1_sync_numbers') and \
        files.has_analysis_data(fp, 'PLU_mrkr_sync_numbers') and \
        files.has_analysis_data(fp, 'psi_mrkr_sync_numbers') and \
        not force_eval:

        return True

    f = h5py.File(fp, 'r')
    sync_times = f['/HH_sync_time-1'].value * 1e-3 # we prefer ns over ps
    sync_numbers = f['/HH_sync_number-1'].value
    channels = f['/HH_channel-1'].value
    ad1_reps = f['{}/adwin_lt1_data/completed_reps'.format(fname)].value
    ad1_ssro1 = f['{}/adwin_lt1_data/SSRO1_results'.format(fname)].value
    ad1_ssro2 = f['{}/adwin_lt1_data/SSRO2_results'.format(fname)].value
    ad2_reps = f['{}/adwin_lt2_data/completed_reps'.format(fname)].value
    ad2_ssro = f['{}/adwin_lt2_data/SSRO_lt2_data'.format(fname)].value
    f.close()

    ### First the HH events

    is_PLU_mrkr = get_markers(fp, 2)
    has_same_sync_number_as_PLU_mrkr = filter_marker(fp, 2)
    is_psiplus_mrkr = get_markers(fp, 3)
    is_psiminus_mrkr = get_markers(fp, 4)
    has_same_sync_number_as_psiminus_mrkr = filter_marker(fp, 4)


    is_ph_ch0, is_ph_ch1 = get_photons(fp)
    is_ph = is_ph_ch0 | is_ph_ch1


    is_ph_with_PLU_mrkr = is_ph & filter_marker(fp, 2)
    is_ph0_with_PLU_mrkr = is_ph_ch0 & filter_marker(fp, 2)
    is_ph1_with_PLU_mrkr = is_ph_ch1 & filter_marker(fp, 2)
    is_ph0_with_psiminus_mrkr = is_ph_ch0 & filter_marker(fp, 4)
    is_ph1_with_psiminus_mrkr = is_ph_ch1 & filter_marker(fp, 4)

    print
    print fp
    print
    print 'HH data'
    print '-------'
    print 'Number of PLU markers received:', len(np.where(is_PLU_mrkr)[0])
    print 'Number of Psi+ markers received:', len(np.where(is_psiplus_mrkr)[0])
    print 'Number of Psi- markers received:', len(np.where(is_psiminus_mrkr)[0])
    print 'Number of photons with the same sync as a PLU marker:', len(np.where(is_ph_with_PLU_mrkr)[0])
    print '    ch0:', len(np.where(is_ph0_with_PLU_mrkr)[0])
    print '    ch1:', len(np.where(is_ph1_with_PLU_mrkr)[0])

    is_ph0_with_PLU_mrkr_in_first_window = is_ph0_with_PLU_mrkr & \
        filter_synctime(sync_times, settings.CH0_START,settings.CH0_STOP)
    is_ph1_with_PLU_mrkr_in_first_window = is_ph1_with_PLU_mrkr & \
        filter_synctime(sync_times, settings.CH1_START, settings.CH1_STOP)
    is_ph_with_PLU_mrkr_in_first_window = is_ph0_with_PLU_mrkr_in_first_window | \
        is_ph1_with_PLU_mrkr_in_first_window

    print
    print 'Number of PLU-marked photons in first window:', \
        len(np.where(is_ph_with_PLU_mrkr_in_first_window)[0])

    is_ph0_with_PLU_mrkr_in_second_window = is_ph0_with_PLU_mrkr & filter_synctime(sync_times, 
        settings.CH0_START+settings.PIPULSESEP, 
        settings.CH0_STOP+settings.PIPULSESEP)
    is_ph1_with_PLU_mrkr_in_second_window = is_ph1_with_PLU_mrkr & filter_synctime(sync_times, 
        settings.CH1_START+settings.PIPULSESEP, 
        settings.CH1_STOP+settings.PIPULSESEP)
    is_ph_with_PLU_mrkr_in_second_window = is_ph0_with_PLU_mrkr_in_second_window | \
        is_ph1_with_PLU_mrkr_in_second_window

    print 'Number of PLU-marked photons in second window:', \
        len(np.where(is_ph_with_PLU_mrkr_in_second_window)[0])

    LDE_ph0_sync_times = sync_times[is_ph0_with_PLU_mrkr]
    LDE_ph0_sync_numbers = sync_numbers[is_ph0_with_PLU_mrkr]
    LDE_ph1_sync_times = sync_times[is_ph1_with_PLU_mrkr]
    LDE_ph1_sync_numbers = sync_numbers[is_ph1_with_PLU_mrkr]
    PLU_mrkr_sync_numbers = sync_numbers[is_PLU_mrkr]
    psi_mrkr_sync_numbers = sync_numbers[is_psiminus_mrkr]

    ### Adwin SSROs

    print 
    print 'Adwin LT1'
    print '---------'
    print 'Number of events:', ad1_reps,
    if ad1_reps != len(np.where(is_PLU_mrkr)[0]):
        raise Exception(fp + '\nAdwin LT1 events do not match HH events - data set seems faulty :(')
    else:
        print 'OK :)'
        
    print 
    print 'Adwin LT2'
    print '---------'
    print 'Number of events:', ad2_reps,
    if ad2_reps != len(np.where(is_PLU_mrkr)[0]):
        raise Exception(fp+'\nAdwin LT2 events do not match HH events - data set seems faulty :(')
    else:
        print 'OK :)'
        
    files.set_analysis_data(fp, 'LDE_ph0_sync_times', LDE_ph0_sync_times)
    files.set_analysis_data(fp, 'LDE_ph0_sync_numbers', LDE_ph0_sync_numbers)
    files.set_analysis_data(fp, 'LDE_ph1_sync_times', LDE_ph1_sync_times)
    files.set_analysis_data(fp, 'LDE_ph1_sync_numbers', LDE_ph1_sync_numbers)
    files.set_analysis_data(fp, 'PLU_mrkr_sync_numbers', PLU_mrkr_sync_numbers)
    files.set_analysis_data(fp, 'psi_mrkr_sync_numbers', psi_mrkr_sync_numbers)
    files.set_analysis_data(fp, 'adwin1_readouts', ad1_reps)
    files.set_analysis_data(fp, 'adwin2_readouts', ad2_reps)

    return True

def filter_LDE_events(folder, **kw):
    ch0_start = kw.pop('ch0_start', settings.CH0_START)
    ch1_start = kw.pop('ch1_start', settings.CH1_START)
    window_length = kw.pop('window_length', settings.WINDOW_LENGTH)
    dt_max = kw.pop('dt_max', settings.DT_MAX)

    ch0_stop = ch0_start + window_length
    ch1_stop = ch1_start + window_length

    fp = files.get_msmt_fp(folder)
    fname = files.get_msmt_name(fp)

    if not files.has_analysis_data(fp, 'LDE_ph0_sync_times') or \
        not files.has_analysis_data(fp, 'LDE_ph1_sync_times') or \
        not files.has_analysis_data(fp, 'LDE_ph0_sync_numbers') or \
        not files.has_analysis_data(fp, 'LDE_ph1_sync_numbers') or \
        not files.has_analysis_data(fp, 'PLU_mrkr_sync_numbers') or \
        not files.has_analysis_data(fp, 'psi_mrkr_sync_numbers'):

        raise Exception(fp+'\nLDE events not extracted yet! Use the get_teleportation_events method for that.')

    LDE_ph0_sync_times, _a = files.get_analysis_data(fp, 'LDE_ph0_sync_times')
    LDE_ph1_sync_times, _a = files.get_analysis_data(fp, 'LDE_ph1_sync_times')
    LDE_ph0_sync_numbers, _a = files.get_analysis_data(fp, 'LDE_ph0_sync_numbers')
    LDE_ph1_sync_numbers, _a = files.get_analysis_data(fp, 'LDE_ph1_sync_numbers')
    PLU_mrkr_sync_numbers, _a = files.get_analysis_data(fp, 'PLU_mrkr_sync_numbers')
    psi_mrkr_sync_numbers, _a = files.get_analysis_data(fp, 'psi_mrkr_sync_numbers')

    ph0_w1_filter = filter_synctime(LDE_ph0_sync_times,
                                    ch0_start, ch0_stop)
    ph1_w1_filter = filter_synctime(LDE_ph1_sync_times,
                                    ch1_start, ch1_stop)
    ph0_w2_filter = filter_synctime(LDE_ph0_sync_times,
                                    ch0_start + settings.PIPULSESEP, 
                                    ch0_stop + settings.PIPULSESEP)
    ph1_w2_filter = filter_synctime(LDE_ph1_sync_times,
                                    ch1_start + settings.PIPULSESEP, 
                                    ch1_stop + settings.PIPULSESEP)

    filters = [ph0_w1_filter, ph1_w1_filter, ph0_w2_filter, ph1_w2_filter]
    photon_sync_numbers = [LDE_ph0_sync_numbers, LDE_ph1_sync_numbers, LDE_ph0_sync_numbers, LDE_ph1_sync_numbers]
    photon_sync_times = [LDE_ph0_sync_times, LDE_ph1_sync_times, LDE_ph0_sync_times, LDE_ph1_sync_times]

    # encode events as follows: [ ch0_w1, ch1_w1, ch0_w2, ch1_w2, psi- ]; 0 = False, 1 = True
    LDE_events = np.zeros((len(PLU_mrkr_sync_numbers), 5))
    dts = np.ones(len(PLU_mrkr_sync_numbers)) * (-1)
    
    for i,s in enumerate(PLU_mrkr_sync_numbers):
        ts = np.array([])
        for j,f,p,t in zip(range(4), filters, photon_sync_numbers, photon_sync_times):
            LDE_events[i,j] += len(p[f][p[f]==s])
            if len(p[f][p[f]==s]) == 1:
                ts = np.append(ts, t[f][p[f]==s])
        
        if len(ts) == 2:
            ts[np.argmax(ts)] -= settings.PIPULSESEP
            dts[i] = abs(ts[0]-ts[1])

        LDE_events[i,4] += len(psi_mrkr_sync_numbers[psi_mrkr_sync_numbers==s])

    is_psiminus_event = (((LDE_events[:,0] > 0) & (LDE_events[:,3] > 0)) | \
                        ((LDE_events[:,1] > 0) & (LDE_events[:,2] > 0))) & \
                        ((dts > 0) & (dts < dt_max))
    is_psiminus_event_with_psi_mrkr = is_psiminus_event & (LDE_events[:,4] > 0)
    is_psiplus_event = (((LDE_events[:,0] > 0) & (LDE_events[:,2] > 0)) | \
                       ((LDE_events[:,1] > 0) & (LDE_events[:,3] > 0))) & \
                       ((dts > 0) & (dts < dt_max))

    print
    print fp
    print '----'
    print 'Psi- events: ', len(np.where(is_psiminus_event)[0])
    print 'Psi- events (incl mrkr): ', len(np.where(is_psiminus_event_with_psi_mrkr)[0])
    print 'Psi+ events: ', len(np.where(is_psiplus_event)[0])

    return is_psiminus_event_with_psi_mrkr, is_psiplus_event

def get_adwin_readouts(folder):
    fp = files.get_msmt_fp(folder)
    fname = files.get_msmt_name(fp)

    f = h5py.File(fp, 'r')
    ad1_reps = f['{}/adwin_lt1_data/completed_reps'.format(fname)].value
    ad1_ssro1 = f['{}/adwin_lt1_data/SSRO1_results'.format(fname)].value
    ad1_ssro2 = f['{}/adwin_lt1_data/SSRO2_results'.format(fname)].value
    ad2_reps = f['{}/adwin_lt2_data/completed_reps'.format(fname)].value
    ad2_ssro = f['{}/adwin_lt2_data/SSRO_lt2_data'.format(fname)].value
    f.close()

    BSM_outcomes = (2**1 * ad1_ssro2) + (2**0 * ad1_ssro1) # binary encoding; i just made it painfully obvious :)

    return BSM_outcomes, ad2_ssro

def get_readouts(folder, **kw):
    fp = files.get_msmt_fp(folder)
    fname = files.get_msmt_name(fp)

    is_psiminus_event, is_psiplus_event = filter_LDE_events(folder)

    f = h5py.File(fp, 'r')
    ad1_reps = f['{}/adwin_lt1_data/completed_reps'.format(fname)].value
    ad1_ssro1 = f['{}/adwin_lt1_data/SSRO1_results'.format(fname)].value
    ad1_ssro2 = f['{}/adwin_lt1_data/SSRO2_results'.format(fname)].value
    ad2_reps = f['{}/adwin_lt2_data/completed_reps'.format(fname)].value
    ad2_ssro = f['{}/adwin_lt2_data/SSRO_lt2_data'.format(fname)].value
    f.close()

    # note that for the BSM outcome the first qubit is the nitrogen
    BSM_outcome_names = ['00', '01', '10', '11']
    BSM_outcomes = (2**1 * ad1_ssro2) + (2**0 * ad1_ssro1) # binary encoding; i just made it painfully obvious :)
    Bob_outcomes_psiminus = np.empty((0,2))
    Bob_outcomes_psiplus = np.empty((0,2))
    BSM_outcomes_psiminus = np.array([])
    BSM_outcomes_psiplus = np.array([])

    for outcome in range(4):
        ssro_results_psiminus = ad2_ssro[((BSM_outcomes==outcome) & is_psiminus_event)]   
        BSM_outcomes_psiminus = np.append(BSM_outcomes_psiminus,
            len(ssro_results_psiminus))
        Bob_outcomes_psiminus = np.vstack((Bob_outcomes_psiminus,
            np.array([len(ssro_results_psiminus[ssro_results_psiminus==0]),
                      len(ssro_results_psiminus[ssro_results_psiminus==1])])))

        ssro_results_psiplus = ad2_ssro[((BSM_outcomes==outcome) & is_psiplus_event)]   
        BSM_outcomes_psiplus = np.append(BSM_outcomes_psiplus,
            len(ssro_results_psiplus))
        Bob_outcomes_psiplus = np.vstack((Bob_outcomes_psiplus,
            np.array([len(ssro_results_psiplus[ssro_results_psiplus==0]),
                      len(ssro_results_psiplus[ssro_results_psiplus==1])])))

    return BSM_outcomes_psiminus, BSM_outcomes_psiplus, \
        Bob_outcomes_psiminus, Bob_outcomes_psiplus


    





