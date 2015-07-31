import os
import datetime
import numpy as np
import h5py
import analysis.lib.tools.toolbox as tb
from analysis.lib.bell import bell_events as be
import shutil

def get_latest_analysis_fp(folder, pattern='total_events'):
    fns_srt=np.sort(os.listdir(folder))[::-1]
    for fn in fns_srt:
        if (pattern in fn) and (fn[-5:]=='.hdf5'):
            return os.path.join(folder,fn)

def get_lt_fps(fps_bs, lt3_folder, lt4_folder):
    max_measurement_delay = datetime.timedelta(minutes=2,seconds=10)
    fps_lt3 = []
    fps_lt4 = []
    for fp_bs in fps_bs:
        bs_m_folder = os.path.split(fp_bs)[0] 
        bs_m_time = tb.get_datetime_from_folder(bs_m_folder)
        bs_m_name = tb.get_measurement_name_from_folder(bs_m_folder)[8:]
        min_timestamp = tb.timestamp_from_datetime(bs_m_time - max_measurement_delay)
        max_timestamp = tb.timestamp_from_datetime(bs_m_time + max_measurement_delay)
        
        lt3_m_folder = tb.latest_data(contains = bs_m_name, folder = lt3_folder, older_than = max_timestamp, newer_than = min_timestamp)
        fps_lt3.append(tb.get_msmt_fp(lt3_m_folder))
        lt4_m_folder = tb.latest_data(contains = bs_m_name, folder = lt4_folder, older_than = max_timestamp, newer_than = min_timestamp)
        fps_lt4.append(tb.get_msmt_fp(lt4_m_folder))

    return fps_lt3, fps_lt4

def collect_lt_data(bs_folder, lt3_source_folder, lt4_source_folder, measurement_pattern):
    fps_bs = tb.get_all_msmt_filepaths(bs_folder, pattern=measurement_pattern)
    print 'Found {} filepaths'.format(len(fps_bs))
    target_folder = os.path.split(bs_folder)[0]
    fps_lt3, fps_lt4 = get_lt_fps(fps_bs, lt3_source_folder, lt4_source_folder)
    skipped_fps=[]
    for target_name, fps in zip(['LT3','LT4'],[fps_lt3,fps_lt4]):
        print 'copying {} files'.format(target_name)    
        for i,fp in enumerate(fps):
            print i,
            src_folder = os.path.split(fp)[0]
            folder,dst_folder = os.path.split(src_folder)
            date_folder = os.path.split(folder)[1]
            try:
                shutil.copytree(src_folder,os.path.join(target_folder,target_name,date_folder,dst_folder))
            except WindowsError:
                print 'Error copying fp {}'.format(fp)
                skipped_fps.append(fp)

    return skipped_fps

def process_bell_data(bs_folder, lt3_folder, lt4_folder, measurement_pattern, bs_params, lt_params,
                       analysis_fp=None, update_previous_analysis_fp=None, ignore_unequal_markers=False, 
                       process_lt3=True, process_lt4=True, VERBOSE=False):
    fps_bs = tb.get_all_msmt_filepaths(bs_folder, pattern=measurement_pattern)
    print 'Found {} filepaths'.format(len(fps_bs))
    
    if update_previous_analysis_fp!=None:
        f = h5py.File(update_previous_analysis_fp,'r')
        total_ent_events     = f['analysis']['total_ent_events'].value
        total_lt3_ssro       = f['analysis']['total_lt3_ssro'].value
        total_lt4_ssro       = f['analysis']['total_lt4_ssro'].value
        total_ent_events_fps = f['analysis']['total_ent_events_fps'].value
        total_lt3_ssro_fps   = f['analysis']['total_lt3_ssro_fps'].value
        total_lt4_ssro_fps   = f['analysis']['total_lt4_ssro_fps'].value
        fps_bs = np.setdiff1d(fps_bs,total_ent_events_fps)
        f.close()
        print '{} filepaths left to analyze'.format(len(fps_bs))
    else:
        total_ent_events     = np.empty((0,be._bs_noof_columns), dtype=np.uint64)
        total_lt3_ssro       = np.empty((0,be._lt_noof_columns), dtype=np.uint64)
        total_lt4_ssro       = np.empty((0,be._lt_noof_columns), dtype=np.uint64)
        total_ent_events_fps = []
        total_lt3_ssro_fps   = []
        total_lt4_ssro_fps   = []

    fps_lt3, fps_lt4 = get_lt_fps(fps_bs, lt3_folder, lt4_folder)

    p_bs = bs_params
    p_lt = lt_params

    for i,fp_bs,fp_lt3,fp_lt4 in zip(range(len(fps_bs)),fps_bs,fps_lt3,fps_lt4):
        print i,

        ent_event_list = be.get_entanglement_event_list(fp_bs,
                                                        st_start_ch0   = p_bs['st_start_ch0'],   st_start_ch1 = p_bs['st_start_ch1'], st_len = p_bs['st_len'], pulse_sep = p_bs['pulse_sep'],
                                                        st_pulse_start_ch0 = p_bs['st_pulse_start_ch0'],st_pulse_start_ch1 = p_bs['st_pulse_start_ch1'], st_pulse_len = p_bs['st_pulse_len'], pulse_max_sn_diff = p_bs['pulse_max_sn_diff'],
                                                        ent_marker_channel_bs = p_bs['ent_marker_channel_bs'],
                                                        VERBOSE=VERBOSE)

        lt3_ssro_list  = be.get_ssro_result_list(fp_lt3,
                                                 ro_start  = p_lt['ro_start'],  ro_length  = p_lt['ro_length'],  ro_channel  = p_lt['ro_channel'],
                                                 rnd_start = p_lt['rnd_start'], rnd_length = p_lt['rnd_length'], rnd_channel = p_lt['rnd_channel'], 
                                                 rnd_0_channel = p_lt['rnd_0_channel'], rnd_1_channel = p_lt['rnd_1_channel'],
                                                 psb_tail_start = p_lt['psb_tail_start_lt3'], psb_tail_len = p_lt['psb_tail_len'], pulse_sep = p_bs['pulse_sep'], 
                                                 ent_marker_channel_lt = p_lt['ent_marker_channel_lt3'], ent_marker_lt_timebin_limit = p_lt['ent_marker_lt_timebin_limit'],
                                                 sn_diff_marker_ent_early = p_lt['sn_diff_marker_ent_early'], sn_diff_marker_ent_late = p_lt['sn_diff_marker_ent_late'],
                                                 invalid_marker_channel_lt = p_lt['invalid_marker_channel_lt'],invalid_marker_max_sn_diff = p_lt['invalid_marker_max_sn_diff'],
                                                 VERBOSE = VERBOSE) if process_lt3 else np.zeros((len(ent_event_list),be._lt_noof_columns))
        #lt3_ssro_list  = be.get_ssro_result_list_adwin(fp_lt3, ssro_result_list=lt3_ssro_list)
        lt4_ssro_list  = be.get_ssro_result_list(fp_lt4,
                                                 ro_start  = p_lt['ro_start'],  ro_length  = p_lt['ro_length'],  ro_channel  = p_lt['ro_channel'],
                                                 rnd_start = p_lt['rnd_start'], rnd_length = p_lt['rnd_length'], rnd_channel = p_lt['rnd_channel'], 
                                                 rnd_0_channel = p_lt['rnd_0_channel'], rnd_1_channel = p_lt['rnd_1_channel'],
                                                 psb_tail_start = p_lt['psb_tail_start_lt4'], psb_tail_len = p_lt['psb_tail_len'], pulse_sep = p_bs['pulse_sep'],
                                                 ent_marker_channel_lt = p_lt['ent_marker_channel_lt4'], ent_marker_lt_timebin_limit = p_lt['ent_marker_lt_timebin_limit'],
                                                 sn_diff_marker_ent_early = p_lt['sn_diff_marker_ent_early'], sn_diff_marker_ent_late = p_lt['sn_diff_marker_ent_late'],
                                                 invalid_marker_channel_lt = p_lt['invalid_marker_channel_lt'], invalid_marker_max_sn_diff = p_lt['invalid_marker_max_sn_diff'],
                                                 VERBOSE = VERBOSE) if process_lt4 else np.zeros((len(ent_event_list),be._lt_noof_columns))
        #lt4_ssro_list = be.get_ssro_result_list_adwin(fp_lt4, ssro_result_list=lt4_ssro_list)
        
        if (len(ent_event_list) != len(lt3_ssro_list)) or (len(ent_event_list) != len(lt4_ssro_list)):
            print 'WARNING: measurement with filepath {}: Number of markers is unequal'.format(fp_bs)
            print 'BS markers: {}, LT3 markers: {}, LT4 markers: {}'.format(len(ent_event_list),len(lt3_ssro_list),len(lt4_ssro_list))
            minlen=min((len(ent_event_list),len(lt3_ssro_list),len(lt4_ssro_list)))
            maxlen=max((len(ent_event_list),len(lt3_ssro_list),len(lt4_ssro_list)))
            if VERBOSE:
                print ent_event_list[:,be._cl_sn], lt3_ssro_list[:,be._cl_sn_ma], lt4_ssro_list[:,be._cl_sn_ma]
            if not ignore_unequal_markers:
                print 'File ignored'
                continue
            elif ignore_unequal_markers == 'fix_last':
                print 'trying to fix unequal markers by discarding last'
                
                if maxlen - minlen >1:
                    print 'Fix failed, number of markers too different: {}!'.format(maxlen - minlen)
                    continue
                if minlen == 0:
                    continue
                ent_event_list = ent_event_list[:minlen,:]
                lt3_ssro_list = lt3_ssro_list[:minlen,:]
                lt4_ssro_list = lt4_ssro_list[:minlen,:]
                if abs(int(ent_event_list[-1,be._cl_sn])-int(lt4_ssro_list[-1,be._cl_sn_ma]))>250 or \
                    abs(int(ent_event_list[-1,be._cl_sn])-int(lt3_ssro_list[-1,be._cl_sn_ma]/250.*251))>250:
                    print 'Fix failed, last marker sync number too different:\n \
                                            BS: {}, LT3 (/250*251): {}. LT4: {}'.format(ent_event_list[-1,be._cl_sn],
                                                                                        int(lt3_ssro_list[-1,be._cl_sn_ma]/250.*251),
                                                                                        lt4_ssro_list[-1,be._cl_sn_ma])
                print 'Fix suceeded'
            elif ignore_unequal_markers == 'append_zeros':
                print 'Appending empty events for missing markers'
                ent_event_list=  np.vstack((ent_event_list, np.zeros((maxlen-len(ent_event_list),be._bs_noof_columns), dtype=be._bs_dtype)))
                lt3_ssro_list =  np.vstack((lt3_ssro_list, np.zeros((maxlen-len(lt3_ssro_list),be._lt_noof_columns), dtype=be._lt_dtype)))
                lt4_ssro_list =  np.vstack((lt4_ssro_list, np.zeros((maxlen-len(lt4_ssro_list),be._lt_noof_columns), dtype=be._lt_dtype)))

        if len(ent_event_list) > 50:
            print 'Measurement with filepath {} has more than 100 events'.format(fp_bs)
        total_ent_events = np.vstack((total_ent_events, ent_event_list))
        total_lt3_ssro   = np.vstack((total_lt3_ssro, lt3_ssro_list))
        total_lt4_ssro   = np.vstack((total_lt4_ssro, lt4_ssro_list))
        for j in range(len(ent_event_list)):
            total_ent_events_fps = np.append(total_ent_events_fps,fp_bs)
            total_lt3_ssro_fps   = np.append(total_lt3_ssro_fps,fp_lt3)
            total_lt4_ssro_fps   = np.append(total_lt4_ssro_fps,fp_lt4)
        
    print 'Done, total_events:',len(total_ent_events)

    if analysis_fp!=None:
        tb.set_analysis_data(analysis_fp,'total_ent_events',     data=total_ent_events,     attributes=[], permissions='a')
        tb.set_analysis_data(analysis_fp,'total_lt3_ssro',       data=total_lt3_ssro,       attributes=[])
        tb.set_analysis_data(analysis_fp,'total_lt4_ssro',       data=total_lt4_ssro,       attributes=[])
        tb.set_analysis_data(analysis_fp,'total_ent_events_fps', data=total_ent_events_fps, attributes=[])
        tb.set_analysis_data(analysis_fp,'total_lt3_ssro_fps',   data=total_lt3_ssro_fps,   attributes=[])
        tb.set_analysis_data(analysis_fp,'total_lt4_ssro_fps',   data=total_lt4_ssro_fps,   attributes=[])

    return total_ent_events, total_lt3_ssro, total_lt4_ssro, total_ent_events_fps, total_lt3_ssro_fps, total_lt3_ssro_fps, total_lt4_ssro_fps

def process_tpqi_data(fps_bs, bs_params, analysis_fp=None, update_previous_analysis_fp=None):
    p_bs = bs_params

    if update_previous_analysis_fp!=None:
        f = h5py.File(update_previous_analysis_fp,'r')
        all_coincidences = f['analysis']['tpqi'].value
        total_ent_events_fps = f['analysis']['total_ent_events_fps'].value
        fps_bs = np.setdiff1d(fps_bs,total_ent_events_fps)
        f.close()
    else:
        all_coincidences = np.empty((0,be._tpqi_noof_columns), dtype=be._tpqi_dtype)

    print 'Processing TPQI {} files'.format(len(fps_bs))
    
    for i,fp in enumerate(fps_bs):
        print i,
        coincidences = be.get_coincidences(fp,st_start_ch0 = p_bs['st_start_ch0'], st_start_ch1 = p_bs['st_start_ch1'], st_len = p_bs['st_len'], pulse_sep = p_bs['pulse_sep'])
        all_coincidences = np.vstack((all_coincidences, coincidences))
    
    if analysis_fp!=None:
        tb.set_analysis_data(analysis_fp, 'tpqi', data=all_coincidences, attributes=[], permissions='a')
    print 'Done!'
    return all_coincidences

def process_lt_stats(fps_lt, lt_params, lt3, analysis_fp=None,  update_previous_analysis_fp=None):

    p_lt = lt_params
    psb_tail_start = p_lt['psb_tail_start_lt3'] if lt3 else p_lt['psb_tail_start_lt4']
    psb_prepulse_start = p_lt['psb_prepulse_start_lt3'] if lt3 else p_lt['psb_prepulse_start_lt4']

    if update_previous_analysis_fp!=None:
        f = h5py.File(update_previous_analysis_fp,'r')
        all_lt_stats = f['analysis']['lt3_stats'].value if lt3 else f['analysis']['lt4_stats'].value
        total_lt_fps = f['analysis']['total_lt3_ssro_fps'].value if lt3 else  f['analysis']['total_lt4_ssro_fps'].value
        noofs_fps = len(fps_lt)
        fps_lt = np.setdiff1d(fps_lt,total_lt_fps)
        noof_old_fps = noofs_fps - len(fps_lt)
        all_lt_stats = np.resize(all_lt_stats,(noofs_fps, be._lt_stats_noof_columns))
        f.close()
    else:
        all_lt_stats = np.zeros((len(fps_lt),be._lt_stats_noof_columns), dtype=be._lt_stats_dtype)
        noof_old_fps = 0

    print 'Processing LT {} files'.format(len(fps_lt))   

    for i,fp in enumerate(fps_lt):
        print i,
        all_lt_stats[i+noof_old_fps,:]=  be.get_lt_stats(fp, ro_start  = p_lt['ro_start'],  ro_length  = p_lt['ro_length'],  ro_channel  = p_lt['ro_channel'],
                                       rnd_start = p_lt['rnd_start'], rnd_length = p_lt['rnd_length'], 
                                       rnd_0_channel = p_lt['rnd_0_channel'], rnd_1_channel = p_lt['rnd_1_channel'], 
                                       psb_tail_start = psb_tail_start, psb_tail_len = p_lt['psb_tail_len'], pulse_sep = p_lt['pulse_sep'],
                                       psb_prepulse_start = psb_prepulse_start, psb_prepulse_len =  p_lt['psb_prepulse_len'])

    if analysis_fp!=None:
        tb.set_analysis_data(analysis_fp, 'lt3_stats' if lt3 else 'lt4_stats', data=all_lt_stats, attributes=[], permissions='a')
    print 'Done!'
    return all_lt_stats

def process_bs_hist_stats(fps_bs, bs_params, analysis_fp=None, update_previous_analysis_fp=None):
    p_bs = bs_params

    if update_previous_analysis_fp!=None:
        f = h5py.File(update_previous_analysis_fp,'r')
        hist_stats = f['analysis']['bs_hist_stats'].value
        total_hist = f['analysis']['bs_total_hist'].value.T 
        total_ent_events_fps = f['analysis']['total_ent_events_fps'].value
        noofs_fps = len(fps_bs)
        fps_bs = np.setdiff1d(fps_bs,total_ent_events_fps)
        noof_old_fps = noofs_fps - len(fps_bs)
        hist_stats = np.resize(hist_stats,(noofs_fps, be._bs_hist_stats_noof_columns))
        f.close()
    else:
        hist_stats = np.zeros((len(fps_bs), be._bs_hist_stats_noof_columns), dtype=be._bs_hist_stats_dtype)
        total_hist = None
        noof_old_fps = 0

    print 'Processing tail {} files'.format(len(fps_bs))

    for i,fp in enumerate(fps_bs):
        print i,
        hi_0, hi_1, hist_stats[i+noof_old_fps,:] = be.get_bs_hist_stats(fp, st_start_ch0   = p_bs['st_start_ch0'],   st_start_ch1 = p_bs['st_start_ch1'], st_len = p_bs['st_len'], pulse_sep = p_bs['pulse_sep'],
                                              st_pulse_start = p_bs['st_pulse_start'], st_pulse_len = p_bs['st_pulse_len'], hist_binsize_ps =p_bs['hist_binsize_ps'] )
        if total_hist == None:
            total_hist = np.vstack((hi_0,hi_1))
        else: 
            total_hist += np.vstack((hi_0,hi_1))
    
    if analysis_fp!=None:
        tb.set_analysis_data(analysis_fp, 'bs_hist_stats', data=hist_stats, attributes=[], permissions='a')
        tb.set_analysis_data(analysis_fp, 'bs_total_hist', data=total_hist.T, attributes=[], permissions='a')
    print 'Done!'
    return total_hist.T,hist_stats

def process_afterpulsing_data(fps_bs, afterpulsing_params, analysis_fp=None, update_previous_analysis_fp=None):
    p_ap = afterpulsing_params
    
    if update_previous_analysis_fp!=None:
        f = h5py.File(update_previous_analysis_fp,'r')
        all_afterpulsing = f['analysis']['afterpulsing'].value
        total_ent_events_fps = f['analysis']['total_ent_events_fps'].value
        noofs_fps = len(fps_bs)
        fps_bs = np.setdiff1d(fps_bs,total_ent_events_fps)
        noof_old_fps = noofs_fps - len(fps_bs)
        all_afterpulsing = np.resize(all_afterpulsing,(noofs_fps, be._bs_afterpulsing_noof_columns))
        f.close()
    else:
        all_afterpulsing = np.zeros((len(fps_bs),be._bs_afterpulsing_noof_columns), dtype=be._bs_afterpulsing_dtype)
        noof_old_fps = 0
    
    print 'Processing afterpulsing {} files'.format(len(fps_bs))

    for i,fp in enumerate(fps_bs):
        print i,
        all_afterpulsing[i] += be.get_bs_afterpulsing(fp, first_st_start = p_ap['first_st_start'], first_st_len = p_ap['first_st_len'], after_st_start = p_ap['after_pulse_st_start'], after_st_len = p_ap['after_pulse_st_len'])
    
    if analysis_fp!=None:
        tb.set_analysis_data(analysis_fp, 'afterpulsing', data=all_afterpulsing, attributes=[], permissions='a')
    print 'Done!'
    return all_afterpulsing

def get_unique_bell_fps_from_analysis_file(analysis_fp):
    f = h5py.File(analysis_fp,'r')
    bs_fps  = np.unique(f['analysis']['total_ent_events_fps'].value)
    lt3_fps = np.unique(f['analysis']['total_lt3_ssro_fps'].value)
    lt4_fps = np.unique(f['analysis']['total_lt4_ssro_fps'].value)
    f.close()
    return bs_fps, lt3_fps, lt4_fps