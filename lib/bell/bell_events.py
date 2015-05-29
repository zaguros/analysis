import numpy as np
import h5py
import os
from analysis.lib.m2 import m2

#cl is short for column
_bs_noof_columns = 11
_cl_sn          = 0 
_cl_type        = 1 #type of entanglement event: 1 is valid double click event, 2 is w1 click event, 3 is w2 click event, 4 is invalid event
_cl_ch_w1       = 2 
_cl_ch_w2       = 3 
_cl_st_w1       = 4 
_cl_st_w2       = 5 
_cl_tt_w1       = 6
_cl_tt_w2       = 7
_cl_psi_min_ma  = 8
_cl_psi_plus_ma = 9
_cl_pulse_cts   = 10

_lt_noof_columns = 16
_cl_sn_ma       = 0
_cl_sn_ro       = 1
_cl_noof_ph_ro  = 2
_cl_st_ma       = 3
_cl_noof_rnd    = 4
_cl_noof_rnd_0  = 5
_cl_noof_rnd_1  = 6
_cl_cr_after    = 7
_cl_ro_after    = 8
_cl_inv_mrkr    = 9
_cl_tt_ma = 10
_cl_tt_rnd = 11
_cl_first_ph_st = 12
_cl_noof_ph_tail = 13
_cl_noof_rnd_0_prev  = 14
_cl_noof_rnd_1_prev  = 15
#_cl_noof_ph_ro_prev  = 16

def elements_unique(arr):
    return np.array_equal(arr, np.unique(arr))


def get_ssro_result_list_adwin(fp_lt, ssro_result_list=None):
    a= m2.M2Analysis(os.path.split(fp_lt)[0])
    ssro_g = a.g['ssro']
    ro_data = ssro_g['RO_data'].value
    cr_after = ssro_g['CR_after'].value
    if ssro_result_list == None:
        ssro_result_list=np.zeros((len(ro_data),_lt_noof_columns), dtype=np.uint64)
        ssro_result_list[:,_cl_noof_ph_ro]  = ro_data
    elif len(ro_data) == len(ssro_result_list):
        ssro_g_result_list[:,_cl_ro_after]  = ro_data
    else:
        print 'WARNING: {} has adwin RO length {:d}, but number of TH readouts is {:d} not match'.format(fp_lt,len(ro_data), len(ssro_result_list))
        print 'ignoring adwin file'
        return ssro_result_list
    ssro_result_list[:,_cl_cr_after]  = cr_after
    a.finish()
    return ssro_result_list

def get_ssro_result_list(fp_lt,
                            ro_start, ro_length, ro_channel,
                            rnd_start, rnd_length, rnd_channel, rnd_0_channel, rnd_1_channel,
                            psb_tail_start,psb_tail_len,
                            ent_marker_channel_lt, ent_marker_lt_timebin_limit, sn_diff_marker_ent_early, sn_diff_marker_ent_late,
                            invalid_marker_channel_lt, invalid_marker_max_sn_diff,
                            VERBOSE=False):

    pqf_lt  = h5py.File(fp_lt,  'r')
    i=1
    while ('PQ_special-'+str(i)) in pqf_lt.keys():
        sp_lt = pqf_lt['/PQ_special-'+str(i)].value      
        ch_lt = pqf_lt['/PQ_channel-'+str(i)].value
        sn_lt = pqf_lt['/PQ_sync_number-'+str(i)].value
        st_lt = pqf_lt['/PQ_sync_time-'+str(i)].value
        tt_lt = pqf_lt['/PQ_time-'+str(i)].value 
        
        is_ent_marker =     ((sp_lt==1) & (ch_lt & ent_marker_channel_lt     == ent_marker_channel_lt))
        is_invalid_marker = ((sp_lt==1) & (ch_lt & invalid_marker_channel_lt == invalid_marker_channel_lt))
        marker_sn = sn_lt[is_ent_marker]
        marker_or_ent_sn      = np.hstack((marker_sn,marker_sn-1,marker_sn-2))
        marker_or_ent_sn_fltr = np.in1d(sn_lt,marker_or_ent_sn)
        
        invalid_marker_sn = sn_lt[is_invalid_marker] if i==1 else np.hstack((invalid_marker_sn,sn_lt[is_invalid_marker]))
        marker_sns = marker_sn if i==1 else np.hstack((marker_sns,marker_sn))
        sp = sp_lt[marker_or_ent_sn_fltr] if i==1 else np.hstack((sp,sp_lt[marker_or_ent_sn_fltr]))
        ch = ch_lt[marker_or_ent_sn_fltr] if i==1 else np.hstack((ch,ch_lt[marker_or_ent_sn_fltr]))
        sn = sn_lt[marker_or_ent_sn_fltr] if i==1 else np.hstack((sn,sn_lt[marker_or_ent_sn_fltr]))
        st = st_lt[marker_or_ent_sn_fltr] if i==1 else np.hstack((st,st_lt[marker_or_ent_sn_fltr]))
        tt = tt_lt[marker_or_ent_sn_fltr] if i==1 else np.hstack((tt,tt_lt[marker_or_ent_sn_fltr]))
        
        i+=1
    pqf_lt.close()

    if not elements_unique(marker_sns):
        raise Exception('File {} has multiple entanglement markers in one sync'.format(fp_lt))
    if VERBOSE:
        print 'Found {} entanglement markers in {}-{}'.format(len(marker_sns),fp_lt, i)

    #print 'made subset of data'
    fltr_ro         = (sp == 0) & (ch == ro_channel)  & (st > ro_start)  & (st < (ro_start  + ro_length))
    fltr_rnd        = (sp == 0) & (ch == rnd_channel) & (st > rnd_start) & (st < (rnd_start + rnd_length))
    fltr_tail       = (sp == 0) & (ch == ro_channel)  & (st > psb_tail_start)  & (st < (psb_tail_start  + psb_tail_len))
    fltr_rnd0       = (sp == 1) & (ch == rnd_0_channel)
    fltr_rnd1       = (sp == 1) & (ch == rnd_1_channel)
    fltr_ent_marker = (sp == 1) & (ch == ent_marker_channel_lt)

    ssro_result_list=np.zeros((len(marker_sns),_lt_noof_columns), dtype=np.uint64)
    for i,cur_sn in enumerate(marker_sns):

        is_cur_ent_marker = fltr_ent_marker & (sn == cur_sn)
        cur_ent_marker_st = st[is_cur_ent_marker][0]
        cur_ent_marker_tt = tt[is_cur_ent_marker][0]
        if cur_ent_marker_st <= ent_marker_lt_timebin_limit:
            ent_sn = cur_sn+sn_diff_marker_ent_early
        elif cur_ent_marker_st > ent_marker_lt_timebin_limit:
            ent_sn = cur_sn+sn_diff_marker_ent_late
        else:
            raise Exception('Unable to compare marker synctime for sync {} in file {}'.format(sn,fp_lt))
        
        fltr_ent =( sn == ent_sn)
        filter_lt4 = ( sn == ent_sn-1)
        fltr_ro_photon = fltr_ent & fltr_ro
        fltr_tail_photon = fltr_ent & fltr_tail
        fltr_rnd_click = fltr_ent & fltr_rnd
        fltr_rnd0_ev = fltr_ent & fltr_rnd0
        fltr_rnd1_ev = fltr_ent & fltr_rnd1
        diff_invalid_ev = np.array(ent_sn - invalid_marker_sn.astype(np.int64), dtype=np.int64)
        fltr_rnd0_prev = ( sn == ent_sn-1) & fltr_rnd0
        fltr_rnd1_prev = ( sn == ent_sn-1) & fltr_rnd1
        fltr_ro_photon_prev = ( sn == ent_sn-1) & fltr_ro
        
        if np.sum(fltr_rnd0_ev | fltr_rnd1_ev )==1:
            rnd_ma_tt =  tt[fltr_rnd0_ev | fltr_rnd1_ev]
        else:
            rnd_ma_tt = 0

        if np.sum(fltr_ro_photon)>0:
            first_ph_st = st[fltr_ro_photon][0]
        else:
            first_ph_st = 0

        
        fltr_invalid_ev = (diff_invalid_ev > 0) & (diff_invalid_ev <= invalid_marker_max_sn_diff)

        ssro_result_list[i,_cl_sn_ma]       = cur_sn
        ssro_result_list[i,_cl_sn_ro]       = ent_sn
        ssro_result_list[i,_cl_noof_ph_ro]  = np.sum(fltr_ro_photon)
        ssro_result_list[i,_cl_st_ma]       = cur_ent_marker_st
        ssro_result_list[i,_cl_noof_rnd]    = np.sum(fltr_rnd_click)
        ssro_result_list[i,_cl_noof_rnd_0]  = np.sum(fltr_rnd0_ev)
        ssro_result_list[i,_cl_noof_rnd_1]  = np.sum(fltr_rnd1_ev)
        ssro_result_list[i,_cl_inv_mrkr]    = np.sum(fltr_invalid_ev)
        ssro_result_list[i,_cl_tt_ma]       = cur_ent_marker_tt
        ssro_result_list[i,_cl_tt_rnd]      = rnd_ma_tt
        ssro_result_list[i,_cl_first_ph_st] = first_ph_st
        ssro_result_list[i,_cl_noof_ph_tail]= np.sum(fltr_tail_photon)
        ssro_result_list[i,_cl_noof_rnd_0_prev]= np.sum(fltr_rnd0_prev)
        ssro_result_list[i,_cl_noof_rnd_1_prev]= np.sum(fltr_rnd1_prev)
        #ssro_result_list[i,_cl_noof_ph_ro_prev]= np.sum(fltr_ro_photon_prev)
        

    return ssro_result_list


def get_entanglement_event_list(fp_bs,
                                st_start_ch0, st_start_ch1, st_len, pulse_sep,
                                st_pulse_start, st_pulse_len, pulse_max_sn_diff,
                                ent_marker_channel_bs,psi_min_marker_bs, psi_plus_marker_bs,
                                VERBOSE=False):
    
    pqf_bs  = h5py.File(fp_bs,  'r')
    sp_bs = pqf_bs['/PQ_special-1'].value      
    ch_bs = pqf_bs['/PQ_channel-1'].value
    sn_bs = pqf_bs['/PQ_sync_number-1'].value
    st_bs = pqf_bs['/PQ_sync_time-1'].value
    tt_bs = pqf_bs['/PQ_time-1'].value 
    pqf_bs.close()
    
    is_ent_marker=((sp_bs==1) & ((ch_bs & ent_marker_channel_bs)==ent_marker_channel_bs))
    ent_sn = sn_bs[is_ent_marker]
    if not elements_unique(ent_sn):
        raise Exception('File {} has multiple entanglement markers in one sync'.format(fp_bs))
    if VERBOSE:
        print 'Found {} entanglement markers in {}'.format(len(ent_sn),fp_bs)

    marker_sn_fltr = np.in1d(sn_bs,ent_sn)
    pulse_st_fltr = (st_pulse_start <= st_bs ) & (st_bs < st_pulse_start+st_pulse_len) & (sp_bs == 0)
    sn = sn_bs[marker_sn_fltr]
    st = st_bs[marker_sn_fltr]
    ch = ch_bs[marker_sn_fltr]
    sp = sp_bs[marker_sn_fltr]
    tt = tt_bs[marker_sn_fltr]
    #print 'made subset of data'
    fltr_w1_ch0 = (((st_start_ch0 <= st)  & (st < (st_start_ch0 + st_len))) & (ch == 0) & (sp == 0))
    fltr_w1_ch1 = (((st_start_ch1 <= st)  & (st < (st_start_ch1 + st_len))) & (ch == 1) & (sp == 0)) 
    fltr_w2_ch0 = (((st_start_ch0 +  pulse_sep <= st) & (st < (st_start_ch0 + pulse_sep + st_len))) & (ch == 0) & (sp == 0)) 
    fltr_w2_ch1 = (((st_start_ch1 +  pulse_sep <= st) & (st < (st_start_ch1 + pulse_sep + st_len))) & (ch == 1) & (sp == 0)) 
    fltr_w1 = fltr_w1_ch0 | fltr_w1_ch1
    fltr_w2 = fltr_w2_ch0 | fltr_w2_ch1

    fltr_psi_min_marker  = ((sp==1) & ((ch & psi_min_marker_bs)==psi_min_marker_bs))
    fltr_psi_plus_marker = ((sp==1) & ((ch & psi_plus_marker_bs)==psi_plus_marker_bs))

    ent_event_list=np.zeros((len(ent_sn),_bs_noof_columns),dtype=np.uint64)
    for i,cur_sn in enumerate(ent_sn):
        fltr_ent = (sn==cur_sn)
        
        pulse_sn_diffs = (cur_sn - sn_bs.astype(np.int64)) 
        pulse_sn_fltr = (pulse_sn_diffs > 0) & (pulse_sn_diffs<3000000)#3 million ~ 60 secs
        ent_event_list[i,_cl_pulse_cts]= np.sum(pulse_sn_fltr & pulse_st_fltr)


        ent_event_list[i,_cl_sn]    = cur_sn   
        ent_event_list[i,_cl_psi_min_ma] = np.sum(fltr_ent & fltr_psi_min_marker)
        ent_event_list[i,_cl_psi_plus_ma] = np.sum(fltr_ent & fltr_psi_plus_marker)

        if (np.sum(fltr_ent & fltr_w1) == 1) and (np.sum(fltr_ent & fltr_w2) == 1):
            #a valid entanglement evetn: one photon in window 1, and one in window 2
            ent_event_list[i,_cl_type]  = 1
            ent_event_list[i,_cl_ch_w1] = ch[fltr_ent & fltr_w1]
            ent_event_list[i,_cl_ch_w2] = ch[fltr_ent & fltr_w2]
            ent_event_list[i,_cl_st_w1] = st[fltr_ent & fltr_w1]
            ent_event_list[i,_cl_st_w2] = st[fltr_ent & fltr_w2]
            ent_event_list[i,_cl_tt_w1] = tt[fltr_ent & fltr_w1]
            ent_event_list[i,_cl_tt_w2] = tt[fltr_ent & fltr_w2]

        elif (np.sum(fltr_ent & fltr_w1) == 0) and (np.sum(fltr_ent & fltr_w2) == 1):
            #a valid spin-photon entanglement event: zero photon in window 1, and one in window 2
            ent_event_list[i,_cl_type]  = 2
            ent_event_list[i,_cl_ch_w1] = 0
            ent_event_list[i,_cl_ch_w2] = ch[fltr_ent & fltr_w2]
            ent_event_list[i,_cl_st_w1] = 0
            ent_event_list[i,_cl_st_w2] = st[fltr_ent & fltr_w2]
            ent_event_list[i,_cl_tt_w1] = 0
            ent_event_list[i,_cl_tt_w2] = tt[fltr_ent & fltr_w2]

        elif (np.sum(fltr_ent & fltr_w1) == 1) and (np.sum(fltr_ent & fltr_w2) == 0):
            #a valid spin-photon entanglement event: one photon in window 1, and zero in window 2
            ent_event_list[i,_cl_type]  = 3
            ent_event_list[i,_cl_ch_w1] = ch[fltr_ent & fltr_w1]
            ent_event_list[i,_cl_ch_w2] = 0
            ent_event_list[i,_cl_st_w1] = st[fltr_ent & fltr_w1]
            ent_event_list[i,_cl_st_w2] = 0
            ent_event_list[i,_cl_tt_w1] = tt[fltr_ent & fltr_w1]
            ent_event_list[i,_cl_tt_w2] = 0

        else:
            #anything else is not valid.
            ent_event_list[i,_cl_type]  = 4
            ent_event_list[i,_cl_ch_w1] = 0
            ent_event_list[i,_cl_ch_w2] = 0
            ent_event_list[i,_cl_st_w1] = 0
            ent_event_list[i,_cl_st_w2] = 0
            ent_event_list[i,_cl_tt_w1] = 0
            ent_event_list[i,_cl_tt_w2] = 0 
            
    if VERBOSE:
        print 'Of which {} are valid entanglement events'.format(np.sum(ent_event_list[:,1]==1))
        print 'and of which {} are valid spin-photon events'.format(np.sum(ent_event_list[:,1]==2)+np.sum(ent_event_list[:,1]==3))
    return ent_event_list
    
def get_latest_analysis_fp(folder, pattern='total_events'):
    fns_srt=np.sort(os.listdir(folder))[::-1]
    for fn in fns_srt:
        if (pattern in fn) and (fn[-5:]=='.hdf5'):
            return os.path.join(folder,fn)