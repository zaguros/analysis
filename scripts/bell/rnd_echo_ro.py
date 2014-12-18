import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, sequence, pqsequence
from analysis.lib.m2 import m2
from analysis.lib.pq import pq_tools
from analysis.lib.tools import toolbox
from analysis.lib.math import error

from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.pq import pq_plots
#timestamp='20140521172533'
#folder=tb.data_from_time(timestamp)

def analyse_rnd_ro_bell(folder, save = True,RO_start=10730, **kw):
    a = pqsequence.PQSequenceAnalysis(folder)
    a.reps=a.g.attrs['repetitions']*a.g['joint_params'].attrs['LDE_attempts_before_CR']
    
    pq_binsize_ns=1
    
    RO_length = a.g['joint_params'].attrs['LDE_RO_duration']*1e9
    
    a.plot_histogram(0,start=0, length=3000, hist_binsize=1, save=False, log_plot=False)
    ssro_calib_folder = kw.pop('ssro_calib_folder', toolbox.latest_data('FastSSRO'))
    

    roc = error.SingleQubitROC()
    roc.F0, roc.u_F0, roc.F1, roc.u_F1 =pqsequence.get_analysed_fast_ssro_calibration(ssro_calib_folder, RO_length)
    #(0.9398,0.0034,0.9942,0.0013)
    print pqsequence.get_analysed_fast_ssro_calibration(ssro_calib_folder, RO_length)
    ##ssro.get_SSRO_calibration(ssro_calib_folder, ro_duration)
    
    fig, ax = plt.subplots(1,1, figsize=(4.5,4))
    
    sync_nrs=a.pqf['/PQ_sync_number-1'].value 
    is_marker_1_event=pq_tools.get_markers(a.pqf,1)
    is_marker_2_event=pq_tools.get_markers(a.pqf,2)
    noof_rnd_0_events=np.sum(is_marker_1_event)
    noof_rnd_1_events=np.sum(is_marker_2_event)
    print 'noof_rnd 0/1 events:',noof_rnd_0_events, '/' , noof_rnd_1_events
    print 'bias toward 0 : {:.2f} % '.format(50-float(noof_rnd_0_events)/(noof_rnd_0_events+noof_rnd_1_events)*100),', error : {:.2f} %'.format(1/np.sqrt(len(np.where(is_marker_1_event)[0])+len(np.where(is_marker_2_event)[0]))*100)
    print 'noof syncs:', sync_nrs[-1]
    print 'Detected marker events {} / {}:'.format(noof_rnd_0_events+noof_rnd_1_events, a.reps)
    
    is_photon_0, is_rnd_clk=pq_tools.get_photons(a.pqf)
    sync_time_ns = a.pqf['/PQ_sync_time-1'].value * pq_binsize_ns
    
    
    st_fltr = (RO_start  <= sync_time_ns) &  (sync_time_ns< (RO_start + RO_length))
    is_photon_0_in_ro_window = st_fltr & is_photon_0
    photon_in_0_ro_window_sync_numbers = sync_nrs[is_photon_0_in_ro_window]
    no_photon_in_0_ro_window_sync_numbers = np.setdiff1d(sync_nrs,photon_in_0_ro_window_sync_numbers)
    av_p0=float(len(np.unique(photon_in_0_ro_window_sync_numbers)))/a.reps
    u_av_p0 = np.sqrt(av_p0*(1-av_p0)/a.reps)
    av_F0, u_av_F0 = roc.num_eval(np.array([av_p0]),np.array([u_av_p0]))
    print 'averaged RO results: F0 {:.2f} $\pm$ {:.2f} % '.format(av_F0[0]*100,u_av_F0[0]*100 )
    
    
    noof_reps_wo_rnd_clk=len(np.unique(sync_nrs[is_rnd_clk]))
    print 'number of reps with a random clk', noof_reps_wo_rnd_clk
    print 'syncs without a random click: {} / {} = {:.2f} %'.format(a.reps-noof_reps_wo_rnd_clk, a.reps, float(a.reps-noof_reps_wo_rnd_clk)/a.reps*100.)
    
    is_last_random_click=np.append(np.diff(np.asarray(is_rnd_clk, dtype='int'))==-1,is_rnd_clk[-1])
    start_rnd=np.min(sync_time_ns[is_rnd_clk])-20
    length_rnd=np.max(sync_time_ns[is_rnd_clk])-start_rnd+20
    pq_plots.plot_photon_hist_filter_comparison(a.pqf,is_last_random_click,start = start_rnd, length = length_rnd, hist_binsize = 1, save = False)

    marker_1_sync_numbers= sync_nrs[is_marker_1_event]
    marker_2_sync_numbers= sync_nrs[is_marker_2_event]
    
    marker_1_ro_ms0_events=pq_tools.filter_on_same_sync_number(photon_in_0_ro_window_sync_numbers,marker_1_sync_numbers)
    marker_2_ro_ms0_events=pq_tools.filter_on_same_sync_number(photon_in_0_ro_window_sync_numbers,marker_2_sync_numbers)
    marker_1_ro_ms1_events=pq_tools.filter_on_same_sync_number(no_photon_in_0_ro_window_sync_numbers,marker_1_sync_numbers)#np.invert(marker_1_ro_ms0_events) #this also works.
    marker_2_ro_ms1_events=pq_tools.filter_on_same_sync_number(no_photon_in_0_ro_window_sync_numbers,marker_2_sync_numbers)#np.invert(marker_2_ro_ms0_events)
    noof_marker_1_ro_ms0_events=np.sum(marker_1_ro_ms0_events)
    noof_marker_2_ro_ms0_events=np.sum(marker_2_ro_ms0_events)
    noof_marker_1_ro_ms1_events=np.sum(marker_1_ro_ms1_events)
    noof_marker_2_ro_ms1_events=np.sum(marker_2_ro_ms1_events)
    
    print 'MA1 & RO0: {}, MA1 & RO1: {}, MA2 & RO0: {}, MA2 & RO1: {}'.format(noof_marker_1_ro_ms0_events, noof_marker_1_ro_ms1_events,noof_marker_2_ro_ms0_events, noof_marker_2_ro_ms1_events)
    
    
    
    ma_1_p0=(float(noof_marker_1_ro_ms0_events)/(noof_marker_1_ro_ms1_events+noof_marker_1_ro_ms0_events))
    ma_1_u_p0 = np.sqrt(ma_1_p0*(1-ma_1_p0)/(noof_marker_1_ro_ms1_events+noof_marker_1_ro_ms0_events))
    ma_2_p0=(float(noof_marker_2_ro_ms0_events)/(noof_marker_2_ro_ms1_events+noof_marker_2_ro_ms0_events))
    ma_2_u_p0 = np.sqrt(ma_2_p0*(1-ma_2_p0)/(noof_marker_2_ro_ms1_events+noof_marker_2_ro_ms0_events))        
    
    print 'Uncorrected: RND 0: F0 {:.2f}%, RND 1: F0 {:.2f}%'.format(ma_1_p0*100, ma_2_p0*100)

    p0, u_p0 = roc.num_eval(np.array([ma_1_p0,ma_2_p0]),np.array([ma_1_u_p0,ma_2_u_p0]))
    
    ax.bar( range(2),p0, 
        #color=[settings.COLORS[0], settings.COLORS[1]],
        align='center', yerr=u_p0, 
        ecolor='k', width=0.8)
    ax.set_title(a.timestamp+'\n Corrected RO')
    ax.text(0, -.15, 'Rnd_no = 0',ha='center', va='bottom')
    ax.text(1, -.15, 'Rnd_no = 1',ha='center', va='bottom')
    ax.set_xticks([0,1])
    #ax.text(1, 1.05, '{:.0f}+/-{:.0f} %'.format(p0*100., u_p0*100.),
    #    ha='center', va='bottom', color=settings.COLORS[1])  
    
    ax.text(0, 1.02,'F0: {:.2f} $\pm$ {:.2f} %'.format(p0[0]*100,u_p0[0]*100 ),ha='center', va='bottom')
    ax.text(1, 1.02,'F0: {:.2f} $\pm$ {:.2f} %'.format(p0[1]*100,u_p0[1]*100 ),ha='center', va='bottom')
    ax.set_ylabel('Fidelity ms0')
    ax.set_ylim(0,1.1)
    if save:
        a.save_fig_incremental_filename(fig,'random_mw_correlation_corrected')

if __name__ == '__main__':
    folder= toolbox.latest_data('Bell_RND')
    analyse_rnd_ro_bell(folder)