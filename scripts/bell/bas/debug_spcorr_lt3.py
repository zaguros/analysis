pqf=h5py.File(fn,'r')
reps=19929*250

pq_binsize_ns=1

RO_length = 4000
RO_start = 10700

fig, ax = plt.subplots(1,1, figsize=(4.5,4))

sync_nrs=pqf['/PQ_sync_number-1'].value 
reps=len(np.unique(sync_nrs))
is_photon_0, is_rnd_clk=pq_tools.get_photons(pqf)
sync_time_ns = pqf['/PQ_sync_time-1'].value * pq_binsize_ns

is_marker_1_event=pq_tools.get_markers(pqf,1)
is_marker_2_event=pq_tools.get_markers(pqf,2)
noof_rnd_0_events=np.sum(is_marker_1_event)
noof_rnd_1_events=np.sum(is_marker_2_event)
print 'noof_rnd 0/1 events:',noof_rnd_0_events, '/' , noof_rnd_1_events
print 'bias toward 0 : {:.2f} % '.format(50-float(noof_rnd_0_events)/(noof_rnd_0_events+noof_rnd_1_events)*100),', error : {:.2f} %'.format(1/np.sqrt(len(np.where(is_marker_1_event)[0])+len(np.where(is_marker_2_event)[0]))*100)
print 'noof syncs:', sync_nrs[-1]
print 'Detected marker events {} / {}:'.format(noof_rnd_0_events+noof_rnd_1_events, reps)
noof_reps_wo_rnd_clk=len(np.unique(sync_nrs[is_rnd_clk]))
print 'number of reps with a random clk', noof_reps_wo_rnd_clk
print 'syncs without a random click: {} / {} = {:.2f} %'.format(reps-noof_reps_wo_rnd_clk, reps, float(reps-noof_reps_wo_rnd_clk)/reps*100.)

is_last_random_click=np.append(np.diff(np.asarray(is_rnd_clk, dtype='int'))==-1,is_rnd_clk[-1])
start_rnd=np.min(sync_time_ns[is_rnd_clk])-20
length_rnd=np.max(sync_time_ns[is_rnd_clk])-start_rnd+20
pq_plots.plot_photon_hist_filter_comparison(pqf,is_last_random_click,start = start_rnd, length = length_rnd, hist_binsize = 1, save = False)

marker_1_sync_numbers= sync_nrs[is_marker_1_event]
marker_2_sync_numbers= sync_nrs[is_marker_2_event]

st_fltr = (RO_start  <= sync_time_ns) &  (sync_time_ns< (RO_start + RO_length))
is_photon_0_in_ro_window = st_fltr & is_photon_0
photon_in_0_ro_window_sync_numbers = sync_nrs[is_photon_0_in_ro_window]
no_photon_in_0_ro_window_sync_numbers = np.setdiff1d(sync_nrs,photon_in_0_ro_window_sync_numbers)
av_p0=float(len(np.unique(photon_in_0_ro_window_sync_numbers)))/reps
u_av_p0 = np.sqrt(av_p0*(1-av_p0)/reps)
    

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
plt.close('all')
pqf.close()