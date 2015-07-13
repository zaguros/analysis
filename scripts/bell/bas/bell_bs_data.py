#folder =r'D:\measuring\data\2015-04-14_SPCorrs\BS_lt4'
folder=r'D:\measuring\data\2015-04-21-lhfbt4'
#folder =r'D:\measuring\data\2015-03-03-lhfbt\BS' # r'D:\measuring\data\2014-12-Entanglement_XX_back_to_basics\BS'#
if load_from_analysis_file:
    afp = be.get_latest_analysis_fp(folder, pattern ='total_events.hdf5')
    af = h5py.File(afp,'r')
    fps=np.unique(af['analysis']['total_ent_events_fps'].value)
    af.close()
else:
    fps=tb.get_all_msmt_filepaths(folder, pattern='third_ever')
ttime=0
reps=0
if False:
    for fp in fps:
        try:
            f=h5py.File(fp,'r')
            ttime += f['/PQ_time-1'][-1]/1e12
            reps  += f['/PQ_sync_number-1'][-1]
        except ValueError:
            print fp
        finally:
            f.close()
    print 'Total measurement time: {:.2f} hours'.format(float(ttime)/60/60)
    print 'Total reps: {}'.format(str(reps))

#PQ_sync_number = np.empty((0,), dtype = np.uint32) 
#PQ_special = np.empty((0,), dtype = np.uint32)         
#PQ_sync_time = np.empty((0,), dtype = np.uint64)
#PQ_time = np.empty((0,), dtype = np.uint64)      
#PQ_channel = np.empty((0,), dtype = np.uint32)
st_start_ch0 = 5443500
st_len   = 50000 #50 ns
ch0_ch1_diff = 600 #1 ns
st_start_ch1=st_start_ch0 + ch0_ch1_diff
p_sep = 350000 #600 ns
tot_time=0
#load all BS data
if True:
    noof_markers=0
    psi_plus_0_sns=array([],dtype='uint32')
    psi_plus_1_sns=array([],dtype='uint32')
    psi_min_0_sns =array([],dtype='uint32')
    psi_min_1_sns =array([],dtype='uint32')
    marker_sns=array([],dtype='uint32')
    bins= np.linspace(5400,5600,2001)
    hists_0=np.empty((len(fps),len(bins)-1), dtype='int' )
    hists_1=np.empty((len(fps),len(bins)-1), dtype='int' )
    plu_misseds = np.zeros(len(fps))
    for i,fp in enumerate(fps):
        f=h5py.File(fp,'r')
        print 'Running{}/{}'.format(i,len(fps))
        if True:
            sn = f['/PQ_sync_number-1'].value 
            sp = f['/PQ_special-1'].value      
            st = f['/PQ_sync_time-1'].value
            #_PQ_time = f['/PQ_time'].value   
            ch = f['/PQ_channel-1'].value
        if False:
            try:
                tot_time += (f['/PQ_time-1'][-1] - f['/PQ_time-1'][0])/1e12
            except ValueError:
                continue
        f.close()
        if True:
            noof_markers+=np.sum(sp)
            marker_sn=sn[(sp==1) &(ch&1==1)]
            marker_sns=np.hstack((marker_sns,marker_sn))
            st_fltr_w1_ch0 = (((st_start_ch0<=st)  & (st<(st_start_ch0+st_len))) & (ch == 0) & (sp==0))
            st_fltr_w2_ch0 = (((st_start_ch0+p_sep<=st) & (st<(st_start_ch0+p_sep+st_len))) & (ch == 0) & (sp==0)) 
            st_fltr_w1_ch1 = (((st_start_ch1<=st)  & (st<(st_start_ch1+st_len)))  & (ch == 1) & (sp==0)) 
            st_fltr_w2_ch1 = (((st_start_ch1+p_sep<=st) & (st<(st_start_ch1+p_sep+st_len))) & (ch == 1) & (sp==0)) 

            w1_ch0_sn = sn[st_fltr_w1_ch0] 
            w2_ch0_sn = sn[st_fltr_w2_ch0] 
            w1_ch1_sn = sn[st_fltr_w1_ch1]
            w2_ch1_sn = sn[st_fltr_w2_ch1]
            psi_plus_0_sn=w1_ch0_sn[np.in1d(w1_ch0_sn,w2_ch0_sn)]
            psi_plus_1_sn=w1_ch1_sn[np.in1d(w1_ch1_sn,w2_ch1_sn)]
            psi_min_0_sn =w1_ch0_sn[np.in1d(w1_ch0_sn,w2_ch1_sn)]
            psi_min_1_sn =w1_ch1_sn[np.in1d(w1_ch1_sn,w2_ch0_sn)]
            all_ent=np.hstack((psi_plus_0_sn,psi_plus_1_sn,psi_min_0_sn,psi_min_1_sn))
            plu_misseds[i] = len(np.setdiff1d(all_ent,marker_sn))

            psi_plus_0_sns=np.hstack((psi_plus_0_sns,psi_plus_0_sn))
            psi_plus_1_sns=np.hstack((psi_plus_1_sns,psi_plus_1_sn))
            psi_min_0_sns =np.hstack((psi_min_0_sns, psi_min_0_sn))
            psi_min_1_sns =np.hstack((psi_min_1_sns, psi_min_1_sn))
        if False:
            hist_0,_tmp = np.histogram(st[(ch==0) & (sp==0)]/1000., bins=bins)
            hist_1,_tmp = np.histogram(st[(ch==1) & (sp==0)]/1000., bins=bins)
            hists_0[i]=hist_0
            hists_1[i]=hist_1

#compares the PLU entanglement markers, and the raw entanglement events
if True:
    all_ents=np.hstack((psi_min_0_sns,psi_min_1_sns,psi_plus_0_sns,psi_plus_1_sns))
    plu_missed = np.setdiff1d(all_ent,marker_sns)
    plu_wrong = np.setdiff1d(marker_sns,all_ent)

#semilogy(bins[:-1],np.sum(hists_0, axis=0))

#plot the position of the maximum count in the tail+pulse for each file.
if False:
    bins= np.linspace(5400,5600,2000)
    maxs=np.zeros(len(ch0))
    fn_ids=np.arange(len(fps))
    for i,y in enumerate(ch0):
        if np.sum(y)>1000:
            maxs[i] = bins[np.argmax(y)]

    plot(fn_ids[maxs>0], maxs[maxs>0]-5440)

#some histogram plotting
if False:
    rmin=370
    rmax=430
    for i in [8]:
        semilogy(bins[rmin:rmax],np.sum(ch0,axis=0)[rmin:rmax])

#tail counts, pulse counts and rejection
if False:
    bins= np.linspace(5400,5600,2000)
    tails=np.zeros(len(ch0))
    pulses=np.zeros(len(ch0))
    rmin=370
    rmax=430
    for i,y in enumerate(ch0):
        if np.sum(y)>1000:
            pulses[i]=np.sum(y[rmin:rmax])
            tails[i]=np.sum(y[rmax+20:])
    ax=plt.subplot(111)
    ax.plot(fn_ids[pulses>0],pulses[pulses>0])
    ax.plot(fn_ids[pulses>0],tails[pulses>0])
    ax2=ax.twinx()
    ax2.plot(fn_ids[pulses>0],pulses[pulses>0]/tails[pulses>0], color='r')