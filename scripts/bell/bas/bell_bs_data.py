folder = r'D:\bjhensen\data\Marked Bell Data\BS'



#PQ_sync_number = np.empty((0,), dtype = np.uint32) 
#PQ_special = np.empty((0,), dtype = np.uint32)         
#PQ_sync_time = np.empty((0,), dtype = np.uint64)
#PQ_time = np.empty((0,), dtype = np.uint64)      
#PQ_channel = np.empty((0,), dtype = np.uint32)
st_start_ch0 = 5442000
st_len   = 50000 #50 ns
ch0_ch1_diff = 1000 #1 ns
st_start_ch1=st_start_ch0 + ch0_ch1_diff
p_sep = 600000 #600 ns

#load all BS data
if False:
    fps=tb.get_all_msmt_filepaths(folder)
    noof_markers=0
    psi_plus_0_sns=array([],dtype='uint32')
    psi_plus_1_sns=array([],dtype='uint32')
    psi_min_0_sns =array([],dtype='uint32')
    psi_min_1_sns =array([],dtype='uint32')
    marker_sns=array([],dtype='uint32')
    bins= np.linspace(5400,5600,2001)
    hists_0=np.empty((len(fps),len(bins)-1), dtype='int' )
    hists_1=np.empty((len(fps),len(bins)-1), dtype='int' )
    for i,fp in enumerate(fps):
        f=h5py.File(fp,'r')
        print 'Running{}/{}'.format(i,len(fps))
        sn = f['/PQ_sync_number-1'].value 
        sp = f['/PQ_special-1'].value      
        st = f['/PQ_sync_time-1'].value
        #_PQ_time = f['/PQ_time'].value   
        ch = f['/PQ_channel-1'].value
        f.close()
        if False:
            noof_markers+=np.sum(sp)
            marker_sn=sn[sp==1]
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
            psi_plus_0_sns=np.hstack((psi_plus_0_sns,psi_plus_0_sn))
            psi_plus_1_sns=np.hstack((psi_plus_1_sns,psi_plus_1_sn))
            psi_min_0_sns =np.hstack((psi_min_0_sns, psi_min_0_sn))
            psi_min_1_sns =np.hstack((psi_min_1_sns, psi_min_1_sn))
        if True:
            hist_0,_tmp = np.histogram(st[(ch==0) & (sp==0)]/1000., bins=bins)
            hist_1,_tmp = np.histogram(st[(ch==1) & (sp==0)]/1000., bins=bins)
            hists_0[i]=hist_0
            hists_1[i]=hist_1

#compares the PLU entanglement markers, and the raw entanglement events
if False:
    all_ent=np.hstack((psi_min_0_sns,psi_min_1_sns,psi_plus_0_sns,psi_plus_1_sns))
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
if True:
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