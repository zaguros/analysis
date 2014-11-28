folder = r'K:\ns\qt\Diamond\Projects\Bell\Data\Marked Bell Data\LT4'

fps=tb.get_all_msmt_filepaths(folder)

#load all data
if False:
    PQ_sync_number = np.empty((0,), dtype = np.uint32) 
    PQ_special = np.empty((0,), dtype = np.uint32)         
    PQ_sync_time = np.empty((0,), dtype = np.uint64)
    PQ_time = np.empty((0,), dtype = np.uint64)      
    PQ_channel = np.empty((0,), dtype = np.uint32)

    for i,fp in enumerate(fps):
        f=h5py.File(fp,'r')
        print 'PROGRESS {}/{}'.format(i,len(fps))
        _PQ_sync_number = f['/PQ_sync_number'].value 
        _PQ_special = f['/PQ_special'].value      
        _PQ_sync_time = f['/PQ_sync_time'].value
        _PQ_time = f['/PQ_time'].value   
        _PQ_channel = f['/PQ_channel'].value
        f.close()

        PQ_sync_number = np.hstack((PQ_sync_number,_PQ_sync_number)) 
        PQ_special = np.hstack((PQ_special, _PQ_special))         
        PQ_sync_time = np.hstack((PQ_sync_time, _PQ_sync_time)) 
        PQ_time = np.hstack((PQ_time, _PQ_time))      
        PQ_channel = np.hstack((PQ_channel, _PQ_channel))

#plot a histogram of the varous events
if False:
    figure()
    bins=np.linspace(0,15000,1500)
    hist(PQ_sync_time[PQ_channel==0], bins=bins,color='r', histtype='step')
    hist(PQ_sync_time[PQ_channel==1], bins=bins,color='b', histtype='step')
    hist(PQ_sync_time[PQ_channel==2], bins=bins,color='g', histtype='step')
    hist(PQ_sync_time[PQ_channel==4], bins=bins,color='k', histtype='step')

#plot the time difference between the arrival of the entanglement markers 
#and the RND marker in the previous sync.
if False:
    m_sn =PQ_sync_number[(PQ_channel==4)&(PQ_sync_time<15000)]
    noof_m = len(m_sn)
    e_sn=m_sn-1
    m_t =PQ_time[(PQ_channel==4)&(PQ_sync_time<15000)]
    r_t=np.empty(noof_m)
    for i,s in enumerate(e_sn):
        r_t[i] = PQ_time[(PQ_sync_number==s) & (PQ_special==1) & (PQ_channel!=4)]
    dt_m_r=m_t-r_t
    figure()
    hist(dt_m_r)

if True:
    m_sn =PQ_sync_number[(PQ_channel==4)&(PQ_sync_time<15000)]
    noof_m = len(m_sn)
    e_sn=m_sn-1
    noof_0_markers=np.empty(noof_m)
    noof_1_markers=np.empty(noof_m)
    for i,s in enumerate(e_sn):
        noof_0_markers[i] = np.sum((PQ_sync_number==s) & (PQ_special==1) & (PQ_channel==1))
        noof_1_markers[i] = np.sum((PQ_sync_number==s) & (PQ_special==1) & (PQ_channel==2))
    print 'too many markers:',np.sum((noof_0_markers > 0) & (noof_1_markers > 0))