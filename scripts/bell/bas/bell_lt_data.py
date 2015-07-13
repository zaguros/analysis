from analysis.scripts.bell.bas import bell_events as be

def sp_corr(sn_lt,sp_lt,ch_lt,st_lt, lt3):
    ro_channel    = 0
    ro_start      = 10450
    ro_length     = 3700

    st_start = 7480 if lt3 else 5350
    st_len       = 200 #50 ns
    st_len_w2    = 200
    p_sep        = 350 #600 ns

    st_fltr_w1 = (sp_lt == 0) & (ch_lt == ro_channel)  & (st_lt > st_start)  & (st_lt < (st_start  + st_len)) 
    st_fltr_w2 = (sp_lt == 0) & (ch_lt == ro_channel)  & (st_lt > st_start + p_sep)  & (st_lt < (st_start + p_sep + st_len_w2)) 
    fltr_ro    = (sp_lt == 0) & (ch_lt == ro_channel)  & (st_lt > ro_start)  & (st_lt < (ro_start  + ro_length))   
    ro0_sn = sn_lt[fltr_ro]
    noof_ro0 = np.sum(fltr_ro)
    w1_sn = sn_lt[st_fltr_w1]
    w2_sn = sn_lt[st_fltr_w2]
    noof_w1_ro0 = np.sum(np.in1d(w1_sn,ro0_sn))
    noof_w2_ro0 = np.sum(np.in1d(w2_sn,ro0_sn))
    noof_w1 = np.sum(st_fltr_w1)
    noof_w2 = np.sum(st_fltr_w2)


    p_w1_ro0 = float(noof_w1_ro0)/noof_w1 if noof_w1 > 0 else 0.
    #u_p_w1_ro0 = np.sqrt(p_w1_ro0 * (1.-p_w1_ro0)/noof_w1)

    p_w2_ro0 = float(noof_w2_ro0)/noof_w2 if noof_w2 > 0 else 0.
    #u_p_w2_ro0 = np.sqrt(p_w2_ro0 * (1.-p_w2_ro0)/noof_w2)

    return p_w1_ro0, p_w2_ro0, noof_w1, noof_w2

def sp_corr_rnd(rnd_channel,sn,sp,ch,st, lt3):
    filter_rnd    = True
    rnd_start     = 10000
    rnd_length    = 1000

    fltr_rnd   = (sp == 1) & (ch == rnd_channel)&(st > rnd_start) & (st < (rnd_start  + rnd_length))
    sn_rnd = sn[fltr_rnd]
    sn_rnd_fltr = np.in1d(sn,sn_rnd)
    sp_lt = sp[sn_rnd_fltr]
    ch_lt = ch[sn_rnd_fltr]
    sn_lt = sn[sn_rnd_fltr]
    st_lt = st[sn_rnd_fltr]
    return sp_corr(sn_lt,sp_lt,ch_lt,st_lt, lt3)

def sp_corr_double_rnd(rnd_channel,rnd_channel_prev,sn,sp,ch,st, lt3):
    filter_rnd    = True
    rnd_start     = 10000
    rnd_length    = 1000

    fltr_rnd   = (sp == 1) & (ch == rnd_channel)&(st > rnd_start) & (st < (rnd_start  + rnd_length))
    sn_rnd = sn[fltr_rnd]
    fltr_rnd_prev  = (sp == 1) & (ch == rnd_channel_prev)&(st > rnd_start) & (st < (rnd_start  + rnd_length))
    sn_rnd_prev = sn[fltr_rnd_prev]
    sn_rnd_fltr = np.in1d(sn,sn_rnd) & np.in1d(sn,sn_rnd_prev+1)
    sp_lt = sp[sn_rnd_fltr]
    ch_lt = ch[sn_rnd_fltr]
    sn_lt = sn[sn_rnd_fltr]
    st_lt = st[sn_rnd_fltr]
    return sp_corr(sn_lt,sp_lt,ch_lt,st_lt, lt3)

load_from_analysis_file = True

#folder = r'K:\ns\qt\Diamond\Projects\Bell\Data\Marked Bell Data\LT4'
folders = [r'D:\measuring\data\2015-04-21-lhfbt4']#r'D:\measuring\data\2015-03-18-lhfbt2',r'D:\measuring\data\2015-04-13-lhfbt3'
#folder = r'D:\measuring\data\2015-04-21-lhfbt4'
#
if load_from_analysis_file:
    afname = 'total_lt3_ssro_fps'
    lt3 = 'lt3' in afname
    fps=np.empty((0,),dtype='str')
    for folder in folders:
        afp = be.get_latest_analysis_fp(folder, pattern ='total_events.hdf5')
        af = h5py.File(afp,'r')
        fps=np.hstack((fps,np.unique(af['analysis'][afname].value)))
        af.close()
else:
    fps=tb.get_all_msmt_filepaths(folder, pattern='third_ever')
#load all data
if True:
    load_all=False
    if load_all:
        PQ_sync_number = np.empty((0,), dtype = np.uint32) 
        PQ_special = np.empty((0,), dtype = np.uint32)         
        PQ_sync_time = np.empty((0,), dtype = np.uint64)
        PQ_time = np.empty((0,), dtype = np.uint64)      
        PQ_channel = np.empty((0,), dtype = np.uint32)
    rnd_zeros=0
    rnd_ones=0
    prepulses = np.zeros(len(fps))
    tails = np.zeros(len(fps))
    p_double = np.zeros(len(fps))
    sp_corrs = np.zeros((len(fps),4,4))
    bins=np.linspace(0,8000,8000)

    print 'starting data loading'
    for i,fp in enumerate(fps):
        try: 
            f=h5py.File(fp,'r')
            print 'PROGRESS {}/{}'.format(i,len(fps)),
            sn = f['/PQ_sync_number-1'].value 
            sp = f['/PQ_special-1'].value      
            st = f['/PQ_sync_time-1'].value
            #_PQ_time = f['/PQ_time-1'].value   
            ch = f['/PQ_channel-1'].value

            if False:
                st_start = 7480 if lt3 else 5350
                st_len       = 200 #50 ns
                st_len_w2    = 200
                p_sep        = 350 #600 ns
                ro_channel=0
                st_fltr_w1 = (sp == 0) & (ch == ro_channel)  & (st > st_start)  & (st < (st_start  + st_len)) 
                st_fltr_w2 = (sp == 0) & (ch == ro_channel)  & (st > st_start + p_sep)  & (st < (st_start + p_sep + st_len_w2)) 
           
                noof_w1 = np.sum(st_fltr_w1)
                noof_w2 = np.sum(st_fltr_w2)
                w1_sn = sn[st_fltr_w1]
                w2_sn = sn[st_fltr_w2]
                noof_w1w2 = len(np.intersect1d(w1_sn,w2_sn))
                reps = sn[-1]
                p_double[i] = float(noof_w1w2)/(noof_w1+noof_w2)

            if True:
                #rnd_zeros+=np.sum((sp==1)&(ch==1))
                #rnd_ones+=np.sum((sp==1)&(ch==2))
                #for j in range(2):
                #    sp_corrs[i,j] = sp_corr_rnd(j+1,sn,sp,ch,st, lt3)

                for j in range(4):
                    rnd_ma_chs=[[1,1], [1,2], [2,1], [2,2]]
                    sp_corrs[i,j] = sp_corr_double_rnd(rnd_ma_chs[j][0],rnd_ma_chs[j][1],sn,sp,ch,st, lt3)
            
            if load_all:
                PQ_sync_number = np.hstack((PQ_sync_number,sn)) 
                PQ_special = np.hstack((PQ_special, sp))         
                PQ_sync_time = np.hstack((PQ_sync_time, st)) 
                #PQ_time = np.hstack((PQ_time, _PQ_time))      
                PQ_channel = np.hstack((PQ_channel, ch))

            if False:
                photons_st = st[(ch==0) & (sp==0)]
                bins=np.linspace(0,15000,15000)
                #figure()
                #hist(photons_st, bins=bins,color='g', histtype='step', normed=True)
                prepulses[i] = np.sum((photons_st>(7560+350)) & (photons_st<=(7568+350)))if lt3 else np.sum((photons_st>5370) & (photons_st<=5383))
                tails[i] = np.sum((photons_st>(7575+350)) & (photons_st<=(7575+100+350)))if lt3 else np.sum((photons_st>5383) & (photons_st<=5500))
                #hist(photons_st, bins=bins, histtype='step', normed=True)
        except KeyError:
            print 'skipping', fp
        finally:
            f.close()

if False:
    figure()
    x =np.arange(len(tails))
    plt.plot(x[tails>0], prepulses[tails>0].astype(np.float64)/tails[tails>0], 'o')

if False:
    figure()
    x =np.arange(len(fps))
    plt.plot(x,p_double, 'o')

if True:
    figure()
    ax=plt.subplot(111)
    for j in range(4):
        for k in range(2):
            ax.plot(sp_corrs[:,j,k], label = 'F0 rnd{}, w{}'.format(str(rnd_ma_chs[j]),k+1))
    ax.legend()
    ax.set_ylim(0,1)
    figure()
    for j in range(4):
        w1 = sp_corrs[:,j,2]
        w2 = sp_corrs[:,j,3]
        plot(w1[w2>0]/w2[w2>0], label = 'w1/w2 rnd{}'.format(j))
    plt.legend()

#load and plot all the master of space positions of the runs.
if False:
    #zs = np.zeros(len(fps))
    #ys = np.zeros(len(fps))
    #xs = np.zeros(len(fps))
    a_vs = np.zeros(len(fps))
    e_vs = np.zeros(len(fps))
    p_vs = np.zeros(len(fps))
    ro_vs= np.zeros(len(fps))
    for i,fp in enumerate(fps):
        try:
            f=h5py.File(fp,'r')
            for k in f.keys():
                if type(f[k])==h5py.Group:
                    g=f[k]
                #print k   
            #zs[i] = (g['instrument_settings']['master_of_space'].attrs)['z']
            #ys[i] = (g['instrument_settings']['master_of_space'].attrs)['y']
            #xs[i] = (g['instrument_settings']['master_of_space'].attrs)['x']
            e_vs[i] = g.attrs['Ex_CR_voltage']
            a_vs[i] = g.attrs['A_CR_voltage']
            p_vs[i] = g.attrs['aom_amplitude']
            ro_vs[i]= g.attrs['RO_voltage_AWG']
        except KeyError:
            print 'skipping', fp
        finally:
            f.close()
    #figure()
    plot(a_vs,'o', label='A CR voltage')
    plot(e_vs,'o', label='E CR voltage')
    plot(p_vs,'o', label='Pulse voltage')
    plot(ro_vs,'o', label='AWG RO voltage')
    #plot(xs, label='x')
    #plot(ys, label='y')
    #plot(zs, label='z')
    plt.legend()
#plot a histogram of the varous events
if False:
    figure()
    bins=np.linspace(0,15000,15000)
    hist(PQ_sync_time[(PQ_channel==0)&(PQ_special==0)], bins=bins,color='g', histtype='step', normed=True)
    #hist(PQ_sync_time[PQ_channel==1], bins=bins,color='b', histtype='step')
    #hist(PQ_sync_time[PQ_channel==2], bins=bins,color='g', histtype='step')
    #hist(PQ_sync_time[PQ_channel==4], bins=bins,color='b', histtype='step')

#plot the time difference between the arrival of the entanglement markers 
#and the RND marker in the previous sync.
if False:
    m_sn =PQ_sync_number[(PQ_channel==4)&(PQ_sync_time<15000)]
    noof_m = len(m_sn)
    e_sn=m_sn-1
    m_t =PQ_time[(PQ_channel==4)&(PQ_sync_time<15000)]
    r_t=np.empty(noof_m)
    for i,s in enumerate(e_sn):
        try:
            r_t[i] = PQ_time[(PQ_sync_number==s) & (PQ_special==1) & (PQ_channel!=4)]
        except:
            continue
    dt_m_r=m_t-r_t
    figure()
    hist(dt_m_r)

if False:
    m_sn =PQ_sync_number[(PQ_channel==4)&(PQ_sync_time<15000)]
    noof_m = len(m_sn)
    e_sn=m_sn-1
    noof_0_markers=np.empty(noof_m)
    noof_1_markers=np.empty(noof_m)
    for i,s in enumerate(e_sn):
        noof_0_markers[i] = np.sum((PQ_sync_number==s) & (PQ_special==1) & (PQ_channel==1))
        noof_1_markers[i] = np.sum((PQ_sync_number==s) & (PQ_special==1) & (PQ_channel==2))
    print 'too many markers:',np.sum((noof_0_markers > 0) & (noof_1_markers > 0))