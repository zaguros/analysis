import time
folder = r'D:\measuring\data\2015-02-Entanglement_XX_LOTR\BS'#'D:\measuring\data\2014-12-Entanglement_XX_lock\BS'#r'D:\measuring\data\2014-11-Entanglement_XX_data\BS'# r'D:\measuring\data\21411_ZZ\BS'#


filepaths = tb.get_all_msmt_filepaths(folder, pattern='New') 
co = np.ones([1,4])
sync_num_name = 'PQ_sync_number-1'

st_start_ch0 = 5430000
st_len   = 200000 #50 ns
ch0_ch1_diff = 600 #1 ns
st_start_ch1=st_start_ch0 + ch0_ch1_diff
p_sep = 350000 #600 ns

#load all BS data and collect the coincidences
if True:
    for i,f in enumerate(filepaths):
        print 'Processing {}/{}'.format(i, len(filepaths))

        pqf = pq_tools.pqf_from_fp(f, rights = 'r+')
        if sync_num_name in pqf.keys():
            st=pqf['/PQ_sync_time-1'].value
            st_fltr_0 = (((st_start_ch0<=st)  & (st<(st_start_ch0+st_len))) | ((st_start_ch0+p_sep<=st) & (st<(st_start_ch0+p_sep+st_len))) )  
            st_fltr_1 = (((st_start_ch1<=st)  & (st<(st_start_ch1+st_len))) | ((st_start_ch1+p_sep<=st) & (st<(st_start_ch1+p_sep+st_len))) )
            co = pq_tools.get_coincidences(pqf, index = 1, fltr0=st_fltr_0, fltr1=st_fltr_1, force_coincidence_evaluation = True, save = False)

            if i == 0:
                cos=co         
            else:
                cos = np.vstack((cos, co))
        pqf.close()
    fn = '{}_coincidences'.format(time.strftime('%Y%m%d_%H%M'))
    np.savez(os.path.join(folder,fn), cos)

#load coincidences from previously saved file
if False:
    fn='20141126_0954_coincidences.npz'
    fn ='20141203_0918_coincidences.npz'
    fn = '20150219_2251_coincidences.npz'
    cos = np.load(os.path.join(folder,fn))['arr_0']
#load all tpqi so far:
if False:
    folder=r'D:\measuring\data\20141206_all_tpqi_combined'
    fn='20141125_1755_coincidences.npz'
    cos = np.load(os.path.join(folder,fn))['arr_0']
    fn='20141126_0954_coincidences.npz'
    cos = np.vstack((cos,np.load(os.path.join(folder,fn))['arr_0']))
    fn='20141203_0918_coincidences.npz'
    cos = np.vstack((cos,np.load(os.path.join(folder,fn))['arr_0']))
    fn='20141205_1544_coincidences.npz'
    cos = np.vstack((cos,np.load(os.path.join(folder,fn))['arr_0']))

#plot the tpqi visibility vs something (or plot histogram for pts=1)
if True:
    pts=1
    Xs=np.linspace(100,100,pts)
    Vs=zeros(pts)
    Tots=zeros(pts)
    for i,x in enumerate(Xs):
        st_start_ch0 = 5443000.
        st_len   = 50*1000. #50 ns
        ch0_ch1_diff = 1000 #1 ns
        st_start_ch1=st_start_ch0 + ch0_ch1_diff
        p_sep = 350000#+14000 #600 ns
        st0=cos[:,1]
        st1=cos[:,2]
        st_fltr_0 = (((st_start_ch0<=st0)  & (st0<(st_start_ch0+st_len))) | ((st_start_ch0+p_sep<=st0) & (st0<(st_start_ch0+p_sep+st_len))) )  
        st_fltr_1 = (((st_start_ch1<=st1)  & (st1<(st_start_ch1+st_len))) | ((st_start_ch1+p_sep<=st1) & (st1<(st_start_ch1+p_sep+st_len))) )
        cos_fltr=cos[st_fltr_0 & st_fltr_1 ,0]/1000.
        ll=len(cos_fltr)
        dt=80
        ii=np.sum((cos_fltr>-dt) & (cos_fltr<dt))
        #print ii
        oo=np.sum((cos_fltr>(-dt-p_sep/1000)) & (cos_fltr<(dt-p_sep/1000)))+np.sum((cos_fltr>(-dt+p_sep/1000)) & (cos_fltr<(dt+p_sep/1000)))
        #print oo
        Tots[i]=ii+oo
        Vs[i]=1-ii/float(oo)

    if pts>1:
        figure()
        print Xs
        print Vs
        print Tots
        Xs=Xs
        ax=plt.subplot(111)
        ax.plot(Xs,Vs)
        ax.set_xlabel('windows length [ns]')
        ax.set_ylabel('TPQI Visibility')
        ax2=ax.twinx()
        ax2.plot(Xs,Tots, color='g')
        ax2.set_ylabel('Total number of events')
    else:
        figure()
        ax=plt.subplot(111)
        bins=np.linspace(-600,600,30)
        hist_TPQI, bins_TPQI, patches_TPQI = ax.hist(cos_fltr,bins=bins, color='b', cumulative=False, histtype = 'step')
        ax.set_title('TPQI histogram', fontsize =17)
        ax.set_ylabel('Counts',fontsize = 15)
        ax.tick_params(axis = 'y', labelsize = 13)
        ax.set_xlabel('Time [ns]',fontsize = 15)
        ax.tick_params(axis = 'x', labelsize = 13)
        #ax.set_xlim(-1000,30)
#hist(cos[st_fltr_0 & st_fltr_1 ,0]/1000., bins=linspace(-700,700,10))