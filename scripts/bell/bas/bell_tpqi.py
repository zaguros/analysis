import time
folder = r'D:\measuring\data\After2014-11-19analysis\BS\20141109'# r'D:\measuring\data\2014-11-Entanglement_XX_data\BS'

filepaths = tb.get_all_msmt_filepaths(folder) 
co = np.ones([1,4])
sync_num_name = 'PQ_sync_number-1'

st_start_ch0 = 5444000
st_len   = 100000 #50 ns
ch0_ch1_diff = 600 #1 ns
st_start_ch1=st_start_ch0 + ch0_ch1_diff
p_sep = 600000 #600 ns

#load all BS data and collect the coincidences
if False:
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

#plot the tpqi visibility vs something (or plot histogram for pts=1)
if True:
    pts=1
    Xs=np.linspace(1,50,pts)
    Vs=zeros(pts)
    Tots=zeros(pts)
    for i,x in enumerate(Xs):
        st_start_ch0 = 5444500
        st_len   = 40000 #50 ns
        ch0_ch1_diff = 600 #1 ns
        st_start_ch1=st_start_ch0 + ch0_ch1_diff
        p_sep = 600000 #600 ns
        st0=cos[:,1]
        st1=cos[:,2]
        st_fltr_0 = (((st_start_ch0<=st0)  & (st0<(st_start_ch0+st_len))) | ((st_start_ch0+p_sep<=st0) & (st0<(st_start_ch0+p_sep+st_len))) )  
        st_fltr_1 = (((st_start_ch1<=st1)  & (st1<(st_start_ch1+st_len))) | ((st_start_ch1+p_sep<=st1) & (st1<(st_start_ch1+p_sep+st_len))) )
        cos_fltr=cos[st_fltr_0 & st_fltr_1 ,0]/1000.
        ll=len(cos_fltr)
        dt=300
        ii=np.sum((cos_fltr>-dt) & (cos_fltr<dt))
        #print ii
        oo=np.sum((cos_fltr>(-dt-600)) & (cos_fltr<(dt-600)))+np.sum((cos_fltr>(-dt+600)) & (cos_fltr<(dt+600)))
        #print oo
        Tots[i]=ii+oo
        Vs[i]=1-ii/float(oo)

    if pts>1:
        print Xs
        print Vs
        print Tots
        ax=plt.subplot(111)
        ax.plot(Xs,Vs)
        ax2=ax.twinx()
        ax2.plot(Xs,Tots, color='g')
    else:
        figure()
        bins=np.linspace(-1000,1000,200)
        hist(cos_fltr,bins=bins)
#hist(cos[st_fltr_0 & st_fltr_1 ,0]/1000., bins=linspace(-700,700,10))