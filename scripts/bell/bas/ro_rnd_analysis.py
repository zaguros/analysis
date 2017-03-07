f=h5py.File(r'D:\measuring\data\2014-11-Entanglement_XX_data\20141127171518_total_events_m1_m1.hdf5','r')
db=f['analysis']['total_ent_events'].value
d3=f['analysis']['total_lt3_ssro'].value
d4=f['analysis']['total_lt4_ssro'].value
f.close()

sn   		= 0	
valid 		= 1
ch_1 		= 2	
ch_2  		= 3	
st_1 		= 4	
st_2 		= 5	
t_1			= 6
t_2			= 7

sn_ma = 0
sn_ro = 1
noof_ph_ro=2
st_ma = 3
noof_rnd = 4
noof_rnd_0 = 5
noof_rnd_1 = 6

noof_ev=len(db)
print 'noof events {}'.format(noof_ev)
#sanity checks
ro_ms0_lt3=np.sum(d3[:,noof_ph_ro]>0)/float(noof_ev)
ro_ms0_lt4=np.sum(d4[:,noof_ph_ro]>0)/float(noof_ev)
print 'RO ms=0 LT3: {:.2f}%, RO ms=0 LT4: {:.2f}% '.format(ro_ms0_lt3*100,ro_ms0_lt4*100 )

rnd_0_lt3=np.sum(d3[:,noof_rnd_0]==0)/float(noof_ev)
rnd_0_lt4=np.sum(d4[:,noof_rnd_0]==0)/float(noof_ev)
print 'RND0 LT3: {:.2f}%, RND0 LT4: {:.2f}% '.format(rnd_0_lt3*100,rnd_0_lt4*100 )
print 'no of events without RND generated: {} '.format(np.sum(d3[:,noof_rnd]==0))

pts=1
X=np.linspace(1000,100000,pts)
Y=np.zeros((pts,2))
T=np.zeros((pts,2))
for k,x in enumerate(X):
	#tail/laser filter
	st_start_ch0 = 5444500
	st_len   = 50000 #50 ns
	ch0_ch1_diff = 600 #1 ns
	st_start_ch1=st_start_ch0 + ch0_ch1_diff
	p_sep = 600000 #600 ns
	st_fltr_w1 = (((st_start_ch0<=db[:,st_1])       & (db[:,st_1]<(st_start_ch0+st_len)) 	  & (db[:,ch_1]==0)) | ((st_start_ch1<=db[:,st_1]) 	   & (db[:,st_1]<(st_start_ch1+st_len)) 		 & (db[:,ch_1]==1)) )  
	st_fltr_w2=  (((st_start_ch0+p_sep<=db[:,st_2]) & (db[:,st_2]<(st_start_ch0+p_sep+st_len)) & (db[:,ch_2]==0)) | ((st_start_ch1+p_sep<=db[:,st_2]) & (db[:,st_2]<(st_start_ch1+p_sep+st_len)) & (db[:,ch_2]==1)) ) 

	st_fltr = st_fltr_w1 & st_fltr_w2

	#dt filter
	dt = 200000 # 10 ns
	dt_fltr = np.abs(np.array(db[:,st_2] - db[:,st_1], dtype='int')-p_sep) < dt 
	#print np.abs(np.array(d[:,st_2] - d[:,st_1], dtype='int')-p_sep)
	#psi_min_plus_filter
	psi_min_fltr = (db[:,ch_1]==1) & (db[:,ch_2]==0)# db[:,ch_1] != db[:,ch_2] #(db[:,ch_1]==0) & (db[:,ch_2]==) 
	psi_plus_fltr= db[:,ch_1] == db[:,ch_2] #(db[:,ch_1]==1) & (db[:,ch_2]==1) 

	#RND err filter
	rnd_fltr= (d3[:,noof_rnd_0] + d3[:,noof_rnd_1] ==1 ) & (d4[:,noof_rnd_0] + d4[:,noof_rnd_1] ==1 ) #\
	          #& (d3[:,noof_rnd]>0) & (d4[:,noof_rnd]>0) & 

	valid_event_fltr = db[:,valid] == 1
	valid_event_fltr_SP = (db[:,valid] == 2) | (db[:,valid] == 3)

	#filtered data
	for psi,psi_fltr in zip(['psi_minw1','psi_plusw2'], [psi_min_fltr,psi_plus_fltr]): #[st_fltr_w1 ,st_fltr_w2]):#
		fltr = st_fltr &  dt_fltr & valid_event_fltr & rnd_fltr & psi_fltr
		#fltr = valid_event_fltr_SP & rnd_fltr & psi_fltr
		db_fltr = db[fltr]
		d3_fltr = d3[fltr]
		d4_fltr = d4[fltr]
		print '-'*40
		print 'FILTERED EVENTS {} \n  noof events {}\n'.format(psi,len(db_fltr))

		rnd_corr=[[0,0],[0,1],[1,0],[1,1]] #'RND [LT3,LT4]'
		ro_corr =[[1,1],[1,0],[0,1],[0,0]] #'RO [LT3,LT4] ms  00, 01, 10, 11'

		corr_mat=np.zeros((4,4))
		Es=np.zeros(4)
		for i,rnd in enumerate(rnd_corr):
			for j,ro in enumerate(ro_corr):
				corr_mat[i,j] =np.sum( (d3_fltr[:,noof_rnd_0] == rnd[0]) \
									 & (d4_fltr[:,noof_rnd_0] == rnd[1]) \
									 & ((d3_fltr[:,noof_ph_ro]>0) == ro[0]) \
									 & ((d4_fltr[:,noof_ph_ro]>0) == ro[1]))
				#print 'rnd: {}, ro {}, N {}'.format(rnd,ro,corr_mat[i,j])
			Es[i] = (corr_mat[i,0] - corr_mat[i,1] - corr_mat[i,2] + corr_mat[i,3])/float(np.sum(corr_mat[i,:]))

		print 'RO ms  00, 01, 10, 11'
		print 'RND00', corr_mat[0], '  +pi/2, +3pi/4 (?)'
		print 'RND01', corr_mat[1], '  +pi/2, -3pi/4 (?)'
		print 'RND10', corr_mat[2], '  0,     +3pi/4 (?)'
		print 'RND11', corr_mat[3], '  0,     -3pi/4 (?)'

		print '\n E (00     01    10    11 )'
		print '   ({:.2f}, {:.2f}, {:.2f}, {:.2f})'.format(Es[0],Es[1],Es[2],Es[3])
		print 'RO ms  00, 01, 10, 11'
		print '    ', np.sum(corr_mat,axis=0)
		if psi == 'psi_min':
			CHSH  = -Es[0] + Es[1] + Es[2] + Es[3]
			Y[k,0]= Es[0]
			T[k,0]=len(db_fltr)
		elif psi == 'psi_plus':
			CHSH  = Es[0] - Es[1] + Es[2] + Es[3] 
			Y[k,1]= Es[0]
			T[k,1]=len(db_fltr)
		print 'CHSH', CHSH

#figure()
#bins=linspace(5435,5465,300)
#hist(d[d[:,ch_1]==0,st_1]/1000., bins=bins, log=True, color='r', histtype='step')
#hist(d[d[:,ch_1]==1,st_1]/1000., bins=bins, log=True, color='b', histtype='step')
if pts>1:
    figure()
    ax=plt.subplot(111)
    ax.plot(X/1000.,Y[:,0])
    ax.plot(X/1000.,Y[:,1])
    ax2=ax.twinx()
    ax2.plot(X/1000.,T[:,0], color='r')