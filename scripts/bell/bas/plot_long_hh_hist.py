#folder=r'D:\measuring\data\After2014-11-19analysis\BS\20141120\132743_Bell_BS_full_Belllhfbt_day7_run2'
folder=r'D:\measuring\data\After2014-11-19analysis\BS\20141118'

if True:
	fps=tb.get_all_msmt_filepaths(folder)
	first_sts_all = np.empty((0,), dtype = np.uint64)

	for fp in fps:
		f=h5py.File(fp,'r')
		sn = f['/PQ_sync_number-1'].value 
		sp = f['/PQ_special-1'].value      
		st = f['/PQ_sync_time-1'].value
		tt = f['/PQ_time-1'].value   
		ch = f['/PQ_channel-1'].value
		f.close()

		syncs_per_rep = 251
		time_per_sync = 15000000
		 
		rep_indices=np.floor((sn-1)/syncs_per_rep)
		sync_in_rep = np.mod(sn-1,syncs_per_rep)
		reps = int(rep_indices[-1])

		is_first_event_in_rep=np.insert(np.diff(np.asarray(rep_indices, dtype='int'))>0,0,True)
		noof_unique_rep_indeces=np.sum(is_first_event_in_rep)#=np.unique(rep_indices)
		tt_first_sync = tt[is_first_event_in_rep] - sync_in_rep[is_first_event_in_rep]*time_per_sync - st[is_first_event_in_rep] 


		first_event_indices=np.append(np.where(is_first_event_in_rep)[0],len(sn))

		rep_time_0s=np.zeros(len(sn), dtype = np.uint64)
		for i in range(noof_unique_rep_indeces):
			rep_time_0s[first_event_indices[i]:first_event_indices[i+1]] =  tt_first_sync[i]

		first_sts = tt - rep_time_0s

		first_sts_all = np.hstack((first_sts_all,first_sts)) 

if True:
	bins=np.linspace(0,260*15e-6*1e12/1e6,20000)
	figure()
	hist(first_sts_all/1e6, bins=bins,color='r', histtype='step')