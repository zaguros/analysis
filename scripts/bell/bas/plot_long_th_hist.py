#folder=r'D:\measuring\data\After2014-11-19analysis\BS\20141120\132743_Bell_BS_full_Belllhfbt_day7_run2'
folder=r'D:\measuring\data\After2014-11-19analysis\LT4'

if True:
	fps=tb.get_all_msmt_filepaths(folder)
	first_sts_all = np.empty((0,), dtype = np.uint64)

	for fp in fps:
		open_size=25000
		f=h5py.File(fp,'r')
		sn = f['/PQ_sync_number'][open_size].value
		#sp = f['/PQ_special-1'][open_size].value      
		st = f['/PQ_sync_time'][open_size].value
		tt = f['/PQ_time'][open_size].value   
		#ch = f['/PQ_channel-1'][open_size].value
		f.close()

		syncs_per_rep = 251
		time_per_sync = 15000000

		tt_first_sync=tt[0] - sn[0]*time_per_sync - st[0]
		first_sts = tt[sn<=syncs_per_rep] - tt_first_sync
		first_sts_all = np.hstack((first_sts_all,first_sts)) 

if True:
	bins=np.linspace(0,260*15e-6*1e12/1e6,20000)
	figure()
	hist(first_sts_all/1e6, bins=bins,color='r', histtype='step')