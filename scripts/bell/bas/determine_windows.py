import time
from analysis.lib.fitting import fit,common
plt.close('all')
if True:
    folder = r'D:\measuring\data\2015-03-03-lhfbt\BS'#'D:\measuring\data\2014-12-Entanglement_XX_lock\BS'#r'D:\measuring\data\2014-11-Entanglement_XX_data\BS'# r'D:\measuring\data\21411_ZZ\BS'#
    filepaths = tb.get_all_msmt_filepaths(folder, pattern='first_ever') 
    sync_num_name = 'PQ_sync_number-1'
    bins= np.linspace(5400000,5600000,2001)
    hists_0=np.empty((len(filepaths),len(bins)-1), dtype='uint32' )
    hists_1=np.empty((len(filepaths),len(bins)-1), dtype='uint32' )
    reps = np.zeros(len(filepaths))
#load all BS data and collect the coincidences

    for i,fp in enumerate(filepaths):
        print 'Processing {}/{}'.format(i, len(filepaths))

        pqf = pq_tools.pqf_from_fp(fp, rights = 'r+')
        if sync_num_name in pqf.keys():
            st=pqf['/PQ_sync_time-1'].value
            ch = pqf['/PQ_channel-1'].value
            sp = pqf['/PQ_special-1'].value 
            hist_0,_tmp = np.histogram(st[(ch==0) & (sp==0)], bins=bins)
            hist_1,_tmp = np.histogram(st[(ch==1) & (sp==0)], bins=bins)
            hists_0[i]=hist_0
            hists_1[i]=hist_1
            reps[i] = pqf['/PQ_sync_number-1'][-1]
        pqf.close()

#some histogram plotting
if True:
    x = bins[1:]/1000.
    y = np.sum(hists_0,axis=0)
    semilogy(x,y)
    x_range_fltr = (x>=5438) & (x<5442.5)
    xf = x[x_range_fltr]
    yf = y[x_range_fltr]
    #a + A * exp(-(x-x0)**2/(2*sigma**2))
    #['g_a', 'g_A', 'g_x0', 'g_sigma']
    #['g_a1', 'g_A1', 'g_x01', 'g_sigma1', 'g_A2', 'g_x02', 'g_sigma2']
    f = common.fit_2gauss
    x_bg_fltr = (x>=5550) & (x<5600)
    a_bg=np.average(y[x_bg_fltr])
    args=[a_bg,10000,5440.73,2,10000,5440.73,2]
    fitres = fit.fit1d(xf, yf, f, *args, fixed = [0],
                   do_print = True, ret = True)
    yff = fitres['fitfunc'](x)
    semilogy(x,yff)

    #p_ratio = np.array(yff, dtype=np.float)/np.array(y, dtype=np.float)
    #figure()
    #semilogy(x,p_ratio)

    x_start_r = np.arange(5442,5448,0.1)
    x_len_r = np.arange(10,80,1)
    p_ratios=np.zeros((len(x_len_r),len(x_start_r)))
    y_tot = np.zeros((len(x_len_r),len(x_start_r)))
    for i,x_start in enumerate(x_start_r):
        for j,x_len in enumerate(x_len_r):
            x_fltr = (x>=x_start) & (x<(x_start+x_len))
            p_ratios[j,i] = float(np.sum(yff[x_fltr]))/np.sum(y[x_fltr])
            y_tot[j,i] = np.sum(y[x_fltr])/1e6
            
    
    #plt.imshow(p_ratios)
    figure()
    xx,yy = np.meshgrid(x_start_r,x_len_r)
    cs = plt.contour(xx,yy,p_ratios, levels=[0.1,0.08,0.06,0.04,0.02,0.01,0.008,0.006,0.005,0.004])
    plt.clabel(cs, inline=1, fontsize=10)

    figure()
    cs = plt.contour(xx,yy,y_tot)
    plt.clabel(cs, inline=1, fontsize=10)
    p_ratio_range_fltr = (p_ratios<0.005)
    
    ii,jj= numpy.unravel_index(np.argmax(y_tot*p_ratio_range_fltr), np.shape(y_tot))
    print x_start_r[jj], x_len_r[ii]

#tail counts, pulse counts and rejection
if False:
    #np.loadtxt?
    d=np.loadtxt(r'Z:\data\20150310\1540_laser_measurement.dat', skiprows=10)
    #d=np.loadtxt(r'K:\ns\qt\Diamond\Projects\Bell\Data\2014-12-17_laser_lock_spectra\2014-12-18_with_analog_and_pid_locks.csv', skiprows=2, dtype='float', usecols=(0,2), delimiter=',')
    ch0=d[:,0]
    ch1=d[:,1]
    y = ch0
    x=np.arange(5400,5600,0.1280)[:-1]
    print len(x), len(y)
    #plt.semilogy(x,y, '-')
    f = common.fit_2gauss
    args=[0.1,10000,5440.73,2,10000,5440.73,2]
    fitres = fit.fit1d(x, y, f, *args, fixed = [0],
                   do_print = True, ret = True)
    p1 = fitres['params_dict']
    plt.semilogy(x,y, 'o')
    plot_pts=10000
    x_p=np.linspace(min(x),max(x),plot_pts)
    if fitres['success']:
        f_p = fitres['fitfunc'](x_p)
        plt.semilogy(x_p,f_p)
    plt.xlim(5420,5470)
    plt.ylim(1,20000)

if False:
    x = bins[1:]/1000.
    pulse_range_fltr_0 = (x>=5438) & (x<5442.5)
    tail_range_fltr_0 = (x>=5444) & (x<(5444+50))
    pulse_range_fltr_1 = (x>=(5438+0.7)) & (x<(5442.5+0.7))
    tail_range_fltr_1 = (x>=(5444+0.7)) & (x<(5444+50+0.7))
    pts=len(hists_0)
    tails = np.zeros(pts)
    pulses = np.zeros(pts)
    for i in range(pts):
        pulses[i]=float(np.sum(hists_0[i,pulse_range_fltr_0]) + np.sum(hists_1[i,pulse_range_fltr_1]))/(reps[i])*1e4 # factor 2 comes from pi/2 pulse
        tails[i] =float(np.sum(hists_0[i, tail_range_fltr_0]) + np.sum(hists_1[i, tail_range_fltr_1]))/(reps[i])*1e4*2 # factor 2 comes from pi/2 pulse
    ax=plt.subplot(111)
    ax.plot(np.arange(pts)[pulses>0],pulses[pulses>0])
    ax.plot(np.arange(pts)[pulses>0],tails[pulses>0])
    ax2=ax.twinx()
    ax2.plot(np.arange(pts)[pulses>0],pulses[pulses>0]/tails[pulses>0], color='r')
