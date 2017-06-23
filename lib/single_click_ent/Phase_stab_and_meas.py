import os
import numpy as np
from matplotlib import pyplot as plt
import h5py
from analysis.lib.fitting import fit, common
from analysis.lib.nv import detection
from analysis.lib.tools import toolbox as tb
from analysis.lib.m2.ssro import ssro
from analysis.lib.purification import purify_pq as ppq; reload(ppq)
from analysis.lib.tools import plot; reload(plot)
import Analysis_params_SCE as analysis_params; reload(analysis_params)
import scipy.fftpack

def analyze_phase(contains, do_only_meas,do_only_stable, do_both):
    # Import
    lt3_analysis = kw.pop('lt3_analysis', False)

    if not(lt3_analysis):
        folder= tb.latest_data(contains)
    else:
        base_folder_lt4 = analysis_params.data_settings['base_folder_lt4']
        lt4_folder = os.path.join(base_folder_lt4,contains)
        filename_str = kw.pop('filename_str', analysis_params.data_settings['filenames_for_expms'][contains])
        a_list=tb.latest_data(contains = filename_str,folder =lt4_folder,**kw)

    a = ppq.purifyPQAnalysis(folder, hdf5_mode='r')

    # general params
    

    delay_stable = a.g.attrs['count_int_time_stab']
    delay_meas = a.g.attrs['count_int_time_meas']
    pid_cycles = a.g.attrs['pid_points_to_store']
    if a.g.attrs['do_post_ent_phase_msmt'] and do_only_meas:
        sample_cycles = 1
        max_repetitions = a.g['adwindata']['completed_reps'].value
    else:
        sample_cycles = a.g.attrs['sample_points']
        max_repetitions = a.g['adwindata']['store_index_stab'].value/a.g.attrs['pid_points_to_store']

    g_0 = a.g.attrs['Phase_Msmt_g_0']
    visibility = a.g.attrs['Phase_Msmt_Vis']

    if do_both:
        pid_counts_1 = a.g['adwindata']['pid_counts_1'].value
        sample_counts_1 = a.g['adwindata']['sampling_counts_1'].value
        pid_counts_2 = a.g['adwindata']['pid_counts_2'].value
        sample_counts_2 = a.g['adwindata']['sampling_counts_2'].value
        u_1 = []
        u_2 = []
        v_1 = []
        v_2 = []
        angle = []

        for i in xrange(max_repetitions):
            for j in xrange(pid_cycles):
                u_1.append(pid_counts_1[j+i*pid_cycles])
                v_1.append(None)
                u_2.append(pid_counts_2[j+i*pid_cycles])
                v_2.append(None)
                cosval = 2*(float(u_1[-1])/(float(u_1[-1])+float(u_2[-1])*g_0)-0.5)*visibility
                cosval = cosval if np.abs(cosval) < 1 else (1.0 * np.sign(cosval))
                angle.append(180*np.arccos(cosval)/np.pi)   
            for k in xrange(sample_cycles):
                u_1.append(None)
                v_1.append(sample_counts_1[k+i*sample_cycles])
                u_2.append(None)
                v_2.append(sample_counts_2[k+i*sample_cycles])
                cosval = 2*(float(v_1[-1])/(float(v_1[-1])+float(v_2[-1])*g_0)-0.5)*visibility
                cosval = cosval if np.abs(cosval) < 1 else (1.0 * np.sign(cosval)) 
                angle.append(180*np.arccos(cosval)/np.pi)

        t = np.arange(0, (len(v_1)), (1))

        # plot counts
        fig = plt.figure(figsize=(17,6))
        ax = plt.subplot(211)
        plt.plot(t, u_1, 'r')
        plt.plot(t, v_1, 'b')
        plt.title('Counts ZPL detector 1 {0}'.format(folder))
        ax.set_xlabel('elapsed time (a.u.)')
        ax.set_ylabel('counts')
        ax2 = plt.subplot(212)
        plt.plot(t, u_2, 'r')
        plt.plot(t, v_2, 'b')
        plt.title('Counts ZPL detector 2 {0}'.format(folder))
        ax2.set_xlabel('elapsed time (a.u.)')
        ax2.set_ylabel('counts')
        plt.tight_layout()  
        fig.savefig(os.path.join(folder, 'both_trace_counts.png'))

        # zoomed counts
        zoom_1 = 1*(pid_cycles+sample_cycles)
        zoom_2 = 2*(pid_cycles+sample_cycles)
        
        fig = plt.figure(figsize=(17,6))
        ax = plt.subplot(211)
        plt.plot(t[zoom_1:zoom_2], u_1[zoom_1:zoom_2], 'r')
        plt.plot(t[zoom_1:zoom_2], v_1[zoom_1:zoom_2], 'b')
        plt.title('Counts ZPL detector 1'.format(folder))
        ax.set_xlabel('elapsed time (a.u.)')
        ax.set_ylabel('counts')
        ax2 = plt.subplot(212)
        plt.plot(t[zoom_1:zoom_2], u_2[zoom_1:zoom_2], 'r')
        plt.plot(t[zoom_1:zoom_2], v_2[zoom_1:zoom_2], 'b')
        plt.title('Counts ZPL detector 2'.format(folder))
        ax2.set_xlabel('elapsed time (a.u.)')
        ax2.set_ylabel('counts')
        plt.tight_layout()
        fig.savefig(os.path.join(folder, 'both_trace_counts_zoomed.png'))

        # phase
        fig = plt.figure(figsize=(17,6))
        ax = plt.subplot(111)
        plt.plot(t, angle, 'r')
        plt.title('Phase of ZPL photons {0}'.format(folder))
        plt.ylim([0,180])
        ax.set_xlabel('elapsed time (a.u.)')
        ax.set_ylabel('angle ($^o$)')
        fig.savefig(os.path.join(folder, 'both_trace_angle.png'))

        # fft
        # yf = scipy.fftpack.fft(v_1)
        # xf = np.linspace(0, 1.0/(2*delay_meas*1e-6), (pid_cycles)*max_repetitions/2)

        # fig, ax = plt.subplots()
        # ymax = 1.2*np.max(yf[10:-10])
        # plt.ylim([0,ymax])
        # ax.plot(xf, yf[:len(yf)/2])
        # plt.title('FFT {0}'.format(folder))
        # ax.set_xlabel('frequency (Hz)')
        # ax.set_ylabel('Amplitude (a.u.)')
        # fig.savefig(os.path.join(folder, 'both_fft.png'))

    if do_only_meas:
        sample_counts_1 = a.g['adwindata']['sampling_counts_1'].value
        sample_counts_2 = a.g['adwindata']['sampling_counts_2'].value
        v_1 = sample_counts_1
        v_2 = sample_counts_2

        cosvals = [2*(float(n0)/(float(n0)+float(n1)*g_0)-0.5)*visibility for n0,n1 in zip(v_1,v_2)]
        cosvals = [cosval if np.abs(cosval) < 1 else (1.0 * np.sign(cosval)) for cosval in cosvals]
        angle = 180*np.arccos(cosvals)/np.pi

        t = np.arange(0, (len(v_1)*delay_meas/1000), (float(delay_meas)/1000))

        # counts
        fig = plt.figure(figsize=(17,6))
        ax = plt.subplot(211)
        plt.plot(t, v_1, 'b')
        plt.title('Counts ZPL detector 1 {0}'.format(folder))
        ax.set_xlabel('elapsed time (ms)')
        ax.set_ylabel('counts')
        ax2 = plt.subplot(212)
        plt.plot(t, v_2, 'b')
        plt.title('Counts ZPL detector 2 {0}'.format(folder))
        ax2.set_xlabel('elapsed time (ms)')
        ax2.set_ylabel('counts')
        plt.tight_layout()  
        fig.savefig(os.path.join(folder, 'only_meas_trace_counts.png'))

        # zoomed counts
        zoom_1 = 1*(sample_cycles)
        zoom_2 = 2*(sample_cycles)
        
        fig = plt.figure(figsize=(17,6))
        ax = plt.subplot(211)
        plt.plot(t[zoom_1:zoom_2], v_1[zoom_1:zoom_2], 'b')
        plt.title('Counts ZPL detector 1'.format(folder))
        ax.set_xlabel('elapsed time (milliseconds)')
        ax.set_ylabel('counts')
        ax2 = plt.subplot(212)
        plt.plot(t[zoom_1:zoom_2], v_2[zoom_1:zoom_2], 'b')
        plt.title('Counts ZPL detector 2'.format(folder))
        ax2.set_xlabel('elapsed time (ms)')
        ax2.set_ylabel('counts')
        plt.tight_layout()
        fig.savefig(os.path.join(folder, 'only_meas_trace_counts_zoomed.png'))

        # phase
        fig = plt.figure(figsize=(17,6))
        ax = plt.subplot(111)
        plt.plot(t, angle, 'r')
        plt.title('Phase of ZPL photons {0}'.format(folder))
        plt.ylim([0,180])
        ax.set_xlabel('elapsed time (ms)')
        ax.set_ylabel('angle ($^o$)')
        fig.savefig(os.path.join(folder, 'only_meas_trace_angle.png'))

        # fft
        yf = scipy.fftpack.fft(v_1)
        xf = np.linspace(0, 1.0/(2*delay_meas*1e-6), (sample_cycles)*max_repetitions/2)

        fig, ax = plt.subplots()
        ymax = 1.2*np.max(yf[15:-15])
        plt.ylim([0,ymax])
        ax.plot(xf[:len(yf)/2], yf[:len(yf)/2])
        xlim = plt.xlim()
        if (xlim[1]>1000):
            plt.xlim(0,1000)
        plt.title('FFT {0}'.format(folder))
        ax.set_xlabel('frequency (Hz)')
        ax.set_ylabel('Amplitude (a.u.)')
        fig.savefig(os.path.join(folder, 'only_meas_fft.png'))

        # histogram
        fig = plt.figure()
        ax = plt.subplot(111)
        hist, bins = np.histogram(angle,bins= 100,normed = True)
        width = np.diff(bins)
        center = (bins[:-1] + bins[1:]) / 2
        ax.bar(center, hist, align='center', width=width)
        g_a = 0.0
        g_x0 = 90
        g_sigma = 45
        g_A = 1/(np.sqrt(2 * np.pi) * g_sigma)

        p0, fitfunc,fitfunc_str = common.fit_gauss(g_a, g_A, g_x0, g_sigma)
        fit_result = fit.fit1d(center,hist, None, p0=p0, fitfunc=fitfunc,
                             ret=True,fixed=[])
        plot.plot_fit1d(fit_result, np.linspace(center[0],center[-1],201), ax=ax, 
                            plot_data=False,print_info = True)
        fig.savefig(os.path.join(folder, 'only_meas_histogram.png'))
        print 'x0, sigma ', fit_result['params_dict']['x0'] , fit_result['params_dict']['sigma'] 

        # standard dev
        fig = plt.figure()
        ax = plt.subplot(111)
        angle_reshape = (np.reshape(angle,[sample_cycles,-1]))

        binsize = 100.0
        var_array = np.zeros([int(np.floor(sample_cycles/binsize)),1])
        for x in range(int(np.floor(sample_cycles/binsize))):
            var_array[x] = np.sqrt(np.var(angle_reshape[binsize*x:binsize*(x+1),:]))


        plt.plot(binsize*t[0:(int(np.floor(sample_cycles/binsize)))], var_array)
        plt.title('Standard deviation of Phase {0}'.format(folder))
        ax.set_xlabel('time (ms)')
        ax.set_ylabel('std dev ($^o$)')
        fig.savefig(os.path.join(folder, 'only_stab_standard_dev.png'))
            


    if do_only_stable:
        pid_counts_1 = a.g['adwindata']['pid_counts_1'].value
        pid_counts_2 = a.g['adwindata']['pid_counts_2'].value
        u_1 = pid_counts_1
        u_2 = pid_counts_2

        cosvals = [2*(float(n0)/(float(n0)+float(n1)*g_0)-0.5)*visibility for n0,n1 in zip(u_1,u_2)]
        cosvals = [cosval if np.abs(cosval) < 1 else (1.0 * np.sign(cosval)) for cosval in cosvals]
        angle = 180*np.arccos(cosvals)/np.pi

        t = np.arange(0, (len(u_1)*delay_stable/1000), (delay_stable/1000))

        # counts
        fig = plt.figure(figsize=(17,6))
        ax = plt.subplot(211)
        plt.plot(t, u_1, 'r')
        plt.title('Counts ZPL detector 1 {0}'.format(folder))
        ax.set_xlabel('elapsed time (ms)')
        ax.set_ylabel('counts')
        ax2 = plt.subplot(212)
        plt.plot(t, u_2, 'r')
        plt.title('Counts ZPL detector 2 {0}'.format(folder))
        ax2.set_xlabel('elapsed time (ms)')
        ax2.set_ylabel('counts')
        plt.tight_layout()  
        fig.savefig(os.path.join(folder, 'only_stab_trace_counts.png'))

        # zoomed counts
        zoom_1 = 1*(pid_cycles)
        zoom_2 = 2*(pid_cycles)
        
        fig = plt.figure(figsize=(17,6))
        ax = plt.subplot(211)
        plt.plot(t[zoom_1:zoom_2], u_1[zoom_1:zoom_2], 'b')
        plt.title('Counts ZPL detector 1'.format(folder))
        ax.set_xlabel('elapsed time (ms)')
        ax.set_ylabel('counts')
        ax2 = plt.subplot(212)
        plt.plot(t[zoom_1:zoom_2], u_2[zoom_1:zoom_2], 'b')
        plt.title('Counts ZPL detector 2'.format(folder))
        ax2.set_xlabel('elapsed time (ms)')
        ax2.set_ylabel('counts')
        plt.tight_layout()
        fig.savefig(os.path.join(folder, 'only_stab_trace_counts_zoomed.png'))

        # phase
        fig = plt.figure(figsize=(17,6))
        ax = plt.subplot(111)
        plt.plot(t, angle, 'r')
        plt.title('Phase of ZPL photons {0}'.format(folder))
        plt.ylim([0,180])
        ax.set_xlabel('elapsed time (ms)')
        ax.set_ylabel('angle ($^o$)')
        fig.savefig(os.path.join(folder, 'only_stab_trace_angle.png'))

        # fft
        yf = scipy.fftpack.fft(u_1)
        xf = np.linspace(0, 1.0/(2*delay_stable*1e-6), (pid_cycles)*max_repetitions/2)

        fig, ax = plt.subplots()
        ymax = 1.2*np.max(yf[15:-15])
        plt.ylim([0,ymax])
        ax.plot(xf, yf[:len(yf)/2])
        xlim = plt.xlim()
        if (xlim[1]>1000):
            plt.xlim(0,1000)
        plt.title('FFT {0}'.format(folder))
        ax.set_xlabel('frequency (Hz)')
        ax.set_ylabel('Amplitude (a.u.)')
        fig.savefig(os.path.join(folder, 'only_stab_fft.png'))

        # histogram
        fig = plt.figure()
        ax = plt.subplot(111)
        hist, bins = np.histogram(angle,bins= 50,normed = True)
        width = np.diff(bins)
        center = (bins[:-1] + bins[1:]) / 2
        ax.bar(center, hist, align='center', width=width)
        g_a = 0.0
        g_x0 = 90
        g_sigma = 20
        g_A = 1/(np.sqrt(2 * np.pi) * g_sigma)

        p0, fitfunc,fitfunc_str = common.fit_gauss(g_a, g_A, g_x0, g_sigma)
        fit_result = fit.fit1d(center,hist, None, p0=p0, fitfunc=fitfunc,
                             ret=True,fixed=[0])
        plot.plot_fit1d(fit_result, np.linspace(center[0],center[-1],201), ax=ax, 
                            plot_data=False,print_info = True)
        fig.savefig(os.path.join(folder, 'only_stab_histogram.png'))
        print 'x0, sigma ', fit_result['params_dict']['x0'] , fit_result['params_dict']['sigma'] 

        # standard dev
        fig = plt.figure()
        ax = plt.subplot(111)
        plt.plot(np.sqrt(np.var(np.reshape(angle,[pid_cycles,-1]),axis=1)))
        plt.title('Standard deviation of Phase {0}'.format(folder))
        ax.set_xlabel('pid cycles (-)')
        ax.set_ylabel('std dev ($^o$)')
        fig.savefig(os.path.join(folder, 'only_stab_standard_dev.png'))
            

