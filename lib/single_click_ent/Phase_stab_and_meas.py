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

def analyze_phase(contains, mode, plot_zoomed = [], start_rep_no = 1,**kw):
    # Import
    lt3_analysis = kw.pop('lt3_analysis', False)

    if not(lt3_analysis):
        measfile= tb.latest_data(contains)
    else:
        base_folder_lt4 = analysis_params.data_settings['base_folder_lt4']
        folder = os.path.join(base_folder_lt4,contains)
        filename_str = kw.pop('filename_str', analysis_params.data_settings['filenames_for_expms'][contains])
        measfile=tb.latest_data(contains = filename_str,folder =folder,**kw)

    a = ppq.purifyPQAnalysis(measfile, hdf5_mode='r')

    # general params
    

    delay_stable = np.float(a.g.attrs['count_int_time_stab'])
    delay_meas = np.float(a.g.attrs['count_int_time_meas'])
    pid_cycles = a.g.attrs['pid_points_to_store']
    if a.g.attrs['do_post_ent_phase_msmt'] and mode == 'do_only_meas':
        sample_cycles = 1
        print 'ey'
        max_repetitions = a.g['adwindata']['completed_reps'].value
    else:
        sample_cycles = a.g.attrs['sample_points']
        max_repetitions = a.g['adwindata']['store_index_stab'].value/a.g.attrs['pid_points_to_store']

    

    g_0 = a.g.attrs['Phase_Msmt_g_0']
    visibility = a.g.attrs['Phase_Msmt_Vis']

    pid_counts_1 = a.g['adwindata']['pid_counts_1'].value
    sample_counts_1 = a.g['adwindata']['sampling_counts_1'].value
    pid_counts_2 = a.g['adwindata']['pid_counts_2'].value
    sample_counts_2 = a.g['adwindata']['sampling_counts_2'].value

    pid_counts_1 = pid_counts_1[(start_rep_no - 1)*pid_cycles:]
    pid_counts_2 = pid_counts_2[(start_rep_no - 1)*pid_cycles:]
    sample_counts_1 = sample_counts_1[(start_rep_no - 1)*sample_cycles:]
    sample_counts_2 = sample_counts_2[(start_rep_no - 1)*sample_cycles:]

    if mode == 'only_meas':
        
        delay = delay_meas
        total_cycles = sample_cycles

        v_1 = sample_counts_1
        v_2 = sample_counts_2
        t = np.arange(0, (len(v_1)*delay/1000), (float(delay)/1000))



    elif mode == 'only_stab':

        delay = delay_stable
        total_cycles = pid_cycles

        v_1 = pid_counts_1
        v_2 = pid_counts_2
        t = np.arange(0, (len(v_1)*delay/1000), (float(delay)/1000))



    else:
        
        v_1 = []
        v_2 = []
        t = []
        angle = []

        for i in xrange(len(pid_counts_1)/pid_cycles):
            
            v_1.extend(pid_counts_1[i*pid_cycles:((i+1)*pid_cycles)])
            v_2.extend(pid_counts_2[i*pid_cycles:((i+1)*pid_cycles)])
            if len(t) == 0:
                t.extend(delay_stable * (1 + np.arange(pid_cycles)))
            else:
                t.extend(t[-1] + delay_stable * (1 + np.arange(pid_cycles)))

            v_1.extend(sample_counts_1[i*sample_cycles:((i+1)*sample_cycles)])
            v_2.extend(sample_counts_2[i*sample_cycles:((i+1)*sample_cycles)])
            if len(t) == 0:
                t.extend(delay_meas * (1 + np.arange(sample_cycles)))
            else:
                t.extend(t[-1] + delay_meas * (1 + np.arange(sample_cycles)))
        t = np.array(t)
        total_cycles = pid_cycles + sample_cycles

    cosvals = [2*(float(n0)/(float(n0)+float(n1)*g_0)-0.5)*visibility for n0,n1 in zip(v_1,v_2)]
    cosvals = [cosval if np.abs(cosval) < 1 else (1.0 * np.sign(cosval)) for cosval in cosvals]
    angle = 180*np.arccos(cosvals)/np.pi

   


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
    fig.savefig(os.path.join(folder, 'trace_counts.png'))

    if len(plot_zoomed):
        
        fig = plt.figure(figsize=(17,6))
        ax = plt.subplot(111)
        plt.plot(t[total_cycles*plot_zoomed[0]:total_cycles*plot_zoomed[1]], angle[total_cycles*plot_zoomed[0]:total_cycles*plot_zoomed[1]])
        plt.title('Zoomed trace {0}'.format(folder))
        ax.set_xlabel('elapsed time (milliseconds)')
        ax.set_ylabel('Phase')
        plt.tight_layout()
        fig.savefig(os.path.join(folder, 'trace_zoomed.png'))

    # phase
    fig = plt.figure(figsize=(17,6))
    ax = plt.subplot(111)
    plt.plot(t, angle, 'r')
    plt.title('Phase of ZPL photons {0}'.format(folder))
    plt.ylim([0,180])
    ax.set_xlabel('elapsed time (ms)')
    ax.set_ylabel('angle ($^o$)')
    fig.savefig(os.path.join(folder, 'trace_angle.png'))


    # fft

    if mode == 'only_meas' or mode == 'only_stab':

        yf = np.abs(scipy.fftpack.fft(angle))
        xf = np.linspace(0, 1.0/(2*delay*1e-6), len(angle)/2)

        fig, ax = plt.subplots()
        ymax = 1.2*np.max(yf[15:-15])
        plt.ylim([0,ymax])
        ax.plot(xf[:len(yf)], yf[:len(yf)/2])
        xlim = plt.xlim()
        if (xlim[1]>1000):
            plt.xlim(0,1000)
        plt.title('FFT {0}'.format(folder))
        ax.set_xlabel('frequency (Hz)')
        ax.set_ylabel('Amplitude (a.u.)')
        fig.savefig(os.path.join(folder, 'fft.png'))

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

    ax.set_xlabel('Phase')
    fig.savefig(os.path.join(folder, 'histogram.png'))
    print 'x0, sigma ', fit_result['params_dict']['x0'] , fit_result['params_dict']['sigma'] 

    # standard dev
    fig = plt.figure()
    ax = plt.subplot(111)
    angle_reshape = (np.reshape(angle,[total_cycles,-1]))
    binsize = 1.0
    var_array = np.zeros([int(np.floor(total_cycles/binsize)),1])
    for x in range(int(np.floor(total_cycles/binsize))):
        var_array[x] = np.sqrt(np.var(angle_reshape[binsize*x:binsize*(x+1),:]))


    plt.plot(binsize*t[0:(int(np.floor(total_cycles/binsize)))], var_array)
    plt.title('Standard deviation of Phase {0}'.format(folder))
    ax.set_xlabel('time (ms)')
    ax.set_ylabel('std dev ($^o$)')
    # fig.savefig(os.path.join(folder, 'standard_dev.png'))
  
