"""
script to analyse lengthscan data
"""
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from matplotlib import pyplot as plt
import numpy as np
from analysis.scripts.cavity import peakdetect; 
from analysis.scripts.cavity import cavity_general_analysis as cga; reload(cga)




def lengthscan_analysis(folder,**kw):
    do_get_peaks = kw.pop('do_get_peaks',False)


    a = cga.cavity_analysis(folder)
    x,y = a.get_lengthscan_data()

    total_sweep_time = a.f.attrs['nr_steps']*(a.f.attrs['ADC_averaging_cycles']+a.f.attrs['wait_cycles'])/3*10/1e6 # in us (adwin cycle 300 kHz)in seconds
    print 'total sweep time = ', total_sweep_time 

    fig,ax = plt.subplots()
    ax.set_title(folder+'\ntotalsweeptime=%.1fs'%(total_sweep_time))
    ax.set_xlabel('laser tuning voltage (V)')
    ax.set_ylabel('PD signal (V)')
    plt.plot(x,y)

    if do_get_peaks:
        minimum_peak_height = kw.pop('minimum_peak_height',0.4*max(a.y_data))
        minimum_peak_distance = kw.pop('minimum_peak_distance',60)
        a.peak_xs = np.array([])
        a.peak_ys = np.array([])

        indices = peakdetect.detect_peaks(a.y_data,mph=minimum_peak_height,mpd=minimum_peak_distance)
        
        for ii in indices:
            a.peak_xs = np.append(a.peak_xs,a.x_data[ii])
            a.peak_ys = np.append(a.peak_ys,a.y_data[ii])
        
        ax.plot(a.peak_xs,a.peak_ys,'+', mfc=None, mec='r', mew=2, ms=8)

    fig.savefig(a.folder+'/PDsignal_vs_length_tuning_voltage.png')
    plt.show()



    a.finish()

    if do_get_peaks:
        return a.peak_xs

