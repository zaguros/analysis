"""
Script to analyze the data obtained in a sweep of the montana sync delays measurement
it is based on the general fit formulas in analysis_linewidth_fsr 
- SvD 28-1-2016
"""
import os
import h5py
import numpy as np
from analysis.lib.tools import toolbox as tb
from matplotlib import pyplot as plt

from analysis.scripts.cavity import peakdetect as pd; reload(pd)

from analysis.scripts.cavity import cavity_general_analysis as cga; reload(cga)

from analysis.scripts.cavity import cavity_fitting as cf; reload(cf)



def fit_sweep_msync_delays(folder = tb.latest_data('141702'),threshold=0.01,delta=0.05):
    a = cga.cavity_analysis(folder)
    print a.f
    a.get_x_pts() #the x points of the laser of voltage scan
    a.get_sweep_pts() #the sweep points  (here: delays)
    a.get_sweep_data()

    #x_datas contains nr_repetitions copies of a.x_pts
    x_datas = np.full((a.nr_repetitions,len(a.x_pts)),a.x_pts)  

    selected_reps = np.arange(a.nr_repetitions)#np.array([0,5,10,15])#np.arange(a.nr_repetitions)#np.array([0,1,10,11,20,21,30,31,40,41]) #

    avg_lws = np.zeros(len(a.sweep_pts))
    u_avg_lws = np.zeros(len(a.sweep_pts))
    peak_range=np.zeros([len(a.sweep_pts),len(selected_reps)])
    nr_peaks=np.zeros([len(a.sweep_pts),len(selected_reps)])
 
    #new_delays = np.zeros(len(a.sweep_pts)*len(a.nr_scans))
    name = ''#reps_0_1_10_11_20_21_30_31_40_41'


    for i,pt in enumerate(a.sweep_pts):
        for j,rep in enumerate(selected_reps):
            peak_range[i,j],nr_peaks[i,j]=cf.determine_peak_range(folder,a.sweep_data[i,j],a.x_pts,delta,tag=str(pt)+'_'+str(rep),show_plot=False)

    avg_peak_range = np.average(peak_range,axis=1)
    u_avg_peak_range = np.std(peak_range,axis=1)
    avg_nr_peaks = np.average(nr_peaks,axis=1)
    u_avg_nr_peaks = np.std(nr_peaks,axis=1)


    print avg_peak_range
    print u_avg_peak_range

    print avg_nr_peaks
    print u_avg_nr_peaks

    #plot and save peak range vs sync delay
    fig,ax = plt.subplots(figsize=(6,4.7))
    ax.errorbar(a.sweep_pts, avg_peak_range, fmt='o', yerr=u_avg_peak_range)
    # ax.legend()
    ax.set_xlabel('sync delay (ms)',fontsize=14)
    ax.set_ylabel('average peak range (V)',fontsize=14)
    ax.tick_params(axis = 'both', which = 'major',labelsize = 14)
    ax.set_xlim(a.sweep_pts[0],a.sweep_pts[-1])
    ax.set_ylim(0,max(avg_peak_range)+max(u_avg_peak_range))

    ax.set_title(folder +'\n'+'peakrange_vs_delays'+name )

    fig.savefig(folder +'/'+'peakrange_vs_delays'+name+'.png')
    plt.show()
    plt.close()


    #plot and save nr peaks vs sync delay
    fig,ax = plt.subplots(figsize=(6,4.7))
    ax.errorbar(a.sweep_pts, avg_nr_peaks, fmt='o', yerr=u_avg_nr_peaks)
    # ax.legend()
    ax.set_xlabel('sync delay (ms)',fontsize=14)
    ax.set_ylabel('average nr peaks',fontsize=14)
    ax.tick_params(axis = 'both', which = 'major',labelsize = 14)
    ax.set_xlim(a.sweep_pts[0],a.sweep_pts[-1])
    ax.set_ylim(0,max(avg_nr_peaks)+max(u_avg_nr_peaks))

    ax.set_title(folder +'\n'+'nrpeaks_vs_delays'+name)

    fig.savefig(folder +'/'+'nrpeaks_vs_delays'+name+'.png')
    plt.show()
    plt.close()

    #closed the file
    a.finish()

    return avg_peak_range, u_avg_peak_range, avg_nr_peaks, u_avg_nr_peaks

    # for i,pt in enumerate(a.sweep_pts):
    #     avg_lws[i],u_avg_lws[i]=cf.fit_piezo_plots(folder,x_datas,a.sweep_data[i],
    #         show_plots = True, nr_repetitions=a.nr_repetitions,threshold=threshold,tag =str(pt))
    
    # for i,pt in enumerate(a.sweep_pts):
    #     avg_lws[i],u_avg_lws[i]=cf.fit_piezo_plots(folder,a.x_pts,a.avg_sweep_data[i],
    #         show_plots = True, nr_repetitions=1,threshold=threshold,tag =str(pt),averaging = False)
   

    # fig,ax = plt.subplots(figsize=(6,4.7))
    # ax.errorbar(a.sweep_pts, avg_lws, fmt='o', yerr=u_avg_lws)
    # # ax.legend()
    # ax.set_xlabel('sync delay (ms)',fontsize=14)
    # ax.set_ylabel('average linewidth (nm)',fontsize=14)
    # ax.tick_params(axis = 'both', which = 'major',labelsize = 14)
    # ax.set_xlim(a.sweep_pts[0],a.sweep_pts[-1])
    # ax.set_ylim(0,max(avg_lws))

    # ax.set_title(folder )

    # fig.savefig(folder +'/'+"lw_vs_delays.png")
    # plt.show()
    # plt.close()


def fit_sweep_mult_folders(from_folder='20160204145330', until_folder='20160204145336', delta = 0.01):

    avg_lw = np.array([])
    u_avg_lw = 0

    max_nr = 10

    for i in np.arange(max_nr):
        print i
        avg_peak_range, u_avg_peak_range, avg_nr_peaks, u_avg_nr_peaks = fit_sweep_msync_delays(folder = tb.latest_data('piezo_scan', older_than = str(i)),threshold=0.01,delta=0.03)
        avg_lw = np.append(avg_lw,lws)
        print avg_lw
        i = i - 1

    print 'Voila! Hvala!'


if __name__ == '__main__':
    fit_sweep_msync_delays()