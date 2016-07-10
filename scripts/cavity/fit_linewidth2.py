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



#def fit_linewidth(folder = tb.latest_data('piezo_scan'),threshold=0.01,delta=0.03):
def fit_linewidth(folder = tb.latest_data('piezo_scan', older_than = '20160204145336'),threshold=0.01,delta=0.03, **kw):

    print folder
    show_plots = kw.pop('show_plots',True) #if the keyword show_plots exists, set to the given value, otherwise default is True
    a = cga.cavity_analysis(folder)
    print a.f
    a.get_x_pts() #the x points of the laser of voltage scan
    a.get_data()

    #x_datas contains nr_repetitions copies of a.x_pts
    x_datas = np.full((a.nr_scans,len(a.x_pts)),a.x_pts)  

    single_scan_lw = np.zeros(a.nr_scans)
    u_single_scan_lw = np.zeros(a.nr_scans)
 
    #new_delays = np.zeros(len(a.sweep_pts)*len(a.nr_scans))
    name = ''#reps_0_1_10_11_20_21_30_31_40_41'


    for i in np.arange(a.nr_scans):
        single_scan_lw[i],u_single_scan_lw[i]=cf.fit_piezo_plots(folder,x_datas[i],a.data[i],
        show_plots = True, threshold=threshold, averaging = False)

    return single_scan_lw, u_single_scan_lw #SvD: you want to return not the index i, but the whole array

#-------------------------------------------------------------------

def find_average_linewidth(from_folder='20160204145330', until_folder='20160204145336'):

    avg_lw = np.array([])
    u_avg_lw = 0

    i = int(until_folder)+1
    print i

    while i > int(from_folder):
        print i
        lws, u_lws = fit_linewidth(folder = tb.latest_data('piezo_scan', older_than = str(i)),threshold=0.01,delta=0.03)
        avg_lw = np.append(avg_lw,lws)
        print avg_lw
        i = i - 1

    print 'Voila! Hvala!'

#-------------------------------------------------------------------

    # fig,ax = plt.subplots(figsize=(6,4.7))
    # ax.errorbar(np.arange(a.nr_scans), avg_lws, yerr=u_avg_lws)
    # ax.legend()
    # ax.set_xlabel('scan number',fontsize=14)
    # ax.set_ylabel('average linewidth (nm)',fontsize=14)
    # ax.tick_params(axis = 'both', which = 'major',labelsize = 14)
    # ax.set_xlim(1,a.nr_scans)
    # ax.set_ylim(0,max(avg_lws))

    # ax.set_title(folder )

    # fig.savefig(folder +'/'+"lw_vs_scan_number .png")
    # plt.show()
    # plt.close()
    

if __name__ == '__main__':
    find_average_linewidth()