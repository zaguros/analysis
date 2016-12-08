'''
Undocumented...
'''

import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import analysis.lib.QEC.nuclear_spin_characterisation as SC 
import matplotlib.cm as cm
import os
from matplotlib import pyplot as plt
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.tools import analysis_magnet_tools as amt
from analysis.lib.m2.ssro import sequence
from analysis.lib.fitting import dark_esr_auto_analysis; 
import matplotlib.mlab as mlab
from tempfile import TemporaryFile
reload(dark_esr_auto_analysis)
reload(toolbox)
reload(amt)
reload(fit)

ZFS                 = 2.877480e9
g_factor            = 2.8025e6

def fit_B_msmt_loop(older_than = None, newer_than = None):
    ZFS = 2.877480e9
    f0m = []; u_f0m = []; f0p = [] ;u_f0p = []
    Bx_field_measured = []
    Bz_field_measured = []
    f_centre_list = []; f_diff_list=[]
    timestamp_list = []
    f_centre_error_list= []
    it_list = []
    f_diff_error_list = []
    
    #msm
    print 'start'
    older_than_SSRO = older_than
    # print older_than, newer_than
    iteration = 0
    while toolbox.latest_data(contains='msmt_msm_', older_than=older_than, newer_than=newer_than,raise_exc = False)!=False:
        print 'iteration'+str(iteration)
        timestamp,folder = toolbox.latest_data(contains='msmt_msm_', older_than=older_than, newer_than=newer_than,return_timestamp = True)
        print 'm folder '+folder
        ## Data location ##

        ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration', older_than=older_than_SSRO)

        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

        f0m_temp,u_f0m_temp = dark_esr_auto_analysis.analyze_dark_esr(None, 2.196*1e-3, add_folder = folder )
       
        #msp
        folder = toolbox.latest_data(contains='msmt_msp_', older_than=older_than, newer_than=newer_than,)
        print folder
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

        f0p_temp,u_f0p_temp = dark_esr_auto_analysis.analyze_dark_esr(None, 2.196*1e-3, add_folder = folder )
        print f0p_temp >0 and f0m_temp >0
        if f0p_temp >0 and f0m_temp >0:
            Bz_measured, Bx_measured = amt.get_B_field(msm1_freq=f0m_temp*1e9, msp1_freq=f0p_temp*1e9, u_msm1_freq =u_f0m_temp ,u_msp1_freq=u_f0p_temp)
            f_centre = (f0m_temp+f0p_temp)/2
            f_centre_error = np.sqrt(u_f0m_temp**2+u_f0p_temp**2)/2
            f_diff = (f_centre-ZFS*1e-9)*1e6
            f_diff_error = f_centre_error*1e6
            f0m.append(f0m_temp)
            u_f0m.append(u_f0m_temp)
            f0p.append(f0p_temp)
            u_f0p.append(u_f0p_temp)
            f_centre_list.append(f_centre)
            f_centre_error_list.append(f_centre_error)
            f_diff_list.append(f_diff)
            f_diff_error_list.append(f_diff_error)
            Bx_field_measured.append(Bx_measured)
            Bz_field_measured.append(Bz_measured)
            timestamp_list.append(timestamp)
            it_list.append(iteration)

        older_than = str(int(timestamp)-1)

        iteration = iteration+1
        print iteration

    
    it_list = linspace(0,len(it_list)-1,len(it_list))
    outfile = TemporaryFile()
    
    np.savez('test_B_field_meas',f0m=f0m,
            u_f0m=u_f0m,
            f0p=f0p,
            u_f0p=u_f0p,
            f_centre_list=f_centre_list,
            f_centre_error_list=f_centre_error_list,
            f_diff_list=f_diff_list,
            f_diff_error_list=f_diff_error_list,
            Bx_field_measured=Bx_field_measured,
            Bz_field_measured=Bz_field_measured,
            timestamp_list=timestamp_list,
            it_list=it_list)

def plot_meas_B_loop():

    d = np.load('test_B_field_meas.npz')
    f0m = d['f0m']; u_f0m = d['u_f0m']; f0p = d['f0p'] ;u_f0p = d['u_f0p']
    Bx_field_measured = d['Bx_field_measured']
    Bz_field_measured = d['Bz_field_measured']
    f_centre_list = d['f_centre_list']; f_diff_list=d['f_diff_list']
    timestamp_list = d['timestamp_list']
    f_centre_error_list= d['f_centre_error_list']
    it_list = d['it_list']
    f_diff_error_list = d['f_diff_error_list']

    Bz_diff_list = [j-304.21 for j in Bz_field_measured]
    Bz_error = [1/(4.*ZFS*g_factor*Bz_field_measured[j])*
            (f0p[j]**2*u_f0p[j]**2+f0m[j]**2*u_f0m[j]**2)**(1/2.) for j in range(len(Bz_field_measured))]
    Bx_error = [1/Bx_field_measured[j]*(f0m[j]**2*u_f0m[j]**2/g_factor**2
            +(ZFS-g_factor*Bz_field_measured[j])**2*Bz_error[j]**2) for j in range(len(Bx_field_measured))]

    print len(Bz_field_measured)
    total_B = [(Bx_field_measured[j]**2+Bz_field_measured[j]**2)**(1/2.) for j in range(len(Bz_field_measured))]
    total_B_error = [1/total_B[j]*(Bx_field_measured[j]**2*Bx_error[j]**2+Bz_field_measured[j]**2*Bz_error[j]**2)**(1/2.) for j in range(len(Bz_field_measured))]

    mean_f_centre   = np.mean(f_diff_list)
    stdev_f_centre  = np.std(f_diff_list)
    mean_Bx         = np.mean(Bx_field_measured)
    stdev_Bx        = np.std(Bx_field_measured)
    mean_Bz         = np.mean(Bz_diff_list)   
    stdev_Bz        = np.std(Bz_diff_list)

    print np.mean(total_B)
    print np.std(total_B)

    # fig = figure(1,figsize=(18,5))
    # ax = fig.add_subplot(111)
    # ax.errorbar(it_list,f_diff_list,f_diff_error_list)
    # ax.set_xlabel('msmt #')
    # ax.set_ylabel('relative centrer freq (kHz) (offset 2.87748 GHz)')
    # ax.annotate('filled', xy=(301, 10), xycoords='data', xytext=(-50, 110), textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    # ax.annotate('1 day break', xy=(606, 10), xycoords='data', xytext=(-50, 110), textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    # ax.annotate('filled', xy=(612, 10), xycoords='data', xytext=(20, 115), textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    # plt.savefig('freq_vs_time',format='png')
    # figure(2,figsize=(18,5))
    # plt.errorbar(it_list,Bx_field_measured, Bx_error)
    # plt.xlabel('msmt #')
    # plt.ylabel('measured Bx field (G)')
    # plt.savefig('Bx_vs_time',format='png')
    # figure(3,figsize=(18,5))
    # plt.errorbar(it_list,Bz_diff_list,Bz_error)
    # plt.xlabel('msmt #')
    # plt.ylabel('measured relative Bz (G) (offset 304.21 G)')
    # plt.savefig('Bz_vs_time',format='png')
    # figure(4)
    # n, bins, patches = plt.hist(f_diff_list,50,normed = 1)
    # bincenters = 0.5*(bins[1:]+bins[:-1])
    # y = mlab.normpdf( bincenters, mean_f_centre, stdev_f_centre)
    # plt.plot(bincenters, y, 'r--', linewidth=1)
    # plt.xlabel('binned relative center freq (kHz)')
    # plt.title('Mean '+str(mean_f_centre)+' kHz, stdev '+str(stdev_f_centre)+' kHz')
    # plt.savefig('binned_freq',format='png')
    # figure(5)
    # n, bins, patches = plt.hist(Bx_field_measured,50,normed = 1)
    # bincenters = 0.5*(bins[1:]+bins[:-1])
    # y = mlab.normpdf( bincenters, mean_Bx, stdev_Bx)
    # plt.plot(bincenters, y, 'r--', linewidth=1)
    # plt.xlabel('binned Bx (G)')
    # plt.title('Mean '+str(mean_Bx)+' G, stdev '+str(stdev_Bx)+' G')
    # plt.savefig('binned_Bx',format='png')
    # figure(6)
    # n, bins, patches = plt.hist(Bz_diff_list,50,normed = 1)
    # bincenters = 0.5*(bins[1:]+bins[:-1])
    # y = mlab.normpdf( bincenters, mean_Bz, stdev_Bz)
    # plt.plot(bincenters, y, 'r--', linewidth=1)
    # plt.xlabel('binned relative Bz (G) (offset 304.21 G)')
    # plt.title('Mean '+str(mean_Bz+304.21)+' G, stdev '+str(stdev_Bz)+' G')
    # plt.savefig('binned_Bz',format='png')

plot_meas_B_loop()
# fit_B_msmt_loop(older_than = '20140428101106', newer_than = '20140426092819')


