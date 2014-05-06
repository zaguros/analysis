'''
Script to analyze the dynamical decoupling data
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
    f0mc = []; u_f0mc = []; f0pc = [] ;u_f0pc = []
    f0mf = []; u_f0mf = []; f0pf = [] ;u_f0pf = []
    Bx_field_measured = []
    Bz_field_measured = []
    f_centrec_list = []; f_centre_errorc_list= [];f_diffc_list=[];f_diffc_error_list = []
    f_centref_list = []; f_centre_errorf_list= [];f_difff_list=[];f_difff_error_list = []
    timestamp_list = []
    it_list = []
    
    
    #msm
    print 'start'
    older_than_SSRO = older_than
    # print older_than, newer_than
    iteration = 0

    if type == 'course':
        contains_name_m = 'magnet_msm1_coarse'
        contains_name_p = 'magnet_msp1_coarse'
    else:
        contains_name_m = 'magnet_msm1_fine'
        contains_name_p = 'magnet_msp1_fine'

    while toolbox.latest_data(contains='magnet_msm1_coarse', older_than=older_than, newer_than=newer_than,raise_exc = False)!=False:
        print 'iteration'+str(iteration)


        ## Data location ##

        ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration', older_than=older_than_SSRO)
        
        #msm1 coarse
        timestamp,folder = toolbox.latest_data(contains='magnet_msm1_coarse', older_than=older_than, newer_than=newer_than,return_timestamp = True)
        print 'm folder '+folder
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

        f0mc_temp,u_f0mc_temp = dark_esr_auto_analysis.analyze_dark_esr(None, 2.196*1e-3, add_folder = folder)
         
        #msm1 fine
        timestamp,folder = toolbox.latest_data(contains='magnet_msm1_fine', older_than=older_than, newer_than=newer_than,return_timestamp = True)
        print 'm folder '+folder
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

        f0mf_temp,u_f0mf_temp = dark_esr_auto_analysis.analyze_dark_esr_single(2.02)     

        #msp
        #msm1 coarse
        timestamp,folder = toolbox.latest_data(contains='magnet_msp1_coarse', older_than=older_than, newer_than=newer_than,return_timestamp = True)
        print 'p folder '+folder
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

        f0pc_temp,u_f0pc_temp = dark_esr_auto_analysis.analyze_dark_esr(None, 2.196*1e-3, add_folder = folder)
         
        #msm1 fine
        timestamp,folder = toolbox.latest_data(contains='magnet_msp1_fine', older_than=older_than, newer_than=newer_than,return_timestamp = True)
        print 'm folder '+folder
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

        f0pf_temp,u_f0pf_temp = dark_esr_auto_analysis.analyze_dark_esr_single(3.73) 

        print f0p_temp >0 and f0m_temp >0
        if f0p_temp >0 and f0m_temp >0:
            Bz_measured, Bx_measured = amt.get_B_field(msm1_freq=f0mf_temp*1e9, msp1_freq=f0pf_temp*1e9, u_msm1_freq =u_f0m_temp ,u_msp1_freq=u_f0p_temp)
            f_centrec = (f0mc_temp+f0pc_temp)/2
            f_centre_errorc = np.sqrt(u_f0mc_temp**2+u_f0pc_temp**2)/2
            f_diffc = (f_centrec-ZFS*1e-9)*1e6
            f_diff_errorc = f_centre_errorc*1e6
            f_centref = (f0mf_temp+f0pf_temp)/2
            f_centre_errorf = np.sqrt(u_f0mf_temp**2+u_f0pf_temp**2)/2
            f_difff = (f_centref-ZFS*1e-9)*1e6
            f_diff_errorf = f_centre_errorf*1e6
            f0mc.append(f0m_temp)
            u_f0mc.append(u_f0m_temp)
            f0pc.append(f0p_temp)
            u_f0pc.append(u_f0p_temp)
            f0mf.append(f0m_temp)
            u_f0mf.append(u_f0m_temp)
            f0pf.append(f0p_temp)
            u_f0pf.append(u_f0p_temp)
            f_centrec_list.append(f_centrec)
            f_centrec_error_list.append(f_centre_errorc)
            f_diffc_list.append(f_diffc)
            f_diffc_error_list.append(f_diff_errorc)
            f_centref_list.append(f_centref)
            f_centref_error_list.append(f_centre_errorf)
            f_difff_list.append(f_difff)
            f_diff_errorf_list.append(f_diff_errorf)
            Bx_field_measured.append(Bx_measured)
            Bz_field_measured.append(Bz_measured)
            timestamp_list.append(timestamp)
            it_list.append(iteration)

        older_than = str(int(timestamp)-1)

        iteration = iteration+1
        print iteration

    
    it_list = linspace(0,len(it_list)-1,len(it_list))
    
    np.savez('test_B_field_drift',f0m=f0m,
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

    d = np.load('B_field_drift.npz')
    f0mc = d['f0mc']; u_f0mc = d['u_f0mc']; f0pc = d['f0pc'] ;u_f0pc = d['u_f0pc']
    f0mf = d['f0mf']; u_f0mf = d['u_f0mf']; f0pf = d['f0pf'] ;u_f0pf = d['u_f0pf']
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

    fig = figure(1,figsize=(18,5))
    ax = fig.add_subplot(111)
    ax.errorbar(it_list,f_diff_list,f_diff_error_list)
    ax.errorbar(it_list,f_diff_list,f_diff_error_list)
    ax.set_xlabel('msmt #')
    ax.set_ylabel('relative centrer freq (kHz) (offset 2.87748 GHz)')
    figure(2,figsize=(18,5))
    plt.errorbar(it_list,Bx_field_measured, Bx_error)
    plt.xlabel('msmt #')
    plt.ylabel('measured Bx field (G)')
    plt.savefig('Bx_vs_time',format='png')
    figure(3,figsize=(18,5))
    plt.errorbar(it_list,Bz_diff_list,Bz_error)
    plt.xlabel('msmt #')
    plt.ylabel('measured relative Bz (G) (offset 304.21 G)')
    plt.savefig('Bz_vs_time',format='png')
    figure(4)
    n, bins, patches = plt.hist(f_diff_list,50,normed = 1)
    bincenters = 0.5*(bins[1:]+bins[:-1])
    y = mlab.normpdf( bincenters, mean_f_centre, stdev_f_centre)
    plt.plot(bincenters, y, 'r--', linewidth=1)
    plt.xlabel('binned relative center freq (kHz)')
    plt.title('Mean '+str(mean_f_centre)+' kHz, stdev '+str(stdev_f_centre)+' kHz')
    plt.savefig('binned_freq',format='png')
    figure(5)
    n, bins, patches = plt.hist(Bx_field_measured,50,normed = 1)
    bincenters = 0.5*(bins[1:]+bins[:-1])
    y = mlab.normpdf( bincenters, mean_Bx, stdev_Bx)
    plt.plot(bincenters, y, 'r--', linewidth=1)
    plt.xlabel('binned Bx (G)')
    plt.title('Mean '+str(mean_Bx)+' G, stdev '+str(stdev_Bx)+' G')
    plt.savefig('binned_Bx',format='png')
    figure(6)
    n, bins, patches = plt.hist(Bz_diff_list,50,normed = 1)
    bincenters = 0.5*(bins[1:]+bins[:-1])
    y = mlab.normpdf( bincenters, mean_Bz, stdev_Bz)
    plt.plot(bincenters, y, 'r--', linewidth=1)
    plt.xlabel('binned relative Bz (G) (offset 304.21 G)')
    plt.title('Mean '+str(mean_Bz+304.21)+' G, stdev '+str(stdev_Bz)+' G')
    plt.savefig('binned_Bz',format='png')






fit_B_msmt_loop(older_than = '20140428101106', newer_than = '20140426092819')

plot_meas_B_loop()
