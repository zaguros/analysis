'''Written by MAB 10-3-15 for a general coherence msmt with a "free" exponent'''

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(mbi)

def Carbon_T_mult(timestamp=None, older_than =None, posneg = True, folder_name = 'Hahn', measurement_name = 'adwindata', ssro_calib_timestamp =None,
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            decay_constant = 200, 
            exponent = 2, 
            partstr = 'part', plot_fit = True, do_print = True, fixed = [0,2], show_guess = False):
    ''' 
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''
    figsize=(6,4.7)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        timestamp, folder = toolbox.latest_data(folder_name, older_than=older_than, return_timestamp=True)
    try:
        Number_of_pulses = int(folder[folder.rfind('_N')+2:folder.rfind('_N')+4].replace('_',''))
        Adressed_carbon = int(folder[folder.rfind('_C')+2:folder.rfind('_C')+3])
        print 'N=', Number_of_pulses, '   C=', Adressed_carbon
        DD_msmt = True
    except:
        DD_msmt = False
    
    print 'Timestamp', timestamp
    print 'Folder', folder

    if partstr in folder:
        numberstart = folder.find(partstr)+len(partstr)
        numberofparts = int(folder[numberstart:len(folder)])
        basis_str_pos = folder[folder.rfind('\\')+7:numberstart]
        if posneg:
            posneg_str = 'posneg'
            if 'positive' in basis_str_pos:
                basis_str_neg = basis_str_pos.replace('positive', 'negative')
            else:
                basis_str_neg = basis_str_pos
                basis_str_pos = basis_str_neg.replace('negative', 'positive')
        else:
            if 'positive' in basis_str_pos:
                posneg_str = 'positive'
            else:
                posneg_str = 'negative'
    else:
        numberofparts = 1
        if 'positive' in folder:
            posneg_str = 'positive'
            basis_str_pos = folder[folder.rfind('\\')+7:len(folder)]
            basis_str_neg = basis_str_pos.replace('positive', 'negative')
        else:
            posneg_str = 'negative'
            basis_str_neg = folder[folder.rfind('\\')+7:len(folder)]
            basis_str_pos = basis_str_neg.replace('negative', 'positive')


        # basis_str_pos , basis_str_neg = basis_str_neg , basis_str_pos



    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print ssro_calib_folder

    cum_pts = 0
    print numberofparts
    if posneg:
        
        for kk in range(numberofparts):
            if partstr in folder:
                folder_pos = toolbox.latest_data(basis_str_pos+str(kk+1), older_than = older_than)
                folder_neg = toolbox.latest_data(basis_str_neg+str(kk+1), older_than = older_than)
            else:
                folder_pos = toolbox.latest_data(basis_str_pos, older_than = older_than)
                folder_neg = toolbox.latest_data(basis_str_neg, older_than = older_than)
            a = mbi.MBIAnalysis(folder_pos)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder)
            cum_pts += a.pts

            b = mbi.MBIAnalysis(folder_neg)
            b.get_sweep_pts()
            b.get_readout_results(name='adwindata')
            b.get_electron_ROC(ssro_calib_folder)
            print a.p0, b.p0
            if kk == 0:
                cum_sweep_pts = a.sweep_pts
                cum_p0 = (a.p0+(1-b.p0))/2.
                cum_u_p0 = np.sqrt(a.u_p0**2+b.u_p0**2)/2
                reps_per_datapoint = a.reps
            else:
                cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
                cum_p0 = np.concatenate((cum_p0, (a.p0+(1-b.p0))/2))
                cum_u_p0 = np.concatenate((cum_u_p0, np.sqrt(a.u_p0**2+b.u_p0**2)/2))

        a.pts   = cum_pts
        a.sweep_pts = cum_sweep_pts
        a.p0    = cum_p0
        a.u_p0  = cum_u_p0

    else:
        for kk in range(numberofparts):
            folder = toolbox.latest_data(basis_str_pos+str(kk+1), older_than = older_than)
            a = mbi.MBIAnalysis(folder)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder)
            reps_per_datapoint = a.reps
            cum_pts += a.pts

            if kk == 0:
                cum_sweep_pts = a.sweep_pts
                cum_p0 = a.p0
                cum_u_p0 = a.u_p0
            else:
                cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
                cum_p0 = np.concatenate((cum_p0, a.p0))
                cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

        a.pts   = cum_pts
        a.sweep_pts = cum_sweep_pts
        a.p0    = cum_p0
        a.u_p0  = cum_u_p0


    sorting_order=a.sweep_pts.argsort()
    a.sweep_pts.sort()
    a.p0=a.p0[sorting_order]
    print a.p0
    print a.u_p0
    a.u_p0=a.u_p0[sorting_order]

    ax=a.plot_results_vs_sweepparam(ret='ax',ax=None, 
                                    figsize=figsize, 
                                    ylim=(0.0,1.05)
                                    )
    fit_results = []
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    # print np.shape(a.sweep_pts[:]),np.shape(a.p0[:,0]),np.shape(a.u_p0[:])
    # print np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose()
    ax.plot(x,y)
    if DD_msmt:
        savestr = timestamp + '_C' + str(Adressed_carbon) + '_N' + str(Number_of_pulses) + '_' + posneg_str + '_Pts' + str(cum_pts) + '_Reps' + str(reps_per_datapoint) + '.txt'    
        save_folder_str = 'D:/Dropbox/QEC LT/Decoupling memory/XYdata/' + savestr
        np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)

    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
             x0, decay_constant,exponent)

         #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    ## plot data and fit as function of total time
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    fit_results.append(fit_result)

    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')

    return fit_results


def Carbon_T(timestamp=None, folder_name = 'Hahn', measurement_name = 'adwindata', ssro_calib_timestamp =None,
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            decay_constant = 200, 
            exponent = 2, 
            plot_fit = True, do_print = True, fixed = [2], show_guess = False):
    ''' 
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''
    figsize=(6,4.7)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data(folder_name)

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print ssro_calib_folder

    print folder
    fit_results = []
    
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')

    a.get_electron_ROC(ssro_calib_folder)
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]

    ax=a.plot_results_vs_sweepparam(ret='ax',ax=None, 
                                    figsize=figsize, 
                                    ylim=(0.0,1.0)
                                    )

    
    
    ax.plot(x,y)

    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
             x0, decay_constant,exponent)

         #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    ## plot data and fit as function of total time
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    fit_results.append(fit_result)

    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')

    return fit_results

