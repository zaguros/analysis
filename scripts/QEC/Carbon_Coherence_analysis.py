'''Written by MAB 10-3-15 for a general coherence msmt with a "free" exponent'''

import numpy as np
import os, sys
if os.name == 'posix':
    sys.path.append("/Users/"+os.getlogin()+"/Documents/teamdiamond/")
else:
    sys.path.append("/measuring/")
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(mbi)
reload(fit)

def Carbon_T_mult_and_averaging(timestamp=None, older_than =None, posneg = True, folder_name = 'Hahn', measurement_name = 'adwindata',ssro_calib_timestamp=None,
            offset = 0.5,
            save_data = True,
            save_name = None,
            save_with_folder_name=False,
            x0 = 0, 
            fancy_plot = False,
            amplitude = 0.5,
            fit_func = 'exp',  
            decay_constant = 200, 
            exponent = 1.75,
            frequency = 0.8,
            slope = 0.,
            zfill = 1,
            axes = None, return_axes = False, reppartstr = '_reppart',
            FETpartstr = '_FETpart', plot_fit = True, do_print = True, fixed = [0,2], show_guess = False):
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
    
    try:
        

        Tomo = folder[folder.find('Tomo')+4:folder.find('Tomo')+6]
        Adressed_carbon = folder[folder.rfind('_C')+2:folder.rfind('_C')+5]
        print   'C=', Adressed_carbon
        print   'Tomo=', Tomo
        if Tomo == folder[folder.find('Tomo'):folder.find('Tomo')+4]:
            Tomo_msmt = True
        else:
            Tomo_msmt = False
    except:
        Tomo_msmt = False

    numberofparts = 1
    FETnumberofparts = 1
    repnumberofparts = 1
    print 'Timestamp', timestamp
    print 'Folder', folder

    if FETpartstr in folder:
        FETnumberstart = folder.rfind(FETpartstr)+len(FETpartstr)
        FETnumberofparts = int(folder[FETnumberstart:len(folder)])
        numberstart = FETnumberstart

    if reppartstr in folder:
        repnumberstart = folder.rfind(reppartstr)+len(reppartstr)
        numberstart = repnumberstart
        if FETpartstr in folder:
            repnumberofparts = int(folder[repnumberstart:folder.rfind(FETpartstr)])
        else:
            repnumberofparts = int(folder[repnumberstart:len(folder)])

    if reppartstr in folder or reppartstr in folder:
        basis_str_pos = folder[folder.rfind('\\')+7:numberstart]
        # print 'Ja', basis_str_pos
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
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
        print ssro_calib_folder

    cum_pts = 0
    timestamp_list = []
    print numberofparts
    if posneg:
        for kk in range(numberofparts):
            if partstr in folder:
                timestamp, folder_pos = toolbox.latest_data(basis_str_pos+str(kk+1).zfill(zfill), older_than = older_than, return_timestamp=True)
                timestamp_list.append(timestamp)
                # print folder_pos
                timestamp, folder_neg = toolbox.latest_data(basis_str_neg+str(kk+1).zfill(zfill), older_than = older_than, return_timestamp=True)
                timestamp_list.append(timestamp)
                # print folder_neg
            else:
                timestamp, folder_pos = toolbox.latest_data(basis_str_pos, older_than = older_than,  return_timestamp=True)
                timestamp_list.append(timestamp)
                timestamp, folder_neg = toolbox.latest_data(basis_str_neg, older_than = older_than, return_timestamp=True)
                timestamp_list.append(timestamp)
            
            a = mbi.MBIAnalysis(folder_pos)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder)
            cum_pts += a.pts
            print a.sweep_pts
            b = mbi.MBIAnalysis(folder_neg)
            b.get_sweep_pts()
            b.get_readout_results(name='adwindata')
            b.get_electron_ROC(ssro_calib_folder)
            if a.pts == 1:
                a.sweep_pts = [a.sweep_pts[0]*len(a.sweep_pts)]
                b.sweep_pts = [b.sweep_pts[0]*len(b.sweep_pts)]
                print a.sweep_pts
            if False:
                a.p0 = (1.-a.p0)
                b.p0 = (1.-b.p0)
            # print a.p0, b.p0
            if kk == 0:
                cum_sweep_pts = a.sweep_pts
                cum_p0 = (a.p0+(1-b.p0))/2.
                cum_u_p0 = np.sqrt(a.u_p0**2+b.u_p0**2)/2
                reps_per_datapoint = a.reps
            else:
                cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
                cum_p0 = np.concatenate((cum_p0, (a.p0+(1-b.p0))/2))
                cum_u_p0 = np.concatenate((cum_u_p0, np.sqrt(a.u_p0**2+b.u_p0**2)/2))
            # print a.get_sweep_pts
        a.pts   = cum_pts
        a.sweep_pts = cum_sweep_pts
        a.p0    = cum_p0
        # print 
        # print cum_sweep_pts
        # print
        # print cum_p0
        # print
        a.u_p0  = cum_u_p0

    else:
        for kk in range(FETnumberofparts):
            for ll in range(repnumberofparts):
                timestamp, folder = toolbox.latest_data(basis_str_pos+str(ll+1).zfill(zfill)+FETpartstr+ str(kk+1).zfill(zfill), older_than = older_than, return_timestamp=True)
                timestamp_list.append(timestamp)
                a = mbi.MBIAnalysis(folder)
                a.get_sweep_pts()
                a.get_readout_results(name='adwindata')
                print a.u_normalized_ssro
                if ll == 0:
                    reps_per_datapoint = a.reps
                    norm_ssro = a.normalized_ssro
                else:
                    norm_ssro += a.normalized_ssro
                    reps_per_datapoint += a.reps
                
                if ll == repnumberofparts-1:
                    a.reps = reps_per_datapoint
                    a.normalized_ssro = norm_ssro / repnumberofparts
                    a.u_normalized_ssro = (a.normalized_ssro*(1.-a.normalized_ssro)/a.reps)**0.5
                    print 'final'
                    print a.u_normalized_ssro
                    print 'final'
            a.get_electron_ROC(ssro_calib_folder)
            cum_pts += a.pts
            if a.pts == 1:
                a.sweep_pts = [a.sweep_pts[0]*len(a.sweep_pts)]
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
    # print a.p0
    # print a.u_p0
    a.u_p0=a.u_p0[sorting_order]
    print 'plot1'
    if not fancy_plot:
        if axes == None:
            ax=a.plot_results_vs_sweepparam(ret='ax',ax=None, 
                                            figsize=figsize, 
                                            ylim=(0.0,1.05)
                                            ) 
        else:
            ax=a.plot_results_vs_sweepparam(ret='ax',ax=axes, 
                                            figsize=figsize, 
                                            ylim=(0.0,1.05)
                                            )
    x = a.sweep_pts.reshape(-1)[:]
    fit_results = []
    y = a.p0.reshape(-1)[:]
    y_err = a.u_p0.reshape(-1)[:]
    # print np.shape(a.sweep_pts[:]),np.shape(a.p0[:,0]),np.shape(a.u_p0[:])
    # print np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose()
    
    print Tomo_msmt
    print DD_msmt
    if save_data:
        if save_with_folder_name == False:
            if Tomo_msmt and save_data:
                savestr = timestamp + '_Tomo' + Tomo + '_C' + str(Adressed_carbon) + '_N' + str(Number_of_pulses) + '_' + posneg_str + '_Pts' + str(cum_pts) + '_Reps' + str(reps_per_datapoint) + '.txt'    
                save_folder_str = 'D:/Dropbox/QEC LT/Decoupling memory/MultiCarbon_Tomo_msmt/' + savestr
                np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)
            elif DD_msmt and save_data:
                savestr = timestamp + '_C' + str(Adressed_carbon) + '_N' + str(Number_of_pulses) + '_' + posneg_str + '_Pts' + str(cum_pts) + '_Reps' + str(reps_per_datapoint) + '.txt'    
                save_folder_str = 'D:/Dropbox/QEC LT/Decoupling memory/XYdata/' + savestr
                np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)
        else:
            savestr = 'From' + str(min(map(int,timestamp_list))) + '_To' + str(max(map(int,timestamp_list))) + basis_str_pos + '.txt'
            save_folder_str = r'D:/Dropbox/QEC LT/Decoupling memory/General_Data/' + savestr
            print save_folder_str
            np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)



    if fit_func == 'line': 
        fixed = [1]
        p0, fitfunc, fitfunc_str = common.fit_line(offset, slope)
    elif fit_func == 'exp':
        p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
            x0, decay_constant,exponent)
    elif fit_func == 'decaying cosine':
        fixed = [0,2,6]
        phase = 0.
        # frequency = 
        p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset, amplitude, 
            x0, decay_constant, exponent,frequency,phase)

    
    
    # fixed=[]
         #plot the initial guess

    

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    # print 'zeropoint', (0.5-fit_result['params_dict']['a'])/(fit_result['params_dict']['b'])
    ## plot data and fit as function of total time

    if fancy_plot:
        xmax = 60
        if axes==None:
            fig,ax = plt.subplots(figsize=(7,5))
            color='b'
            label = r'$m_s=0$'
        else:
            color='g'
            label = r'$m_s=-1$'
            ax = axes
        plot.plot_fit1d(fit_result, np.linspace(x[0],xmax,1001), ax=ax, plot_data=False,add_txt = False,linewidth =2, linestyle = '-',color = color)
        errlines = ax.errorbar(x,y,yerr = y_err,ls = '', marker = 'o',markersize = 8,capsize=5, lw = 2,color = color, label=label)
        ax.set_ylim([-0.1,1.1])
        xticks = [int(0)]+(np.arange(10,xmax+10,10)).tolist()
        plt.xticks(xticks)
        ax.set_xticklabels(xticks)
        ax.set_xlabel('Relaxation time (s)',fontsize = 15)
        ax.set_ylabel('Fidelity',fontsize = 15)
        #ax.hlines([0.5],x[0]-1e-3,x[-1]+1e3,linestyles='dotted',linewidth = 2)
        ax.hlines([1.],0,xmax,linestyles='dotted',linewidth = 2)
        ax.hlines([0.],0,xmax,linestyles='dotted',linewidth = 2)
        ax.set_xlim([x[0],xmax])
        plt.yticks([0,0.2,0.4,0.6,0.8,1])

        ax.tick_params(axis='x', which='major', labelsize=15)
        ax.tick_params(axis='y', which='major', labelsize=15)
        plt.rcParams['axes.linewidth'] = 2
        ax.tick_params('both', length=4, width=2, which='major')
        plt.legend(loc = (0,0.2),fontsize = 18,numpoints = 1,frameon = False,columnspacing=0.5,handletextpad=0.0)


        plt.savefig(r'D:\Dropbox\QEC LT\Decoupling memory\00_Thesis_plots\electronshort.pdf')        
    else:
        ax.plot(x,y)
        ax.hlines([0.5],x[0],x[-1],linestyles='dotted',linewidth = 2)
        if show_guess:
            ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

        if plot_fit == True:
            print fit_result
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    fit_results.append(fit_result)

    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')
    if return_axes == True:
        return ax, fit_result
    else:
        return fit_results

def Carbon_T_mult(timestamp=None, older_than =None, posneg = True, folder_name = 'Hahn', measurement_name = 'adwindata',ssro_calib_timestamp=None,
            offset = 0.5,
            save_data = True,
            save_name = None,
            save_with_folder_name=False,
            x0 = 0, 
            fancy_plot = False,
            amplitude = 0.5,
            slope=0.,
            fit_func = 'exp',  
            decay_constant = 200, 
            exponent = 1.75,
            frequency = 0.8,
            zfill = 1,
            axes = None, return_axes = False,
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
    
    try:
        

        Tomo = folder[folder.find('Tomo')+4:folder.find('Tomo')+6]
        Adressed_carbon = folder[folder.rfind('_C')+2:folder.rfind('_C')+5]
        print   'C=', Adressed_carbon
        print   'Tomo=', Tomo
        if Tomo == folder[folder.find('Tomo'):folder.find('Tomo')+4]:
            Tomo_msmt = True
        else:
            Tomo_msmt = False
    except:
        Tomo_msmt = False

    print 'Timestamp', timestamp
    print 'Folder', folder

    if partstr in folder:
        numberstart = folder.rfind(partstr)+len(partstr)
        numberofparts = int(folder[numberstart:len(folder)])
        basis_str_pos = folder[folder.rfind('\\')+7:numberstart]
        # print 'Ja', basis_str_pos
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

    # numberofparts = 6

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
        print ssro_calib_folder

    cum_pts = 0
    timestamp_list = []
    print numberofparts
    if posneg:
        for kk in range(numberofparts):
            if partstr in folder:
                timestamp, folder_pos = toolbox.latest_data(basis_str_pos+str(kk+1).zfill(zfill), older_than = older_than, return_timestamp=True)
                timestamp_list.append(timestamp)
                # print folder_pos
                timestamp, folder_neg = toolbox.latest_data(basis_str_neg+str(kk+1).zfill(zfill), older_than = older_than, return_timestamp=True)
                timestamp_list.append(timestamp)
                # print folder_neg
            else:
                timestamp, folder_pos = toolbox.latest_data(basis_str_pos, older_than = older_than,  return_timestamp=True)
                timestamp_list.append(timestamp)
                timestamp, folder_neg = toolbox.latest_data(basis_str_neg, older_than = older_than, return_timestamp=True)
                timestamp_list.append(timestamp)
            
            a = mbi.MBIAnalysis(folder_pos)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder)
            cum_pts += a.pts
            print a.sweep_pts
            b = mbi.MBIAnalysis(folder_neg)
            b.get_sweep_pts()
            b.get_readout_results(name='adwindata')
            b.get_electron_ROC(ssro_calib_folder)
            if a.pts == 1:
                a.sweep_pts = [a.sweep_pts[0]*len(a.sweep_pts)]
                b.sweep_pts = [b.sweep_pts[0]*len(b.sweep_pts)]
                print a.sweep_pts
            print 'Yaha!'
            if False:
                a.p0 = (1.-a.p0)
                b.p0 = (1.-b.p0)
            # print a.p0, b.p0
            if kk == 0:
                cum_sweep_pts = a.sweep_pts
                cum_p0 = (a.p0+(1-b.p0))/2.
                cum_u_p0 = np.sqrt(a.u_p0**2+b.u_p0**2)/2
                reps_per_datapoint = a.reps
            else:
                cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
                cum_p0 = np.concatenate((cum_p0, (a.p0+(1-b.p0))/2))
                cum_u_p0 = np.concatenate((cum_u_p0, np.sqrt(a.u_p0**2+b.u_p0**2)/2))
            # print a.get_sweep_pts
        a.pts   = cum_pts
        a.sweep_pts = cum_sweep_pts
        a.p0    = cum_p0
        # print 
        # print cum_sweep_pts
        # print
        # print cum_p0
        # print
        a.u_p0  = cum_u_p0

    else:
        for kk in range(numberofparts):
            timestamp, folder = toolbox.latest_data(basis_str_pos+str(kk+1).zfill(zfill), older_than = older_than, return_timestamp=True)
            timestamp_list.append(timestamp)
            a = mbi.MBIAnalysis(folder)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder)
            reps_per_datapoint = a.reps
            cum_pts += a.pts
            if a.pts == 1:
                a.sweep_pts = [a.sweep_pts]
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
    print sorting_order
    # print a.p0
    # print a.u_p0
    a.u_p0=a.u_p0[sorting_order]
    print 'plot1'
    if not fancy_plot:
        if axes == None:
            ax=a.plot_results_vs_sweepparam(ret='ax',ax=None, 
                                            figsize=figsize, 
                                            ylim=(0.0,1.05)
                                            ) 
        else:
            ax=a.plot_results_vs_sweepparam(ret='ax',ax=axes, 
                                            figsize=figsize, 
                                            ylim=(0.0,1.05)
                                            )
    x = a.sweep_pts.reshape(-1)[:]
    print x
    fit_results = []
    y = a.p0.reshape(-1)[:]
    y_err = a.u_p0.reshape(-1)[:]
    # firstpterr = np.sqrt(np.sum(y_err[0:4]**2.))/4.
    # firstpv = np.sum(y[0:4])/4.
    # y_err = y_err[3::]
    # y_err[3] = firstpterr
    # x = x[3::]
    # y = y[3::]
    # y[3]=firstpv
    # print x
    # print np.shape(a.sweep_pts[:]),np.shape(a.p0[:,0]),np.shape(a.u_p0[:])
    # print np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose()
    
    print Tomo_msmt
    print DD_msmt
    if save_data:
        if save_with_folder_name == False:
            if Tomo_msmt and save_data:
                savestr = timestamp + '_Tomo' + Tomo + '_C' + str(Adressed_carbon) + '_N' + str(Number_of_pulses) + '_' + posneg_str + '_Pts' + str(cum_pts) + '_Reps' + str(reps_per_datapoint) + '.txt'    
                save_folder_str = 'D:/Dropbox/QEC LT/Decoupling memory/MultiCarbon_Tomo_msmt/' + savestr
                np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)
            elif DD_msmt and save_data:
                savestr = timestamp + '_C' + str(Adressed_carbon) + '_N' + str(Number_of_pulses) + '_' + posneg_str + '_Pts' + str(cum_pts) + '_Reps' + str(reps_per_datapoint) + '.txt'    
                save_folder_str = 'D:/Dropbox/QEC LT/Decoupling memory/XYdata/' + savestr
                np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)
        else:
            savestr = 'From' + str(min(map(int,timestamp_list))) + '_To' + str(max(map(int,timestamp_list))) + basis_str_pos + '.txt'
            save_folder_str = r'D:/Dropbox/QEC LT/Decoupling memory/General_Data/' + savestr
            print save_folder_str
            np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)



    if fit_func == 'line':
        fixed = []
        p0, fitfunc, fitfunc_str = common.fit_line(offset, slope)
    elif fit_func == 'exp':
        p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
            x0, decay_constant,exponent)
    elif fit_func == 'decaying cosine':
        fixed = [0,2,6]
        phase = 0.
        # frequency = 
        p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset, amplitude, 
            x0, decay_constant, exponent,frequency,phase)

    
    
    # fixed=[]
         #plot the initial guess

    

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    # print 'zeropoint', (0.5-fit_result['params_dict']['a'])/(fit_result['params_dict']['b'])
    ## plot data and fit as function of total time

    if fancy_plot:
        # xmax = 0.75
        if axes==None:
            fig,ax = plt.subplots(figsize=(7,5))
            color='b'
        else:
            color='g'
            ax = axes
        plot.plot_fit1d(fit_result, np.linspace(0,60,10001), ax=ax, plot_data=False,add_txt = False,linewidth =2, linestyle = '-',color = color)
        errlines = ax.errorbar(x,y,yerr = y_err,ls = '', marker = 'o',markersize = 8,capsize=5, lw = 2,color=color)
        ax.set_ylim([-0.1,1.1])
        # xticks = [int(0)]+(np.arange(0.25,0.76,0.25)).tolist()
        # plt.xticks([0.5,0.75,1])
        # ax.set_xticklabels(xticks)
        # circle1=plt.Circle((0.305,0.5),.2,color='r')
        # ax.add_artist(circle1)
        # ax.arrow(0.305, 0.75, 0., -0.15, head_width=0.02, head_length=0.04, fc='k', ec='k')
        ax.set_xlabel('Free evolution time (s)',fontsize = 15)
        # ax.set_ylabel(r'$\langle $IX$ \rangle$',fontsize = 20)
        ax.set_ylabel('Fidelity',fontsize = 15)
        x[-1]=0.15
        #ax.hlines([0.5],0,x[-1]+1e3,linestyles='dotted',linewidth = 2)
        ax.hlines([1.],xvals[0]-1,xvals[-1]+1,linestyles='dotted',linewidth = 2)
        ax.hlines([0.],xvals[0]-1,xvals[-1]+1,linestyles='dotted',linewidth = 2)
        ax.set_xlim([0,x[-1]])
        ax.set_ylim([0.45,1.05])
        plt.yticks([0,0.2,0.4,0.6,0.8,1])
        ax.tick_params(axis='x', which='major', labelsize=15)
        ax.tick_params(axis='y', which='major', labelsize=15)
        plt.rcParams['axes.linewidth'] = 2
        ax.tick_params('both', length=4, width=2, which='major')
        # plt.savefig(r'D:\Dropbox\QEC LT\Decoupling memory\00_Thesis_plots\C2T2ms02.eps')

    else:
        ax.plot(x,y)
        ax.hlines([0.5],x[0],x[-1],linestyles='dotted',linewidth = 2)
        ax.set_ylim([-0.1,1.1])
        if show_guess:
            print 'Hello'
            ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)


        if plot_fit == True:
            print fit_result
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    fit_results.append(fit_result)

    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')
    if return_axes == True:
        return ax, fit_result
    else:
        return fit_results

def Carbon_DFS(older_than =None, posneg = True, folder_name = 'NuclearDD_111_1_sil18_sweep_evolution_time_auto_C3&6_LogicpX_TomoZZ_ROpositive_el1_N01', measurement_name = 'adwindata',ssro_calib_timestamp=None,
            offset = 0.5,
            save_data = False,
            save_name = None,
            save_with_folder_name=False,
            x0 = 0, 
            fancy_plot = False,
            amplitude = 0.5,
            fit_func = 'exp',  
            decay_constant = 200, 
            exponent = 1.75,
            frequency = 0.8,
            zfill = 2,
            axes = None, return_axes = False,
            partstr = 'part', plot_fit = True, do_print = True, fixed = [], show_guess = False):
    ''' 
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''
    figsize=(6,4.7)

    
    timestamp, folder = toolbox.latest_data(folder_name, older_than=older_than, return_timestamp=True)
    
   
    print 'Timestamp', timestamp
    print 'Folder', folder

    
    numberstart = folder.rfind(partstr)+len(partstr)
    numberofparts = int(folder[numberstart:len(folder)])
    basis_str_pos = folder[folder.rfind('\\')+7:numberstart]
        # print 'Ja', basis_str_pos
    posneg_str = 'posneg'
    

        # basis_str_pos , basis_str_neg = basis_str_neg , basis_str_pos


    ssro_calib_timestamp = '20150504_130503'
    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
        print ssro_calib_folder

    cum_pts = 0
    timestamp_list = []
    color = ['b','r']
    print numberofparts
    for aaaa, Logic in enumerate(['pX','mX']):
        logicstart = basis_str_pos.rfind('Logic')+len('Logic')
        basis_str_pos = basis_str_pos[0:logicstart] + Logic + basis_str_pos[logicstart+2:] 
        tomo_p0 = []
        tomo_u_p0 = []
        T_s = []
        A_s = []
        n_s = []
        for tomo_basis in ['XX','YY','ZZ']:
            tomostart = basis_str_pos.rfind('Tomo')+len('Tomo')
            basis_str_pos = basis_str_pos[0:tomostart] + tomo_basis + basis_str_pos[tomostart+2:] 
            if 'positive' in basis_str_pos:
                basis_str_neg = basis_str_pos.replace('positive', 'negative')
            else:
                basis_str_neg = basis_str_pos
                basis_str_pos = basis_str_neg.replace('negative', 'positive')

            for kk in range(numberofparts):
                if partstr in folder:
                    timestamp, folder_pos = toolbox.latest_data(basis_str_pos+str(kk+1).zfill(zfill), older_than = older_than, return_timestamp=True)
                    timestamp_list.append(timestamp)
                    print folder_pos
                    timestamp, folder_neg = toolbox.latest_data(basis_str_neg+str(kk+1).zfill(zfill), older_than = older_than, return_timestamp=True)
                    timestamp_list.append(timestamp)
                    # print folder_neg
                else:
                    timestamp, folder_pos = toolbox.latest_data(basis_str_pos, older_than = older_than,  return_timestamp=True)
                    timestamp_list.append(timestamp)
                    timestamp, folder_neg = toolbox.latest_data(basis_str_neg, older_than = older_than, return_timestamp=True)
                    timestamp_list.append(timestamp)
                
                a = mbi.MBIAnalysis(folder_pos)
                a.get_sweep_pts()
                a.get_readout_results(name='adwindata')
                a.get_electron_ROC(ssro_calib_folder)
                cum_pts += a.pts
                # print a.sweep_pts
                b = mbi.MBIAnalysis(folder_neg)
                b.get_sweep_pts()
                b.get_readout_results(name='adwindata')
                b.get_electron_ROC(ssro_calib_folder)
                if a.pts == 1:
                    a.sweep_pts = [a.sweep_pts[0]*len(a.sweep_pts)]
                    b.sweep_pts = [b.sweep_pts[0]*len(b.sweep_pts)]
                    print a.sweep_pts
                if (tomo_basis == 'YY' and Logic == 'pX') or (tomo_basis == 'ZZ' and Logic == 'mX'):
                    a.p0 = (1.-a.p0)
                    b.p0 = (1.-b.p0)

                # print a.p0, b.p0
                if kk == 0:
                    cum_sweep_pts = a.sweep_pts
                    cum_p0 = (a.p0+(1-b.p0))/2.
                    cum_u_p0 = np.sqrt(a.u_p0**2+b.u_p0**2)/2
                    reps_per_datapoint = a.reps
                else:
                    cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
                    cum_p0 = np.concatenate((cum_p0, (a.p0+(1-b.p0))/2))
                    cum_u_p0 = np.concatenate((cum_u_p0, np.sqrt(a.u_p0**2+b.u_p0**2)/2))
                # print a.get_sweep_pts
            p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, x0, decay_constant,exponent)
            if tomo_basis == 'ZZ':
                slope = 0
                p0, fitfunc, fitfunc_str = common.fit_line(offset, slope)

            
            a.pts   = cum_pts
            print 'pts'
            print a.pts
            a.sweep_pts = cum_sweep_pts
            a.p0    = cum_p0
            a.u_p0  = cum_u_p0
            sorting_order=a.sweep_pts.argsort()
            print a.sweep_pts
            a.sweep_pts.sort()
            print a.sweep_pts
            print a.p0
            print sorting_order
            a.p0=a.p0[sorting_order]
            print a.p0
            print a.u_p0
            a.u_p0=a.u_p0[sorting_order]
            tomo_p0.append(a.p0)
            tomo_u_p0.append(a.u_p0)

            x = a.sweep_pts.reshape(-1)[:]
            fit_results = []
            y = a.p0.reshape(-1)[:]
            y_err = a.u_p0.reshape(-1)[:]
            # print tomo_basis
            # print Logic
            # print 'x'
            # print x
            # print 'y'
            # print y
            # print fitfunc
            # print p0
            if tomo_basis == 'ZZ':
                fixed = [1]
            else:
                fixed = [0,2,4]
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
            
            print 'HIERO----------------------------------------'
            if tomo_basis == 'ZZ':
                amp = fit_result['params_dict']['a']
                print amp
            else:
                A_s.append(fit_result['params_dict']['A'])
                print A_s
                T_s.append(fit_result['params_dict']['T'])
                if 4 in fixed:
                    n_s.append(exponent)
                else:
                    n_s.append(fit_result['params_dict']['n'])

        # print 'Tomos'
        # print tomo_p0[0], tomo_p0[0]*2.-1. 
        # print tomo_p0[1]*2.-1.
        # print tomo_p0[2]*2.-1.
        if True:
            cum_p0 = ((1+tomo_p0[0]*2.-1 + tomo_p0[1]*2.-1. + tomo_p0[2]*2.-1)/4.)
            cum_u_p0  = ( np.sqrt( (tomo_u_p0[0]*2.)**2. + (tomo_u_p0[1]*2.)**2. + (tomo_u_p0[2]*2.)**2. ))/4.
            a.p0 = cum_p0
            a.u_p0 = cum_u_p0

        print 'cumulative'
        print cum_p0
        print cum_u_p0
        # a.pts   = cum_pts
        # a.sweep_pts = cum_sweep_pts
        # a.p0    = cum_p0
        # # print 
        # # print cum_sweep_pts
        # # print
        # # print cum_p0
        # # print
        # a.u_p0  = cum_u_p0
        

        # sorting_order=a.sweep_pts.argsort()
        # a.sweep_pts.sort()
        # a.p0=a.p0[sorting_order]
        # print a.p0
        # print a.u_p0
        #a.u_p0=a.u_p0[sorting_order]
        print 'plot1'
        if not fancy_plot:
            if axes == None:
                ax=a.plot_results_vs_sweepparam(ret='ax',ax=None,fmt='o', 
                                                figsize=figsize, 
                                                ylim=(0.0,1.05)
                                                ) 
                axes = ax
            else:
                ax=a.plot_results_vs_sweepparam(ret='ax',ax=axes, 
                                                figsize=figsize, 
                                                ylim=(0.0,1.05)
                                                )
        x = a.sweep_pts.reshape(-1)[:]
        fit_results = []
        y = a.p0.reshape(-1)[:]
        y_err = a.u_p0.reshape(-1)[:]
        # x = a.sweep_pts
        # y = a.p0
        # y_err = a.u_p0
        # print np.shape(a.sweep_pts[:]),np.shape(a.p0[:,0]),np.shape(a.u_p0[:])
        # print np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose()
        
        # print Tomo_msmt
        # print DD_msmt
        if save_data:
            if save_with_folder_name == False:
                if Tomo_msmt and save_data:
                    savestr = timestamp + '_Tomo' + Tomo + '_C' + str(Adressed_carbon) + '_N' + str(Number_of_pulses) + '_' + posneg_str + '_Pts' + str(cum_pts) + '_Reps' + str(reps_per_datapoint) + '.txt'    
                    save_folder_str = 'D:/Dropbox/QEC LT/Decoupling memory/MultiCarbon_Tomo_msmt/' + savestr
                    np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)
                elif DD_msmt and save_data:
                    savestr = timestamp + '_C' + str(Adressed_carbon) + '_N' + str(Number_of_pulses) + '_' + posneg_str + '_Pts' + str(cum_pts) + '_Reps' + str(reps_per_datapoint) + '.txt'    
                    save_folder_str = 'D:/Dropbox/QEC LT/Decoupling memory/XYdata/' + savestr
                    np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)
            else:
                savestr = 'From' + str(min(map(int,timestamp_list))) + '_To' + str(max(map(int,timestamp_list))) + basis_str_pos + '.txt'
                save_folder_str = r'D:/Dropbox/QEC LT/Decoupling memory/General_Data/' + savestr
                print save_folder_str
                np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)



        if fit_func == 'line':
            fixed = [1]
            amplitude = 0.
            p0, fitfunc, fitfunc_str = common.fit_line(offset, amplitude)
        elif fit_func == 'exp':
            p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
                x0, decay_constant,exponent)
        elif fit_func == 'decaying cosine':
            fixed = [0,2,6]
            phase = 0.
            # frequency = 
            p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset, amplitude, 
                x0, decay_constant, exponent,frequency,phase)

        
        
        # fixed=[]
             #plot the initial guess

        

        # fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

        # print 'zeropoint', (0.5-fit_result['params_dict']['a'])/(fit_result['params_dict']['b'])
        ## plot data and fit as function of total time

        if fancy_plot:

            xmax = 0.75
            if axes==None:
                fig,ax = plt.subplots(figsize=(7,5))
                axes = ax
            else:
                ax = axes
            xarray = np.linspace(x[0],x[-1],201)
            def state_F_function(x,A,T,n,a):
                y1 = 0.5+A[0]*np.exp(-(x/T[0])**n[0])
                y2 = 0.5+A[1]*np.exp(-(x/T[1])**n[1])
                print  ((1+y1*2.-1 + y2*2.-1. + a*2.-1)/4.)
                return (1+y1*2.-1 + y2*2.-1. + a*2.-1)/4.
            ax.plot(xarray,state_F_function(xarray,A_s,T_s,n_s,amp),linewidth=2,color=color[aaaa])
            # plot.plot_fit1d(fit_result, np.linspace(x[0],0.75,1001), ax=ax, plot_data=False,add_txt = False,linewidth =2, linestyle = '-',color = '0.25')
            errlines = ax.errorbar(x,y,yerr = y_err, color = color[aaaa],ls = '', marker = 'o',markersize = 8,capsize=5, lw = 2)
            
            # xticks = [int(0)]+(np.arange(0.25,0.76,0.25)).tolist()
            # plt.xticks(xticks)
            # ax.set_xticklabels(xticks)
            ax.set_xlabel('Free evolution time (s)',fontsize = 15)
            ax.set_ylabel('Fidelity',fontsize = 15)
            ax.hlines([0.5],x[0]-1e-3,x[-1]+1e3,linestyles='dotted',linewidth = 2)
            ax.set_xlim([x[0],x[-1]])
            # ax.set_ylim([0,1])
            plt.xticks(np.linspace(0,1,3))
            yticks = np.linspace(0.3,0.7,5)
            plt.yticks(yticks)
            ax.tick_params(axis='x', which='major', labelsize=15)
            ax.tick_params(axis='y', which='major', labelsize=15)
            plt.rcParams['axes.linewidth'] = 2
            ax.tick_params('both', length=4, width=2, which='major')
            ax.set_ylim([0.25,0.75])
            plt.savefig(r'D:\Dropbox\QEC LT\Decoupling memory\00_Thesis_plots\StateF.pdf',bbox_inches='tight')

        else:
            ax.plot(x,y)
            ax.hlines([0.5],x[0],x[-1],linestyles='dotted',linewidth = 2)
            if show_guess:
                ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

            if plot_fit == True:
                print fit_result
                x = np.linspace(x[0],x[-1],201)
                def state_F_function(x,A,T,n,a):
                    y1 = 0.5+A[0]*np.exp(-(x/T[0])**n[0])
                    y2 = 0.5+A[1]*np.exp(-(x/T[1])**n[1])
                    print  ((1+y1*2.-1 + y2*2.-1. + a*2.-1)/4.)
                    return (1+y1*2.-1 + y2*2.-1. + a*2.-1)/4.
                ax.plot(x,state_F_function(x,A_s,T_s,n_s,amp))
                # plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

        fit_results.append(fit_result)

    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')
    if return_axes == True:
        return ax, fit_result
    else:
        return fit_results


def Carbon_T_averaging(timestamp=None, older_than =None, posneg = True, folder_name = 'Hahn', measurement_name = 'adwindata', ssro_calib_timestamp =None,
            offset = 0.3, 
            x0 = 0,  
            amplitude = 0.7,  
            decay_constant = 10, 
            exponent = 1, 
            partstr = 'part', plot_fit = True, do_print = True, fixed = [2,4], show_guess = False):
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

    print 'Timestamp', timestamp
    print 'Folder', folder

    numberstart = folder.rfind('reppart') + 7
    numberofparts = 12
    basis_str_pos = folder[folder.rfind('\\')+7:numberstart]
        
  

        # basis_str_pos , basis_str_neg = basis_str_neg , basis_str_pos



    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print ssro_calib_folder

    cum_pts = 0
    print numberofparts
    
    for kk in range(numberofparts):
        if len(str(kk+1)) == 2:
            strprt = str(kk+1)
        else:
            strprt = str(kk+1) + '_' 
        folder = toolbox.latest_data(basis_str_pos+strprt, older_than = older_than)
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        reps_per_datapoint = a.reps
        cum_pts += a.pts

        if kk == 0:
            # cum_sweep_pts = a.sweep_pts
            cum_p0 = np.zeros((1,len(a.sweep_pts)))
            cum_u_p0 = np.zeros((1, len(a.sweep_pts)))
    
        # cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
        # print a.p0
        # print cum_p0[0,:]
        # print a.u_p0
        cum_p0[0,:] += a.p0.reshape(-1)
        cum_u_p0[0,:] += np.power(a.u_p0.reshape(-1),2)
    
    # sumvar = np.zeros((1,len(a.p0)))
    # for ii in range(numberofparts):
    #     sumvar += np.power(cum_p0[ii,:],2)  
    TEMP = np.sqrt(cum_u_p0)/numberofparts
    a.u_p0 = TEMP.transpose()
    # np.sqrt(a.u_p0**2+b.u_p0**2)/2
    a.pts   = cum_pts
    # a.sweep_pts = cum_sweep_pts
    TEMP = cum_p0 / numberofparts
    a.p0 = TEMP.transpose()
    # a.u_p0  = cum_u_p0


    # sorting_order=a.sweep_pts.argsort()
    # a.sweep_pts.sort()
    # a.p0=a.p0[sorting_order]
    # print a.p0
    # print a.u_p0
    # a.u_p0=a.u_p0[sorting_order]

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


    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
             x0, decay_constant,exponent)

    # p0, fitfunc, fitfunc_str = common.fit_line(offset, amplitude)
    # fixed=[]
         #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)
    ax.hlines([0.5],0,1.4,linestyles='dotted',linewidth = 2)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    # print 'zeropoint', (0.5-fit_result['params_dict']['a'])/(fit_result['params_dict']['b'])
    ## plot data and fit as function of total time
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    fit_results.append(fit_result)

    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')

    return fit_results


def Carbon_T_mult_save(timestamp=None, older_than =None, posneg = True, folder_name = 'Hahn', measurement_name = 'adwindata', ssro_calib_timestamp =None,
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

    Number_of_pulses = int(folder[folder.rfind('_N')+2:folder.rfind('_N')+4].replace('_',''))
    Adressed_carbon = int(folder[folder.rfind('_C')+2:folder.rfind('_C')+3])
    print 'N=', Number_of_pulses, '  C=', Adressed_carbon
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
            basis_str_pos = folder[folder.rfind('\\')+7:len(folder)]
            basis_str_neg = basis_str_pos.replace('positive', 'negative')
        else:
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
            cum_pts += a.pts

            if kk == 0:
                cum_sweep_pts = a.sweep_pts
                cum_p0 = a.p0
                cum_u_p0 = a.u_p0
                reps_per_datapoint = a.reps
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

    # ax=a.plot_results_vs_sweepparam(ret='ax',ax=None, 
    #                                 figsize=figsize, 
    #                                 ylim=(0.0,1.05)
    #                                 )
    # fit_results = []
    # x = a.sweep_pts.reshape(-1)[:]
    # y = a.p0.reshape(-1)[:]
    
    # ax.plot(x,y)

    # p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
    #          x0, decay_constant,exponent)
    savestr = timestamp + '_C' + Adressed_carbon + '_N' + Number_of_pulses + '_' + posneg_str + '_Pts' + cum_pts + '_Reps' + reps_per_datapoint + '.txt'
    
    np.savetxt("/Users/"+os.getlogin()+"/Documents/teamdiamond/XY4data/" + savestr, np.concatenate((a.sweep_pts,a.p0,a.u_p0),comments='#'+savestr))
         #plot the initial guess
    # if show_guess:
    #     ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

    # fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    ## plot data and fit as function of total time
    # if plot_fit == True:
    #     plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    # fit_results.append(fit_result)

    # plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    # format='pdf')
    # plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    # format='png')

    # return fit_results

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
    print a.normalized_ssro
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

def Carbon_assymetric_RO(timestamps=None, timestamps_neg=None ,older_than =None, posneg = False, folder_name = 'Hahn', measurement_name = 'adwindata', ssro_calib_timestamp =None,
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            sigma = 5, 
            partstr = 'part', plot_fit = True, do_print = True, fixed = [0], show_guess = False):
    ''' 
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''
    figsize=(6,4.7)
    Number_of_folders = len(timestamps)
    folders = []
    folders_neg = []
    Tomo_names = []
    for tstmp_nr, timestamp in enumerate(timestamps):
        folder = toolbox.data_from_time(timestamp)
        folders.append(folder)
        Tomo_names.append(folder[folder.rfind('Tomo')+4:folder.rfind('Tomo')+6])
        if timestamps_neg != None:
            folder =  toolbox.data_from_time(timestamps_neg[tstmp_nr])
            folders_neg.append(folder)
        print 'Tomo=', Tomo_names
    DD_msmt = False
    
    print 'Timestamp', timestamp
    print 'Folder', folder

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print ssro_calib_folder
    colors = ['b','g','r']

    fit_results = []

    for kk in range(Number_of_folders):
        
        a = mbi.MBIAnalysis(folders[kk])
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        reps_per_datapoint = a.reps
        
        if timestamps_neg != None:
            b = mbi.MBIAnalysis(folders_neg[kk])
            b.get_sweep_pts()
            b.get_readout_results(name='adwindata')
            b.get_electron_ROC(ssro_calib_folder)

            a.p0    = (a.p0+(1-b.p0))/2.
            a.u_p0  = np.sqrt(a.u_p0**2+b.u_p0**2)/2

        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]

        
        # ax.plot(x,y)
        
        p0, fitfunc, fitfunc_str = common.fit_gauss(offset, amplitude, 
             x0, sigma)
        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
        print 'T2*= %.2f +- %.2f' % (fit_result['params_dict']['sigma'],fit_result['error_dict']['sigma'])
        label = Tomo_names[kk]+', x0= %.2f +- %.2f' % (2.**.5*fit_result['params_dict']['x0'],2.**.5*fit_result['error_dict']['x0'])
        if kk==0:
            ax=a.plot_results_vs_sweepparam(ret='ax',ax=None, 
                                            figsize=figsize, 
                                            ylim=(0.4,1.0),
                                            labels=[label],
                                            )
            ax.set_xlim((-15.,15.))
        else:
           a.plot_results_vs_sweepparam(ax=ax, 
                                       figsize=figsize, 
                                       ylim=(0.4,1.),
                                       labels=[label]
                                       ) 
        
       
    ## plot data and fit as function of total time
        print fit_result['params_dict']['x0']
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False,color=colors[kk], add_txt=False)
        # ax.set_label(Tomo_names[kk])
        fit_results.append(fit_result)
    plt.legend(loc='lower')
    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')

    return fit_results