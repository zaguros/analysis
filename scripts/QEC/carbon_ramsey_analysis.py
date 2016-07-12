import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(mbi)


def Carbon_Ramsey(timestamp=None, carbon=None, transition=None, measurement_name = ['adwindata'], ssro_calib_timestamp =None,
            frequency = 1, 
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            decay_constant = 200, 
            phase =0, 
            exponent = 2, 
            plot_fit = False, do_print = False, fixed = [2], show_guess = True,
            return_phase = False,
            return_freq = False,
            return_amp = False,
            return_results = True,
            close_plot = False,
            title = 'Carbon'):
    ''' 
    Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    elif carbon != None:
        folder = toolbox.latest_data(contains='C'+str(carbon)+'_ms'+str(transition))
    else:
        folder = toolbox.latest_data(title)

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSROCalibration')

    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
        print ssro_calib_folder


    fit_results = []
    for k in range(0,len(measurement_name)):
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        ax = a.plot_results_vs_sweepparam(ret='ax')

        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]

        ax.plot(x,y)
        p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset, amplitude, 
                x0, decay_constant,exponent,frequency ,phase )

        #plot the initial guess
        if show_guess:
            ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

        print 'fitfunction: '+fitfunc_str

        ## plot data and fit as function of total time
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

        fit_results.append(fit_result)
        if title == None:
            title = 'analyzed_result'
        plt.savefig(os.path.join(folder, title + '.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder, title + '.png'),
        format='png')
        if close_plot == True:
            plt.close()

        # for item in fit_result['params_dict']:
        #     print item

        if return_freq == True:
            f0 = fit_result['params_dict']['f']
            u_f0 = fit_result['error_dict']['f']
            return f0, u_f0

        if return_phase == True and return_amp == False:
            phi0 = fit_result['params_dict']['phi']
            u_phi0 = fit_result['error_dict']['phi']
            return phi0, u_phi0

        if return_phase == True and return_amp == True:
            phi0 = fit_result['params_dict']['phi']
            u_phi0 = fit_result['error_dict']['phi']
            A = fit_result['params_dict']['A']
            u_A = fit_result['error_dict']['A']
            return phi0, u_phi0, A, u_A

        if return_amp == True:
            A = fit_result['params_dict']['A']
            u_A = fit_result['error_dict']['A']
            return A, u_A

    if return_results == True:
        return fit_results

def Carbon_Ramsey_mult_msmts(timestamp=None, measurement_name = ['adwindata'], ssro_calib_timestamp =None,
            frequency = 1, 
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            decay_constant = 200, 
            phase =0, 
            exponent = 2, 
            plot_fit = False, do_print = False, fixed = [2], show_guess = True,
            return_phase = False,
            return_freq = False,
            return_amp = False,
            return_results = True,
            close_plot = False,
            partstr = 'part',
            contains=[],
            title = 'Carbon'):
    ''' 
    Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data(title)
        if partstr in folder:
            numberstart = folder.find(partstr)+len(partstr)
            numberofparts = int(folder[numberstart:len(folder)])
            basis_str = folder[folder.rfind('\\')+7:numberstart]
        else:
            numberofparts = 1

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')

    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
        print ssro_calib_folder

    fit_results = []
    #for kk in range(numberofparts):
    for kk,cnts in enumerate(contains):    
        '''
        if partstr in folder:
            folder = toolbox.latest_data(basis_str+str(kk+1))
        else:
            folder = toolbox.latest_data(basis_str)
        '''
        folder = toolbox.latest_data(cnts)    
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        ax = a.plot_results_vs_sweepparam(ret='ax')

        if kk == 0:
            cum_sweep_pts = a.sweep_pts
            cum_p0 = a.p0
            cum_u_p0 = a.u_p0
            cum_pts = a.pts
        else:
            cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
            cum_p0 = np.concatenate((cum_p0, a.p0))
            cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))
            cum_pts += a.pts

    a.pts   = cum_pts
    a.sweep_pts = cum_sweep_pts
    a.p0    = cum_p0
    a.u_p0  = cum_u_p0

    sorting_order=a.sweep_pts.argsort()
    a.sweep_pts.sort()
    a.p0=a.p0[sorting_order]
    a.u_p0=a.u_p0[sorting_order]

    ax=a.plot_results_vs_sweepparam(ret='ax')

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]

    ax.plot(x,y)

    p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset, amplitude, 
            x0, decay_constant,exponent,frequency ,phase )

    #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    ## plot data and fit as function of total time
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    fit_results.append(fit_result)
    if title == None:
        title = 'analyzed_result'
    plt.savefig(os.path.join(folder, title + '.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, title + '.png'),
    format='png')
    if close_plot == True:
        plt.close()

    if return_freq == True:
        f0 = fit_result['params_dict']['f']
        u_f0 = fit_result['error_dict']['f']
        return f0, u_f0

    if return_phase == True and return_amp == False:
        phi0 = fit_result['params_dict']['phi']
        u_phi0 = fit_result['error_dict']['phi']
        return phi0, u_phi0

    if return_phase == True and return_amp == True:
        phi0 = fit_result['params_dict']['phi']
        u_phi0 = fit_result['error_dict']['phi']
        A = fit_result['params_dict']['A']
        u_A = fit_result['error_dict']['A']
        return phi0, u_phi0, A, u_A

    if return_results == True:
        return fit_results

def Carbon_Ramsey_Crosstalk(timestamp=None, measurement_name = ['adwindata'], ssro_calib_timestamp =None,
            frequency = 1, 
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            decay_constant = 200, 
            phase =0, 
            exponent = 2, 
            plot_fit = False, do_print = False, fixed = [2,3,4], show_guess = True,
            return_phase = False,
            return_freq = False,
            return_results = True,
            return_amp = False,
            close_plot = False,
            title = None):
    ''' 
    Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''

    if timestamp != None:
        folder_a = toolbox.data_from_time(timestamp)
    else:
        folder_a, timestamp = toolbox.latest_data('Crosstalk', return_timestamp = True)
        
    folder_b =  toolbox.latest_data('Crosstalk',older_than = timestamp)   

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print ssro_calib_folder

    fit_results = []
    for k in range(0,len(measurement_name)):
        a = mbi.MBIAnalysis(folder_a)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        ax = a.plot_results_vs_sweepparam(ret='ax')

        X_RO_data = 2*(a.p0.reshape(-1)[:])-1
        X_RO_data_u = 2*(a.u_p0.reshape(-1)[:])


        a = mbi.MBIAnalysis(folder_b)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        ax = a.plot_results_vs_sweepparam(ret='ax')

        x = a.sweep_pts.reshape(-1)[:]
        Y_RO_data = 2*(a.p0.reshape(-1)[:])-1
        Y_RO_data_u = 2*(a.u_p0.reshape(-1)[:])

        RO_data     = (X_RO_data**2 + Y_RO_data**2)**0.5
        RO_data_u   = (1./(X_RO_data**2 + Y_RO_data**2)*(X_RO_data**2 * X_RO_data_u**2 + Y_RO_data**2 *Y_RO_data_u**2))**0.5

        fig = a.default_fig(figsize=(7.5,5))
        ax2 = a.default_ax(fig)
        ax2.axhspan(0,1,fill=False,ls='dotted')
        ax2.set_ylim(-1,1)
        ax2.errorbar(x,RO_data,RO_data_u)

        p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset, amplitude, 
        x0, decay_constant,exponent,frequency ,phase )

        #plot the initial guess
        if show_guess:
            ax2.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

        fit_result = fit.fit1d(x,RO_data, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

        print 'fitfunction: '+fitfunc_str

        ## plot data and fit as function of total time
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax2, plot_data=False)

        fit_results.append(fit_result)
        if title == None:
            title = 'analyzed_result'
        plt.savefig(os.path.join(folder_a, title + '.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder_a, title + '.png'),
        format='png')
        if close_plot == True:
            plt.close()

        if return_freq == True:
            f0 = fit_result['params_dict']['f']
            u_f0 = fit_result['error_dict']['f']
            return f0, u_f0

        if return_phase == True:
            phi0 = fit_result['params_dict']['phi']
            u_phi0 = fit_result['error_dict']['phi']
            Amp = fit_result['params_dict']['A']
            if return_amp == True:    
                return phi0, u_phi0, Amp
            else:
                return phi0, u_phi0

    if return_results == True:
        return fit_results



def Carbon_Ramsey_Crosstalk_no_fit(older_than=None, crosstalk = None,measurement_name = ['adwindata'], 
    ssro_calib_folder =None,title = None):
    ''' 
    Function to analyze 
    Crosstalk is for example ['1to2', '1to5', '2to1', '2to5', '5to2', '5to1'] 
    '''
    #1to5_RO_X

    if crosstalk == None:
        crosstalk = ['1to2', '1to5', '2to1', '2to5', '5to2', '5to1'] 

    for kk in crosstalk:
        
        folder_a = toolbox.latest_data(kk + '_RO_X', older_than = older_than)
        folder_b = toolbox.latest_data(kk + '_RO_Y', older_than = older_than)  

        if ssro_calib_folder == None: 
            ssro_calib_folder = toolbox.latest_data('SSRO')

        for k in range(0,len(measurement_name)):
            a = mbi.MBIAnalysis(folder_a)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder)
            # ax = a.plot_results_vs_sweepparam(ret='ax')

            X_RO_data = 2*(a.p0.reshape(-1)[:])-1
            X_RO_data_u = 2*(a.u_p0.reshape(-1)[:])

            a = mbi.MBIAnalysis(folder_b)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder)
            # ax = a.plot_results_vs_sweepparam(ret='ax')

            # x = a.sweep_pts.reshape(-1)[:]

            Y_RO_data = 2*(a.p0.reshape(-1)[:])-1
            Y_RO_data_u = 2*(a.u_p0.reshape(-1)[:])

            RO_data     = (X_RO_data**2 + Y_RO_data**2)**0.5
            RO_data_u   = (1./(X_RO_data**2 + Y_RO_data**2)*(X_RO_data**2 * X_RO_data_u**2 + Y_RO_data**2 *Y_RO_data_u**2))**0.5
            
            x_ticks = a.sweep_pts.reshape(-1)
            x = range(len(RO_data))

            fig = a.default_fig(figsize=(7.5,5))
            ax2 = a.default_ax(fig)
            ax2.axhspan(0,RO_data[0],fill=False)
            ax2.axhspan(RO_data[0]-RO_data_u[0],RO_data[0]+RO_data_u[0],fill=False,ls='dotted')
            ax2.set_ylim(0,2*np.max(RO_data))
            ax2.set_xlim(x[1]-1, x[-1]+1)
            ax2.errorbar(x[1:],RO_data[1:],RO_data_u[1:])
            ax2.xaxis.set_ticks( x[1:] )
            ax2.set_xticklabels(x_ticks, rotation=90)

      

def Carbon_Ramsey_DD_freq(older_than = None,
            transition = None,  
            carbon = 1,
            frequency = 1, 
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            decay_constant = 200, 
            phase =0, 
            exponent = 2, 
            plot_fit = False, do_print = False, fixed = [2], show_guess = True,
            return_phase = False,
            return_freq = True,
            return_results = True,
            close_plot = False):
    ''' 
    '''
    if transition != None:

        folder1 = toolbox.latest_data(contains = 'evo_times_1_C' + str(carbon)+str(transition), older_than = older_than)
        folder2 = toolbox.latest_data(contains = 'evo_times_2_C' + str(carbon)+str(transition), older_than = older_than)
        folder3 = toolbox.latest_data(contains = 'evo_times_3_C' + str(carbon)+str(transition), older_than = older_than)
    else:
        folder1 = toolbox.latest_data(contains = 'evo_times_1_C' + str(carbon), older_than = older_than)
        folder2 = toolbox.latest_data(contains = 'evo_times_2_C' + str(carbon), older_than = older_than)
        folder3 = toolbox.latest_data(contains = 'evo_times_3_C' + str(carbon), older_than = older_than)




    
    
    a1 = mbi.MBIAnalysis(folder1)
    a1.get_sweep_pts()
    a1.get_readout_results(name='adwindata')
    a1.get_electron_ROC()

    a2 = mbi.MBIAnalysis(folder2)
    a2.get_sweep_pts()
    a2.get_readout_results(name='adwindata')
    a2.get_electron_ROC()

    a3 = mbi.MBIAnalysis(folder3)
    a3.get_sweep_pts()
    a3.get_readout_results(name='adwindata')
    a3.get_electron_ROC()
    '''
    folder4 = toolbox.latest_data(contains = 'evo_times_4_C' + str(carbon), older_than = older_than)
    folder5 = toolbox.latest_data(contains = 'evo_times_5_C' + str(carbon), older_than = older_than)
    folder6 = toolbox.latest_data(contains = 'evo_times_6_C' + str(carbon), older_than = older_than)
    folder7 = toolbox.latest_data(contains = 'evo_times_7_C' + str(carbon), older_than = older_than)
    
    a4 = mbi.MBIAnalysis(folder4)
    a4.get_sweep_pts()
    a4.get_readout_results(name='adwindata')
    a4.get_electron_ROC()

    a5 = mbi.MBIAnalysis(folder5)
    a5.get_sweep_pts()
    a5.get_readout_results(name='adwindata')
    a5.get_electron_ROC()

    a6 = mbi.MBIAnalysis(folder6)
    a6.get_sweep_pts()
    a6.get_readout_results(name='adwindata')
    a6.get_electron_ROC()

    a7 = mbi.MBIAnalysis(folder7)
    a7.get_sweep_pts()
    a7.get_readout_results(name='adwindata')
    a7.get_electron_ROC()
    
    a1.p0           = np.r_[a1.p0,a2.p0,a3.p0,a4.p0,a5.p0,a6.p0,a7.p0]
    a1.u_p0         = np.r_[a1.u_p0,a2.u_p0,a3.u_p0,a4.u_p0,a5.u_p0,a6.u_p0,a7.u_p0]
    a1.sweep_pts     = np.r_[a1.sweep_pts,a2.sweep_pts,a3.sweep_pts,a4.sweep_pts,a5.sweep_pts,a6.sweep_pts,a7.sweep_pts]
    '''

    a1.p0           = np.r_[a1.p0,a2.p0,a3.p0]
    a1.u_p0         = np.r_[a1.u_p0,a2.u_p0,a3.u_p0]
    a1.sweep_pts     = np.r_[a1.sweep_pts,a2.sweep_pts,a3.sweep_pts]
    

    x = a1.sweep_pts.reshape(-1)[:]
    y = a1.p0.reshape(-1)[:]
 
    ax = a1.plot_results_vs_sweepparam(ret='ax')

    ax.plot(x,y)
    p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset, amplitude, 
            x0, decay_constant,exponent,frequency ,phase )

    #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    print 'fitfunction: '+fitfunc_str

    ## plot data and fit as function of total time
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

   
    title = 'DD freq ramsey C' +str(carbon)
    plt.savefig(os.path.join(folder1, title + '.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder1, title + '.png'),
    format='png')

