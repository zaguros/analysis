import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import sequence, mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot


guess_offset = 1
guess_x0 = 3.730
guess_splitB = 30.
guess_splitN = 2.18e-3
# guess_splitC = .8e-3 #12.78
guess_width = 0.2e-3
guess_splitB = 30.
guess_splitN = 2.18e-3
guess_sigma = 0.2e-3
guess_amplitude = 0.3

# try fitting
guess_A_min1 = 0.3
guess_A_plus1 = 0.3
guess_A_0 = 0.3
guess_Nsplit = 2.196e-3


def analyze_2_desr(timestamps=[None,None], measurement_name = ['DarkESR'], ssro_calib_timestamps =[None,None],MBI = [False,False],
       center_guess = False, ax=None, ret=None,min_dip_depth = 0.85 ,
        fixed = [],
        plot_fit = False, do_print = False, show_guess = True, print_info = True,
        figsize = (3,2), linewidth = 2, markersize = 2, fontsize =10,
        title =None,savename ='DarkESR_2',**kw):
    '''
    Same basic function as normal analyze_dark_esr except that it makes two plots that share axes.
    '''
    f , (ax1,ax2) = plt.subplots(2, sharex=True, sharey=True,figsize = figsize)
    folder = [None,None]
    ssro_calib_folder = [None,None]
    a =[None,None]
    fit_result = [None,None]
    x = np.arange(5)
    y = x*x

    # ax1.plot(x,y)
    # ax2.errorbar(x,x*y,yerr=x)

    for i in range(len(timestamps)): #Get data folders


        if timestamps[i] == None:
            timestamps[i], folder[i]   = toolbox.latest_data(measurement_name[0],return_timestamp=True)
        else:
            folder[i] = toolbox.data_from_time(timestamps[i])

        if ssro_calib_timestamps[i] == None:
            ssro_calib_folder[i] = toolbox.latest_data('SSRO')
        else:
            ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamps[i])
            ssro_calib_folder[i] = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'

        ssro_calib_folder[i] = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration')

        if MBI[i] == False:
            a[i] = sequence.SequenceAnalysis(folder[i])
            a[i].get_sweep_pts()
            a[i].get_readout_results('ssro')
        else:
            print 'printing folder'
            print folder[i]
            a[i] = mbi.MBIAnalysis(folder[i])
            a[i].get_sweep_pts()
            a[i].sweep_pts = a[i].sweep_pts *1e-3
            a[i].get_readout_results(name = 'adwindata')
        a[i].get_electron_ROC(ssro_calib_folder=ssro_calib_folder[i])

        x = a[i].sweep_pts # convert to MHz
        y = a[i].p0.reshape(-1)[:]

        if center_guess == True:
            guess_ctr = float(raw_input('Center guess?'))
        else:
            j=0
            print min_dip_depth
            print y[21]
            while y[j]>min_dip_depth and j < len(y)-2:  #y[j]>0.93*y[j+1]: # such that we account for noise
                k = j
                j += 1
            #j = len(y)-2
            if k > len(y)-5:
                print 'Could not find dip'
                return
            else:
                print 'k'+str(k)
                print len(y)
                guess_ctr = x[k]+ guess_splitN #convert to GHz and go to middle dip
                print 'guess_ctr= '+str(guess_ctr)

        ### fitfunction
        if MBI ==True:
            ### fitfunction
            A_min1 = fit.Parameter(guess_A_min1, 'A_min1')
            A_plus1 = fit.Parameter(guess_A_plus1, 'A_plus1')
            A_0 = fit.Parameter(guess_A_0, 'A_0')
            o = fit.Parameter(guess_offset, 'o')
            x0 = fit.Parameter(guess_x0, 'x0')
            sigma = fit.Parameter(guess_sigma, 'sigma')
            Nsplit = fit.Parameter(guess_Nsplit, 'Nsplit')

            def fitfunc(x):
                return o() - A_min1()*np.exp(-((x-(x0()-Nsplit()))/sigma())**2) \
                        - A_plus1()*np.exp(-((x-(x0()+Nsplit()))/sigma())**2) \
                        - A_0()*np.exp(-((x-x0())/sigma())**2) \

            fit_result[i] = fit.fit1d(x, y, None, p0 = [A_min1, A_plus1, A_0, sigma, o, x0],
                    fitfunc = fitfunc, do_print=True, ret=True, fixed=[])
            print 'test'

        else: #MBI ==False:


            A_min1 = fit.Parameter(guess_A_min1, 'A_min1')
            A_plus1 = fit.Parameter(guess_A_plus1, 'A_plus1')
            A_0 = fit.Parameter(guess_A_0, 'A_0')
            o = fit.Parameter(guess_offset, 'o')
            x0 = fit.Parameter(guess_x0, 'x0')
            sigma = fit.Parameter(guess_sigma, 'sigma')
            Nsplit = fit.Parameter(guess_Nsplit, 'Nsplit')
            def fitfunc(x):
                return o() - A_min1()*np.exp(-((x-(x0()-Nsplit()))/sigma())**2) \
                        - A_plus1()*np.exp(-((x-(x0()+Nsplit()))/sigma())**2) \
                        - A_0()*np.exp(-((x-x0())/sigma())**2) \

            try:
                fit_result[i] = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
                        guess_amplitude, guess_width, guess_ctr,
                        (3, guess_splitN),
                        do_print=True, ret=True, fixed=[])
                # fit_result[i] = fit.fit1d(x, y, None, p0 = [A_min1, A_plus1, A_0, sigma, o, x0],
                        # fitfunc = fitfunc, do_print=True, ret=True, fixed=[])


            except Exception:
                print 'Exception'
                guess_ctr = float(raw_input('Center guess?'))
                fit_result[i] = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
                        guess_amplitude, guess_width, guess_ctr,
                        (3, guess_splitN),
                        do_print=True, ret=True, fixed=[4])

    #Plotting results
        #Both ax =1 and 2


    a[0].plot_result_vs_sweepparam(ret=ret, name='ssro', ax=ax1)
        # ax1.errorbar(a[i].sweep_pts, a[i].p0, fmt='o',
        #         yerr=a[i].u_p0,markersize=markersize)

    plot.plot_fit1d(fit_result[0], np.linspace(min(x), max(x), 1000), ax=ax1, plot_data=False, print_info = print_info,**kw)

    a[1].plot_results_vs_sweepparam(ret=ret, name='adwindata', ax=ax2)
        # ax2.errorbar(a[i].sweep_pts, a[i].p0, fmt='o',
        #         yerr=a[i].u_p0,markersize=markersize)

    plot.plot_fit1d(fit_result[1], np.linspace(min(x), max(x), 1000), ax=ax2, plot_data=False, print_info = print_info,**kw)


    ax2.set_xlabel('Microwave frequency (GHz)',fontsize = fontsize)
    ax2.set_ylabel(r'$F$ $\left( |0\rangle \right)$', fontsize = fontsize)
    if title ==None:
        ax2.set_title(a[0].timestamp+'\n'+a[0].measurementstring)
    else:
        ax2.set_title(title)


    f.subplots_adjust(hspace=0)
    plt.setp([ax.get_xticklabels() for ax in f.axes[:-1]], visible=False)
    print 'folder again'
    print folder[0]

    print savename+'.pdf'
    plt.savefig(os.path.join(folder[1], savename+'.pdf'),
    format='pdf',bbox_inches='tight')
    print 'Saved file in %s' %folder[1]


    return fit_result
