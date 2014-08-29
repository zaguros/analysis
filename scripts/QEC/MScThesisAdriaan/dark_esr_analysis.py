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

### settings
timestamp =None#'20140710_205010' #' #'114103_PulsarD' #YYYYmmddHHMMSS
#Timestamp = 20140808_110902
#SSRO Timestamp = 20140808_110615


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


def analyze_dark_esr(timestamp=None, measurement_name = ['DarkESR'], ssro_calib_timestamp =None,
           center_guess = False, ax=None, ret=None,min_dip_depth = 0.85 ,
            fixed = [],
            plot_fit = False, do_print = False, show_guess = True, print_info = True,
            figsize = (3,2), linewidth = 2, markersize = 2, fontsize =10,
            title =None,savename ='DarkESR',**kw):

    if timestamp == None:
        timestamp, folder   = toolbox.latest_data(measurement_name[0],return_timestamp =True)
    else:
        folder = toolbox.data_from_time(timestamp)

    if ssro_calib_timestamp == None:
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'

    ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration')
    print ssro_calib_folder
    a = sequence.SequenceAnalysis(folder)
    fig = a.default_fig(figsize=figsize)
    ax = a.default_ax(fig)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)
    x = a.sweep_pts # convert to MHz
    y = a.p0.reshape(-1)[:]
    # ax.plot(x,y)
    a.plot_result_vs_sweepparam(ret=ret, name='ssro', ax=ax)
    #ax.set_ylim(0.1,1.05)


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
        # fit_result = fit.fit1d(x, y, None, p0 = [A_min1, A_plus1, A_0, sigma, o, x0],
        # fitfunc = fitfunc, do_print=True, ret=True, fixed=[])

        fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
                guess_amplitude, guess_width, guess_ctr,
                # (2, guess_splitN),
                # (2, guess_splitC),
                # (2, guess_splitB),
                (3, guess_splitN),
                do_print=True, ret=True, fixed=[])

        plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), ax=ax, plot_data=False, print_info = print_info,**kw)

        Norm=(fit_result['params'][0]+fit_result['params'][1]+fit_result['params'][2])
        Population_left=fit_result['params'][0]/Norm
        Population_middle=fit_result['params'][2]/Norm
        Population_right=fit_result['params'][1]/Norm
        print '############################'
        print 'Population left ' , Population_left
        print 'Population middle ' , Population_middle
        print 'Population right ' , Population_right
        print '#############################'
    except Exception:
        guess_ctr = float(raw_input('Center guess?'))
        fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
                guess_amplitude, guess_width, guess_ctr,
                (3, guess_splitN),
                do_print=True, ret=True, fixed=[4])
        plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), ax=ax, plot_data=False, **kw)

    ax.set_xlabel('Microwave frequency (GHz)',fontsize = fontsize)
    ax.set_ylabel(r'$F$ $\left( |0\rangle \right)$', fontsize = fontsize)
    if title ==None:
        ax.set_title(a.timestamp+'\n'+a.measurementstring)
    else:
        ax.set_title(title)

    plt.savefig(os.path.join(folder, savename+'.pdf'),
    format='pdf',bbox_inches='tight')


    return fit_result
