import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(fit)


def OneQubitTomo(timestamp=None, measurement_name = ['adwindata'], ssro_calib_timestamp =None,
        frequency = [1,1,1], offset =[ 0.5,0.5,0.5], amplitude =[ 0.5,0.5,0.5],  phase =[0.5,0.5,0.5],
        fixed = [],
        plot_fit = False, do_print = False, show_guess = True):
    '''
    Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed')
    [freq, offset, Amplitude, phase]
    '''

    if timestamp != None:
        timestampZ = timestamp
        folderZ = toolbox.data_from_time(timestamp)
    else:
        timestampZ, folderZ   = toolbox.latest_data('CarbonR',return_timestamp =True)
    timestampY, folderY = toolbox.latest_data('CarbonR',older_than = timestampZ, return_timestamp =True)
    timestampX, folderX = toolbox.latest_data('CarbonR',older_than = timestampY, return_timestamp =True)
    folders = [folderX,folderY, folderZ]

    if ssro_calib_timestamp == None:
        ssro_calib_folder = toolbox.latest_data('SSRO',older_than = timestampX)
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'

    order_of_bases = ['X', 'Y','Z']
    x = [None]*len(folders)
    y = [None]*len(folders)
    fit_result = [None]*len(folders)
    for i,folder in enumerate(folders):
        fit_results = []

        print '*'*60
        print order_of_bases[i] + ' Tomography'
        print 'folder %s' %folder


        for k in range(0,len(measurement_name)):
            a = mbi.MBIAnalysis(folder)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder)
            if i ==0:
                ax = a.plot_results_vs_sweepparam(ret='ax',labels = order_of_bases[i])
                # a.plot_results_vs_sweepparam(ax=ax)
            else:
                a.plot_results_vs_sweepparam(ax=ax,labels = order_of_bases[i])
            x[i] = a.sweep_pts.reshape(-1)[:]
            y[i]= a.p0.reshape(-1)[:]


            print frequency[i]
            p0, fitfunc, fitfunc_str = common.fit_cos(frequency[i], offset[i], amplitude [i],phase[i] )
            try:
                fit_result[i] = fit.fit1d(x[i],y[i], None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=fixed)
                if plot_fit == True:
                    plot.plot_fit1d(fit_result[i], np.linspace(x[i][0],x[i][-1],201), ax=ax,
                            plot_data=False,print_info = False)
                fit.write_to_file(fit_result[i],folder,fitname = str(order_of_bases[i])+'-tomography')
            except:
                pass
            if show_guess:
                ax.plot(np.linspace(x[i][0],x[i][-1],201), fitfunc(np.linspace(x[i][0],x[i][-1],201)), ':', lw=2)


    print 'fitfunction: '+fitfunc_str
    if plot_fit ==True:
        ax.legend(('X data','X-fit','Y data','Y-fit','Z data','Z-fit'),fontsize='x-small')
    elif plot_fit == False:
        ax.legend(('X data','Y data','Z data'),fontsize='small')

    ## plot data and fit as function of total time

    fit_results.append(fit_result[i])

    plt.savefig(os.path.join(folder, 'analyzed_tomography_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_tomography_result.png'),
    format='png')

    return fit_results
