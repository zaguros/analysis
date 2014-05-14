'''
Script to analyze the dynamical decoupling data
'''

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import analysis.lib.QEC.nuclear_spin_characterisation as SC #used for simulating FP response of spins
# import magnettools as mt # Does not work atm because of qt lab being imported in MT
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from analysis.lib.tools import plot

from analysis.lib.m2.ssro import sequence
import matplotlib.mlab as mlab
from tempfile import TemporaryFile

reload(toolbox)


def fingerprint_loop(older_than = None, newer_than = None,contains = '',number = 0,filename = 'loop'):
    fig = figure(number,figsize=(10,5))
    ax = fig.add_subplot(111)
    histdata_1521 = []
    histdata_1531 = []
    # ssro_calib_folders = ['d:\\measuring\\data\\20140504\\070533_AdwinSSRO_SSROCalibration_Hans_sil1']

    ## Data location ##
    while toolbox.latest_data(contains=contains, older_than=older_than, newer_than=newer_than,raise_exc = False)!=False:
        timestamp,folder = toolbox.latest_data(contains=contains, older_than=older_than, newer_than=newer_than,return_timestamp = True)
        ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration', older_than='20140502193810')
        ## Data location ##
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
       
    # ############
    # ## Plotting ###
    # # ############
    #     print folder
        x = a.sweep_pts.reshape(-1)[:]
        y = a.p0.reshape(-1)[:]
        y_error = a.u_p0.reshape(-1)[:]
        # ax.errorbar(x,y,y_error)
        older_than = str(int(timestamp)-1)


        idx_1521 = np.argmin(abs(x-15.21))
        idx_1531 = np.argmin(abs(x-15.31))

        if y[idx_1521]<0.8:
          histdata_1521.append(y[idx_1521])
          histdata_1531.append(y[idx_1531])
          ax.errorbar(x,y,y_error)
        print y_error
#     print os.path.join(folder, 'fingerprint_'+filename+'.pdf')
#     plt.savefig(os.path.join(folder, 'fingerprint_'+filename+'.pdf'),
#         format='pdf')
#     plt.savefig(os.path.join(folder, 'fingerprint_'+filename+'.png'),
#         format='png')

#     mean_hist_1531 = np.mean(histdata_1531)
#     stdev_histdata_1531 = np.std(histdata_1531)
#     print 'dip 15.31'
#     print mean_hist_1531
#     print stdev_histdata_1531

#     figure(number*10)
#     n, bins, patches = plt.hist(histdata_1531,50,normed = 1)
#     bincenters = 0.5*(bins[1:]+bins[:-1])
#     y = mlab.normpdf( bincenters, mean_hist_1531, stdev_histdata_1531)
#     plt.plot(bincenters, y, 'r--', linewidth=1)
#     plt.xlabel('binned relative depth of dip')
#     plt.title('Histogram of dip at tau = 15.31 us')
#     plt.savefig('binned_data'+filename,format='png')

#     mean_hist_1521 = np.mean(histdata_1521)
#     stdev_histdata_1521 = np.std(histdata_1521)
#     print 'dip 15.21'
#     print mean_hist_1521
#     print stdev_histdata_1521

#     figure(number*11)
#     n, bins, patches = plt.hist(histdata_1521,50,normed = 1)
#     bincenters = 0.5*(bins[1:]+bins[:-1])
#     y = mlab.normpdf( bincenters, mean_hist_1521, stdev_histdata_1521)
#     plt.plot(bincenters, y, 'r--', linewidth=1)
#     plt.xlabel('binned relative depth of dip')
#     plt.title('Histogram of dip at tau = 15.21 us')
#     plt.savefig('binned_data15_2'+filename,format='png')

# # fingerprint_loop(older_than = '20140505090835', newer_than = '20140502190202',contains='around18',number = 1,filename = 'loop_total')

# fingerprint_loop(older_than = '20140505090835', newer_than = '20140503215120',contains='around18',number = 2,filename = 'loop_1')

# fingerprint_loop(older_than = '20140503215120', newer_than = '20140502190202',contains='around18',number = 3,filename = 'loop_1')

fingerprint_loop(older_than = '20140505090835', newer_than = '201405005090000',contains='around15',number = 1,filename = 'loop_total')

# fingerprint_loop(older_than = '20140505090835', newer_than = '20140503215120',contains='around15',number = 2,filename = 'loop_1')

# fingerprint_loop(older_than = '20140503215120', newer_than = '20140502190202',contains='around15',number = 3,filename = 'loop_2')





