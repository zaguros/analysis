import os, sys
import numpy as np
import h5py
import logging
import sympy

from matplotlib import pyplot as plt

from analysis.lib import fitting; 
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, esr
from analysis.lib.tools import plot
from analysis.lib.math import error


def convert_tstamp_to_minutes(tstamp):
    ### there are probably better ways to do this....

    d = int(tstamp[6:8])*60*24
    h = int(tstamp[8:10])*60
    mi = int(tstamp[10:12])
    s = float(tstamp[12:14])/60.
    # print y,m,d,h,mi,s
    return d+h+mi+s

### settings
timestamp = '20160110_113400'#'20150609_231010' ### one can use this timestamp as starting point.
ssro_calib_timestamp = '20160109_102710'
analysis_length = 804#11#804 # number of dark esr measurements under consideration.


average_data = True

show_guess = False


ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Pippin_SIL1'

folder_list = []
tstamp_list = []

if timestamp == None:
    timestamp = '20800101_101010' ## some future date... only younger data is considered.

for i in range(analysis_length):
    timestamp, folder = toolbox.latest_data(contains='DarkESR',
                                        return_timestamp =True,
                                        older_than=timestamp,
                                        raise_exc=False)
    folder_list.append(folder)
    tstamp_list.append(timestamp)



# print folder_list

guess_offset = 1.00
guess_A_0 = 0.08
guess_x0 = 1713.468427
guess_sigma = 250e-3
guess_splitC = 0.4#2.182 #


###################################################

# get data

###################################################

cum_pts = 0
cum_sweep_pts       = np.empty(0)
cum_p0              = np.empty(0)
cum_u_p0              = np.empty(0)
fits = []

cum_normalized_ssro = np.empty(0)

for kk,folder in enumerate(folder_list):

    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)
    # cum_pts += a.pts
    a.sweep_pts = a.sweep_pts*1e3
    x = a.sweep_pts # convert to MHz
    y = a.p0.reshape(-1)[:]

    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_A_0, guess_sigma, guess_x0,
            # (2, guess_splitN),
             (2, guess_splitC),
            do_print=False, ret=True, fixed=[])

    fits.append(fit_result)

    # print fits


    if kk == 0:
        cum_sweep_pts = a.sweep_pts
        if average_data:
            cum_p0 = a.p0/float(analysis_length)
            cum_u_p0 = a.u_p0**2
        else:
            cum_p0 = a.p0
            cum_u_p0 = a.u_p0

        # reps_per_datapoint = a.reps
    else:
        if average_data:
            cum_p0 = a.p0/float(analysis_length) + cum_p0
            cum_u_p0 = a.u_p0**2+ cum_u_p0
        else:
            cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
            cum_p0 = np.concatenate((cum_p0, a.p0))
            cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))
        
if average_data:
    cum_u_p0 = np.sqrt(cum_u_p0)/float(analysis_length)
sorting_order=cum_sweep_pts.argsort()
cum_sweep_pts.sort()
cum_p0=cum_p0[sorting_order]
cum_u_p0=cum_u_p0[sorting_order]

## print the result of fit to the averaged data
fit_result = fit.fit1d(x, cum_p0[sorting_order], esr.fit_ESR_gauss, guess_offset,
            guess_A_0, guess_sigma, guess_x0,
            # (2, guess_splitN),
             (2, guess_splitC),
            do_print=True, ret=True, fixed=[])


# a.pts   = cum_pts
a.sweep_pts = cum_sweep_pts
a.p0    = cum_p0
a.u_p0  = cum_u_p0

minute_list = []
for tstamp in tstamp_list:
    minute_list.append(convert_tstamp_to_minutes(tstamp))

minute_list = np.array(minute_list)
minute_list = minute_list-minute_list[-1]

### extract duration. assumes that all measurements are spaced equally
x =minute_list
y = []
y_err = []

for res in fits:
    y.append(res['params_dict']['x0'])
    y_err.append(res['error_dict']['x0'])

y = np.array(y)
y_err = np.array(y_err)
# #########################################################

# # plot

# #########################################################



ax3 = a.plot_result_vs_sweepparam(ret='ax',name='ssro',fmt='-')
plot.plot_fit1d(fit_result, np.linspace(1713,1714,1001), ax=ax3, plot_data=False,add_txt = False)

fig = plt.figure()
ax = plt.subplot()

# n, bins, patches = plt.hist(y,13, facecolor='g')
plt.errorbar(x,y,y_err)
centre = 1.713e3+0.42
plt.ylim([centre-0.05,centre+0.05])
plt.xlabel('Time (minutes)')
plt.ylabel('Central frequency (MHz)')


# ax2 = plt.subplot()
# p0, fitfunc, fitfunc_str = common.fit_gauss(0.0, 14., 0,20)
# print len(n),len(bins)
# bins_rescaled = [(x-1.713e3+0.42)*1e3 for x in bins]
# fit_result = fit.fit1d(bins_rescaled[:-1],n, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0])

# plot.plot_fit1d(fit_result, np.linspace(-35,20,101), ax=ax, plot_data=True,color = 'r',add_txt = True, lw = 1)
# plt.xlabel("Detuning from central frequency (kHz)")
# plt.ylabel("# of occurences")

# plt.close('all')
plt.show()
plt.close('all')
