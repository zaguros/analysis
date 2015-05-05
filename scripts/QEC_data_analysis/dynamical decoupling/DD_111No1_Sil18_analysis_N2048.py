'''Script to analyze the dynamical decoupling data
by THT'''

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)


### Inputs ###

## Data location ##
measurement_name = ['adwindata']

timestamp = []

#### loop the get the timestmaps for 2048. this can also be done with a list of newer and older than timestamps. (if more measurements are to be used.)

new_tsmp = '20150430_104800' ## newer than
old_tsmp = '20150430_121207' ## older than
search_string = '_DecouplingSequence_111_1_sil18sweep_tau_N_2048'

while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)
    timestamp.append(old_tsmp)


## fit paramaters initial values (general exponential function) ##
offset      = 0.5
amplitude   = 0.45
position    = 0
T2          = 80
power       = 1.5

## other settings ##
plot_fit    = True
show_guess  = False

fit_results = []

cum_pts = 0
cum_sweep_pts       = np.empty(0)
cum_p0              = np.empty(0)
cum_u_p0              = np.empty(0)

cum_normalized_ssro = np.empty(0)
#for k in range(0,len(measurement_name)):

for kk in range(len(timestamp)):
    folder = toolbox.data_from_time(timestamp[kk])

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
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

ax = a.plot_results_vs_sweepparam(ret='ax')

x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)

#plot the initial guess
if show_guess:
    ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), '-', lw=2)

fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])

## plot data and fit as function of total time
if plot_fit == True:
    plot.plot_fit1d(fit_result, np.linspace(0,np.amax(x),201), ax=ax, plot_data=False)

fit_results.append(fit_result)

print folder
plt.savefig(os.path.join(folder, 'combined_result.pdf'),
format='pdf')
plt.savefig(os.path.join(folder, 'combined_result.png'),
format='png')
