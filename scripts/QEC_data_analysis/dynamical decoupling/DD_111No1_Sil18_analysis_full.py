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

timestamp128 = ['20150429_110632', '20150429_110915', '20150429_111413','20150430_103829','20150430_104617']

timestamp256 = ['20150429_113126', '20150429_113332', '20150429_113704','20150429_114159',
                '20150430_095434','20150430_100149']

timestamp512 = ['20150429_115420', '20150429_115558', '20150429_115839',
                '20150429_120219','20150429_120701','20150430_100255','20150430_100834','20150430_101511']

timestamp1024 = ['20150429_120802', '20150429_120904', '20150429_121029','20150429_121220',
                '20150429_121427','20150429_121659','20150429_121945','20150430_101738',
                '20150430_102105','20150430_102455','20150430_102908','20150430_103340']

timestamp2048 = []

timestamplist = [timestamp128,timestamp256,timestamp512,timestamp1024]
labels = ['128','256','512','1024']
color_list =  ['b','r','g','c']
## fit paramaters initial values (general exponential function) ##
offset      = 0.5
amplitude   = 0.45
position    = 0
T2          = 60
power       = 2.

## other settings ##
plot_fit    = True
show_guess  = False

fit_results = []

cum_pts = 0
cum_sweep_pts       = np.empty(0)
cum_p0              = np.empty(0)
cum_u_p0              = np.empty(0)

cum_normalized_ssro = np.empty(0)

# fig=plt.figure()
# ax=plt.subplot()
#for k in range(0,len(measurement_name)):
for ii, timestamp in enumerate(timestamplist):
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
    if ii == 0:
        ax = a.plot_results_vs_sweepparam(ret='ax',labels=[labels[ii]])
    else:
        plt.errorbar(cum_sweep_pts,cum_p0,cum_u_p0,fmt='o',label=labels[ii],color=color_list[ii])

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    y_u = cum_u_p0
    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)

    #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), '-', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2,4])

    ## plot data and fit as function of total time
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax,add_txt = False, plot_data=False,color=color_list[ii])

    fit_results.append(fit_result)

    cum_pts = 0
    cum_sweep_pts       = np.empty(0)
    cum_p0              = np.empty(0)
    cum_u_p0              = np.empty(0)

    cum_normalized_ssro = np.empty(0)

print folder

ax.set_ylim(0.4,1.02)
ax.set_xlim(1.,410.)
ax.set_xscale('log')
plt.legend()


plt.savefig(os.path.join(folder, 'combined_result.pdf'),
format='pdf')
plt.savefig(os.path.join(folder, 'combined_result.png'),
format='png')
