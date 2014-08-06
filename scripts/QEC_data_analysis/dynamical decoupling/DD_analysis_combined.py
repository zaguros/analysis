'''Script to analyze the dynamical decoupling data
by THT'''

import numpy as np
import os
from measurement.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt

reload(common)


show_fit = True

### N=1024 ###

## Data location ##
measurement_name = ['adwindata']

timestamp = ['20140405_185814', '20140405_190116', '20140405_190604', '20140405_191312']
cum_pts = 0

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

fig = a.default_fig(figsize=(6,4))
ax = a.default_ax(fig)
ax.plot(np.log10(a.sweep_pts), a.p0, '.b', lw=1)
#ax.plot((a.sweep_pts), a.p0, '.', lw=1)


#fit
x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

offset      = 0.5
amplitude   = 0.4
position    = 0
T2          = 80000
power       = 1.5

p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)
fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])
if show_fit:
    ax.plot(np.log10(np.linspace(10**(0.5),x[-1]+50e4,2001)), fit_result['fitfunc'](np.linspace(0,x[-1]+50e4,2001)), '-b', lw=1)


### N=256 ###

## Data location ##
measurement_name = ['adwindata']

timestamp = ['20140405_145055', '20140405_150942', '20140405_152120', '20140405_170830', '20140405_171629', '20140405_175233', '20140405_180505']
cum_pts = 0

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
#ax = a.plot_results_vs_sweepparam(ret='ax')
ax.plot(np.log10(a.sweep_pts), a.p0, '.r', lw=1)

#fit
x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

offset      = 0.5
amplitude   = 0.45
position    = 0
T2          = 45000
power       = 1.5

p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)
fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])
if show_fit:
    ax.plot(np.log10(np.linspace(10**(0.5),x[-1]+50e4,2001)), fit_result['fitfunc'](np.linspace(0,x[-1]+50e4,2001)), '-r', lw=1)




### N=64###

timestamp = ['20140405_135543', '20140405_142700', '20140405_143211', '20140405_143754']

cum_pts = 0

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
ax.plot(np.log10(a.sweep_pts), a.p0, '.g', lw=1)

#fit
x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

offset      = 0.5
amplitude   = 0.45
position    = 0
T2          = 19000
power       = 3

p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)
fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])
if show_fit:
    ax.plot(np.log10(np.linspace(10**(0.5),x[-1]+50e4,2001)), fit_result['fitfunc'](np.linspace(0,x[-1]+50e4,2001)), '-g', lw=1)









### N=8###

timestamp = ['20140405_133301']

cum_pts = 0

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
ax.plot(np.log10(a.sweep_pts), a.p0, '.y', lw=1)

#fit
x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

offset      = 0.5
amplitude   = 0.45
position    = 0
T2          = 4500
power       = 3


p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)
fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])
if show_fit:
    ax.plot(np.log10(np.linspace(10**(0.5),x[-1]+50e4,2001)), fit_result['fitfunc'](np.linspace(0,x[-1]+50e4,2001)), '-y', lw=1)













### N=4###

timestamp = ['20140405_130516', '20140405_131954']

cum_pts = 0

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
ax.plot(np.log10(a.sweep_pts), a.p0, '.m', lw=1)

#fit
x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

offset      = 0.5
amplitude   = 0.48
position    = 0
T2          = 2400
power       = 2


p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)
fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])
if show_fit:
    ax.plot(np.log10(np.linspace(10**(0.5),x[-1]+50e4,2001)), fit_result['fitfunc'](np.linspace(0,x[-1]+50e4,2001)), '-m', lw=1)











### N=1###

timestamp = ['20140405_123712']

cum_pts = 0

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
ax.plot(np.log10(a.sweep_pts), a.p0, '.c', lw=1)

#fit
x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

offset      = 0.5
amplitude   = 0.45
position    = 0
T2          = 1000
power       = 3

p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)
fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])
if show_fit:
    ax.plot(np.log10(np.linspace(10**(0.5),x[-1]+50e4,2001)), fit_result['fitfunc'](np.linspace(0,x[-1]+50e4,2001)), '-c', lw=1)










# ## plot data and fit as function of total time
# if plot_fit == True:
# plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)














## fit paramaters initial values (general exponential function) ##
# offset      = 0.5
# amplitude   = 0.45
# position    = 0
# T2          = 18000
# power       = 3

## other settings ##
# plot_fit    = True
# show_guess  = True
# fit_results = []





#p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)
#plot the initial guess
#if show_guess:
#    ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), '-', lw=1)

#fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])

## plot data and fit as function of total time
#if plot_fit == True:
#    plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)

#fit_results.append(fit_result)





print folder
plt.savefig(os.path.join(folder, 'combined_result.pdf'),
format='pdf')
plt.savefig(os.path.join(folder, 'combined_result.png'),
format='png')
