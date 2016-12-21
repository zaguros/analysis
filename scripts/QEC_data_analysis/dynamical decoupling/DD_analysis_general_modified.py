'''Script to analyze the dynamical decoupling data
by THT''' 

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt

from scipy.stats.stats import pearsonr
reload(common)


### Inputs ###

## Data location ##
measurement_name = ['adwindata']
timestamp = []
temperature_list = []
absolute_time_list = []

new_tsmp = '20161101_131800' ## newer than
old_tsmp = '20161101_140000' ## older than

search_string = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY16sweep_tau_N_1024_part1'
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

    absolute_time = int(old_tsmp[-2:]) + int(old_tsmp[-4:-2])*60 + int(old_tsmp[-6:-4])*60**2 + int(old_tsmp[6:8])*(60**2)*24
    absolute_time_list.append(absolute_time)
    timestamp.append(old_tsmp)


timestamp = timestamp[::-1]
absolute_time_list = absolute_time_list[::-1]

print 'time stamp is    ' +str(timestamp)

offset      = 0.5
amplitude   = 0.45
position    = 0
T2          = 150
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
    temperature = (a.g.attrs['temp']-100)/0.385
    temperature_list.append(temperature)
    

    if kk == 0:
        cum_sweep_pts = a.sweep_pts
        cum_p0 = a.p0
        cum_u_p0 = a.u_p0
    else:
        cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
        cum_p0 = np.concatenate((cum_p0, a.p0))
        cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))


#print cum_p0
#print cum_u_p0
relative_time   = (np.array(absolute_time_list)-absolute_time_list[0])/60.


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
    plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)

fit_results.append(fit_result)

print folder
plt.savefig(os.path.join(folder, 'combined_result.pdf'),
format='pdf')
plt.savefig(os.path.join(folder, 'combined_result.png'),
format='png')



#plot fidelity vs time

fig = plt.figure(2,figsize=(18,5))
ax = fig.add_subplot(111)
ax.errorbar(relative_time,cum_p0,cum_u_p0)
ax.set_xlabel('minutes from ' + timestamp[-1])
ax.set_ylabel('fidelity')
# splt.ylim(48s,500)
# plt.ylim(np.max(relative_freq)-400,np.max(relative_freq)+50)
# plt.xlim(-300,-100)
plt.title('fidelity_vs_time')
plt.savefig(folder+'/fidelity_vs_time.png',format='png')

# plot temperature vs time
fig = plt.figure(3,figsize=(18,5))
plt.plot(relative_time,temperature_list)
plt.xlabel('minutes from ' + timestamp[-1])
plt.ylabel('Temperature')
plt.title('Stability measurement M1 - temperature')
plt.savefig(folder+'/temperature_vs_time.png',format='png')




# plot fidelity distribution
mean_fidelity=np.mean(cum_p0)
stdev_fidelity=np.std(cum_p0) 

plt.figure(4)
n, bins, patches = plt.hist(cum_p0,50,normed = 1)
bincenters = 0.5*(bins[1:]+bins[:-1])
# y = mlab.normpdf( bincenters, 1e6*(mean_fms0-1.746), 1e6*stdev_fms0)
# plt.plot(bincenters, y, 'r--', linewidth=1)
plt.xlabel('binned fidelity')
plt.title('Mean '+str(round(mean_fidelity,6))+' stdev '+str(round(stdev_fidelity,2)))
plt.savefig(folder+'/binned_fidelity.png',format='png')

#pearsonr(np.array(temperature_list),cum_p0)



#plot fidelity vs temperature

fig = plt.figure(5,figsize=(10,5))
ax = fig.add_subplot(111)
ax.scatter(temperature_list, cum_p0)
ax.set_xlabel('temperature')
ax.set_ylabel('fidelity')
plt.title('fidelity_vs_temperature')
plt.savefig(folder+'/fidelity_vs_temperature.png',format='png')
plt.show()