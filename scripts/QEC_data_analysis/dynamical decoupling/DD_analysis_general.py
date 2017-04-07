'''Script to analyze the dynamical decoupling data
by THT''' 

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt




from analysis.lib.m2.ssro import mbi
import analysis.lib.QEC.nuclear_spin_characterisation as SC 
import matplotlib.cm as cm
import os



from analysis.lib.tools import analysis_magnet_tools as amt
from analysis.lib.m2.ssro import sequence
from analysis.lib.fitting import dark_esr_auto_analysis; 
import matplotlib.mlab as mlab
from tempfile import TemporaryFile
reload(common)


### Inputs ###

## Data location ##
measurement_name = ['adwindata']
timestamp_list = []
temperature_list = []
absolute_time_list = []

newer_than = '20161028_000500' ## newer than
older_than = '20161028_003800' ## older than

search = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY8sweep_tau_N'

iteration = 0
# while toolbox.latest_data(contains=contains, older_than=older_than, newer_than=newer_than,raise_exc = False)!=False:
#     print 'ITERATION '+str(iteration)
while toolbox.latest_data(contains=contains,
                                        return_timestamp_list =True,
                                        older_than=older_than,
                                        newer_than=newer_than,
                                        raise_exc=False) != False:

    timestamp, folder = toolbox.latest_data(contains=search,
                                        return_timestamp_list =True,
                                        older_than=older_than,
                                        newer_than=newer_than,
                                        raise_exc=False)

    absolute_time = int(timestamp[-2:]) + int(timestamp[-4:-2])*60 + int(timestamp[-6:-4])*60**2 + int(timestamp[6:8])*(60**2)*24
    absolute_time_list.append(absolute_time)


    # timestamp,folder = toolbox.latest_data(contains=contains, older_than=older_than, newer_than=newer_than,return_timestamp_list = True)
    # print 'm folder '+folder

    timestamp_list.append(timestamp_list)
    iteration = iteration+1
    


timestamp_list = timestamp_list[::-1]
print 'time stamp is    ' +str(timestamp_list)

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

for kk in range(len(timestamp_list)):
    folder = toolbox.data_from_time(timestamp_list[kk])

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    cum_pts += a.pts
    temperature = (a.g.attrs['temp']-100)/0.385
    temperature_list.append(temperature)
    absolute_time = int(timestamp_list[-2:]) + int(timestamp_list[-4:-2])*60 + int(timestamp_list[-6:-4])*60**2 + int(timestamp_list[6:8])*(60**2)*24
    absolute_time_list.append(absolute_time)


    if kk == 0:
        cum_sweep_pts = a.sweep_pts
        cum_p0 = a.p0
        cum_u_p0 = a.u_p0
    else:
        cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
        cum_p0 = np.concatenate((cum_p0, a.p0))
        cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))


print cum_p0
relative_time   = (absolute_time_list-abs_time[0])/60.


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


fig = plt.figure(1,figsize=(18,5))
ax = fig.add_subplot(111)
ax.errorbar(relative_time,cum_p0,cum_u_p0*1e3)
ax.set_xlabel('minutes from ' + timestamp_list[-1])
ax.set_ylabel('fidelity')
# splt.ylim(48s,500)
# plt.ylim(np.max(relative_freq)-400,np.max(relative_freq)+50)
# plt.xlim(-300,-100)
plt.title('fidelity_vs_time')
plt.savefig(folder+'/fidelity_vs_time.png',format='png')



