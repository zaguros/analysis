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


new_tsmp = '20170216_160000' ## newer than
old_tsmp = '20170216_232000' ## older than

new_tsmp = '20140404_152200' ## newer than  Hans N=16
old_tsmp = '20140404_161400' ## older than

new_tsmp = '20140404_163520' ## newer than  Hans N=32
old_tsmp = '20140404_163600' ## older than

search_string = 'N32'
#search_string = 'DecouplingSequence_111_1'
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
T2          = 5000
power       = 1.5

## other settings ##
plot_fit    = True
show_guess  = False

fit_results = []

cum_pts = 0
cum_sweep_pts       = np.empty(0)
cum_p0              = np.empty(0)
cu                  =np.empty(0) 
cu_u                =np.empty(0) 
cum_u_p0            = np.empty(0)
cum_tau_list            = np.empty(0)

cum_normalized_ssro = np.empty(0)
#for k in range(0,len(measurement_name)):

for kk in range(len(timestamp)):
    folder = toolbox.data_from_time(timestamp[kk])

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    cum_pts += a.pts
    # temperature = (a.g.attrs['temp']-100)/0.385
    # temperature_list.append(temperature)


    if kk == 0:
        cum_sweep_pts = a.sweep_pts
        cum_p0 = a.p0
        cum_u_p0 = a.u_p0
        #cum_tau_list = a.tau_list
    #elif kk in [5,10,15,21,26,31]: 
    else:
        cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
        cum_p0 = np.concatenate((cum_p0, a.p0))
        cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))
        #cum_tau_list = np.concatenate((cum_tau_list, a.tau_list))

#print cum_sweep_pts
#print cum_p0
#print cum_u_p0
relative_time   = (np.array(absolute_time_list)-absolute_time_list[0])/60.




a.pts   = cum_pts 
a.sweep_pts = cum_sweep_pts #/(2*4.371)
a.p0    = cum_p0
a.u_p0  = cum_u_p0

print cum_u_p0
#a.tau_list = cum_tau_list

# a.p0 = np.array(cum_p0[0],cum_p0[5],cum_p0[11],cum_p0[17],cum_p0[23],cum_p0[29])
# a.sweep_pts = np.array(cum_sweep_pts[0],cum_sweep_pts[5],cum_sweep_pts[11],cum_sweep_pts[17],cum_sweep_pts[23],cum_sweep_pts[29])

#ax = a.plot_results_vs_sweepparam(ret='ax',fmt='-o',figsize=(20,8))

x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

index_list=[]
# To avoid dips in the analysis
for jj in range (1,2):
    index_list=[]
    for ii in range(1,len(x)-40):
        #if y[ii+2] > 0.65:
        for kk in range (ii+1, len(x)-1):
            if y[kk] > 1.03*y[ii] :
                index_list.append(ii)
    x=np.delete(x,[index_list])
    y=np.delete(y,[index_list])

fig = plt.figure(2,figsize=(10,5))
ax = fig.add_subplot(111)

ax.plot(x,y,'o')

#x=np.delete(x,[28,29,30,31,32,33,34,50,51,52,53,54])
#y=np.delete(y,[28,29,30,31,32,33,34,50,51,52,53,54])

p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)

print fitfunc_str

#plot the initial guess
if show_guess:
    ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), '-', lw=2)

fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])

## plot data and fit as function of total time
if plot_fit == True:
    plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)

fit_results.append(fit_result)


#ax.set_xlim(6.0,7.0)
# ax.set_ylim(0.2,1.03)
#ax.set_xscale('log')
ax.set_xlabel('Total evolution time (ms)')
#plt.axvline(x=52*2.3147, ymin=0, ymax=1, linewidth=0.5)
#ticks= range(40,80,2)
#ticks=np.linspace(2.6,4.6,10)
# ticks=[2.4739662673301623, 2.5550799154393475, 2.6361935635485336, 2.7173072116577188, 2.7984208597669049, 2.8795345078760901, 2.9606481559852762, 3.0417618040944614, 3.1228754522036475, 3.2039891003128327, 3.2851027484220188, 3.366216396531204, 3.4473300446403901, 3.5284436927495753, 3.6095573408587613, 3.6906709889679465, 3.7717846370771326, 3.8528982851863178, 3.9340119332955039, 4.0151255814046891, 4.0962392295138752, 4.1773528776230604, 4.2584665257322465, 4.3395801738414317, 4.4206938219506178, 4.501807470059803, 4.5829211181689891, 4.6640347662781743, 4.7451484143873603, 4.8262620624965455, 4.9073757106057316]
# ticks2=[2.4813501439668433, 2.5618263648522546, 2.6423025857376654, 2.7227788066230767, 2.8032550275084875, 2.8837312483938988, 2.9642074692793101, 3.044683690164721, 3.1251599110501322, 3.2056361319355431, 3.2861123528209544, 3.3665885737063657, 3.4470647945917765, 3.5275410154771882, 3.6080172363625991, 3.6884934572480104, 3.7689696781334217, 3.8494458990188325, 3.9299221199042433, 4.0103983407896546, 4.0908745616750659, 4.1713507825604772, 4.2518270034458876, 4.3323032243312989, 4.4127794452167102, 4.4932556661021215, 4.5737318869875327, 4.654208107872944, 4.7346843287583544, 4.8151605496437657, 4.895636770529177]
#ax.xaxis.set_ticks(ticks)
# for xc in ticks:
#     plt.axvline(x=xc,color='g')

# for xc in ticks2:
#     plt.axvline(x=xc,color='r')

ax.grid()

print folder
plt.savefig(os.path.join(folder, 'combined_result.pdf'),
format='pdf')
plt.savefig(os.path.join(folder, 'combined_result.png'),
format='png')



# #plot fidelity vs time
# # This is to pick one of the data points
# point_number=3
# for i in range (0,len(cum_p0)):
#     if cum_sweep_pts[i]==cum_sweep_pts[point_number]:
#         cu=np.concatenate((cu, cum_p0[i]))
#         cu_u=np.concatenate((cu_u, cum_u_p0[i]))


# #print 'cu is '+str(cu)        

# fig = plt.figure(2,figsize=(18,5))
# ax = fig.add_subplot(111)
# ax.errorbar(relative_time,cu,cu_u)
# ax.set_xlabel('minutes from ' + timestamp[-1])
# ax.set_ylabel('fidelity')
# # splt.ylim(48s,500)
# # plt.ylim(np.max(relative_freq)-400,np.max(relative_freq)+50)
# # plt.xlim(-300,-100)
# plt.title('fidelity_vs_time')
# plt.savefig(folder+'/fidelity_vs_time.png',format='png')

# # plot temperature vs time
# fig = plt.figure(3,figsize=(18,5))
# plt.plot(relative_time,temperature_list)
# plt.xlabel('minutes from ' + timestamp[-1])
# plt.ylabel('Temperature')
# plt.title('Stability measurement M1 - temperature')
# plt.savefig(folder+'/temperature_vs_time.png',format='png')




# # plot fidelity distribution
# mean_fidelity=np.mean(cu)
# stdev_fidelity=np.std(cu) 

# plt.figure(4)
# n, bins, patches = plt.hist(cu,50,normed = 1)
# bincenters = 0.5*(bins[1:]+bins[:-1])
# # y = mlab.normpdf( bincenters, 1e6*(mean_fms0-1.746), 1e6*stdev_fms0)
# # plt.plot(bincenters, y, 'r--', linewidth=1)
# plt.xlabel('binned fidelity')
# plt.title('Mean '+str(round(mean_fidelity,6))+' stdev '+str(round(stdev_fidelity,2)))
# plt.savefig(folder+'/binned_fidelity.png',format='png')

# # pearsonr(np.array(temperature_list),cum_p0)



# #plot fidelity vs temperature

# fig = plt.figure(5,figsize=(10,5))
# ax = fig.add_subplot(111)
# ax.scatter(temperature_list, cu)
# ax.set_xlabel('temperature')
# ax.set_ylabel('fidelity')
# plt.title('fidelity_vs_temperature')
# plt.savefig(folder+'/fidelity_vs_temperature.png',format='png')
# plt.show()



# # # To plot average only
# mean_f  = []


# for j in range (0,7):
#     cu1     = np.empty(0)
#     point_number=j
#     for i in range (0,len(cum_p0)):
#         if cum_sweep_pts[i]==cum_sweep_pts[point_number]:
#             cu1=np.concatenate((cu1, cum_p0[i]))
            

#     mean_fidelity=np.mean(cu1)
#     mean_f.append(mean_fidelity)
# print mean_f
# print cum_sweep_pts

