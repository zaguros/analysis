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


timestamp2048RepT = []
timestamp2048NoShutter = []
timestamp2048Block = []
tst2048comparison = [] ### for freshly taken data with 2048 pulses.


#### loop the get the timestmaps for 2048. this can also be done with a list of newer and older than timestamps. (if more measurements are to be used.)

# new_tsmp = '20150430_104800' ## newer than
# old_tsmp = '20150430_121207' ## older than
new_tsmp = '20150508_172620' ## newer than
old_tsmp = '20150508_173208' ## older than
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
    timestamp2048RepT.append(old_tsmp)


new_tsmp = '20150507_144811'
old_tsmp = '20150507_150042'
search_string = 'sil18_Single_Block_ShutterNO_XY8_sweep_tau_N_2048'
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
    timestamp2048NoShutter.append(old_tsmp)


new_tsmp = '20150507_150237'
old_tsmp = '20150507_151715'
search_string = 'sil18_Single_Block_ShutterYES_XY8_sweep_tau_N_2048'
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
    timestamp2048Block.append(old_tsmp)


new_tsmp = '20150508_172136'
old_tsmp = '20150508_172517'
search_string = '_Single_Block_ShutterNO_XY8_sweep_tau_N_2048'
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
    tst2048comparison.append(old_tsmp)
print tst2048comparison


'''
b: blue
g: green
r: red
c: cyan
m: magenta
y: yellow
k: black
w: white
'''

timestamplist = [timestamp2048RepT,timestamp2048NoShutter,timestamp2048Block,tst2048comparison]
labels = ['repeating T','Block, no shutter','Shutter + Block','rerun repeating T']
timestamplist = [timestamp2048RepT,timestamp2048Block]
labels = ['repeating T','Three Blocks']
color_list =  ['b','r','g','c','y','k','m']
# N_list = [128,256,512,1024,2048,4096]
# XY_list = [8]*6



## fit paramaters initial values (general exponential function) ##
offset      = 0.5
amplitude   = 0.4
position    = 0
T2          = 400
# T2          = 30
power       = 3.

## other settings ##
plot_fit    = True
show_guess  = False
save_data = False
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
            reps_per_datapoint = a.reps
        else:
            cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
            cum_p0 = np.concatenate((cum_p0, a.p0))
            cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

    sorting_order=cum_sweep_pts.argsort()
    cum_sweep_pts.sort()
    cum_p0=cum_p0[sorting_order]
    cum_u_p0=cum_u_p0[sorting_order]

    a.pts   = cum_pts
    a.sweep_pts = cum_sweep_pts
    a.p0    = cum_p0
    a.u_p0  = cum_u_p0
    # print a.p0
    # print a.u_p0
    if ii == 0:
        ax = a.plot_results_vs_sweepparam(ret='ax',labels=[labels[ii]])
        x_max = max(cum_sweep_pts)
        x_min = min(cum_sweep_pts)    
    else:
        cum_sweep_pts = np.array(cum_sweep_pts)
        cum_p0 = np.array(cum_p0)
        # cum_u_p0 = np.array(cum_u_p0)
        x_max = max(x_max,max(cum_sweep_pts))
        x_min = min(x_min,min(cum_sweep_pts))
        cum_u_p0_ar = cum_u_p0
        cum_u_p0_ar = np.reshape(cum_u_p0_ar, (np.shape(cum_u_p0_ar)[0], ))
        # print np.shape(cum_u_p0)
        plt.errorbar(cum_sweep_pts.flatten(),cum_p0.flatten(),yerr=cum_u_p0_ar,fmt='o',label=labels[ii],color=color_list[ii])

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    y_u = cum_u_p0
    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)

    #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), '-', lw=2)
    
    # if save_data:
    #     print timestamp[-1]
    #     print cum_pts
    #     print reps_per_datapoint
    #     savestr = 'FirstTS_' + timestamp[0].replace(' ','') +'_LastTS_'+timestamp[-1].replace(' ','')+ '_Electron_DD_N' + str(labels[ii]) + '_XY' + str(XY_list[ii]) + '_Pts' + str(cum_pts) + '_Reps' + str(reps_per_datapoint) + '.txt'    
    #     save_folder_str = 'D:/Dropbox/QEC LT/Decoupling memory/Electron_DD_Data/' + savestr
    #     np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)
    # print
    # print np.shape(x)
    # print np.shape(y)
    # print 
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2,4])

    ## plot data and fit as function of total time
    print 
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(0,np.amax(x)+10.,201), ax=ax,add_txt = False, plot_data=False,color=color_list[ii])

    fit_results.append(fit_result)

    cum_pts = 0
    cum_sweep_pts       = np.empty(0)
    cum_p0              = np.empty(0)
    cum_u_p0              = np.empty(0)

    cum_normalized_ssro = np.empty(0)

print folder

ax.set_ylim(0.4,1.02)
ax.set_xlim(0.,x_max)
# ax.set_xlim(140.,230.)
# ax.set_xscale('log')
plt.legend()


plt.savefig(os.path.join(folder, 'combined_result.pdf'),
format='pdf')
plt.savefig(os.path.join(folder, 'combined_result.png'),
format='png')
