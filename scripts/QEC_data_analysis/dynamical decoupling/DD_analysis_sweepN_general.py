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

timestamp1 = []

timestamp2 = []

timestamp4 = []

timestamp8 = [] 

timestamp16 = [] 

timestamp32 = [] 

timestamp64 = []

timestamp128 = []

timestamp128shutter = []

timestamp256 = []

timestamp512 = []

timestamp1024 = []

timestamp2048 = []

timestamp4096 = []

timestamp8192 = []

timestamp16384 = []

timestamp32768 = []

timestamp65536old = []

timestamp65536 = []


#### loop the get the timestmaps for 2048. this can also be done with a list of newer and older than timestamps. (if more measurements are to be used.)

# new_tsmp = '20150603_222730' ## newer than
# old_tsmp = '20150604_030652' ## older than

# search_string = '_DecouplingSequence_111_1_sil18_RepT_ShutterYES_XY8sweep_tau_N_256'
# while toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False) != False:
#     old_tsmp, folder = toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False)
#     timestamp256.append(old_tsmp)

# new_tsmp = '20150604_205612' ## newer than
# old_tsmp = '20150604_222819' ## older than

# search_string = '_DecouplingSequence_111_1_sil18_RepT_ShutterYES_XY8sweep_tau_N_128'
# while toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False) != False:
#     old_tsmp, folder = toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False)
#     timestamp128.append(old_tsmp)

# new_tsmp = '20150604_155513' ## newer than
# old_tsmp = '20150604_155515' ## older than

# search_string = '_DecouplingSequence_111_1_sil18_RepT_ShutterYES_XY8sweep_tau_N_128'
# while toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False) != False:
#     old_tsmp, folder = toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False)
#     timestamp128shutter.append(old_tsmp)

# new_tsmp = '20150604_030651' ## newer than
# old_tsmp = '20150604_083412' ## older than


# search_string = '_DecouplingSequence_111_1_sil18_RepT_ShutterYES_XY8sweep_tau_N_512'
# while toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False) != False:
#     old_tsmp, folder = toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False)
#     timestamp512.append(old_tsmp)

new_tsmp = '20161025_181000' ## newer than
old_tsmp = '20161025_183000'

search_string = 'DecouplingSequence_111_1_sil18_RepT4sweep_N_on_tau_L5'
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
    timestamp1024.append(old_tsmp)

# new_tsmp = '20150606_045358' ## newer than
# old_tsmp = '20150606_112435'

# search_string = '_DecouplingSequence_111_1_sil18_RepT_ShutterYES_XY8sweep_tau_N_2048'
# while toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False) != False:
#     old_tsmp, folder = toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False)
#     timestamp2048.append(old_tsmp)
# search_string = '_DecouplingSequence_111_1_sil18sweep_tau_N_2048'
# while toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False) != False:
#     print 'I WAS HERE'
#     old_tsmp, folder = toolbox.latest_data(contains=search_string,
#                                         return_timestamp =True,
#                                         older_than=old_tsmp,
#                                         newer_than=new_tsmp,
#                                         raise_exc=False)
#     timestamp2048.append(old_tsmp)

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

timestamplist = [timestamp1,timestamp2,timestamp4,timestamp8,timestamp16,timestamp32,timestamp64,timestamp128,timestamp256,timestamp512,timestamp1024,timestamp2048]
labels = ['1','2','4','8','16','32','64','128','256','512','1024','2048','4096']
color_list =  ['b','r','g','c','y','k','m','b','r','g','c','y']
# timestamplist = [timestamp1,timestamp2,timestamp4,timestamp8,timestamp16,timestamp32,timestamp64,timestamp128,timestamp256,timestamp512,timestamp1024,timestamp2048]
N_list = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096]
# N_list = [1]
XY_list = [8]*len(timestamplist)

print timestamplist
#timestampnr = 0
# timestamplist = [timestamplist[timestampnr]]
# labels = [labels[timestampnr]]

#color_list = [color_list[timestampnr]]
## fit paramaters initial values (general exponential function) ##
offset      = 0.5
amplitude   = 0.4
position    = 0
T2          = 200
# T2          = 30
power       = 3.

## other settings ##
plot_fit    = True
show_guess  = False
save_data = True
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
            if ii==0:
                cum_p0= 1.-a.p0

            cum_u_p0 = a.u_p0
            reps_per_datapoint = a.reps
        else:
            cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
            cum_p0 = np.concatenate((cum_p0, a.p0))
            cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

    #sorting_order=cum_sweep_pts.argsort()
    #cum_sweep_pts.sort()
    #cum_p0=cum_p0[sorting_order]
    #cum_u_p0=cum_u_p0[sorting_order]
    print type(cum_sweep_pts)
    print cum_sweep_pts
    a.pts   = cum_pts
    #cum_sweep_pts = np.array(range(len(cum_sweep_pts))[::-1])
    a.sweep_pts = cum_sweep_pts
    a.p0    = cum_p0
    a.u_p0  = cum_u_p0

    if ii == 0:
        ax = a.plot_results_vs_sweepparam(ret='ax',labels=[labels[ii]])
        # x_max = max(cum_sweep_pts)
        # x_min = min(cum_sweep_pts)    
    else:
        cum_sweep_pts = np.array(cum_sweep_pts)
        cum_p0 = np.array(cum_p0)
        # cum_u_p0 = np.array(cum_u_p0)
        # x_max = max(x_max,max(cum_sweep_pts))
        # x_min = min(x_min,min(cum_sweep_pts))
        cum_u_p0_ar = cum_u_p0
        cum_u_p0_ar = np.reshape(cum_u_p0_ar, (np.shape(cum_u_p0_ar)[0], ))
        # print np.shape(cum_u_p0)
        plt.errorbar(cum_sweep_pts.flatten(),cum_p0.flatten(),yerr=cum_u_p0_ar,fmt='o',label=labels[ii],color=color_list[ii])

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    y_u = cum_u_p0
    T2 = 1.2*N_list[ii]**(2./3)
    print 'Guess T2', T2
    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)

    #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), '-', lw=2)
    
    if save_data:
        print timestamp[-1]
        print N_list[ii]
        print XY_list[ii]
        print cum_pts
        print reps_per_datapoint
        savestr = 'FirstTS_' + timestamp[0].replace(' ','') +'_LastTS_'+timestamp[-1].replace(' ','')+ '_Electron_DD_N' + str(labels[ii]) + '_XY' + str(XY_list[ii]) + '_Pts' + str(cum_pts) + '_Reps' + str(reps_per_datapoint) + '.txt'    
        save_folder_str = 'D:/Dropbox/QEC LT/Decoupling memory/Electron_DD_Data_NEW/' + savestr
        np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)
    # print
    # print np.shape(x)
    # print np.shape(y)
    # print 
    if True:
        print N_list[ii]
        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])

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
x_max = 600
ax.set_xlim(0.1,x_max)
# ax.set_xlim(140.,230.)
ax.set_xscale('log')
plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))


plt.savefig(os.path.join(folder, 'combined_result.pdf'),
format='pdf')
plt.savefig(os.path.join(folder, 'combined_result.png'),
format='png')
