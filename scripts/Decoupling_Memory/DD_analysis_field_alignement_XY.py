'''Script to analyze the dynamical decoupling data for eigenstate
by MA''' 

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
timestamp1 = []
timestamp2 = []
timestamp3 = []
timestamp4 = []
timestamp5 = []
timestamp6 = []
timestamp7 = []
timestamp8 = []
timestamp9 = []
timestamp10 = []
timestamp11 = []
timestamp12 = []
timestamp13 = []
timestamp14 = []
timestamp15 = []
timestamp16 = []
timestamp17 = []
timestamp18 = []
timestamp19 = []
timestamp20 = []
timestamp21 = []

timestamp64 = []
timestamp128=[]
timestamp256=[]
timestamp512=[]
timestamp1024=[]
timestamp16_m7=[]
timestamp16_p4=[]
timestamp16_m2=[]
temperature_list = []
absolute_time_list = []

# ##############################################
# new_tsmp = '20161222_233800' ## newer than
# old_tsmp = '20161223_010000' ## older than
# search_string = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY8sweep_tau_N_64'
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


#     timestamp64.append(old_tsmp)

# timestamp64 = timestamp64[::-1]
# # # print timestamp64


#####################################################
new_tsmp = '20170213_234000' ## newer than
old_tsmp = '20170214_004000' ## older than

search_string = '_DecouplingSequence_111_1'
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


    timestamp1.append(old_tsmp)
    
timestamp1 = timestamp1[::-1]
print timestamp1


#########################################################
new_tsmp = '20170214_141600' ## newer than
old_tsmp = '20170214_144000' ## older than

search_string = '_DecouplingSequence_111_1'
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


    timestamp2.append(old_tsmp)
    
timestamp2 = timestamp2[::-1]

# ###########################################
new_tsmp = '20170214_144700' ## newer than
old_tsmp = '20170214_183500' ## older than

search_string = '_DecouplingSequence_111_1'
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


    timestamp3.append(old_tsmp)
    
timestamp3 = timestamp3[::-1]

# #################################

new_tsmp = '20170214_183800' ## newer than
old_tsmp = '20170214_190800' ## older than

search_string = '_DecouplingSequence_111_1'
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


    timestamp4.append(old_tsmp)
    
timestamp4 = timestamp4[::-1]

# #################################

new_tsmp = '20170214_190800' ## newer than
old_tsmp = '20170214_193300' ## older than

search_string = '_DecouplingSequence_111_1_'
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


    timestamp5.append(old_tsmp)
    
timestamp5 = timestamp5[::-1]

# #################################

new_tsmp = '20170214_193600' ## newer than
old_tsmp = '20170214_200000' ## older than

search_string = '_DecouplingSequence_111_1'
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


    timestamp6.append(old_tsmp)
    
timestamp6 = timestamp6[::-1]

# #################################

new_tsmp = '20170214_200000' ## newer than
old_tsmp = '20170214_203000' ## older than

search_string = '_DecouplingSequence_111_'
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


    timestamp7.append(old_tsmp)
    
timestamp7 = timestamp7[::-1]
#################################

# #################################

new_tsmp = '20170210_221000' ## newer than
old_tsmp = '20170210_224200' ## older than

search_string = '_DecouplingSequence_111'
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


    timestamp8.append(old_tsmp)
    
timestamp8 = timestamp8[::-1]

# #################################

new_tsmp = '20170210_224200' ## newer than
old_tsmp = '20170210_231400' ## older than

search_string = '_DecouplingSequence_111_'
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


    timestamp9.append(old_tsmp)
    
timestamp9 = timestamp9[::-1]

# #################################

new_tsmp = '20170210_231800' ## newer than
old_tsmp = '20170210_234700' ## older than

search_string = 'DecouplingSequence_111_1_'
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


    timestamp10.append(old_tsmp)
    
timestamp10 = timestamp10[::-1]

# #################################

new_tsmp = '20170210_234700' ## newer than
old_tsmp = '20170211_001900' ## older than

search_string = 'DecouplingSequence_111_1_'
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


    timestamp11.append(old_tsmp)
    
timestamp11 = timestamp11[::-1]

# #################################

new_tsmp = '20170211_181700' ## newer than
old_tsmp = '20170211_200000' ## older than

search_string = '_DecouplingSequence_111'
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


    timestamp12.append(old_tsmp)
    
timestamp12 = timestamp12[::-1]

# #################################

new_tsmp = '20170209_165900' ## newer than
old_tsmp = '20170209_173200' ## older than

search_string = '_DecouplingSequence_111_1_'
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


    timestamp13.append(old_tsmp)
    
timestamp13 = timestamp13[::-1]

# #################################

new_tsmp = '20170213_113000' ## newer than
old_tsmp = '20170213_124000' ## older than

search_string = '_DecouplingSequence_111_1_'
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


    timestamp14.append(old_tsmp)
    
timestamp14 = timestamp14[::-1]

# #################################

new_tsmp = '20170213_214000' ## newer than
old_tsmp = '20170213_223000' ## older than

search_string = '_DecouplingSequence_111'
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


    timestamp15.append(old_tsmp)
    
timestamp15 = timestamp15[::-1]

# #################################

new_tsmp = '20170204_035400' ## newer than
old_tsmp = '20170204_101100' ## older than

search_string = '_DecouplingSequence_111_1_sil18_S'
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


    timestamp16.append(old_tsmp)
    
timestamp16 = timestamp16[::-1]

# #################################

new_tsmp = '20170204_035400' ## newer than
old_tsmp = '20170204_100000' ## older than

search_string = '_DecouplingSequence_111_1'
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


    timestamp17.append(old_tsmp)
    
timestamp17 = timestamp17[::-1]


timestamp_list=[timestamp1,timestamp2,timestamp3,timestamp4,timestamp5,timestamp6,timestamp7]#,timestamp8,timestamp9,timestamp10,timestamp11,timestamp12,timestamp13,timestamp14]#,timestamp14,timestamp15,timestamp16,timestamp17]
labels = ['4.54','4.58','4.62','4.66','4.70','4.74','4.78','4.68','4.7','4.72','4.74','4.4-dZFS 346 kHz','older_data_4.665','newer_Data_4.665']#['64-Eigenstate-optimal','64-Eigenstate-overrotation','64-Superposition-optimal','256-Eigenstate-optimal','256-Eigenstate-over','256-superposition','-0.8','-0.6','-0.4','-0.2','0','+0.2','+0.4','0.6','+0.8','+1','1.2','1.4','1.6']
# timestamp_list=[timestamp13,timestamp14,timestamp15,timestamp7]
# labels = ['old_4.665-superposition','new_4.665_eigenstate','new_4.665_superposition','4.66_eigenstate']

N = [3072,3072,3072,3072,3072,3072,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,]
color_list =  ['b','r','g','c','y','k','m','b','r','g','c','y','k','m','b','r','g','c','y','k','m']
fmt_label = ['-o','-o','-o','-o','-o','-o','-o','--o','--o','--o','--o','--o','--o','--o','--o','--o']


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
cu                  =np.empty(0) 
cu_u                =np.empty(0) 
cum_u_p0            = np.empty(0)

cum_normalized_ssro = np.empty(0)
#for k in range(0,len(measurement_name)):

for ii, timestamp in enumerate(timestamp_list):
    print timestamp

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

    #print cum_sweep_pts
    #print cum_p0
    #print cum_u_p0
    #relative_time   = (np.array(absolute_time_list)-absolute_time_list[0])/60.


    a.pts   = cum_pts
    a.sweep_pts = cum_sweep_pts
    a.p0    = cum_p0
    a.u_p0  = 0*cum_u_p0


    # if ii == 0:
    #     ax = a.plot_results_vs_sweepparam(ret='ax',labels=[labels[ii]],figsize=(10,6))
    #     x_max = max(cum_sweep_pts)
    #     x_min = min(cum_sweep_pts)    
    # else:
    #     cum_sweep_pts = np.array(cum_sweep_pts)
    #     cum_p0 = np.array(cum_p0)
    #     # cum_u_p0 = np.array(cum_u_p0)
    #     x_max = max(x_max,max(cum_sweep_pts))
    #     x_min = min(x_min,min(cum_sweep_pts))
    #     cum_u_p0_ar = cum_u_p0
    #     cum_u_p0_ar = np.reshape(cum_u_p0_ar, (np.shape(cum_u_p0_ar)[0], ))
    #     # print np.shape(cum_u_p0)
    #     plt.errorbar(cum_sweep_pts.flatten(),cum_p0.flatten(),yerr=cum_u_p0_ar,fmt='o',label=labels[ii],color=color_list[ii])

   

    # x = a.sweep_pts.reshape(-1)[:]
    # y = a.p0.reshape(-1)[:]

    # p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)

    # # plt.plot(cum_sweep_pts,cum_p0,'o')

    # #plot the initial guess
    # # if show_guess:
    # #     ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), '-', lw=2)

    # fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])

    # ## plot data and fit as function of total time
     
    # if plot_fit == True:
    #     plot.plot_fit1d(fit_result, np.linspace(0,np.amax(x)+10.,201), ax=ax,add_txt = False, plot_data=False,color=color_list[ii])

    # fit_results.append(fit_result)
    # cum_pts = 0
    # cum_sweep_pts       = np.empty(0)
    # cum_p0              = np.empty(0)
    # cum_u_p0              = np.empty(0)



###################################################### plot against tau (in tau larmor)

    a.sweep_pts = 1000*cum_sweep_pts/(2*3072)

    if ii == 0:
        ax = a.plot_results_vs_sweepparam(ret='ax',name='fidelity_vs_tau',fmt='-o',labels=[labels[ii]],figsize=(40,8))
        x_max = max(1000*cum_sweep_pts/(2*2.3158*N[ii]**1))
        x_min = min(1000*cum_sweep_pts/(2*2.3158*N[ii]**1))    
    else:
        cum_sweep_pts = np.array(cum_sweep_pts)
        cum_p0 = np.array(cum_p0)
        # cum_u_p0 = np.array(cum_u_p0)
        x_max = max(x_max,max(1000*cum_sweep_pts/(2*2.3158*N[ii]**1.)))
        x_min = min(x_min,min(1000*cum_sweep_pts/(2*2.3158*N[ii]**1.)))
        cum_u_p0_ar = cum_u_p0
        cum_u_p0_ar = np.reshape(cum_u_p0_ar, (np.shape(cum_u_p0_ar)[0], ))
        # print np.shape(cum_u_p0)
        plt.errorbar(cum_sweep_pts.flatten()/(1e-3*2*3072),cum_p0.flatten(),yerr=0*cum_u_p0_ar,fmt=fmt_label[ii],label=labels[ii],color=color_list[ii])

   

    # x = a.sweep_pts.reshape(-1)[:]
    # y = a.p0.reshape(-1)[:]

    # p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)

  

    # # plot the initial guess
    # # if show_guess:
    # #     ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), '-', lw=2)

    # # fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])

    # ## plot data and fit as function of total time
     
    # # if plot_fit == True:
    # #     plot.plot_fit1d(fit_result, np.linspace(0,np.amax(x)+10.,201), ax=ax,add_txt = False, plot_data=False,color=color_list[ii])

    # # fit_results.append(fit_result)
    cum_pts = 0
    cum_sweep_pts       = np.empty(0)
    cum_p0              = np.empty(0)
    cum_u_p0              = np.empty(0)
#############################################################################################

print folder
ax.set_ylim(0.3,1.02)
# ax.set_xlim(0.,x_max)
ax.set_xlabel('Tau (us) ')
plt.axvline(x=2.314, ymin=0, ymax=1, linewidth=0.5)
plt.grid()

ax.set_xlim(2.18,2.45)
# ax.set_xscale('log')
plt.legend()
plt.savefig(os.path.join(folder, 'combined_result.pdf'),
format='pdf')
plt.savefig(os.path.join(folder, 'combined_result.png'),
format='png')



