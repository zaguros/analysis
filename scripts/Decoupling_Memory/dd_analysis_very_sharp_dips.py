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
timestamp64 = []
timestamp128=[]
timestamp256=[]
timestamp256_old=[]
timestamp256_2=[]
timestamp256_22=[]
timestamp256_3=[]
timestamp256_4=[]
timestamp256_5=[]
timestamp256_6=[]
timestamp256_7=[]
timestamp512=[]
timestamp1024=[]
temperature_list = []
absolute_time_list = []

##############################################
# new_tsmp = '20161212_233000' ## newer than
# old_tsmp = '20161213_233558' ## older than

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
# # print timestamp64


# #####################################################
# new_tsmp = '20161213_114000' ## newer than
# old_tsmp = '20161213_114558' ## older than

# search_string = '_DecouplingSequence_111_1_sil18_RepT'
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
    
# timestamp128 = timestamp128[::-1]
# # print timestamp128

#########################################################
# new_tsmp = '20161214_122000' ## newer than
# old_tsmp = '20161214_131000' ## older than

# search_string = '_DecouplingSequence_111_1_sil18_RepT'
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


#     timestamp256_old.append(old_tsmp)
    
# timestamp256_old = timestamp256_old[::-1]


#########################################################
new_tsmp = '20170209_214100' ## newer than
old_tsmp = '20170210_042800' ## older than


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


    timestamp256.append(old_tsmp)
    
timestamp256 = timestamp256[::-1]


#########################################################
new_tsmp = '20170207_192300' ## newer than
old_tsmp = '20170207_202800' ## older than

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


    timestamp256_2.append(old_tsmp)
    
timestamp256_2 = timestamp256_2[::-1]

#########################################################
new_tsmp = '20170207_182014'## newer than
old_tsmp = '20170207_192016' ## older than

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


    timestamp256_22.append(old_tsmp)
    
timestamp256_22 = timestamp256_22[::-1]



#########################################################
new_tsmp = '20161217_164000' ## newer than
old_tsmp = '20161217_164200' ## older than

search_string = '_DecouplingSequence_111_1_sil18_RepT'
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


    timestamp256_3.append(old_tsmp)
    
timestamp256_3 = timestamp256_3[::-1]
#########################################################
new_tsmp = '20161217_165900' ## newer than
old_tsmp = '20161217_166000' ## older than

search_string = '_DecouplingSequence_111_1_sil18_RepT'
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


    timestamp256_4.append(old_tsmp)
    
timestamp256_4 = timestamp256_4[::-1]



#########################################################
new_tsmp = '20161217_171000' ## newer than
old_tsmp = '20161217_173000' ## older than

search_string = '_DecouplingSequence_111_1_sil18_RepT'
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


    timestamp256_5.append(old_tsmp)
    
timestamp256_5 = timestamp256_5[::-1]


#########################################################
new_tsmp = '20161217_173700' ## newer than
old_tsmp = '20161217_174400' ## older than

search_string = '_DecouplingSequence_111_1_sil18_RepT'
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


    timestamp256_6.append(old_tsmp)
    
timestamp256_6 = timestamp256_6[::-1]

#########################################################
new_tsmp = '20161217_181800' ## newer than
old_tsmp = '20161217_183300' ## older than

search_string = '_DecouplingSequence_111_1_sil18_RepT'
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


    timestamp256_7.append(old_tsmp)
    
timestamp256_7 = timestamp256_7[::-1]
###########################################
# new_tsmp = '20161129_015400' ## newer than
# old_tsmp = '20161129_022058' ## older than

# search_string = '_DecouplingSequence_111_1_sil18_RepT'
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
    
# timestamp512 = timestamp512[::-1]

# #################################

# new_tsmp = '20161213_025900' ## newer than
# old_tsmp = '20161213_052058' ## older than

# search_string = '_DecouplingSequence_111_1_sil18_RepT'
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


#     timestamp1024.append(old_tsmp)
    
# timestamp1024 = timestamp1024[::-1]

#################################

timestamp_list=[timestamp256,timestamp256_2,timestamp256_22]#,timestamp256_3,timestamp256_4,timestamp256_5,timestamp256_6]
labels = ['3072','1024','512','64','96','128','160']
N = [3072,1024,512]
color_list =  ['b','r','g','c','y','k','m','b']


print 'time stamp is    ' +str(timestamp_list)

offset      = 0.5
amplitude   = 0.45
position    = 0
T2          = 150
power       = 1.5

## other settings ##
plot_fit    = False
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
    a.sweep_pts = 1.0e3*cum_sweep_pts/(2*N[ii])
    a.p0    = cum_p0
    a.u_p0  = 0*cum_u_p0


    if ii == 0:
        ax = a.plot_results_vs_sweepparam(ret='ax',fmt='-o',labels=[labels[ii]],figsize=(36,12))
        #x_max = max(cum_sweep_pts)
        #x_min = min(cum_sweep_pts)    
    else:
        cum_sweep_pts = np.array(cum_sweep_pts)
        cum_p0 = np.array(cum_p0)
        # cum_u_p0 = np.array(cum_u_p0)
        #x_max = max(x_max,max(cum_sweep_pts))
        #x_min = min(x_min,min(cum_sweep_pts))
        cum_u_p0_ar = cum_u_p0
        cum_u_p0_ar = np.reshape(cum_u_p0_ar, (np.shape(cum_u_p0_ar)[0], ))
        # print np.shape(cum_u_p0)
        plt.errorbar(cum_sweep_pts.flatten()/(2*N[ii]*1.0e-3),cum_p0.flatten(),yerr=0*cum_u_p0_ar,fmt='-o',label=labels[ii],color=color_list[ii])

   

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

    # a.sweep_pts = cum_sweep_pts/(2*N[ii]**1.)

    # if ii == 0:
    #     ax = a.plot_results_vs_sweepparam(ret='ax',name='fidelity_vs_tau',labels=[labels[ii]],figsize=(10,6))
    #     x_max = max(cum_sweep_pts/(2*N[ii]**1))
    #     x_min = min(cum_sweep_pts/(2*N[ii]**1))    
    # else:
    #     cum_sweep_pts = np.array(cum_sweep_pts)
    #     cum_p0 = np.array(cum_p0)
    #     # cum_u_p0 = np.array(cum_u_p0)
    #     x_max = max(x_max,max(cum_sweep_pts/(2*N[ii]**1.)))
    #     x_min = min(x_min,min(cum_sweep_pts/(2*N[ii]**1.)))
    #     cum_u_p0_ar = cum_u_p0
    #     cum_u_p0_ar = np.reshape(cum_u_p0_ar, (np.shape(cum_u_p0_ar)[0], ))
    #     # print np.shape(cum_u_p0)
    #     plt.errorbar(cum_sweep_pts.flatten()/(2*N[ii]**1.),cum_p0.flatten(),yerr=cum_u_p0_ar,fmt='o',label=labels[ii],color=color_list[ii])

   

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
    # cum_pts = 0
    # cum_sweep_pts       = np.empty(0)
    # cum_p0              = np.empty(0)
    # cum_u_p0              = np.empty(0)
#############################################################################################

# mean_switching=[0.96663071353427199, 0.94961382027607688, 0.88815650869974549, 0.82289096680456064, 0.74745290559569089, 0.7453359313305451, 0.69172356255372258]
# evolution_time=[ 7.110656,  18.962432,  30.81216,   42.663936,  54.515712,  66.367488,
#   78.217216]
# plt.plot(evolution_time, mean_switching, 'mo', label='switch_no_attenuator')
print folder
#ax.set_ylim(0.3,1.03)
#ax.set_xlim(.1,x_max)
ticks=[2.4739662673301623, 2.5550799154393475, 2.6361935635485336, 2.7173072116577188, 2.7984208597669049, 2.8795345078760901, 2.9606481559852762, 3.0417618040944614, 3.1228754522036475, 3.2039891003128327, 3.2851027484220188, 3.366216396531204, 3.4473300446403901, 3.5284436927495753, 3.6095573408587613, 3.6906709889679465, 3.7717846370771326, 3.8528982851863178, 3.9340119332955039, 4.0151255814046891, 4.0962392295138752, 4.1773528776230604, 4.2584665257322465, 4.3395801738414317, 4.4206938219506178, 4.501807470059803, 4.5829211181689891, 4.6640347662781743, 4.7451484143873603, 4.8262620624965455, 4.9073757106057316]
ticks2=[2.4813501439668433, 2.5618263648522546, 2.6423025857376654, 2.7227788066230767, 2.8032550275084875, 2.8837312483938988, 2.9642074692793101, 3.044683690164721, 3.1251599110501322, 3.2056361319355431, 3.2861123528209544, 3.3665885737063657, 3.4470647945917765, 3.5275410154771882, 3.6080172363625991, 3.6884934572480104, 3.7689696781334217, 3.8494458990188325, 3.9299221199042433, 4.0103983407896546, 4.0908745616750659, 4.1713507825604772, 4.2518270034458876, 4.3323032243312989, 4.4127794452167102, 4.4932556661021215, 4.5737318869875327, 4.654208107872944, 4.7346843287583544, 4.8151605496437657, 4.895636770529177]
ticks=[2.4,2.42,2.44,2.46,2.48,2.5,2.52,2.54,2.56,2.58,2.6]
ax.xaxis.set_ticks(ticks)
# for xc in ticks:
#     plt.axvline(x=xc,color='g')

# for xc in ticks2:
#     plt.axvline(x=xc,color='r')

ax.set_xlabel('Tau (us) ')
ax.set_xlim(2.4,2.6)
#ax.grid()
# ax.set_xlim(140.,230.)
#ax.set_xscale('log')
plt.legend(loc='lower left')
plt.savefig(os.path.join(folder, 'combined_result.pdf'),
format='pdf')
plt.savefig(os.path.join(folder, 'combined_result.png'),
format='png')



