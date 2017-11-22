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
timestamp22 = []
timestamp23 = []
timestamp24 = []
timestamp25 = []
timestamp26 = []

timestamp30 = []
timestamp31 = []
timestamp32 = []
timestamp33 = []
timestamp34 = []
timestamp35 = []
timestamp36 = []
timestamp37 = []
timestamp38 = []
timestamp39 = []
timestamp40 = []
timestamp41 = []
timestamp42 = []

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
temp_stamp1 = []
temp_stamp2 = []
temp_stamp3 = []
temp_stamp11 = []
temp_stamp12 = []
x_list = []
y_list = []
T_list = []
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


##################################################### Data for x 32 to 8192 #########
new_tsmp = '20170315_235400' ## newer than
old_tsmp = '20170316_001100' ## older than

search_string = 'DecouplingSequence_111_1_sil18_'
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


    temp_stamp1.append(old_tsmp)
temp_stamp1 =temp_stamp1[::-1]

new_tsmp = '20170316_105500' ## newer than
old_tsmp = '20170316_105700' ## older than

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

    temp_stamp2.append(old_tsmp)
temp_stamp2 =temp_stamp2[::-1]
 
new_tsmp = '20170316_140500' ## newer than
old_tsmp = '20170316_140600' ## older than

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

    temp_stamp3.append(old_tsmp)
temp_stamp3 =temp_stamp3[::-1]

timestamp1 = temp_stamp1+temp_stamp2+temp_stamp3

temp_stamp1=[]
temp_stamp2=[]
temp_stamp3=[]



#########################################################  Data for y 32 to 8192
new_tsmp = '20170316_120300' ## newer than
old_tsmp = '20170316_122500' ## older than

search_string = 'DecouplingSequence_111_1_sil18_'
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


    temp_stamp1.append(old_tsmp)
temp_stamp1 =temp_stamp1[::-1]

new_tsmp = '20170316_111900' ## newer than
old_tsmp = '20170316_112000' ## older than

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

    temp_stamp2.append(old_tsmp)
temp_stamp2 =temp_stamp2[::-1]

new_tsmp = '20170316_162600' ## newer than
old_tsmp = '20170316_163000' ## older than

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

    temp_stamp3.append(old_tsmp)
temp_stamp3 =temp_stamp3[::-1]

timestamp2 = temp_stamp1+temp_stamp2+temp_stamp3

temp_stamp1=[]
temp_stamp2=[]
temp_stamp3=[]

# ###########################################  Data for -x 32 to 8192
new_tsmp = '20170316_004200' ## newer than
old_tsmp = '20170316_010100' ## older than

search_string = 'DecouplingSequence_111_1_sil18_'
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


    temp_stamp1.append(old_tsmp)
temp_stamp1 =temp_stamp1[::-1]

new_tsmp = '20170316_111200' ## newer than
old_tsmp = '20170316_111300' ## older than

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

    temp_stamp2.append(old_tsmp)
temp_stamp2 =temp_stamp2[::-1]


new_tsmp = '20170316_164000' ## newer than
old_tsmp = '20170316_165000' ## older than

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

    temp_stamp3.append(old_tsmp)
temp_stamp3 =temp_stamp3[::-1]

timestamp3 = temp_stamp1+temp_stamp2+temp_stamp3

temp_stamp1=[]
temp_stamp2=[]
temp_stamp3=[]
############
# ###########################################  Data for -y 32 to 8192
new_tsmp = '20170316_183100' ## newer than
old_tsmp = '20170316_184500' ## older than

search_string = 'DecouplingSequence_111_1_sil18_'
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


    temp_stamp1.append(old_tsmp)
temp_stamp1 =temp_stamp1[::-1]

new_tsmp = '20170316_113500' ## newer than
old_tsmp = '20170316_113700' ## older than

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

    temp_stamp2.append(old_tsmp)
temp_stamp2 =temp_stamp2[::-1]


new_tsmp = '20170316_150700' ## newer than
old_tsmp = '20170316_150800' ## older than

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

    temp_stamp3.append(old_tsmp)
temp_stamp3 =temp_stamp3[::-1]

timestamp4 = temp_stamp1+temp_stamp2+temp_stamp3

temp_stamp1=[]
temp_stamp2=[]
temp_stamp3=[]





# ###########################################  Data for eigenstate 0 (32 to 8192)
new_tsmp = '20170316_011000' ## newer than
old_tsmp = '20170316_012800' ## older than

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


    temp_stamp1.append(old_tsmp)
temp_stamp1 =temp_stamp1[::-1]

new_tsmp = '20170316_114800' ## newer than
old_tsmp = '20170316_114900' ## older than

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

    temp_stamp2.append(old_tsmp)
temp_stamp2 =temp_stamp2[::-1]


new_tsmp = '20170316_165000' ## newer than
old_tsmp = '20170316_165300' ## older than

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

    temp_stamp3.append(old_tsmp)
temp_stamp3 =temp_stamp3[::-1]

timestamp5 = temp_stamp1+temp_stamp2+temp_stamp3

temp_stamp1=[]
temp_stamp2=[]
temp_stamp3=[]

# ###########################################  Data for eigenstate 1 (32 to 8192)
new_tsmp = '20170316_204000' ## newer than
old_tsmp = '20170316_222000' ## older than

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


    timestamp6.append(old_tsmp)
timestamp6 =timestamp6[::-1]



temp_stamp1=[]
temp_stamp2=[]
temp_stamp3=[]




'''
# Data for N= 3072 #######################################################################

# ################################# here data were flipped

new_tsmp = '20170220_131000' ## newer than
old_tsmp = '20170220_143100' ## older than

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

    temp_stamp1.append(old_tsmp)
temp_stamp1 =temp_stamp1[::-1]
#timestamp4.append(temp_stamp)

new_tsmp = '20170221_020200' ## newer than
old_tsmp = '20170221_024500' ## older than

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

    temp_stamp2.append(old_tsmp)
temp_stamp2 =temp_stamp2[::-1]
#timestamp4.append(temp_stamp)

new_tsmp = '20170220_143100' ## newer than
old_tsmp = '20170220_180000' ## older than


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

    temp_stamp3.append(old_tsmp)
temp_stamp3 =temp_stamp3[::-1]
timestamp4=temp_stamp1+temp_stamp2+temp_stamp3
    
#timestamp4 = timestamp4[::-1]

# ################################# ############### ################


# #################################

new_tsmp = '20170221_202500' ## newer than
old_tsmp = '20170222_150900' ## older than

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

'''
########################################################

# timestamp_list=[timestamp1,timestamp2,timestamp3,timestamp4,timestamp5,timestamp6,timestamp7,timestamp8,timestamp9,timestamp10,timestamp11,timestamp12,timestamp13,timestamp14]#,timestamp14,timestamp15,timestamp16,timestamp17]
# labels = ['4.54','4.56','4.58','4.6','4.62','4.64','4.66','4.68','4.7','4.72','4.74','4.4-dZFS 346 kHz','older_data_4.665','newer_Data_4.665']#['64-Eigenstate-optimal','64-Eigenstate-overrotation','64-Superposition-optimal','256-Eigenstate-optimal','256-Eigenstate-over','256-superposition','-0.8','-0.6','-0.4','-0.2','0','+0.2','+0.4','0.6','+0.8','+1','1.2','1.4','1.6']


timestamp_list=[timestamp1, timestamp2, timestamp3,timestamp4,timestamp5,timestamp6]


labels = ['|-i>','|+>','|+i>','|->','|o>','|1>']
labels_N = ['32','64','128','512','1024','2048','3072','4096','6144','10240','20480']
color_list =  ['black','darkred','r','orangered','gold','darkgreen','g','b','darkblue','indigo','violet','c','k','m','coral','olive']#,'deeppink','b','r','g','c','k','m','gold','b','r','g','c','k','m']
#color_list= [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
fmt_label = ['o','s','^','o','s','^','o','s','^','o','s','^','o','s','^','^','^','^','^','s','s','s','s','s']


print 'time stamp is    ' +str(timestamp4)

offset      = 0.5
amplitude   = 0.45
position    = 0
T2          = 1
power       = 1.5

## other settings ##
plot_fit    = False
show_guess  = False

fit_results = []

cum_pts = 0
cum_sweep_pts       =np.empty(0)
cum_p0              =np.empty(0)
cu                  =np.empty(0) 
cu_u                =np.empty(0) 
cum_u_p0            =np.empty(0)

cum_normalized_ssro = np.empty(0)
e_list=[]
# cum1 = np.empty(0)
# cum2 = np.empty(0)
# cum3 = np.empty(0)
# cum33 = np.empty(0)
# cum4= np.empty(0)
#for k in range(0,len(measurement_name)):


for ii, timestamp in enumerate(timestamp_list):
    #print timestamp

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

    
    a.pts   = cum_pts
    a.sweep_pts = cum_sweep_pts
    a.p0    = cum_p0
    a.u_p0  = cum_u_p0
       
    x_list.append(list(cum_sweep_pts))
    y_list.append(list(cum_p0))
    e_list.append(list(cum_u_p0))

   



################################################## 
    cum_pts = 0
    cum_sweep_pts       = np.empty(0)
    cum_p0              = np.empty(0)
    cum_u_p0              = np.empty(0)




 
    cum_pts = 0
    cum_sweep_pts       = np.empty(0)
    cum_p0              = np.empty(0)
    cum_u_p0            = np.empty(0)
    cum_s=np.empty(0)
    cum1=np.empty(0)
    cum2=np.empty(0)
    cum3=np.empty(0)
    cum4=np.empty(0)
    cum5=np.empty(0)
    cum6=np.empty(0)
    cum7=np.empty(0)

    
######################################### plotting data ################
fig = plt.figure(2,figsize=(5,3))
ax1 = fig.add_subplot(111)
ax1.set_color_cycle(color_list)



for jj in range(0,len(x_list)):

    
    ##############  Fitting the data #########
    x= np.array(x_list[jj]) #/(2*N[jj])
    y= np.array(y_list[jj])
    e= np.array(e_list[jj])

    x = x.reshape(-1)[:]
    y = y.reshape(-1)[:]
    e = e.reshape(-1)[:]
    
    if jj==0:
        y_avg=y
        e_avg=e**2
    else:
        y_avg=y_avg+y
        e_avg=e_avg**2+e**2


    #########################################################
    T2= 20

    ax1.errorbar(x.flatten(),y.flatten(),yerr=e[jj],fmt=fmt_label[jj],label=labels[jj],color=color_list[jj],markersize=3)

    
    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])

    # plot data and fit as function of total time
    plot_fit= False

    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(0,np.amax(x)+400.,1001), ax=ax1,add_txt = False, plot_data=False,color=color_list[jj],lw=1)

    fit_results.append(fit_result)
    T_list.append(fit_result['params_dict']['T'])
y_avg=y_avg/6 
e_avg=(e_avg**0.5)/6 

  
##################################
ax1.set_xlabel('Total free evolution time (ms)')
ax1.set_ylabel('State fidelity')
ax1.set_xscale('log')
#ax1.set_xlim(0,80)
#ax1.grid()
ax1.legend(bbox_to_anchor=(1.25,0.5),loc='center right',numpoints=1,fontsize=8,frameon=False)

folder='C:\Users\TUD277931\Dropbox\TaminiauLab\Projects\Coherence in multi-qubit systems\Paper\Figures'
plt.savefig(os.path.join(folder,'Process_fidelity_separate_inputs.pdf'),
format='pdf',bbox_inches='tight',pad_inches=0.2,transparent=True)

# folder='C:\Users\TUD277931\Dropbox\TaminiauLab\Projects\Coherence in multi-qubit systems\Paper\Figures'
# plt.savefig(os.path.join(folder,'Process_fidelity_detailed.pdf'),
# format='pdf',bbox_inches='tight',pad_inches=0.2,transparent=True)

# plt.savefig(os.path.join(folder, 'selected_max_vs_tot_time_log.pdf'),
# format='pdf')
# plt.savefig(os.path.join(folder, 'selected_max_vs_tot_time_log.png'),
# format='png')

# Plot the average




fig = plt.figure(3,figsize=(4,2.5))
ax2 = fig.add_subplot(111)
ax2.set_color_cycle(color_list)
for ii in range (0,len(x)):
    ax2.errorbar(x[ii],y_avg[ii],yerr=e_avg[ii], fmt=fmt_label[ii],color=color_list[ii],label=labels_N[ii],markersize=3) 


## fitting the average

y=y_avg
y = y.reshape(-1)[:]
print (x)
p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)

fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])

# plot data and fit as function of total time
plot_fit= True

if plot_fit == True:
    plot.plot_fit1d(fit_result, np.linspace(0,np.amax(x)+400.,1001), ax=ax2,add_txt = False, plot_data=False,color='green',lw=1)

fit_results.append(fit_result)
T_list.append(fit_result['params_dict']['T'])

####################################################




ax2.set_xlabel('Total evolution time (ms)')
ax2.set_ylabel('Average state fidelity')
#ax1.set_xlabel('Tau (us)')

ax2.set_xscale('log')
plt.axhline(y=0.6667, xmin= 0.55, xmax=1.e3, linewidth=0.5,color='gray')
#ax2.set_xlim(0.5e0,1.e5)
ax2.set_ylim(0.45,1.08)

# get handles
#handles, labels = ax2.get_legend_handles_labels()
# remove the errorbars
#handles = [h[0] for h in handles]
# use them in the legend
#ax2.legend(handles, labels, loc='upper left',numpoints=1)
#ax2.grid()
#ax2.legend(bbox_to_anchor=(1.25,0.5),loc='center right',ncol=1,fontsize=9)
ax2.legend(loc='bottom left',ncol=2,numpoints=1,fontsize=8,frameon=False)

folder='C:\Users\TUD277931\Dropbox\TaminiauLab\Projects\Coherence in multi-qubit systems\Paper\Figures'
plt.savefig(os.path.join(folder,'Process_fidelity.pdf'),
format='pdf',bbox_inches='tight',pad_inches=0.2,transparent=True)



        






