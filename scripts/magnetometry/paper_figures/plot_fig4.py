#######################################
###   PLOT Fig 4                    ###
#######################################

import matplotlib as mpl
#mpl.rcParams['text.usetex']=True
#mpl.rcParams['text.latex.unicode']=True
execfile('D:/measuring/analysis/scripts/setup_analysis.py')
from matplotlib import rc, cm
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
reload(magnetometry)
load_data=False

def analyze_saved_simulations (timestamp,error_bars=False):
    mgnt_exp = magnetometry.AdaptiveMagnetometry(N=14, tau0=20e-9)
    mgnt_exp.error_bars=error_bars
    mgnt_exp.load_analysis (timestamp=timestamp)
    mgnt_exp.calculate_scaling()
    return mgnt_exp


def add_scaling_plot(timestamp, ax, exp_data, label, marker_settings, color):
    #adds a scaling plot to axis 'ax', loading from analyzed data with a given 'timestamp'
    #exp_data=boolena, if 'True' then data is plotted with markers and errorbars are calculated, 
    #otherwise it is considered a simulation, and plotted like a line
    #label, string for legend
    data_file = analyze_saved_simulations (timestamp=timestamp, error_bars=exp_data)

    ax.plot (data_file.total_time, data_file.sensitivity, marker_settings+color, label=label)
    if exp_data: 
        ax.fill_between (data_file.total_time, data_file.sensitivity-data_file.err_sensitivity, data_file.sensitivity+data_file.err_sensitivity, color=color, alpha=0.1)
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.ylabel ('$V_{H}$ T ')
    #plt.show()
    return ax

def fig4a_plot():
    f = plt.figure(figsize=(8,6))
    ax = f.add_subplot(1,1,1)
    ax = add_scaling_plot (timestamp = '20141125_131323', exp_data=True, ax=ax, label = '', marker_settings='o', color='b')
    ax = add_scaling_plot (timestamp = '20141126_103521', exp_data=True, ax=ax, label = '', marker_settings='^', color='r')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel ('total ramsey time [$\mu$s]', fontsize=15)
    ax.set_ylabel ('sensitivity [$\mu$T$^2$*Hz$^{-1}$]', fontsize=15)

    plt.show()

fig4a_plot()



'''
#mgnt_G5F7_N14_sim_sweep_phase_nonadapt=asim.analyze_saved_simulations('20141118_161545')
#mgnt_G5F7_N14_sim_sweep_phase_adapt=asim.analyze_saved_simulations('20141118_163111')
#mgnt_G6F2_N14_sim=asim.analyze_saved_simulations('20141118_165516')
#mgnt_G2F2_N14_sim=asim.analyze_saved_simulations('20141118_172359')
#mgnt_G5F7_N14_sim=asim.analyze_saved_simulations('20141112_182801')
#mgnt_G5F2_N9_sim_sweep_phase_nonadapt=asim.analyze_saved_simulations('20141118_192523')
mgnt_G5F2_N11_sim=analyze_saved_simulations('20141118_195110')
#mgnt_G5F7_N14_data_1period=asim.analyze_saved_simulations('20141112_130924')
#mgnt_G5F7_N14_sim_sweep_phase_nonadapt_OH=asim.analyze_saved_simulations('20141118_161545')
#mgnt_G5F2_N11_sim_OH=asim.analyze_saved_simulations('20141118_195110')

mgnt_G5F2_N14_data=analyze_saved_simulations('20141120_115830')
#mgnt_G5F2_N14_data_OH=asim.analyze_saved_simulations('20141120_115830')
mgnt_G5F2_N14_data_sweep_phase_adapt=analyze_saved_simulations('20141125_131323')
mgnt_G5F7_N14_data_sweep_phase_adapt=analyze_saved_simulations('20141121_143256')

if load_data:
    mgnt_G5F7_N14_sim=asim.analyze_saved_simulations('20141112_182801')
    mgnt_G5F7_N14_sim_2=asim.analyze_saved_simulations('20141117_134837')
    mgnt_G5F7_N14_data=asim.analyze_saved_simulations('20141113_205258')
    mgnt_G5F7_N14_data_1period=asim.analyze_saved_simulations('20141112_130924')
    mgnt_G2F1_N14_data=asim.analyze_saved_simulations('20141117_155419')
    mgnt_G2F1_N14_data_3=asim.analyze_saved_simulations('20141114_080410')
    mgnt_G2F1_N14_sim_again=asim.analyze_saved_simulations('20141117_132520')
    mgnt_G2F1_N14_sim_again_2=asim.analyze_saved_simulations('20141117_133843')
    mgnt_G2F1_N14_sim_3=asim.analyze_saved_simulations('20141117_142604')
    
    
    mgnt_G2F1_N14_lowerfid=asim.analyze_saved_simulations('20141117_181925')
    mgnt_G5F7_N14_sim_sweep_phase_nonadapt=asim.analyze_saved_simulations('20141118_135425')
    mgnt_G5F7_N14_sim_sweep_phase_adapt=asim.analyze_saved_simulations('20141118_134405')    
    mgnt_G5F2_N14_data_OH=asim.analyze_saved_simulations('20141120_115830')

mgnt_G5F2_N9_sim_sweep_phase_nonadapt=asim.analyze_saved_simulations('20141118_192523')   
fig = plt.figure(figsize=(7,3))
#fig.clear()
plt.clf()

##################
###  G=5, F=7  ###
##################

#mgnt=mgnt_G5F7_N14_sim
#20141112_182801
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,color='DarkGreen',label='   G=5 F=7 Sim',linewidth=2)

#mgnt=mgnt_G5F7_N14_sim_2
#20141117_134837
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,color='k',label='   G=5 F=7 Sim2')
 
#mgnt=mgnt_G5F7_N14_data
#20141113_205258
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'o',color='DarkGreen',label='   G=5 F=7 Data 4 periods',markersize=10)

#mgnt=mgnt_G5F7_N14_data_1period
#20141113_205258
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'o',color='DarkGreen',label='   G=5 F=7 Data',markersize=10)

#mgnt=mgnt_G5F7_N14_sim_sweep_phase_nonadapt
#20141118_135425
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'x',color='Crimson',label='   G=5 F=7 Sim Berry non adapt',markersize=10,mew=2)

#mgnt=mgnt_G5F7_N14_sim_sweep_phase_adapt
#20141118_134405
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'-',color='RoyalBlue',label='   G=5 F=7 Sim Berry adapt',linewidth=2)


#mgnt=mgnt_G5F7_N14_data_sweep_phase_adapt
#20141121_143256
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'+',color='RoyalBlue',label='   G=5 F=7 Data Berry adapt',markersize=10,mew=2)



##################
###  G=2, F=1  ###
##################

#mgnt=mgnt_G2F1_N14_sim
#20141114_112636
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'-',color='r',label='   G=2 F=1 Sim',markersize=10)

#mgnt=mgnt_G2F1_N14_data
#20141117_155419
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'*',color='Crimson',label='   G=2 F=1 Data',markersize=10)

#mgnt=mgnt_G2F1_N14_lowerfid
#20141117_181925
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'-',color='RoyalBlue',label='   G=2 F=1 Sim reduced Fid',linewidth=2)

#mgnt=mgnt_G2F1_N14_data_3
#20141114_080410
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'s',color='Crimson',label='   G=2 F=1 Data v3',markersize=10)


#mgnt=mgnt_G2F1_N14_sim_again
#20141117_132520
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'-',color='r',label='   G=2 F=1 Sim_v1',markersize=10)

#mgnt=mgnt_G2F1_N14_sim_again_2
#20141117_133843
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'-',color='r',label='   G=2 F=1 Sim_v2',markersize=10)
#
#mgnt=mgnt_G2F1_N14_sim_3
#20141117_142604
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'-',color='r',label='   G=2 F=1 Sim_v3',markersize=10)

##################
###  G=6, F=2  ###
##################

#mgnt=mgnt_G6F2_N14_sim
#
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'-',color='DarkGreen',label='   G=5 F=2 Sim',markersize=10)

#mgnt=mgnt_G2F2_N14_sim
#
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'.',color='DarkGreen',label='   G=2 F=2 Sim',markersize=10)

##################
###  G=5, F=2  ###
##################

mgnt=mgnt_G5F2_N9_sim_sweep_phase_nonadapt
#'20141118_192523'
plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'-',color='Crimson',label='   G=5 F=2 Sim Berry non adapt',linewidth=2)

mgnt=mgnt_G5F2_N14_data_sweep_phase_adapt
#'20141121_142425'
plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'x',color='Crimson',label='   G=5 F=2 Data Berry',markersize=10,mew=2)

mgnt=mgnt_G5F2_N11_sim
#'20141118_195110'
plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'-',color='DarkGreen',label='   G=5 F=2 Sim',linewidth=2)

mgnt=mgnt_G5F2_N14_data
#'20141120_115830'
plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'o',color='DarkGreen',label='   G=5 F=2 Data',markersize=10)



### WITH OVERHEAD###

#mgnt=mgnt_G5F7_N14_sim_sweep_phase_nonadapt_OH
#'20141118_161545'
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'x',color='Crimson',label='   G=5 F=7 Sim Berry non adapt',markersize=10,mew=2)

#mgnt=mgnt_G5F2_N11_sim_OH
#'20141118_195110'
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'-',color='DarkGreen',label='   G=5 F=2 Sim',markersize=10)

#mgnt=mgnt_G5F2_N14_data_OH
#'20141120_115830'
#plt.loglog(mgnt.total_time/mgnt.t0,mgnt.scaling_variance*mgnt.total_time/mgnt.t0,'o',color='DarkGreen',label='   G=5 F=2 Data',markersize=10)

plt.legend(loc=3,prop={'size':12})
plt.tick_params(axis='x', labelsize=14)
plt.tick_params(axis='y', labelsize=14)

plt.xlabel (r'$\frac{T}{\tau_{0}}$', fontsize=24)
plt.ylabel (r'$\frac{V_{H}T}{\tau_{0}}$',fontsize=24)\

#plt.xlabel (r'$\frac{T+T_{overhead}}{\tau_{0}}$', fontsize=24)
#plt.ylabel (r'$\frac{V_{H}(T+T_{overhead})}{\tau_{0}}$',fontsize=24)
plt.xlim([10**0,10**5])
plt.ylim([10**-2,10**1.5])
#plt.ylabel(r'$\frac{1}{2}$')
plt.show()

#fig.savefig(r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\analyzed data\20141125_124500_Berry_and_Cappellaro_G5F2_sim_and_data.pdf', bbox_inches='tight')
plt.close()
'''