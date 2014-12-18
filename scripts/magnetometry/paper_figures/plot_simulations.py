#######################################
###   PLOT Fig 4                    ###
#######################################

import matplotlib
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

    ax.plot (data_file.total_time, data_file.sensitivity, marker_settings,color=color, label=label)
    if exp_data: 
        ax.fill_between (data_file.total_time, data_file.sensitivity-data_file.err_sensitivity, data_file.sensitivity+data_file.err_sensitivity, color=color, alpha=0.1)
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.ylabel ('$V_{H}$ T ')
    #plt.show()
    return ax

def compare_2plots(timestamp1, timestamp2, title, label1 = 'adapt', label2='non-adapt'):
    f = plt.figure(figsize=(8,6))
    ax = f.add_subplot(1,1,1)
    ax = add_scaling_plot (timestamp = timestamp1, exp_data=True, ax=ax, label = label1, marker_settings='o', color='b')
    ax = add_scaling_plot (timestamp = timestamp2, exp_data=True, ax=ax, label = label2, marker_settings='^', color='r')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel ('total ramsey time [$\mu$s]', fontsize=15)
    ax.set_ylabel ('$V_{H}T$', fontsize=15)
    plt.title (title, fontsize=20)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(15)
    plt.legend()
    plt.show()

def compare_multiple_plots(timestamps, labels, title):
    #MAX 10 plots!!!! Then no more markers!
    f = plt.figure(figsize=(8,6))
    ax = f.add_subplot(1,1,1)
    markers = ['o', '+', '^', 'v', '>', '<', 's', '*', 'D', '|']
    ccc = np.linspace(0,1,len(timestamps))
    for i,k in enumerate(timestamps):
        ax = add_scaling_plot (timestamp = k, exp_data=True, ax=ax, label = labels[i], marker_settings=markers[i], color=cm.Set1(ccc[i]))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel ('total ramsey time [$\mu$s]', fontsize=15)
    ax.set_ylabel ('$V_{H}T$', fontsize=15)
    plt.title (title, fontsize=20)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(15)
    plt.legend()
    plt.show()


def compare_G3_variableF_imperfRO():
    #compares adaptive (Cappellaro only) to non-adaptive (phase_update=True) for G=3
    compare_2plots (timestamp1= '20141215_152517', timestamp2='20141215_152604', title = 'G=3, F=0, RO_fid = 0.87')
    compare_2plots (timestamp1= '20141215_152820', timestamp2='20141215_152939', title = 'G=3, F=1, RO_fid = 0.87')
    compare_2plots (timestamp1= '20141215_153251', timestamp2='20141215_153454', title = 'G=3, F=2, RO_fid = 0.87')
    compare_2plots (timestamp1= '20141215_153913', timestamp2='20141215_154152', title = 'G=3, F=3, RO_fid = 0.87')
    compare_2plots (timestamp1= '20141215_154728', timestamp2='20141215_155054', title = 'G=3, F=4, RO_fid = 0.87')
    #compare_2plots (timestamp1= '20141215_135229', timestamp2='20141215_135352', title = 'G=3, F=5, RO_fid = 0.87')
    #compare_2plots (timestamp1= '20141215_135556', timestamp2='20141215_135723', title = 'G=3, F=6, RO_fid = 0.87')


def compare_G3_variableF_imperfRO_recalcPhase():
    #compares adaptive (calculating optimal phase for each msmnt) and normal Cappellaro
    compare_2plots (timestamp1= '20141215_152517', timestamp2='20141215_160631', title = 'G=3, F=0, RO_fid = 0.87', label1='Cappellaro', label2='always optimal')
    compare_2plots (timestamp1= '20141215_152820', timestamp2='20141215_160844', title = 'G=3, F=1, RO_fid = 0.87', label1='Cappellaro', label2='always optimal')
    compare_2plots (timestamp1= '20141215_153251', timestamp2='20141215_161225', title = 'G=3, F=2, RO_fid = 0.87', label1='Cappellaro', label2='always optimal')
    compare_2plots (timestamp1= '20141215_153913', timestamp2='20141215_161702', title = 'G=3, F=3, RO_fid = 0.87', label1='Cappellaro', label2='always optimal')
    compare_2plots (timestamp1= '20141215_154728', timestamp2='20141215_162257', title = 'G=3, F=4, RO_fid = 0.87', label1='Cappellaro', label2='always optimal')

def compare_G3_variableF_imperfRO_optCapp_vs_nonadaptive():
    #compares adaptive (Cappellaro optimal phase for each msmnt) to non-adaptive (phase_update=True) for G=3
    compare_2plots (timestamp1= '20141215_160631', timestamp2='20141215_152604', title = 'G=3, F=0, RO_fid = 0.87', label1='opt Cappellaro', label2='non-adaptive')
    compare_2plots (timestamp1= '20141215_160844', timestamp2='20141215_152939', title = 'G=3, F=1, RO_fid = 0.87', label1='opt Cappellaro', label2='non-adaptive')
    compare_2plots (timestamp1= '20141215_161225', timestamp2='20141215_153454', title = 'G=3, F=2, RO_fid = 0.87', label1='opt Cappellaro', label2='non-adaptive')
    compare_2plots (timestamp1= '20141215_161702', timestamp2='20141215_154152', title = 'G=3, F=3, RO_fid = 0.87', label1='opt Cappellaro', label2='non-adaptive')
    compare_2plots (timestamp1= '20141215_162257', timestamp2='20141215_155054', title = 'G=3, F=4, RO_fid = 0.87', label1='opt Cappellaro', label2='non-adaptive')

def check_recalcPhase_N1sweep():
    #compares adaptive (Cappellaro optimal phase for each msmnt) to the same, but sweeping phase of the first msmnts (N=1) according to m*pi/G
    compare_2plots (timestamp1= '20141215_172959', timestamp2='20141215_160631', title = 'G=3, F=0, RO_fid = 0.87', label1='N1_sweep', label2='always optimal')
    compare_2plots (timestamp1= '20141215_173214', timestamp2='20141215_160844', title = 'G=3, F=1, RO_fid = 0.87', label1='N1_sweep', label2='always optimal')
    compare_2plots (timestamp1= '20141215_173528', timestamp2='20141215_161225', title = 'G=3, F=2, RO_fid = 0.87', label1='N1_sweep', label2='always optimal')
    compare_2plots (timestamp1= '20141215_173929', timestamp2='20141215_161702', title = 'G=3, F=3, RO_fid = 0.87', label1='N1_sweep', label2='always optimal')
    compare_2plots (timestamp1= '20141215_174445', timestamp2='20141215_162257', title = 'G=3, F=4, RO_fid = 0.87', label1='N1_sweep', label2='always optimal')


def compare_G5_variableF_imperfRO():
    compare_2plots (timestamp1= '20141215_134142', timestamp2='20141215_134205', title = 'G=5, F=0, RO_fid = 0.87')
    compare_2plots (timestamp1= '20141215_134300', timestamp2='20141215_134334', title = 'G=5, F=1, RO_fid = 0.87')
    compare_2plots (timestamp1= '20141215_134443', timestamp2='20141215_134529', title = 'G=5, F=2, RO_fid = 0.87')
    compare_2plots (timestamp1= '20141215_134649', timestamp2='20141215_134741', title = 'G=5, F=3, RO_fid = 0.87')
    compare_2plots (timestamp1= '20141215_134924', timestamp2='20141215_135030', title = 'G=5, F=4, RO_fid = 0.87')
    compare_2plots (timestamp1= '20141215_135229', timestamp2='20141215_135352', title = 'G=5, F=5, RO_fid = 0.87')
    compare_2plots (timestamp1= '20141215_135556', timestamp2='20141215_135723', title = 'G=5, F=6, RO_fid = 0.87')

def compare_G5_variableF_perfRO():
    compare_2plots (timestamp1= '20141215_140947', timestamp2='20141215_141011', title = 'G=5, F=0, RO_fid = 1.00')
    compare_2plots (timestamp1= '20141215_141102', timestamp2='20141215_141139', title = 'G=5, F=1, RO_fid = 1.00')
    compare_2plots (timestamp1= '20141215_141245', timestamp2='20141215_141333', title = 'G=5, F=2, RO_fid = 1.00')
    compare_2plots (timestamp1= '20141215_141503', timestamp2='20141215_141601', title = 'G=5, F=3, RO_fid = 1.00')
    compare_2plots (timestamp1= '20141215_141745', timestamp2='20141215_141851', title = 'G=5, F=4, RO_fid = 1.00')
    compare_2plots (timestamp1= '20141215_142057', timestamp2='20141215_142214', title = 'G=5, F=5, RO_fid = 1.00')

def compare_cappellaro_varG ():
    t_stamps = ['20141215_143011','20141215_143043','20141215_143125','20141215_143221','20141215_143328','20141215_143447']
    labels = ['G=0', 'G=1', 'G=2', 'G=3', 'G=4', 'G=5']
    compare_multiple_plots (timestamps=t_stamps, labels=labels, title = 'Cappellaro protocol (F=0, RO_fid=1)')


#compare_G5_variableF_imperfRO()
#compare_G5_variableF_perfRO()
#compare_cappellaro_varG()
#compare_G3_variableF_imperfRO_recalcPhase()
#compare_G3_variableF_imperfRO_optCapp_vs_nonadaptive()
check_recalcPhase_N1sweep()

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