#######################################
###   PLOT Fig 4                    ###
#######################################

import matplotlib as mpl
#mpl.rcParams['text.usetex']=True
#mpl.rcParams['text.latex.unicode']=True
#execfile('D:/measuring/analysis/scripts/setup_analysis.py')
from matplotlib import rc, cm
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
from analysis.scripts.magnetometry import adaptive_magnetometry_analysis as data_analysis
from analysis.scripts.magnetometry.paper_figures import plots
reload(plots)
reload(magnetometry)
load_data=False
gamma_e = 28e9
def analyze_data():
    ###############
    ### Fig 4 a ###
    ### G=5 F=7 ###
    ###############

    # Cappellaro
    # in analyzed data folder: 20141215_183540_adaptive_magnetometry_N=14G=5F=7_fid0=0.87.hdf5
    N14_G5_F7=data_analysis.analyze_sweep_field(F=7,G=5,nr_periods=4,phase_update=False,CR_after_threshold=6,newer_than='20141111_172700',older_than='20141112_161817')

    #Berry adapt
    # in analyzed data folder: 20141215_192100_adaptive_magnetometry_N=14G=5F=7_fid0=0.87.hdf5
    N14_G5_F7_adapt=data_analysis.analyze_sweep_field(F=7,G=5,nr_periods=1,phase_update=True,CR_after_threshold=7,newer_than='20141121_131000',older_than='20141121_161500')
    #Berry non-adapt
    # in analyzed data folder: 20141215_194843_adaptive_magnetometry_N=14G=5F=7_fid0=0.87.hdf5
    N14_G5_F7_non_adapt=data_analysis.analyze_sweep_field(F=7,G=5,nr_periods=1,phase_update=True,CR_after_threshold=7,newer_than='20141126_103400',older_than='20141126_131200')


    ###############
    ### Fig 4 b ###
    ### G=5 F=2 ###
    ###############

    # Cappellaro
    # in analyzed data folder: 20141215_200053_adaptive_magnetometry_N=14G=5F=2_fid0=0.87.hdf5
    N14_G5_F2=data_analysis.analyze_sweep_field(F=2,G=5,nr_periods=1,phase_update=False,CR_after_threshold=7,newer_than='20141119_234700',older_than='20141120_025300')
    #Berry adapt
    # in analyzed data folder: 20141215_201910_adaptive_magnetometry_N=14G=5F=2_fid0=0.87.hdf5
    N14_G5_F2_adapt=data_analysis.analyze_sweep_field(F=2,G=5,nr_periods=1,phase_update=True,CR_after_threshold=7,newer_than='20141121_101700',older_than='20141121_131000')
    #Berry non-adapt
    # in analyzed data folder: 20141215_203241_adaptive_magnetometry_N=14G=5F=2_fid0=0.87.hdf5
    N14_G5_F2_non_adapt=data_analysis.analyze_sweep_field(F=2,G=5,nr_periods=1,phase_update=True,CR_after_threshold=7,newer_than='20141125_105600',older_than='20141125_133500')
    #Swarm
    # in analyzed data folder: 20150106_221250_adaptive_magnetometry_N=14G=5F=2_fid0=0.87.hdf5
    N14_G5_F2_swarm=analyze_sweep_field(F=2,G=5,nr_periods=1,phase_update=False,swarm_opt = True, CR_after_threshold=10,newer_than='20150106_145501',older_than='20150106_213944')

def analyze_saved_simulations (timestamp,error_bars=False,include_overhead=False):
    mgnt_exp = magnetometry.AdaptiveMagnetometry(N=14, tau0=20e-9)
    mgnt_exp.error_bars=error_bars
    mgnt_exp.load_analysis (timestamp=timestamp)
    mgnt_exp.calculate_scaling(include_overhead=include_overhead)
    return mgnt_exp


def add_scaling_plot(timestamp, ax, exp_data, label, marker_settings, color_nr,timestamp_sim=None,N_max=None,include_overhead=False):
    #adds a scaling plot to axis 'ax', loading from analyzed data with a given 'timestamp'
    #exp_data=boolena, if 'True' then data is plotted with markers and errorbars are calculated, 
    #otherwise it is considered a simulation, and plotted like a line
    #label, string for legend
    data_file = analyze_saved_simulations (timestamp=timestamp, error_bars=exp_data,include_overhead=include_overhead)
    
    if N_max:
        Nm=N_max
    else:
        Nm=len(data_file.total_time)
    print data_file.err_sensitivity[Nm-1:Nm]
    #ax.errorbar(data_file.total_time[Nm-1:Nm]*1e6, data_file.sensitivity[Nm-1:Nm],yerr=data_file.err_sensitivity[Nm-1:Nm],fmt=marker_settings,color=plots.colors[color_nr], label=label)    
        
    if exp_data: 
        #ax.fill_between (data_file.total_time[0:Nm], data_file.sensitivity[0:Nm]-data_file.err_sensitivity[0:Nm], data_file.sensitivity[0:Nm]+data_file.err_sensitivity[0:Nm], rasterized=True,color=plots.colors[color_nr+1], alpha=0.5)
        if label:
            (_, caps, _)=ax.errorbar(data_file.total_time[0:Nm], data_file.sensitivity[0:Nm], yerr=data_file.err_sensitivity[0:Nm],fmt=marker_settings,elinewidth=1.5,color=plots.colors[color_nr], label=label)
        else:
            (_, caps, _)=ax.errorbar(data_file.total_time[0:Nm], data_file.sensitivity[0:Nm], yerr=data_file.err_sensitivity[0:Nm], fmt=marker_settings,elinewidth=1.5,color=plots.colors[color_nr])    
    for cap in caps:
        cap.set_markeredgewidth(1.5)
    
    if timestamp_sim:
        sim_file = analyze_saved_simulations (timestamp=timestamp_sim, error_bars=False,include_overhead=include_overhead)
        if label:
            ax.plot (sim_file.total_time[0:Nm], sim_file.sensitivity[0:Nm], '-',color=plots.colors[color_nr], label=label,linewidth=1.5)
        else:
            ax.plot (sim_file.total_time[0:Nm], sim_file.sensitivity[0:Nm], '-',color=plots.colors[color_nr],linewidth=1.5)
    
    tau0=data_file.t0
    ymin=1e18*min(data_file.sensitivity[0:Nm])/((2*np.pi*gamma_e*tau0)**2)
    ymin_err=1e18*min(data_file.err_sensitivity[0:Nm])/((2*np.pi*gamma_e*tau0)**2)
    print '######################'
    print label
    print 'Minimum Sensitivity^2 is ', ymin,'+-', ymin_err,' nT^2/Hz'
    #print 'full y', 1e18*data_file.sensitivity[0:Nm]/((2*np.pi*gamma_e*tau0)**2)
    #print 'y err', 1e18*data_file.err_sensitivity[0:Nm]/((2*np.pi*gamma_e*tau0)**2)
    #print '1/x', 1/data_file.total_time[0:Nm]
    frq=1/data_file.total_time[0:Nm]
    frq21Hz=np.abs(frq - 21).argmin()
    #print 'index of min 21 Hz', frq21Hz
    S21Hz=1e18*data_file.sensitivity[frq21Hz]/((2*np.pi*gamma_e*tau0)**2)
    errS21Hz=1e18*data_file.err_sensitivity[frq21Hz]/((2*np.pi*gamma_e*tau0)**2)
    print '@ rep rate', frq[frq21Hz] ,' Hz' , '(Closest to 21 Hz)'
    print 'Sensitivity^2 is ',S21Hz   , '+-', errS21Hz,'nT^2/Hz'
    print 'Sensitivity^2 is ', 1e-6*S21Hz  , '+-', 1e-6*errS21Hz,'uT^2/Hz'
    print 'or'
    print 'Sensitivity is ',sqrt(S21Hz)   , '+-', sqrt(errS21Hz),'nT/sqrt(Hz)'
    print 'Sensitivity is ', sqrt(1e-6*S21Hz)  , '+-', sqrt(1e-6*errS21Hz),'uT/sqrt(Hz)'
    print '######################'


    
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.ylabel ('$V_{H}$ T ')
    #plt.show()
    return ax,tau0
def add_sensitivity_plot(timestamp, ax, exp_data, marker_settings, timestamp_sim=None,color_nr=2,N_max=None,include_overhead=False):
    #adds a scaling plot to axis 'ax', loading from analyzed data with a given 'timestamp' and in units of nT/sqrt(Hz)
    #exp_data=boolena, if 'True' then data is plotted with markers and errorbars are calculated, 
    #otherwise it is considered a simulation, and plotted like a line
    #label, string for legend
    data_file = analyze_saved_simulations (timestamp=timestamp, error_bars=exp_data,include_overhead=include_overhead)
    
    if N_max:
        Nm=N_max
    else:
        Nm=len(data_file.total_time)
    tau0=data_file.t0    
    data=1e9*data_file.sensitivity[0:Nm] /((2*np.pi*gamma_e*tau0)**2)
    ax.plot (1/data_file.total_time[0:Nm], data, marker_settings,color=plots.colors[color_nr])    
    

    
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.ylabel ('$V_{H}$ T ')
    #plt.show()
    return ax,tau0
def fig4a_plot(do_save=True):
    ################
    ###  G=5 F=2 ###
    ################
    f = plt.figure(figsize=(3.5,2))
    ax1 = f.add_subplot(1,1,1)
    
    # V2-2015-01-11
    #ax1,tau0 = add_scaling_plot (timestamp = '20141215_200053', exp_data=True, ax=ax1, label = 'Cappellaro', marker_settings='o', color_nr=0,timestamp_sim='20141118_195110' )
    #ax1,tau0 = add_scaling_plot (timestamp = '20141215_201910', exp_data=True, ax=ax1, label = 'Berry adaptive', marker_settings='^', color_nr=2,timestamp_sim=None,N_max=9)
    #ax1,tau0 = add_scaling_plot (timestamp = '20141215_203241', exp_data=True, ax=ax1, label = 'Berry non-adaptive', marker_settings='v', color_nr=4,timestamp_sim='20141118_192523',N_max=9)
    

    # Cappellaro G=5 F=7
    ax1,tau0 = add_scaling_plot (timestamp = '20141215_200053', exp_data=True, ax=ax1, label = 'Limited-adaptive (G=5 F=2)', marker_settings='o', color_nr=2,timestamp_sim='20150420_195041')#'20150420_171446')#'20141118_195110' )

    # Berry non adapt G5 F2
    ax1,tau0 = add_scaling_plot (timestamp = '20141215_203241', exp_data=True, ax=ax1, label = 'Non-adaptive (G=5,F=2)', marker_settings='v', color_nr=4,N_max=9,timestamp_sim='20141118_192523')
    # Berry non adapt G5 F7
    ax1,tau0 = add_scaling_plot (timestamp = '20141215_194843', exp_data=True, ax=ax1, label = 'Non-adaptive (G=5 F=7)', marker_settings='v', color_nr=0,timestamp_sim='20141118_135425',N_max=12)
    # Swarm G5 F2
    ax1,tau0 = add_scaling_plot (timestamp = '20150106_221250', exp_data=True, ax=ax1, label = 'Adaptive (G=5,F=2)', marker_settings='^',color_nr=6,timestamp_sim=None)
    T=np.linspace(1e-7,1e-3,2001)/tau0
    ax1.plot(T*tau0,tau0*np.pi**2/T,'k--')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel ('Total phase accumulation time T (s)')
    #ax.set_ylabel ('sensitivity [$\mu$T$^2$*Hz$^{-1}$]', fontsize=15)
    ax1.set_ylabel ('$V_H$ T (rad$^2$ Hz$^{-1}$)')
    ax1.legend(loc=3)

    y1_max=1e-6
    y1_min=1e-10
    ax1.set_ylim([y1_min,y1_max])

    
    # Double y axis: Sensitivity nT sqrtHz
    ax2 = ax1.twinx()
    ax2.set_yscale('log')
    #ax2,tau0 = add_sensitivity_plot (timestamp = '20141215_203241', exp_data=True, ax=ax2,  marker_settings='.', color_nr=2,N_max=1)#,timestamp_sim='20141118_192523')
    ax2.set_ylim([1e18*y1_min/((2*np.pi*gamma_e*tau0)**2),1e18*y1_max /((2*np.pi*gamma_e*tau0)**2)])
    ax2.set_ylabel('Sensitiviy$^2$ ($nT^2$ $Hz^{-1}$)')


    plt.show()
    if do_save:
        f.savefig(r'M:\tnw\ns\qt\Diamond\Eigenpapers\15-Adaptive DC Magnetometry\Figures\Fig4a.pdf', bbox_inches='tight')

def fig4b_plot(do_save=True):
    ################
    ###  G=5 F=7 ###
    ################
    f = plt.figure(figsize=(3.5,2))
    ax1 = f.add_subplot(1,1,1)
    
    # V2-2015-01-11
    #ax1,tau0 = add_scaling_plot (timestamp = '20141215_183540', exp_data=True, ax=ax1, label = 'Cappellaro', marker_settings='o', color_nr=0,timestamp_sim='20141112_182801' )
    #ax1,tau0 = add_scaling_plot (timestamp = '20141215_192100', exp_data=True, ax=ax1, label = 'Berry adaptive', marker_settings='^', color_nr=2,timestamp_sim='20141118_134405')
    #ax1,tau0 = add_scaling_plot (timestamp = '20141215_194843', exp_data=True, ax=ax1, label = 'Berry non-adaptive', marker_settings='v', color_nr=4,timestamp_sim='20141118_135425')
    #ax1,tau0 = add_scaling_plot (timestamp = '20150106_221250', exp_data=True, ax=ax1, label = 'Swarm (G=5,F=2)', marker_settings='^',color_nr=6,timestamp_sim=None)

    

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel ('Total phase accumulation time T ($\mu$s)')
    #ax.set_ylabel ('sensitivity [$\mu$T$^2$*Hz$^{-1}$]')
    ax1.set_ylabel ('$V_H$ T')



    # Double x axis: Longest interaction time tau
    tempax = ax1.twinx()
    ax2 = tempax.twiny()
    ax2.set_xscale('log')

    n=np.array([2,7,12])
    F=5
    G=7
    T=[]
    L=[]
    t0=20e-9
    T.append(t0*(G*(2**(n)-1)+F*(2**(n)-1-n)))
    for i in np.arange(len(n)):
        #print i
        L.append(str(2**n[i]*20*1e-3))
    print L

    ax2.set_xlabel(r'Longest interaction time $\tau^N$ ($\mu s$) ')
    ax2.set_xlim(1e-1, 1e4)
    # NOTE: check this conversion, should depend on G and F!!!!
    ax2.set_xticks([  5.20000000e-01 ,  2.97800000e1, 9.81600000e2])
    ax2.set_xticklabels(L)
    ax2.axvspan(1.0817e3, 1e4, alpha=0.5, color='LightGrey')
    
    # Double y axis: Sensitivity nT sqrtHz
    #ax3 = ax1.twiny()
    tempax.set_yscale('log')
    
    tempax.set_ylim(1e-9, 1e-6)
    ax1Ys = ax1.get_yticks()
    print tau0*gamma_e
    ax2Ys = []
    for Y in ax1Ys:
        ax2Ys.append(int(1e9*Y /((2*np.pi*gamma_e*tau0)**2)))
    
    tempax.set_yticks(ax1Ys)
    tempax.set_ybound(ax1.get_ybound())
    tempax.set_yticklabels(ax2Ys)
    tempax.set_ylabel('Sensitiviy ($nT$ $Hz^{-0.5}$)')
    ax1.legend()

    plt.show()
    if do_save:
        f.savefig(r'M:\tnw\ns\qt\Diamond\Eigenpapers\15-Adaptive DC Magnetometry\Figures\Fig4b.pdf', bbox_inches='tight')



def fig4c_plot(do_save=True):
    f = plt.figure(figsize=(3.5,2))
    ax1 = f.add_subplot(1,1,1)
    x1_max=1e0
    x1_min=1e-3
    y1_max=1e-2
    y1_min=1e-8
    
    # For zoom:
    #ax1.set_xlim([4e4,2e5])
    #y1_max=3.2e-8
    #y1_min=2.6e-8
    
    
    #ax1.set_xlim([1e7,2e8])
    #y1_max=3e-5
    #y1_min=2.5e-5
    
    ax1.set_ylim([y1_min,y1_max])
    ax1.set_xlim([x1_min,x1_max])
    # V2-2015-01-11
    #ax1,tau0 = add_scaling_plot (timestamp = '20141215_200053', exp_data=True, ax=ax1, label = 'Cappellaro', marker_settings='o', color_nr=0,timestamp_sim='20141118_195110',include_overhead=True )

    # Berry non adapt G5 F7
    ax1,tau0 = add_scaling_plot (timestamp = '20141215_194843', exp_data=True, ax=ax1, label = 'Non-adaptive (G=5, F=7)', marker_settings='v', color_nr=0,timestamp_sim='20141118_135425',include_overhead=True,N_max=12)
    # Swarm G5 F2
    ax1,tau0 = add_scaling_plot (timestamp = '20150106_221250', exp_data=True, ax=ax1, label = 'Adaptive (G=5, F=2)', marker_settings='^', color_nr=6,timestamp_sim=None,include_overhead=True)

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel ('Total phase accumulation time T + $T_{overhead}$ (s)')

    #ax1.set_xlabel ('Repetition rate [Hz]')
    #ax.set_ylabel ('sensitivity [$\mu$T$^2$*Hz$^{-1}$]', fontsize=15)
    ax1.set_ylabel ('$V_H$ (T + $T_{overhead}$) (rad$^2$ Hz$^{-1}$)')

    
    
    # Double y axis: Sensitivity nT sqrtHz
    ax2 = ax1.twinx()
    ax2.set_yscale('log')
    #ax2,tau0 = add_sensitivity_plot (timestamp = '20141215_203241', exp_data=True, ax=ax2,  marker_settings='.', color_nr=2,N_max=1)#,timestamp_sim='20141118_192523')
    ax2.set_ylim([1e18*y1_min/((2*np.pi*gamma_e*tau0)**2),1e18*y1_max /((2*np.pi*gamma_e*tau0)**2)])
    ax2.set_ylabel('Sensitiviy$^2$ ($nT^2$ $Hz^{-1}$)')
    
    ax3 = ax1.twiny()
    ax3.set_xscale('log')
    #ax2,tau0 = add_sensitivity_plot (timestamp = '20141215_203241', exp_data=True, ax=ax2,  marker_settings='.', color_nr=2,N_max=1)#,timestamp_sim='20141118_192523')
    ax3.set_xlim([1/float(x1_min),1/float(x1_max)])
    ax3.set_xlabel('Repetition rate (Hz)')
    ax1.tick_params(axis='x')
    ax1.tick_params(axis='y')
    ax2.tick_params(axis='y')
    ax3.tick_params(axis='x')
    #ax1.legend(loc=3)

    plt.show()
    if do_save:
        f.savefig(r'M:\tnw\ns\qt\Diamond\Eigenpapers\15-Adaptive DC Magnetometry\Figures\Fig4c.pdf', bbox_inches='tight')    

def failed_CR_checks(do_save=True):
    perc_cr_failed=np.arange(0)
    perc_cr_failed_averaged_per_datapoint  = np.arange(0)
    timestamps=['20141215_203241','20141215_194843','20150106_221250']
    for t in timestamps:
        testdata=analyze_saved_simulations(timestamp=t)
        for i in np.arange(testdata.N-2):
            perc_cr_failed=np.append(perc_cr_failed,100*testdata.results_dict[str(i+2)]['nr_discarded_elements']/float(testdata.repetitions))
            perc_cr_failed_averaged_per_datapoint=np.append(perc_cr_failed_averaged_per_datapoint,[mean(100*testdata.results_dict[str(i+2)]['nr_discarded_elements']/float(testdata.repetitions))])
        
        #print perc_cr_failed
        #print 100*testdata.results_dict[str(i+2)]['nr_discarded_elements']/float(testdata.repetitions)
        #print testdata.N    
        #print perc_cr_failed
    
    f = plt.figure(figsize=(3.5,2))
    ax = f.add_subplot(1,1,1)

    ax.hist(perc_cr_failed,bins=10,color='RoyalBlue')
    ax.set_yticks([0,25,50,75,100])
    ax.set_xticks([0,25,50,75,100])
    ax.set_ylabel('Occurence')
    ax.set_xlabel('Percentage of rejected runs')
    if do_save:
        f.savefig(r'M:\tnw\ns\qt\Diamond\Eigenpapers\15-Adaptive DC Magnetometry\Figures\Failed_CR_checks.pdf', bbox_inches='tight')    
    
    f = plt.figure(figsize=(3.5,2))
    ax = f.add_subplot(1,1,1)
  
    ax.hist(perc_cr_failed_averaged_per_datapoint,bins=10,color='RoyalBlue')
    #ax.set_yticks([0,25,50,75,100])
    ax.set_xticks([0,25,50,75,100])
    ax.set_ylabel('Occurence')
    ax.set_xlabel('Percentage of rejected runs')
    if do_save:
        f.savefig(r'M:\tnw\ns\qt\Diamond\Eigenpapers\15-Adaptive DC Magnetometry\Figures\Failed_CR_checks_averaged_per_datapoint.pdf', bbox_inches='tight')   
    
fig4a_plot(do_save=False)
#fig4b_plot(do_save=False)
#fig4c_plot(do_save=True)
#failed_CR_checks(do_save=True)


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