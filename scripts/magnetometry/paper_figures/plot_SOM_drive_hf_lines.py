from analysis.lib.m2.ssro import sequence
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
reload(sequence)
from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, ramsey,rabi
from analysis.lib.tools import plot
from analysis.lib.math import error
from matplotlib import rc, cm

###########################
##### Plot Fig 1 Ramsey ###
###########################

#folder=r'K:\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20140918\222845_PulsarMBIElectronRamsey_gretel-sil10mbi_eramsey_pump_pm1_check_leakage'
foldermIm1=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141127\170640_adptv_estimation_rabi_awg_pulse_msm1_lines_AWG'
foldermI0=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141127\170751_adptv_estimation_rabi_awg_pulse_ms0_lines_AWG'
foldermIp1=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141127\171050_adptv_estimation_rabi_awg_pulse_msp1_lines_AWG'
folder_all_lines=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141127\170247_adptv_estimation_rabi_awg_pulse_all_lines_AWG'
#folder20kHz=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141219\110821_adptv_estimation_ramsey_sweep_timedet_20kHz'
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
  
for i in range(len(tableau20)):  
    r, g, b = tableau20[i]  
    tableau20[i] = (r / 255., g / 255., b / 255.) 

colors=[tableau20[0],tableau20[2],tableau20[4],tableau20[6]]

def plot_single_rabi_curve(folder,ax,color,label,tot,utot,tot_mean):
    a=sequence.MagnetometrySequenceAnalysis(folder)
    RO,uRO=a.get_magnetometry_phase_calibration(name='adwindata',ssro_calib_folder=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141127\162130_AdwinSSRO_SSROCalibration_Gretel_sil10')
    
    guess_f1 = 20e-3 #in GHz
    guess_A1 = 0.5
    guess_phi1 = np.pi
    guess_tau = 96
    guess_a = 0.5
    
    a.sweep_pts=a.sweep_pts
    print a.sweep_pts
    #ax.set_ylim([0.0,1.05])
    amp = (max(a.p0)-min(a.p0))/2.
    offset = (max(a.p0)+min(a.p0))/2.
    freq = 1/8000.
    decay = 10000.
    fit_result = fit.fit1d(a.sweep_pts, a.p0, rabi.fit_rabi_damped_exp,
        freq,amp, offset,decay,
         fixed=[],
        do_print=False, ret=True)
    A=fit_result['params_dict']['A']*2
    F=fit_result['params_dict']['f']
    eF=fit_result['error_dict']['f']
    #print fit_result
    print min(a.p0),A
    print 'frq',F*1e6, '+-',eF,'kHz'
    x_fit=np.linspace(a.sweep_pts[0],a.sweep_pts[-1],501)
    y_fit=fit_result['fitfunc'](x_fit)
    ax.plot(x_fit*1e-3,y_fit,'-',color='Grey',linewidth=2)
    ax.errorbar(a.sweep_pts*1e-3,a.p0,yerr=a.u_p0,color=color,fmt='o',label=label)

    if len(tot)==0:
        tot=(a.p0-mean(a.p0))
        utot=a.u_p0
        tot_mean=mean(a.p0)
    else:
        tot=tot+(a.p0-mean(a.p0))
        utot=utot+a.u_p0
        tot_mean=tot_mean+mean(a.p0)
    print max(a.u_p0)    
    return tot,a.sweep_pts,utot,tot_mean,A
def plot_hf_lines(do_save=False):


    folders=[foldermIm1,foldermI0,foldermIp1]
    f_names=['m1','zero','p1']

    labels=['mI=-1','mI=0','mI=+1']

    tot=[]
    utot=[]
    tot_mean=0
    A_tot=0
    for i in np.arange(len(folders)):
        fig2 = plt.figure(figsize=(3,1.25))
        fig2.clf()
        ax2 = fig2.add_subplot(111)
        tot,sweep_pts,utot,tot_mean,A=plot_single_rabi_curve(folders[i],ax2,colors[i],labels[i],tot,utot,tot_mean)
        A_tot+=A
        ax2.set_xlabel('Pulse duration (us)')#,fontsize=24)
        ax2.set_ylabel('P($m_s =0$)')#,fontsize=24)
        ax2.tick_params(axis='x')#, labelsize=18)
        ax2.tick_params(axis='y')#, labelsize=18)
        ax2.yaxis.set_ticks([0,0.5,1])
        ax2.set_xlim([0,14])
        ax2.set_ylim([0,1.1])
        ax2.legend(loc=4)
        if do_save:
            folder_name=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\analyzed data\hf_lines_'+str(f_names[i])+'.pdf'
            fig2.savefig(folder_name, bbox_inches='tight')

    print 'A_tot',A_tot
    fig2 = plt.figure(figsize=(6,4))
    fig2.clf()
    ax2 = fig2.add_subplot(111)
    plot_single_rabi_curve(folder_all_lines,ax2,colors[-1],'All_lines',tot,utot,tot_mean)
    #ax2.errorbar(sweep_pts,tot+0.5,yerr=utot,color='k',fmt='o',label='Total')
    ax2.set_xlabel('Pulse duration (us)')#,fontsize=24)
    ax2.set_ylabel('P($m_s =0$)')#,fontsize=24)
    ax2.tick_params(axis='x')#, labelsize=18)
    ax2.tick_params(axis='y')#, labelsize=18)
    ax2.set_xlim([0,14])
    ax2.set_ylim([0,1.1])

    ax2.legend()
    if do_save:
        fig2.savefig(r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\analyzed data\all_hf_lines.pdf', bbox_inches='tight')

def plot_desr(do_save=False):
    folder=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141014\130500_PulsarDarkESR_FM_Gretel_sil10'

    ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration')
    print ssro_calib_folder
    
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

    x = a.sweep_pts # convert to MHz
    y = a.p0.reshape(-1)[:]
    # ax.plot(x,y)
    fig2 = plt.figure(figsize=(2,1.25))
    fig2.clf()
    ax2 = fig2.add_subplot(111) 
    #a.plot_result_vs_sweepparam(ret=True, name='ssro',ax=ax2)
    ax2.set_xlim([-3,3])
    ax2.set_ylim(0.5,1)
    ax2.yaxis.set_ticks([0.5,0.75,1])
    #ax2.xaxis.set_ticks([0.5,0.75,1])
    guess_x0= 2.845334
    guess_offset = 1
    guess_Nsplit = 2.18e-3
    guess_sigma = 0.0001

    guess_A_min1 = 0.4
    guess_A_plus1 = 0.25
    guess_A_0 = 0.3
 

    A_min1 = fit.Parameter(guess_A_min1, 'A_min1')
    A_plus1 = fit.Parameter(guess_A_plus1, 'A_plus1')
    A_0 = fit.Parameter(guess_A_0, 'A_0')
    o = fit.Parameter(guess_offset, 'o')
    x0 = fit.Parameter(guess_x0, 'x0')
    sigma = fit.Parameter(guess_sigma, 'sigma')
    Nsplit = fit.Parameter(guess_Nsplit, 'Nsplit')
    def fitfunc(x):

        return o() - np.abs(A_min1())*np.exp(-((x-(x0()-Nsplit()))/sigma())**2) \
                - np.abs(A_plus1())*np.exp(-((x-(x0()+Nsplit()))/sigma())**2) \
                - np.abs(A_0())*np.exp(-((x-x0())/sigma())**2) \
    

    print x
    fit_result = fit.fit1d(x, y, None, p0 = [A_min1, A_plus1, A_0, sigma, o, x0, Nsplit],
        fitfunc = fitfunc, do_print=True, ret=True, fixed=[])  
    print fitfunc(x)
    x_fit=np.linspace(a.sweep_pts[0],a.sweep_pts[-1],2001)
    y_fit=fit_result['fitfunc'](x_fit)
    
    
    ax2.plot(1e3*(x_fit-2.845334),y_fit,color='Grey',linewidth=2)
    ax2.errorbar(1e3*(a.sweep_pts-2.845334),a.p0,yerr=a.u_p0,color=colors[3],fmt='o')
    ax2.set_ylabel('P($m_s =0$)')
    ax2.set_xlabel('Pulse detuning (MHz)')
    if do_save:
        fig2.savefig(r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\analyzed data\desr.pdf', bbox_inches='tight')

## Interesting DESR measurements
#M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141014\130500_PulsarDarkESR_FM_Gretel_sil10
plot_hf_lines(do_save=False)
#plot_desr(do_save=True)