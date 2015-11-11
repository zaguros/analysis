from analysis.lib.m2.ssro import sequence
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
reload(sequence)
from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, ramsey
from analysis.lib.tools import plot
from analysis.lib.math import error
from matplotlib import rc, cm

###########################
##### Plot Fig 1 Ramsey ###
###########################

#folder=r'K:\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20140918\222845_PulsarMBIElectronRamsey_gretel-sil10mbi_eramsey_pump_pm1_check_leakage'
folder0kHz=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141218\161012_adptv_estimation_ramsey_sweep_timedet_0kHz'
folder20kHz=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141218\161315_adptv_estimation_ramsey_sweep_timedet_20kHz'
#folder20kHz=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141219\110821_adptv_estimation_ramsey_sweep_timedet_20kHz'

def plot_single_ramsey_curve(folder,ax,color,label):
    a=sequence.MagnetometrySequenceAnalysis(folder)
    RO,uRO=a.get_magnetometry_phase_calibration(name='adwindata')
    
    guess_f1 = 20e-3 #in GHz
    guess_A1 = 0.5
    guess_phi1 = np.pi
    guess_tau = 96
    guess_a = 0.5
    
    a.sweep_pts=a.sweep_pts*1e6
    #ax.set_ylim([0.0,1.05])

    fit_result = fit.fit1d(a.sweep_pts, a.p0, ramsey.fit_ramsey_gaussian_decay,
        guess_tau, guess_a, (guess_f1, guess_A1, guess_phi1),
        #(guess_f3, guess_A3, guess_phi3),
         fixed=[],
        do_print=False, ret=True)
    
    ax.errorbar(a.sweep_pts,a.p0,yerr=a.u_p0,color=color,fmt='o',label=label)
    x_fit=np.linspace(a.sweep_pts[0],a.sweep_pts[-1],501)
    y_fit=fit_result['fitfunc'](x_fit)
    ax.plot(x_fit,y_fit,'-',color='Grey',linewidth=2)

def plot_ramsey(do_save=False):


    folders=[folder0kHz,folder20kHz]
    colors=['RoyalBlue','Crimson']
    labels=['0 kHz','20 kHz']
    fig2 = plt.figure(figsize=(6,4))
    fig2.clf()
    ax2 = fig2.add_subplot(111)

    for i in np.arange(len(folders)):
        plot_single_ramsey_curve(folders[i],ax2,colors[i],labels[i])
    ax2.set_xlabel('Free evolution time (us)')#,fontsize=24)
    ax2.set_ylabel('P($m_s =0$)')#,fontsize=24)
    ax2.tick_params(axis='x')#, labelsize=18)
    ax2.tick_params(axis='y')#, labelsize=18)
    #ax2.set_xlim([0,200])
    ax2.set_ylim([0,1])
    ax2.legend()
    if do_save:
        fig2.savefig(r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\analyzed data\Ramsey.pdf', bbox_inches='tight')



###########################
##### Plot Fig 1 SSRO   ###
###########################
import h5py
from matplotlib import rc

def plot_ssro(do_save=False):
    f=h5py.File(r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141119\160309_AdwinSSRO_SSROCalibration_Gretel_sil10\analysis.hdf5')
    #print f.attrs
    cpsh_ms0= f['cpsh/ms0'].value
    cpsh_ms1= f['cpsh/ms1'].value
    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111)


    ax.hist(cpsh_ms0, np.arange(max(cpsh_ms0)+2)-0.5, align='mid',
                        normed=True,facecolor='Grey',edgecolor='None',alpha=0.3,linewidth=0) # , stacked=True)
    ax.hist(cpsh_ms1, np.arange(max(cpsh_ms1)+2)-0.5, align='mid', 
                        normed=True,facecolor='Grey',edgecolor='None',alpha=0.3,linewidth=0) # , stacked=True)  

    #ax.hist(cpsh_ms1, np.arange(max(cpsh_ms1)+2)-0.5, align='mid', label='$m_s = \pm 1$\n F=0.99',
    #                    normed=True,facecolor='None',edgecolor='Crimson',alpha=1,linewidth=2) # , stacked=True)    
    #ax.hist(cpsh_ms0, np.arange(max(cpsh_ms0)+2)-0.5, align='mid', label='$m_s = 0$ \n F=0.88' ,
    #                    normed=True,facecolor='None',edgecolor='RoyalBlue',alpha=1,linewidth=2) # , stacked=True)    

    ax.hist(cpsh_ms1, np.arange(max(cpsh_ms1)+2)-0.5, align='mid', histtype='step',label='$m_s = \pm 1$\n F=0.99',
                        normed=True,facecolor='None',edgecolor='Crimson',alpha=1,linewidth=2) # , stacked=True)    
    ax.hist(cpsh_ms0, np.arange(max(cpsh_ms0)+2)-0.5, align='mid', histtype='step',label='$m_s = 0$ \n F=0.88' ,
                        normed=True,facecolor='None',edgecolor='RoyalBlue',alpha=1,linewidth=2) # , stacked=True)    

    ax.set_xlabel('Photon number',fontsize=24)
    ax.set_ylabel('Probability',fontsize=24)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    #ax.set_title(self.default_plot_title + title_suffix)
    ax.set_xlim(-0.75, max(cpsh_ms1-9)+0.5)
    ax.legend(prop={'size':18})
    f.close()
    if do_save:
        fig.savefig(r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\analyzed data\SSRO.pdf', bbox_inches='tight')

########################
### Inset: Linescan ####
########################

def plot_linescan(do_save=False):
    spctr_data=np.loadtxt(r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\raw_data\20141113\133506_LaserScansGreenRepump_Gretel_SIL10_LT1_Red_Yellow\133506_LaserScansGreenRepump_Gretel_SIL10_LT1_Red_Yellow.dat')
    #plt.plot(spctr_data[:,1],spctr_data[:,2])
    #print shape(spctr_data[:,1])
    j=0
    linescan_nr=18
    counts=[]
    freq=[]
    print len(spctr_data[:,3])
    for i in np.arange(len(spctr_data[:,3])):
        if int(spctr_data[i,3]) == linescan_nr:
            counts.append(spctr_data[i,2])
            freq.append(spctr_data[i,1])
            j+=1
    fig = plt.figure(figsize=(5,3))
    ax = fig.add_subplot(111)        
    ax.plot(freq,counts, color='Royalblue',linewidth=2)
    plt.xlabel('Laser detuning (GHz)')
    plt.ylabel('Intensity (kcnt/s)')
    plt.ylim([0,6000])
    plt.xlim([58,65])
    if do_save:
        fig.savefig(r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\analyzed data\linescan.pdf', bbox_inches='tight')


plot_ssro(do_save=False)
plot_linescan(do_save=False)
plot_ramsey(do_save=False)