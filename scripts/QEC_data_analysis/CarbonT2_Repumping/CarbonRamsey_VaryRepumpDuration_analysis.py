import numpy as np
import os,sys

sys.path.append("/measuring/")
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(fit) 
reload(toolbox)

"""
this script analyzes a data set in which the combined system of electronic system and a single carbon was initialized in the state |-1>|-x> via MBI.
We then apply a repumping laser which transfers the electronic spin back to |0> for a varying duration.
The carbon spin in this experiment has a measured hyperfine coupling of -36.8kHz (see also the QEC labbook/PPt slides, Carbon 1)

Norbert Kalb
2014
"""

def CosineSum_MBI_data(timestamp=None, measurement_name = ['adwindata'], ssro_calib_timestamp =None,
        frequency = [1,1], offset =0.5, amplitude =[ 0.5,0.5],  phase =[0,0], 
        fixed = [], 
        plot_fit = False, do_print = False, show_guess = True):
    ''' 
    Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    [freq, offset, Amplitude, phase] 
    '''


    if timestamp == None:
        timestamp, folder   = toolbox.latest_data('CarbonR',return_timestamp =True)
    else: 
        folder = toolbox.data_from_time(timestamp) 

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'

    fit_result = [None]

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder)
    ax = a.plot_results_vs_sweepparam(ret='ax')
    x = a.sweep_pts.reshape(-1)[:]
    y= a.p0.reshape(-1)[:]

    p0, fitfunc, fitfunc_str = common.fit_sum_2cos(offset,amplitude[0],frequency[0],phase[0],amplitude[1],frequency[1],phase[1]) 
    if show_guess:
        ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=fixed)
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax, 
                plot_data=False,print_info = True)
    
    fit.write_to_file(fit_result,folder,fitname = 'Sum of cosine fit') 


    ## plot data and fit as function of total time

    plt.savefig(os.path.join(folder, 'CosineSumFit.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'CosineSumFit.png'),
    format='png')
    return fit_result,folder

    # return fit_result

def fit_decaying_cos(g_f, g_a, g_A, g_phi,g_t, *arg):
    fitfunc_str = 'A *exp(-x/t) *cos(2pi * (f*x + phi/360) ) + (1-exp(-x/t)a)'

    f = fit.Parameter(g_f, 'f')
    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    phi = fit.Parameter(g_phi, 'phi')
    t   = fit.Parameter(g_t, 't')
    print 'guessed frequency is '+str(g_f)
    p0 = [f, a, A,phi,t]

    def fitfunc(x):
        return (1-np.exp(-x/t()))*a() + A()*np.exp(-x/t()) * np.cos(2*np.pi*( f()*x + phi()/360.))

    return p0, fitfunc, fitfunc_str


def FitSpinPumpTimes(fit_res_list,spin_pump_list,folder):
    """
    this function takes the fitted amplitudes for the population in ms=-1 and ms=0 and returns them as a plot.
    """

    spin_pump_list
    x=np.ones(len(spin_pump_list))
    #intialize array for ms=0
    y0=np.ones(len(spin_pump_list))
    y0_u=np.ones(len(spin_pump_list))
    #intialize array for ms=-1
    y1=np.ones(len(spin_pump_list))
    y1_u=np.ones(len(spin_pump_list))


    #extract relevant values from the fit results.
    for i,sp in enumerate(spin_pump_list):
        x[i]=sp
        fit_res=fit_res_list[i]
        y0[i]=np.absolute(fit_res['params_dict']['A'])
        y0_u[i]=np.absolute(fit_res['error'][1])
        y1[i]=np.absolute(fit_res['params_dict']['B'])
        y1_u[i]=np.absolute(fit_res['error'][4])

    print x,y0
    fig=plt.figure(10)
    ax=plt.subplot(111)
    #work around, because the function 'errorbar is a bit bugged...'
    plt.rc('lines',**{'linestyle':'None'})
    plt.errorbar(x,y0,y0_u,color='blue',marker='o',fmt='')
    plt.errorbar(x,y1,y1_u,color='red',marker='o',fmt='')
    plt.axis([0,60,0,1])

    p0, fitfunc, fitfunc_str = fit_decaying_cos(1/20.,0.,0.02,-90.,1./0.037)


    #fit an exponential starting at 0 repumping time and an exponent of 1.
    fit_result=fit.fit1d(x,y0,None,p0=p0,fitfunc=fitfunc,do_print=True,ret=True,print_info=False,fixed=[3])
    plot.plot_fit1d(fit_result,np.linspace(x[0],x[-1],1001),ax=ax,print_info=False,plot_data=False,linestyle='-b')


    p0, fitfunc, fitfunc_str = common.fit_general_exponential(0.,0.5,0.,20.,1.)
    print fitfunc_str

    fit_result=fit.fit1d(x,y1,None,p0=p0,fitfunc=fitfunc,do_print=True,print_info=False,ret=True,fixed=[0,2,4])
    plot.plot_fit1d(fit_result,np.linspace(x[0],x[-1],1001),ax=ax,print_info=False,plot_data=False,linestyle='-r')


    #fit with a double exponential:
    p0, fitfunc, fitfunc_str = common.fit_double_exp_decay_with_offset(-0.2,0.5,5.,0.5,20.)
    print fitfunc_str
    fit_result=fit.fit1d(x,y1,None,p0=p0,fitfunc=fitfunc,do_print=True,print_info=False,ret=True,fixed=[0])
    plot.plot_fit1d(fit_result,np.linspace(x[0],x[-1],1001),ax=ax,print_info=False,plot_data=False,linestyle='-g')


    plt.title('Sample_111_No1_C13_C1_Ramsey_contrast_over_spinpumping_duration')
    plt.xlabel('Spin pumping duration (us)')
    plt.ylabel('Amplitude')
    plt.legend(['f_0','f_1'])

    plt.savefig(os.path.join(folder,'SpinPumpData.pdf'),format='pdf')
    plt.savefig(os.path.join(folder,'SpinPumpData.png'),format='png')








timestamp_list=['20141202_211228','20141202_213444','20141202_215711','20141202_221923','20141202_224216','20141202_230523','20141202_232856','20141202_235306','20141203_001820']
ssro_stamp='20141202_205422'
spin_pump_list=[0,2,5,7,10,20,30,40,50]

fit_result_list=[]
for i,timestamp in enumerate(timestamp_list):
    cur_Fit_result,last_folder=CosineSum_MBI_data(timestamp=timestamp,
        ssro_calib_timestamp=ssro_stamp,
        frequency=[430e3,468e3],
        offset=0.5,
        amplitude=[0.2,0.3],
        fixed=[],
        plot_fit=True,
        do_print=False,
        show_guess=False)
    fit_result_list.append(cur_Fit_result)

FitSpinPumpTimes(fit_result_list,spin_pump_list,last_folder)