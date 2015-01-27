import numpy as np
import os
from analysis.lib.tools import toolbox; reload(toolbox)
from analysis.lib.tools import plot; reload(plot)
from analysis.lib.fitting import fit, common, ramsey;reload(common); reload(fit) 
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt

''' Script to analyze and fit the electron Ramsey data with 13C initialization, 
    as well as the simultaneously taken tomography data THT 141201 ''' 

''' Plan: 
- We then want to fit this data with a model that gives us the initialization fidelity
'''
### Data folders
    ### First data set
# timestamp_list_noinit = ['20141201_192629', '20141201_194015', '20141201_195413', '20141201_201324', '20141201_202722', '20141201_204121',
#                         '20141201_211401', '20141201_212807'] 
# timestamp_list_up =     ['20141201_191556', '20141201_192955', '20141201_194338', '20141201_200247', '20141201_201647', '20141201_203051',
#                         '20141201_210314', '20141201_211722'] 
# timestamp_list_down =   ['20141201_192115', '20141201_193458',  '20141201_194859', '20141201_200759', '20141201_202211', '20141201_203606',
#                         '20141201_210831', '20141201_212240'] 

    ### Second data set


def get_and_plot_data(timestamp_list, ssro_calib_folder):

    for ii, timestamp in enumerate(timestamp_list):

        folder = toolbox.data_from_time(timestamp)
        
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        
        if ii==0:
            ax = a.plot_results_vs_sweepparam(ret='ax', markersize = 4, save=False, fmt = 'o-')
            ax.set_ylim(-0.05,1.05)
            ax.set_xlim(a.sweep_pts[0],a.sweep_pts[-1])
            ax.axhspan(0,1,fill=False,ls='dotted')
            ax.axhspan(0,0.5,fill=False,ls='dotted')
        else: 
            a.plot_results_vs_sweepparam(ax= ax,  markersize = 4, save=False, fmt = 'o-')
        
        if ii==0:
            p0_sum     = a.p0
            u_p0_sum   = a.u_p0**2 
        else:
            p0_sum     = p0_sum + a.p0
            u_p0_sum   = u_p0_sum + a.u_p0**2
    
    a.p0   = p0_sum/len(timestamp_list) 
    a.u_p0 = (u_p0_sum**0.5) / len(timestamp_list) 
    
    a.x    = a.sweep_pts.reshape(-1)[:]
    a.p0   = a.p0.reshape(-1)[:]
    a.u_p0 = a.u_p0.reshape(-1)[:] 

    return a, folder

def get_and_plot_data_nuc_RO(timestamp_list_pos, timestamp_list_neg, ssro_calib_folder):

    for ii, timestamp in enumerate(timestamp_list_pos):

        folder = toolbox.data_from_time(timestamp)
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        a.x    = a.sweep_pts.reshape(-1)[:]
        a.sweep_pts = range(len(a.p0))
        
        folder = toolbox.data_from_time(timestamp_list_neg[ii])
        b = mbi.MBIAnalysis(folder)
        b.get_sweep_pts()
        b.get_readout_results(name='adwindata')
        b.get_electron_ROC(ssro_calib_folder)
        
        ### Combine positive/negative data
        a.p0 = ( (2*a.p0-1) - (2*b.p0-1))/2.
        a.u_p0 =  (a.u_p0**2 + a.u_p0**2)**0.5 

        # if ii==0:
        #     ax = a.plot_results_vs_sweepparam(ret='ax', markersize = 4, save=False, fmt = 'o')
        #     ax.set_ylim(-1.05,1.05)
        #     ax.set_xlim(a.sweep_pts[0]-0.1,a.sweep_pts[-1]+0.1)
        #     ax.axhspan(0,1,fill=False,ls='dotted')
        #     ax.axhspan(-1,0,fill=False,ls='dotted')
        # else: 
        #     a.plot_results_vs_sweepparam(ax= ax,  markersize = 4, save=False, fmt = 'o')
        #     ax.set_ylim(-1.05,1.05)

        if ii==0:
            p0_sum     = a.p0
            u_p0_sum   = a.u_p0**2 
        else:
            p0_sum     = p0_sum + a.p0
            u_p0_sum   = u_p0_sum + a.u_p0**2
    
    a.p0   = p0_sum/len(timestamp_list_pos)
    a.u_p0 = (u_p0_sum**0.5) / len(timestamp_list_pos) 
    
    a.p0   = a.p0.reshape(-1)[:]
    a.u_p0 = a.u_p0.reshape(-1)[:] 
    return a, folder

def get_data_timestamps(older_than, nr_of_repetitions, carbon = 'C1'):

    ### Initialize data lists
    data_types = [carbon+'_noInit', carbon+'_down', carbon+'_up',
                  '_down_'+carbon+'_neg', '_down_'+carbon+'_pos', '_up_'+carbon+'_neg', '_up_'+carbon+'_pos',
                  '_no_'+carbon+'_neg', '_no_'+carbon+'_pos']
    timestamps = {}
    folders    = {}

    for kk in range(len(data_types)):
        temp_time   = [0]*nr_of_repetitions
        temp_folder = [0]*nr_of_repetitions
        temp_older_than = older_than

        for ii in range(nr_of_repetitions):     
            temp_time[ii], temp_folder[ii] = toolbox.latest_data(contains = data_types[kk], older_than = temp_older_than, return_timestamp = True)
            temp_older_than = temp_time[ii]

        if data_types[kk][0:5] == '_down':    
            timestamps[data_types[kk][1:5]+data_types[kk][8:12]]  = temp_time  
            folders[data_types[kk][1:5]+data_types[kk][8:12]]     = temp_folder
        elif data_types[kk][0:3] == '_up':    
            timestamps[data_types[kk][1:3]+data_types[kk][6:10]]  = temp_time  
            folders[data_types[kk][1:3]+data_types[kk][6:10]]     = temp_folder
        elif data_types[kk][0:3] == '_no':    
            timestamps[data_types[kk][1:3]+data_types[kk][6:10]]  = temp_time  
            folders[data_types[kk][1:3]+data_types[kk][6:10]]     = temp_folder
        else:    
            timestamps[data_types[kk][3:]] = temp_time  
            folders[data_types[kk][3:]]    = temp_folder
        
    return timestamps, folders


### the fit function for the Ramsey fit
def fit_ramsey_gaussian_decay(g_tau, g_a, *arg):
    """
    fitfunction for a ramsey modulation, with gaussian decay,
        y(x) = a + A*exp(-(x/tau)**2) * mod,

        where:
        mod = sum_i(cos(2pi*f_i*x +phi) - 1)

    Initial guesses (in this order):
        g_tau : decay const
        g_A : Amplitude
        g_a : offset

        For the modulation:
        an arbitrary no of tuples, in the form
        (g_f, g_A)[i] = (frequency, Amplitude, phase)[i]
    """
    fitfunc_str = 'y(x) = a + exp(-(x/tau)**2)*('
    no_frqs = len(arg)
    if no_frqs == 0:
        print 'no modulation frqs supplied'
        return False
    
    tau = fit.Parameter(g_tau, 'tau')
    # A = fit.Parameter(g_A, 'A')
    a = fit.Parameter(g_a, 'a')
    p0 = [tau, a]

    print 'fitting with %d modulation frequencies' % no_frqs

    frqs = []
    amplitudes = []
    phases = []
    print arg
    for i, m in enumerate(arg):
        fitfunc_str += 'A%d*cos(2pi*f%d*x+phi%d)' % (i, i, i)
        frqs.append(fit.Parameter(m[0], 'f%d'%i))
        phases.append(fit.Parameter(m[2], 'phi%d'%i))
        amplitudes.append(fit.Parameter(m[1], 'A%d'%i))
        p0.append(frqs[i])
        p0.append(amplitudes[i])
        p0.append(phases[i])
    fitfunc_str += ')'

    def fitfunc(x):
        prd = exp(-(x/tau())**2)
        mod = 0
        for i in range(no_frqs):
            mod += amplitudes[i]() * (cos(2*pi*frqs[i]()*x+phases[i]()))
        return a() + prd*mod

    return p0, fitfunc, fitfunc_str

### the fit function for the initialization fit
def fit_initialization_fidelity(g_F1, g_A1, g_F2,  g_A2, g_f, g_a, g_A, g_t):
    '''Fit function for an electron Ramsey with/without C13 initialization
    g_F1 -  Fidelity of the strongly coupled C13
    g_A1 -  Hyperfine of the strongly coupled C13
    g_F2 -  Fidelity of the initialized C13
    g_A2 -  Hyperfine of the initialized C13
    g_f  -  Detuning (carrier frequency)
    g_a  -  Offset, maybe fix to 0.5
    g_A  -  Amplitude maybe fix to 0.5
    g_t  -  Decay of Ramsey
    
    The function should take the complete data set (3 measurements) at once and fit it
    '''
 
    fitfunc_str = '''a + A * exp(-(x/t)^2) * [ 
                      F1*F2         *cos[2pi*(f + (A1+A2)/2 )*x] 
                    + F1*(1-F2)     *cos[2pi*(f + (A1-A2)/2 )*x] 
                    + (1-F1)*F2     *cos[2pi*(f + (-A1+A2)/2 )*x] 
                    + (1-F1)*(1-F2) *cos[2pi*(f + (-A1-A2)/2 )*x] 
                    ]'''

    F1  = fit.Parameter(g_F1, 'F1')
    A1  = fit.Parameter(g_A1, 'A1')
    F2  = fit.Parameter(g_F2, 'F2')
    A2  = fit.Parameter(g_A2, 'A2')

    f   = fit.Parameter(g_f, 'f')
    a   = fit.Parameter(g_a, 'a')
    A   = fit.Parameter(g_A, 'A')
    t   = fit.Parameter(g_t, 't')

    p0 = [F1, A1, F2, A2, f, a, A, t]

    def fitfunc(x):
        return (a() + A()*np.exp(-(x/t())**2) * (
                          F1()*F2()         *np.cos(2*np.pi*(f() + (A1()+A2())/2)*x) +
                          F1()*(1-F2())     *np.cos(2*np.pi*(f() + (A1()-A2())/2)*x) +
                          (1-F1())*F2()     *np.cos(2*np.pi*(f() + (-A1()+A2())/2)*x) +
                          (1-F1())*(1-F2()) *np.cos(2*np.pi*(f() + (-A1()-A2())/2)*x))
                          )
    return p0, fitfunc, fitfunc_str

def fit_initialization_fidelity_all(g_F1, g_A1, g_F2,  g_A2, g_f, g_a, g_A, g_t):
    '''Fit function for an electron Ramsey with/without C13 initialization
    g_F1 -  Fidelity of the strongly coupled C13
    g_A1 -  Hyperfine of the strongly coupled C13
    g_F2 -  Fidelity of the initialized C13
    g_A2 -  Hyperfine of the initialized C13
    g_f  -  Detuning (carrier frequency)
    g_a  -  Offset, maybe fix to 0.5
    g_A  -  Amplitude maybe fix to 0.5
    g_t  -  Decay of Ramsey
    
    The function should take the complete data set (3 measurements) at once and fit it
    '''
 
    fitfunc_str = '''a + A * exp(-(x/t)^2) * [ 
                      F1*F2         *cos[2pi*(f + (A1+A2)/2 )*x] 
                    + F1*(1-F2)     *cos[2pi*(f + (A1-A2)/2 )*x] 
                    + (1-F1)*F2     *cos[2pi*(f + (-A1+A2)/2 )*x] 
                    + (1-F1)*(1-F2) *cos[2pi*(f + (-A1-A2)/2 )*x] 
                    ]'''

    F1  = fit.Parameter(g_F1, 'F1')
    A1  = fit.Parameter(g_A1, 'A1')
    F2  = fit.Parameter(g_F2, 'F2')
    A2  = fit.Parameter(g_A2, 'A2')

    f   = fit.Parameter(g_f, 'f')
    a   = fit.Parameter(g_a, 'a')
    A   = fit.Parameter(g_A, 'A')
    t   = fit.Parameter(g_t, 't')

    p0 = [F1, A1, F2, A2, f, a, A, t]

    def fitfunc(x):

        L = len(x)/3
               

        x_no_init   = x[0:L]
        x_up        = x[L:2*L]
        x_down      = x[2*L:3*L]


        S_no_init =       (a() + A()*np.exp(-(x_no_init/t())**2) * (
                          0.5*F1()          *np.cos(2*np.pi*(f() + (A1()+A2())/2)*x_no_init) +
                          0.5*F1()          *np.cos(2*np.pi*(f() + (A1()-A2())/2)*x_no_init) +
                          0.5*(1-F1())      *np.cos(2*np.pi*(f() + (-A1()+A2())/2)*x_no_init) +
                          0.5*(1-F1())      *np.cos(2*np.pi*(f() + (-A1()-A2())/2)*x_no_init))
                          )

        S_up =       (a() + A()*np.exp(-(x_up/t())**2) * (
                          F1()*F2()         *np.cos(2*np.pi*(f() + (A1()+A2())/2)*x_up) +
                          F1()*(1-F2())     *np.cos(2*np.pi*(f() + (A1()-A2())/2)*x_up) +
                          (1-F1())*F2()     *np.cos(2*np.pi*(f() + (-A1()+A2())/2)*x_up) +
                          (1-F1())*(1-F2()) *np.cos(2*np.pi*(f() + (-A1()-A2())/2)*x_up))
                          )

        S_down =       (a() + A()*np.exp(-(x_down/t())**2) * (
                          F1()*F2()         *np.cos(2*np.pi*(f() + (A1()-A2())/2)*x_down) +
                          F1()*(1-F2())     *np.cos(2*np.pi*(f() + (A1()+A2())/2)*x_down) +
                          (1-F1())*F2()     *np.cos(2*np.pi*(f() + (-A1()-A2())/2)*x_down) +
                          (1-F1())*(1-F2()) *np.cos(2*np.pi*(f() + (-A1()+A2())/2)*x_down))
                          )

        S_tot = r_[S_no_init, S_up, S_down]

        return S_tot

    return p0, fitfunc, fitfunc_str

### get the data locations
''' Data set 1: older than 20141202_092002, 10 measurements, C1, C2, C5
    Data set 2: older than 20141203_212603, 20 measurements, C2
    Second data set for C1: older than 20150103_100036, 39 msmnts '''

# ssro_calib_folder = 'D:\\measuring\\data\\20141201\\114408_AdwinSSRO_SSROCalibration_111_1_sil18'

ssro_calib_folder = 'D:\\measuring\\data\\20150108\\065633_AdwinSSRO_SSROCalibration_111_1_sil18'

carbon_nr = 'C2'

# timestamp_dict, folders_dict = get_data_timestamps('20141202_092002', 10, carbon = carbon_nr)
# timestamp_dict, folders_dict = get_data_timestamps('20141203_212603', 20, carbon = carbon_nr)
if carbon_nr == 'C1':
    timestamp_dict, folders_dict = get_data_timestamps('20150103_100036', 38, carbon = carbon_nr)

elif carbon_nr == 'C5':
    timestamp_dict, folders_dict = get_data_timestamps('20150108_084312', 29, carbon = carbon_nr)

elif  carbon_nr == 'C2':
    timestamp_dict, folders_dict = get_data_timestamps('20150109_090120', 29, carbon = carbon_nr)



#######################
### Nuclear RO data ###
#######################

### Get the data
# a_C13_RO_up, folder_up = get_and_plot_data_nuc_RO(timestamp_dict['up_pos'],timestamp_dict['up_neg'],ssro_calib_folder)
# a_C13_RO_down, folder_down = get_and_plot_data_nuc_RO(timestamp_dict['down_pos'],timestamp_dict['down_neg'],ssro_calib_folder)
# a_C13_RO_no, folder_no = get_and_plot_data_nuc_RO(timestamp_dict['no_pos'],timestamp_dict['no_neg'],ssro_calib_folder)


if 0:

    ### plot for up initialization 
    fig,ax = plt.subplots(figsize=(5,5)) 
    rects = ax.bar(a_C13_RO_up.sweep_pts,a_C13_RO_up.p0,yerr=a_C13_RO_up.u_p0,align ='center',ecolor = 'k' )
    ax.set_xticks(a_C13_RO_up.sweep_pts)
    ax.set_xticklabels(a_C13_RO_up.x)
    ax.set_ylim(-1.1,1.1)
    ax.set_title('Up initialization integrated'+'---'+ str(folder_up))
    ax.hlines([-1,0,1],-0.1,2.1,linestyles='dotted')
    for ii,rect in enumerate(rects):
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(a_C13_RO_up.p0[ii],2)) +'('+ str(int(round(a_C13_RO_up.u_p0[ii]*100))) +')',
            ha='center', va='bottom')
    fig.savefig(os.path.join(folder_up,'tomo.png'))

    ### plot for down initialization 
    fig,ax = plt.subplots(figsize=(5,5)) 
    rects = ax.bar(a_C13_RO_down.sweep_pts,a_C13_RO_down.p0,yerr=a_C13_RO_down.u_p0,align ='center',ecolor = 'k' )
    ax.set_xticks(a_C13_RO_down.sweep_pts)
    ax.set_xticklabels(a_C13_RO_down.x)
    ax.set_ylim(-1.1,1.1)
    ax.set_title('down initialization integrated'+'---'+ str(folder_down))
    ax.hlines([-1,0,1],-0.1,2.1,linestyles='dotted')
    for ii,rect in enumerate(rects):
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(a_C13_RO_down.p0[ii],2)) +'('+ str(int(round(a_C13_RO_down.u_p0[ii]*100))) +')',
            ha='center', va='bottom')
    fig.savefig(os.path.join(folder_down,'tomo.png'))

    ### plot for no initialization 
    fig,ax = plt.subplots(figsize=(5,5)) 
    rects = ax.bar(a_C13_RO_no.sweep_pts,a_C13_RO_no.p0,yerr=a_C13_RO_no.u_p0,align ='center',ecolor = 'k' )
    ax.set_xticks(a_C13_RO_no.sweep_pts)
    ax.set_xticklabels(a_C13_RO_no.x)
    ax.set_ylim(-1.1,1.1)
    ax.set_title('no initialization integrated'+'---'+ str(folder_no))
    ax.hlines([-1,0,1],-0.1,2.1,linestyles='dotted')
    for ii,rect in enumerate(rects):
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(a_C13_RO_no.p0[ii],2)) +'('+ str(int(round(a_C13_RO_no.u_p0[ii]*100))) +')',
            ha='center', va='bottom')
    fig.savefig(os.path.join(folder_no,'tomo.png'))

    ### Final initialization results
    average_initialization  = (a_C13_RO_up.p0[2] - a_C13_RO_down.p0[2])/2
    average_initialization_u= ((a_C13_RO_up.u_p0[2]**2 + a_C13_RO_down.u_p0[2]**2)**0.5)/2

    print 'Initialization by nuclear RO'
    print 'Up = '   + str(round(a_C13_RO_up.p0[2],3)) + '+/-' + str(round(a_C13_RO_up.u_p0[2],3))
    print 'Down = ' + str(round(a_C13_RO_down.p0[2],3)) + '+/-' + str(round(a_C13_RO_down.u_p0[2],3))
    print 'Average = ' + str(round(average_initialization,3)) + '+/-' + str(round(average_initialization_u,3))

#########################
### Electron Ramsey data ###
#########################
### load electron ramsey data (and plot individual plots)
# a_noinit, folder = get_and_plot_data(timestamp_dict['noInit'], ssro_calib_folder)
# a_up, folder     = get_and_plot_data(timestamp_dict['up'], ssro_calib_folder)
# a_down, folder   = get_and_plot_data(timestamp_dict['down'], ssro_calib_folder)

    ###########################
    ###  fit data seperately ###
    ###########################

if 0:
    ### Plot the data all together
    fig = a_up.default_fig(figsize=(10,5))
    ax  = a_up.default_ax(fig)
    ax.set_xlim(a_up.x[0]-1,a_up.x[-1]+1)

    ax.errorbar(a_noinit.x,    a_noinit.p0, a_noinit.u_p0,0*np.ones(len(a_noinit.u_p0)), '.b', markersize = 2, label = 'no_init') 
    ax.errorbar(a_up.x,        a_up.p0,     a_up.u_p0,    0*np.ones(len(a_up.u_p0)),     '.r', markersize = 2, label = 'up') 
    ax.errorbar(a_down.x,      a_down.p0,   a_down.u_p0,  0*np.ones(len(a_down.u_p0)),   '.k', markersize = 2, label = 'down') 
    ax.legend()

    ### Fitting the data seperately to Ramsey decays
    guess_f1    = 0.6e-3 #in GHz
    guess_A1    = 0.5
    guess_phi1  = 0.
    guess_f2    = 0.40e-3 #in GHz
    guess_A2    = 0.5
    guess_phi2  = 0.
    guess_tau   = 4600
    guess_a     = 0.5

    fit_result1 = fit.fit1d(a_noinit.x, a_noinit.p0, fit_ramsey_gaussian_decay,
            guess_tau, guess_a, (guess_f1, guess_A1, guess_phi1),
            (guess_f2, guess_A2, guess_phi2),
             fixed=[1,4,7],
            do_print=True, ret=True)
    fit_result2 = fit.fit1d(a_up.x, a_up.p0, fit_ramsey_gaussian_decay,
            guess_tau, guess_a, (guess_f1, guess_A1, guess_phi1),
            (guess_f2, guess_A2, guess_phi2),
             fixed=[1,4,7],
            do_print=True, ret=True)
    fit_result3 = fit.fit1d(a_down.x, a_down.p0, fit_ramsey_gaussian_decay,
            guess_tau, guess_a, (guess_f1, guess_A1, guess_phi1),
            (guess_f2, guess_A2, guess_phi2),
             fixed=[1,4,7],
            do_print=True, ret=True)

    plot.plot_fit1d(fit_result1, np.linspace(0,a_noinit.x[-1],201), ax=ax,
            plot_data=False, linestyle = '-b')
    plot.plot_fit1d(fit_result2, np.linspace(0,a_noinit.x[-1],201), ax=ax,
            plot_data=False, linestyle = '-r')
    plot.plot_fit1d(fit_result3, np.linspace(0,a_noinit.x[-1],201), ax=ax,
            plot_data=False, linestyle = '-k')

    plt.savefig(os.path.join(folder, 'electronramsey_analysis.pdf'),
            format='pdf')
    plt.savefig(os.path.join(folder, 'electronramsey_analysis.png'),
                format='png')

    ###########################
    ### to fit all at once: ###
    ###########################

if 0:
    fig = a_up.default_fig(figsize=(10,5))
    ax  = a_up.default_ax(fig)
    ax.set_xlim(a_up.x[0]-1,a_up.x[-1]+1)

    ax.errorbar(a_noinit.x,    a_noinit.p0, a_noinit.u_p0,0*np.ones(len(a_noinit.u_p0)), 'ob', markersize = 6, label = 'no_init') 
    ax.errorbar(a_up.x,        a_up.p0,     a_up.u_p0,    0*np.ones(len(a_up.u_p0)),     'or', markersize = 6, label = 'up') 
    ax.errorbar(a_down.x,      a_down.p0,   a_down.u_p0,  0*np.ones(len(a_down.u_p0)),   'ok', markersize = 6, label = 'down') 
    ax.legend()

    x_all       =  r_[ a_noinit.x, a_up.x, a_down.x]
    y_all       =  r_[ a_noinit.p0, a_up.p0, a_down.p0]
    y_all_error =  r_[ a_noinit.u_p0, a_up.u_p0, a_down.u_p0] 

        
    
    if carbon_nr == 'C1': ### Guess parameters, for C1 - Frequencies in kHz
        guess_F1    =  0.372863
        guess_A1    =  0.184e-3
        guess_F2    =  0.93
        guess_A2    =  -0.0364e-3     #Carbon 1:  -36.4

        guess_f     = 0.525e-3 
        guess_a     = 0.5
        guess_A     = 0.48
        guess_tau   = 5314

    elif carbon_nr == 'C5':         ### Guess parameters, for C5 - Frequencies in kHz

           
        guess_F1    =  0.372863
        guess_A1    =  0.184e-3
        guess_F2    =  0.93
        guess_A2    =  0.0242e-3     #Carbon 1:  -36.4 - C2: 20.6  -C5: 24.4

        guess_f     = 0.525e-3 
        guess_a     = 0.5
        guess_A     = 0.48
        guess_tau   = 5314


    elif carbon_nr == 'C2':         ### Guess parameters, for C2 - Frequencies in kHz


        guess_F1    =  0.372863
        guess_A1    =  0.184e-3
        guess_F2    =  0.93
        guess_A2    =  0.0176e-3     #Carbon 1:  -36.4 - C2: 20.6  -C5: 24.4

        guess_f     = 0.525e-3 
        guess_a     = 0.5
        guess_A     = 0.48
        guess_tau   = 5314

    p0, fitfunc, fitfunc_str = fit_initialization_fidelity_all(guess_F1, guess_A1, guess_F2, guess_A2,guess_f, guess_a, guess_A, guess_tau)

    fit_result1 = fit.fit1d(x_all, y_all, fit_initialization_fidelity_all,
            guess_F1, guess_A1, guess_F2, guess_A2,
            guess_f, guess_a, guess_A, guess_tau,
             fixed=[3],
            do_print=True, ret=True)

    ### show guess ###
    if 0:
        x_temp = linspace(a_noinit.x[0],a_noinit.x[-1],200)
        x_temp_list = r_[x_temp,x_temp,x_temp]
        L = len(x_temp_list)/3
        y_temp = fitfunc(x_temp_list)

        guess_curves_no_init =  y_temp[0:L]
        guess_curves_up      =  y_temp[L:2*L]
        guess_curves_down    =  y_temp[2*L:3*L]

        ax.plot(x_temp, guess_curves_no_init, 'b', lw=2)
        ax.plot(x_temp, guess_curves_up, 'r', lw=2)
        ax.plot(x_temp, guess_curves_down, 'k', lw=2)


    print fit_result1['params_dict']

    ### plot fit results
    p0_2, fitfunc_2, fitfunc_str_2 = fit_initialization_fidelity_all(
                                    fit_result1['params_dict']['F1'], 
                                    fit_result1['params_dict']['A1'], 
                                    fit_result1['params_dict']['F2'], 
                                    # 1,
                                    # fit_result1['params_dict']['A2'], 
                                    guess_A2,
                                    fit_result1['params_dict']['f'], 
                                    fit_result1['params_dict']['a'], 
                                    fit_result1['params_dict']['A'], 
                                    fit_result1['params_dict']['t'])

    x_temp = linspace(a_noinit.x[0],a_noinit.x[-1],200)
    x_temp_list = r_[x_temp,x_temp,x_temp]
    L = len(x_temp_list)/3
    y_temp = fitfunc_2(x_temp_list)

    guess_curves_no_init =  y_temp[0:L]
    guess_curves_up      =  y_temp[L:2*L]
    guess_curves_down    =  y_temp[2*L:3*L]

    ax.plot(x_temp, guess_curves_no_init, 'b', lw=2)
    ax.plot(x_temp, guess_curves_up, 'r', lw=2)
    ax.plot(x_temp, guess_curves_down, 'k', lw=2)




##################################
'''Fitting to a numerical model'''
################################## 

### THT: this seems an ideal job for using Qutip. Will do this on my laptop...
# from scipy.linalg import expm, sinm, cosm
# def fit_init_fidelity_numerical(g_F1, g_A_par_1, g_A_perp_1,  g_F2, g_A_par_2, g_A_perp_2, g_det, g_t):
#     '''Fit function for an electron Ramsey with/without C13 initialization
#     g_F1 -  Fidelity of the strongly coupled C13
#     g_A1 -  Hyperfine of the strongly coupled C13
#     g_F2 -  Fidelity of the initialized C13
#     g_A2 -  Hyperfine of the initialized C13
#     g_f  -  Detuning (carrier frequency)
#     g_a  -  Offset, maybe fix to 0.5
#     g_A  -  Amplitude maybe fix to 0.5
#     g_t  -  Decay of Ramsey
    
#     The function should take the complete data set (3 measurements) at once and fit it
#     '''

#     fitfunc_str = '''numerical solution'''

#     ### Parameters
#         ### Spin 1
#     F1          = fit.Parameter(g_F1, 'F1')
#     A_par_1     = fit.Parameter(g_A_par_1, 'A_par_1')
#     A_perp_1    = fit.Parameter(g_A_perp_1, 'A_perp_1')
        
#         ### Spin 2
#     F2          = fit.Parameter(g_F2, 'F2')
#     A_par_2     = fit.Parameter(g_A_par_2, 'A_par_2')
#     A_perp_2    = fit.Parameter(g_A_perp_2, 'A_perp_2')

#         ### Ramsey
#     det = fit.Parameter(g_det, 'det')
#     t   = fit.Parameter(g_t, 't')

#     p0 = [F1, A_par_1, A_perp_1, F2, A_par_2, A_perp_2, det, t]

#     def fitfunc(x):

#         ### pauli matrices and basis states
#         Id = np.matrix([[1, 0], [0, 1]])
#         sz = np.matrix([[1, 0], [0, -1]])
#         Iz = np.matrix([[1, 0], [0, -1]])
#         Ix = np.matrix([[0, 1], [1,  0]])

#         rho_x = (Id+Ix)/2.
#         rho_Z = (Id+sz)/2.
#         rho_mZ = (Id-sz)/2.

#         ### Initial states 
#         rho_e    = (Id+Ix)/2.
#         rho_c1   = F1*rho_Z + (1-F1)*rho_mZ 
#         rho_c2   = F2*rho_Z + (1-F2)*rho_mZ 
#         rho = np.kron(np.kron(rho_e,rho_c1),rho_c2) 

#         ### Hamiltonian
#         H  =  (det*np.kron(sz,Id,Id) +
#               A_par_1*np.kron(sz,Iz,Id) + A_par_2*np.kron(sz,Ix,Id) +
#               A_par_2*np.kron(sz,Id,Iz) + A_perp_2*np.kron(sz,Id,Ix) )

#         for tau in x:
 
#             U  =  expm(-1j*H*tau)
            
#             ### State evolution
#             rho_final = U * rho * U
#             print rho_final

#             ### Final fidelity

#     return p0, fitfunc, fitfunc_str

# data_F1          = 1
# data_A_par_1     = 1e3
# data_A_perp_1    = 0
# data_F2          = 1
# data_A_par_2     = 1e3
# data_A_perp_2    = 0
# data_det         = 400e3
# data_t           = 0

# p0, fitfunc, fitfunc_str = fit_init_fidelity_numerical(data_F1, data_A_par_1, data_A_perp_1,
#                                                        data_F2, data_A_par_2, data_A_perp_2,
#                                                        data_det, data_t)

# xdata = np.linspace(0,6e-6,2)
# ydata = fitfunc(xdata)

# plt.figure()
# plt.plot(xdata,ydata,'bo')
