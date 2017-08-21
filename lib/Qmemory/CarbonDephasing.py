import numpy as np
import os
import h5py
from analysis.lib.tools import toolbox; reload(toolbox)
from analysis.lib.m2.ssro import mbi
reload(mbi)
from matplotlib import pyplot as plt
import matplotlib as mpl
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
reload(fit)
reload(common)
import string
import analysis.lib.QEC.hyperfine_params as hf ### used for perp_coupling vs ZZ
import csv
import copy as cp
import matplotlib.cm as cm

mpl.rc('pdf', fonttype=42)
pdf_with_rc_fonts = {
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'font.sans-serif': ['Arial'],
}
mpl.rcParams.update(pdf_with_rc_fonts)

from matplotlib import rcParams
rcParams['xtick.major.size'] = 4
#xtick.minor.size
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['lines.linewidth'] = 1
rcParams['axes.linewidth'] = 0.6


VERBOSE = False

color_list = ['b','g','y','r','brown','m','c']
linewidth = 1
errorbar_width = 2
figwidthPRL=3.+3./8.
golden_ratio = 1.62
figsize=(figwidthPRL,figwidthPRL/1.4)
axeslabel_fontsize = 7
ticklabel_fontsize = 7
fignumber_fontsize=11
legend_fontsize = 7
markersize = 3
majorticklength = 3
minorticklength = 1.5
tickwidth = 0.6
axeswidth = 0.6
save_figure_to = 'D:\measuring\QMem_plots'


CR_after_check = True ### discard events with ionization for data analysis? (this relies on the CR check after the SSRO.)

def get_from_hdf5(folder,key_list):
    # gets a msmt_parameter from an hdf5 file
    Datafile=h5py.File(folder+folder[26:] + '.hdf5','r') 

    first_dict_layer = Datafile[Datafile.keys()[0]]
    return_list = []

    for key in key_list:
        return_list = return_list+[first_dict_layer.attrs[key]]
    return return_list

def get_tstamp_from_folder(folder):
    return folder[18:18+15]

def get_dephasing_data(folder_dict,ssro_calib_folder,**kw):

    tomos = kw.pop('tomos',['X','Y'])
    ## contains tomo values
    data_dict = {   
    'sweep_pts': [],
    }


    for t in ['X','Y','XX','XY','YX','YY', 'Z', 'ZZ']:
        data_dict.update({t:[]})
        data_dict.update({t+'_u':[]})
    if VERBOSE:
        print 'tomos:',  tomos
        print 'I got folder dict ', folder_dict
    for t in tomos:
        for i,f in enumerate(folder_dict[t]):
            a = mbi.MBIAnalysis(f)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata',CR_after_check = CR_after_check)
            a.get_electron_ROC(ssro_calib_folder)
            a.get_sequence_length()
            x_labels = a.sweep_pts.reshape(-1)
            if len(x_labels) != len(a.p0.reshape(-1)) :
                print 'Warning: X and y axis have different size. Check mm setting pts and sweep_pts'
            if i == 0:
                data_dict[t] = ((a.p0.reshape(-1))-0.5)*2
                data_dict[t+'_u'] = 2*a.u_p0.reshape(-1)
            else:
                y = ((a.p0.reshape(-1))-0.5)*2
                y_u = 2*a.u_p0.reshape(-1)
                data_dict[t] = [y0/2-y[ii]/2 for ii,y0 in enumerate(data_dict[t])]
                data_dict[t+'_u'] = [np.sqrt(y0**2+y_u[ii]**2)/2 for ii,y0 in enumerate(data_dict[t+'_u'])]
        #if VERBOSE:
        #    print data_dict
    ## one carbon experiment
    if len(tomos[0]) ==1:
        npY = np.array(data_dict['Y'])
        npX = np.array(data_dict['X'])
        npZ = np.array(data_dict['Z'])
        npY_u = np.array(data_dict['Y_u'])
        npX_u = np.array(data_dict['X_u'])
        npZ_u = np.array(data_dict['Z_u'])

        is_X_measurement = kw.pop('is_X_measurement', True)
        if VERBOSE:
            print 'is_X_measurement in get_dephasing_data: ', is_X_measurement
        if not is_X_measurement:
            if VERBOSE:
                print 'get_dephasing_data will return Z'
            npY = npZ
            npY_u = npZ_u
            npX = npZ
            npX_u = npZ_u
        if len(x_labels) != len(npX):
            print 'Warning'
        if kw.pop('do_get_sequence_length', False):
            if VERBOSE:
                print 'get_dephasing_data also returns the sequence length'
            return x_labels,npX,npY,npX_u,npY_u, 2*a.repump_wait + a.fast_repump_duration - a.AOM_delay + a.initial_wait - a.avg_repump_time # XXXXX check
        else:
            if VERBOSE:
                print 'get_dephasing_data does not return the sequence length'
            return x_labels,npX,npY,npX_u,npY_u

    elif len(tomos[0]) == 2:
        data_dict['sweep_pts'] = x_labels
        if kw.pop('do_get_sequence_length', False):
            return data_dict, 2*a.repump_wait + a.fast_repump_duration - a.AOM_delay + a.initial_wait - a.avg_repump_time # XXXXX check
        else:
            return data_dict

def extract_data_from_sweep(older_than = None,
        folder_name ='Repetitions_',
        carbon = '2',
        ssro_calib_timestamp = None, 
        do_T2correct=False, **kw) :

    '''
    searches for all necessary files, extracts the data, does T2* correction 
    and returns the resulting data in folder_dict
    '''

    folder_dict = {
    'sweep_pts': [],   # what has been swept
    'res' : [],         #Bloch vector length
    'res_u' : []
    }

    is_X_measurement = kw.get('is_X_measurement', True)
    logicstate = kw.get('logicstate', 'X')
    tau_larmor = kw.get('tau_larmor', None)

    if VERBOSE:
        print 'extract data: is_X_measurement ', is_X_measurement, ' , logicstate', logicstate
    
    for t in ['X','Y','XX','XY','YX','YY', 'Z', 'ZZ']:
        folder_dict.update({t:[]})

    ### search data
    if len(carbon) ==1:
        for ro in ['positive','negative']:
            if is_X_measurement:
                single_tomos = ['X','Y']
            else:
                single_tomos = ['Z']
            for t in single_tomos:
                if tau_larmor == None:
                    search_string = ro+'_Tomo_'+t+'_'+'C'+carbon
                else:
                    search_string = ro+'_Tomo_'+t+'_'+'C'+carbon+'tLarmor'+str(tau_larmor)
                if VERBOSE:
                    print 'search string is', search_string
                folder_dict[t].append(toolbox.latest_data(contains = search_string,
                    older_than=older_than, raise_exc = False, VERBOSE=False))

    ### two carbons were involved
    elif len(carbon) == 2:
        for ro in ['positive','negative']:
            if is_X_measurement:
                double_tomos =['XX','YY','XY','YX']
            else:
                double_tomos =['ZZ']
            for t in double_tomos:
                search_string = ro+'_state'+logicstate+'_Tomo_'+t+'_'+'C'+carbon
                folder_dict[t].append(toolbox.latest_data(contains = search_string,
                    older_than = older_than,raise_exc = False))

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSROCalibration', older_than = older_than)
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.data_from_time(ssro_calib_timestamp)

    if len(carbon) == 1:
        if VERBOSE:
            print 'extracting data for single carbon number', carbon
        x_labels,npX,npY,npX_u,npY_u, seq_length = get_dephasing_data(folder_dict, ssro_calib_folder, do_get_sequence_length=True, tomos=single_tomos, **kw)
        if VERBOSE:
            print 'x_labels ', x_labels, ' X_readout (npx) ', npX
        if is_X_measurement:
            folder_dict['res'] = np.sqrt(npY**2+npX**2)
            folder_dict['res_u'] = np.sqrt((npX*npX_u)**2+(npY*npY_u)**2)/np.sqrt((npX**2+npY**2))
        else: # Give Bloch vector length Z coordinates, which come from get_dephasing data as 'X'
            folder_dict['res'] = npX
            folder_dict['res_u'] = npX_u
        folder_dict['sweep_pts'] = x_labels


    elif len(carbon) == 2:
        if VERBOSE:
            print 'extracting data for carbons', carbon
        tomo_dict, seq_length = get_dephasing_data(folder_dict,ssro_calib_folder, tomos = double_tomos, do_get_sequence_length=True)
        x_labels = tomo_dict['sweep_pts']
        XX,XX_u = np.array(tomo_dict['XX']),np.array(tomo_dict['XX_u'])
        YY,YY_u = np.array(tomo_dict['YY']),np.array(tomo_dict['YY_u'])
        XY,XY_u = np.array(tomo_dict['XY']),np.array(tomo_dict['XY_u'])
        YX,YX_u = np.array(tomo_dict['YX']),np.array(tomo_dict['YX_u'])

        folder_dict['res'] = np.sqrt(XX**2+YY**2+XY**2+YX**2)/np.sqrt(2) ### geometric sum of possible correlations divided by maximum --> sqrt(2)
        ### calculate error bar in the same fashion
        folder_dict['res_u'] = (np.sqrt((XX*XX_u)**2+(YY*YY_u)**2+(YX*YX_u)**2+(XY*XY_u)**2)/np.sqrt(XX**2+YY**2+XY**2+YX**2))/np.sqrt(2)
        folder_dict['sweep_pts'] = x_labels
        if VERBOSE:
            print 'res in extract_data_from_sweep: ', folder_dict['res']

    if do_T2correct:
        if VERBOSE:
            print 'Im doing T2 correction in get_data_from_sweep '
        folder_dict = do_T2_correction(
                carbon=carbon, folder_dict=folder_dict, sequence_duration_us=seq_length*10**6)

    return folder_dict

def do_T2_correction(carbon='2', folder_dict=[], sequence_duration_us=10):
    if VERBOSE:
        print 'correcting for T2* decay before fitting. Are you sure this is what you want?'
    T2star_us = {'2': 12600, '1': 10150, '5': 4300, '3': 6350, '6': 4150 }

    n_of_reps = np.array(folder_dict['sweep_pts'])
    bloch_vec_len = np.array(folder_dict['res'])
    if VERBOSE:
        print 'folder_dict in do_T2_correction: ', folder_dict
    y_u = np.array(folder_dict['res_u'])

    extrapolatedT2star={}
    for first in T2star_us:
        for second in T2star_us:
            extrapolatedT2star[first+second]= 1./ np.sqrt( (1./T2star_us[first])**2+(1./T2star_us[second])**2)
    T2star_us.update( extrapolatedT2star)
    datapoints = len(n_of_reps)
    T2_Factors=np.zeros(datapoints)

    cut_index = len(n_of_reps)
    for count in np.arange(datapoints):
        T2_Factors[count] = np.exp(-((float(n_of_reps[count])*sequence_duration_us)**2)/(T2star_us[carbon]**2))
        if T2_Factors[count] < 0.2:
            cut_index = count
            break
            
    #### cut away unreliable datapoints due to T2*
    bloch_vec_len = bloch_vec_len[:cut_index]
    T2_Factors = T2_Factors[:cut_index]
    n_of_reps = n_of_reps[:cut_index]
    y_u = y_u[:cut_index]

    print sequence_duration_us
    if 0. in T2_Factors:
        print 'Warning: devision by zero would be required. I take uncorrected data '
        return folder_dict
    else:
        folder_dict['sweep_pts'] = n_of_reps
        folder_dict['res'] = np.divide( bloch_vec_len, T2_Factors)
        folder_dict['res_u'] = np.divide( y_u, T2_Factors) 
        return folder_dict


def create_plot(folder_name, folder_dict, carbon, fit_result):
    if folder_name == 'Memory_Sweep_repump_time_':
        plot_Xlabel = 'average repump time (us)'
    elif folder_name == 'Repetitions_':
        plot_Xlabel = 'Number of repetitions'
    elif folder_name == 'Memory_Sweep_repump_duration':
        plot_Xlabel = 'Repump duration'
    elif folder_name == 'Memory_sweep_timing_':
        plot_Xlabel = 'Waiting time (us)'
    else:
        plot_Xlabel = folder_name

    fig = plt.figure()
    ax = plt.subplot()
    plt.errorbar(folder_dict['sweep_pts'],folder_dict['res'], folder_dict['res_u'], marker='.', label='C'+carbon)
    plot.plot_fit1d(fit_result, np.linspace(folder_dict['sweep_pts'][0],folder_dict['sweep_pts'][-1],1001),ax=ax, plot_data=False,add_txt=False, lw = 2)

    plt.xlabel(plot_Xlabel)
    plt.ylabel('Bloch vector length')
    plt.title('Dephasing for C'+carbon+' '+ get_tstamp_from_folder(folder_dict['X'][0]))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=legend_fontsize)
    plt.savefig(os.path.join(folder_dict['X'][0],'CarbonDephasing.pdf'),format='pdf')
    plt.savefig(os.path.join(folder_dict['X'][0],'CarbonDephasing.png'),format='png')
    plt.show()
    plt.close('all')

def extract_coupling_strength(folder_dict):
     ### after fitting: get the coupling strength. Here we need logic state ms = 0 and ms = +-1 frequencies.
        coupling = 0
        #extract one complete path
        for key in ['XX','X', 'Z', 'ZZ']:
            if folder_dict[key] != []:
                folder = folder_dict[key][0]

        carbon_list,logic_state = get_from_hdf5(folder,['carbon_list','2qb_logical_state'])
        ### extract the expected coupling strength (works for 1 and 2 carbons)
        for ii, c in enumerate(carbon_list):
            if ii == 0:
                [coupling] = get_from_hdf5(folder,['C'+str(c)+'_freq_1_m1'])
                f_m1,f_0 = get_from_hdf5(folder,['C'+str(c)+'_freq_1_m1','C'+str(c)+'_freq_0'])
            else:
                f_m1,f_0 = get_from_hdf5(folder,['C'+str(c)+'_freq_1_m1','C'+str(c)+'_freq_0'])
                if logic_state == 'X':
                    coupling += f_m1-f_0
                    if VERBOSE:
                        print 'X ', coupling
                elif logic_state == 'mX':
                    coupling += - f_0 - f_m1
                    if VERBOSE:
                        print 'mX ', coupling
                else:
                    print 'Warning: Logical state not properly read out'
                    coupling = 0
        
        #calculate Delta_f w.r.t. to f_0
        if ii ==0:
            coupling = np.abs(f_0-f_m1)
        else:
            coupling = np.abs(np.abs(coupling)-f_0)       

        return folder_dict, coupling, folder

def Sweep_repetitions(older_than = None,
        folder_name ='Repetitions_',
        carbon = '1',
        ssro_calib_timestamp =None, 
        plot_result= True,
        fit_result = True,
        do_T2correct=False, **kw) :

    folder_dict =  extract_data_from_sweep(older_than = older_than, folder_name =folder_name, carbon = carbon,
        ssro_calib_timestamp =ssro_calib_timestamp, do_T2correct=do_T2correct, **kw)

    if folder_name == 'Memory_sweep_timing_':
        folder_dict['sweep_pts'] = folder_dict['sweep_pts']*1e6
    
    ### fit a function to the extracted data. return the fit results and coupling strength (extract from hdf5 file).
    A0 = np.max(folder_dict['res'])
    offset = 0
    decay, x0 = 50, 0

    fitGauss = kw.get('fitGauss', False)
    if fitGauss:
        offset = folder_dict['sweep_pts'][np.argmin(folder_dict['res'])]
        x0 = folder_dict['sweep_pts'][np.argmax(folder_dict['res'])]
        p0,fitfunc,fitfunc_str = common.fit_gauss(offset, A0, x0, 1)
        fixed = []
    else:
        p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)
        fixed = [0,3]

    fit_result = fit.fit1d(folder_dict['sweep_pts'],folder_dict['res'],None,
            p0 = p0, fitfunc = fitfunc, do_print = False, ret = True, fixed = fixed)

    folder_dict, coupling, folder = extract_coupling_strength(folder_dict)

    if plot_result:
        create_plot(folder_name, folder_dict, carbon, fit_result)

    if not fitGauss:
        return coupling, fit_result['params_dict']['tau'], fit_result['error_dict']['tau'], folder, fit_result
    else:
        print fit_result['params_dict'], fit_result['error_dict']

def Sweep_Rep_List( carbons = ['1'],
        older_than = None, 
        do_T2correct = False,
        folder_name = 'Repetitions_', 
        ssro_calib_timestamp = None,**kw):

    fit_results = kw.get('fit_results',True)
    sequence_length = kw.get('sequence_length',None)
    logicstate_list = kw.get('logicstate_list',len(carbons)*['X']) ## can be list such as ['X','mX'], used for DFS measurements.
    colors = kw.get('colors', color_list)
    log_plot = kw.get('log_plot', False)
    fitGauss = kw.get('fitGauss', False)
    return_fits = kw.get('return_fits', False)

    do_plot_results = kw.pop('do_plot_results', True)

    x_arr = []
    y_arr = []
    y_u_arr = []
    for c,logicstate in zip(carbons,logicstate_list):
        folder_dict= extract_data_from_sweep(
            older_than = older_than,
            folder_name =folder_name, carbon = c,
            ssro_calib_timestamp =ssro_calib_timestamp,
            logicstate = logicstate,
            do_T2correct=do_T2correct, **kw)
        if VERBOSE:
            print 'folder_dict in Sweep_Rep_List: ', folder_dict
        x_arr.append(folder_dict['sweep_pts'])
        y_arr.append(folder_dict['res'])
        y_u_arr.append(folder_dict['res_u'])

    ### convert to time instead of repetitions:
    if sequence_length != None:
        x_arr = [x*sequence_length for x in x_arr]
    
    if do_plot_results:
        fig = plt.figure()
        ax = plt.subplot()
    is_X_measurement = kw.get('is_X_measurement', False)
    print is_X_measurement
    if is_X_measurement:
        for key in ['XX','X']:
            if folder_dict[key] != []:
                folder = folder_dict[key][0]
    else:
        for key in ['ZZ','Z']:
            if folder_dict[key] != []:
                folder = folder_dict[key][0]

    for x,y,y_u,carbon,logicstate,jj in zip(x_arr,y_arr,y_u_arr,carbons,logicstate_list,range(len(x_arr))):

        if fit_results:
            if sequence_length != None:
                decay = decay*sequence_length
            if fitGauss:
                
                offset = min(y)/4
                A0 = max(y)/4
                x0 = x[np.argmax(y)]
                decay = np.abs( x[np.argmax(y)]-x[np.argmin(y)] ) /3
                if VERBOSE:
                    print 'Fit guess A0', A0, 'offset', offset, 'x0', x0, 'decay', decay
                p0,fitfunc,fitfunc_str = common.fit_gauss(offset,A0,x0,decay)
                fixed=[]
            else:
                A0 = y[0]
                x0 = 0
                offset = 0
                decay = 500
                p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)
                fixed = [0,3]

            fit_result = fit.fit1d(x,y,None,p0 = p0, fitfunc = fitfunc, do_print = not return_fits, ret = True, fixed = fixed)
            
            if do_plot_results:
                plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax,color = colors[jj], plot_data=False,add_txt = False, lw = 2)

        label_txt = 'C'+carbon
        if len(carbon)!=1:
            label_txt = label_txt+'_'+logicstate
        
        if do_plot_results:
            plt.errorbar(x,y,y_u,marker='.',color = colors[jj],label=label_txt)

    if do_plot_results:
        plt.xlabel('Number of repetitions')
        if sequence_length != None:
            plt.xlabel('elapsed time (us)')
        elif folder_name == 'Memory_Sweep_repump_time_':
            plt.xlabel('average repump time (us)')

        plt.ylabel('Bloch vector length')

        if log_plot:
            ax.set_yscale("log", nonposy='clip')
        y_min=kw.get('ymin', 0.8 * min(y))    
        plt.ylim(y_min,1.0)

        plt.title(get_tstamp_from_folder(folder) + ' Dephasing for C'+carbon)
        plt.legend()#bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(os.path.join(folder,'CarbonDephasing.pdf'),format='pdf')
        plt.savefig(os.path.join(folder,'CarbonDephasing.png'),format='png')
        #if not return_fits:
        plt.show()
        plt.close('all')

    if return_fits:
        return folder, fit_result

def sweep_avg_repump_and_tau_larmor( carbons = ['1','2'],
        older_than = None, 
        do_T2correct = False,
        folder_name = 'Repetitions_', 
        ssro_calib_timestamp = None,
        tau_larmor_list = [],
        **kw):
    is_X_measurement = kw.get('is_X_measurement', True)
    A_list = [[] for i in enumerate(carbons)]

    A_u_list    = cp.deepcopy(A_list)
    x0_list     = cp.deepcopy(A_list)
    x0_u_list   = cp.deepcopy(A_list)


    for ii,c in enumerate(carbons):
        for tau_larmor in tau_larmor_list:
            folder, fit_result = Sweep_Rep_List( carbons = [c],
                older_than = older_than, 
                do_T2correct = do_T2correct,
                folder_name = 'Repetitions_', 
                ssro_calib_timestamp = None,
                tau_larmor = tau_larmor,
                return_fits = True,
                fit_result = True,
                fitGauss = True,
                **kw
                )

            A_list[ii].append(fit_result['params_dict']['A'])
            A_u_list[ii].append(fit_result['error_dict']['A'])
            x0_list[ii].append(fit_result['params_dict']['x0'])
            x0_u_list[ii].append(fit_result['error_dict']['x0'])
    # print fit_result

    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    plot_title = 'X Dephasing ' if is_X_measurement else 'Z decay '
    
    for ii, c in enumerate(carbons):
        ax.errorbar(tau_larmor_list,A_list[ii],A_u_list[ii], marker='.',label='C'+str(c))
        ax2.errorbar(tau_larmor_list,x0_list[ii],x0_u_list[ii], marker='.',label='C'+str(c))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_xlabel('t (us)')
    ax2.set_xlabel('t (us)')
    ax.set_ylabel('fitted Amplitude')
    ax.set_ylim([0.3,0.9])
    ax2.set_xlim([tau_larmor_list[0]-0.1,tau_larmor_list[-1]+0.1])
    ax.set_xlim([tau_larmor_list[0]-0.1,tau_larmor_list[-1]+0.1])
    ax2.set_ylabel('fitted tau (us)')
    ax.set_title(plot_title + get_tstamp_from_folder(folder))
    plt.show()
    plt.close('all')

def coupling_vs_repetitions(c_identifiers,**kw):

    is_X_measurement = kw.get('is_X_measurement', True)
    folder = None

    older_than = kw.get('older_than',None)

    x, y, y_u = [],[],[]
    ## acquire data
    threshold = 0.
    for carbon in c_identifiers:
        if len(carbon) > 1:
            for logicstate in ['X','mX']:
                x_temp,y_temp,y_u_temp,folder,fit_result = Sweep_repetitions(carbon = carbon,logicstate = logicstate,return_fits=True,plot_result = False,**kw)
                if fit_result['params_dict']['A'] > threshold:
                    x.append(x_temp)
                    y.append(y_temp)
                    y_u.append(y_u_temp)
                else:
                    print 'Low initial amplitude ', folder, 'A0 is ', fit_result['params_dict']['A']
        else:
            x_temp,y_temp,y_u_temp,folder,fit_result = Sweep_repetitions(carbon = carbon,return_fits=True,plot_result = False,**kw)
            x.append(x_temp)
            y.append(y_temp)
            y_u.append(y_u_temp)
            if (fit_result['reduced_chisq'] > 0.01 or y_u_temp/y_temp > 0.3):
                print 'fit failed: Chisq is', fit_result['reduced_chisq'], '. One should discard the data point ', carbon, '. y_u is ', y_u_temp, ' y is ', y_temp

    if VERBOSE:
        print 'folder in coupling vs rep: ', folder

    return 1e-3*np.array(np.sort(x)), np.array(y)[np.argsort(np.array(x))], np.array(y_u)[np.argsort(np.array(x))], folder

def repump_power_vs_repetitions(c_identifier, repump_powers=[0], **kw):

    older_than = kw.get('older_than',None)
    p = len(repump_powers)
    x = np.zeros(p)
    y = np.zeros(p)
    y_u = np.zeros(p)

    ## acquire data
    ii = 0
    for p in repump_powers:
        #print 'identifier: ', c_identifier

        ## necessary to find the correct power. Is used as an intermediate older_than timestamp
        tst = get_tstamp_from_folder(toolbox.latest_data(contains = str(p), older_than = older_than))
        ### compiles a timestamp that is 1 second older than the previously found timestamp
        folder_bookmark = str(int(tst[:8]))+str(int(tst[-6:])+1)

        if len(c_identifier) > 1:
            for logicstate in ['X','mX']:
                x[ii],y[ii],y_u[ii],folder,fit_result = Sweep_repetitions(carbon = c_identifier,
                    logicstate = logicstate,return_fits=True, plot_result = False ,**kw)
                ii +=1
        else:
            x[ii],y[ii],y_u[ii],folder,fit_result = Sweep_repetitions(carbon = c_identifier,
                older_than = folder_bookmark,return_fits=True, plot_result = False ,**kw)
            ii +=1

        if fit_result['reduced_chisq']>0.005:
                    print 'Bad fit in ', folder, 'Xi_square is ', fit_result['reduced_chisq']

    if VERBOSE:
        print 'folder ', folder
    return x,y,y_u,folder
    

def Osci_period(carbon = '1',older_than = None,ssro_calib_timestamp = None, do_print=False, add_txt=True, **kw):

    fit_results = kw.get('fit_results',True)
    folder_name = kw.get('folder_name','Memory_NoOf_Repetitions_')

    ### fit parameters
    freq = kw.pop('freq',1/170.)
    offset = kw.pop('offfset',0.)
    decay = kw.get('decay',200)
    fixed = kw.pop('fixed',[1])
    show_guess = kw.pop('show_guess',False)

    auto_analysis = kw.get('auto_analysis',False)

    if auto_analysis:
        do_print  = False
        add_txt = True


    folder_dict = {
    'X' : [],
    'Y' : [],
    'resX' : [],
    'resX_u' : [],
    'resY' : [],
    'resY_u': [],
    'sweep_pts': [],
    'res' : [],
    'res_u' : []
    }

    ### search data
    if auto_analysis:
        tomos = ['X']
    else:
        tomos = ['X','Y']
    for ro in ['positive','negative']:
        for t in tomos:
            search_string = folder_name+ro+'_Tomo_'+t+'_'+'C'+carbon
            folder_dict[t].append(toolbox.latest_data(contains = search_string,older_than = older_than,raise_exc = False))
    print folder_dict[t]

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSROCalibration', older_than = older_than)
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.data_from_time(ssro_calib_timestamp)
        if VERBOSE:
            print 'ssro Folder: ', ssro_calib_folder

        ### extract data
    x_labels,npX,npY,npX_u,npY_u = get_dephasing_data(folder_dict,ssro_calib_folder,tomos=tomos)

    fig = plt.figure()
    ax = plt.subplot()
    if auto_analysis:
        res,res_u,rng = [npX],[npX_u],range(1)
    else:
        res,res_u,rng = [npX,npY],[npX_u,npY_u],range(2)
    for y,y_u,jj in zip(res,res_u,rng):
        if fit_results:
            A0 = max(y)
            phi0 = 0
            p0,fitfunc,fitfunc_str = common.fit_decaying_cos(freq,offset,A0,phi0,decay)

            fit_result = fit.fit1d(x_labels,y,None,p0 = p0, fitfunc = fitfunc, do_print = do_print, ret = True, fixed = fixed)
            plot.plot_fit1d(fit_result, np.linspace(x_labels[0],x_labels[-1],1001), ax=ax,color = color_list[jj], plot_data=False,add_txt = add_txt, lw = 2)

            if show_guess:
                print decay
                ax.plot(np.linspace(x_labels[0],x_labels[-1],201), fitfunc(np.linspace(x_labels[0],x_labels[-1],201)), ':', lw=2)
        
        plt.errorbar(x_labels,y,y_u,marker='.',color = color_list[jj],label='C'+carbon+['X','Y'][jj])

    ## define folder for data saving
    folder = folder_dict[t][0]

    plt.xlabel('Repump repetitions')
    plt.ylabel('Contrast')
    plt.title(get_tstamp_from_folder(folder) + ' Dephasing for C'+carbon)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(os.path.join(folder,'CarbonDephasing_osci.pdf'),format='pdf')
    plt.savefig(os.path.join(folder,'CarbonDephasing_osci.png'),format='png')
    plt.show()
    # if not auto_analysis:
    #     plt.show()
    plt.close('all')

    print 'Results are saved in ', folder[18:18+15]
    if auto_analysis:
        return fit_result

# all these functions are a mess and a lot is hardcoded, gonna write my own stuff. SK.
def fit_sin_pos_neg_data(folder_name = '', measurement_name = ['adwindata'],
            offset=[0], amplitude = [0.5], center = [0], decay_constant = [200], exp_power = [0],
            frequency = [1], phase =[0],
            fixed = [], ylim = [-0.5, 1.05],
            plot_fit = False, do_print = False, show_guess = True, **kw):
    ''' Function to fit positive and negative mbi type data 
    with exponential and sinusoidal functions or combinations thereof. ~ SK
    folder_name: ro + foldername is what it searches for
    measurement_name: list of measurement names

    '''
    ax = kw.pop('ax',None)
    label = kw.pop('label','')
    older_than = kw.pop('older_than',None)

    fit_results = []
    fig = plt.figure()
    ax = plt.subplot()
    for k in range(0,len(measurement_name)):

        x,y,y_u,f = get_PosNeg_data(folder_name, older_than = older_than, **kw)
        print f
        plt.errorbar(x,y,y_u,label = label,fmt='o',**kw)
        

        if ylim != None:
            ax.set_ylim(ylim[0],ylim[1])

        ### fit depending on the number of the frequencies
        if len(frequency) == 1:
            print exp_power[0]
            p0, fitfunc, fitfunc_str = common.fit_exp_cos(offset[0],
                    amplitude[0], center[0], decay_constant[0], exp_power[0],
                    frequency[0], phase[0])
            if show_guess:
                plt.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(0,x[-1],201)), lw=2)
            print 'starting fit.fit1d'
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
        elif len(frequency) == 2:
            p0, fitfunc, fitfunc_str = common.fit_gaussian_decaying_2cos(offset[0],amplitude[0],decay_constant[0],amplitude[0],
                frequency[0],  phase[0], amplitude[1], frequency[1],  phase[1])
            if show_guess:
                plt.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

        ## plot fit
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax, plot_data=False)
        
        print 'folder combo: ' + str(f)
        print fit_result['params_dict']
        # print fit_results

        # save as pdf and png
        plt.savefig(os.path.join(f, 'analyzed_result.pdf'),         format='pdf')
        plt.savefig(os.path.join(f, 'analyzed_result.png'),    format='png')

        plt.show()
        
    return fit_result


def attempt_decay(folder_name = '', 
            tomo_basis = ['X','Y','Z'],
            offset=[0], amplitude = [0.5], center = [0], decay_constant = [200], exp_power = [0],frequency = [1], phase =[0],
            fixed = [], ylim = [-0.5, 1.05],
            plot_fit = False, do_print = False, show_guess = True, **kw):
    ''' Function to fit positive and negative mbi type data 
    with exponential and sinusoidal functions or combinations thereof. ~ SK
    folder_name: ro + foldername is what it searches for
    measurement_name: list of measurement names

    '''
    label = kw.pop('label','')
    older_than = kw.pop('older_than',None)

    fit_results = []
    fig = plt.figure()
    ax = plt.subplot()


    if len(amplitude) == 1:
        amplitude =  amplitude*len(tomo_basis)


    
    for ii, t in enumerate(tomo_basis):
        ### supply it data to fit. Can be lists of data with list of amplitudes to do multi plot
        x,y,y_u,f = get_PosNeg_data('_Tomo_'+ t + folder_name, older_than = older_than, **kw)


        plt.errorbar(x,y,y_u,label = label,fmt='.',**kw)

        #### Fitting
        ## plot guess?
        if show_guess:
            plt.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(0,x[-1],201)), lw=2)
        print 'starting fit.fit1d'

        ## do fit?
        if plot_fit:
            p0, fitfunc, fitfunc_str = common.fit_exp_cos(offset[0],
                amplitude[ii], center[0], decay_constant[0], exp_power[0],
                frequency[0], phase[0])
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
            # print fit_result

            ## plot fit
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax, plot_data=False)
            # print 'folder combo: ' + str(f)
            # print fit_result['params_dict']
    
    # figure properties
    if ylim != None:
        ax.set_ylim(ylim[0],ylim[1])



    #save
    plt.savefig(os.path.join(f, 'analyzed_result.pdf'),         format='pdf')
    plt.savefig(os.path.join(f, 'analyzed_result.png'),    format='png')

    plt.title(folder_name)
    plt.show()


### BULK ANALYSIS FUCNTIONS SK 1-5
def attempt_decay_all_data(carbon = 1, el_bases =['X','mX','Y','mY','Z','mZ'],
            tomo_bases = ['X','Y','Z'],**kw):
    ''' Function to fit positive and negative mbi type data 
    with exponential and sinusoidal functions or combinations thereof. ~ SK
    folder_name: ro + foldername is what it searches for
    measurement_name: list of measurement names

    '''
    older_than = kw.pop('older_than',None)
    x_list = []
    y_list = []
    y_u_list = []
    f_list = []
    labels = []

    # brute way to get 6x3 list of lists, would a dictionary be more useful?
    for e in el_bases:
        # x_temp = []
        # y_temp = []
        # y_u_temp =[]
        # f_temp = []
        for t in tomo_bases:
 
            ### supply it data to fit. Can be lists of data with list of amplitudes to do multi plot
            x,y,y_u,f = get_PosNeg_data('_Tomo_'+ t + '_elState_' + e + '_C' + str(carbon), 
                older_than = older_than, **kw)
            # x_temp.append(x)
            # y_temp.append(y)
            # y_u_temp.append(y_u)
            # f_temp.append(f)
            labels.append('e in '+ e +', tomo in ' + t)
            x_list.append(x)
            y_list.append(y)
            y_u_list.append(y_u)
            f_list.append(f)


    return x_list, y_list, y_u_list,f_list, labels
        

def errorplot_data(x_input,y_input,y_u_input,ax,**kw):
    
    labels = kw.pop('labels',['']*len(x_input))
    y_lim = kw.pop('y_lim',None)
    colors = cm.rainbow(np.linspace(0, 1, len(x_input)))
    for ii,(x,y,y_u,l) in enumerate(zip(x_input,y_input,y_u_input,labels)):
        ax.errorbar(x,y,yerr=y_u,fmt='.',label = l,color = colors[ii])
    

    #plot options
    if labels != ['']*len(x_input):
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if y_lim != None:
        ax.set_ylim(y_lim[0],y_lim[1])

# Not doing anything at the moment 29-04 SK
def fit_and_plot_exp_sin(
    x,
    y,
    ax,
    offset=[0], 
    amplitudes = [0.5], 
    center = [0], 
    decay_constant = [200], 
    exp_power = [0],
    frequency = [1], 
    phase =[0],
    fixed = [], 
    plot_fit = True, 
    do_print = False, 
    show_guess = False, 
    **kw):
    
    labels = kw.pop('labels',[None]*len(x))
    #too lazy for multiple amplitudes?
    if len(amplitudes)==1:
        amplitudes = amplitudes*len(x)
    colors = cm.rainbow(np.linspace(0, 1, len(x)))
    ### fit depending on the number of the frequencies
    for ii, (a,x,y) in enumerate(zip(amplitudes,x,y)):
        # print offset[0], a, center[0], decay_constant[0], exp_power[0], frequency[0], phase[0]
        p0, fitfunc, fitfunc_str = common.fit_exp_cos(offset[0],
            a, center[0], decay_constant[0], exp_power[0],
            frequency[0], phase[0])
        if show_guess:
            plt.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(0,x[-1],201)), lw=2)



        if plot_fit:
            fit_result = fit.fit1d(x,y, None, 
                p0=p0, fitfunc=fitfunc, 
                do_print=do_print, 
                ret=True,
                fixed=fixed,
                VERBOSE=False)
            # print fit_result
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), 
                ax=ax, 
                plot_data=False,
                print_info = True,
                add_txt = False,
                color = colors[ii],
                label = labels[ii] + ', Decay constant = ' + str(fit_result['params_dict']['T']))

        if labels != [None]*len(x):
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
### END BULK ANALYSIS FUNCTIONS


# Nice function!
def get_PosNeg_data(name,**kw):

    """
    The idea is that we have two folders with the same ending: name
    searches for one folder positive_+name and for negative_+name 
    averages the data and returns the sweep points, measured contrast, uncertainty and the positive folder.
    """
    ssro_calib_timestamp = kw.pop('ssro_calib_timestamp',None)
    older_than = kw.pop('older_than',None)

    data_dict = {
    'folders' : [],
    'sweep_pts': [],
    'res' : [],
    'res_u' : []
    }


    for ro in ['positive','negative']:
        search_string = ro+name
        if VERBOSE:
            print 'looking for ', search_string, '  older than is ', older_than
        data_dict['folders'].append(toolbox.latest_data(contains = search_string, older_than = older_than,raise_exc = True))


    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSROCalibration', older_than = older_than)
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.data_from_time(ssro_calib_timestamp)
    if VERBOSE:
        print 'folders: ',data_dict['folders']
    for i,f in enumerate(data_dict['folders']):
        if VERBOSE:
            print 'doing mbi analysis for', f
        a = mbi.MBIAnalysis(f)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata',CR_after_check = CR_after_check)
        a.get_electron_ROC(ssro_calib_folder)



        x_labels = a.sweep_pts.reshape(-1)
        if i == 0:
            # data_dict['res'] = abs(a.p0.reshape(-1)) # fidelity
            data_dict['res'] = ((a.p0.reshape(-1))-0.5)*2
            # data_dict['res_u'] = a.u_p0.reshape(-1) # fidelity
            data_dict['res_u'] = 2*a.u_p0.reshape(-1)
        else:
            y = ((a.p0.reshape(-1))-0.5)*2 # Contrast
            # y = abs(a.p0.reshape(-1)) # fidelity
            y_u = 2*a.u_p0.reshape(-1) # contrast
            # y_u = a.u_p0.reshape(-1) # fidelity
            data_dict['res'] = [y0/2-y[ii]/2 for ii,y0 in enumerate(data_dict['res'])]
            data_dict['res_u'] = [np.sqrt(y0**2+y_u[ii]**2)/2 for ii,y0 in enumerate(data_dict['res_u'])]

    return x_labels,data_dict['res'],data_dict['res_u'],data_dict['folders'][1]

def plot_data(x,y,**kw):
    label = kw.get('label',None)
    y_u = kw.pop('y_u',None)
    if y_u != None:
        plt.errorbar(x,y,y_u,label = label,**kw)
    else: plt.plot(x,y)

def fit_exp_pos_neg_data(folder_name,**kw):
    ax = kw.pop('ax',None)
    label = kw.pop('label','')
    older_than = kw.pop('older_than',None)

    x,y1,y_u = np.array([]),np.array([]),np.array([])

    loop_bit = True
    i = 0
    while loop_bit:

        x2,y2,y_u2,f = get_PosNeg_data(folder_name,label = label,older_than = older_than, **kw)
        print '*'*10
        print f
        # new reference
        older_than = get_tstamp_from_folder(f)

        ## add data to arrays
        if len(x) == 0:
            x = np.array(x2)
            y1 = np.array(y2)
            y_u = np.array(y_u2)
        else:
            x = np.append(x,x2); y1 = np.append(y1,y2); y_u = np.append(y_u,y_u2)

        ### if there is a single repetition then we stop looping
        i +=1
        if VERBOSE:
            print 'x2 ', x2
        if 1 in x2 or 0 in x2 or 100 in x2:
            if VERBOSE:
                print 'x2 ', x2, ' f ', f
            loop_bit = False

    ### sort arrays after adding them up
    y1 = y1[np.argsort(x)]
    y_u = y_u[np.argsort(x)]
    x = np.sort(x)


    plot_data(x,y1,y_u=y_u,**kw)

    offset,A0,decay,x0 = 0,np.amax(y1),400,0
    p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)
    fixed = [0,3]
    fit_result = fit.fit1d(x,y1,None,p0 = p0, fitfunc = fitfunc, do_print = True, ret = True, fixed = fixed)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001),ax=ax, plot_data=False,add_txt=False, lw = 2,**kw)
    return f

def Z_decay_vs_coupling(c_idents, perpendicular = True, **kw):

    ### perpendicular hyperfine copulings in kHz (from fingerprint fits).
    hyperfine_params = {}
    hyperfine_params['C1']  = 25.0
    hyperfine_params['C2']  = 43.0
    hyperfine_params['C3']  = 55.0
    hyperfine_params['C4']  = 21.0
    hyperfine_params['C5']  = 26.0
    hyperfine_params['C6']  = 12.0

    folder = kw.pop('folder',None)
    label = kw.pop('label','')
    older_than0 = kw.pop('older_than', None)
    if VERBOSE:
        print 'older than ', older_than0
    return_vals = kw.pop('return_vals',False)

    xarray,zarray,zarray_u, f = [],[],[], []
    coupling_minus = np.array([])
    coupling_plus = np.array([])
    coupling_parallel_kHz = np.array([])
    fitted_decays = np.array([])
    fitted_decays_u = np.array([])

    ### coolect data from different measurements and fuse them. (needed for ZZ decays)
    for Cs in c_idents:
        print 'Im investigating carbon(s) ', Cs
        older_than = older_than0
        for state in ['mX','X']: ### need this loop for DFS configurations
            if len(Cs) == 2:
                folder_name = '_state'+state+'_Tomo_ZZ_C'+Cs
            else:
                folder_name = '_Tomo_Z_C'+Cs

            ### reset values before new round
            loop_bit = True
            
            x, z, z_u = np.array([]),np.array([]),np.array([])
            loop_counter = 0

            while loop_bit:
                if VERBOSE:
                    print 'analysing folder: ', folder_name, 'older than is: ', older_than

                xarray, zarray, zarray_u, f = get_PosNeg_data(folder_name,label = label, older_than = older_than, **kw)
                if VERBOSE:
                    print 'x and z array:', xarray, zarray, 'folder f', f, 'folder_name is', folder_name
                folder_dict = extract_data_from_sweep(older_than = older_than,
                    folder_name = folder_name, carbon = Cs, logicstate = state,
                    do_T2correct = False, **kw)
                if VERBOSE:
                    print 'folder dict is ', folder_dict
                    print 'folder is now: ', f, ', older than is ', older_than
                older_than = get_tstamp_from_folder(f) # new reference

                ## add data to arrays
                if len(x) == 0:
                    x = np.array(xarray)
                    z = np.array(zarray)
                    z_u = np.array(zarray_u)
                else:
                    x = np.append(x, xarray)
                    z = np.append(z, zarray)
                    z_u = np.append(z_u, zarray_u)

                ### if there is a single repetition then we stop looping
                if VERBOSE:
                    print 'x array is ', x
                if 1 in x or 0 in x : #or loop_counter>0: #at most two files

                    loop_bit = False

                    if len(Cs) == 2:
                        [C1_perp, C2_perp] = [hyperfine_params['C'+Cs[0]],hyperfine_params['C'+Cs[1]]]
                        cplus = np.sqrt(C1_perp**2 + C2_perp**2)/np.sqrt(2.)
                        cminus = np.sqrt(abs(C1_perp**2 - C2_perp**2))/np.sqrt(2.)
                        coupling_minus = np.append(coupling_minus,[cminus])
                        coupling_plus = np.append(coupling_plus,[cplus])
                        folder_dict, coupling, folder = extract_coupling_strength(folder_dict)
                    else:
                        C1_perp = hyperfine_params['C'+Cs[0]]
                        coupling_minus = np.append(coupling_minus,[C1_perp])
                        coupling_plus = np.append(coupling_plus,[C1_perp])
                        folder_dict, coupling, folder = extract_coupling_strength(folder_dict)
                loop_counter +=1


            ### sort arrays after adding them up
            z = z[np.argsort(x)]
            z_u = z_u[np.argsort(x)]
            x = np.sort(x)
            ### fit exponential to the data
            if False:
                fig = plt.figure()
                ax = plt.subplot()
                plt.errorbar(x, z, yerr = z_u,fmt = '.',color = 'b', label = str(Cs))
                plt.legend()
                plt.close('all')
            offset, A0, decay, x0 = 0, np.amax(zarray), 400, 0
            p0, fitfunc, fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)
            fixed = [0,3]
            fit_result = fit.fit1d(x,z,None,p0 = p0, fitfunc = fitfunc, do_print = False, ret = True, fixed = fixed)
            if VERBOSE:
                print 'fit_result', fit_result
            fitted_decays = np.append(fitted_decays,[fit_result['params_dict']['tau']])
            fitted_decays_u = np.append(fitted_decays_u,[fit_result['error_dict']['tau']])
            coupling_parallel_kHz = np.append(coupling_parallel_kHz, coupling / 1000)

            # print fitted_decays
            ### break the state loop for a single carbon in order to not take the same dat twice.
            if len(Cs) == 1:
                break
    ### plot results:
    if return_vals:
        if VERBOSE:
            print 'returning ', coupling_parallel_kHz, fitted_decays,fitted_decays_u, f
        return coupling_parallel_kHz, fitted_decays,fitted_decays_u, f

    fig = plt.figure()
    ax = plt.subplot()

    plt.xlabel('Perpendicular hyperfine (kHz)')
    plt.ylabel('Fitted Z decay constant')

    if True: #Log-Plot
        plt.ylim(50,10000)
        ax.set_yscale("log", nonposy='clip')
    else:
        plt.ylim(0,4000)
    plt.errorbar(coupling_plus, fitted_decays, yerr = fitted_decays_u,fmt = '.',color = 'b', label = 'sum')
    plt.errorbar(coupling_minus, fitted_decays, yerr = fitted_decays_u,fmt = '.',color = 'g', label = 'difference')
    plt.legend(numpoints=1)
    if VERBOSE:
        print 'saving to ', f
    plt.savefig(os.path.join(f,'perp_coupling_vs_repetitions.pdf'),format='pdf')
    plt.savefig(os.path.join(f,'perp_coupling_vs_repeititons.png'),format='png')
    plt.show()
    plt.close('all')

def repump_speed(timestamp=None, ssro_calib_timestamp =None, older_than=None, powers = [0],
            exclude_first_n_points = 0., log_plot = False,
            amplitude =0.2, decay_constant_one = 20., decay_constant_two = 40.,
            x_offs = 0, offset=0.01, fixed = [0], 
            do_plot = True, do_fit = False, print_fit = True, 
            plot_fit=False,  plot_fit_guess = True, 
            init_states = ['m1','p1'], ro_states=['0','m1','p1']):
   
    ''' 
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''
    fitted_tau, fitted_tau2, fitted_tau_err, fitted_tau2_err = [],[], [], []
    fig = plt.figure()
    ax = plt.subplot()
    plt.xlabel('time (ns)')
    plt.ylabel('p')
    plt.tight_layout()
    fit_results = []

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)

    if ssro_calib_timestamp == None:
        ssro_calib_folder = toolbox.latest_data('SSROCalibration', older_than=older_than)
        if VERBOSE:
            print 'Using SSRO timestamp ', ssro_calib_folder
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.data_from_time(ssro_calib_timestamp)
        if VERBOSE:
            print 'Using SSRO timestamp ', ssro_calib_folder
    for power_elem in powers:
        for init_element in init_states:
            for ro_element in ro_states:
                if timestamp != None:
                    folder = toolbox.data_from_time(timestamp)
                elif len(powers) >1:
                    folder = toolbox.latest_data( \
                        'Repump_'+str(power_elem)+'nW_' \
                        +str(ro_element)+'RO_'+str(init_element)+'init',
                        older_than=older_than)
                else:
                    #folder = toolbox.latest_data('nW_'+ro_element+'RO_'+init_element+'init',
                     #   older_than=older_than)
                    folder = toolbox.latest_data('repump',
                        older_than=older_than)

                print 'folder is ', folder
                plt.title(folder)
                a = mbi.MBIAnalysis(folder)
                a.get_sweep_pts()
                CR_after_check = None
                a.get_readout_results(name='adwindata',CR_after_check = CR_after_check)
                a.get_electron_ROC(ssro_calib_folder)

                x = a.sweep_pts.reshape(-1)[exclude_first_n_points:]
                if ro_element == '0':
                    y = np.array(1.) - a.p0.reshape(-1)[exclude_first_n_points:]
                    #y = a.p0.reshape(-1)[exclude_first_n_points:]
                else:
                    y = a.p0.reshape(-1)[exclude_first_n_points:]
                y_u = a.u_p0.reshape(-1)[exclude_first_n_points:]

                fmt = '.-' if init_element == '0' else 'o-' if init_element == 'm1' else 'x--' 
                color = 'k' if ro_element == '0' else 'r' if ro_element == 'm1' else 'b' 
                if do_plot:
                    if log_plot:
                        ax.set_yscale("log", nonposy='clip')
                        plt.ylim(0.0001,1.05)
                        plt.xlim(-10,np.amax(x))
                    else:
                        plt.ylim(0.0,1.05)
                        plt.xlim(-10,np.amax(x))
                    plt.errorbar(x,y, yerr = y_u, fmt = fmt, color = color, \
                        label = init_element+'init_'+ro_element + 'RO' )
                    plt.legend(numpoints=1, fontsize=legend_fontsize, loc=1,
                        frameon=False, labelspacing=-0.15, borderpad=.5, handletextpad=0, borderaxespad=0)

                #fitfunction: y(x) = A * exp(-x/tau)+ A2 * exp(-x/tau2) + a
                p0, fitfunc, fitfunc_str = common.fit_repumping( offset, amplitude, decay_constant_one,
                        decay_constant_two, x_offs )
                if do_fit:
                    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=print_fit, ret=True,fixed=fixed)
                    if plot_fit == True:
                        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)
                if plot_fit_guess:
                    ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

        
    if do_plot:
        plt.savefig(os.path.join(folder, 'analyzed_result.pdf'), format='pdf')
        plt.savefig(os.path.join(folder, 'analyzed_result.png'), format='png')
        plt.show()

    plt.close('all')


def bin_data(x=[], y=[], y_u=[], binwidth_ns = None):
    
    #Sort data
    sortedx = np.argsort(x)
    y = y[sortedx]
    y_u = y_u[sortedx]
    x = x[sortedx]

    binned_x, binned_y, binned_yu, temp_y, temp_yu= [],[],[],[],[]

    if binwidth_ns != None:
        last_x = 0
        #Calculate mx number of bins
        max_number_of_bins = np.floor(x[-1]/binwidth_ns)+1
        print 'last x:', x[-1]
        # go over bins
        for bin_no in np.arange(max_number_of_bins): 
            #print 'last_x ', last_x, ' bin_no ', bin_no, ' x value: ', bin_no * binwidth_ns
            while True:
                #go over data and pull all data that belongs to a bin into the temp_y list
                if last_x >= len(x):
                    print 'Binning done.'
                    break # quit loop.
                if np.floor(x[last_x] / binwidth_ns) > bin_no: #it's time for the next bin; write values into output array and clear temp array
                    if len(temp_y) > 0:
                        binned_x.append(bin_no * binwidth_ns)
                        binned_y.append(np.sum(temp_y)/len(temp_y))
                        temp_y = []
                    if len(temp_yu) > 0:
                        binned_yu.append(np.sqrt(np.sum(np.square(temp_yu)))/len(temp_yu))
                        temp_yu = []
                    break # quit loop.
                else:
                    temp_y.append(y[last_x])
                    temp_yu.append(y_u[last_x])
                    last_x += 1

    else:
        print 'no binwidth specified'
        binned_x, binned_y, binned_yu = x, y, y_u
    return binned_x, binned_y, binned_yu

def repump_speed_paper_plot(timestamp=None, measurement_name = 'adwindata', ssro_calib_timestamp =None,
            exclude_first_n_points = 0.,
            offset = 0., x0 = 0, older_than= None, newer_than=None, binwidth_ns = None,
            amplitude = 0.8, decay_constant_one = 0.2, decay_constant_two = 0.6, x_offs = 0,
            plot_results = True, do_T2correct=False, labels = [],
            plot_fit = True, do_print = False, fixed = [2], show_guess = True,
            powers = [0], colors=['b'], cutoff_ns = [3500]):
   

    fitted_tau, fitted_tau2, fitted_tau_err, fitted_tau2_err = [],[],[],[]
    CR_after_check = True

    #p0, fitfunc, fitfunc_str = [], [], []
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot()
    fig.text(0.04,0.92, '(a)', fontsize=fignumber_fontsize)
    plt.xlabel('Repump time (ns)', size=axeslabel_fontsize)
    plt.ylabel(r'$1-p_{\left| 0 \right\rangle}$', size=axeslabel_fontsize)
    plt.ylim(0.01,1.05)
    ax.set_yscale("log", nonposy='clip')
    plt.xlim(0.,1000)
    plt.tick_params(pad = 4, axis='both', which='major', labelsize=ticklabel_fontsize, width = tickwidth, length=majorticklength)
    plt.tick_params(pad = 4, axis='both', which='minor', labelsize=ticklabel_fontsize, width = tickwidth, length=minorticklength)
    
    plt.tight_layout()
    

    for count in np.arange(len(older_than)):
        x, y, y_u = [], [], []
        print 'older than ', older_than[count], ' and newer_than ', newer_than[count]
        fit_results = []
        #folder = toolbox.data_from_time(timestamp[count])
        folder_list = toolbox.latest_data('Repump', 
            folder = 'D:\measuring\data',
            older_than=older_than[count],
            newer_than=newer_than[count],
            return_all=True)
        print 'I found ', len(folder_list), ' files:'
        #print 'FolderList: ', folder_list
        folder_list_ext = folder_list
        # for av_elem in folder_list:
        #     folder_list_ext.append(toolbox.latest_data(av_elem))
        #print 'folder is ', folder_list_ext
        
        if ssro_calib_timestamp == None :
            ssro_calib_folder = toolbox.latest_data('SSROCalibration', older_than=older_than[count])
        else:
            ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp[count])
            ssro_calib_folder = toolbox.data_from_time(ssro_calib_timestamp[count])
        if VERBOSE:
            print 'Using SSRO timestamp ', ssro_calib_folder
            
        for elem in np.arange(len(folder_list_ext)):
            #print folder_list_ext[elem]
            a = mbi.MBIAnalysis(folder_list_ext[elem])
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata', CR_after_check = CR_after_check)
            a.get_electron_ROC(ssro_calib_folder)

            x_list = 1000*a.sweep_pts.reshape(-1)[exclude_first_n_points[count]:]
            x = np.append(x, x_list - x_list[0])  #shifts the graph such that it starts at t=0
            y = np.append(y,1-a.p0.reshape(-1)[exclude_first_n_points[count]:])
            y_u = np.append(y_u,a.u_p0.reshape(-1)[exclude_first_n_points[count]:])
            #print 'lengths are: x ', len(x), ' y ', len(y), ' y_u ', len(y_u)

        binned_x, binned_y, binned_yu = bin_data(x, y, y_u, binwidth_ns[count])

        #plt.errorbar(a.sweep_pts[exclude_first_n_points[count]:], 1-a.p0[exclude_first_n_points[count]:,0], yerr = a.u_p0[exclude_first_n_points[count]:,0], fmt = '.',color = colors[count], label = '')
        #print x, y, y_u
        #print 'x ', binned_x, ' cutoff ', cutoff_ns[count]
        elements_to_cut = np.nonzero( binned_x > np.array(cutoff_ns[count]) )[0]
        if elements_to_cut != []:
            print 'cut index ', elements_to_cut[0]
            binned_x = binned_x[0:elements_to_cut[0]]
            binned_y = binned_y[0:elements_to_cut[0]]
            binned_yu = binned_yu[0:elements_to_cut[0]]

        plt.errorbar(binned_x, binned_y, zorder = 500-count, capsize= errorbar_width, yerr = binned_yu, fmt = '', ls= '', color = colors[count], label = labels[count], elinewidth=linewidth)

        #fitfunction: y(x) = A * exp(-x/tau)+ A2 * exp(-x/tau2) + a
        p0, fitfunc, fitfunc_str = common.fit_repumping( 
            offset[count], amplitude[count], decay_constant_one[count],
            decay_constant_two[count], x_offs[count] )
        # p0, fitfunc, fitfunc_str = common.fit_double_exp_decay_with_x_offset(np.array(offset[count]), np.array(amplitude[count]), 
        #          np.array(decay_constant_one[count]), np.array(decay_constant_two[count]), np.array(x_offs[count]))

        if plot_results and show_guess:
            ax.plot(np.linspace(binned_x[0],binned_x[-1],201), fitfunc(np.linspace(binned_x[0],binned_x[-1],201)), ':', lw=2)

        fit_result = fit.fit1d( binned_x, binned_y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, label = labels[count], ret=True, fixed=fixed[count])

        ## plot data and fit as function of total time
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001),color = colors[count],log=True, 
                ax=ax, plot_data=False, legend=False, add_txt=False, lw= linewidth)
        
        plt.legend(numpoints=1, fontsize=legend_fontsize, loc=3, frameon=False, labelspacing=-0.15, borderpad=.5, handletextpad=0, borderaxespad=0)
        
        fit_results.append(fit_result['params_dict'])
        fitted_tau.append(fit_result['params_dict']['tau'])
        fitted_tau_err.append(fit_result['error_dict']['tau'])
        fitted_tau2.append(fit_result['params_dict']['tau2'])
        fitted_tau2_err.append(fit_result['error_dict']['tau2'])
        
        print fit_result['params_dict'], fit_result['error_dict']


        #fitted_tau.append(fit_result['params_dict']['tau'])
        #fitted_tau_err.append(fit_result['error_dict']['tau'])

    if plot_results:
        save_figure_to = 'D:\measuring\QMem_plots'
        #save_figure_to = 'K:\ns\qt\Diamond\Eigenpapers\15-WeaklyCoupledQuantumMemory\Figures'
        print 'saving to: ', save_figure_to
        plt.savefig(os.path.join(save_figure_to, 'Fig2.pdf'), format='pdf')
        plt.savefig(os.path.join(save_figure_to, 'Fig2.png'), format='png')

    print 'tau ', fitted_tau, ' tau err ', fitted_tau_err, ' tau2 ', fitted_tau2, ' tau2 err ', fitted_tau2_err

# def analyze_avg_repump_time(carbons = ['2','3'],folder_name = 'Memory_Sweep_repump_time_',fit_results = False,older_than = older_than):
#     CD.Sweep_Rep_List(carbons = carbons, folder_name = folder_name, fit_results = False, older_than = older_than)

def rate_equations_pm0S(t, transition = '_E'):   #paste results from mathematica here
    if transition == '_E':
        pop_0, pop_m1, pop_p1, pop_S =  \
            0.05833 *np.exp(-0.358328 *t) + 0.823026 *np.exp(-0.01777 *t) + \
             0.118644 *np.exp(-0.0020322 *t), \
            0.0411599 *np.exp(-0.358328 *t) + 0.0407096 *np.exp(-0.35748 *t) + \
             0.45929 *np.exp(-0.0183774 *t) + 0.399694 *np.exp(-0.01777 *t) + \
             0.0591464 *np.exp(-0.0020322 *t), \
            0.0411599 *np.exp(-0.358328 *t) - 0.0407096 *np.exp(-0.35748 *t) - \
             0.45929 *np.exp(-0.0183774 *t) + 0.399694 *np.exp(-0.01777 *t) + \
             0.0591464 *np.exp(-0.0020322 *t), \
            0.00962525 *np.exp(-0.358328 *t) - 0.217911 *np.exp(-0.01777 *t) + \
             0.208285 *np.exp(-0.0020322 *t)

    else:
        pop_0, pop_m1, pop_p1, pop_S =  \
            0.0364711 *np.exp(-0.294573 *t) - 0.000531361 *np.exp(-0.201993 *t) + \
             0.617408 *np.exp(-0.0438856 *t) + 0.211492 *np.exp(-0.011164 *t) + \
             0.13516 *np.exp(-0.00198961 *t), \
            0.0744674 *np.exp(-0.294573 *t) + 0.011632 *np.exp(-0.201993 *t) + \
             0.765363 *np.exp(-0.0438856 *t) + 0.0826231 *np.exp(-0.011164 *t) + \
             0.0659143 *np.exp(-0.00198961 *t), \
            0.00274379 *np.exp(-0.294573 *t) - 0.0110715 *np.exp(-0.201993 *t) - \
             0.186566 *np.exp(-0.0438856 *t) + 0.125965 *np.exp(-0.011164 *t) + \
             0.0689292 *np.exp(-0.00198961 *t), \
            0.027869 *np.exp(-0.294573 *t) + 0.00111591 *np.exp(-0.201993 *t) - \
             0.185802 *np.exp(-0.0438856 *t) - 0.0646157 *np.exp(-0.011164 *t) + \
             0.221433 *np.exp(-0.00198961 *t)
    return pop_p1, pop_m1, pop_0, pop_S

def repump_speed_pm1_paper_plot(timestamp=None, measurement_name = 'adwindata', ssro_calib_timestamp =None,
            exclude_first_n_points_A = 0., exclude_first_n_points_E = 0.,
            offset = 0., x0 = 0,tstamps_A= None, tstamps_E= None, binwidth_ns = None,
            amplitude = 0.8, decay_constant_one = 0.2, decay_constant_two = 0.6, x_offs = 0,
            plot_results = True, do_T2correct=False, marker = '.',
            plot_fit = True, do_print = False, fixed = [2], show_guess = True,
            powers = [0], colors=['b'], cutoff_ns = [3500], labels = []):
   
    fig, axs = plt.subplots(1, 2, sharey=True)

    fig.set_size_inches(figsize)
    fig.subplots_adjust(wspace=0.1)
    E_or_A = ('_E', '_A')
    fig.text(0.005,0.91, '(b)', fontsize=fignumber_fontsize)
    axs[0].set_title(r'$E_{1,2}$ transitions', size=8)
    axs[1].set_title(r'$A_{1,2}$ transitions', size=8)
    fitted_tau, fitted_tau2, fitted_tau_err, fitted_tau2_err = [],[],[],[]
    plt.ylim(0.01,1.01)
    plt.tick_params(pad = 4, axis='both', which='major', labelsize=ticklabel_fontsize, width = tickwidth, length=majorticklength)
    plt.tick_params(pad = 4, axis='both', which='minor', labelsize=ticklabel_fontsize, width = tickwidth, length=minorticklength)
    plt.tight_layout()


    for panel_no in range(2):


        axs[panel_no].tick_params(axis='both', which='major', labelsize=ticklabel_fontsize, width = tickwidth, length = majorticklength)
        axs[panel_no].set_yscale("log", nonposy='clip')
        axs[panel_no].set_xlabel('Repump time (ns)', size=axeslabel_fontsize)
        axs[panel_no].xaxis.set_ticks(np.arange(0, 801,200))
    

        print 'Panel No ', panel_no
        if panel_no == 0:
            tstamps = tstamps_E
            exclude_first_n_points = exclude_first_n_points_E
        else:
            tstamps = tstamps_A
            exclude_first_n_points = exclude_first_n_points_A

        for count in np.arange(len(tstamps)):
            x, y, y_u = [], [], []
            folder = toolbox.data_from_time(tstamps[count])
            if VERBOSE:
                print folder
            if ssro_calib_timestamp == None :
                ssro_calib_folder = toolbox.latest_data('SSROCalibration', older_than=tstamps[count])
            else:
                ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp[count])
                ssro_calib_folder = toolbox.data_from_time(ssro_calib_timestamp)
            if VERBOSE:
                print 'Using SSRO timestamp ', ssro_calib_folder
                
            a = mbi.MBIAnalysis(folder)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata', CR_after_check = CR_after_check)
            a.get_electron_ROC(ssro_calib_folder)

            x_list = 1000*a.sweep_pts.reshape(-1)[exclude_first_n_points[count]:]
            x = np.append(x, x_list - x_list[0])  #shifts the graph such that it starts at t=0
            if '0RO' in folder:
                print 'inverting'
                y = np.append(y,1-a.p0.reshape(-1)[exclude_first_n_points[count]:])
            else:
                print 'not inverting'
                y = np.append(y,a.p0.reshape(-1)[exclude_first_n_points[count]:])
            y_u = np.append(y_u,a.u_p0.reshape(-1)[exclude_first_n_points[count]:])
            #print 'lengths are: x ', len(x), ' y ', len(y), ' y_u ', len(y_u)

            binned_x, binned_y, binned_yu = bin_data(x, y, y_u, binwidth_ns[count])

            elements_to_cut = np.nonzero( binned_x > np.array(cutoff_ns[count]) )[0]
            if elements_to_cut != []:
                print 'cut index ', elements_to_cut[0]
                binned_x = binned_x[0 : elements_to_cut[0]]
                binned_y = binned_y[0 : elements_to_cut[0]]
                binned_yu = binned_yu[0 : elements_to_cut[0]]

            #axs[panel_no].errorbar(binned_x, binned_y, 
            #    yerr = binned_yu, fmt = '', ls= '', 
            #    color = colors[count], label = labels[count], elinewidth=2)
            axs[panel_no].plot(binned_x, binned_y, marker[count], color = colors[count], label = labels[count], markersize = markersize)

            if 'p1RO' in folder:
                print 'fitting p1'
                p0, fitfunc, fitfunc_str = common.fit_repumping_p1( 
                    offset[count], amplitude[count],amplitude[count], decay_constant_one[count],
                    decay_constant_two[count], x_offs[count] )
            else:
                p0, fitfunc, fitfunc_str = common.fit_repumping( 
                    offset[count], amplitude[count], decay_constant_one[count],
                    decay_constant_two[count], x_offs[count] )

            ## plot Model data

            if False:
                plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001),color = colors[count],log=True, 
                    ax=ax, plot_data=False, legend=None, add_txt=False)
            
            
            try:
                fit_results.append(fit_result['params_dict'])
                fitted_tau.append(fit_result['params_dict']['tau'])
                fitted_tau_err.append(fit_result['error_dict']['tau'])
                fitted_tau2.append(fit_result['params_dict']['tau2'])
                fitted_tau2_err.append(fit_result['error_dict']['tau'])
                
                print fit_result['params_dict'], fit_result['error_dict']
            except:
                print 'Fit didnt fit.'

            with open(os.path.join(save_figure_to,str(tstamps[count])+'data_2b.csv'), 'w') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(x)
                writer.writerow(y)
                writer.writerow(y_u)
            #fitted_tau.append(fit_result['params_dict']['tau'])
            #fitted_tau_err.append(fit_result['error_dict']['tau'])
        model_x = np.arange(cutoff_ns[count])
        model_pop_p, model_pop_m, model_pop_0, model_pop_S = [],[],[],[]
        for time in model_x:
            pop_p,pop_m,pop_0,pop_S=rate_equations_pm0S(time, E_or_A[panel_no])
            model_pop_p.append(pop_p)
            model_pop_m.append(pop_m)
            model_pop_0.append(pop_0)
            model_pop_S.append(pop_S)
        axs[panel_no].plot(model_x, model_pop_p, color = colors[0], lw=linewidth, zorder = 1)
        axs[panel_no].plot(model_x, model_pop_m, color = colors[1], lw=linewidth, zorder = 1)
        axs[panel_no].plot(model_x, model_pop_0, color = colors[2], lw=linewidth, zorder = 1)
        axs[panel_no].plot(model_x, model_pop_S, ':', color = colors[3], lw=linewidth, zorder = 1, label = r'$| S \rangle$')
        axs[panel_no].xaxis.set_tick_params(width=tickwidth, length=majorticklength)
        axs[panel_no].yaxis.set_tick_params(width=tickwidth, length=majorticklength)
        axs[panel_no].yaxis.set_tick_params('minor', width=tickwidth, length=minorticklength)
        axs[panel_no].legend(numpoints=1, fontsize=legend_fontsize,
            borderpad=.5, handletextpad=0, borderaxespad=0.5, loc=1, frameon=False,  labelspacing=0.)

        if plot_results:
            print 'saving to: ', save_figure_to
            plt.savefig(os.path.join(save_figure_to, 'Fig2b.pdf'), format='pdf')
            plt.savefig(os.path.join(save_figure_to, 'Fig2b.png'), format='png')
        print 'tau ', fitted_tau, ' tau err ', fitted_tau_err, ' tau2 ', fitted_tau2, ' tau2 err ', fitted_tau2_err
    plt.close('all')    

def ionization_paper_plot(timestamps=None, measurement_name = 'adwindata', ssro_calib_timestamp =None,
            exclude_first_n_points = 0., binwidth_ns = None,
            amplitude = 1, decay_constant_guess = 3000, decay_constant_two = 0.6, x_offs = 0,
            plot_results = True, do_T2correct=False, labels = [''], log_plot = True,
            do_fit = True, do_print = False, fixed = [2], show_guess = True,save_fig = True,
            colors=['b']):
   

    fitted_tau, fitted_tau_err = [],[]
    CR_after_check = True


    fig = plt.figure(figsize)
    ax = plt.subplot()
    # fig.text(0.04,0.92, '(a)', fontsize=fignumber_fontsize)
    plt.xlabel('Repetitions', size=axeslabel_fontsize)
    plt.ylabel(r'$p_{NV^-}$', size=axeslabel_fontsize)
    # ax.set_yscale("log", nonposy='clip')

    plt.tick_params(pad = 4, axis='both', which='major', labelsize=ticklabel_fontsize, width = tickwidth, length=majorticklength)
    plt.tick_params(pad = 4, axis='both', which='minor', labelsize=ticklabel_fontsize, width = tickwidth, length=minorticklength)
    
    plt.tight_layout()
    
    for count in np.arange(len(timestamps)):
        x, y, y_u, fit_results = [], [], [], []

        folder = toolbox.data_from_time(timestamps[count])

        if ssro_calib_timestamp == None :
            ssro_calib_folder = toolbox.latest_data('SSROCalibration', older_than=timestamps[count])
        else:
            ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp[count])
            ssro_calib_folder = toolbox.data_from_time(ssro_calib_timestamp)
        if VERBOSE:
            print 'Using SSRO timestamp ', ssro_calib_folder
            
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata', CR_after_check = CR_after_check)
        a.get_electron_ROC(ssro_calib_folder)

        x_list = a.sweep_pts.reshape(-1)[exclude_first_n_points[count]:]
        x = np.append(x, x_list - x_list[0])  #shifts the graph such that it starts at t=0
        y = np.append(y, a.p0.reshape(-1)[exclude_first_n_points[count]:])
        y_u = np.append(y_u,a.u_p0.reshape(-1)[exclude_first_n_points[count]:])

        if binwidth_ns != None:
            binned_x, binned_y, binned_yu = bin_data(x, y, y_u, binwidth_ns[count])
        else:
            binned_x, binned_y, binned_yu = x, y, y_u
        if VERBOSE:
            print 'binned x and y: ', binned_x, binned_y

        plt.ylim(0.1,1.0)
        plt.xlim(0.,1.05*max(binned_x))

        plt.errorbar(binned_x, binned_y, zorder = 500-count, yerr = binned_yu, capsize= errorbar_width, fmt = '.', ls= '', color = colors[count], label = labels[count], elinewidth=linewidth)

        if do_fit:
            p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(0.,1.,decay_constant_guess[count])
            print  p0, fitfunc, fitfunc_str 
            fit_result = fit.fit1d( binned_x, binned_y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True, fixed = fixed)
        if plot_results and show_guess:
            ax.plot(np.linspace(binned_x[0],binned_x[-1],201), fitfunc(np.linspace(binned_x[0],binned_x[-1],201)), ':', lw=2)

        if do_fit:  ## plot data and fit
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001),color = colors[count],log=True, 
                ax=ax, plot_data=False, legend=False, add_txt=False, lw= linewidth)
        
        plt.legend(numpoints=1, fontsize=legend_fontsize, loc=3, frameon=False, labelspacing=0.3, borderpad=.5, handletextpad=0, borderaxespad=0)

    if save_fig:
        save_figure_to = 'D:\measuring\QMem_plots'
        #save_figure_to = 'K:\ns\qt\Diamond\Eigenpapers\15-WeaklyCoupledQuantumMemory\Figures'
        print 'saving to: ', save_figure_to
        plt.savefig(os.path.join(save_figure_to, 'Supp_Ionization.pdf'), format='pdf')
        plt.savefig(os.path.join(save_figure_to, 'Supp_Ionization.png'), format='png')

def coupling_vs_rep_paper_plot(c_idents = ['1'], do_Z=False, older_than_list= None, older_than_Z=None, 
    labels = [], styles=['ko'], fit_colors = [], LogPlot=True, fixed = [0],
    update_data=True, lastdata = None, do_T2correct=False, **kw):
    
    arraylength=0
    for c_ident_count in c_idents:  # add one for singles and two (pm configuration) for DFS
        arraylength += len(c_ident_count)

    singles_idents = []
    for singels_elem in c_idents:
        if len(singels_elem)==1:
            singles_idents.append(singels_elem)

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot()
   
    if not do_Z:
        fig.text(0.02,0.9, '(b)', fontsize=fignumber_fontsize)
    plt.tick_params(pad = 4, axis='both', which='major', labelsize=ticklabel_fontsize, width = tickwidth, length=majorticklength)
    plt.tick_params(pad = 4, axis='both', which='minor', labelsize=ticklabel_fontsize, width = tickwidth, length=minorticklength)
    plt.xlabel(r'Coupling strength $|\Delta \omega|/(2\pi)$ (kHz)', size=axeslabel_fontsize)   
    plt.ylabel(r'Bloch vector decay constant $N_{1/e}$', size=axeslabel_fontsize, labelpad =2)
    plt.title('')

    plt.tight_layout()

    x, y, y_u, z, z_u= np.zeros( (len(older_than_list),arraylength) ), np.zeros( (len(older_than_list),arraylength) ), np.zeros( (len(older_than_list),arraylength) ), np.zeros( (len(older_than_list),arraylength) ), np.zeros( (len(older_than_list),arraylength) )
    singles_x, singles_y, singles_y_u = np.zeros( (len(older_than_list),len(singles_idents)) ), np.zeros( (len(older_than_list),len(singles_idents)) ), np.zeros( (len(older_than_list),len(singles_idents)) )
    if LogPlot:
        if older_than_Z == None:
            plt.ylim([10,4000])
        else:
            plt.ylim([10,50000])
        plt.xlim([0,90])
        #ax.set_xscale("log", nonposy='clip')
        ax.set_yscale("log", nonposy='clip')
    else:
        plt.ylim([0,1000])

    do_T2correct_list = kw.pop('do_T2correct_list', None)

    for count in np.arange(len(older_than_list)):
        if do_T2correct_list != None:
            do_T2correct = do_T2correct_list[count]
        if update_data or lastdata == None:
            print 'processing timestamps before ', older_than_list[count]
            x[count], y[count], y_u[count], folder_temp = coupling_vs_rep_update_data(c_idents=c_idents, older_than = older_than_list[count] , do_Z = False, do_T2correct=do_T2correct)
            singles_x[count], singles_y[count], singles_y_u[count], folder_temp = coupling_vs_rep_update_data(c_idents=singles_idents, older_than = older_than_list[count], do_Z = False, do_T2correct=do_T2correct)
        else:
            x, y, y_u, singles_x, singles_y, singles_y_u = lastdata
        
        guess_values= kw.get('fit_guess',[(1,0.000435,19)])
        p0,fitfunc,fitfunc_str = common.fit_dephasing_coupl(guess_values[count][0],guess_values[count][1],guess_values[count][2])

        if count < len(fit_colors):  #do fit
            x_tempfit, y_tempfit, y_u_tempfit  = [],[],[]
            for fit_count in np.arange(len(x[count])):
                if y_u[count][fit_count]/y[count][fit_count] < 0.9: # exclude data points with too large error
                    x_tempfit.append(x[count][fit_count])
                    y_tempfit.append(y[count][fit_count])
                    y_u_tempfit.append(y_u[count][fit_count])
                else:
                    print 'excluding from fit: ', x_tempfit
            if False:  # Show fit guess
                ax.plot(np.linspace(2,90,201), fitfunc(np.linspace(2,90,201)), ':', lw=linewidth, color = 'r')
            fit_result = fit.fit1d(x_tempfit, y_tempfit, None, p0 = p0, fitfunc = fitfunc, do_print = True, ret = True, fixed = fixed[count])
            print 'fit_result ', fit_result 
            plot.plot_fit1d(fit_result, np.linspace(0,100,1001), label = labels[count], ax=ax, color = fit_colors[count], plot_data=False, add_txt = False, lw = linewidth)

        
    if do_Z:
        print 'processing Z, timestamps before ', older_than_Z
        x_zmeas, z, z_u, folder_temp = coupling_vs_rep_update_data(c_idents=c_idents, do_Z=True, older_than = older_than_Z , do_T2correct=False)
        singles_x_zmeas, singles_z, singles_z_u, folder_temp = coupling_vs_rep_update_data(c_idents=singles_idents, older_than = older_than_Z, do_Z = True, do_T2correct=False)

    
    if True:  # Show Guess from repump time
        if VERBOSE:
            print 'Plotting Guess'
        fit_guess_art= kw.get('fit_guess_art',[(1,0.000435,19)])
        if fit_guess_art != []:
            p0,fitfunc,fitfunc_str = common.fit_dephasing_coupl(fit_guess_art[0],fit_guess_art[1],fit_guess_art[2])
            ax.plot(np.linspace(2,90,201), fitfunc(np.linspace(2,90,201)), ':', lw=linewidth, color = 'r')
    if do_Z:  #plot Z
        ax.errorbar(x_zmeas,z,z_u, fmt = '.', color='r',
            zorder = 500, label = 'Z decay', capsize= errorbar_width, linewidth=linewidth, markeredgewidth = 1, markeredgecolor = 'r', markersize=2)
        ax.plot(singles_x_zmeas,singles_z, 'o', mfc='white', color = 'r',
            zorder = 500, markeredgecolor = 'r', markersize=2)
        
    for count in np.arange(len(older_than_list)):   
        ax.errorbar(x[count],y[count],y_u[count], fmt = '.', color=fit_colors[count],
            zorder = 500-count, label = labels[count], capsize= errorbar_width, linewidth=linewidth, markeredgewidth = 1, markeredgecolor = fit_colors[count], markersize=2)
        ax.plot(singles_x[count], singles_y[count], 'o', mfc='white', color = fit_colors[count],
            zorder = 500-count, markeredgecolor = fit_colors[count], markersize=2)
        leg=ax.legend(ncol=1, loc= kw.get('legend_pos', 1), fontsize=legend_fontsize, labelspacing=0.3, 
            borderpad=.5, handletextpad=0, borderaxespad=0,
            numpoints=1, title = ' ', frameon=False)
        
        # swap_and_right_align_legend
        vp = leg._legend_box._children[-1]._children[0]
        for c in vp._children:
            c._children.reverse()
        vp.align="right" 
        
    ax.xaxis.set_tick_params(width=1, length=2)
    ax.yaxis.set_tick_params(width=1, length=2)
    ax.yaxis.set_tick_params('minor', width=1, length=1)

    plt.xticks(np.arange(0,90,20))
    

    #plt.legend(loc=1)#bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    save_figure_to = 'D:\measuring\QMem_plots'
    #save_figure_to = 'K:\ns\qt\Diamond\Eigenpapers\15-WeaklyCoupledQuantumMemory\Figures'
    print 'saving to: ', save_figure_to
    if older_than_Z == None:
        plt.savefig(os.path.join(save_figure_to, 'Fig4b.pdf'), format='pdf')
        plt.savefig(os.path.join(save_figure_to, 'Fig4b.png'), format='png')
    else:
        plt.savefig(os.path.join(save_figure_to, 'Supp_Zdecay.pdf'), format='pdf')
        plt.savefig(os.path.join(save_figure_to, 'Supp_Zdecay.png'), format='png')
    
    plt.show()
    plt.close('all')

    return x, y, y_u, singles_x, singles_y, singles_y_u

def decay_vs_rep_paper_plot(DPS = False, carbons = ['1','2'],
        older_than = None, 
        do_T2correct = False,
        folder_name = 'Repetitions_', plotclassical = False,
        ssro_calib_timestamp = None, plotlabel='', **kw):

    ## other key word arguments
    x_max = kw.pop('x_max', None)
    fig_name = kw.pop('fig_name', 'Fig')
    tick_spacing = kw.pop('tick_spacing', 500)
    logicstate_list = kw.pop('logicstate_list',len(carbons)*['X']) ## can be list such as ['X','mX'], used for DFS measurements.
    colors = kw.pop('colors', color_list)
    is_X_measurement = kw.get('is_X_measurement', True)
    labels=[]

    x_arr, y_arr, y_u_arr = [],[],[]
    for c,logicstate in zip(carbons,logicstate_list):
        folder_dict= extract_data_from_sweep(older_than = older_than,
            folder_name =folder_name, carbon = c,
            ssro_calib_timestamp =ssro_calib_timestamp,
            logicstate = logicstate,
            do_T2correct=do_T2correct, **kw)
        x_arr.append(folder_dict['sweep_pts'])
        y_arr.append(folder_dict['res'])
        y_u_arr.append(folder_dict['res_u'])
        folder_dict, coupling, folder = extract_coupling_strength(folder_dict)
        labels.append(np.round(coupling/100.)/10)
    
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot()
    fig.text(0.03,0.9, plotlabel, fontsize=fignumber_fontsize)
    ax.tick_params(pad=4)
    plt.tick_params(pad = 4, axis='both', which='major', labelsize=ticklabel_fontsize, width = tickwidth, length=majorticklength)
    plt.tick_params(pad = 4, axis='both', which='minor', labelsize=ticklabel_fontsize, width = tickwidth, length=minorticklength)



    if is_X_measurement:
        for key in ['XX','X']:
            if folder_dict[key] != []:
                folder = folder_dict[key][0]
    else:
        for key in ['ZZ','Z']:
            if folder_dict[key] != []:
                folder = folder_dict[key][0]
                print 'measuring Z'

    for x,y,y_u,carbon,logicstate,jj in zip(x_arr,y_arr,y_u_arr,carbons,logicstate_list,range(len(x_arr))):
        A0, offset, decay, x0 = y[0], 0, 50, 0
        p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)
        fixed = [0,3]
        fit_result = fit.fit1d(x,y,None,p0 = p0, fitfunc = fitfunc, do_print = False, ret = True, fixed = fixed)
        if x_max == None:
            x_max = 2 * x[-1]
        if carbon =='25' and logicstate == 'X':
            fit_range = 400
        else:
            fit_range = x_max
        plot.plot_fit1d(fit_result, np.linspace( 0*x[0], fit_range, 1001), zorder = 400, ax=ax,color = colors[jj], plot_data=False,add_txt = False, lw = linewidth)
        plt.errorbar(x,y,y_u, zorder = 500, capsize= errorbar_width, fmt='.',lw=1,color = colors[jj], label=labels[jj], elinewidth=linewidth, markersize=markersize)

    plt.ylim(0.07,1.0)
    ax.set_yscale("log", nonposy='clip')
    plt.errorbar(np.linspace(0, x_max, num=5), np.linspace(1./3., 1./3., num=5), None, 
        zorder = 1, capsize= errorbar_width,elinewidth=linewidth,markersize=markersize, fmt = '--r', lw=linewidth)

    plt.xlabel(r'Number of repetitions $N$', size=axeslabel_fontsize, labelpad =2)
    plt.xlim(0, x_max)
    plt.ylabel('Bloch vector length', size=axeslabel_fontsize, labelpad =2)
    plt.title('')
    leg=plt.legend(ncol=2, loc= 1 if DPS else 1, fontsize=legend_fontsize, labelspacing=0.3, 
        borderpad=0.5, handletextpad=0, borderaxespad=0.2,
        numpoints=1, title = r'$|\Delta \omega |/ (2 \pi)$(kHz)', frameon=False)
    leg.get_title().set_fontsize(legend_fontsize)
    plt.tight_layout()
    plt.xticks(np.arange(0, x_max+1, tick_spacing))
    plt.minorticks_on() 
    if plotclassical:
        fig.text(0.8,0.58, 'classical', fontsize=legend_fontsize, color='r')
        fig.text(0.8,0.63, 'quantum', fontsize=legend_fontsize, color='r')
    if DPS:
        fig.text(0.67,0.53, r'$(\left |\uparrow \downarrow \right \rangle + \left | \downarrow \uparrow \right \rangle)/\sqrt{2}$', fontsize=legend_fontsize, color='r')
        fig.text(0.16,0.23, r'$(\left |\uparrow \uparrow \right \rangle + \left | \downarrow \downarrow \right \rangle)/\sqrt{2}$', fontsize=legend_fontsize, color='m')

    save_figure_to = 'D:\measuring\QMem_plots'
    #save_figure_to = 'K:\ns\qt\Diamond\Eigenpapers\15-WeaklyCoupledQuantumMemory\Figures'
    print 'saving to: ', save_figure_to
    plt.savefig(os.path.join(save_figure_to, fig_name+'.pdf'), format='pdf')
    plt.savefig(os.path.join(save_figure_to, fig_name+'.png'), format='png')
    plt.show()
    plt.close('all')

def t_sweep_paper_plot( carbons = ['1','2'],
        older_than = None, 
        do_T2correct = False,
        folder_name = 'Repetitions_', 
        ssro_calib_timestamp = None,**kw):

    ## other key word arguments
    x_max = kw.pop('x_max', None)
    labels = kw.get('labels',[])
    fig_name = kw.pop('fig_name', 'Fig')
    fit_results = kw.pop('fit_results',True)
    sequence_length = kw.pop('sequence_length',None)
    colors = kw.pop('colors', color_list)
    x_arr, y_arr, y_u_arr = [], [], []

    for c in carbons:
        folder_dict= extract_data_from_sweep(older_than = older_than,
            folder_name = folder_name, carbon = c, logicstate = 'Z',
            ssro_calib_timestamp =ssro_calib_timestamp,
            do_T2correct=do_T2correct, **kw)
        if VERBOSE:
            print 'folder_dict ', folder_dict
        x_arr.append(folder_dict['sweep_pts'])
        y_arr.append(folder_dict['res'])
        y_u_arr.append(folder_dict['res_u'])
        labels.append(str(c))
    
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot()
    #fig.text(0.03,0.9, '(a)', fontsize=fignumber_fontsize)
    ax.tick_params(pad=4)
    plt.tick_params(pad = 4, axis='both', which='major', labelsize=ticklabel_fontsize, width = tickwidth, length=majorticklength)
    plt.tick_params(pad = 4, axis='both', which='minor', labelsize=ticklabel_fontsize, width = tickwidth, length=minorticklength)
    plt.xticks(np.arange(0, 3001, 1000))
    plt.minorticks_on()
    plt.yticks(np.arange(0,1.001,0.2))

    is_X_measurement = kw.get('is_X_measurement', True)
    if is_X_measurement:
        for key in ['XX','X']:
            if folder_dict[key] != []:
                folder = folder_dict[key][0]
    else:
        for key in ['ZZ','Z']:
            if folder_dict[key] != []:
                folder = folder_dict[key][0]
    if VERBOSE:
        print 'y_arr ', y_arr
    for x,y,y_u,carbon,jj in zip(x_arr,y_arr,y_u_arr,carbons,range(len(x_arr))):
        x = x*1e9 #convert to ns
        if fit_results:
            A0, offset, decay, x0 = 1, 0., 500, 500
            if sequence_length != None:
                decay = decay*sequence_length
            #p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)
            p0,fitfunc,fitfunc_str = common.fit_gauss(offset,A0,x0,decay)
            fixed = [0]

            fit_result = fit.fit1d(x,y,None,p0 = p0, fitfunc = fitfunc, do_print = False, ret = True, fixed = fixed)
            print fit_result['params_dict']['x0'], ' +- ', fit_result['error_dict']['x0'] 
            if x_max == None:
                x_max = 1.1* x[-1]
            plot.plot_fit1d(fit_result, np.linspace( 1.05*x[0], 1.05*x[-1], 1001), ax=ax,color = colors[jj], plot_data=False,add_txt = False, lw = linewidth)

        if VERBOSE:
            print 'xy ', x, y
        plt.plot((2300,2300),(0,1), ls=':',color = 'r') # optimized time
        plt.errorbar(x,abs(y),yerr=y_u, zorder = 500, capsize= errorbar_width, fmt='o-',lw=linewidth,color = colors[jj],label=labels[jj], elinewidth = linewidth, markersize=markersize)

    plt.ylim(0.0,1)
    plt.xlim(0, 1.05*x[-1])

    #plt.errorbar( np.linspace(0, x_max, num=5),np.linspace(1./3., 1./3., num=5), None, fmt = '--y', lw=2)


    plt.xlabel(r't (ns)', size=axeslabel_fontsize, labelpad =2)
    #plt.xlim(0, x_max)


    plt.ylabel('Bloch vector length',size=axeslabel_fontsize, labelpad =2)
    plt.title('')
    leg=plt.legend(loc=(0.2,0.65), fontsize=legend_fontsize, labelspacing=0.3, 
        borderpad=.5, handletextpad=0.2, borderaxespad=0.2,
        numpoints=1, title = r'$\Delta \omega / (2 \pi)$ (kHz)', frameon=False)
    leg.get_title().set_fontsize(legend_fontsize)
    plt.tight_layout()
    
    print 'saving to: ', save_figure_to
    plt.savefig(os.path.join(save_figure_to, fig_name+'.pdf'), format='pdf')
    plt.savefig(os.path.join(save_figure_to, fig_name+'.png'), format='png')
    plt.show()
    plt.close('all')

def avg_repump_time_paper_plot( carbons = ['1','2'],
        older_than = None, 
        do_T2correct = False,
        folder_name = 'Repetitions_', 
        ssro_calib_timestamp = None,**kw):

    ## other key word arguments
    x_max = kw.pop('x_max', None)
    fig_name = kw.pop('fig_name', 'Fig')
    fit_results = kw.pop('fit_results',True)
    sequence_length = kw.pop('sequence_length',None)
    logicstate_list = kw.pop('logicstate_list',len(carbons)*['X']) ## can be list such as ['X','mX'], used for DFS measurements.
    colors = kw.pop('colors', color_list)

    x_arr, y_arr, y_u_arr = [], [], []
    labels = []

    for c,logicstate in zip(carbons,logicstate_list):
        folder_dict= extract_data_from_sweep(older_than = older_than,
            folder_name =folder_name, carbon = c,
            ssro_calib_timestamp =ssro_calib_timestamp,
            logicstate = logicstate,
            do_T2correct=do_T2correct, **kw)
        x_arr.append(folder_dict['sweep_pts'])
        y_arr.append(folder_dict['res'])
        y_u_arr.append(folder_dict['res_u'])
        folder_dict, coupling, folder = extract_coupling_strength(folder_dict)
        labels.append(np.round(coupling/100.)/10)
    ### convert to time instead of repetitions:
    if sequence_length != None:
        x_arr = [x*sequence_length for x in x_arr]
    
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot()
    fig.text(0.03,0.9, '(a)', fontsize=fignumber_fontsize)
    ax.tick_params(pad=4)
    plt.tick_params(pad = 4, axis='both', which='major', labelsize=ticklabel_fontsize, width = tickwidth, length=majorticklength)
    plt.tick_params(pad = 4, axis='both', which='minor', labelsize=ticklabel_fontsize, width = tickwidth, length=minorticklength)
    plt.xticks(np.arange(-1000, 2001, 1000))
    plt.minorticks_on()
    plt.yticks(np.arange(0,0.7,0.2))

    is_X_measurement = kw.get('is_X_measurement', True)
    if is_X_measurement:
        for key in ['XX','X']:
            if folder_dict[key] != []:
                folder = folder_dict[key][0]
    else:
        for key in ['ZZ','Z']:
            if folder_dict[key] != []:
                folder = folder_dict[key][0]

    for x,y,y_u,carbon,logicstate,jj in zip(x_arr,y_arr,y_u_arr,carbons,logicstate_list,range(len(x_arr))):
        x = x*1000 #convert to ns
        if fit_results:
            A0, offset, decay, x0 = 1, 0., 500, 500
            if sequence_length != None:
                decay = decay*sequence_length
            #p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)
            p0,fitfunc,fitfunc_str = common.fit_gauss(offset,A0,x0,decay)
            fixed = [0]

            fit_result = fit.fit1d(x,y,None,p0 = p0, fitfunc = fitfunc, do_print = False, ret = True, fixed = fixed)
            print fit_result['params_dict']['x0'], ' +- ', fit_result['error_dict']['x0'] 
            if x_max == None:
                x_max = 1.1* x[-1]
            plot.plot_fit1d(fit_result, np.linspace( 1.05*x[0], 1.05*x[-1], 1001), ax=ax,color = colors[jj], plot_data=False,add_txt = False, lw = linewidth)


        if len(carbon)!=1:
            label_txt = label_txt+'_'+logicstate

        plt.errorbar(x,y,y_u, zorder = 500, capsize= errorbar_width, fmt='.',lw=linewidth,color = colors[jj],label=labels[jj], elinewidth = linewidth, markersize=markersize)

    plt.ylim(0.0,0.7)
    plt.xlim(1.05*x[0], 1.05*x[-1])

    #plt.errorbar( np.linspace(0, x_max, num=5),np.linspace(1./3., 1./3., num=5), None, fmt = '--y', lw=2)


    plt.xlabel(r'$\tau$ (ns)', size=axeslabel_fontsize, labelpad =2)
    #plt.xlim(0, x_max)


    plt.ylabel('Bloch vector length',size=axeslabel_fontsize, labelpad =2)
    plt.title('')
    leg=plt.legend(loc=2, fontsize=legend_fontsize, labelspacing=0.3, 
        borderpad=.5, handletextpad=0, borderaxespad=0,
        numpoints=1, title = r'$|\Delta \omega|  / (2 \pi)$ (kHz)', frameon=False)
    leg.get_title().set_fontsize(legend_fontsize)
    plt.tight_layout()
    
    print 'saving to: ', save_figure_to
    plt.savefig(os.path.join(save_figure_to, fig_name+'.pdf'), format='pdf')
    plt.savefig(os.path.join(save_figure_to, fig_name+'.png'), format='png')
    plt.show()
    plt.close('all')


def coupling_vs_rep_update_data( c_idents=['1'], older_than=None, do_Z = False, do_T2correct=False):

    if do_Z:
        x_zmeas,z,z_u, folder = Z_decay_vs_coupling(c_idents, perpendicular = False, older_than=older_than, return_vals = True, is_X_measurement=False)
        #x,y,y_u,folder = coupling_vs_repetitions(c_idents,older_than = older_than,
        #    folder_name = 'Memory_NoOfRepetitions_', do_T2correct=do_T2correct)
        return x_zmeas, z, z_u, folder
    else:    
        x,y,y_u,folder = coupling_vs_repetitions(c_idents,older_than = older_than,
            folder_name = 'Memory_NoOfRepetitions_', do_T2correct=do_T2correct)
        return x, y, y_u, folder



def ionization_probabilities(timestamps=None, measurement_name = 'adwindata', ssro_calib_timestamp =None,
            exclude_first_n_points = 0., binwidth_ns = None,
            amplitude = 1, decay_constant_guess = 3000, decay_constant_two = 0.6, x_offs = 0,
            plot_results = True, do_T2correct=False, labels = [''], log_plot = True,
            do_fit = True, do_print = False, fixed = [2], show_guess = True,save_fig = True,
            colors=['b'],**kw):
   

    fitted_tau, fitted_tau_err = [],[]
    CR_after_check = kw.pop('CR_after_check',False)

    fig = plt.figure()
    ax = plt.subplot()
    # fig.text(0.04,0.92, '(a)', fontsize=fignumber_fontsize)
    plt.xlabel('Repetitions')
    plt.ylabel(r'$p_{NV^-}$')
    # ax.set_yscale("log", nonposy='clip')

    # plt.tick_params(pad = 4, axis='both', which='major', labelsize=ticklabel_fontsize, width = tickwidth, length=majorticklength)
    # plt.tick_params(pad = 4, axis='both', which='minor', labelsize=ticklabel_fontsize, width = tickwidth, length=minorticklength)
    
    plt.tight_layout()
    
    for count in np.arange(len(timestamps)):
        x, y, y_u, fit_results = [], [], [], []

        folder = toolbox.data_from_time(timestamps[count])

        if ssro_calib_timestamp == None :
            ssro_calib_folder = toolbox.latest_data('SSROCalibration', older_than=timestamps[count])
        else:
            ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp[count])
            ssro_calib_folder = toolbox.data_from_time(ssro_calib_timestamp)
        if VERBOSE:
            print 'Using SSRO timestamp ', ssro_calib_folder
            
        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata', CR_after_check = CR_after_check)
        a.get_electron_ROC(ssro_calib_folder)

        x_list = a.sweep_pts.reshape(-1)[exclude_first_n_points[count]:]
        x = np.append(x, x_list - x_list[0])  #shifts the graph such that it starts at t=0
        y = np.append(y, a.p0.reshape(-1)[exclude_first_n_points[count]:])
        y_u = np.append(y_u,a.u_p0.reshape(-1)[exclude_first_n_points[count]:])

        if binwidth_ns != None:
            binned_x, binned_y, binned_yu = bin_data(x, y, y_u, binwidth_ns[count])
        else:
            binned_x, binned_y, binned_yu = x, y, y_u
        if VERBOSE:
            print 'binned x and y: ', binned_x, binned_y

        plt.ylim(0.1,1.0)
        plt.xlim(0.,1.05*max(binned_x))

        plt.errorbar(binned_x, binned_y, zorder = 500-count, yerr = binned_yu,  fmt = 'o', ls= '', color = colors[count], label = labels[count])

        if do_fit:
            p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(0.,1.,decay_constant_guess[count])
            print  p0, fitfunc, fitfunc_str 
            fit_result = fit.fit1d( binned_x, binned_y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True, fixed = fixed)
        if plot_results and show_guess:
            ax.plot(np.linspace(binned_x[0],binned_x[-1],201), fitfunc(np.linspace(binned_x[0],binned_x[-1],201)), ':', lw=2)

        if do_fit:  ## plot data and fit
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001),color = colors[count],log=log_plot, 
                ax=ax, plot_data=False, legend=False, add_txt=False)
        
        # plt.legend(numpoints=1, fontsize=legend_fontsize, loc=3, frameon=False, labelspacing=0.3, borderpad=.5, handletextpad=0, borderaxespad=0)
