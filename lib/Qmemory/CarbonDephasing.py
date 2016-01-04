import numpy as np
import os,h5py
from analysis.lib.tools import toolbox; reload(toolbox)
from analysis.lib.m2.ssro import mbi
reload(mbi)
from matplotlib import pyplot as plt
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
reload(fit)
reload(common)
import string
import analysis.lib.QEC.hyperfine_params as hf ### used for perp_coupling vs ZZ

color_list = ['b','g','y','r','brown','m','c']
CR_after_check = False ### discard events with ionization for data analysis? (this relies on the CR check after the SSRO.)

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

    # print tomos
    for t in ['X','Y','XX','XY','YX','YY']:
        data_dict.update({t:[]})
        data_dict.update({t+'_u':[]})

    for t in tomos:
        for i,f in enumerate(folder_dict[t]):
            a = mbi.MBIAnalysis(f)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata',CR_after_check = False)
            a.get_electron_ROC(ssro_calib_folder)
            a.get_sequence_length()

            
            x_labels = a.sweep_pts.reshape(-1)


            if i == 0:
                data_dict[t] = ((a.p0.reshape(-1))-0.5)*2
                data_dict[t+'_u'] = 2*a.u_p0.reshape(-1)
            else:
                y = ((a.p0.reshape(-1))-0.5)*2
                y_u = 2*a.u_p0.reshape(-1)
                data_dict[t] = [y0/2-y[ii]/2 for ii,y0 in enumerate(data_dict[t])]
                data_dict[t+'_u'] = [np.sqrt(y0**2+y_u[ii]**2)/2 for ii,y0 in enumerate(data_dict[t+'_u'])]


    
    ## one carbon experiment
    if len(tomos[0]) ==1:
        npY = np.array(data_dict['Y'])
        npX = np.array(data_dict['X'])
        npY_u = np.array(data_dict['Y_u'])
        npX_u = np.array(data_dict['X_u'])
        if kw.pop('do_get_sequence_length', False):
            return x_labels,npX,npY,npX_u,npY_u, 2*a.repump_wait[1]+a.fast_repump_duration[1] # XXXXX check
        else:
            return x_labels,npX,npY,npX_u,npY_u

    elif len(tomos[0]) == 2:
        data_dict['sweep_pts'] = x_labels
        if kw.pop('do_get_sequence_length', False):
            return data_dict, 2*a.repump_wait[1]+a.fast_repump_duration[1]# XXXXX check
        else:
            return data_dict

def extract_data_from_sweep(older_than = None,
        folder_name ='Repetitions_',
        carbon = '2',
        ssro_calib_timestamp =None, 
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

    for t in ['X','Y','XX','XY','YX','YY', 'Z', 'ZZ']:
        folder_dict.update({t:[]})

    logicstate = kw.get('logicstate',None)

    ### search data
    if len(carbon) ==1:
        for ro in ['positive','negative']:
            for t in ['X','Y']:
                search_string = folder_name+ro+'_Tomo_'+t+'_'+'C'+carbon
                #print search_string
                folder_dict[t].append(toolbox.latest_data(contains = search_string,older_than = older_than,raise_exc = False))
    ### two carbons were involved
    elif len(carbon) == 2:
        for ro in ['positive','negative']:
            for t in ['XX','YY','XY','YX']:
                search_string = folder_name+ro+'_state'+logicstate+'_Tomo_'+t+'_'+'C'+carbon
                folder_dict[t].append(toolbox.latest_data(contains = search_string,older_than = older_than,raise_exc = False))

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil8'

    if len(carbon) == 1:
        #print folder_dict
        x_labels,npX,npY,npX_u,npY_u, seq_length = get_dephasing_data(folder_dict,ssro_calib_folder, do_get_sequence_length=True)

        folder_dict['res'] = np.sqrt(npY**2+npX**2)
        folder_dict['res_u'] = np.sqrt((npX*npX_u)**2+(npY*npY_u)**2)/np.sqrt((npX**2+npY**2))
        folder_dict['sweep_pts'] = x_labels


    elif len(carbon) == 2:

        tomo_dict, seq_length = get_dephasing_data(folder_dict,ssro_calib_folder, tomos = ['XX','YY','XY','YX'], do_get_sequence_length=True)
        x_labels = tomo_dict['sweep_pts']
        XX,XX_u = np.array(tomo_dict['XX']),np.array(tomo_dict['XX_u'])
        YY,YY_u = np.array(tomo_dict['YY']),np.array(tomo_dict['YY_u'])
        XY,XY_u = np.array(tomo_dict['XY']),np.array(tomo_dict['XY_u'])
        YX,YX_u = np.array(tomo_dict['YX']),np.array(tomo_dict['YX_u'])

        folder_dict['res'] = np.sqrt(XX**2+YY**2+XY**2+YX**2)/np.sqrt(2) ### geometric sum of possible correlations divided by maximum --> sqrt(2)
        ### calculate error bar in the same fashion
        folder_dict['res_u'] = (np.sqrt((XX*XX_u)**2+(YY*YY_u)**2+(YX*YX_u)**2+(XY*XY_u)**2)/np.sqrt(XX**2+YY**2+XY**2+YX**2))/np.sqrt(2)
        folder_dict['sweep_pts'] = x_labels

    if do_T2correct:
        folder_dict = do_T2_correction(
                carbon=carbon, folder_dict=folder_dict, sequence_duration_us=seq_length*10**6)

    return folder_dict

def do_T2_correction(carbon='2', folder_dict=[], sequence_duration_us=10):
    #print 'correcting for T2* decay before fitting. Are you sure this is what you want?'
    T2star_us = {'2': 12600, '1': 10150, '5': 18550, '3': 6350, '6': 4150 }

    n_of_reps = np.array(folder_dict['sweep_pts'])
    bloch_vec_len = np.array(folder_dict['res'])
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
        T2_Factors[count] = np.exp(-((float(n_of_reps[count])*sequence_duration_us)**2)/(2*T2star_us[carbon]**2))
        if T2_Factors[count] < 0.25:
            # print count
            cut_index = count
            # print T2_Factors[count]
            break
            
    #### cut away unreliable datapoints due to T2*
    bloch_vec_len = bloch_vec_len[:cut_index]
    T2_Factors = T2_Factors[:cut_index]
    n_of_reps = n_of_reps[:cut_index]
    y_u = y_u[:cut_index]

    if 0. in T2_Factors:
        print 'Warning: devision by zero would be required. I take uncorrected data '
        return folder_dict
    else:
        folder_dict['sweep_pts'] = n_of_reps
        folder_dict['res'] = np.divide( bloch_vec_len, T2_Factors)
        folder_dict['res_u'] = y_u
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
    plt.errorbar(folder_dict['sweep_pts'],folder_dict['res'], folder_dict['res_u'], marker='o', label='C'+carbon)
    plot.plot_fit1d(fit_result, np.linspace(folder_dict['sweep_pts'][0],folder_dict['sweep_pts'][-1],1001),ax=ax, plot_data=False,add_txt=False, lw = 2)

    plt.xlabel(plot_Xlabel)
    plt.ylabel('Bloch vector length')
    plt.title('Dephasing for C'+carbon+' '+ get_tstamp_from_folder(folder_dict['X'][0]))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(os.path.join(folder_dict['X'][0],'CarbonDephasing.pdf'),format='pdf')
    plt.savefig(os.path.join(folder_dict['X'][0],'CarbonDephasing.png'),format='png')
    plt.show()
    plt.close('all')

def extract_coupling_strength(folder_dict):
     ### after fitting: get the coupling strength. (one could make this a subroutine at some point.)
        # for the coupling strength we need logic state ms = 0 and ms = +-1 frequencies.
        coupling = 0
        #extract one complete path
        for key in ['XX','X']:
            if folder_dict[key] != []:
                folder = folder_dict[key][0]
            # else:
            #     print 'Warning'

        carbon_list,logic_state = get_from_hdf5(folder,['carbon_list','2qb_logical_state'])
        ### extract the expected coupling strength (works for 1 and 2 carbons)
        for ii,c in enumerate(carbon_list):
            if ii == 0:
                [coupling] = get_from_hdf5(folder,['C'+str(c)+'_freq_1_m1'])
                f_m1,f_0 = get_from_hdf5(folder,['C'+str(c)+'_freq_1_m1','C'+str(c)+'_freq_0'])
            else:
                f_m1,f_0 = get_from_hdf5(folder,['C'+str(c)+'_freq_1_m1','C'+str(c)+'_freq_0'])
                if logic_state == 'X':
                    # print 'carbon freqs',f_m1,f_0
                    coupling += (f_m1-f_0)
                else: 
                    # print 'mX'
                    # print 'carbon freqs',f_m1,f_0
                    coupling +=  - f_m1-f_0

        #calculate Delta_f w.r.t. to f_0

        coupling = np.abs(np.abs(coupling)-f_0)
        print 'Coupling ', coupling

        return folder_dict, coupling, folder

def Sweep_repetitions(older_than = None,
        folder_name ='Repetitions_',
        carbon = '2',
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
        offset = folder_dict['sweep_pts'][np.argmax(folder_dict['res'])]
    print 'A0 ', A0, ' offset ', offset
    
    if fitGauss:
        #a() + A() * np.exp(-(x-x0())**2/(2*sigma()**2))
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
        #print fit_result
        return coupling, fit_result['params_dict']['tau'], fit_result['error_dict']['tau'], folder, fit_result
    else:
        print fit_result['params_dict'], fit_result['error_dict']

def Sweep_Rep_List( carbons = ['1','2'],
        older_than = None, 
        do_T2correct = False,
        folder_name = 'Repetitions_', 
        ssro_calib_timestamp = None,**kw):

    ## other key word arguments
    fit_results = kw.pop('fit_results',True)
    sequence_length = kw.pop('sequence_length',None)
    logicstate_list = kw.pop('logicstate_list',len(carbons)*['X']) ## can be list such as ['X','mX'], used for DFS measurements.

    x_arr = []
    y_arr = []
    y_u_arr = []
    for c,logicstate in zip(carbons,logicstate_list):

        folder_dict= extract_data_from_sweep(older_than = older_than,
            folder_name =folder_name, carbon = c,
            ssro_calib_timestamp =ssro_calib_timestamp,
            logicstate = logicstate,
            do_T2correct=do_T2correct, **kw)
        x_arr.append(folder_dict['sweep_pts'])
        y_arr.append(folder_dict['res'])
        y_u_arr.append(folder_dict['res_u'])
    ### convert to time instead of repetitions:
    if sequence_length != None:
        x_arr = [x*sequence_length for x in x_arr]
    
    fig = plt.figure()
    ax = plt.subplot()
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

        if fit_results:
            A0 = y[0]
            offset = 0
            decay = 50
            if sequence_length != None:
                decay = decay*sequence_length
            x0 = 0
            p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)
            # p0,fitfunc,fitfunc_str = common.fit_gauss(offset,A0,x0,decay)
            fixed = [0,3]

            fit_result = fit.fit1d(x,y,None,p0 = p0, fitfunc = fitfunc, do_print = True, ret = True, fixed = fixed)
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax,color = color_list[jj], plot_data=False,add_txt = False, lw = 2)

        label_txt = 'C'+carbon
        if len(carbon)!=1:
            label_txt = label_txt+'_'+logicstate

        plt.errorbar(x,y,y_u,marker='o',color = color_list[jj],label=label_txt)


    plt.xlabel('Number of repetitions')

    if sequence_length != None:
        plt.xlabel('elapsed time (us)')
    elif folder_name == 'Memory_Sweep_repump_time_':
        plt.xlabel('average repump time (us)')

    plt.ylabel('Bloch vector length')
    plt.title(get_tstamp_from_folder(folder) + ' Dephasing for C'+carbon)
    plt.legend()#bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(os.path.join(folder,'CarbonDephasing.pdf'),format='pdf')
    plt.savefig(os.path.join(folder,'CarbonDephasing.png'),format='png')
    plt.show()
    plt.close('all')

def coupling_vs_repetitions(c_identifiers,**kw):

    is_X_measurement = kw.get('is_X_measurement', True)
    older_than = kw.get('older_than',None)
    s = 0
    for c in c_identifiers:
        s += len(c)
    x = np.zeros(s)
    y = np.zeros(s)
    y_u = np.zeros(s)

    ## acquire data
    ii = 0
    for carbon in c_identifiers:
        if len(carbon) > 1:
            for logicstate in ['X','mX']:
                x[ii],y[ii],y_u[ii],folder,fit_result = Sweep_repetitions(carbon = carbon,logicstate = logicstate,return_fits=True,plot_result = False,**kw)
                ii +=1
        else:
            x[ii],y[ii],y_u[ii],folder,fit_result = Sweep_repetitions(carbon = carbon,return_fits=True,plot_result = False,**kw)
            ii +=1
        if fit_result['reduced_chisq']>0.005:
            print 'Bad fit in ', folder, 'Xi_square is ', fit_result['reduced_chisq']

    return x,y,y_u,folder

def repump_power_vs_repetitions(c_identifier, repump_powers=[0], **kw):

    older_than = kw.pop('older_than',None)
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

    print folder
    return x,y,y_u,folder
    

def Osci_period(carbon = '1',older_than = None,ssro_calib_timestamp = None, do_print=False, add_txt=True, **kw):

    fit_results = kw.pop('fit_results',True)
    folder_name = kw.pop('folder_name','Memory_NoOf_Repetitions_')

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

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil8'
        # print ssro_calib_folder

        ### extract data
    x_labels,npX,npY,npX_u,npY_u = get_dephasing_data(folder_dict,ssro_calib_folder)

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
                print 'i was here'
                print decay
                ax.plot(np.linspace(x_labels[0],x_labels[-1],201), fitfunc(np.linspace(x_labels[0],x_labels[-1],201)), ':', lw=2)
        
        plt.errorbar(x_labels,y,y_u,marker='o',color = color_list[jj],label='C'+carbon+['X','Y'][jj])

    ## define folder for data saving
    folder = folder_dict[t][0]

    plt.xlabel('Repump repetitions')
    plt.ylabel('Contrast')
    plt.title(get_tstamp_from_folder(folder) + ' Dephasing for C'+carbon)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(os.path.join(folder,'CarbonDephasing_osci.pdf'),format='pdf')
    plt.savefig(os.path.join(folder,'CarbonDephasing_osci.png'),format='png')
    if not auto_analysis:
        plt.show()
    plt.close('all')

    print 'Results are saved in ', folder[18:18+15]
    if auto_analysis:
        return fit_result

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
        data_dict['folders'].append(toolbox.latest_data(contains = search_string,older_than = older_than,raise_exc = False))


    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil8'
        # print ssro_calib_folder
    for i,f in enumerate(data_dict['folders']):
        a = mbi.MBIAnalysis(f)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata',CR_after_check = CR_after_check)
        a.get_electron_ROC(ssro_calib_folder)



        x_labels = a.sweep_pts.reshape(-1)
        if i == 0:
            data_dict['res'] = ((a.p0.reshape(-1))-0.5)*2
            data_dict['res_u'] = 2*a.u_p0.reshape(-1)
        else:
            y = ((a.p0.reshape(-1))-0.5)*2
            y_u = 2*a.u_p0.reshape(-1)
            data_dict['res'] = [y0/2-y[ii]/2 for ii,y0 in enumerate(data_dict['res'])]
            data_dict['res_u'] = [np.sqrt(y0**2+y_u[ii]**2)/2 for ii,y0 in enumerate(data_dict['res_u'])]

    return x_labels,data_dict['res'],data_dict['res_u'],data_dict['folders'][0]

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
        print x2
        if 1 in x2 or 0 in x2 or 100 in x2:
            print x2
            print f
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

def Z_decay_vs_perp_coupling(c_idents,**kw):

    ### perpendicular hyperfine copulings in kHz (from fingerprint fits).
    hyperfine_params = {}
    hyperfine_params['C1']  = 25.0
    hyperfine_params['C2']  = 43.0
    hyperfine_params['C3']  = 55.0
    hyperfine_params['C4']  = 21.0
    hyperfine_params['C5']  = 26.0
    hyperfine_params['C6']  = 12.0



    label = kw.pop('label','')
    older_than0 = kw.pop('older_than',None)
    return_vals = kw.pop('return_vals',False)

    x,y1,y_u = np.array([]),np.array([]),np.array([])
    coupling_minus = np.array([])
    coupling_plus =np.array([])
    fitted_decays = np.array([])
    fitted_decays_u = np.array([])
    loop_bit = True

    ### coolect data from different measurements and fuse them. (needed for ZZ decays)
    for Cs in c_idents:
        print Cs
        for state in ['mX','X']: ### need this loop for DFS configurations
            if len(Cs) == 2:
                folder_name = '_state'+state+'_Tomo_ZZ_C'+Cs
            else:
                folder_name = '_Tomo_Z_C'+Cs


            ### reset vlaues before new round
            loop_bit = True
            older_than = older_than0
            x,y1,y_u = np.array([]),np.array([]),np.array([])

            while loop_bit:

                x2,y2,y_u2,f = get_PosNeg_data(folder_name,label = label,older_than = older_than, **kw)

                # print f
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
                if 1 in x2 or 0 in x2:

                    loop_bit = False

                    # extract perpendicular coupling strength.
                    if len(Cs) == 2:
                        

                        [C1_perp,C2_perp] = [hyperfine_params['C'+Cs[0]],hyperfine_params['C'+Cs[1]]]
                        cplus = np.sqrt(C1_perp**2 + C2_perp**2)/np.sqrt(2.)
                        cminus = np.sqrt(abs(C1_perp**2 - C2_perp**2))/np.sqrt(2.)
                        coupling_minus = np.append(coupling_minus,[cminus])
                        coupling_plus = np.append(coupling_plus,[cplus])
                    else:
                        C1_perp = hyperfine_params['C'+Cs[0]]
                        coupling_minus = np.append(coupling_minus,[C1_perp])
                        coupling_plus = np.append(coupling_plus,[C1_perp])
                # print coupling_minus

            ### sort arrays after adding them up
            y1 = y1[np.argsort(x)]
            y_u = y_u[np.argsort(x)]
            x = np.sort(x)
            ### fit exponential to the data

            offset,A0,decay,x0 = 0,np.amax(y1),400,0
            p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)
            fixed = [0,3]
            fit_result = fit.fit1d(x,y1,None,p0 = p0, fitfunc = fitfunc, do_print = False, ret = True, fixed = fixed)
            # print fit_result['params_dict'].keys()
            fitted_decays = np.append(fitted_decays,[fit_result['params_dict']['tau']])
            fitted_decays_u = np.append(fitted_decays_u,[fit_result['error_dict']['tau']])

            # print fitted_decays
            ### break the state loop for a single carbon in roder to not take the same dat twice.
            if len(Cs) == 1:
                break
    ### plot results:
    if return_vals:
        return fitted_decays,fitted_decays_u

    fig = plt.figure()
    ax = plt.subplot()

    plt.xlabel('Perpendicular hyperfine (kHz)')
    plt.ylabel('Fitted Z decay constant')

    if True: #Log-Plot
        plt.ylim(50,10000)
        ax.set_yscale("log", nonposy='clip')
    else:
        plt.ylim(0,3000)
    plt.errorbar(coupling_plus, fitted_decays, yerr = fitted_decays_u,fmt = 'o',color = 'b', label = 'sum')
    plt.errorbar(coupling_minus, fitted_decays, yerr = fitted_decays_u,fmt = 'o',color = 'g', label = 'difference')
    plt.legend()
    print 'saving to ', f
    plt.savefig(os.path.join(f,'perp_coupling_vs_repetitions.pdf'),format='pdf')
    plt.savefig(os.path.join(f,'perp_coupling_vs_repeititons.png'),format='png')
    plt.show()
    plt.close('all')

def repump_speed_doubleExp(timestamp=None, measurement_name = 'adwindata', ssro_calib_timestamp =None,
            exclude_first_n_points = 0.,
            offset = 0., x0 = 0, older_than= None,
            amplitude = 1., x_offs=0,   
            decay_constant_one = 0.2, decay_constant_two = 0.6, 
            plot_results = True, do_T2correct=False, log_plot=True,
            plot_fit = True, do_print = False, fixed = [2], show_guess = True, **kw):
   
    ''' 
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''
    figsize=(6,4.7)
    fitted_tau, fitted_tau2 = [],[] 
    fitted_tau_err, fitted_tau2_err = [],[]
    fit_results = []

    nf_powers = kw.get('nf_powers', [-1])
    repumper_powers = kw.get('repumper_powers', [-1])
    invert = kw.get('invert', 'check_fn')
    readout_bases = kw.get('ro_bases', [0,-1])
    
    if timestamp != None:  #evaluate particular file
        folder = toolbox.data_from_time(timestamp)
        nf_powers=[-1]

    if ssro_calib_timestamp == None:
        ssro_calib_folder = toolbox.latest_data('SSRO', older_than=older_than)
        print 'Using SSRO timestamp ', ssro_calib_folder
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print 'Using SSRO timestamp ', ssro_calib_folder

    for sweep_elem in range(len(nf_powers)):
        fig = plt.figure()
        ax = plt.subplot()
        plt.xlabel('time (ns)')
        plt.ylabel('p(m$_s$=$\pm1$)')
        folderlist = [] 
        for readout_elem in readout_bases:  
            if timestamp != None:
                folderlist = [toolbox.data_from_time(timestamp)]
            elif nf_powers == [-1]:
                folderlist = [toolbox.latest_data('ElectronRepump', older_than=older_than)]
            else:
                if readout_elem==-1:
                    add_to_msmt_name = 'm1RO_'
                elif readout_elem==1:
                    add_to_msmt_name = 'p1RO_'
                else:
                    add_to_msmt_name = '0RO_'
                searchstring = 'ElectronRepump_'+str(nf_powers[sweep_elem]*1e9)+'nW_'+str(repumper_powers[sweep_elem]*1e9)+'nW_'+add_to_msmt_name
                try:
                    new_folder = toolbox.latest_data(searchstring, older_than=older_than)
                    folderlist.append(new_folder)
                except:
                    print 'No file found for readout in ', readout_elem
        print folderlist

        for folder in folderlist:
            print 'folder is ', folder
            plt.title(folder) 
            a = mbi.MBIAnalysis(folder)
            a.get_sweep_pts()
            CR_after_check = None
            a.get_readout_results(name='adwindata',CR_after_check = CR_after_check)
            a.get_electron_ROC(ssro_calib_folder)
            x = 1000* a.sweep_pts.reshape(-1)[exclude_first_n_points:]

            if invert == 'check_fn':
                do_invert = True if '0RO' in folder else False 
            else:
                do_invert = invert        
            if do_invert:
                print 'inverting ', folder 
                y = np.array(1.) - a.p0.reshape(-1)[exclude_first_n_points:]
            else: 
                y = a.p0.reshape(-1)[exclude_first_n_points:]
            y_u = a.u_p0.reshape(-1)[exclude_first_n_points:]

            if log_plot:
                 # 0])+str(c_list[1]+'_repump_power'+str(sweep      plt.ylim(0.005,1.05)
                ax.set_yscale("log", nonposy='clip')
                plt.ylim(0.0001,1.05)
                plt.xlim(-10,np.amax(x))
            else:
                plt.ylim(0.0,1.05)
                plt.xlim(-10,np.amax(x))

            plot_color = 'r' if '0RO' in folder else 'b' if 'p1RO' in folder else 'k'
            plt.errorbar(x,y, yerr = y_u, fmt = 'o', color = plot_color , label = 'x')
            #plt.legend(bbox_to_anchor=(1.05, 1), loc=3, borderaxespad=0.)
            #fitfunction: y(x) = A * exp(-x/tau)+ A2 * exp(-x/tau2) + a
            p0, fitfunc, fitfunc_str = common.fit_repumping( offset, amplitude, decay_constant_one,
                decay_constant_two, x_offs )

            if plot_results and show_guess:
                ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True, fixed=fixed)

            ## plot data and fit as function of total time
            if plot_fit == True:
                plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), color = plot_color, ax=ax, plot_data=False)
            
            try:
                fitted_tau.append(fit_result['params_dict']['tau']) 
                fitted_tau2.append(fit_result['params_dict']['tau2']) 
                fitted_tau_err.append(fit_result['error_dict']['tau'])
                fitted_tau2_err.append(fit_result['error_dict']['tau2'])
            except:
                print 'fit didnt succeed'

            if plot_results:
                plt.savefig(os.path.join(folder, 'analyzed_result.pdf'), format='pdf')
                plt.savefig(os.path.join(folder, 'analyzed_result.png'), format='png')
    return fitted_tau, fitted_tau2, fitted_tau_err, fitted_tau2_err

def repump_speed_paper_plot(timestamp=None, measurement_name = 'adwindata', ssro_calib_timestamp =None,
            exclude_first_n_points = 0.,
            offset = 0., x0 = 0, older_than= None, newer_than=0, binwidth_ns = None,
            amplitude = 0.8, decay_constant_one = 0.2, decay_constant_two = 0.6, x_offs = 0,
            plot_results = True, do_T2correct=False,
            plot_fit = True, do_print = False, fixed = [2], show_guess = True, invert = [True],
            powers = [0], colors=['b']):
   
    ''' 
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''
    figwidthPRL=3.+3./8.
    golden_ratio = 1.62
    figsize=(figwidthPRL,figwidthPRL/golden_ratio)

    fitted_tau, fitted_tau2 = [],[]
    fitted_tau_err, fitted_tau2_err = [],[]

    #p0, fitfunc, fitfunc_str = [], [], []
    fig = plt.figure()
    ax = plt.subplot()
    plt.xlabel('time (ns)')
    plt.ylabel('p(m$_s$=$\pm1$)')
    plt.ylim(0.01,1.05)
    ax.set_yscale("log", nonposy='clip')
    plt.xlim(0.,2500)
    

    for count in np.arange(len(older_than)):
        x, y, y_u = [], [], []
        print 'older than ', older_than[count], ' and newer_than ', newer_than[count]
        fit_results = []
        #folder = toolbox.data_from_time(timestamp[count])
        folder_list = toolbox.latest_data('ElectronRepump', older_than=older_than[count], newer_than=newer_than[count], return_all=True)
        folder_list_ext = []
        for av_elem in folder_list[2]:
            folder_list_ext.append(toolbox.latest_data(av_elem))
        print 'folder is ', folder_list_ext
        
        folder = folder_list_ext[0]

        if ssro_calib_timestamp == None :
            ssro_calib_folder = toolbox.latest_data('SSRO', older_than=older_than[count])
        else:
            ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp[count])
            ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print 'Using SSRO timestamp ', ssro_calib_folder
            
        for elem in np.arange(len(folder_list_ext)):
            #print folder_list_ext[elem]
            a = mbi.MBIAnalysis(folder_list_ext[elem])
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata', CR_after_check = CR_after_check)
            a.get_electron_ROC(ssro_calib_folder)

            x_list = 1000*a.sweep_pts.reshape(-1)[exclude_first_n_points[count]:]
            x = np.append(x, x_list - x_list[0])  #shifts the graph such that it starts at t=0
            if invert[count]:
                y = np.append(y,1-a.p0.reshape(-1)[exclude_first_n_points[count]:])
            else:
                y = np.append(y,a.p0.reshape(-1)[exclude_first_n_points[count]:])
            y_u = np.append(y_u,a.u_p0.reshape(-1)[exclude_first_n_points[count]:])
            #print 'lengths are: x ', len(x), ' y ', len(y), ' y_u ', len(y_u)

        sortedx = np.argsort(x)
        y = y[sortedx]
        y_u = y_u[sortedx]
        x = x[sortedx]

        binned_x, binned_y, binned_yu, temp_y, temp_yu= [],[],[],[],[]

        if binwidth_ns[count] != None:
            last_x = 0
            for x_count in np.arange(len(x)):    
                if np.floor(x[x_count] / binwidth_ns[count]) > last_x and len(temp_y)>0:
                    binned_x.append(last_x * binwidth_ns[count])
                    binned_y.append(np.sum(temp_y)/len(temp_y))
                    binned_yu.append(np.sqrt(np.sum(np.square(temp_yu)))/len(temp_yu))
                    last_x += 1
                    temp_y, temp_yu = [],[]
                else: 
                    temp_y.append(y[x_count])
                    temp_yu.append(y_u[x_count])
        else:
            binned_x = x
            binned_y = y
            binned_yu = y_u

        #plt.errorbar(a.sweep_pts[exclude_first_n_points[count]:], 1-a.p0[exclude_first_n_points[count]:,0], yerr = a.u_p0[exclude_first_n_points[count]:,0], fmt = 'o',color = colors[count], label = '')
        #print x, y, y_u
        plt.errorbar(binned_x, binned_y, yerr = binned_yu, fmt = 'o', color = colors[count], label = '')

        #fitfunction: y(x) = A * exp(-x/tau)+ A2 * exp(-x/tau2) + a
        p0, fitfunc, fitfunc_str = common.fit_repumping( 
            offset[count], amplitude[count], decay_constant_one[count],
            decay_constant_two[count], x_offs[count] )
        # p0, fitfunc, fitfunc_str = common.fit_double_exp_decay_with_x_offset(np.array(offset[count]), np.array(amplitude[count]), 
        #          np.array(decay_constant_one[count]), np.array(decay_constant_two[count]), np.array(x_offs[count]))

        if plot_results and show_guess:
            ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

        fit_result = fit.fit1d( x, y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True, fixed=fixed[count])

        ## plot data and fit as function of total time
        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001),color=colors[count],log=True, 
                ax=ax, plot_data=False, legend=None, add_txt=False)

        fit_results.append(fit_result['params_dict']['tau'])
        #print fit_result['params_dict'], fit_result['error_dict']

        #fitted_tau.append(fit_result['params_dict']['tau'])
        #fitted_tau_err.append(fit_result['error_dict']['tau'])

    if plot_results:
        save_figure_to = 'D:\measuring\QMem_plots'
        #save_figure_to = 'K:\ns\qt\Diamond\Eigenpapers\15-WeaklyCoupledQuantumMemory\Figures'
        print 'saving to: ', save_figure_to
        plt.savefig(os.path.join(save_figure_to, 'Fig2.pdf'), format='pdf')
        plt.savefig(os.path.join(save_figure_to, 'Fig2.png'), format='png')

    return fit_results

# def analyze_avg_repump_time(carbons = ['2','3'],folder_name = 'Memory_Sweep_repump_time_',fit_results = False,older_than = older_than):
#     CD.Sweep_Rep_List(carbons = carbons, folder_name = folder_name, fit_results = False, older_than = older_than)

