import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from analysis.lib.QEC import ConditionalParity as CP
from analysis.lib.fitting import fit, common, ramsey;reload(common); reload(fit)

reload (CP)
import h5py
import csv

RO_corr_1qb = 1.
RO_corr_3qb = 1.


from matplotlib import pyplot as plt
script_name = 'three_qubit_QEC_analysis.py'

''' These functions are old but can be useful '''

def Contrast_Plot_QEC_full(timestamp = None, measurement_name = ['adwindata'],folder_name ='RO1',
        ssro_calib_timestamp =None, save = True,
        do_plot  = True):

    ''' this function is currently not used anymore '''

        ### SSRO calibration

    for k in range(3):
        print k 
        timestamp_pos, folder_a = toolbox.latest_data(contains = 'positive_'+folder_name+ '_k'+ str(k), older_than = timestamp,return_timestamp = True)
        timestamp_neg, folder_b = toolbox.latest_data(contains = 'negative_'+folder_name+ '_k'+ str(k), older_than = timestamp,return_timestamp = True)
        
        x_t, y_t, y_err_t, y_00_t, y_00_err_t, p00_avg_t, y_01_t, y_01_err_t, p01_avg_t, y_10_t, y_10_err_t, p10_avg_t, y_11_t, y_11_err_t, p11_avg_t  = Contrast_Plot_QEC(timestamps=[timestamp_pos, timestamp_neg], 
                                measurement_name = ['adwindata'],folder_name =folder_name,
                                post_select_QEC = True, ssro_calib_timestamp =ssro_calib_timestamp, do_plot = False, return_data = True)
        if k == 0:
            x = list(x_t) 
            y = list(y_t)
            y_err = list(y_err_t)
            y_00 = list(y_00_t) 
            y_00_err = list(y_00_err_t) 
            p00_avg = list(p00_avg_t) 
            y_01 = list(y_01_t) 
            y_01_err = list(y_01_err_t) 
            p01_avg = list(p01_avg_t) 
            y_10 = list(y_10_t) 
            y_10_err = list(y_10_err_t) 
            p10_avg = list(p10_avg_t) 
            y_11 = list(y_11_t) 
            y_11_err = list(y_11_err_t) 
            p11_avg = list(p11_avg_t) 
        else:
            x.extend(list(x_t) )
            y.extend(list(y_t))
            y_err.extend(list(y_err_t))
            y_00.extend(list(y_00_t) )
            y_00_err.extend(list(y_00_err_t) )
            p00_avg.extend(list(p00_avg_t) )
            y_01.extend(list(y_01_t) )
            y_01_err.extend(list(y_01_err_t) )
            p01_avg.extend(list(p01_avg_t) )
            y_10.extend(list(y_10_t) )
            y_10_err.extend(list(y_10_err_t) )
            p10_avg.extend(list(p10_avg_t) )
            y_11.extend(list(y_11_t) )
            y_11_err.extend(list(y_11_err_t) )
            p11_avg.extend(list(p11_avg_t) )


    if do_plot == True:
        fig,ax = plt.subplots() 
        rects = ax.errorbar(x,y,yerr=y_err,color = 'k' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title(str(folder_a)+'/'+'\n' + script_name)
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')

        try:
            fig.savefig(
                os.path.join(folder_a,'QEC_full.png'))
        except:
            print 'Figure has not been saved.'

        fig,ax = plt.subplots() 
        rects = ax.errorbar(x,y_00,yerr=y_00_err,color = 'k' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title(str(folder_a)+'/'+'\n postselect_00')
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')

        try:
            fig.savefig(
                os.path.join(folder_a,'QEC_ps_00_full.png'))
        except:
            print 'Figure has not been saved.'

        fig,ax = plt.subplots() 
        rects = ax.errorbar(x,y_01,yerr=y_01_err,color = 'k' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title(str(folder_a)+'/'+ '\n postselect_01')
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')

        try:
            fig.savefig(
                os.path.join(folder_a,'QEC_ps_01_full.png'))
        except:
            print 'Figure has not been saved.'

        fig,ax = plt.subplots() 
        rects = ax.errorbar(x,y_10,yerr=y_10_err,color = 'k' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)

        ax.set_title(str(folder_a)+'/'+'\n postselect_10')
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')

        try:
            fig.savefig(
                os.path.join(folder_a,'QEC_ps_10_full.png'))
        except:
            print 'Figure has not been saved.'

        fig,ax = plt.subplots() 
        rects = ax.errorbar(x,y_11,yerr=y_11_err,color = 'k' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title(str(folder_a)+'/'+ '\n postselect_11')
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')

        try:
            fig.savefig(
                os.path.join(folder_a,'QEC_ps_11_full.png'))
        except:
            print 'Figure has not been saved.'


        fig,ax = plt.subplots()
        ax.set_title(str(folder_a)+'/'+ '\n probabilities')
        ax.plot(x,p00_avg, 'c', label = 'p00')
        ax.plot(x,p01_avg, 'k', label = 'p01')
        ax.plot(x,p10_avg, 'm', label = 'p10')
        ax.plot(x,p11_avg, 'b', label = 'p11')
        plt.legend()
        ax.set_xlabel('error probability')
        ax.set_ylabel('outcome probability')

        try:
            fig.savefig(
                os.path.join(folder_a,'QEC_probabilities_full.png'))
        except:
            print 'Figure has not been saved.'

def Plot_errorcurve_no_QEC(timestamp = None, measurement_name = ['adwindata'],folder_name ='QEC',
        ssro_calib_timestamp =None, save = True,
        plot_fit = True, return_data = False) :
    ''' Currently not used '''

    ### SSRO calibration
    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'

    ### Obtain and analyze data
        ### postive RO data
    timestamp_pos, folder_a = toolbox.latest_data(contains = 'positive_'+folder_name, older_than = timestamp,return_timestamp = True)
    timestamp_neg, folder_b = toolbox.latest_data(contains = 'negative_'+folder_name, older_than = timestamp,return_timestamp = True)
    # print folder_a
    # print folder_b
    
    a = mbi.MBIAnalysis(folder_a)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder)
    y_a= ((a.p0.reshape(-1)[:])-0.5)*2
    y_err_a = 2*a.u_p0.reshape(-1)[:] 

    b = mbi.MBIAnalysis(folder_b)
    b.get_sweep_pts()
    b.get_readout_results(name='adwindata')
    b.get_electron_ROC(ssro_calib_folder)
    y_b= ((b.p0.reshape(-1)[:])-0.5)*2
    y_err_b = 2*b.u_p0.reshape(-1)[:] 

    x = a.sweep_pts.reshape(-1)[:]
    # x = range(len(y_a)) 


    
    ### Combine data
    y = (y_a - y_b)/2.
    y_err =  1./2*(y_err_a**2 + y_err_b**2)**0.5 
    


    if plot_fit ==True: 
        fig,ax = plt.subplots() 
        ax.errorbar(x,y,yerr=y_err,color = 'k' )
        ax.set_ylim(-1.1,1.1)
        ax.set_title(str(folder_a)+'/'+str(timestamp_pos))
        ax.set_xticks(x)
        ax.set_xlim([-0.1,1.1])
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')


    if save and ax != None:
        try:
            fig.savefig(
                os.path.join(folder_a,'QEC.png'))
        except:
            print 'Figure has not been saved.'

    if return_data == True:
        return x, y, y_err

''' Used basic functions '''

def load_QEC_data(folder, ssro_calib_folder, post_select = True):
    ''' Loads a QEC measurment and returns all 
    the results in a dictionairy'''

    a = CP.ConditionalParityAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata', post_select_QEC = False)
    # print ssro_calib_folder
    a.get_electron_ROC(ssro_calib_folder)

    x = a.sweep_pts.reshape(-1)
    c0, c0_u = a.convert_fidelity_to_contrast(a.p0,a.u_p0)
  
    if post_select:
        a = CP.ConditionalParityAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata', post_select_QEC = True)
        a.get_electron_ROC(ssro_calib_folder  , post_select_QEC = True)

        c0_00,c0_00_u =  a.convert_fidelity_to_contrast(a.p0_00,a.u_p0_00)
        c0_01,c0_01_u =  a.convert_fidelity_to_contrast(a.p0_01,a.u_p0_01)
        c0_10,c0_10_u =  a.convert_fidelity_to_contrast(a.p0_10,a.u_p0_10)
        c0_11,c0_11_u =  a.convert_fidelity_to_contrast(a.p0_11,a.u_p0_11)

    data_dict = {}
    # data_dict['a'] = a
    data_dict['x']          = x
    data_dict['c0']         = c0 
    data_dict['c0_u']       = c0_u 

    if post_select:
        data_dict['c0_00']      = c0_00 
        data_dict['c0_00_u']    = c0_00_u 
        data_dict['c0_01']      = c0_01 
        data_dict['c0_01_u']    = c0_01_u 
        data_dict['c0_10']      = c0_10 
        data_dict['c0_10_u']    = c0_10_u 
        data_dict['c0_11']      = c0_11 
        data_dict['c0_11_u']    = c0_11_u 
        data_dict['p00']        = a.p00
        data_dict['p01']        = a.p01
        data_dict['p10']        = a.p10
        data_dict['p11']        = a.p11

    return data_dict

def get_folder(timestamp, folder_name):
    if timestamp == None:
        timestamp, folder   = toolbox.latest_data(folder_name, return_timestamp =True)
    else:
        folder = toolbox.data_from_time(timestamp)
    return timestamp, folder    

def get_ssro_folder(ssro_calib_timestamp):
    if ssro_calib_timestamp == None:
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
    return ssro_calib_folder

''' These functions both load and plot individual measurements TODO THT, seperate loading from plotting'''

def plot_single_QEC_result(timestamps = [None], folder_name ='QEC', ssro_calib_timestamp = None,
        post_select = False, save = True, title = None, fontsize = 10) :
    '''
    Plots the results of a single QEC/Encoding measurement, from a raw data folder.
    post_select must be false if no parity measurements
    Length of the timestamp gives eaither a single measurement or a positive/negative one
    '''
    ### Timestamps and folders
    timestamp, folder = get_folder(timestamps[0], folder_name)
    if len(timestamps) == 2:
        timestamp2, folder2 = get_folder(timestamps[1], folder_name)
    ssro_calib_folder = get_ssro_folder(ssro_calib_timestamp)

    ### Get the data
    data = load_QEC_data(folder = folder, ssro_calib_folder = ssro_calib_folder, post_select=post_select)
    if len(timestamps) == 2:
        data2 = load_QEC_data(folder = folder2, ssro_calib_folder = ssro_calib_folder, post_select=post_select)

    ### Combine the data in case of positive/negative measurement    
    if len(timestamps) == 2:
        
        data['c0'] = (data['c0'] - data2['c0'])/2.
        data['c0_u'] =  1./2*(data['c0_u']**2 + data2['c0_u']**2)**0.5      

    ### Plots
    plt.rc('font', size=fontsize)
    
    fig,ax = plt.subplots()
    
    ax.errorbar(data['x'], data['c0'],yerr=data['c0_u'], color = 'b' )
    
    if title == None:
        ax.set_title(str(folder)+'/'+str(timestamp))
    else:
        ax.set_title(title)
    
    ax.set_ylim(-1,1)
    ax.set_xlim(-0.05,1.05)
    ax.hlines([-1,0,1],data['x'][0]-1,data['x'][-1]+1,linestyles='dotted')

    fig.savefig(os.path.join(folder,'QEC_single_measurement.png'))
    fig.savefig(os.path.join(folder, 'QEC_single_measurement.pdf'),
            format='pdf',bbox_inches='tight')


    ### Post selected data
    if post_select ==True:
        
        ### Combine negative/positive post selected data
        if len(timestamps) ==2:
            data['c0_00'] = (data['c0_00'] - data2['c0_00'])/2.
            data['c0_00_u'] =  1./2*(data['c0_00_u']**2 + data2['c0_00_u']**2)**0.5   
            data['c0_01'] = (data['c0_01'] - data2['c0_01'])/2.
            data['c0_01_u'] =  1./2*(data['c0_01_u']**2 + data2['c0_01_u']**2)**0.5
            data['c0_10'] = (data['c0_10'] - data2['c0_10'])/2.
            data['c0_10_u'] =  1./2*(data['c0_10_u']**2 + data2['c0_10_u']**2)**0.5
            data['c0_11'] = (data['c0_11'] - data2['c0_11'])/2.
            data['c0_11_u'] =  1./2*(data['c0_11_u']**2 + data2['c0_11_u']**2)**0.5

        ### Plot postselected data
        fig,ax = plt.subplots()
        ax.errorbar(data['x'], data['c0_00'], yerr=data['c0_00_u'], label = '00',color = 'k' )
        ax.errorbar(data['x'], data['c0_01'], yerr=data['c0_01_u'], label = '01',color = 'c' )
        ax.errorbar(data['x'], data['c0_10'], yerr=data['c0_10_u'], label = '10',color = 'g' )
        ax.errorbar(data['x'], data['c0_11'], yerr=data['c0_11_u'], label = '11',color = 'r' )
        ax.set_xlim(-0.2,1.2)
        ax.legend()
 
        if title == None:
            ax.set_title(str(folder)+'/'+str(timestamp))
        else:
            ax.set_title(title)
        ax.hlines([-1,0,1],data['x'][0]-1,data['x'][-1]+1,linestyles='dotted')

        fig.savefig(os.path.join(folder,'QEC_single_measurment_ps.png'))
        fig.savefig(os.path.join(folder,'QEC_single_measurment_ps.pdf'),
                format='pdf',bbox_inches='tight')

        
        ### Cobine neagtive/postitive outcome probabilities
        if len(timestamps) == 2:
            data['p00'] = (data['p00'] + data2['p00'])/2
            data['p01'] = (data['p01'] + data2['p01'])/2
            data['p10'] = (data['p10'] + data2['p10'])/2
            data['p11'] = (data['p11'] + data2['p11'])/2

        ### Outcome probabilities
        fig,ax = plt.subplots()
        ax.set_title(str(folder)+'/'+ '\n probabilities')
        ax.plot(data['x'],data['p00'], 'co', label = 'p00')
        ax.plot(data['x'],data['p01'], 'ko', label = 'p01')
        ax.plot(data['x'],data['p10'], 'mo', label = 'p10')
        ax.plot(data['x'],data['p11'], 'bo', label = 'p11')
        ax.plot(data['x'],data['p00']+data['p01']+ data['p10']+data['p11'], 'go', label = 'sum')
        plt.legend()
        ax.set_xlim(-0.2,1.2)
        ax.set_xlabel('error probability')
        ax.set_ylabel('outcome probability')  
        ax.set_title(str(folder)+'/'+str(timestamp) + '_QEC_probs')                

        print data['p00'] + data['p01'] + data['p10'] +data['p11']


        fig.savefig(os.path.join(folder,'QEC_probs'+'.png'))
       

def QEC_create_data_dict(older_than = None, RO = 0, state = 'Z', len_k = 6, sym = '11'):
    QEC_dict = {}
    k_dict = {}
    
    for error_sign in [1,-1]:
        # print 'sign_'+str(error_sign)
        QEC_dict[str(error_sign)] ={}
        for direction in ['positive','negative']:
            QEC_dict[str(error_sign)][direction] = {}

            for k in range(len_k):
                # print 'k_'+str(k)
                
                print '----'
                timestamp, folder = toolbox.latest_data(contains = sym +'_'+direction+'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_sign)
                                                        +'_'+state, older_than = older_than,return_timestamp = True)
                # print folder
                SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
                print SSRO_folder
                print '----'
                k_dict['k_'+str(k)] ={}
                k_dict['k_'+str(k)] = load_QEC_data(folder, SSRO_folder, post_select = True) 
                          
            for item in k_dict['k_0']:
                if len_k == 4:
                    QEC_dict[str(error_sign)][direction][item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item], k_dict['k_3'][item]), axis=0)
                elif len_k == 6:
                    QEC_dict[str(error_sign)][direction][item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item], k_dict['k_3'][item], k_dict['k_4'][item], k_dict['k_5'][item]), axis=0)

    return QEC_dict,folder


def no_QEC_create_data_dict(older_than = None, RO = 0, state = 'Z'):
    QEC_dict = {}
    k_dict = {}
    
    for error_sign in [1,-1]:
        # print 'sign_'+str(error_sign)
        QEC_dict[str(error_sign)] ={}
        for direction in ['positive','negative']:
            QEC_dict[str(error_sign)][direction] = {}

            for k in range(2):
                # print 'k_'+str(k)
                
                timestamp, folder = toolbox.latest_data(contains = 'no_corr_'+direction+'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_sign)+'_'+state, older_than = older_than,return_timestamp = True)

                SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
                # print SSRO_timestamp
                k_dict['k_'+str(k)] ={}
                k_dict['k_'+str(k)], folder = Plot_QEC(timestamp = timestamp, folder_name = folder,
                    ssro_calib_timestamp = SSRO_timestamp, return_raw = False, return_dict = True, post_select_QEC = False) 
                
                
            for item in k_dict['k_0']:
                QEC_dict[str(error_sign)][direction][item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item]), axis=0)

    return QEC_dict,folder


def QEC_create_data_dict_single_error_single_elRO(older_than = None, RO = 0, state = 'Z', len_k = 6, sym = '11',error_sign = 1, el_RO = 'positive',sweep_time = False):
    QEC_dict = {}
    k_dict = {}
    
    for k in range(len_k):
        # print 'k_'+str(k)
        if sweep_time == True:
            timestamp, folder = toolbox.latest_data(contains = sym +'_sweep_time_'+el_RO+'_RO'+str(RO)+'_k'+str(k)
                                                    +'_'+state, older_than = older_than,return_timestamp = True)
        elif sweep_time == False:
            timestamp, folder = toolbox.latest_data(contains = sym +'_'+el_RO+'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_sign)
                                                    +'_'+state, older_than = older_than,return_timestamp = True)
        SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
        print SSRO_folder
        k_dict['k_'+str(k)] ={}
        k_dict['k_'+str(k)] = load_QEC_data(folder, SSRO_folder, post_select = True) 
                  
    for item in k_dict['k_0']:
        if len_k == 4:
            QEC_dict[item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item], k_dict['k_3'][item]), axis=0)
        elif len_k == 6:
            QEC_dict[item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item], k_dict['k_3'][item], k_dict['k_4'][item], k_dict['k_5'][item]), axis=0)

    return QEC_dict,folder


def no_QEC_create_data_dict_single_error_single_elRO(idle = False,sweep_time = False, older_than = None, RO = 0, state = 'Z', error_sign = 1, el_RO = 'positive'):
    QEC_dict = {}
    k_dict = {}
    k_length = 3
    if sweep_time == True:
        k_length =4

    for k in range(k_length):
        # print 'k_'+str(k)
        if idle == False and sweep_time == False:
            timestamp, folder = toolbox.latest_data(contains = 'no_correct' +'_'+el_RO+'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_sign)
                                                    +'_'+state, older_than = older_than,return_timestamp = True)
        elif idle == True:
            timestamp, folder = toolbox.latest_data(contains = 'no_correct_idle' +'_'+el_RO+'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_sign)
                                                    +'_'+state, older_than = older_than,return_timestamp = True)
        elif sweep_time == True:
            print('no_correct_sweep_time' +'_'+el_RO+'_RO'+str(RO)+'_k'+str(k)
                                                    +'_'+state) == 'no_correct_sweep_time_negative_RO6_k5_mZ'
            timestamp, folder = toolbox.latest_data(contains = 'no_correct_sweep_time' +'_'+el_RO+'_RO'+str(RO)+'_k'+str(k)
                                                    +'_'+state, older_than = older_than,return_timestamp = True)
            print folder
        SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
        print SSRO_folder
        k_dict['k_'+str(k)] ={}
        k_dict['k_'+str(k)] = load_QEC_data(folder, SSRO_folder, post_select = False) 
                  
    for item in k_dict['k_0']:
        print item
        if k_length ==3:
            QEC_dict[item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item]), axis=0)
        if k_length ==6:
            QEC_dict[item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item],k_dict['k_3'][item],k_dict['k_4'][item],k_dict['k_5'][item]), axis=0)
        if k_length == 4:
            QEC_dict[item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item],k_dict['k_3'][item]), axis=0)


    return QEC_dict,folder


def single_Qubit_QEC_create_data_dict_single_error_single_elRO(older_than = None, sweep_time = False,Qubit = 1, state = 'Z', error_sign = 1, el_RO = 'positive'):
    QEC_dict = {}
    k_dict = {}
    print Qubit

    if sweep_time == False:
        len_k =2
    elif sweep_time == True:
        len_k = 4

    if Qubit == 1:
        carbon = 'C1'
    elif Qubit == 2:
        carbon = 'C5'
    elif Qubit == 3:
        carbon = 'C2'
    
    if state == 'X' or state == 'mX':
        RO = 2+(Qubit-1)*3
    elif state == 'Y' or state == 'mY':
        RO = 1+(Qubit-1)*3
    elif state == 'Z' or state == 'mZ': 
        RO = 0+(Qubit-1)*3
        if sweep_time == True:
            RO = Qubit-1

    for k in range(len_k+1):
        print 'k = '+str(k)
        print 'single_Qbit_sweep_time' +'_'+el_RO+'_'+ carbon +'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(-1)+'_'+state
        if sweep_time == False:
            timestamp, folder = toolbox.latest_data(contains = 'single_Qbit' +'_'+el_RO+'_'+ carbon +'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_sign)+'_'+state, older_than = older_than,return_timestamp = True)
        elif sweep_time == True and k !=4:
            print 'ok'
            timestamp, folder = toolbox.latest_data(contains = 'single_Qbit_sweep_time' +'_'+el_RO+'_'+ carbon +'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(-1)+'_'+state, older_than = '20150105_142000',return_timestamp = True)
        if sweep_time == True and  k == 4:
            timestamp, folder = toolbox.latest_data(contains = 'single_Qbit_sweep_time' +'_'+el_RO+'_'+ carbon +'_RO'+str(RO)+'_k'+str(0)+'_sign'+ str(-1)+'_'+state, older_than = '20150106_100000',return_timestamp = True)
            print 'k =4, timestamp = '+str(timestamp)
        SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
        print timestamp
        print SSRO_folder
        k_dict['k_'+str(k)] ={}
        k_dict['k_'+str(k)] = load_QEC_data(folder, SSRO_folder, post_select = False) 
                  
    for item in k_dict['k_0']:
        if len_k ==2:
            QEC_dict[item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item]), axis=0)
        elif len_k == 4:
            QEC_dict[item] = np.concatenate((k_dict['k_4'][item],k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item],k_dict['k_3'][item]), axis=0)
    return QEC_dict,folder

''' simple plotting QEC data without loading/saving '''

def QEC_data_single_state_RO(older_than = None,state = 'Z',RO = 0, sym = '00'):

    QEC_data_dict = {}
    u_list = ['c0_u', 'c0_00_u','c0_01_u','c0_10_u','c0_11_u']
    c_list = ['c0', 'c0_00','c0_01','c0_10','c0_11']
    p_list = ['p00','p01','p10','p11']
    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']


    QEC_dict, folder = QEC_create_data_dict(older_than = older_than, RO = RO, state = state, sym = sym)
    for v in range(5):
        QEC_data_dict[y_list[v]] = {}
        QEC_data_dict[y_err_list[v]] = {}


        QEC_data_dict[y_list[v]] = (QEC_dict[str(-1)]['positive'][c_list[v]]+
                                                        QEC_dict[str(1)]['positive'][c_list[v]]-
                                                        QEC_dict[str(-1)]['negative'][c_list[v]]-
                                                        QEC_dict[str(1)]['negative'][c_list[v]])/4       
        
        QEC_data_dict[y_err_list[v]] = ((QEC_dict[str(-1)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(1)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(-1)]['negative'][u_list[v]]**2+
                                                        QEC_dict[str(1)]['negative'][u_list[v]]**2)**0.5)/4
    
    QEC_data_dict['RO_contrast_pos_error'] =        (QEC_dict[str(1)]['positive']['c0']+
                                                    QEC_dict[str(1)]['negative']['c0'])/2     
    
    QEC_data_dict['RO_contrast_neg_error'] =        (QEC_dict[str(-1)]['positive']['c0']+
                                                    QEC_dict[str(-1)]['negative']['c0'])/2     

    QEC_data_dict['u_RO_contrast_pos_error'] =        (QEC_dict[str(1)]['positive']['c0_u']**2+
                                                    QEC_dict[str(1)]['negative']['c0_u']**2)**0.5/2     
    
    QEC_data_dict['u_RO_contrast_neg_error'] =        (QEC_dict[str(-1)]['positive']['c0_u']**2+
                                                    QEC_dict[str(-1)]['negative']['c0_u']**2)**0.5/2     

    for p in range(4):
            QEC_data_dict[p_list[p]] = {}
            
            QEC_data_dict[p_list[p]] =          (QEC_dict[str(-1)]['positive'][p_list[p]]+
                                                            QEC_dict[str(1)]['positive'][p_list[p]]+
                                                            QEC_dict[str(-1)]['negative'][p_list[p]]+
                                                            QEC_dict[str(1)]['negative'][p_list[p]])/4

    QEC_data_dict['x'] = QEC_dict[str(1)]['positive']['x']
    QEC_data_dict['folder'] = folder

    return QEC_data_dict, folder


def QEC_data_single_state_RO_single_error_sign(older_than = None,state = 'Z',RO = 0, sym = '00',e_sign = 1):

    QEC_data_dict = {}
    u_list = ['c0_u', 'c0_00_u','c0_01_u','c0_10_u','c0_11_u']
    c_list = ['c0', 'c0_00','c0_01','c0_10','c0_11']
    p_list = ['p00','p01','p10','p11']
    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']


    QEC_dict, folder = QEC_create_data_dict(older_than = older_than, RO = RO, state = state, sym = sym)

    for v in range(5):
        QEC_data_dict[y_list[v]] = {}
        QEC_data_dict[y_err_list[v]] = {}


        QEC_data_dict[y_list[v]] = (QEC_dict[str(e_sign)]['positive'][c_list[v]]-
                                                        QEC_dict[str(e_sign)]['negative'][c_list[v]])/2
        
        
        QEC_data_dict[y_err_list[v]] = ((QEC_dict[str(e_sign)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(e_sign)]['negative'][u_list[v]]**2)**0.5)/2
    for p in range(4):
            QEC_data_dict[p_list[p]] = {}
            
            QEC_data_dict[p_list[p]] = (    QEC_dict[str(e_sign)]['positive'][p_list[p]]+
                                                                QEC_dict[str(e_sign)]['negative'][p_list[p]])/2


    QEC_data_dict['x'] = QEC_dict[str(e_sign)]['positive']['x']
    QEC_data_dict['folder'] = folder

    return QEC_data_dict, folder

def QEC_plot_single_state_RO(older_than = None, no_error = '00',state = 'Z',RO = 0, e_sign = None,plot_guide = True):        
    

    if e_sign == None:
        QEC_data_dict, folder =  QEC_data_single_state_RO(older_than = older_than,state = state,RO = RO, sym = no_error)
    else:
        QEC_data_dict, folder =  QEC_data_single_state_RO_single_error_sign(older_than = older_than,state = state,RO = RO, sym = no_error,e_sign = e_sign)
    folder  = r'D:\measuring\data\QEC_data\figs'

    x = QEC_data_dict['x']
    y = QEC_data_dict['y']
    y_00 = QEC_data_dict['y_00']
    y_01 = QEC_data_dict['y_01']
    y_10 = QEC_data_dict['y_10']
    y_11 = QEC_data_dict['y_11']

    y_err = QEC_data_dict['y_err']
    y_err_00 = QEC_data_dict['y_err_00']
    y_err_01 = QEC_data_dict['y_err_01']
    y_err_10 = QEC_data_dict['y_err_10']
    y_err_11 = QEC_data_dict['y_err_11']

    p_00 = QEC_data_dict['p00']
    p_01 = QEC_data_dict['p01']
    p_10 = QEC_data_dict['p10']
    p_11 = QEC_data_dict['p11']

    x_g = [x[0],x[-1]]
    y_g = [y[0],y[-1]]

    fig,ax = plt.subplots() 
    ax.errorbar(x,y,yerr=y_err,color = 'k' )
    ax.plot(x_g,y_g,color = 'g' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_QEC')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')
    try:
        fig.savefig(
            os.path.join(folder,'errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_all'+'.png'))
    except:
        print 'Figure has not been saved.'

    fig,ax = plt.subplots() 
    ax.errorbar(x,y_00,yerr=y_err_00,color = 'c', label = 'y_00' )
    ax.errorbar(x,y_01,yerr=y_err_01,color = 'k', label = 'y_01' )
    ax.errorbar(x,y_10,yerr=y_err_10,color = 'm', label = 'y_10' )
    ax.errorbar(x,y_11,yerr=y_err_11,color = 'b', label = 'y_11' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    plt.legend()
    ax.set_title('errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_QEC_PS')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')

    try:
        fig.savefig(
            os.path.join(folder,'errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_ps'+'.png'))
    except:
        print 'Figure has not been saved.'




    fig,ax = plt.subplots()
    ax.set_title(str(folder)+'/'+ '\n probabilities')
    ax.plot(x,p_00, 'c', label = 'p00')
    ax.plot(x,p_01, 'k', label = 'p01')
    ax.plot(x,p_10, 'm', label = 'p10')
    ax.plot(x,p_11, 'b', label = 'p11')
    # ax.plot(x,p_00+p_01+p_10+p_11, 'g', label = 'sum')
    plt.legend()
    ax.set_xlabel('error probability')
    ax.set_ylabel('outcome probability')  
    ax.set_title('errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_QEC_probs')                
    
    try:
        fig.savefig(
            os.path.join(folder,'errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_probs'+'.png'))
    except:
        print 'Figure has not been saved.'

    return QEC_data_dict, folder

''' fitting the data '''

### Fitfunction

def fit_QEC(g_O, g_A, g_p):
    '''Fit function for QEC process fidelity data
    g_O -  Offset, given by the fidelity of the state that is insensitive to errors
    g_A -  Amplitude, g_iven by the fidelity of the states that are sensitive
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''test'''

    O   = fit.Parameter(g_O , 'O')
    A   = fit.Parameter(g_A, 'A')
    p   = fit.Parameter(g_p, 'p')

    p0 = [O, A, p]

    def fitfunc(x):
        '''test'''
        return (O() + A()*(  1-3*x+3*x**2-2*x**3 + 3*(2*p()-1)*(x-3*x**2+2*x**3)))

    return p0, fitfunc, fitfunc_str

def fit_QEC_curve(x,y):
        
        guess_O = 0.5
        guess_A = 0.5
        guess_p = 1
        p0, fitfunc, fitfunc_str = fit_QEC(guess_O, guess_A, guess_p)

        fit_result = fit.fit1d(x, y, fit_QEC,
                guess_O, guess_A, guess_p,
                fixed=[],
                do_print=True, ret=True)
        
        p02, fitfunc2, fitfunc_str2 = fit_QEC(fit_result['params_dict']['O'], fit_result['params_dict']['A'], fit_result['params_dict']['p'])

        x_temp      = np.linspace(x[0],x[-1],200)
        y_temp      =  fitfunc2(x_temp)

        return x_temp, y_temp,fit_result['params_dict']['p']

def fit_general_exponential(g_a, g_A, g_x0, g_T, g_n):
    fitfunc_str = 'a + A * exp(-((x-x0)/T )**n)'

    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    x0 = fit.Parameter(g_x0, 'x0')
    T = fit.Parameter(g_T, 'T')
    n = fit.Parameter(g_n, 'n')

    p0 = [a, A, x0, T, n]

    def fitfunc(x):
        return a() + A() * np.exp(-(x-x0())**n()/(T()**n()))
    return p0, fitfunc, fitfunc_str

''' these functions are used to open/close save/load from and to HDF5 files '''

def openfile(name = ''):
   datafile = h5py.File(os.path.join(r'K:\ns\qt\Diamond\Projects\QEC LT\QEC data10', name)) 
   return datafile

def closefile(datafile):
    datafile.close()

def openfile_single_state_RO_run(sym = '00', RO = 1, state = 'Z', error_sign = 1,el_RO = 'positive', run = 1):
    name = 'run_'+str(run)+ '_error_sym_' + sym + '_state_'+state + '_RO_'+str(RO)+ '_sign'+str(error_sign)+'_el_RO_'+el_RO+'.hdf5'
    print name
    datafile = h5py.File(os.path.join(r'D:\measuring\data\QEC_data\data', name)) 
    return datafile    

def save_data_single_state_RO_hdf5file(datafile, data_dict):
    f = datafile
    f_grp = f.create_group('data')

    for item in data_dict:
        f.attrs [item] = data_dict[item]
        f_grp.create_dataset (item, data = data_dict[item])

def load_single_data_hdf5file(datafile):
    f = datafile
    f_grp = f['/'+'data']

    data_dict = {}
    for item in f_grp.keys():
        data_dict[item] = f_grp[item].value
    return data_dict

''' here you save new data '''

def save_QEC_dataset_single_sign_single_elRO(older_than = None,sym = '11', RO = 1, state = 'Z', error_sign = 1, el_RO = 'positive', run =1, sweep_time =False ):
    

    QEC_temp_dict, folder = QEC_create_data_dict_single_error_single_elRO(older_than = older_than, RO = RO, state = state, 
                                                                        len_k = 6, sym = sym,error_sign = error_sign, el_RO = el_RO,sweep_time = sweep_time)

    if sweep_time == False:
        print 'yup'
        datafile = openfile_single_state_RO_run(sym = sym, RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)
    elif sweep_time ==True:
        datafile = openfile_single_state_RO_run(sym = sym  +'sweep_time_', RO = RO, state = state, error_sign = 0,el_RO = el_RO, run = run)
    
    save_data_single_state_RO_hdf5file(datafile, QEC_temp_dict)
       
    closefile(datafile)

def save_no_QEC_dataset_single_sign_single_elRO(idle = False, sweep_time = False, older_than = None, RO = 1, state = 'Z', error_sign = 1, el_RO = 'positive'):
    

    QEC_temp_dict, folder = no_QEC_create_data_dict_single_error_single_elRO(older_than = older_than, RO = RO, idle = idle, sweep_time = sweep_time,
                             state = state, error_sign = error_sign, el_RO = el_RO)
    if idle == True:
        datafile = openfile_single_state_RO_run(sym = 'no_correction_idle', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = 0)
    elif sweep_time == True:
        datafile = openfile_single_state_RO_run(sym = 'no_correction_sweep_time', RO = RO, state = state, error_sign = 0,el_RO = el_RO, run = 1)
    else:
        datafile = openfile_single_state_RO_run(sym = 'no_correction', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = 0)
    save_data_single_state_RO_hdf5file(datafile, QEC_temp_dict)
       
    closefile(datafile)

def save_single_Qubit_QEC_dataset_single_sign_single_elRO(older_than = None, Qubit =1, sweep_time = False,state = 'Z', error_sign = 1, el_RO = 'positive'):
    
    if state == 'X' or state == 'mX':
        RO = 2+(Qubit-1)*3
    elif state == 'Y' or state == 'mY':
        RO = 1+(Qubit-1)*3
    elif state == 'Z' or state == 'mZ': 
        RO = 0+(Qubit-1)*3
        if sweep_time == True:
            RO = Qubit-1

    QEC_temp_dict, folder = single_Qubit_QEC_create_data_dict_single_error_single_elRO(older_than = older_than, Qubit = Qubit, 
                                                            state = state, error_sign = error_sign, el_RO = el_RO, sweep_time =sweep_time)
    if sweep_time == False:
        datafile = openfile_single_state_RO_run(sym = 'single_qubit', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = 0)
    elif sweep_time == True:
        datafile = openfile_single_state_RO_run(sym = 'single_qubit_sweep_time', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = 0)
    save_data_single_state_RO_hdf5file(datafile, QEC_temp_dict)
       
    closefile(datafile)

def load_QEC_dataset_single_sign_single_elRO(sym = '11', RO = 1, state = 'Z', error_sign = 1, el_RO = 'positive', run =1,sweep_time = False):
    
    datafile = openfile_single_state_RO_run(sym = sym, RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)

    if sweep_time == False:
        datafile = openfile_single_state_RO_run(sym = sym, RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)
    elif sweep_time ==True:
        datafile = openfile_single_state_RO_run(sym = sym  +'sweep_time_' , RO = RO, state = state, error_sign = 0,el_RO = el_RO, run = run)

    QEC_temp_dict = load_single_data_hdf5file(datafile)
       
    closefile(datafile)

    return QEC_temp_dict

def load_no_QEC_dataset_single_sign_single_elRO(idle = False, sweep_time = False,RO = 1, state = 'Z', error_sign = 1, el_RO = 'positive'):
    
    if idle == False and sweep_time == False:
        datafile = openfile_single_state_RO_run(sym = 'no_correction', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = 0)
    elif sweep_time == True:

        datafile = openfile_single_state_RO_run(sym = 'no_correction_sweep_time', RO = RO, state = state, error_sign = 0,el_RO = el_RO, run = 1)
    else:
        datafile = openfile_single_state_RO_run(sym = 'no_correction_idle', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = 0)

    QEC_temp_dict = load_single_data_hdf5file(datafile)
       
    closefile(datafile)

    return QEC_temp_dict

def load_single_Qubit_QEC_dataset_single_sign_single_elRO(Qubit =1, state = 'Z', error_sign = 1,sweep_time = False, el_RO = 'positive'):
    
    if state == 'X' or state == 'mX':
        RO = 2+(Qubit-1)*3
    elif state == 'Y' or state == 'mY':
        RO = 1+(Qubit-1)*3
    elif state == 'Z' or state == 'mZ': 
        RO = 0+(Qubit-1)*3
        if sweep_time == True:
            RO = Qubit-1
    
    if sweep_time == False:
        datafile = openfile_single_state_RO_run(sym = 'single_qubit', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = 0)
    elif sweep_time == True:
        datafile = openfile_single_state_RO_run(sym = 'single_qubit_sweep_time', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = 0)

    QEC_temp_dict = load_single_data_hdf5file(datafile)
       
    closefile(datafile)

    return QEC_temp_dict

''' from here you can plot data taken from an existing HDF5 file '''

def QEC_sum_data_single_state_RO(run = 1, no_error = '00',state = 'Z',RO = 0):

    QEC_data_dict = {}
    u_list = ['c0_u', 'c0_00_u','c0_01_u','c0_10_u','c0_11_u']
    c_list = ['c0', 'c0_00','c0_01','c0_10','c0_11']
    p_list = ['p00','p01','p10','p11']
    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']
    
    if RO < 3:
        RO_C = 1/RO_corr_1qb
    else:
        RO_C = 1/RO_corr_3qb
    
    QEC_dict = {}
    
    for error_sign in [-1,1]:
        QEC_dict[str(error_sign)] = {}
        for el_RO in ['positive','negative']:
            QEC_dict[str(error_sign)][el_RO] = {}
            QEC_dict[str(error_sign)][el_RO] = load_QEC_dataset_single_sign_single_elRO(sym = no_error, RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)
    
    for v in range(5):
        QEC_data_dict[y_list[v]] = {}
        QEC_data_dict[y_err_list[v]] = {}


        QEC_data_dict[y_list[v]] = RO_C*(QEC_dict[str(-1)]['positive'][c_list[v]]+
                                                        QEC_dict[str(1)]['positive'][c_list[v]]-
                                                        QEC_dict[str(-1)]['negative'][c_list[v]]-
                                                        QEC_dict[str(1)]['negative'][c_list[v]])/4
        
        
        QEC_data_dict[y_err_list[v]] = RO_C*((QEC_dict[str(-1)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(1)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(-1)]['negative'][u_list[v]]**2+
                                                        QEC_dict[str(1)]['negative'][u_list[v]]**2)**0.5)/4
    
    QEC_data_dict['RO_contrast_pos_error'] =        (QEC_dict[str(1)]['positive']['c0']+
                                                    QEC_dict[str(1)]['negative']['c0'])/2     
    
    QEC_data_dict['RO_contrast_neg_error'] =        (QEC_dict[str(-1)]['positive']['c0']+
                                                    QEC_dict[str(-1)]['negative']['c0'])/2     

    QEC_data_dict['u_RO_contrast_pos_error'] =        (QEC_dict[str(1)]['positive']['c0_u']**2+
                                                    QEC_dict[str(1)]['negative']['c0_u']**2)**0.5/2     
    
    QEC_data_dict['u_RO_contrast_neg_error'] =        (QEC_dict[str(-1)]['positive']['c0_u']**2+
                                                    QEC_dict[str(-1)]['negative']['c0_u']**2)**0.5/2     

    for p in range(4):
            QEC_data_dict[p_list[p]] = {}
            
            QEC_data_dict[p_list[p]] = (QEC_dict[str(-1)]['positive'][p_list[p]]+
                                                            QEC_dict[str(1)]['positive'][p_list[p]]+
                                                            QEC_dict[str(-1)]['negative'][p_list[p]]+
                                                            QEC_dict[str(1)]['negative'][p_list[p]])/4

    QEC_data_dict['x'] = QEC_dict[str(1)]['positive']['x']

    return QEC_data_dict

def QEC_sum_data_single_state_RO_single_error_sign(run = 1, no_error = '00',state = 'Z',RO = 0,error_sign = 1, load_set = True, older_than = None,sweep_time = False):

    QEC_data_dict = {}
    u_list = ['c0_u', 'c0_00_u','c0_01_u','c0_10_u','c0_11_u']
    c_list = ['c0', 'c0_00','c0_01','c0_10','c0_11']
    p_list = ['p00','p01','p10','p11']
    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']
    QEC_dict = {}

    if RO < 3:
        RO_C = 1/RO_corr_1qb
    else:
        RO_C = 1/RO_corr_3qb

    QEC_dict[str(error_sign)] = {}
    for el_RO in ['positive','negative']:
        QEC_dict[str(error_sign)][el_RO] = {}
        if load_set == True:
            QEC_dict[str(error_sign)][el_RO] = load_QEC_dataset_single_sign_single_elRO(sym = no_error, RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run,sweep_time = sweep_time)
        elif load_set == False:
            QEC_dict[str(error_sign)][el_RO] , folder = QEC_create_data_dict_single_error_single_elRO(older_than = older_than, RO = RO, state = state, 
                                                                                                len_k = 6, sym = no_error,error_sign = error_sign, el_RO = el_RO, sweep_time = sweep_time)

    for v in range(5):
        QEC_data_dict[y_list[v]] = {}
        QEC_data_dict[y_err_list[v]] = {}


        QEC_data_dict[y_list[v]] = RO_C*(QEC_dict[str(error_sign)]['positive'][c_list[v]]-
                                                        QEC_dict[str(error_sign)]['negative'][c_list[v]])/2
        
        
        QEC_data_dict[y_err_list[v]] = RO_C*((QEC_dict[str(error_sign)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(error_sign)]['negative'][u_list[v]]**2)**0.5)/2
    for p in range(4):
            QEC_data_dict[p_list[p]] = {}
            
            QEC_data_dict[p_list[p]] = (    QEC_dict[str(error_sign)]['positive'][p_list[p]]+
                                                                QEC_dict[str(error_sign)]['negative'][p_list[p]])/2


    QEC_data_dict['x'] = QEC_dict[str(error_sign)]['positive']['x']

    return QEC_data_dict

def no_QEC_data_single_state_RO(idle = False, sweep_time = False, older_than = None,state = 'Z',RO = 0, load_set = False):

    QEC_data_dict = {}
    u_list = ['c0_u']
    c_list = ['c0']
    y_list = ['y']
    y_err_list = ['y_err']

    if RO < 3:
        RO_C = 1/RO_corr_1qb
    else:
        RO_C = 1/RO_corr_3qb


    QEC_dict = {}
    for error_sign in [-1,1]:
        QEC_dict[str(error_sign)] = {}
        for el_RO in ['positive','negative']:
            QEC_dict[str(error_sign)][el_RO] = {}
            if load_set == False:
                print 'ok'
                QEC_dict[str(error_sign)][el_RO], folder = no_QEC_create_data_dict_single_error_single_elRO(idle = idle, sweep_time = sweep_time, older_than = older_than, RO = RO, state = state, error_sign = error_sign, el_RO = el_RO)
            else:
                QEC_dict[str(error_sign)][el_RO]= load_no_QEC_dataset_single_sign_single_elRO(idle = idle, sweep_time = sweep_time, RO = RO, state = state, error_sign = error_sign,el_RO = el_RO)
    for v in range(1):
        QEC_data_dict[y_list[v]] = {}
        QEC_data_dict[y_err_list[v]] = {}
        QEC_data_dict[y_list[v]] = RO_C*(QEC_dict[str(-1)]['positive'][c_list[v]]+
                                                        QEC_dict[str(1)]['positive'][c_list[v]]-
                                                        QEC_dict[str(-1)]['negative'][c_list[v]]-
                                                        QEC_dict[str(1)]['negative'][c_list[v]])/4

        
        
        QEC_data_dict[y_err_list[v]] = RO_C* (QEC_dict[str(-1)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(1)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(-1)]['negative'][u_list[v]]**2+
                                                        QEC_dict[str(1)]['negative'][u_list[v]]**2)**0.5/4


    QEC_data_dict['RO_contrast_pos_error'] =        (QEC_dict[str(1)]['positive']['c0']+
                                                    QEC_dict[str(1)]['negative']['c0'])/2     
    
    QEC_data_dict['RO_contrast_neg_error'] =        (QEC_dict[str(-1)]['positive']['c0']+
                                                    QEC_dict[str(-1)]['negative']['c0'])/2     

    QEC_data_dict['u_RO_contrast_pos_error'] =        (QEC_dict[str(1)]['positive']['c0_u']**2+
                                                    QEC_dict[str(1)]['negative']['c0_u']**2)**0.5/2     
    
    QEC_data_dict['u_RO_contrast_neg_error'] =        (QEC_dict[str(-1)]['positive']['c0_u']**2+
                                                    QEC_dict[str(-1)]['negative']['c0_u']**2)**0.5/2    

    QEC_data_dict['x'] = QEC_dict[str(1)]['positive']['x']
    # QEC_data_dict['folder'] = folder

    return QEC_data_dict


def no_QEC_data_single_state_RO_single_error_sign(idle = False, sweep_time = False, older_than = None,state = 'Z',RO = 0,error_sign = 0, load_set = False):

    QEC_data_dict = {}
    u_list = ['c0_u']
    c_list = ['c0']
    y_list = ['y']
    y_err_list = ['y_err']

    QEC_dict = {}

    if RO < 3:
        RO_C = 1/RO_corr_1qb
    else:
        RO_C = 1/RO_corr_3qb

    QEC_dict[str(error_sign)] = {}
    for el_RO in ['positive','negative']:
        QEC_dict[str(error_sign)][el_RO] = {}
        if load_set == False:

            QEC_dict[str(error_sign)][el_RO], folder = no_QEC_create_data_dict_single_error_single_elRO(idle = idle, sweep_time = sweep_time, older_than = older_than, RO = RO, state = state, error_sign = error_sign, el_RO = el_RO)
        else:

            QEC_dict[str(error_sign)][el_RO]= load_no_QEC_dataset_single_sign_single_elRO(idle = idle, sweep_time = sweep_time, RO = RO, state = state, error_sign = error_sign,el_RO = el_RO)
    for v in range(1):
        QEC_data_dict[y_list[v]] = {}
        QEC_data_dict[y_err_list[v]] = {}
        QEC_data_dict[y_list[v]] = RO_C*(QEC_dict[str(error_sign)]['positive'][c_list[v]]-
                                                        QEC_dict[str(error_sign)]['negative'][c_list[v]])/2

        
        
        QEC_data_dict[y_err_list[v]] = RO_C*(QEC_dict[str(error_sign)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(error_sign)]['negative'][u_list[v]]**2)**0.5/2


    QEC_data_dict['x'] = QEC_dict[str(error_sign)]['positive']['x']
    return QEC_data_dict

def single_qubit_no_QEC_data_single_state_RO(older_than = None,state = 'Z',Qubit = 1,sweep_time =False, load_set = False):

    if state == 'X' or state == 'mX':
        RO = 2+(Qubit-1)*3
    elif state == 'Y' or state == 'mY':
        RO = 1+(Qubit-1)*3
    elif state == 'Z' or state == 'mZ': 
        RO = 0+(Qubit-1)*3


    RO_C = 1/RO_corr_1qb


    QEC_data_dict = {}
    u_list = ['c0_u']
    c_list = ['c0']
    y_list = ['y']
    y_err_list = ['y_err']

    QEC_dict = {}
    for error_sign in [-1,1]:
        QEC_dict[str(error_sign)] = {}
        for el_RO in ['positive','negative']:
            QEC_dict[str(error_sign)][el_RO] = {}
            if load_set == False:
                QEC_dict[str(error_sign)][el_RO], folder = single_Qubit_QEC_create_data_dict_single_error_single_elRO(older_than = older_than, Qubit = Qubit,sweep_time=sweep_time, state = state, error_sign = error_sign, el_RO = el_RO)
            else:
                QEC_dict[str(error_sign)][el_RO]= load_single_Qubit_QEC_dataset_single_sign_single_elRO(Qubit = Qubit, state = state, error_sign = error_sign,el_RO = el_RO)
    for v in range(1):
        QEC_data_dict[y_list[v]] = {}
        QEC_data_dict[y_err_list[v]] = {}
        QEC_data_dict[y_list[v]] = RO_C*(QEC_dict[str(-1)]['positive'][c_list[v]]+
                                                        QEC_dict[str(1)]['positive'][c_list[v]]-
                                                        QEC_dict[str(-1)]['negative'][c_list[v]]-
                                                        QEC_dict[str(1)]['negative'][c_list[v]])/4

        
        
        QEC_data_dict[y_err_list[v]] = RO_C* (QEC_dict[str(-1)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(1)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(-1)]['negative'][u_list[v]]**2+
                                                        QEC_dict[str(1)]['negative'][u_list[v]]**2)**0.5/4


    QEC_data_dict['RO_contrast_pos_error'] =        (QEC_dict[str(1)]['positive']['c0']+
                                                    QEC_dict[str(1)]['negative']['c0'])/2     
    
    QEC_data_dict['RO_contrast_neg_error'] =        (QEC_dict[str(-1)]['positive']['c0']+
                                                    QEC_dict[str(-1)]['negative']['c0'])/2     

    QEC_data_dict['u_RO_contrast_pos_error'] =        (QEC_dict[str(1)]['positive']['c0_u']**2+
                                                    QEC_dict[str(1)]['negative']['c0_u']**2)**0.5/2     
    
    QEC_data_dict['u_RO_contrast_neg_error'] =        (QEC_dict[str(-1)]['positive']['c0_u']**2+
                                                    QEC_dict[str(-1)]['negative']['c0_u']**2)**0.5/2    

    QEC_data_dict['x'] = QEC_dict[str(1)]['positive']['x']
    # QEC_data_dict['folder'] = folder

    return QEC_data_dict

def single_qubit_no_QEC_data_single_state_RO_single_error_sign(older_than = None,state = 'Z',Qubit = 1,error_sign = -1, sweep_time =False, load_set = False):

    if state == 'X' or state == 'mX':
        RO = 2+(Qubit-1)*3
    elif state == 'Y' or state == 'mY':
        RO = 1+(Qubit-1)*3
    elif state == 'Z' or state == 'mZ': 
        RO = [Qubit-1]

    QEC_data_dict = {}
    u_list = ['c0_u']
    c_list = ['c0']
    y_list = ['y']
    y_err_list = ['y_err']


    RO_C = 1/RO_corr_1qb

    QEC_dict = {}

    QEC_dict[str(error_sign)] = {}
    for el_RO in ['positive','negative']:
        QEC_dict[str(error_sign)][el_RO] = {}
        if load_set == False:
            QEC_dict[str(error_sign)][el_RO], folder = single_Qubit_QEC_create_data_dict_single_error_single_elRO(older_than = older_than, Qubit = Qubit,sweep_time=sweep_time, state = state, error_sign = error_sign, el_RO = el_RO)
        else:
            QEC_dict[str(error_sign)][el_RO]= load_single_Qubit_QEC_dataset_single_sign_single_elRO(Qubit = Qubit, state = state, sweep_time = sweep_time, error_sign = error_sign,el_RO = el_RO)
    for v in range(1):
        QEC_data_dict[y_list[v]] = {}
        QEC_data_dict[y_err_list[v]] = {}
        QEC_data_dict[y_list[v]] = RO_C*(QEC_dict[str(error_sign)]['positive'][c_list[v]]-
                                                       QEC_dict[str(error_sign)]['negative'][c_list[v]])/2

        
        
        QEC_data_dict[y_err_list[v]] = RO_C*(QEC_dict[str(error_sign)]['positive'][u_list[v]]**2+
                                                       QEC_dict[str(error_sign)]['negative'][u_list[v]]**2)**0.5/2


    QEC_data_dict['x'] = QEC_dict[str(error_sign)]['positive']['x']
    # QEC_data_dict['folder'] = folder

    return QEC_data_dict

def undo_correction_single_state_RO(run = 1, no_error = '00',state = 'Z',RO = 0):

    dataset_dict_full = QEC_sum_data_single_state_RO(run = run, no_error = no_error,state = state,RO = RO)

    p_list = ['p00','p01','p10','p11']
    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']

    no_error_list = [int(no_error[0]),int(no_error[1])]
    # determine error on QB, give example for 00
    no = no_error #no error, detected by 00
    Q1 = str((no_error_list[0]+1)%2)+str((no_error_list[1]+1)%2) # error on Q1, detected by 11, Carbon 1
    Q2 = str(no_error_list[0])+str((no_error_list[1]+1)%2) # error on Q2, detected by 01
    Q3 = str((no_error_list[0]+1)%2)+str(no_error_list[1])# error on Q3, detected by 10

    RO_state = state+str(RO)
    if RO_state in ['Y6','Z1','mY6','mZ1']:
        y_new = np.zeros(len(dataset_dict_full['p'+no]))
        for i in range(len(dataset_dict_full['p'+no])):
            y_new[i] = (dataset_dict_full['p'+no][i]*dataset_dict_full['y_'+no][i]
                    + dataset_dict_full['p'+Q1][i]*dataset_dict_full['y_'+Q1][i]
                    -dataset_dict_full['p'+Q2][i]*dataset_dict_full['y_'+Q2][i]
                    +dataset_dict_full['p'+Q3][i]*dataset_dict_full['y_'+Q3][i])
    elif RO_state in ['Y4','Z2','mY4','mZ2']:
        y_new = np.zeros(len(dataset_dict_full['p'+no]))
        for i in range(len(dataset_dict_full['p'+no])):
            y_new[i] = (dataset_dict_full['p'+no][i]*dataset_dict_full['y_'+no][i]
                    + dataset_dict_full['p'+Q1][i]*dataset_dict_full['y_'+Q1][i]
                    +dataset_dict_full['p'+Q2][i]*dataset_dict_full['y_'+Q2][i]
                    -dataset_dict_full['p'+Q3][i]*dataset_dict_full['y_'+Q3][i])
    elif RO_state in ['Y5','Z0','mY5','mZ0']:
        y_new = np.zeros(len(dataset_dict_full['p'+no]))
        for i in range(len(dataset_dict_full['p'+no])):
            y_new[i] = (dataset_dict_full['p'+no][i]*dataset_dict_full['y_'+no][i]
                    - dataset_dict_full['p'+Q1][i]*dataset_dict_full['y_'+Q1][i]
                    +dataset_dict_full['p'+Q2][i]*dataset_dict_full['y_'+Q2][i]
                    +dataset_dict_full['p'+Q3][i]*dataset_dict_full['y_'+Q3][i])
    else:
        y_new = dataset_dict_full['y']

    

    return y_new

def undo_correction_single_state_RO_error_sign(run = 1, no_error = '00',state = 'Z',RO = 0,error_sign = 1,sweep_time=False):
    if sweep_time == False:
        dataset_dict_full = QEC_sum_data_single_state_RO_single_error_sign(run = run, no_error = no_error,state = state,RO = RO,error_sign = error_sign, sweep_time=False)
    elif sweep_time == True:
        dataset_dict_full = QEC_sum_data_single_state_RO_single_error_sign(no_error = no_error,state = state,RO = RO,load_set = True,sweep_time = True)

    p_list = ['p00','p01','p10','p11']
    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']

    no_error_list = [int(no_error[0]),int(no_error[1])]
    # determine error on QB, give example for 00
    no = no_error #no error, detected by 00
    Q1 = str((no_error_list[0]+1)%2)+str((no_error_list[1]+1)%2) # error on Q1, detected by 11, Carbon 1
    Q2 = str(no_error_list[0])+str((no_error_list[1]+1)%2) # error on Q2, detected by 01
    Q3 = str((no_error_list[0]+1)%2)+str(no_error_list[1])# error on Q3, detected by 10

    RO_state = state+str(RO)
    if RO_state in ['Y6','Z1','mY6','mZ1']:
        y_new = np.zeros(len(dataset_dict_full['p'+no]))
        for i in range(len(dataset_dict_full['p'+no])):
            y_new[i] = (dataset_dict_full['p'+no][i]*dataset_dict_full['y_'+no][i]
                    + dataset_dict_full['p'+Q1][i]*dataset_dict_full['y_'+Q1][i]
                    -dataset_dict_full['p'+Q2][i]*dataset_dict_full['y_'+Q2][i]
                    +dataset_dict_full['p'+Q3][i]*dataset_dict_full['y_'+Q3][i])
    elif RO_state in ['Y4','Z2','mY4','mZ2']:
        y_new = np.zeros(len(dataset_dict_full['p'+no]))
        for i in range(len(dataset_dict_full['p'+no])):
            y_new[i] = (dataset_dict_full['p'+no][i]*dataset_dict_full['y_'+no][i]
                    + dataset_dict_full['p'+Q1][i]*dataset_dict_full['y_'+Q1][i]
                    +dataset_dict_full['p'+Q2][i]*dataset_dict_full['y_'+Q2][i]
                    -dataset_dict_full['p'+Q3][i]*dataset_dict_full['y_'+Q3][i])
    elif RO_state in ['Y5','Z0','mY5','mZ0']:
        y_new = np.zeros(len(dataset_dict_full['p'+no]))
        for i in range(len(dataset_dict_full['p'+no])):
            y_new[i] = (dataset_dict_full['p'+no][i]*dataset_dict_full['y_'+no][i]
                    - dataset_dict_full['p'+Q1][i]*dataset_dict_full['y_'+Q1][i]
                    +dataset_dict_full['p'+Q2][i]*dataset_dict_full['y_'+Q2][i]
                    +dataset_dict_full['p'+Q3][i]*dataset_dict_full['y_'+Q3][i])
    else:
        y_new = dataset_dict_full['y']

    

    return y_new

def QEC_state_sum_all(state = 'Z', RO = 0,run_list_00 = [1,2,3],run_list_01 = [1,2,3],run_list_10 = [1,2],run_list_11 = [1,2,3]):

    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']
    state_dict_single = {}
    state_dict = {}

    for j, run in enumerate(['00','01','10','11']):
        if run == '00':
            run_list = run_list_00
        elif run == '01':
            run_list = run_list_01
        elif run == '10':
            run_list = run_list_10
        elif run == '11':
            run_list = run_list_11
        
        #sum over runs
        for i in range(len(run_list)):
            state_dict_single[run] = {}
            state_dict_single[run] = QEC_sum_data_single_state_RO(run = run_list[i], no_error = run,state = state,RO = RO)
            # print state_dict_single
            if i ==0:
                
                state_dict[run] = {}
                for k, yy in enumerate(y_list):
                    state_dict[run][yy] = 1/float(len(run_list))*state_dict_single[run][ yy]
                    state_dict[run]['temp'+ y_err_list[k]] = state_dict_single[run][ y_err_list[k]]**2
            else:
                for k, yy in enumerate(y_list):
                    state_dict[run][ yy] += 1/float(len(run_list))*state_dict_single[run][ yy]
                    state_dict[run]['temp'+  y_err_list[k]] += state_dict_single[run][ y_err_list[k]]**2
                
            for k, yy in enumerate(y_list):
                state_dict[run][ y_err_list[k]] = 1/float(len(run_list))*(state_dict[run]['temp'+ y_err_list[k]])**0.5

        

    state_dict['x'] = state_dict_single[run]['x']


        # sum over all
    for j, run in enumerate(['00','01','10','11']):
        state_dict_single = state_dict[run]

        if j ==0:
            summed_dict = {}
            
            for k, yy in enumerate(y_list):
                summed_dict[yy] = 1/4.*state_dict_single[yy]
                summed_dict['temp'+y_err_list[k]] = state_dict_single[y_err_list[k]]**2
        else:
        
            for k, yy in enumerate(y_list):
                summed_dict[yy] += 1/4.*state_dict_single[yy]
                summed_dict['temp'+ y_err_list[k]] += state_dict_single[y_err_list[k]]**2
                        
        for k, yy in enumerate(y_list):
            summed_dict[ y_err_list[k]] = 1/4.*(summed_dict['temp'+ y_err_list[k]])**0.5





    summed_dict['x'] = state_dict['x']
    
    return summed_dict

''' plot single QEC / no QEC lines '''

def QEC_plot_single_state_RO_saved_data(run = 1, no_error = '00',state = 'Z',RO = 0, plot_separate = False,plot_guide = True,plot_no_correct = False):        
    
    
    dataset_dict_full = QEC_sum_data_single_state_RO(run = run, no_error = no_error,state = state,RO = RO)
    QEC_data_dict  = dataset_dict_full
        
    folder  = r'D:\measuring\data\QEC_data\figs'

    x = QEC_data_dict['x']
    y = QEC_data_dict['y']
    y_00 = QEC_data_dict['y_00']
    y_01 = QEC_data_dict['y_01']
    y_10 = QEC_data_dict['y_10']
    y_11 = QEC_data_dict['y_11']

    y_err = QEC_data_dict['y_err']
    y_err_00 = QEC_data_dict['y_err_00']
    y_err_01 = QEC_data_dict['y_err_01']
    y_err_10 = QEC_data_dict['y_err_10']
    y_err_11 = QEC_data_dict['y_err_11']

    p_00 = QEC_data_dict['p00']
    p_01 = QEC_data_dict['p01']
    p_10 = QEC_data_dict['p10']
    p_11 = QEC_data_dict['p11']

    x_g = [x[0],x[-1]]
    y_g = [y[0],y[-1]]

    fig,ax = plt.subplots() 
    ax.errorbar(x,y,yerr=y_err,color = 'k' )
    if plot_guide == True:
        ax.plot(x_g,y_g,color = 'g' )
    if plot_no_correct == True:
        y_no_corr = undo_correction_single_state_RO(run = run, no_error = no_error,state = state,RO = RO)
        print y_no_corr
        ax.errorbar(x,y_no_corr,yerr = y_err, color = 'm')
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')

    if plot_no_correct == True:
        try:
            fig.savefig(
                os.path.join(folder,'undo_correct_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_all'+'.png'))
        except:
            print 'Figure has not been saved.'
    else:
        try:
            fig.savefig(
                os.path.join(folder,'error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_all'+'.png'))
        except:
            print 'Figure has not been saved.'

    fig,ax = plt.subplots() 
    ax.errorbar(x,y_00,yerr=y_err_00,color = 'c', label = 'y_00' )
    ax.errorbar(x,y_01,yerr=y_err_01,color = 'k', label = 'y_01' )
    ax.errorbar(x,y_10,yerr=y_err_10,color = 'm', label = 'y_10' )
    ax.errorbar(x,y_11,yerr=y_err_11,color = 'b', label = 'y_11' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    plt.legend()
    ax.set_title('error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC_PS')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')

    try:
        fig.savefig(
            os.path.join(folder,'error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_ps'+'.png'))
    except:
        print 'Figure has not been saved.'




    fig,ax = plt.subplots()
    ax.set_title(str(folder)+'/'+ '\n probabilities')
    ax.plot(x,p_00, 'c', label = 'p00')
    ax.plot(x,p_01, 'k', label = 'p01')
    ax.plot(x,p_10, 'm', label = 'p10')
    ax.plot(x,p_11, 'b', label = 'p11')
    # ax.plot(x,p_00+p_01+p_10+p_11, 'g', label = 'sum')
    plt.legend()
    ax.set_xlabel('error probability')
    ax.set_ylabel('outcome probability')  
    ax.set_title('error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC_probs')                
    
    try:
        fig.savefig(
            os.path.join(folder,'error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_probs'+'.png'))
    except:
        print 'Figure has not been saved.'

    
    if plot_separate == True:
        QEC_data_dict_n = QEC_sum_data_single_state_RO_single_error_sign(run = run, no_error = no_error,state = state,RO = RO,error_sign = -1)
        QEC_data_dict_p = QEC_sum_data_single_state_RO_single_error_sign(run = run, no_error = no_error,state = state,RO = RO,error_sign = 1)

        QEC_data_dict_pp = load_QEC_dataset_single_sign_single_elRO(sym = no_error, RO = RO, state = state, error_sign = 1,el_RO = 'positive', run = run)
        QEC_data_dict_np = load_QEC_dataset_single_sign_single_elRO(sym = no_error, RO = RO, state = state, error_sign = -1,el_RO = 'positive', run = run)
        QEC_data_dict_pn = load_QEC_dataset_single_sign_single_elRO(sym = no_error, RO = RO, state = state, error_sign = 1,el_RO = 'negative', run = run)
        QEC_data_dict_nn = load_QEC_dataset_single_sign_single_elRO(sym = no_error, RO = RO, state = state, error_sign = -1,el_RO = 'negative', run = run)


        linestyle = [':', '--']
        label = ['negative error','positive error']
        # for QEC_data_dict, ii in enumerate[QEC_data_dict_pp,QEC_data_dict_pn, QEC_data_dict_np, QEC_data_dict_nn]:
        for ii in range(2):
            if ii ==0:
                QEC_data_dict = QEC_data_dict_n
            if ii ==1:
                QEC_data_dict = QEC_data_dict_p



            x = QEC_data_dict['x']
            y = QEC_data_dict['y']
            y_00 = QEC_data_dict['y_00']
            y_01 = QEC_data_dict['y_01']
            y_10 = QEC_data_dict['y_10']
            y_11 = QEC_data_dict['y_11']

            y_err = QEC_data_dict['y_err']
            y_err_00 = QEC_data_dict['y_err_00']
            y_err_01 = QEC_data_dict['y_err_01']
            y_err_10 = QEC_data_dict['y_err_10']
            y_err_11 = QEC_data_dict['y_err_11']

            p_00 = QEC_data_dict['p00']
            p_01 = QEC_data_dict['p01']
            p_10 = QEC_data_dict['p10']
            p_11 = QEC_data_dict['p11']



            # plt.figure(10)
            if ii == 0:
                fig1,ax1 = plt.subplots() 
            ax1.errorbar(x,y,yerr=y_err, ls = linestyle[ii], color = 'b', label = label[ii] )
            if ii ==1:  
                lgd = ax1.legend(loc = 2, bbox_to_anchor = (1,1))
                ax1.set_ylim(-1.1,1.1)
                ax1.set_xlim(-0.1,1.1)
                ax1.set_title('single_error_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC')
                ax1.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
                ax1.set_xlabel('error probability')
                ax1.set_ylabel('Contrast')

                try:
                    fig1.savefig(
                        os.path.join(folder,'single_error_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_all'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
                except:
                    print 'Figure has not been saved.'

            # plt.figure(11)
            if ii == 0:
                fig2,ax2 = plt.subplots() 
            ax2.errorbar(x,y_00,yerr=y_err_00,color = 'c', label = 'y_00_'+label[ii],ls = linestyle[ii] )
            ax2.errorbar(x,y_01,yerr=y_err_01,color = 'k', label = 'y_01_'+label[ii],ls = linestyle[ii] )
            ax2.errorbar(x,y_10,yerr=y_err_10,color = 'm', label = 'y_10_'+label[ii],ls = linestyle[ii] )
            ax2.errorbar(x,y_11,yerr=y_err_11,color = 'b', label = 'y_11_'+label[ii],ls = linestyle[ii] )
            ax2.set_ylim(-1.1,1.1)
            ax2.set_xlim(-0.1,1.1)
            if ii ==1:  
                lgd = ax2.legend(loc = 2, bbox_to_anchor = (1,1))
                ax2.set_title('single_error_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC_PS')
                ax2.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
                ax2.set_xlabel('error probability')
                ax2.set_ylabel('Contrast')

                try:
                    fig2.savefig(
                        os.path.join(folder,'single_error_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_ps'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
                except:
                    print 'Figure has not been saved.'



            # plt.figure(12)
            if ii == 0:
                fig3,ax3 = plt.subplots()
            ax3.set_title(str(folder)+'/'+ '\n probabilities')
            ax3.plot(x,p_00, linestyle[ii], color = 'c', label = 'p00_'+label[ii])
            ax3.plot(x,p_01, linestyle[ii], color = 'k', label = 'p01_'+label[ii])
            ax3.plot(x,p_10, linestyle[ii], color = 'm', label = 'p10_'+label[ii])
            ax3.plot(x,p_11, linestyle[ii], color = 'b', label = 'p11_'+label[ii])
            # ax3.plot(x,p_00+p_01+p_10+p_11, 'g', label = 'sum')
            if ii ==1:  
                lgd = ax3.legend(loc = 2, bbox_to_anchor = (1,1))
                ax3.set_xlabel('error probability')
                ax3.set_ylabel('outcome probability')  
                ax3.set_title('single_error_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC_probs')                
                
                try:
                    fig3.savefig(
                        os.path.join(folder,'single_error_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_probs'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
                except:
                    print 'Figure has not been saved.'



        linestyle = ['-.', '-', '--', ':']
        label = ['pos_RO_pos_err','pos_RO_neg_err','neg_RO_pos_err','neg_RO_neg_err']
        # for QEC_data_dict, ii in enumerate[QEC_data_dict_pp,QEC_data_dict_pn, QEC_data_dict_np, QEC_data_dict_nn]:
        for ii in range(4):
            if ii ==0:
                QEC_data_dict = QEC_data_dict_pp
                sign = 1
            if ii ==1:
                QEC_data_dict = QEC_data_dict_np
                sign = 1
            if ii ==2:
                QEC_data_dict = QEC_data_dict_pn
                sign = -1
            if ii ==3:
                QEC_data_dict = QEC_data_dict_nn
                sign = -1
            # print QEC_data_dict
                # print item
            x = QEC_data_dict['x']
            y = sign*QEC_data_dict['c0']
            y_00 = sign*QEC_data_dict['c0_00']
            y_01 = sign*QEC_data_dict['c0_01']
            y_10 = sign*QEC_data_dict['c0_10']
            y_11 = sign*QEC_data_dict['c0_11']

            y_err = QEC_data_dict['c0_u']
            y_err_00 = QEC_data_dict['c0_00_u']
            y_err_01 = QEC_data_dict['c0_01_u']
            y_err_10 = QEC_data_dict['c0_10_u']
            y_err_11 = QEC_data_dict['c0_11_u']

            p_00 = QEC_data_dict['p00']
            p_01 = QEC_data_dict['p01']
            p_10 = QEC_data_dict['p10']
            p_11 = QEC_data_dict['p11']



            # plt.figure(10)
            if ii == 0:
                fig1,ax1 = plt.subplots() 
            ax1.errorbar(x,y,yerr=y_err, ls = linestyle[ii], color = 'b', label = label[ii] )
            if ii ==3:  
                lgd = ax1.legend(loc = 2, bbox_to_anchor = (1,1))
                ax1.set_ylim(-1.1,1.1)
                ax1.set_xlim(-0.1,1.1)
                ax1.set_title('single_lines_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC')
                ax1.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
                ax1.set_xlabel('error probability')
                ax1.set_ylabel('Contrast')

                try:
                    fig1.savefig(
                        os.path.join(folder,'single_lines_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_all'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
                except:
                    print 'Figure has not been saved.'

            # plt.figure(11)
            if ii == 0:
                fig2,ax2 = plt.subplots() 
            ax2.errorbar(x,y_00,yerr=y_err_00,color = 'c', label = 'y_00_'+label[ii],ls = linestyle[ii] )
            ax2.errorbar(x,y_01,yerr=y_err_01,color = 'k', label = 'y_01_'+label[ii],ls = linestyle[ii] )
            ax2.errorbar(x,y_10,yerr=y_err_10,color = 'm', label = 'y_10_'+label[ii],ls = linestyle[ii] )
            ax2.errorbar(x,y_11,yerr=y_err_11,color = 'b', label = 'y_11_'+label[ii],ls = linestyle[ii] )
            ax2.set_ylim(-1.1,1.1)
            ax2.set_xlim(-0.1,1.1)
            if ii ==3:  
                lgd = ax2.legend(loc = 2, bbox_to_anchor = (1,1))
                ax2.set_title('single_lines_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC_PS')
                ax2.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
                ax2.set_xlabel('error probability')
                ax2.set_ylabel('Contrast')

                try:
                    fig2.savefig(
                        os.path.join(folder,'single_lines_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_ps'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
                except:
                    print 'Figure has not been saved.'



            # plt.figure(12)
            if ii == 0:
                fig3,ax3 = plt.subplots()
            ax3.set_title(str(folder)+'/'+ '\n probabilities')
            ax3.plot(x,p_00, linestyle[ii], color = 'c', label = 'p00_'+label[ii])
            ax3.plot(x,p_01, linestyle[ii], color = 'k', label = 'p01_'+label[ii])
            ax3.plot(x,p_10, linestyle[ii], color = 'm', label = 'p10_'+label[ii])
            ax3.plot(x,p_11, linestyle[ii], color = 'b', label = 'p11_'+label[ii])
            # ax3.plot(x,p_00+p_01+p_10+p_11, 'g', label = 'sum')
            if ii ==3:  
                lgd = ax3.legend(loc = 2, bbox_to_anchor = (1,1))
                ax3.set_xlabel('error probability')
                ax3.set_ylabel('outcome probability')  
                ax3.set_title('single_lines_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC_probs')                
                
                try:
                    fig3.savefig(
                        os.path.join(folder,'single_lines_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_probs'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
                except:
                    print 'Figure has not been saved.'



    return QEC_data_dict, folder

def QEC_sweep_time_sum_error_syns(RO  = 0, state = 'Z',run_list = ['00','11','01']):
    
    p_list = ['p00','p01','p10','p11']
    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']
    dataset_dict_full= {}
    
    for i, run in enumerate(run_list):

            print 'test'
            print run
            dataset_dict_full[run] = QEC_sum_data_single_state_RO_single_error_sign(no_error = run,state = state,RO = RO,load_set = True,sweep_time = True)
            
            if i ==0:
                print float(len(run_list))
                summed_dict = {}
                for k, yy in enumerate(y_list):
                    summed_dict[yy] = 1/float(len(run_list))*dataset_dict_full[run][yy]
                    summed_dict['temp'+ y_err_list[k]] = dataset_dict_full[run][y_err_list[k]]**2
                for jj, p in enumerate(p_list):
                        summed_dict[p_list[jj]]=1/float(len(run_list))*dataset_dict_full[run][p]                    
            else:
                for k, yy in enumerate(y_list):
                    summed_dict[yy] += 1/float(len(run_list))*dataset_dict_full[run][yy]
                    summed_dict['temp'+ y_err_list[k]] += dataset_dict_full[run][y_err_list[k]]**2
                else:
                    summed_dict[p_list[jj]]+=1/float(len(run_list))*dataset_dict_full[run][p]                          
            
            for k, yy in enumerate(y_list):
                summed_dict[y_err_list[k]] = 1/float(len(run_list))*(summed_dict['temp'+ y_err_list[k]])**0.5




    summed_dict['x'] = dataset_dict_full[run]['x']
    
    return summed_dict

def QEC_sweep_time_sum_states(RO  = 0):
    
    p_list = ['p00','p01','p10','p11']
    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']
    dataset_dict_full= {}

    dataset_dict_full['Z'] = QEC_sweep_time_sum_error_syns(RO  = RO, state = 'Z',run_list = ['00','11','01'])
    dataset_dict_full['mZ'] = QEC_sweep_time_sum_error_syns(RO  = RO, state = 'Z',run_list = ['00','11','01'])


    summed_dict = {}
    for k, yy in enumerate(y_list):
        summed_dict[yy] = 1/2.*(dataset_dict_full['Z'][yy]-dataset_dict_full['mZ'][yy])
        summed_dict[y_err_list[k]] = 1/2.*(dataset_dict_full['Z'][y_err_list[k]]**2+dataset_dict_full['Z'][y_err_list[k]]**2)**0.5
    for jj, p in enumerate(p_list):
        summed_dict[p_list[jj]]=1/2.*(dataset_dict_full['Z'][p]+dataset_dict_full['mZ'][p])                    

    summed_dict['x'] = dataset_dict_full['Z']['x']
    
    return summed_dict

def QEC_plot_single_state_sweep_time(older_than = '20150107_090000',run = 1, no_error = '',no_error_list = [],state = 'Z',add_encode = False, add_single =False, plot_no_correct = False,load_set = False):        
    fig1,ax1 = plt.subplots() 
    fig2,ax2 = plt.subplots() 
    fig3,ax3 = plt.subplots() 
    color = ['r','g','b']
    dataset_dict_full = {}
    no_QEC_data_dict = {}
    QEC_single_data_dict = {}

    for RO in [0,1,2]:
        if no_error_list != []:
            dataset_dict_full[RO] = QEC_sweep_time_sum_error_syns(state = state,RO = RO,run_list = no_error_list)
        else:
            dataset_dict_full[RO] = QEC_sum_data_single_state_RO_single_error_sign(no_error = no_error,state = state,RO = RO,load_set = load_set, older_than = older_than,sweep_time = True)

        QEC_data_dict  = dataset_dict_full[RO]
            
        folder  = r'D:\measuring\data\QEC_data\figs\timesweep'

        parity_time = 2*(4.996e-6*34 +11.312e-6*48) +2*(13.616e-6*34+4.996e-6*34) + 2* 150e-6
        x = QEC_data_dict['x']+ np.ones(len(QEC_data_dict['x']))*parity_time
        
        y = QEC_data_dict['y']
        y_00 = QEC_data_dict['y_00']
        y_01 = QEC_data_dict['y_01']
        y_10 = QEC_data_dict['y_10']
        y_11 = QEC_data_dict['y_11']

        y_err = QEC_data_dict['y_err']
        y_err_00 = QEC_data_dict['y_err_00']
        y_err_01 = QEC_data_dict['y_err_01']
        y_err_10 = QEC_data_dict['y_err_10']
        y_err_11 = QEC_data_dict['y_err_11']

        p_00 = QEC_data_dict['p00']
        p_01 = QEC_data_dict['p01']
        p_10 = QEC_data_dict['p10']
        p_11 = QEC_data_dict['p11']

        
        ax1.errorbar(x,y,yerr=y_err,color = color[RO], label = 'QEC, decode to Qubit '+str(RO+1))

        if plot_no_correct == True:
            y_no_corr = undo_correction_single_state_RO_error_sign(run = 1, no_error = no_error,state = state,RO = RO,error_sign = 1,sweep_time=True)
            dataset_dict_full[RO]['y_no_corr'] = y_no_corr
            ax1.errorbar(x,y_no_corr,yerr = y_err, color = color[RO],ls = '-.', label = 'undo correction, Q'+str(RO+1))



        ax2.errorbar(x,y_00,yerr=y_err_00,color = 'c', label = 'y_00' )
        ax2.errorbar(x,y_01,yerr=y_err_01,color = 'k', label = 'y_01' )
        ax2.errorbar(x,y_10,yerr=y_err_10,color = 'm', label = 'y_10' )
        ax2.errorbar(x,y_11,yerr=y_err_11,color = 'b', label = 'y_11' )
        ax2.set_ylim(-1.1,1.1)
        ax2.set_xlim(-1e-3,35e-3)
        ax2.legend()
        ax2.set_title('sweep_time_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC_PS')
        ax2.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax2.set_xlabel('time (s)')
        ax2.set_ylabel('Contrast')

        ax3.set_title(str(folder)+'/'+ '\n probabilities')
        ax3.plot(x,p_00, 'c', label = 'p00')
        ax3.plot(x,p_01, 'k', label = 'p01')
        ax3.plot(x,p_10, 'm', label = 'p10')
        ax3.plot(x,p_11, 'b', label = 'p11')
        ax3.set_xlim(-1e-3,35e-3)
        # ax3.legend()
        ax3.set_xlabel('time (s)')
        ax3.set_ylabel('outcome probability')  
        ax3.set_title('sweep_time_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC_probs')                
        
        if add_encode == True:
            no_QEC_data_dict[RO] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = state,RO = RO, load_set = load_set,error_sign = 0)
            if RO == 1:
                x = no_QEC_data_dict[RO]['x']
                y = no_QEC_data_dict[RO]['y']
                y_err = no_QEC_data_dict[RO]['y_err']
                ax1.errorbar(x,y,yerr=y_err,color = color[RO],ls = '--', label = 'Encoding, sweep time, Q' + str(RO+1) )

        if add_single == True and RO == 1:
            QEC_single_data_dict[RO] =  single_qubit_no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,state = state,sweep_time = True, error_sign = -1, Qubit = RO+1, load_set = True)    
            x_single = QEC_single_data_dict[RO]['x']
            y_single = QEC_single_data_dict[RO]['y']
            y_single_err = QEC_single_data_dict[RO]['y_err']
          
            ax1.errorbar(x_single,y_single,yerr=y_single_err,color = color[RO],ls = ':', label = 'Single Qubit, sweep time, Q' + str(RO+1) )


    if no_error_list != []:
        dataset_dict_full[6] = QEC_sweep_time_sum_error_syns(state = state,RO = 6,run_list = no_error_list)
    else:
        dataset_dict_full[6] = QEC_sum_data_single_state_RO_single_error_sign(no_error = no_error,state = state,RO = 6,load_set = load_set, older_than = older_than,sweep_time = True)
    
    y_toff = 1/2.*(dataset_dict_full[0]['y']+dataset_dict_full[1]['y']+dataset_dict_full[2]['y']-dataset_dict_full[6]['y'])
    y_toff_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
    x = dataset_dict_full[6]['x']+ np.ones(len(dataset_dict_full[6]['x']))*parity_time
    ax1.errorbar(x,y_toff,yerr=y_toff_err,color = 'k', label = 'QEC+ majority vote' )
    if plot_no_correct == True:
        y_no_corr = undo_correction_single_state_RO_error_sign(run = 1, no_error = no_error,state = state,RO = RO,error_sign = 1,sweep_time=True)
        dataset_dict_full[6]['y_no_corr'] = y_no_corr
        y_toff = 1/2.*(dataset_dict_full[0]['y_no_corr']+dataset_dict_full[1]['y_no_corr']+dataset_dict_full[2]['y_no_corr']-dataset_dict_full[6]['y_no_corr'])
        y_toff_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
        x = dataset_dict_full[6]['x']+ np.ones(len(dataset_dict_full[6]['x']))*parity_time
        ax1.errorbar(x,y_toff,yerr=y_toff_err,color = 'k',ls = '-.', label = 'undo QEC+ toffoli' )    
    if add_encode == True:
        no_QEC_data_dict[6] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = state,RO = 6, load_set = load_set,error_sign = 0)
        
        y_toff = 1/2.*(no_QEC_data_dict[0]['y']+no_QEC_data_dict[1]['y']+no_QEC_data_dict[2]['y']-no_QEC_data_dict[6]['y'])
        y_toff_err = 1/2.*(no_QEC_data_dict[0]['y_err']**2+no_QEC_data_dict[1]['y_err']**2+no_QEC_data_dict[2]['y_err']**2+no_QEC_data_dict[6]['y_err']**2)**0.5
        x = no_QEC_data_dict[0]['x']
        ax1.errorbar(x,y_toff,yerr=y_toff_err,color = 'k',ls ='--', label = 'Encoding, sweep time, majority vote' )

    ax1.set_ylim(-1.1,1.1)
    ax1.set_xlim(-1e-3,35e-3)
    ax1.set_title('sweep_time_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC')
    ax1.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('Contrast')
    lgd = ax1.legend(loc = 2, bbox_to_anchor = (1,1))



    if plot_no_correct == True:
        try:
            fig1.savefig(
                os.path.join(folder,'sweep_time_undo_correct_error_syn_'+no_error+'_state_'+state+'_all'+'.png'))
        except:
            print 'Figure has not been saved.'
    else:
        try:
            fig1.savefig(
                os.path.join(folder,'sweep_time_error_syn_'+no_error+'_state_'+state+'_all'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
        except:
            print 'Figure has not been saved.'

    try:
        fig2.savefig(
            os.path.join(folder,'sweep_time_error_syn_'+no_error+'_state_'+state+'_ps'+'.png'))
    except:
        print 'Figure has not been saved.'

    try:
        fig3.savefig(
            os.path.join(folder,'sweep_time_error_syn_'+no_error+'_state_'+state+'_probs'+'.png'))
    except:
        print 'Figure has not been saved.'

# def QEC_plot_all_summed_sweep_time(add_encode = False, add_single =False, plot_no_correct = False,load_set = True):        
#     fig1,ax1 = plt.subplots() 
#     fig2,ax2 = plt.subplots() 
#     fig3,ax3 = plt.subplots() 
#     color = ['r','g','b']
#     dataset_dict_full = {}
#     no_QEC_data_dict = {}
#     QEC_single_data_dict = {}

#     for RO in [0,1,2]:

#         dataset_dict_full[RO] = QEC_sweep_time_sum_states(RO  = RO)
#         QEC_data_dict  = dataset_dict_full[RO]
            
#         folder  = r'D:\measuring\data\QEC_data\figs\timesweep'

#         parity_time = 2*(4.996e-6*34 +11.312e-6*48) +2*(13.616e-6*34+4.996e-6*34) + 2* 150e-6
#         x = QEC_data_dict['x']+ np.ones(len(QEC_data_dict['x']))*parity_time
        
#         y = QEC_data_dict['y']
#         y_00 = QEC_data_dict['y_00']
#         y_01 = QEC_data_dict['y_01']
#         y_10 = QEC_data_dict['y_10']
#         y_11 = QEC_data_dict['y_11']

#         y_err = QEC_data_dict['y_err']
#         y_err_00 = QEC_data_dict['y_err_00']
#         y_err_01 = QEC_data_dict['y_err_01']
#         y_err_10 = QEC_data_dict['y_err_10']
#         y_err_11 = QEC_data_dict['y_err_11']

#         p_00 = QEC_data_dict['p00']
#         p_01 = QEC_data_dict['p01']
#         p_10 = QEC_data_dict['p10']
#         p_11 = QEC_data_dict['p11']

        
#         ax1.errorbar(x,y,yerr=y_err,color = color[RO], label = 'QEC, decode to Qubit '+str(RO+1))

#         if plot_no_correct == True:
#             y_no_corr = undo_correction_single_state_RO_error_sign(run = 1, no_error = no_error,state = state,RO = RO,error_sign = 1,sweep_time=True)
#             dataset_dict_full[RO]['y_no_corr'] = y_no_corr
#             ax1.errorbar(x,y_no_corr,yerr = y_err, color = color[RO],ls = '-.', label = 'undo correction, Q'+str(RO+1))



#         ax2.errorbar(x,y_00,yerr=y_err_00,color = 'c', label = 'y_00' )
#         ax2.errorbar(x,y_01,yerr=y_err_01,color = 'k', label = 'y_01' )
#         ax2.errorbar(x,y_10,yerr=y_err_10,color = 'm', label = 'y_10' )
#         ax2.errorbar(x,y_11,yerr=y_err_11,color = 'b', label = 'y_11' )
#         ax2.set_ylim(-1.1,1.1)
#         ax2.set_xlim(-1e-3,35e-3)
#         ax2.legend()
#         ax2.set_title('sweep_time_error_syn_'+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC_PS')
#         ax2.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
#         ax2.set_xlabel('time (s)')
#         ax2.set_ylabel('Contrast')

#         ax3.set_title(str(folder)+'/'+ '\n probabilities')
#         ax3.plot(x,p_00, 'c', label = 'p00')
#         ax3.plot(x,p_01, 'k', label = 'p01')
#         ax3.plot(x,p_10, 'm', label = 'p10')
#         ax3.plot(x,p_11, 'b', label = 'p11')
#         ax3.set_xlim(-1e-3,35e-3)
#         # ax3.legend()
#         ax3.set_xlabel('time (s)')
#         ax3.set_ylabel('outcome probability')  
#         ax3.set_title('sweep_time_error_syn_'+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC_probs')                
        
#         if add_encode == True:
#             no_QEC_data_dict[RO] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = state,RO = RO, load_set = load_set,error_sign = 0)

#             x = no_QEC_data_dict[RO]['x']
#             y = no_QEC_data_dict[RO]['y']
#             y_err = no_QEC_data_dict[RO]['y_err']
#             ax1.errorbar(x,y,yerr=y_err,color = color[RO],ls = '--', label = 'Encoding, sweep time, Q' + str(RO+1) )

#         if add_single == True:
#             QEC_single_data_dict[RO] =  single_qubit_no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,state = state,sweep_time = True, error_sign = -1, Qubit = RO+1, load_set = True)    
#             x_single = QEC_single_data_dict[RO]['x']
#             y_single = QEC_single_data_dict[RO]['y']
#             y_single_err = QEC_single_data_dict[RO]['y_err']
          
#             ax1.errorbar(x_single,y_single,yerr=y_single_err,color = color[RO],ls = ':', label = 'Single Qubit, sweep time, Q' + str(RO+1) )


    
#     dataset_dict_full[6] = QEC_sweep_time_sum_states(RO  = 6)
    
#     y_toff = 1/2.*(dataset_dict_full[0]['y']+dataset_dict_full[1]['y']+dataset_dict_full[2]['y']-dataset_dict_full[6]['y'])
#     y_toff_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
#     x = dataset_dict_full[6]['x']+ np.ones(len(dataset_dict_full[6]['x']))*parity_time
#     ax1.errorbar(x,y_toff,yerr=y_toff_err,color = 'k', label = 'QEC+ toffoli' )
#     if plot_no_correct == True:
#         y_no_corr = undo_correction_single_state_RO_error_sign(run = 1, = state = state,RO = RO,error_sign = 1,sweep_time=True)
#         dataset_dict_full[6]['y_no_corr'] = y_no_corr
#         y_toff = 1/2.*(dataset_dict_full[0]['y_no_corr']+dataset_dict_full[1]['y_no_corr']+dataset_dict_full[2]['y_no_corr']-dataset_dict_full[6]['y_no_corr'])
#         y_toff_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
#         x = dataset_dict_full[6]['x']+ np.ones(len(dataset_dict_full[6]['x']))*parity_time
#         ax1.errorbar(x,y_toff,yerr=y_toff_err,color = 'k',ls = '-.', label = 'undo QEC+ toffoli' )    
#     if add_encode == True:
#         no_QEC_data_dict[6] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = state,RO = 6, load_set = load_set,error_sign = 0)
        
#         y_toff = 1/2.*(no_QEC_data_dict[0]['y']+no_QEC_data_dict[1]['y']+no_QEC_data_dict[2]['y']-no_QEC_data_dict[6]['y'])
#         y_toff_err = 1/2.*(no_QEC_data_dict[0]['y_err']**2+no_QEC_data_dict[1]['y_err']**2+no_QEC_data_dict[2]['y_err']**2+no_QEC_data_dict[6]['y_err']**2)**0.5
#         x = no_QEC_data_dict[0]['x']
#         ax1.errorbar(x,y_toff,yerr=y_toff_err,color = 'k',ls ='--', label = 'Encoding, sweep time, toff' )

#     ax1.set_ylim(-1.1,1.1)
#     ax1.set_xlim(-1e-3,35e-3)
#     ax1.set_title('sweep_time_error_syn_'+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC')
#     ax1.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
#     ax1.set_xlabel('time (s)')
#     ax1.set_ylabel('Contrast')
#     lgd = ax1.legend(loc = 2, bbox_to_anchor = (1,1))



#     if plot_no_correct == True:
#         try:
#             fig1.savefig(
#                 os.path.join(folder,'sweep_time_undo_correct_all'+'.pdf'))
#         except:
#             print 'Figure has not been saved.'
#     else:
#         try:
#             fig1.savefig(
#                 os.path.join(folder,'sweep_time_error_all'+'.pdf'),bbox_extra_artists = (lgd,),bbox_inches='tight')
#         except:
#             print 'Figure has not been saved.'

#     try:
#         fig2.savefig(
#             os.path.join(folder,'sweep_time_error_syn_'+no_error+'_state_'+state+'_ps'+'.pdf'))
#     except:
#         print 'Figure has not been saved.'

#     try:
#         fig3.savefig(
#             os.path.join(folder,'sweep_time_error_syn_'+no_error+'_state_'+state+'_probs'+'.pdf'))
#     except:
#         print 'Figure has not been saved.'

def no_QEC_plot_single_state_RO(state = 'Z',RO = 0, load_set = False, older_than = None):        
    

    QEC_data_dict =  no_QEC_data_single_state_RO(older_than = older_than,state = state,RO = RO, load_set = load_set)
    QEC_idle_data_dict =  no_QEC_data_single_state_RO(idle = True,older_than = older_than,state = state,RO = RO, load_set = load_set)

    
    folder  = r'D:\measuring\data\QEC_data\figs\Encoding'   


    x = QEC_data_dict['x']
    y = QEC_data_dict['y']
    y_idle = QEC_idle_data_dict['y']


    y_err = QEC_data_dict['y_err']
    y_idle_err = QEC_idle_data_dict['y_err']

    fig,ax = plt.subplots() 
    ax.errorbar(x,y,yerr=y_err,color = 'k', label = 'no idling' )
    ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'b', label = 'idling for QEC time' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('state_'+state+'_RO_'+str(RO)+'_noQEC')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')
    ax.legend()
    try:
        fig.savefig(
            os.path.join(folder,'no_QEC'+'_state_'+state+'_RO_'+str(RO)+'.png'))
    except:
        print 'Figure has not been saved.'

def no_QEC_sweep_time_plot_single_state(state = 'Z',load_set = False, older_than = None):        
    fig,ax = plt.subplots() 
    color = ['r','g','b']
    QEC_data_dict = {}
    QEC_single_data_dict = {}

    for  RO in [0,1,2]:
        QEC_data_dict[RO] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = state,RO = RO, load_set = load_set,error_sign = 0)
        QEC_single_data_dict[RO] =  single_qubit_no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,state = state,sweep_time = True, error_sign = -1, Qubit = RO+1, load_set = True)    
        
        folder  = r'D:\measuring\data\QEC_data\figs\Encoding'   


        x = QEC_data_dict[RO]['x']
        y = QEC_data_dict[RO]['y']
        y_err = QEC_data_dict[RO]['y_err']


        x_single = QEC_single_data_dict[RO]['x']
        y_single = QEC_single_data_dict[RO]['y']
        y_single_err = QEC_single_data_dict[RO]['y_err']

        if RO == 2:
            extra_time = 2*(4.996e-6*34 +11.312e-6*48) # 2*(13.616e-6*34)+116e-6
        elif RO == 0:
             extra_time = 2*(13.616e-6*34+11.312e-6*48)# 2*(4.996e-6*34)+116e-6
        elif RO == 1:
            extra_time =  2*(4.996e-6*34 +13.616e-6*34)#2*(11.312e-6*48)+116e-6

        x_single = x_single #+ np.ones(len(x_single))* extra_time

        ax.errorbar(x,y,yerr=y_err,color = color[RO], label = 'Encoding, sweep time, Q' + str(RO+1) )
        ax.errorbar(x_single,y_single,yerr=y_single_err,color = color[RO],ls = ':', label = 'Single Qubit, sweep time, Q' + str(RO+1) )
    
    QEC_data_dict[6] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = state,RO = 6, load_set = load_set,error_sign = 0)
    
    y_toff = 1/2.*(QEC_data_dict[0]['y']+QEC_data_dict[1]['y']+QEC_data_dict[2]['y']-QEC_data_dict[6]['y'])
    y_toff_err = 1/2.*(QEC_data_dict[0]['y_err']**2+QEC_data_dict[1]['y_err']**2+QEC_data_dict[2]['y_err']**2+QEC_data_dict[6]['y_err']**2)**0.5

    ax.errorbar(x,y_toff,yerr=y_toff_err,color = 'k', label = 'Encoding, sweep time, toff' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.001,0.035)
    # ax.set_title('state_'+state+'_RO_'+str(RO)+'_noQEC_sweep_time')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('time (s)')
    ax.set_ylabel('Contrast')
    if state == 'Z':
        location = 4
    elif state == 'mZ':
        location =1


    legend = ax.legend(loc = location)

    ax.set_title('no_QEC_sweep_time_single_vs_encode'+'_state_'+state)
    for label in legend.get_texts():
        label.set_fontsize('x-small')

    try:
        fig.savefig(
            os.path.join(folder,'no_QEC_sweep_time_single_vs_encode'+'_state_'+state+'.png'))
    except:
        print 'Figure has not been saved.'

def single_Qubit_no_QEC_plot_single_state_RO(state = 'Z',Qubit = 1, load_set = False, older_than = None):        
    

    QEC_data_dict =  single_qubit_no_QEC_data_single_state_RO(older_than = older_than,state = state,Qubit = Qubit, load_set = load_set)

    
    folder  = r'D:\measuring\data\QEC_data\figs\Encoding'   


    x = QEC_data_dict['x']
    y = QEC_data_dict['y']


    y_err = QEC_data_dict['y_err']

    fig,ax = plt.subplots() 
    ax.errorbar(x,y,yerr=y_err,color = 'k' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('state_'+state+'Qubit_'+str(Qubit)+'_single_qubit')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')
    try:
        fig.savefig(
            os.path.join(folder,'single_qubit'+'_state_'+state+'_Qubit_'+str(Qubit)+'.png'))
    except:
        print 'Figure has not been saved.'

def single_Qubit_sweep_time_no_QEC_plot_single_state(state = 'Z', load_set = False, older_than = None):        
    fig,ax = plt.subplots() 
    color = ['r','g','b']
    for Qubit in [1,2,3]:
        
        QEC_data_dict =  single_qubit_no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,state = state,sweep_time = True, error_sign = -1, Qubit = Qubit, load_set = load_set)

        
        folder  = r'D:\measuring\data\QEC_data\figs\Encoding'   


        x = QEC_data_dict['x']

        y = QEC_data_dict['y']


        y_err = QEC_data_dict['y_err']


        ax.errorbar(x,y,yerr=y_err,color = color[Qubit-1] , label = 'Qubit_'+str(Qubit))
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(0,35e-3)
    ax.set_title('state_'+state+'_single_qubit, sweep time')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Contrast')
    plt.legend()
    try:
        fig.savefig(
            os.path.join(folder,'single_qubit'+'_state_'+state+'_sweep_time'+'.png'))
    except:
        print 'Figure has not been saved.'

def plot_only_one(state = 'Z'):
    sign = 1
    if state == 'Y' or state == 'mZ':
        sign = -1
    for Qubit in [1,2,3]:
        if state == 'X' or state == 'mX':
            RO_list = [6]
        if state == 'Y' or state == 'mY':
            RO_list = [4,5,6]
        if state == 'Z' or state == 'mZ':
            RO_list = [0,1,2]
        folder  = r'D:\measuring\data\QEC_data\figs'

        for ii,RO in enumerate(RO_list):

            fig,ax = plt.subplots()
            ax.set_title('Qubit '+ str(Qubit) +' state '+state)

            QEC_data_dict = QEC_state_sum_all(state = state, RO = RO)
            y = QEC_data_dict['y']
            y_err = QEC_data_dict['y_err']
            x = QEC_data_dict['x']
            ax.errorbar(x,sign*y,yerr=y_err,color = 'r', label =  'QEC')

            no_QEC_data_dict =  no_QEC_data_single_state_RO(state = state,RO = RO, load_set = True)
            y = no_QEC_data_dict['y']
            y_err = no_QEC_data_dict['y_err']
            x = no_QEC_data_dict['x']
            ax.errorbar(x,sign*y,yerr=y_err,color = 'b',ls =':', label =  'Encoding')
            
            toff_process_dict = no_QEC_toffoli_fids()
            y = toff_process_dict['toff_'+state+'y']
            y_err = toff_process_dict['toff_'+state+'y_err']
            x = toff_process_dict['x']
            ax.errorbar(x,y,yerr=y_err,color = 'b',ls = '--', label =  'Encode toffoli')

            no_QEC_idle_data_dict =  no_QEC_data_single_state_RO(idle = True,state = state,RO = RO, load_set = True)
            y = no_QEC_idle_data_dict['y']
            y_err = no_QEC_idle_data_dict['y_err']
            x = no_QEC_idle_data_dict['x']
            ax.errorbar(x,sign*y,yerr=y_err,color = 'k',ls =':', label =  'Encode+idle')
            
            toff_process_dict_idle = no_QEC_toffoli_fids(idle = True)
            y = toff_process_dict_idle['toff_'+state+'y']
            y_err = toff_process_dict_idle['toff_'+state+'y_err']
            x = toff_process_dict_idle['x']
            ax.errorbar(x,y,yerr=y_err,color = 'k', ls = '--', label =  'Encode toffoli +idle')
            
            single_no_QEC_data_dict =  single_qubit_no_QEC_data_single_state_RO(state = state,Qubit = Qubit, load_set = True)
            y = single_no_QEC_data_dict['y']
            y_err = single_no_QEC_data_dict['y_err']
            x = single_no_QEC_data_dict['x']
            ax.errorbar(x,sign*y,yerr=y_err,color = 'g',ls =':', label =  'single qubit')

            
            ax.set_ylim(-1.1,1.1)
            ax.set_xlim(-0.1,1.1)
                   
            ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
            # ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = 'b')
            ax.vlines([0.5],-1.1,1.1,linestyles='dotted', color = 'b')
            ax.set_xlabel('error probability')
            ax.set_ylabel('Contrast')
            plt.legend()
            try:
                fig.savefig(
                    os.path.join(folder,'contrast_summed_qubit_'+str(Qubit)+'_state_'+ state+'.png'))
            except:
                print 'Figure has not been saved.'


''' Calculate and plot process fidelities '''


def QEC_process_fids(run = 1 ,no_error = '00'):

    dataset_dict = {}
    for state in ['Z','mZ','Y','mY', 'X','mX']:
        dataset_dict[state] = {}
        if state == 'X' or state == 'mX':
            RO_list = [6]
        if state == 'Y' or state == 'mY':
            RO_list = [4,5,6]
        if state == 'Z' or state == 'mZ':
            RO_list = [0,1,2]
        for RO in RO_list:
            dataset_dict[state]['Tomo_'+str(RO)] = {}

            dataset_dict[state]['Tomo_'+str(RO)] = QEC_sum_data_single_state_RO(run = run, no_error = no_error,state = state,RO = RO)
            dataset_dict[state]['Tomo_'+str(RO)]['y_new'] = {}
            dataset_dict[state]['Tomo_'+str(RO)]['y_new'] = undo_correction_single_state_RO(run = run, no_error = no_error,state = state,RO = RO)
        # print dataset_dict
    process_dict = {}

    y_list = ['y','y_00','y_01','y_10','y_11','y_new']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11','y_err']

    for v in range(6):
        # print v
        process_dict['dec_1_'+y_list[v]] = {}
        process_dict['dec_2_'+y_list[v]] = {}
        process_dict['dec_3_'+y_list[v]] = {}
        # process_dict['dec_toff_'+y_list[v]] = {}
        process_dict['dec_avg_'+y_list[v]] = {}
        process_dict['dec_1_'+y_err_list[v]] = {}
        process_dict['dec_2_'+y_err_list[v]] = {}
        process_dict['dec_3_'+y_err_list[v]] = {}
        # process_dict['toff_'+y_err_list[v]] = {}
        process_dict['dec_avg_'+y_err_list[v]] = {}

        # print dataset_dict['Z']
        process_dict['dec_1_'+y_list[v]] = 1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(0)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(0)][y_list[v]]
                        - dataset_dict['Y']['Tomo_'+str(5)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(5)][y_list[v]]
                        + dataset_dict['X']['Tomo_'+str(6)][y_list[v]] - dataset_dict['mX']['Tomo_'+str(6)][y_list[v]])

        process_dict['dec_1_'+y_err_list[v]] = 1/8.*(dataset_dict['Z']['Tomo_'+str(0)][y_err_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(0)][y_err_list[v]]**2 
                        + dataset_dict['Y']['Tomo_'+str(5)][y_err_list[v]]**2 +  dataset_dict['mY']['Tomo_'+str(5)][y_err_list[v]]**2 
                        + dataset_dict['X']['Tomo_'+str(6)][y_err_list[v]]**2 + dataset_dict['mX']['Tomo_'+str(6)][y_err_list[v]]**2 )**0.5

        process_dict['dec_2_'+y_list[v]] = 1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(1)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(1)][y_list[v]]
                - dataset_dict['Y']['Tomo_'+str(6)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(6)][y_list[v]]
                + dataset_dict['X']['Tomo_'+str(6)][y_list[v]] - dataset_dict['mX']['Tomo_'+str(6)][y_list[v]])

        process_dict['dec_2_'+y_err_list[v]] = 1/8.*(dataset_dict['Z']['Tomo_'+str(1)][y_err_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(1)][y_err_list[v]]**2 
                        + dataset_dict['Y']['Tomo_'+str(6)][y_err_list[v]]**2 +  dataset_dict['mY']['Tomo_'+str(6)][y_err_list[v]]**2 
                        + dataset_dict['X']['Tomo_'+str(6)][y_err_list[v]]**2 + dataset_dict['mX']['Tomo_'+str(6)][y_err_list[v]]**2 )**0.5

        process_dict['dec_3_'+y_list[v]] = 1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(2)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(2)][y_list[v]]
                        - dataset_dict['Y']['Tomo_'+str(4)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(4)][y_list[v]]
                        + dataset_dict['X']['Tomo_'+str(6)][y_list[v]] - dataset_dict['mX']['Tomo_'+str(6)][y_list[v]])

        process_dict['dec_3_'+y_err_list[v]] = 1/8.*(dataset_dict['Z']['Tomo_'+str(2)][y_err_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(2)][y_err_list[v]]**2 
                        + dataset_dict['Y']['Tomo_'+str(4)][y_err_list[v]]**2 +  dataset_dict['mY']['Tomo_'+str(4)][y_err_list[v]]**2 
                        + dataset_dict['X']['Tomo_'+str(6)][y_err_list[v]]**2 + dataset_dict['mX']['Tomo_'+str(6)][y_err_list[v]]**2 )**0.5


        process_dict['dec_avg_'+y_list[v]] = 1/3.* (1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(0)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(0)][y_list[v]]
                        - dataset_dict['Y']['Tomo_'+str(5)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(5)][y_list[v]]
                        + dataset_dict['X']['Tomo_'+str(6)][y_list[v]] - dataset_dict['mX']['Tomo_'+str(6)][y_list[v]])+
                            1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(1)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(1)][y_list[v]]
                                            - dataset_dict['Y']['Tomo_'+str(6)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(6)][y_list[v]]
                                            + dataset_dict['X']['Tomo_'+str(6)][y_list[v]] - dataset_dict['mX']['Tomo_'+str(6)][y_list[v]])+
                        1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(2)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(2)][y_list[v]]
                                                - dataset_dict['Y']['Tomo_'+str(4)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(4)][y_list[v]]
                                                + dataset_dict['X']['Tomo_'+str(6)][y_list[v]] - dataset_dict['mX']['Tomo_'+str(6)][y_list[v]])      )                      

        process_dict['dec_avg_'+y_err_list[v]] = 1/3.*1/8.*((dataset_dict['Z']['Tomo_'+str(0)][y_err_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(0)][y_err_list[v]]**2 
                        + dataset_dict['Y']['Tomo_'+str(5)][y_err_list[v]]**2 +  dataset_dict['mY']['Tomo_'+str(5)][y_err_list[v]]**2 
                        + dataset_dict['X']['Tomo_'+str(6)][y_err_list[v]]**2 + dataset_dict['mX']['Tomo_'+str(6)][y_err_list[v]]**2 ) +
                                    (dataset_dict['Z']['Tomo_'+str(1)][y_err_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(1)][y_err_list[v]]**2 
                                            + dataset_dict['Y']['Tomo_'+str(6)][y_err_list[v]]**2 +  dataset_dict['mY']['Tomo_'+str(6)][y_err_list[v]]**2 
                                            + dataset_dict['X']['Tomo_'+str(6)][y_err_list[v]]**2 + dataset_dict['mX']['Tomo_'+str(6)][y_err_list[v]]**2 )+
                                    (dataset_dict['Z']['Tomo_'+str(2)][y_err_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(2)][y_err_list[v]]**2 
                                            + dataset_dict['Y']['Tomo_'+str(4)][y_err_list[v]]**2 +  dataset_dict['mY']['Tomo_'+str(4)][y_err_list[v]]**2 
                                            + dataset_dict['X']['Tomo_'+str(6)][y_err_list[v]]**2 + dataset_dict['mX']['Tomo_'+str(6)][y_err_list[v]]**2 ))**0.5

    process_dict['x'] = dataset_dict['Z']['Tomo_'+str(0)]['x']
    return process_dict


def QEC_process_fids_sum_runs(run_list = [1,2],no_error = '00'):
    dec_list = ['dec_1_','dec_2_','dec_3_','dec_avg_']
    y_list = ['y','y_00','y_01','y_10','y_11','y_new']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11','y_err']


    for i in range(len(run_list)):
        process_dict_single = {}
        process_dict_single = QEC_process_fids(run = run_list[i], no_error = no_error)
        # print process_dict_single
        if i ==0:
            process_dict = {}
            # print process_dict_single[item+ yy]
            for ii, item in enumerate(dec_list):
                for k, yy in enumerate(y_list):
                    process_dict[item+ yy] = 1/float(len(run_list))*process_dict_single[item+ yy]
                    process_dict['temp'+item+ y_err_list[k]] = process_dict_single[item+ y_err_list[k]]**2
        else:
            for ii, item in enumerate(dec_list):
                for k, yy in enumerate(y_list):
                    process_dict[item+ yy] += 1/float(len(run_list))*process_dict_single[item+ yy]
                    process_dict['temp'+item+ y_err_list[k]] += process_dict_single[item+ y_err_list[k]]**2
        for ii, item in enumerate(dec_list):                    
            for k, yy in enumerate(y_list):
                process_dict[item+ y_err_list[k]] = 1/float(len(run_list))*(process_dict['temp'+item+ y_err_list[k]])**0.5
    
    process_dict['x'] = process_dict_single['x']
    
    return process_dict

def QEC_process_fids_sum_all(run_list_00 = [1,2,3],run_list_01 = [1,2,3],run_list_10 = [2],run_list_11 = [3]):
    dec_list = ['dec_1_','dec_2_','dec_3_','dec_avg_']
    y_list = ['y','y_00','y_01','y_10','y_11','y_new']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11','y_err']

    for i, run in enumerate(['00','01','10','11']):
            if run == '00':
                run_list = run_list_00
            elif run == '01':
                run_list = run_list_01
            elif run == '10':
                run_list = run_list_10
            elif run == '11':
                run_list = run_list_11

            process_dict_single= {}
            process_dict_single = QEC_process_fids_sum_runs(run_list = run_list,no_error = run)
            if i ==0:
                process_dict = {}
                for ii, item in enumerate(dec_list):
                    for k, yy in enumerate(y_list):
                        process_dict[item+ yy] = 1/4.*process_dict_single[item+ yy]
                        process_dict['temp'+item+ y_err_list[k]] = process_dict_single[item+ y_err_list[k]]**2
            else:
                for ii, item in enumerate(dec_list):
                    for k, yy in enumerate(y_list):
                        process_dict[item+ yy] += 1/4.*process_dict_single[item+ yy]
                        process_dict['temp'+item+ y_err_list[k]] += process_dict_single[item+ y_err_list[k]]**2
            for ii, item in enumerate(dec_list):                    
                for k, yy in enumerate(y_list):
                    process_dict[item+ y_err_list[k]] = 1/4.*(process_dict['temp'+item+ y_err_list[k]])**0.5
    
    process_dict['x'] = process_dict_single['x']
    
    return process_dict

def no_QEC_process_fids(idle = False):

    dataset_dict = {}
    for state in ['Z','mZ','Y','mY', 'X','mX']:
        dataset_dict[state] = {}
        if state == 'X' or state == 'mX':
            RO_list = [6]
        if state == 'Y' or state == 'mY':
            RO_list = [3,4,5,6]
        if state == 'Z' or state == 'mZ':
            RO_list = [6,0,1,2]
        for RO in RO_list:
            dataset_dict[state]['Tomo_'+str(RO)] = {}

            dataset_dict[state]['Tomo_'+str(RO)] = no_QEC_data_single_state_RO(idle = idle,state = state,RO = RO,load_set = True)
        # print dataset_dict
    process_dict = {}

    process_dict['dec_1_'+'y'] = {}
    process_dict['dec_2_'+'y'] = {}
    process_dict['dec_3_'+'y'] = {}
    # process_dict['dec_toff_'+'y'] = {}
    process_dict['dec_avg_'+'y'] = {}
    process_dict['dec_1_'+'y_err'] = {}
    process_dict['dec_2_'+'y_err'] = {}
    process_dict['dec_3_'+'y_err'] = {}
    # process_dict['toff_'+'y_err'] = {}
    process_dict['dec_avg_'+'y_err'] = {}

    # print dataset_dict['Z']

    process_dict['dec_1_'+'y'] = 1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(0)]['y'] - dataset_dict['mZ']['Tomo_'+str(0)]['y']
                    - dataset_dict['Y']['Tomo_'+str(5)]['y'] + dataset_dict['mY']['Tomo_'+str(5)]['y']
                    + dataset_dict['X']['Tomo_'+str(6)]['y'] - dataset_dict['mX']['Tomo_'+str(6)]['y'])

    process_dict['dec_1_'+'y_err'] = 1/8.*(dataset_dict['Z']['Tomo_'+str(0)]['y_err']**2 + dataset_dict['mZ']['Tomo_'+str(0)]['y_err']**2 
                    + dataset_dict['Y']['Tomo_'+str(5)]['y_err']**2 +  dataset_dict['mY']['Tomo_'+str(5)]['y_err']**2 
                    + dataset_dict['X']['Tomo_'+str(6)]['y_err']**2 + dataset_dict['mX']['Tomo_'+str(6)]['y_err']**2 )**0.5

    process_dict['dec_2_'+'y'] = 1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(1)]['y'] - dataset_dict['mZ']['Tomo_'+str(1)]['y']
            - dataset_dict['Y']['Tomo_'+str(6)]['y'] + dataset_dict['mY']['Tomo_'+str(6)]['y']
            + dataset_dict['X']['Tomo_'+str(6)]['y'] - dataset_dict['mX']['Tomo_'+str(6)]['y'])

    process_dict['dec_2_'+'y_err'] = 1/8.*(dataset_dict['Z']['Tomo_'+str(1)]['y_err']**2 + dataset_dict['mZ']['Tomo_'+str(1)]['y_err']**2 
                    + dataset_dict['Y']['Tomo_'+str(6)]['y_err']**2 +  dataset_dict['mY']['Tomo_'+str(6)]['y_err']**2 
                    + dataset_dict['X']['Tomo_'+str(6)]['y_err']**2 + dataset_dict['mX']['Tomo_'+str(6)]['y_err']**2 )**0.5

    process_dict['dec_3_'+'y'] = 1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(2)]['y'] - dataset_dict['mZ']['Tomo_'+str(2)]['y']
                    - dataset_dict['Y']['Tomo_'+str(4)]['y'] + dataset_dict['mY']['Tomo_'+str(4)]['y']
                    + dataset_dict['X']['Tomo_'+str(6)]['y'] - dataset_dict['mX']['Tomo_'+str(6)]['y'])

    process_dict['dec_3_'+'y_err'] = 1/8.*(dataset_dict['Z']['Tomo_'+str(2)]['y_err']**2 + dataset_dict['mZ']['Tomo_'+str(2)]['y_err']**2 
                    + dataset_dict['Y']['Tomo_'+str(4)]['y_err']**2 +  dataset_dict['mY']['Tomo_'+str(4)]['y_err']**2 
                    + dataset_dict['X']['Tomo_'+str(6)]['y_err']**2 + dataset_dict['mX']['Tomo_'+str(6)]['y_err']**2 )**0.5


    process_dict['dec_avg_'+'y'] = 1/3.* (1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(0)]['y'] - dataset_dict['mZ']['Tomo_'+str(0)]['y']
                    - dataset_dict['Y']['Tomo_'+str(5)]['y'] + dataset_dict['mY']['Tomo_'+str(5)]['y']
                    + dataset_dict['X']['Tomo_'+str(6)]['y'] - dataset_dict['mX']['Tomo_'+str(6)]['y'])+
                        1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(1)]['y'] - dataset_dict['mZ']['Tomo_'+str(1)]['y']
                                        - dataset_dict['Y']['Tomo_'+str(6)]['y'] + dataset_dict['mY']['Tomo_'+str(6)]['y']
                                        + dataset_dict['X']['Tomo_'+str(6)]['y'] - dataset_dict['mX']['Tomo_'+str(6)]['y'])+
                    1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(2)]['y'] - dataset_dict['mZ']['Tomo_'+str(2)]['y']
                                            - dataset_dict['Y']['Tomo_'+str(4)]['y'] + dataset_dict['mY']['Tomo_'+str(4)]['y']
                                            + dataset_dict['X']['Tomo_'+str(6)]['y'] - dataset_dict['mX']['Tomo_'+str(6)]['y'])      )                      

    process_dict['dec_avg_'+'y_err'] = 1/3.*(1/8.*(dataset_dict['Z']['Tomo_'+str(0)]['y_err']**2 + dataset_dict['mZ']['Tomo_'+str(0)]['y_err']**2 
                    + dataset_dict['Y']['Tomo_'+str(5)]['y_err']**2 +  dataset_dict['mY']['Tomo_'+str(5)]['y_err']**2 
                    + dataset_dict['X']['Tomo_'+str(6)]['y_err']**2 + dataset_dict['mX']['Tomo_'+str(6)]['y_err']**2 ) +
                                1/8.*(dataset_dict['Z']['Tomo_'+str(1)]['y_err']**2 + dataset_dict['mZ']['Tomo_'+str(1)]['y_err']**2 
                                        + dataset_dict['Y']['Tomo_'+str(6)]['y_err']**2 +  dataset_dict['mY']['Tomo_'+str(6)]['y_err']**2 
                                        + dataset_dict['X']['Tomo_'+str(6)]['y_err']**2 + dataset_dict['mX']['Tomo_'+str(6)]['y_err']**2 )+
                                1/8.*(dataset_dict['Z']['Tomo_'+str(2)]['y_err']**2 + dataset_dict['mZ']['Tomo_'+str(2)]['y_err']**2 
                                        + dataset_dict['Y']['Tomo_'+str(4)]['y_err']**2 +  dataset_dict['mY']['Tomo_'+str(4)]['y_err']**2 
                                        + dataset_dict['X']['Tomo_'+str(6)]['y_err']**2 + dataset_dict['mX']['Tomo_'+str(6)]['y_err']**2 ))**0.5


    process_dict['x'] = dataset_dict['Z']['Tomo_'+str(0)]['x']
    return process_dict

def no_QEC_toffoli_fids(idle = False):

    dataset_dict = {}
    toff_dict = {}

    for state in ['Z','mZ','Y','mY', 'X','mX']:
        dataset_dict[state] = {}
        if state == 'X' or state == 'mX':
            RO_list = [6]
        if state == 'Y' or state == 'mY':
            RO_list = [3,4,5,6]
        if state == 'Z' or state == 'mZ':
            RO_list = [6,0,1,2]
        for RO in RO_list:
            dataset_dict[state]['Tomo_'+str(RO)] = {}

            dataset_dict[state]['Tomo_'+str(RO)] = no_QEC_data_single_state_RO(idle = idle,state = state,RO = RO,load_set = True)
        # print dataset_dict


        toff_dict['toff_'+state+'y'] = {}

        toff_dict['toff_'+state+'y_err'] = {}

        if state == 'X':
            toff_dict['toff_'+state+'y'] = (dataset_dict[state]['Tomo_'+str(6)]['y'])
            toff_dict['toff_'+state+'y_err'] = (dataset_dict[state]['Tomo_'+str(6)]['y_err'])
        if state == 'mX':
            toff_dict['toff_'+state+'y'] = (dataset_dict[state]['Tomo_'+str(6)]['y'])
            toff_dict['toff_'+state+'y_err'] = (dataset_dict[state]['Tomo_'+str(6)]['y_err'])
        if state == 'Y':
            toff_dict['toff_'+state+'y'] = 1/2.*(-dataset_dict[state]['Tomo_'+str(3)]['y']-dataset_dict[state]['Tomo_'+str(4)]['y']-dataset_dict[state]['Tomo_'+str(5)]['y']-dataset_dict[state]['Tomo_'+str(6)]['y'])
            toff_dict['toff_'+state+'y_err'] = 1/2.*(+dataset_dict[state]['Tomo_'+str(3)]['y_err']**2+dataset_dict[state]['Tomo_'+str(4)]['y_err']**2+dataset_dict[state]['Tomo_'+str(5)]['y_err']**2+dataset_dict[state]['Tomo_'+str(6)]['y_err']**2)**0.5
        if state == 'mY':
            toff_dict['toff_'+state+'y'] = 1/2.*(-dataset_dict[state]['Tomo_'+str(3)]['y']-dataset_dict[state]['Tomo_'+str(4)]['y']-dataset_dict[state]['Tomo_'+str(5)]['y']-dataset_dict[state]['Tomo_'+str(6)]['y'])
            toff_dict['toff_'+state+'y_err'] = 1/2.*(+dataset_dict[state]['Tomo_'+str(3)]['y_err']**2+dataset_dict[state]['Tomo_'+str(4)]['y_err']**2+dataset_dict[state]['Tomo_'+str(5)]['y_err']**2+dataset_dict[state]['Tomo_'+str(6)]['y_err']**2)**0.5
        if state == 'Z':
            toff_dict['toff_'+state+'y'] = 1/2.*(-dataset_dict[state]['Tomo_'+str(6)]['y']+dataset_dict[state]['Tomo_'+str(0)]['y']+dataset_dict[state]['Tomo_'+str(1)]['y']+dataset_dict[state]['Tomo_'+str(2)]['y'])
            toff_dict['toff_'+state+'y_err'] = 1/2.*(+dataset_dict[state]['Tomo_'+str(6)]['y_err']**2+dataset_dict[state]['Tomo_'+str(0)]['y_err']**2+dataset_dict[state]['Tomo_'+str(1)]['y_err']**2+dataset_dict[state]['Tomo_'+str(2)]['y_err']**2)**0.5
        if state == 'mZ':
            toff_dict['toff_'+state+'y'] = 1/2.*(-dataset_dict[state]['Tomo_'+str(6)]['y']+dataset_dict[state]['Tomo_'+str(0)]['y']+dataset_dict[state]['Tomo_'+str(1)]['y']+dataset_dict[state]['Tomo_'+str(2)]['y'])
            toff_dict['toff_'+state+'y_err'] = 1/2.*(+dataset_dict[state]['Tomo_'+str(6)]['y_err']**2+dataset_dict[state]['Tomo_'+str(0)]['y_err']**2+dataset_dict[state]['Tomo_'+str(1)]['y_err']**2+dataset_dict[state]['Tomo_'+str(2)]['y_err']**2)**0.5

    toff_dict['toff_process_y'] = 1/4.+1/8.*(toff_dict['toff_'+'Z'+'y']-toff_dict['toff_'+'mZ'+'y']+toff_dict['toff_'+'Y'+'y']-toff_dict['toff_'+'mY'+'y']+toff_dict['toff_'+'X'+'y']-toff_dict['toff_'+'mX'+'y'])
    toff_dict['toff_process_y_err'] = 1/8.*(toff_dict['toff_'+'Z'+'y_err']**2+toff_dict['toff_'+'mZ'+'y_err']**2+toff_dict['toff_'+'Y'+'y_err']**2+toff_dict['toff_'+'mY'+'y_err']**2+toff_dict['toff_'+'X'+'y_err']**2+toff_dict['toff_'+'mX'+'y_err']**2)**0.5
    # print toff_dict['toff_process_y'][0]

    toff_dict['x'] = dataset_dict['Z']['Tomo_'+str(0)]['x']
    return toff_dict

def single_Qubit_no_QEC_process_fids():
    dataset_dict = {}
    for Qubit in [1,2,3]:
        dataset_dict['Q'+str(Qubit)] = {}
        for state in ['Z','mZ','Y','mY', 'X','mX']:
            dataset_dict['Q'+str(Qubit)][state] = {}
            
            if state == 'X' or state == 'mX':
                RO = 2+(Qubit-1)*3
            elif state == 'Y' or state == 'mY':
                RO = 1+(Qubit-1)*3
            elif state == 'Z' or state == 'mZ': 
                RO = 0+(Qubit-1)*3

            dataset_dict['Q'+str(Qubit)][state] = single_qubit_no_QEC_data_single_state_RO(state = state,Qubit = Qubit, load_set = True)
        # print dataset_dict
    process_dict = {}

    process_dict['dec_1_'+'y'] = {}
    process_dict['dec_2_'+'y'] = {}
    process_dict['dec_3_'+'y'] = {}
    # process_dict['dec_toff_'+'y'] = {}
    process_dict['dec_avg_'+'y'] = {}
    process_dict['dec_1_'+'y_err'] = {}
    process_dict['dec_2_'+'y_err'] = {}
    process_dict['dec_3_'+'y_err'] = {}
    # process_dict['toff_'+'y_err'] = {}
    process_dict['dec_avg_'+'y_err'] = {}

    # print dataset_dict['Z']

    process_dict['dec_1_'+'y'] = 1/4. + 1/8.*(dataset_dict['Q1']['Z']['y'] - dataset_dict['Q1']['mZ']['y']
                    - dataset_dict['Q1']['Y']['y'] + dataset_dict['Q1']['mY']['y']
                    + dataset_dict['Q1']['X']['y'] - dataset_dict['Q1']['mX']['y'])

    process_dict['dec_1_'+'y_err'] = 1/8.*(dataset_dict['Q1']['Z']['y_err']**2 + dataset_dict['Q1']['mZ']['y_err']**2 
                    + dataset_dict['Q1']['Y']['y_err']**2 +  dataset_dict['Q1']['mY']['y_err']**2 
                    + dataset_dict['Q1']['X']['y_err']**2 + dataset_dict['Q1']['mX']['y_err']**2 )**0.5

    process_dict['dec_2_'+'y'] = 1/4. + 1/8.*(dataset_dict['Q2']['Z']['y'] - dataset_dict['Q2']['mZ']['y']
                    - dataset_dict['Q2']['Y']['y'] + dataset_dict['Q2']['mY']['y']
                    + dataset_dict['Q2']['X']['y'] - dataset_dict['Q2']['mX']['y'])

    process_dict['dec_2_'+'y_err'] = 1/8.*(dataset_dict['Q2']['Z']['y_err']**2 + dataset_dict['Q2']['mZ']['y_err']**2 
                    + dataset_dict['Q2']['Y']['y_err']**2 +  dataset_dict['Q2']['mY']['y_err']**2 
                    + dataset_dict['Q2']['X']['y_err']**2 + dataset_dict['Q2']['mX']['y_err']**2 )**0.5

    process_dict['dec_3_'+'y'] = 1/4. + 1/8.*(dataset_dict['Q3']['Z']['y'] - dataset_dict['Q3']['mZ']['y']
                    - dataset_dict['Q3']['Y']['y'] + dataset_dict['Q3']['mY']['y']
                    + dataset_dict['Q3']['X']['y'] - dataset_dict['Q3']['mX']['y'])

    process_dict['dec_3_'+'y_err'] = 1/8.*(dataset_dict['Q3']['Z']['y_err']**2 + dataset_dict['Q3']['mZ']['y_err']**2 
                    + dataset_dict['Q3']['Y']['y_err']**2 +  dataset_dict['Q3']['mY']['y_err']**2 
                    + dataset_dict['Q3']['X']['y_err']**2 + dataset_dict['Q3']['mX']['y_err']**2 )**0.5


    process_dict['dec_avg_'+'y'] = 1/3.* (process_dict['dec_1_'+'y']+process_dict['dec_2_'+'y']+process_dict['dec_3_'+'y'])                  

    process_dict['dec_avg_'+'y_err'] = 1/3.*(process_dict['dec_1_'+'y_err']**2+process_dict['dec_2_'+'y_err']**2+process_dict['dec_3_'+'y_err']**2)**0.5


    process_dict['x'] = dataset_dict['Q1']['Z']['x']
    return process_dict

def QEC_plot_process_fids(run_list = [1],no_error = '00',append_no_QEC = False, append_undo_corr = False):

    process_dict = QEC_process_fids_sum_runs(run_list = run_list,no_error = no_error)

    x = process_dict['x']
    folder  = r'D:\measuring\data\QEC_data\figs\process fidelities'

    t_list = ['1','2','3','avg']
    color_list = ['c','r','b','g']

    if append_no_QEC == True:
       no_process_dict = no_QEC_process_fids() 


    fig,ax = plt.subplots()

    for i in range(4):
        y = process_dict['dec_'+t_list[i]+'_y']
        y_err = process_dict['dec_'+t_list[i]+'_y_err']
        if i == 3:
            x_g = [x[0],x[-1]]
            y_g = [y[0],y[-1]]  
            # ax.plot(x_g,y_g,color = 'k', linestyle = ':' )          
        ax.errorbar(x,y,yerr=y_err,color = color_list[i], label =  'decode to '+ t_list[i])
        if append_no_QEC == True:
            y = no_process_dict['dec_'+t_list[i]+'_y']
            y_err = no_process_dict['dec_'+t_list[i]+'_y_err']
            ax.errorbar(x,y,yerr=y_err,color = color_list[i], ls = ':', label =  'no QEC, decode to '+ t_list[i])
        # if append_undo_corr == True:
        #     y = process_dict['dec_'+t_list[i]+'_y_new']
        #     y_err = process_dict['dec_'+t_list[i]+'_y_err']
        #     ax.errorbar(x,y,yerr=y_err,color = color_list[i], ls = '--', label =  'undo correction, decode to '+ t_list[i])            
    ax.set_ylim(-0.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('error_syn_'+no_error+'_run_'+str(run_list)+'_process_fidelity'+'.png')                
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = 'b')
    ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = 'b')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Process Fidelity')
    plt.legend()
    
    try:
        fig.savefig(
            os.path.join(folder,'error_syn_'+no_error+'_run_'+str(run_list)+'_process_fidelity'+'.png'))
    except:
        print 'Figure has not been saved.'

    if append_undo_corr == True:
            

        for i in range(4):
            fig,ax = plt.subplots()
            y = process_dict['dec_'+t_list[i]+'_y']
            y_err = process_dict['dec_'+t_list[i]+'_y_err']     
            ax.errorbar(x,y,yerr=y_err,color = color_list[i], label =  'decode to '+ t_list[i])
            
            y = no_process_dict['dec_'+t_list[i]+'_y']
            y_err = no_process_dict['dec_'+t_list[i]+'_y_err']
            ax.errorbar(x,y,yerr=y_err,color = color_list[i], ls = ':', label =  'no QEC, decode to '+ t_list[i])

            y = process_dict['dec_'+t_list[i]+'_y_new']
            y_err = process_dict['dec_'+t_list[i]+'_y_err']
            ax.errorbar(x,y,yerr=y_err,color = color_list[i], ls = '--', label =  'undo correction, decode to '+ t_list[i])            
            ax.set_ylim(-0.1,1.1)
            ax.set_xlim(-0.1,1.1)
            ax.set_title('error_syn_'+no_error+'_run_'+str(run_list)+'_process_fidelity dec to ' +t_list[i]+'.png')                
            ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
            ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = 'b')
            ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = 'b')
            ax.set_xlabel('error probability')
            ax.set_ylabel('Process Fidelity')
            plt.legend()
            
            try:
                fig.savefig(
                    os.path.join(folder,'error_syn_'+no_error+'_run_'+str(run_list)+'_process_fidelity dec to ' +t_list[i]+'.png'))
            except:
                print 'Figure has not been saved.'


def QEC_plot_process_fids_sum(append_no_QEC =True):

    process_dict = QEC_process_fids_sum_all()
    toff_process_dict = no_QEC_toffoli_fids()
    x = process_dict['x']
    folder  = r'D:\measuring\data\QEC_data\figs\process fidelities'

    t_list = ['1','2','3','avg']
    color_list = ['c','r','b','g']

    if append_no_QEC == True:
       no_process_dict = no_QEC_process_fids() 
       single_process_dict = single_Qubit_no_QEC_process_fids()
       process_dict_idle = no_QEC_process_fids(idle = True)

       toff_dict_idle = no_QEC_toffoli_fids(idle = True)

    
    fig,ax = plt.subplots()
    for i in range(4):

        y = process_dict['dec_'+t_list[i]+'_y']
        y_err = process_dict['dec_'+t_list[i]+'_y_err']

        x_fit, y_fit, p_err = fit_QEC_curve(x,y)

        if i == 3:
            x_g = [x[0],x[-1]]
            y_g = [y[0],y[-1]]  
            # ax.plot(x_g,y_g,color = 'k', linestyle = ':' )          
        ax.errorbar(x,y,yerr=y_err,color = color_list[i],ls = '', marker = 'o', ms = 8)
        ax.plot(x_fit, y_fit, color_list[i], lw=1, label =  'decode to '+ t_list[i]+ ', p_correction: '+str(int(p_err*100)/100.))
        
        if append_no_QEC == True:
            y = no_process_dict['dec_'+t_list[i]+'_y']
            y_err = no_process_dict['dec_'+t_list[i]+'_y_err']
            y_idle = process_dict_idle['dec_'+t_list[i]+'_y']
            y_idle_err = process_dict_idle['dec_'+t_list[i]+'_y_err']

            x_fit, y_fit, p_err = fit_QEC_curve(x,y)
            ax.errorbar(x,y,yerr=y_err,color = color_list[i], ls = '',marker = '*')
            ax.plot(x_fit, y_fit, color_list[i], lw=1, ls = ':',  label =  'no QEC, decode to '+ t_list[i]+ ', p_correction: '+str(int(p_err*100)/100.))
            
            x_fit, y_idle_fit, p_err = fit_QEC_curve(x,y_idle)
            ax.errorbar(x,y_idle,yerr=y_idle_err,color = color_list[i], ls = '', marker = '*')
            ax.plot(x_fit, y_idle_fit, color_list[i], lw=1, ls = '--', label =  'with idling, no QEC, decode to '+ t_list[i]+ ', p_correction: '+str(int(p_err*100)/100.))
            

            if i ==3:                
                y = toff_process_dict['toff_process_y']
                y_err = toff_process_dict['toff_process_y_err']
                x_fit, y_fit, p_err = fit_QEC_curve(x,y)
                ax.errorbar(x,y,yerr=y_err,color = 'k', ls = '',marker = '*')
                ax.plot(x_fit, y_fit, 'k', lw=1, ls = ':',label =  'toffoli decoded'+ ', p_correction: '+str(int(p_err*100)/100.))

                y_idle = toff_dict_idle['toff_process_y']
                y_idle_err = toff_dict_idle['toff_process_y_err']  
                x_fit, y_idle_fit, p_err = fit_QEC_curve(x,y_idle)
                ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'k',  ls = '',marker = '*')

                ax.plot(x_fit, y_idle_fit, 'k', lw=1, ls = '--', label =  'with idling, toffoli decoded'+ ', p_correction: '+str(int(p_err*100)/100.))      

    ax.set_ylim(-0.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('all_summed_process_fidelity dec to '+t_list[i]+'.png')                
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = 'b')
    ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = 'b')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Process Fidelity')
    lgd = ax.legend(loc = 2, bbox_to_anchor = (1,1))
    
    try:
        fig.savefig(
            os.path.join(folder,'all_summed_process_fidelity'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
    except:
        print 'Figure has not been saved.'


        

    for i in range(4):
        fig,ax = plt.subplots()
        y = process_dict['dec_'+t_list[i]+'_y']
        y_err = process_dict['dec_'+t_list[i]+'_y_err'] 
        x_fit, y_fit, p_err = fit_QEC_curve(x,y)    
        ax.plot(x_fit, y_fit, 'r', lw=1, label =  'decode to '+ t_list[i]+ ', p_correction: '+str(int(p_err*100)/100.))
        ax.errorbar(x,y,yerr=y_err,color = 'r',ls = '',marker = 'o')
        # x_g = [x[0],x[-1]]
        # y_g = [y[0],y[-1]]
        # ax.plot(x_g,y_g,color = 'k', linestyle = ':' )

        y = process_dict['dec_'+t_list[i]+'_y_new']
        y_err = process_dict['dec_'+t_list[i]+'_y_err']
        x_fit, y_fit, p_err = fit_QEC_curve(x,y)    
        ax.plot(x_fit, y_fit, 'r',ls = ':', lw=1, label =  'undo correction, decode to '+ t_list[i]+ ', p_correction: '+str(int(p_err*100)/100.))  
        ax.errorbar(x,y,yerr=y_err,color = 'r', ls = '',marker = '.')
        
        y = no_process_dict['dec_'+t_list[i]+'_y']
        y_err = no_process_dict['dec_'+t_list[i]+'_y_err']
        x_fit, y_fit, p_err = fit_QEC_curve(x,y)    
        ax.plot(x_fit, y_fit, 'b',ls = '--', lw=1, label =  'no QEC, decode to '+ t_list[i]+ ', p_correction: '+str(int(p_err*100)/100.))     
        ax.errorbar(x,y,yerr=y_err,color = 'b', ls = '',marker ='.')
        
        y = toff_process_dict['toff_process_y']
        y_err = toff_process_dict['toff_process_y_err']
        x_fit, y_fit, p_err = fit_QEC_curve(x,y)    
        ax.plot(x_fit, y_fit, 'b',ls = '-', lw=1,label =  'toffoli decoded'+ ', p_correction: '+str(int(p_err*100)/100.)) 
        ax.errorbar(x,y,yerr=y_err,color = 'b', ls = '',marker = '*')
    
        y_idle = process_dict_idle['dec_'+t_list[i]+'_y']
        y_idle_err = process_dict_idle['dec_'+t_list[i]+'_y_err']
        x_fit, y_fit, p_err = fit_QEC_curve(x,y_idle)    
        ax.plot(x_fit, y_fit, 'k',ls = '--', lw=1, label =  'with idling, no QEC, decode to '+ t_list[i]+ ', p_correction: '+str(int(p_err*100)/100.)) 
        ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'k', ls = '',marker = '.')

        y_idle = toff_dict_idle['toff_process_y']
        y_idle_err = toff_dict_idle['toff_process_y_err']  
        x_fit, y_fit, p_err = fit_QEC_curve(x,y_idle)    
        ax.plot(x_fit, y_fit, 'k',ls = '-', lw=1, label =  'with idling, toffoli decoded'+ ', p_correction: '+str(int(p_err*100)/100.))   
        ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'k',  ls = '', marker = '*')

        y = single_process_dict['dec_'+t_list[i]+'_y']
        y_err = process_dict['dec_'+t_list[i]+'_y_err']
        x_fit, y_fit, p_err = fit_QEC_curve(x,y)    
        ax.plot(x_fit, y_fit, 'g',ls = '--', lw=1, label =  'single Qubit: '+ t_list[i]+ ', p_correction: '+str(int(p_err*100)/100.))     
        ax.errorbar(x,y,yerr=y_err,color = 'g', ls = '',marker = 'o')

        ax.set_ylim(-0.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title('all_summed_process_fidelity dec to ' +t_list[i]+'.png')                
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = 'b')
        ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = 'b')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Process Fidelity')
        lgd = ax.legend(loc = 2, bbox_to_anchor = (1,1))
        
        try:
            fig.savefig(
                os.path.join(folder,'all_summed_process_fidelity dec to ' +t_list[i]+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
        except:
            print 'Figure has not been saved.'


def no_QEC_plot_toff_fids():
    toff_dict = no_QEC_toffoli_fids()
    toff_dict_idle = no_QEC_toffoli_fids(idle = True)

    x = toff_dict['x']
    folder  = r'D:\measuring\data\QEC_data\figs\Encoding'

    
    for state in ['Z','mZ','Y','mY', 'X','mX']:
        fig,ax = plt.subplots()
        y = toff_dict['toff_'+state+'y']
        y_err = toff_dict['toff_'+state+'y_err']
        ax.errorbar(x,y,yerr = y_err, label = 'toffoli state '+state)


        ax.set_title('toffoli_decoded_state_'+state+'.png')                

        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')
    
        try:
            fig.savefig(
                os.path.join(folder,'toffoli_decoded_state_'+state+'.png'))
        except:
            print 'Figure has not been saved.'


def no_QEC_plot_process_fids():

    process_dict = no_QEC_process_fids()
    process_dict_idle = no_QEC_process_fids(idle = True)

    single_process_dict = single_Qubit_no_QEC_process_fids()

    toff_process_dict = no_QEC_toffoli_fids()
    toff_dict_idle = no_QEC_toffoli_fids(idle = True)

    x = process_dict['x']
    folder  = r'D:\measuring\data\QEC_data\figs\process fidelities'

    t_list = ['1','2','3','avg']
    color_list = ['c','r','b','g']



    fig,ax = plt.subplots()

    for i in range(4):
        y = process_dict['dec_'+t_list[i]+'_y']
        y_err = process_dict['dec_'+t_list[i]+'_y_err']
        y_idle = process_dict_idle['dec_'+t_list[i]+'_y']
        y_idle_err = process_dict_idle['dec_'+t_list[i]+'_y_err']      
        ax.errorbar(x,y,yerr=y_err,color = color_list[i], label =  'decode to '+ t_list[i])
        ax.errorbar(x,y_idle,yerr=y_idle_err,color = color_list[i], ls = ':')#, label =  'with idling, decode to '+ t_list[i])

    y = toff_process_dict['toff_process_y']
    y_err = toff_process_dict['toff_process_y_err']
    y_idle = toff_dict_idle['toff_process_y']
    y_idle_err = toff_dict_idle['toff_process_y_err']  
    ax.errorbar(x,y,yerr=y_err,color = 'k', label =  'toffoli decoded')
    ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'k',  ls = ':', label =  'with idling, toffoli decoded')
    ax.set_ylim(-0.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('no_correction_process_fidelity'+'.png')                
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = 'b')
    ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = 'b')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Process Fidelity')
    plt.legend()
    
    try:
        fig.savefig(
            os.path.join(folder,'no_correction_process_fidelity'+'.png'))
    except:
        print 'Figure has not been saved.'


def single_Qubit_no_QEC_plot_process_fids():

    process_dict = single_Qubit_no_QEC_process_fids()

    x = process_dict['x']
    folder  = r'D:\measuring\data\QEC_data\figs\process fidelities'

    t_list = ['1','2','3','avg']
    color_list = ['c','r','b','g']



    fig,ax = plt.subplots()

    for i in range(4):
        y = process_dict['dec_'+t_list[i]+'_y']
        y_err = process_dict['dec_'+t_list[i]+'_y_err']
      
        ax.errorbar(x,y,yerr=y_err,color = color_list[i], label =  'Qubit to '+ t_list[i])
    ax.set_ylim(-0.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('single_qubit_process_fidelity'+'.png')                
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = 'b')
    ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = 'b')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Process Fidelity')
    plt.legend()
    
    try:
        fig.savefig(
            os.path.join(folder,'single_qubit_process_fidelity'+'.png'))
    except:
        print 'Figure has not been saved.'


''' old scripts and test run '''


def QEC_plot_RO_average(run = 1, no_error = '00'):
    color = ['k','b','m']    
    dataset_dict = {}
    folder  = r'D:\measuring\data\QEC_data\figs\el RO avg'
    for state in ['Z','mZ','Y','mY', 'X','mX']:
        dataset_dict[state] = {}
        if state == 'X' or state == 'mX':
            RO_list = [6]
        if state == 'Y' or state == 'mY':
            RO_list = [4,5,6]
        if state == 'Z' or state == 'mZ':
            RO_list = [0,1,2]

        fig,ax = plt.subplots() 
        print RO_list
        for i, RO in enumerate(RO_list):
            dataset_dict[state]['Tomo_'+str(RO)] = {}
            print RO
            print state
            dataset_dict[state]['Tomo_'+str(RO)] = QEC_sum_data_single_state_RO(run = run, no_error = no_error,state = state,RO = RO)

            x = dataset_dict[state]['Tomo_'+str(RO)]['x']
            y1 = dataset_dict[state]['Tomo_'+str(RO)]['RO_contrast_pos_error']
            y2 = dataset_dict[state]['Tomo_'+str(RO)]['RO_contrast_neg_error']
            y1_err = dataset_dict[state]['Tomo_'+str(RO)]['u_RO_contrast_pos_error']
            y2_err = dataset_dict[state]['Tomo_'+str(RO)]['u_RO_contrast_neg_error']
            ax.set_title('state '+ state + ' el RO average')
            ax.errorbar(x,y1,yerr=y1_err,color = color[i] , label = 'positive error RO '+ str(RO))
            ax.errorbar(x,y2,yerr=y2_err,color = color[i], ls = '--', label = 'negative error RO '+ str(RO))
        ax.legend()

        try:
            fig.savefig(
                os.path.join(folder,'error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'.png'))
        except:
            print 'Figure has not been saved.'

def no_QEC_plot_RO_average():
    color = ['k','b','m']    
    dataset_dict = {}
    folder  = r'D:\measuring\data\QEC_data\figs\el RO avg'
    for state in ['Z','mZ','Y','mY', 'X','mX']:
        dataset_dict[state] = {}
        if state == 'X' or state == 'mX':
            RO_list = [6]
        if state == 'Y' or state == 'mY':
            RO_list = [4,5,6]
        if state == 'Z' or state == 'mZ':
            RO_list = [0,1,2]

        fig,ax = plt.subplots() 
        print RO_list
        for i, RO in enumerate(RO_list):
            dataset_dict[state]['Tomo_'+str(RO)] = {}
            print RO
            print state
            dataset_dict[state]['Tomo_'+str(RO)] = no_QEC_data_single_state_RO(state = state,RO = RO,load_set = True)

            x = dataset_dict[state]['Tomo_'+str(RO)]['x']
            y1 = dataset_dict[state]['Tomo_'+str(RO)]['RO_contrast_pos_error']
            y2 = dataset_dict[state]['Tomo_'+str(RO)]['RO_contrast_neg_error']
            y1_err = dataset_dict[state]['Tomo_'+str(RO)]['u_RO_contrast_pos_error']
            y2_err = dataset_dict[state]['Tomo_'+str(RO)]['u_RO_contrast_neg_error']
            ax.set_title('no QEC state '+ state + ' el RO average')
            ax.errorbar(x,y1,yerr=y1_err,color = color[i] , label = 'positive error RO '+ str(RO))
            ax.errorbar(x,y2,yerr=y2_err,color = color[i], ls = '--', label = 'negative error RO '+ str(RO))
        ax.legend()

        try:
            fig.savefig(
                os.path.join(folder,'error_syn_'+'no_correction'+'_run_'+str(0)+'_state_'+state+'.png'))
        except:
            print 'Figure has not been saved.'


def plot_test_run_QEC(older_than = None,state_RO_list = ['X6','Y4','Y5','Y6','Z0','Z1','Z2','Z6'],ssro_calib_timestamp = None):

    ssro_calib_folder = get_ssro_folder(ssro_calib_timestamp)
    data = {}

    y   = []
    u_y = []

    y_00 = []
    y_01 = []
    y_10 = []
    y_11 = []

    u_y_00 = []
    u_y_01 = []
    u_y_10 = []
    u_y_11 = []

    p_00 = []
    p_01 = []
    p_10 = []
    p_11 = []

    for state_RO in state_RO_list:
        # folder_p = toolbox.latest_data(contains = 'positive_test_RO'+state_RO[1]+'_k1_sign1_'+state_RO[0], older_than = older_than)
        # folder_n = toolbox.latest_data(contains = 'negative_test_RO'+state_RO[1]+'_k1_sign1_'+state_RO[0], older_than = older_than)


        folder_p = toolbox.latest_data(contains = 'positive_test_RO'+state_RO[1], older_than = older_than)
        folder_n = toolbox.latest_data(contains = 'negative_test_RO'+state_RO[1], older_than = older_than)
        

        data_p = load_QEC_data(folder_p, ssro_calib_folder, post_select = True)
        data_n = load_QEC_data(folder_n, ssro_calib_folder, post_select = True)

        # print folder_p
        # print folder_n
        

        y = y + [abs(((data_p['c0'] - data_n['c0'])/2.)[0])]
        u_y = u_y + [abs( (1./2*(data_p['c0_u']**2 + data_n['c0_u']**2)**0.5   )[0])]

        y_00 = y_00 + [(((data_p['c0_00'] - data_n['c0_00'])/2.)[0])]
        u_y_00 = u_y_00 + [( ( 1./2*(data_p['c0_00_u']**2 + data_n['c0_00_u']**2)**0.5   )[0])]
        y_01 = y_01 + [(((data_p['c0_01'] - data_n['c0_01'])/2.)[0])]
        u_y_01 = u_y_01 + [( ( 1./2*(data_p['c0_01_u']**2 + data_n['c0_01_u']**2)**0.5)[0])]
        y_10 = y_10 + [(((data_p['c0_10'] - data_n['c0_10'])/2.)[0])]
        u_y_10 = u_y_10 + [( ( 1./2*(data_p['c0_10_u']**2 + data_n['c0_10_u']**2)**0.5)[0])]
        y_11 = y_11 + [(((data_p['c0_11'] - data_n['c0_11'])/2.)[0])]
        u_y_11 = u_y_11 + [( ( 1./2*(data_p['c0_11_u']**2 + data_n['c0_11_u']**2)**0.5)[0])]

        p_00 = p_00+[((data_p['p00'] + data_n['p00'])/2)[0]]
        p_01 = p_01+[((data_p['p01'] + data_n['p01'])/2)[0]]
        p_10 = p_10+[((data_p['p10'] + data_n['p10'])/2)[0]]
        p_11 = p_11+[((data_p['p11'] + data_n['p11'])/2)[0]]


    x_ticks = state_RO_list
    x = range(len(state_RO_list))
    print x
    fig,ax = plt.subplots()
    ax.set_title(str(folder_p)+'/'+ '\n expectation values')
    rects = ax.bar(x,y,yerr=u_y,align ='center',ecolor = 'k' )
    ax.set_xticks(x)
    ax.set_xticklabels(x_ticks)
    ax.set_xlim(x[0]-0.5,x[-1]+0.5)
    ax.set_ylim(0,1)
    def autolabel(rects):
        for ii,rect in enumerate(rects):
            height = rect.get_height()
            plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2)) +'('+ str(int(round(u_y[ii]*100))) +')',
                ha='center', va='bottom')
    
    autolabel(rects)

    fig,ax = plt.subplots()
    ax.set_title(str(folder_p)+'/'+ '\n post_selected')
    ax.errorbar(x, y_00, yerr=u_y_00, label = '00',color = 'k', fmt='o' )
    ax.errorbar(x, y_01, yerr=u_y_01, label = '01',color = 'c', fmt='o' )
    ax.errorbar(x, y_10, yerr=u_y_10, label = '10',color = 'g', fmt='o' )
    ax.errorbar(x, y_11, yerr=u_y_11, label = '11',color = 'r', fmt='o' )
    ax.set_xlim(x[0]-0.5,x[-1]+0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(x_ticks)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    fig,ax = plt.subplots()
    ax.set_title(str(folder_p)+'/'+ '\n probabilities')
    ax.plot(x,p_00, 'co', label = 'p00')
    ax.plot(x,p_01, 'ko', label = 'p01')
    ax.plot(x,p_10, 'mo', label = 'p10')
    ax.plot(x,p_11, 'bo', label = 'p11')
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_xlabel('error probability')
    ax.set_ylabel('outcome probability')      
    ax.set_xlim(x[0]-0.5,x[-1]+0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(x_ticks)

    print p_00[0][0] + p_01[0][0] + p_10[0][0] +p_11[0][0]
    print p_00[1][0] + p_01[1][0] + p_10[1][0] +p_11[1][0]
    print p_00[2][0] + p_01[2][0] + p_10[2][0] +p_11[2][0]
    # print p_00[3][0] + p_01[3][0] + p_10[3][0] +p_11[3][0]

# def QEC_sum_fidelities(date = None, state  = 'Z', no_error = '11'):      

#     sum_list = {}
#     sum_list['Z_3qb'] = [1,1,1,1,1,1,1]
#     sum_list['mZ_3qb'] = [-1,-1,-1,1,1,1,-1]
#     sum_list['Y_3qb'] = [1,1,1,1,-1,-1,-1]
#     sum_list['mY_3qb'] = [1,1,1,-1,1,1,1]
#     sum_list['X_3qb'] = [1,1,1,-1,-1,-1,1]
#     sum_list['mX_3qb'] = [1,1,1,1,1,1,-1]

#     sum_list['Z_toff'] = [1,1,1,0,0,0,-1]
#     sum_list['mZ_toff'] = [-1,-1,-1,0,0,0,1]
#     sum_list['Y_toff'] = [0,0,0,-1,-1,-1,-1]
#     sum_list['mY_toff'] = [0,0,0,1,1,1,1]
#     sum_list['X_toff'] = [0,0,0,0,0,0,1]
#     sum_list['mX_toff'] = [0,0,0,0,0,0,-1]


#     p_list = ['p00','p01','p10','p11']
#     y_list = ['y','y_00','y_01','y_10','y_11']
#     y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']

#     QEC_state_dict = {}
#     QEC_temp_dict = {}


#     dataset_dict = load_QEC_dataset(date = date, no_error = no_error)
#     QEC_temp_dict  = dataset_dict[state]

#     sum_type = '_3qb'
#     QEC_state_dict[state+sum_type] = {}
#     QEC_state_dict[state+sum_type]['x'] = QEC_temp_dict['Tomo_0']['x']
#     for v in range(5):
#         for RO in range(7):
#             if RO ==0:
#                 QEC_state_dict[state+sum_type][y_list[v]] = 1/8.*(1+sum_list[state+sum_type][RO]*QEC_temp_dict['Tomo_'+str(RO)][y_list[v]])
#                 QEC_state_dict[state+sum_type][y_err_list[v]] = QEC_temp_dict['Tomo_'+str(RO)][y_err_list[v]]**2
#             else:
#                 QEC_state_dict[state+sum_type][y_list[v]] = QEC_state_dict[state+sum_type][y_list[v]]+1/8.*(sum_list[state+sum_type][RO]*QEC_temp_dict['Tomo_'+str(RO)][y_list[v]])
#                 QEC_state_dict[state+sum_type][y_err_list[v]] = QEC_state_dict[state+sum_type][y_err_list[v]]+ QEC_temp_dict['Tomo_'+str(RO)][y_err_list[v]]**2
#             if RO == 6:
#                 QEC_state_dict[state+sum_type][y_err_list[v]] = 1/8.*QEC_state_dict[state+sum_type][y_err_list[v]]**0.5

#     for v in range(4):
#         for RO in range(7):
#             if RO ==0:
#                 QEC_state_dict[state+sum_type][p_list[v]] = 1/7.*(QEC_temp_dict['Tomo_'+str(RO)][p_list[v]])
#             else:
#                 QEC_state_dict[state+sum_type][p_list[v]] = QEC_state_dict[state+sum_type][p_list[v]]+ 1/7.*QEC_temp_dict['Tomo_'+str(RO)][p_list[v]]
    
#     if state in ['Z','mZ','Y','mY']:    
#         sum_type = '_toff'
#         QEC_state_dict[state+sum_type] = {}
#         QEC_state_dict[state+sum_type]['x'] = QEC_temp_dict['Tomo_0']['x']
#         for v in range(5):
#             for RO in range(7):
#                 if RO ==0:
#                     QEC_state_dict[state+sum_type][y_list[v]] = 1/2.*(sum_list[state+sum_type][RO]*QEC_temp_dict['Tomo_'+str(RO)][y_list[v]])
#                     QEC_state_dict[state+sum_type][y_err_list[v]] = abs(sum_list[state+sum_type][RO])*QEC_temp_dict['Tomo_'+str(RO)][y_err_list[v]]**2
#                 else:
#                     QEC_state_dict[state+sum_type][y_list[v]] = QEC_state_dict[state+sum_type][y_list[v]]+1/2.*(sum_list[state+sum_type][RO]*QEC_temp_dict['Tomo_'+str(RO)][y_list[v]])
#                     QEC_state_dict[state+sum_type][y_err_list[v]] = QEC_state_dict[state+sum_type][y_err_list[v]]+ abs(sum_list[state+sum_type][RO])*QEC_temp_dict['Tomo_'+str(RO)][y_err_list[v]]**2
#                 if RO == 6:
#                     QEC_state_dict[state+sum_type][y_err_list[v]] = 1/4.*QEC_state_dict[state+sum_type][y_err_list[v]]**0.5
#                     QEC_state_dict[state+sum_type][y_list[v]] = 1/2.*(1+QEC_state_dict[state+sum_type][y_list[v]]) # make it fidelity

#         for v in range(4):
#             for RO in range(7):
#                 if RO ==0:
#                     QEC_state_dict[state+sum_type][p_list[v]] = 1/4.*(abs(sum_list[state+sum_type][RO])*QEC_temp_dict['Tomo_'+str(RO)][p_list[v]])
#                 else:
#                     QEC_state_dict[state+sum_type][p_list[v]] = QEC_state_dict[state+sum_type][p_list[v]]+ 1/4.*abs(sum_list[state+sum_type][RO])*QEC_temp_dict['Tomo_'+str(RO)][p_list[v]]
    
#     else: # state X and -X need only one exp value ZZZ
#         sum_type = '_toff'
#         QEC_state_dict[state+sum_type] = {}
#         QEC_state_dict[state+sum_type]['x'] = QEC_temp_dict['Tomo_0']['x']
#         for v in range(5):
#             QEC_state_dict[state+sum_type][y_list[v]] = (1+sum_list[state+sum_type][6]*QEC_temp_dict['Tomo_'+str(6)][y_list[v]])/2.
#             QEC_state_dict[state+sum_type][y_err_list[v]] = QEC_temp_dict['Tomo_'+str(6)][y_err_list[v]]/2
#         for v in range(4):
#             QEC_state_dict[state+sum_type][p_list[v]] = QEC_temp_dict['Tomo_'+str(6)][p_list[v]]        

#     return QEC_state_dict
    
# def plot_QEC_sum_fidelities(date = None,state = 'Z',no_error = '00'):

#     QEC_state_dict = QEC_sum_fidelities(date = date, state  = state,no_error = no_error)
#     folder  = r'K:\ns\qt\Diamond\Projects\QEC LT\QEC data'  
#     for sum_type in ['_3qb', '_toff']:

#         x = QEC_state_dict[state+sum_type]['x']
#         y = QEC_state_dict[state+sum_type]['y']
#         y_00 = QEC_state_dict[state+sum_type]['y_00']
#         y_01 = QEC_state_dict[state+sum_type]['y_01']
#         y_10 = QEC_state_dict[state+sum_type]['y_10']
#         y_11 = QEC_state_dict[state+sum_type]['y_11']

#         y_err = QEC_state_dict[state+sum_type]['y_err']
#         y_err_00 = QEC_state_dict[state+sum_type]['y_err_00']
#         y_err_01 = QEC_state_dict[state+sum_type]['y_err_01']
#         y_err_10 = QEC_state_dict[state+sum_type]['y_err_10']
#         y_err_11 = QEC_state_dict[state+sum_type]['y_err_11']

#         p_00 = QEC_state_dict[state+sum_type]['p00']
#         p_01 = QEC_state_dict[state+sum_type]['p01']
#         p_10 = QEC_state_dict[state+sum_type]['p10']
#         p_11 = QEC_state_dict[state+sum_type]['p11']

#         fig,ax = plt.subplots() 
#         ax.errorbar(x,y,yerr=y_err)
#         ax.set_ylim(-0.1,1.1)
#         ax.set_xlim(-0.1,1.1)
#         ax.set_title('errorsyn_'+no_error+'_state_'+state+'QEC_'+sum_type)
#         ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
#         ax.set_xlabel('error probability')
#         ax.set_ylabel('Fidelity'+sum_type)
#         try:
#             fig.savefig(
#                 os.path.join(folder,date+'_errorsyn_'+no_error+'_state_'+state+'QEC_'+sum_type+'.png'))
#         except:
#             print 'Figure has not been saved.'

#         fig,ax = plt.subplots() 
#         ax.errorbar(x,y_00,yerr=y_err_00,color = 'c', label = 'y_00' )
#         ax.errorbar(x,y_01,yerr=y_err_01,color = 'k', label = 'y_01' )
#         ax.errorbar(x,y_10,yerr=y_err_10,color = 'm', label = 'y_10' )
#         ax.errorbar(x,y_11,yerr=y_err_11,color = 'b', label = 'y_11' )
#         ax.set_ylim(-0.1,1.1)
#         ax.set_xlim(-0.1,1.1)
#         plt.legend()
#         ax.set_title(date+'_errorsyn_'+no_error+'_state_'+state+'postselect_'+sum_type)
#         ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
#         ax.set_xlabel('error probability')
#         ax.set_ylabel('Fidelity'+sum_type)

#         try:
#             fig.savefig(
#                 os.path.join(folder,date+'_errorsyn_'+no_error+'_state_'+state+'postselect_'+sum_type+'.png'))
#         except:
#             print 'Figure has not been saved.'




#         fig,ax = plt.subplots()
#         ax.set_title(str(folder)+'/'+ '\n probabilities')
#         ax.plot(x,p_00, 'c', label = 'p00')
#         ax.plot(x,p_01, 'k', label = 'p01')
#         ax.plot(x,p_10, 'm', label = 'p10')
#         ax.plot(x,p_11, 'b', label = 'p11')
#         ax.set_ylim(-0.1,1.1)
#         ax.set_xlim(-0.1,1.1)
#         ax.plot(x,p_00+p_01+p_10+p_11,label = 'sum')
#         plt.legend()
#         ax.set_xlabel('error probability')
#         ax.set_ylabel('outcome probability'+sum_type)   
#         ax.set_title(date+'_errorsyn_'+no_error+'_state_'+state+'postselect_'+sum_type)               
        
#         try:
#             fig.savefig(
#                 os.path.join(folder,date+'_errorsyn_'+no_error+'_state_'+state+'probabilities_'+sum_type+'.png'))
#         except:
#             print 'Figure has not been saved.'

# def QEC_process_fids(date = None,no_error = '00'):

#     dataset_dict = load_QEC_dataset(date = date, no_error = no_error)
#     process_dict = {}

#     y_list = ['y','y_00','y_01','y_10','y_11']
#     y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']

#     for v in range(5):
#         print v
#         process_dict['dec_1_'+y_list[v]] = {}
#         process_dict['dec_2_'+y_list[v]] = {}
#         process_dict['dec_3_'+y_list[v]] = {}
#         # process_dict['dec_toff_'+y_list[v]] = {}
#         process_dict['dec_avg_'+y_list[v]] = {}
#         process_dict['dec_1_'+y_err_list[v]] = {}
#         process_dict['dec_2_'+y_err_list[v]] = {}
#         process_dict['dec_3_'+y_err_list[v]] = {}
#         # process_dict['toff_'+y_err_list[v]] = {}
#         process_dict['dec_avg_'+y_err_list[v]] = {}


#         process_dict['dec_1_'+y_list[v]] = 1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(0)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(0)][y_list[v]]
#                         - dataset_dict['Y']['Tomo_'+str(5)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(5)][y_list[v]]
#                         + dataset_dict['X']['Tomo_'+str(6)][y_list[v]] - dataset_dict['mX']['Tomo_'+str(6)][y_list[v]])

#         process_dict['dec_1_'+y_err_list[v]] = 1/8.*(dataset_dict['Z']['Tomo_'+str(0)][y_err_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(0)][y_err_list[v]]**2 
#                         + dataset_dict['Y']['Tomo_'+str(5)][y_err_list[v]]**2 +  dataset_dict['mY']['Tomo_'+str(5)][y_err_list[v]]**2 
#                         + dataset_dict['X']['Tomo_'+str(6)][y_err_list[v]]**2 + dataset_dict['mX']['Tomo_'+str(6)][y_err_list[v]]**2 )**0.5

#         process_dict['dec_2_'+y_list[v]] = 1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(1)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(1)][y_list[v]]
#                 - dataset_dict['Y']['Tomo_'+str(6)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(6)][y_list[v]]
#                 + dataset_dict['X']['Tomo_'+str(6)][y_list[v]] - dataset_dict['mX']['Tomo_'+str(6)][y_list[v]])

#         process_dict['dec_2_'+y_err_list[v]] = 1/8.*(dataset_dict['Z']['Tomo_'+str(1)][y_err_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(1)][y_err_list[v]]**2 
#                         + dataset_dict['Y']['Tomo_'+str(6)][y_err_list[v]]**2 +  dataset_dict['mY']['Tomo_'+str(6)][y_err_list[v]]**2 
#                         + dataset_dict['X']['Tomo_'+str(6)][y_err_list[v]]**2 + dataset_dict['mX']['Tomo_'+str(6)][y_err_list[v]]**2 )**0.5

#         process_dict['dec_3_'+y_list[v]] = 1/4. + 1/8.*(dataset_dict['Z']['Tomo_'+str(2)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(2)][y_list[v]]
#                         - dataset_dict['Y']['Tomo_'+str(4)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(4)][y_list[v]]
#                         + dataset_dict['X']['Tomo_'+str(6)][y_list[v]] - dataset_dict['mX']['Tomo_'+str(6)][y_list[v]])

#         process_dict['dec_3_'+y_err_list[v]] = 1/8.*(dataset_dict['Z']['Tomo_'+str(3)][y_err_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(3)][y_err_list[v]]**2 
#                         + dataset_dict['Y']['Tomo_'+str(4)][y_err_list[v]]**2 +  dataset_dict['mY']['Tomo_'+str(4)][y_err_list[v]]**2 
#                         + dataset_dict['X']['Tomo_'+str(6)][y_err_list[v]]**2 + dataset_dict['mX']['Tomo_'+str(6)][y_err_list[v]]**2 )**0.5


#         process_dict['dec_avg_'+y_list[v]] = 1/3.*(process_dict['dec_1_'+y_list[v]]+process_dict['dec_2_'+y_list[v]]+process_dict['dec_3_'+y_list[v]])
#         process_dict['dec_avg_'+y_err_list[v]] = 1/3.*(process_dict['dec_1_'+y_err_list[v]]**2+process_dict['dec_2_'+y_err_list[v]]**2+process_dict['dec_3_'+y_err_list[v]]**2)**0.5


#     process_dict['x'] = dataset_dict['Z']['Tomo_'+str(0)]['x']
#     return process_dict

# def QEC_plot_process_fids(date = None,no_error = '00'):

#     process_dict = QEC_process_fids(date = date, no_error = no_error)

#     x = process_dict['x']
#     folder  = r'K:\ns\qt\Diamond\Projects\QEC LT\QEC data'  

#     t_list = ['1','2','3','avg']
#     color_list = ['c','r','b','g']

#     fig,ax = plt.subplots() 
#     for i in range(4):
#         y = process_dict['dec_'+t_list[i]+'_y']
#         y_err = process_dict['dec_'+t_list[i]+'_y_err']
#         ax.errorbar(x,y,yerr=y_err,color = color_list[i], label =  'decode to '+ t_list[i])
#     ax.set_ylim(-0.1,1.1)
#     ax.set_xlim(-0.1,1.1)
#     ax.set_title(date+'_errorsyn_'+no_error+'_process_fids.png')                
#     ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
#     ax.set_xlabel('error probability')
#     ax.set_ylabel('Process Fidelity')
#     plt.legend()
    
#     try:
#         fig.savefig(
#             os.path.join(folder,date+'_errorsyn_'+no_error+'process_fids'+'.png'))
#     except:
#         print 'Figure has not been saved.'
