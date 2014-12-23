import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from analysis.lib.QEC import ConditionalParity as CP
reload (CP)
import h5py
import csv

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


''' New functions THT '''

def load_QEC_data(folder, ssro_calib_folder, post_select = True):
    ''' Loads a QEC measurment and returns all 
    the results in a dictionairy'''

    a = CP.ConditionalParityAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata', post_select_QEC = False)
    print ssro_calib_folder
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
                print folder
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


''' these functions are used to open/close save/load from and to HDF5 files '''

def openfile_single_state_RO(sym = '00', RO = 1, state = 'Z'):
    name = 'error_sym_' + sym + '_state_'+state + '_RO_'+str(RO)+ '.hdf5'
    datafile = h5py.File(os.path.join(r'D:\measuring\data\QEC_data\QEC_data', name)) 
    return datafile


def openfile(name = ''):
   datafile = h5py.File(os.path.join(r'K:\ns\qt\Diamond\Projects\QEC LT\QEC data10', name)) 
   return datafile

def closefile(datafile):
    datafile.close()

def save_data_hdf5file(datafile, data_dict, state, RO):
    f = datafile
    f_grp = f.create_group('state_'+state+'_RO_'+str(RO))

    for item in data_dict:
        f.attrs [item] = data_dict[item]
        f_grp.create_dataset (item, data = data_dict[item])

def save_data_single_state_RO_hdf5file(datafile, data_dict):
    f = datafile
    f_grp = f.create_group('data')

    for item in data_dict:
        f.attrs [item] = data_dict[item]
        f_grp.create_dataset (item, data = data_dict[item])

def load_data_hdf5file(datafile,state, RO):
    f = datafile
    f_grp = f['/'+'state_'+state+'_RO_'+str(RO)]

    data_dict = {}
    for item in f_grp.keys():
        data_dict[item] = f_grp[item].value

    return data_dict

def load_single_data_hdf5file(datafile):
    f = datafile
    f_grp = f['/'+'data']

    data_dict = {}
    for item in f_grp.keys():
        data_dict[item] = f_grp[item].value
    return data_dict


''' here you save new data '''

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
    for p in range(4):
            QEC_data_dict[p_list[p]] = {}
            
            QEC_data_dict[p_list[p]] = (QEC_dict[str(-1)]['positive'][p_list[p]]+
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

def no_QEC_data_single_state_RO(older_than = None,state = 'Z',RO = 0):

    QEC_data_dict = {}
    u_list = ['u_c0']
    c_list = ['c0']
    y_list = ['y']
    y_err_list = ['y_err']


    QEC_dict, folder = no_QEC_create_data_dict(older_than = older_than, RO = RO, state = state)
    for v in range(1):
        QEC_data_dict[y_list[v]] = {}
        QEC_data_dict[y_err_list[v]] = {}


        QEC_data_dict[y_list[v]] = (QEC_dict[str(-1)]['positive'][c_list[v]]+
                                                        QEC_dict[str(1)]['positive'][c_list[v]]-
                                                        QEC_dict[str(-1)]['negative'][c_list[v]]-
                                                        QEC_dict[str(1)]['negative'][c_list[v]])/4

        
        
        QEC_data_dict[y_err_list[v]] = (QEC_dict[str(-1)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(1)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(-1)]['negative'][u_list[v]]**2+
                                                        QEC_dict[str(1)]['negative'][u_list[v]]**2)**0.5/4


    QEC_data_dict['x'] = QEC_dict[str(1)]['positive']['x']
    QEC_data_dict['folder'] = folder

    return QEC_data_dict, folder

def save_QEC_dataset_single(older_than = None,sym = '11', RO = 1, state = 'Z'):
    
    datafile = openfile_single_state_RO(sym = sym, RO = RO, state = state)

    QEC_temp_dict, folder = QEC_data_single_state_RO(older_than = older_than, RO = RO, state = state, sym = sym)

    save_data_single_state_RO_hdf5file(datafile, QEC_temp_dict)
       
    closefile(datafile)

def save_QEC_dataset(older_than = None, no_error = '11_1'):
    
    datafile = openfile(name = 'QEC_'+older_than[0:8]+'_error_syn_'+no_error+'.hdf5')
    QEC_temp_dict = {}

    for state in ['mY','Y','Z','mZ','X','mX']:
        QEC_temp_dict[state] = {}
        for RO  in range(7):
            QEC_temp_dict[state]['Tomo_'+str(RO)] = {}
            QEC_temp_dict[state]['Tomo_'+str(RO)], folder = QEC_data_single_state_RO(older_than = older_than, RO = RO, state = state)
            save_data_hdf5file(datafile, state, RO)
       
    closefile(datafile)


def load_QEC_dataset(date = None, no_error = '00'):

    datafile = openfile(name = 'QEC_'+date+'_error_syn_'+no_error+'.hdf5')
    QEC_temp_dict = {}
    for state in ['mY','Y','Z','mZ','X','mX']:
        QEC_temp_dict[state] = {}
        for RO  in range(7):
            QEC_temp_dict[state]['Tomo_'+str(RO)] = {}
            QEC_temp_dict[state]['Tomo_'+str(RO)] = load_data_hdf5file(datafile,state, RO)

    closefile(datafile)

    return QEC_temp_dict

def load_QEC_dataset_single(sym = '11', RO = 1, state = 'Z'):
    
    datafile = openfile_single_state_RO(sym = sym, RO = RO, state = state)

    QEC_temp_dict = load_single_data_hdf5file(datafile)
    
    closefile(datafile)
    return QEC_temp_dict

''' from here you can plot data taken from an existing HDF5 file '''

def QEC_plot_single_state_RO(date = '20141120', no_error = '00',state = 'Z',RO = 0, load_set = True, older_than = None, e_sign = None,plot_guide = True):        
    
    
    if load_set == True:
        QEC_data_dict = {}   
        dataset_dict = load_QEC_dataset_single(sym = no_error, state = state, RO = RO)
        QEC_data_dict  = dataset_dict
        
    else:
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

def no_QEC_plot_single_state_RO(date = '20141120',state = 'Z',RO = 0, load_set = False, older_than = None):        
    
    
    if load_set == True:
        QEC_data_dict = {}   
        dataset_dict = load_QEC_dataset(date = date, no_error = no_error)
        QEC_data_dict  = dataset_dict[state]['Tomo_'+str(RO)]
        
    else:
        QEC_data_dict, folder =  no_QEC_data_single_state_RO(older_than = older_than,state = state,RO = RO)
    
    
    folder  = r'K:\ns\qt\Diamond\Projects\QEC LT\QEC data'   


    x = QEC_data_dict['x']
    y = QEC_data_dict['y']


    y_err = QEC_data_dict['y_err']

    fig,ax = plt.subplots() 
    ax.errorbar(x,y,yerr=y_err,color = 'k' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('state_'+state+'_RO_'+str(RO)+'_noQEC')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')
    try:
        fig.savefig(
            os.path.join(folder,date+'no_QEC'+'_state_'+state+'_RO_'+str(RO)+'.png'))
    except:
        print 'Figure has not been saved.'


def QEC_sum_fidelities(date = None, state  = 'Z', no_error = '11'):      

    sum_list = {}
    sum_list['Z_3qb'] = [1,1,1,1,1,1,1]
    sum_list['mZ_3qb'] = [-1,-1,-1,1,1,1,-1]
    sum_list['Y_3qb'] = [1,1,1,1,-1,-1,-1]
    sum_list['mY_3qb'] = [1,1,1,-1,1,1,1]
    sum_list['X_3qb'] = [1,1,1,-1,-1,-1,1]
    sum_list['mX_3qb'] = [1,1,1,1,1,1,-1]

    sum_list['Z_toff'] = [1,1,1,0,0,0,-1]
    sum_list['mZ_toff'] = [-1,-1,-1,0,0,0,1]
    sum_list['Y_toff'] = [0,0,0,-1,-1,-1,-1]
    sum_list['mY_toff'] = [0,0,0,1,1,1,1]
    sum_list['X_toff'] = [0,0,0,0,0,0,1]
    sum_list['mX_toff'] = [0,0,0,0,0,0,-1]


    p_list = ['p00','p01','p10','p11']
    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']

    QEC_state_dict = {}
    QEC_temp_dict = {}


    dataset_dict = load_QEC_dataset(date = date, no_error = no_error)
    QEC_temp_dict  = dataset_dict[state]

    sum_type = '_3qb'
    QEC_state_dict[state+sum_type] = {}
    QEC_state_dict[state+sum_type]['x'] = QEC_temp_dict['Tomo_0']['x']
    for v in range(5):
        for RO in range(7):
            if RO ==0:
                QEC_state_dict[state+sum_type][y_list[v]] = 1/8.*(1+sum_list[state+sum_type][RO]*QEC_temp_dict['Tomo_'+str(RO)][y_list[v]])
                QEC_state_dict[state+sum_type][y_err_list[v]] = QEC_temp_dict['Tomo_'+str(RO)][y_err_list[v]]**2
            else:
                QEC_state_dict[state+sum_type][y_list[v]] = QEC_state_dict[state+sum_type][y_list[v]]+1/8.*(sum_list[state+sum_type][RO]*QEC_temp_dict['Tomo_'+str(RO)][y_list[v]])
                QEC_state_dict[state+sum_type][y_err_list[v]] = QEC_state_dict[state+sum_type][y_err_list[v]]+ QEC_temp_dict['Tomo_'+str(RO)][y_err_list[v]]**2
            if RO == 6:
                QEC_state_dict[state+sum_type][y_err_list[v]] = 1/8.*QEC_state_dict[state+sum_type][y_err_list[v]]**0.5

    for v in range(4):
        for RO in range(7):
            if RO ==0:
                QEC_state_dict[state+sum_type][p_list[v]] = 1/7.*(QEC_temp_dict['Tomo_'+str(RO)][p_list[v]])
            else:
                QEC_state_dict[state+sum_type][p_list[v]] = QEC_state_dict[state+sum_type][p_list[v]]+ 1/7.*QEC_temp_dict['Tomo_'+str(RO)][p_list[v]]
    
    if state in ['Z','mZ','Y','mY']:    
        sum_type = '_toff'
        QEC_state_dict[state+sum_type] = {}
        QEC_state_dict[state+sum_type]['x'] = QEC_temp_dict['Tomo_0']['x']
        for v in range(5):
            for RO in range(7):
                if RO ==0:
                    QEC_state_dict[state+sum_type][y_list[v]] = 1/2.*(sum_list[state+sum_type][RO]*QEC_temp_dict['Tomo_'+str(RO)][y_list[v]])
                    QEC_state_dict[state+sum_type][y_err_list[v]] = abs(sum_list[state+sum_type][RO])*QEC_temp_dict['Tomo_'+str(RO)][y_err_list[v]]**2
                else:
                    QEC_state_dict[state+sum_type][y_list[v]] = QEC_state_dict[state+sum_type][y_list[v]]+1/2.*(sum_list[state+sum_type][RO]*QEC_temp_dict['Tomo_'+str(RO)][y_list[v]])
                    QEC_state_dict[state+sum_type][y_err_list[v]] = QEC_state_dict[state+sum_type][y_err_list[v]]+ abs(sum_list[state+sum_type][RO])*QEC_temp_dict['Tomo_'+str(RO)][y_err_list[v]]**2
                if RO == 6:
                    QEC_state_dict[state+sum_type][y_err_list[v]] = 1/4.*QEC_state_dict[state+sum_type][y_err_list[v]]**0.5
                    QEC_state_dict[state+sum_type][y_list[v]] = 1/2.*(1+QEC_state_dict[state+sum_type][y_list[v]]) # make it fidelity

        for v in range(4):
            for RO in range(7):
                if RO ==0:
                    QEC_state_dict[state+sum_type][p_list[v]] = 1/4.*(abs(sum_list[state+sum_type][RO])*QEC_temp_dict['Tomo_'+str(RO)][p_list[v]])
                else:
                    QEC_state_dict[state+sum_type][p_list[v]] = QEC_state_dict[state+sum_type][p_list[v]]+ 1/4.*abs(sum_list[state+sum_type][RO])*QEC_temp_dict['Tomo_'+str(RO)][p_list[v]]
    
    else: # state X and -X need only one exp value ZZZ
        sum_type = '_toff'
        QEC_state_dict[state+sum_type] = {}
        QEC_state_dict[state+sum_type]['x'] = QEC_temp_dict['Tomo_0']['x']
        for v in range(5):
            QEC_state_dict[state+sum_type][y_list[v]] = (1+sum_list[state+sum_type][6]*QEC_temp_dict['Tomo_'+str(6)][y_list[v]])/2.
            QEC_state_dict[state+sum_type][y_err_list[v]] = QEC_temp_dict['Tomo_'+str(6)][y_err_list[v]]/2
        for v in range(4):
            QEC_state_dict[state+sum_type][p_list[v]] = QEC_temp_dict['Tomo_'+str(6)][p_list[v]]        

    return QEC_state_dict
    

def plot_QEC_sum_fidelities(date = None,state = 'Z',no_error = '00'):

    QEC_state_dict = QEC_sum_fidelities(date = date, state  = state,no_error = no_error)
    folder  = r'K:\ns\qt\Diamond\Projects\QEC LT\QEC data'  
    for sum_type in ['_3qb', '_toff']:

        x = QEC_state_dict[state+sum_type]['x']
        y = QEC_state_dict[state+sum_type]['y']
        y_00 = QEC_state_dict[state+sum_type]['y_00']
        y_01 = QEC_state_dict[state+sum_type]['y_01']
        y_10 = QEC_state_dict[state+sum_type]['y_10']
        y_11 = QEC_state_dict[state+sum_type]['y_11']

        y_err = QEC_state_dict[state+sum_type]['y_err']
        y_err_00 = QEC_state_dict[state+sum_type]['y_err_00']
        y_err_01 = QEC_state_dict[state+sum_type]['y_err_01']
        y_err_10 = QEC_state_dict[state+sum_type]['y_err_10']
        y_err_11 = QEC_state_dict[state+sum_type]['y_err_11']

        p_00 = QEC_state_dict[state+sum_type]['p00']
        p_01 = QEC_state_dict[state+sum_type]['p01']
        p_10 = QEC_state_dict[state+sum_type]['p10']
        p_11 = QEC_state_dict[state+sum_type]['p11']

        fig,ax = plt.subplots() 
        ax.errorbar(x,y,yerr=y_err)
        ax.set_ylim(-0.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title('errorsyn_'+no_error+'_state_'+state+'QEC_'+sum_type)
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Fidelity'+sum_type)
        try:
            fig.savefig(
                os.path.join(folder,date+'_errorsyn_'+no_error+'_state_'+state+'QEC_'+sum_type+'.png'))
        except:
            print 'Figure has not been saved.'

        fig,ax = plt.subplots() 
        ax.errorbar(x,y_00,yerr=y_err_00,color = 'c', label = 'y_00' )
        ax.errorbar(x,y_01,yerr=y_err_01,color = 'k', label = 'y_01' )
        ax.errorbar(x,y_10,yerr=y_err_10,color = 'm', label = 'y_10' )
        ax.errorbar(x,y_11,yerr=y_err_11,color = 'b', label = 'y_11' )
        ax.set_ylim(-0.1,1.1)
        ax.set_xlim(-0.1,1.1)
        plt.legend()
        ax.set_title(date+'_errorsyn_'+no_error+'_state_'+state+'postselect_'+sum_type)
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Fidelity'+sum_type)

        try:
            fig.savefig(
                os.path.join(folder,date+'_errorsyn_'+no_error+'_state_'+state+'postselect_'+sum_type+'.png'))
        except:
            print 'Figure has not been saved.'




        fig,ax = plt.subplots()
        ax.set_title(str(folder)+'/'+ '\n probabilities')
        ax.plot(x,p_00, 'c', label = 'p00')
        ax.plot(x,p_01, 'k', label = 'p01')
        ax.plot(x,p_10, 'm', label = 'p10')
        ax.plot(x,p_11, 'b', label = 'p11')
        ax.set_ylim(-0.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.plot(x,p_00+p_01+p_10+p_11,label = 'sum')
        plt.legend()
        ax.set_xlabel('error probability')
        ax.set_ylabel('outcome probability'+sum_type)   
        ax.set_title(date+'_errorsyn_'+no_error+'_state_'+state+'postselect_'+sum_type)               
        
        try:
            fig.savefig(
                os.path.join(folder,date+'_errorsyn_'+no_error+'_state_'+state+'probabilities_'+sum_type+'.png'))
        except:
            print 'Figure has not been saved.'


def QEC_process_fids(date = None,no_error = '00'):

    dataset_dict = load_QEC_dataset(date = date, no_error = no_error)
    process_dict = {}

    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']

    for v in range(5):
        print v
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

        process_dict['dec_3_'+y_err_list[v]] = 1/8.*(dataset_dict['Z']['Tomo_'+str(3)][y_err_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(3)][y_err_list[v]]**2 
                        + dataset_dict['Y']['Tomo_'+str(4)][y_err_list[v]]**2 +  dataset_dict['mY']['Tomo_'+str(4)][y_err_list[v]]**2 
                        + dataset_dict['X']['Tomo_'+str(6)][y_err_list[v]]**2 + dataset_dict['mX']['Tomo_'+str(6)][y_err_list[v]]**2 )**0.5


        process_dict['dec_avg_'+y_list[v]] = 1/3.*(process_dict['dec_1_'+y_list[v]]+process_dict['dec_2_'+y_list[v]]+process_dict['dec_3_'+y_list[v]])
        process_dict['dec_avg_'+y_err_list[v]] = 1/3.*(process_dict['dec_1_'+y_err_list[v]]**2+process_dict['dec_2_'+y_err_list[v]]**2+process_dict['dec_3_'+y_err_list[v]]**2)**0.5


    process_dict['x'] = dataset_dict['Z']['Tomo_'+str(0)]['x']
    return process_dict

def QEC_plot_process_fids(date = None,no_error = '00'):

    process_dict = QEC_process_fids(date = date, no_error = no_error)

    x = process_dict['x']
    folder  = r'K:\ns\qt\Diamond\Projects\QEC LT\QEC data'  

    t_list = ['1','2','3','avg']
    color_list = ['c','r','b','g']

    fig,ax = plt.subplots() 
    for i in range(4):
        y = process_dict['dec_'+t_list[i]+'_y']
        y_err = process_dict['dec_'+t_list[i]+'_y_err']
        ax.errorbar(x,y,yerr=y_err,color = color_list[i], label =  'decode to '+ t_list[i])
    ax.set_ylim(-0.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title(date+'_errorsyn_'+no_error+'_process_fids.png')                
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Process Fidelity')
    plt.legend()
    
    try:
        fig.savefig(
            os.path.join(folder,date+'_errorsyn_'+no_error+'process_fids'+'.png'))
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


