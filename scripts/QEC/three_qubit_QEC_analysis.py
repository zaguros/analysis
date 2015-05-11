import numpy as np
import os
import pickle
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from analysis.lib.QEC import ConditionalParity as CP
from analysis.lib.fitting import fit, common, ramsey;reload(common); reload(fit)
import matplotlib.cm as cm
import matplotlib as mpl; reload(mpl)
from pylab import *

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

''' Basic functions '''

def load_QEC_data(folder, ssro_calib_folder, post_select = True, nr_of_parity_msmnts = 2):
    ''' Loads a QEC type measurment and returns all
    the results in a dictionairy'''

    a = CP.ConditionalParityAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata', post_select_QEC = False)
    print ssro_calib_folder
    a.get_electron_ROC(ssro_calib_folder)

    x = a.sweep_pts.reshape(-1)
    c0, c0_u = a.convert_fidelity_to_contrast(a.p0,a.u_p0)

    data_dict = {}
    # data_dict['a'] = a
    data_dict['x']          = x
    data_dict['c0']         = c0
    data_dict['c0_u']       = c0_u

    if post_select:
        if nr_of_parity_msmnts==2:
            a = CP.ConditionalParityAnalysis(folder)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata', post_select_QEC = True)
            a.get_electron_ROC(ssro_calib_folder  , post_select_QEC = True)

            c0_00,c0_00_u =  a.convert_fidelity_to_contrast(a.p0_00,a.u_p0_00)
            c0_01,c0_01_u =  a.convert_fidelity_to_contrast(a.p0_01,a.u_p0_01)
            c0_10,c0_10_u =  a.convert_fidelity_to_contrast(a.p0_10,a.u_p0_10)
            c0_11,c0_11_u =  a.convert_fidelity_to_contrast(a.p0_11,a.u_p0_11)

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

        elif nr_of_parity_msmnts==4:

            a = CP.ConditionalParityAnalysis(folder)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata', post_select_multiple_rounds = True)
            a.get_electron_ROC(ssro_calib_folder  , post_select_multiple_rounds = True)

            c0_0000,c0_0000_u =  a.convert_fidelity_to_contrast(a.p0_0000,a.u_p0_0000)
            c0_0100,c0_0100_u =  a.convert_fidelity_to_contrast(a.p0_0100,a.u_p0_0100)
            c0_1000,c0_1000_u =  a.convert_fidelity_to_contrast(a.p0_1000,a.u_p0_1000)
            c0_1100,c0_1100_u =  a.convert_fidelity_to_contrast(a.p0_1100,a.u_p0_1100)

            c0_0010,c0_0010_u =  a.convert_fidelity_to_contrast(a.p0_0010,a.u_p0_0010)
            c0_0110,c0_0110_u =  a.convert_fidelity_to_contrast(a.p0_0110,a.u_p0_0110)
            c0_1010,c0_1010_u =  a.convert_fidelity_to_contrast(a.p0_1010,a.u_p0_1010)
            c0_1110,c0_1110_u =  a.convert_fidelity_to_contrast(a.p0_1110,a.u_p0_1110)

            c0_0001,c0_0001_u =  a.convert_fidelity_to_contrast(a.p0_0001,a.u_p0_0001)
            c0_0101,c0_0101_u =  a.convert_fidelity_to_contrast(a.p0_0101,a.u_p0_0101)
            c0_1001,c0_1001_u =  a.convert_fidelity_to_contrast(a.p0_1001,a.u_p0_1001)
            c0_1101,c0_1101_u =  a.convert_fidelity_to_contrast(a.p0_1101,a.u_p0_1101)

            c0_0011,c0_0011_u =  a.convert_fidelity_to_contrast(a.p0_0011,a.u_p0_0011)
            c0_0111,c0_0111_u =  a.convert_fidelity_to_contrast(a.p0_0111,a.u_p0_0111)
            c0_1011,c0_1011_u =  a.convert_fidelity_to_contrast(a.p0_1011,a.u_p0_1011)
            c0_1111,c0_1111_u =  a.convert_fidelity_to_contrast(a.p0_1111,a.u_p0_1111)

            data_dict['c0_0000']      = c0_0000; data_dict['c0_0000_u']    = c0_0000_u
            data_dict['c0_0100']      = c0_0100; data_dict['c0_0100_u']    = c0_0100_u
            data_dict['c0_1000']      = c0_1000; data_dict['c0_1000_u']    = c0_1000_u
            data_dict['c0_1100']      = c0_1100; data_dict['c0_1100_u']    = c0_1100_u

            data_dict['c0_0010']      = c0_0010; data_dict['c0_0010_u']    = c0_0010_u
            data_dict['c0_0110']      = c0_0110; data_dict['c0_0110_u']    = c0_0110_u
            data_dict['c0_1010']      = c0_1010; data_dict['c0_1010_u']    = c0_1010_u
            data_dict['c0_1110']      = c0_1110; data_dict['c0_1110_u']    = c0_1110_u

            data_dict['c0_0001']      = c0_0001; data_dict['c0_0001_u']    = c0_0001_u
            data_dict['c0_0101']      = c0_0101; data_dict['c0_0101_u']    = c0_0101_u
            data_dict['c0_1001']      = c0_1001; data_dict['c0_1001_u']    = c0_1001_u
            data_dict['c0_1101']      = c0_1101; data_dict['c0_1101_u']    = c0_1101_u

            data_dict['c0_0011']      = c0_0011; data_dict['c0_0011_u']    = c0_0011_u
            data_dict['c0_0111']      = c0_0111; data_dict['c0_0111_u']    = c0_0111_u
            data_dict['c0_1011']      = c0_1011; data_dict['c0_1011_u']    = c0_1011_u
            data_dict['c0_1111']      = c0_1111; data_dict['c0_1111_u']    = c0_1111_u

            data_dict['p0000'] = a.p0000; data_dict['p0100'] = a.p0100; data_dict['p1000'] = a.p1000; data_dict['p1100'] = a.p1100
            data_dict['p0010'] = a.p0010; data_dict['p0110'] = a.p0110; data_dict['p1010'] = a.p1010; data_dict['p1110'] = a.p1110
            data_dict['p0001'] = a.p0001; data_dict['p0101'] = a.p0101; data_dict['p1001'] = a.p1001; data_dict['p1101'] = a.p1101
            data_dict['p0011'] = a.p0011; data_dict['p0111'] = a.p0111; data_dict['p1011'] = a.p1011; data_dict['p1111'] = a.p1111

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

''' These functions both load and plot individual measurements'''

def plot_single_QEC_result(timestamps = [None], folder_name ='QEC', ssro_calib_timestamp = None,
        post_select = False, save = True, title = None, fontsize = 10) :
    '''
    Plots the results of a single QEC/Encoding measurement, from a raw data folder.
    post_select must be false if no parity measurements
    Length of the timestamp gives eaither a single measurement or a positive/negative one
    '''
    ### Timestamps and folders
    timestamp, folder = get_folder(timestamps[0], folder_name)
    # print folder
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
                elif len_k ==2:
                    QEC_dict[str(error_sign)][direction][item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item]), axis=0)

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


def QEC_create_data_dict_single_error_single_elRO(older_than = None, RO = 0, do_p = False,run = 1, state = 'Z', len_k = 6, sym = '11',error_sign = 1, el_RO = 'positive',sweep_time = False):
    QEC_dict = {}
    k_dict = {}

    if do_p != True:
        for k in range(len_k):
            # print 'k_'+str(k)
            if sweep_time == True:
                print sym +'_sweep_time_'+el_RO+'_RO'+str(RO)+'_k'+str(k)+'_'+state
                print older_than
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
    else:
        if run > 4:
            len_p = 3
        else:
            len_p =2
        for p in range(len_p):
            print older_than
            # print sym +'_'+el_RO+'_RO'+str(RO)+'_p_'+str(p)+'_sign'+ str(error_sign)+'_'+state =='11_positive_RO0_p_0_sign-1_Z'
            timestamp, folder = toolbox.latest_data(contains = sym +'_'+el_RO+'_RO'+str(RO)+'_p_'+str(p)+'_sign'+ str(error_sign)
                                                    +'_'+state, older_than = older_than,return_timestamp = True)
            SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
            print SSRO_folder
            k_dict['p_'+str(p)] ={}
            k_dict['p_'+str(p)] = load_QEC_data(folder, SSRO_folder, post_select = True)

        for item in k_dict['p_0']:
            if len_p == 2:
                QEC_dict[item] = np.concatenate((k_dict['p_0'][item],k_dict['p_1'][item]), axis=0)
            else:
                QEC_dict[item] = np.concatenate((k_dict['p_0'][item],k_dict['p_1'][item],k_dict['p_2'][item]), axis=0)
    return QEC_dict,folder


def no_QEC_create_data_dict_single_error_single_elRO(idle = False,do_p = False,sweep_time = False, older_than = None, RO = 0, state = 'Z',run = 1, error_sign = 1, el_RO = 'positive'):
    QEC_dict = {}
    k_dict = {}
    k_length = 3

    if sweep_time == True:
        k_length =4

    if do_p == True:
        print 'ok'
        print older_than
        print 'no_correct_idle' +'_'+el_RO+'_RO'+str(RO)+'_p_'+str(0)+'_sign'+ str(error_sign)+'_'+state
        timestamp, folder = toolbox.latest_data(contains = 'no_correct_idle' +'_'+el_RO+'_RO'+str(RO)+'_p_'+str(0)+'_sign'+ str(error_sign)
                                                                +'_'+state, older_than = older_than,return_timestamp = True)
        SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
        QEC_dict = load_QEC_data(folder, SSRO_folder, post_select = False)

    elif do_p != True:
        for k in range(k_length):
            # print 'k_'+str(k)
            if run == 3:
                k = k+4
            if idle == False and sweep_time == False:
                # print 'ok'
                print 'no_correct' +'_'+el_RO+'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_sign)+'_'+state
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
                # print folder
            SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
            print SSRO_folder
            k_dict['k_'+str(k)] ={}
            k_dict['k_'+str(k)] = load_QEC_data(folder, SSRO_folder, post_select = False)

        for item in k_dict['k_'+str(k)]:
            print item
            if k_length ==3 and run !=3:
                QEC_dict[item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item]), axis=0)
            if k_length ==3 and run ==3:
                QEC_dict[item] = np.concatenate((k_dict['k_4'][item],k_dict['k_5'][item],k_dict['k_6'][item]), axis=0)            
            if k_length ==6:
                QEC_dict[item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item],k_dict['k_3'][item],k_dict['k_4'][item],k_dict['k_5'][item]), axis=0)
            if k_length == 4:
                QEC_dict[item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item],k_dict['k_3'][item]), axis=0)


    return QEC_dict,folder


def single_Qubit_QEC_create_data_dict_single_error_single_elRO(older_than = None, sweep_time = False,Qubit = 1, state = 'Z', error_sign = 1, el_RO = 'positive',run = 0):
    QEC_dict = {}
    k_dict = {}
    print Qubit

    if sweep_time == False:
        len_k =2
    elif sweep_time == True:
        len_k = 4
    print len_k

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

    for k in range(len_k):
        if run !=0:
            k = k+2
        print 'k = '+str(k)
        print 'single_Qbit_sweep_time' +'_'+el_RO+'_'+ carbon +'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(-1)+'_'+state
        if sweep_time == False:
            timestamp, folder = toolbox.latest_data(contains = 'single_Qbit' +'_'+el_RO+'_'+ carbon +'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_sign)+'_'+state, older_than = older_than,return_timestamp = True)
        elif sweep_time == True and k !=4:
            # print 'ok'
            timestamp, folder = toolbox.latest_data(contains = 'single_Qbit_sweep_time' +'_'+el_RO+'_'+ carbon +'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(-1)+'_'+state, older_than = '20150105_142000',return_timestamp = True)
        if sweep_time == True and  k == 4:
            timestamp, folder = toolbox.latest_data(contains = 'single_Qbit_sweep_time' +'_'+el_RO+'_'+ carbon +'_RO'+str(RO)+'_k'+str(0)+'_sign'+ str(-1)+'_'+state, older_than = '20150106_100000',return_timestamp = True)
            print 'k =4, timestamp = '+str(timestamp)
        SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
        print timestamp
        print SSRO_folder
        k_dict['k_'+str(k)] ={}
        k_dict['k_'+str(k)] = load_QEC_data(folder, SSRO_folder, post_select = False)

    for item in k_dict['k_'+str(k)]:
        if len_k ==2 and run ==0:
            QEC_dict[item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item]), axis=0)
        if len_k ==2 and run !=0:
            QEC_dict[item] = np.concatenate((k_dict['k_2'][item],k_dict['k_3'][item]), axis=0)            
        elif len_k == 4 and run !=0:
            QEC_dict[item] = np.concatenate((k_dict['k_4'][item],k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item],k_dict['k_3'][item]), axis=0)
        elif len_k == 4:
            QEC_dict[item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item],k_dict['k_3'][item]), axis=0)
    return QEC_dict,folder

''' simple plotting QEC data without loading/saving '''

def QEC_data_single_state_RO(older_than = None,state = 'Z',RO = 0, sym = '00', len_k = 6):

    QEC_data_dict = {}
    u_list = ['c0_u', 'c0_00_u','c0_01_u','c0_10_u','c0_11_u']
    c_list = ['c0', 'c0_00','c0_01','c0_10','c0_11']
    p_list = ['p00','p01','p10','p11']
    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']


    QEC_dict, folder = QEC_create_data_dict(older_than = older_than, RO = RO, state = state, sym = sym, len_k = len_k)
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

    QEC_data_dict['y_pos']      = (QEC_dict[str(1)]['positive']['c0']-QEC_dict[str(1)]['negative']['c0'])/2
    QEC_data_dict['y_pos_err']  = (QEC_dict[str(1)]['positive']['c0_u']**2+QEC_dict[str(1)]['negative']['c0_u']**2)**0.5/2

    QEC_data_dict['y_neg']      = (QEC_dict[str(-1)]['positive']['c0']-QEC_dict[str(-1)]['negative']['c0'])/2
    QEC_data_dict['y_neg_err']  = (QEC_dict[str(-1)]['positive']['c0_u']**2+QEC_dict[str(-1)]['negative']['c0_u']**2)**0.5/2


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

def QEC_plot_single_state_RO(older_than = None, no_error = '00',state = 'Z',RO = 0, e_sign = None,plot_guide = True, len_k = 6):


    if e_sign == None:
        QEC_data_dict, folder =  QEC_data_single_state_RO(older_than = older_than,state = state,RO = RO, sym = no_error, len_k=len_k)
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

    plt.show()

    return QEC_data_dict, folder


''' fitting the data '''

### Fitfunctions

def fit_QEC(g_A, g_pc,g_O):
    '''Fit function for QEC process fidelity data
    g_O -  Offset, given by the fidelity of the state that is insensitive to errors
    g_A -  Amplitude, g_iven by the fidelity of the states that are sensitive
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''A*[pc(1-6p**2+4p**3)+(1-pc)(1-2p)]'''

    A   = fit.Parameter(g_A, 'A')
    O = fit.Parameter(g_O,'O')
    pc   = fit.Parameter(g_pc, 'pc')

    p0 = [A, pc,O]

    def fitfunc(p):
        '''test'''
        return (O()+A()*(pc()*(1-6*p**2+4*p**3)+(1-pc())*(1-2*p)))

    return p0, fitfunc, fitfunc_str

def fit_QEC_curve(x,y, return_errorbar=False,return_guess = False):

    guess_A = 0.85
    guess_pc = 0.5

    p0, fitfunc, fitfunc_str = fit_QEC( guess_A, guess_pc,0)

    if return_guess == True:
        x_temp      = np.linspace(x[0],x[-1],200)
        y_temp      =  fitfunc(x_temp) 
        return x_temp, y_temp

    else:

        fit_result = fit.fit1d(x, y, fit_QEC,
                guess_A, guess_pc,0,
                fixed=[2],
                do_print=True, ret=True)

        p02, fitfunc2, fitfunc_str2 = fit_QEC(fit_result['params_dict']['A'], fit_result['params_dict']['pc'],0)

        x_temp      = np.linspace(min(x),max(x),200)
        y_temp      =  fitfunc2(x_temp)

        print 'A'
        print fit_result['params_dict']['A']
        if return_errorbar == False:
            return x_temp, y_temp, fit_result['params_dict']['pc']
        else:
            return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['error_dict']['pc']

def fit_QEC_process_curve(x,y, A = 0.5, pc = 1, O = 0.5, return_errorbar=False,return_guess = False):

    guess_O = O
    guess_A = A
    guess_p = pc

    p0, fitfunc, fitfunc_str = fit_QEC(guess_A, guess_p,guess_O)

    if return_guess == True:
        x_temp      = np.linspace(x[0],x[-1],200)
        y_temp      =  fitfunc(x_temp) 
        return x_temp, y_temp

    else:
        fit_result = fit.fit1d(x, y, fit_QEC,
             guess_A, guess_p,guess_O,
                fixed=[],
                do_print=True, ret=True)

        p02, fitfunc2, fitfunc_str2 = fit_QEC( fit_result['params_dict']['A'], fit_result['params_dict']['pc'],g_O = fit_result['params_dict']['O'])

        x_temp      = np.linspace(min(x),max(x),200)
        y_temp      =  fitfunc2(x_temp)

        if return_errorbar == False:
            return x_temp, y_temp,fit_result['params_dict']['pc']
        elif return_errorbar == 'all':
            return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['params_dict']['A'],fit_result['params_dict']['O']
        else:
            return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['error_dict']['pc']



def fit_QEC_11(g_A, g_pc,g_O,F0=0.890,F1 = 0.988):
    '''Fit function for QEC process fidelity data
    g_O -  Offset, given by the fidelity of the state that is insensitive to errors
    g_A -  Amplitude, g_iven by the fidelity of the states that are sensitive
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''A*[pc(1-6p**2+4p**3)+(1-pc)(1-2p)]'''

    A   = fit.Parameter(g_A, 'A')
    O = fit.Parameter(g_O,'O')
    pc   = fit.Parameter(g_pc, 'pc')

    p0 = [A, pc,O]

    def fitfunc(p):
        '''test'''
        return (O()+A()*(pc()*(1-6*p**2+4*p**3)+(1-pc())*(1-2*p))*(F1**2-(3*F1**2-F0**2-2*F0*F1)*(p-p**2)))

    return p0, fitfunc, fitfunc_str


def fit_QEC_process_curve_11(x,y, A=0.5,pc=1,O=0.5,return_errorbar=False,return_guess=False):

    guess_O = O
    guess_A = A
    guess_p = pc
    p0, fitfunc, fitfunc_str = fit_QEC_11(guess_A, guess_p,guess_O)

    if return_guess == True:
        x_temp      = np.linspace(x[0],x[-1],200)
        y_temp      =  fitfunc(x_temp) 
        return x_temp, y_temp

    else:
        fit_result = fit.fit1d(x, y, fit_QEC_11,
             guess_A, guess_p,guess_O,
                fixed=[],
                do_print=True, ret=True)

        p02, fitfunc2, fitfunc_str2 = fit_QEC_11( fit_result['params_dict']['A'], fit_result['params_dict']['pc'],g_O = fit_result['params_dict']['O'])

        x_temp      = np.linspace(min(x),max(x),200)
        y_temp      =  fitfunc2(x_temp)

        if return_errorbar == False:
            return x_temp, y_temp,fit_result['params_dict']['pc']
        elif return_errorbar == 'all':
            return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['params_dict']['A'],fit_result['params_dict']['O']
        else:
            return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['error_dict']['pc']

def fit_QEC_curve_11(x,y, A=0.5,pc=1,O=0.5,return_errorbar=False,return_guess=False):

    guess_O = O
    guess_A = A
    guess_p = pc
    p0, fitfunc, fitfunc_str = fit_QEC_11(guess_A, guess_p,guess_O)

    if return_guess == True:
        x_temp      = np.linspace(x[0],x[-1],200)
        y_temp      =  fitfunc(x_temp) 
        return x_temp, y_temp

    else:
        fit_result = fit.fit1d(x, y, fit_QEC_11,
             guess_A, guess_p,0,
                fixed=[2],
                do_print=True, ret=True)

        p02, fitfunc2, fitfunc_str2 = fit_QEC_11( fit_result['params_dict']['A'], fit_result['params_dict']['pc'],g_O = 0)

        x_temp      = np.linspace(min(x),max(x),200)
        y_temp      =  fitfunc2(x_temp)

        if return_errorbar == False:
            return x_temp, y_temp,fit_result['params_dict']['pc']
        elif return_errorbar == 'all':
            return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['params_dict']['A'],fit_result['params_dict']['O']
        else:
            return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['error_dict']['pc']

def fit_QEC_11_2_alt(g_A, g_pc,g_O,F0=0.890,F1 = 0.988):
    '''Fit function for QEC process fidelity data
    g_O -  Offset, given by the fidelity of the state that is insensitive to errors
    g_A -  Amplitude, g_iven by the fidelity of the states that are sensitive
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''A*[pc(1-6p**2+4p**3)+(1-pc)(1-2p)]'''

    A   = fit.Parameter(g_A, 'A')
    O = fit.Parameter(g_O,'O')
    pc   = fit.Parameter(g_pc, 'pc')

    p0 = [A, pc,O]

    def fitfunc(p):
        '''test'''
        p2 = 1/2.*(1-np.sqrt(1-2*p))
        return (O()+A()*(pc()*(1-6*p**2+4*p**3)+(1-pc())*(1-2*p))*(F1**2-(3*F1**2-F0**2-2*F0*F1)*(p2-p2**2)))

    return p0, fitfunc, fitfunc_str

def fit_QEC_curve_11_2(x,y, A=0.5,pc=1,O=0.5,return_errorbar=False,return_guess=False):

    guess_O = O
    guess_A = A
    guess_p = pc
    p0, fitfunc, fitfunc_str = fit_QEC_11_2_alt(guess_A, guess_p,guess_O)

    if return_guess == True:
        x_temp      = np.linspace(x[0],x[-1],200)
        y_temp      =  fitfunc(x_temp) 
        return x_temp, y_temp

    else:
        fit_result = fit.fit1d(x, y, fit_QEC_11_2_alt,
             guess_A, guess_p,0,
                fixed=[2],
                do_print=True, ret=True)

        p02, fitfunc2, fitfunc_str2 = fit_QEC_11_2_alt( fit_result['params_dict']['A'], fit_result['params_dict']['pc'],g_O = 0)

        x_temp      = np.linspace(min(x),max(x),200)
        y_temp      =  fitfunc2(x_temp)

        if return_errorbar == False:
            return x_temp, y_temp,fit_result['params_dict']['pc']
        elif return_errorbar == 'all':
            return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['params_dict']['A'],fit_result['params_dict']['O']
        else:
            return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['error_dict']['pc']

def fit_QEC_11_3_alt(g_A, g_pc,g_O,F0=0.890,F1 = 0.988):
    '''Fit function for QEC process fidelity data
    g_O -  Offset, given by the fidelity of the state that is insensitive to errors
    g_A -  Amplitude, g_iven by the fidelity of the states that are sensitive
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''A*[pc(1-6p**2+4p**3)+(1-pc)(1-2p)]'''

    A   = fit.Parameter(g_A, 'A')
    O = fit.Parameter(g_O,'O')
    pc   = fit.Parameter(g_pc, 'pc')

    p0 = [A, pc,O]

    def fitfunc(p):
        '''test'''
        p3 = 1/2.*(1-(1-2*p)**(1/3.))
        return (O()+A()*(pc()*(1-6*p**2+4*p**3)+(1-pc())*(1-2*p))*(F1**2-(3*F1**2-F0**2-2*F0*F1)*(p3-p3**2))**2)

    return p0, fitfunc, fitfunc_str

def fit_QEC_curve_11_3(x,y, A=0.5,pc=1,O=0.5,return_errorbar=False,return_guess=False):

    guess_O = O
    guess_A = A
    guess_p = pc
    p0, fitfunc, fitfunc_str = fit_QEC_11_3_alt(guess_A, guess_p,guess_O)

    if return_guess == True:
        x_temp      = np.linspace(x[0],x[-1],200)
        y_temp      =  fitfunc(x_temp) 
        return x_temp, y_temp

    else:
        fit_result = fit.fit1d(x, y, fit_QEC_11_3_alt,
             guess_A, guess_p,0,
                fixed=[2],
                do_print=True, ret=True)

        p02, fitfunc2, fitfunc_str2 = fit_QEC_11_3_alt( fit_result['params_dict']['A'], fit_result['params_dict']['pc'],g_O = 0)

        x_temp      = np.linspace(min(x),max(x),200)
        y_temp      =  fitfunc2(x_temp)

        if return_errorbar == False:
            return x_temp, y_temp,fit_result['params_dict']['pc']
        elif return_errorbar == 'all':
            return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['params_dict']['A'],fit_result['params_dict']['O']
        else:
            return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['error_dict']['pc']
def fit_QEC_00(g_A, g_pc,g_O,F0=0.890,F1 = 0.988):
    '''Fit function for QEC process fidelity data
    g_O -  Offset, given by the fidelity of the state that is insensitive to errors
    g_A -  Amplitude, g_iven by the fidelity of the states that are sensitive
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''A*[pc(1-6p**2+4p**3)+(1-pc)(1-2p)]'''

    A   = fit.Parameter(g_A, 'A')
    O = fit.Parameter(g_O,'O')
    pc   = fit.Parameter(g_pc, 'pc')

    p0 = [A, pc,O]

    def fitfunc(p):
        '''test'''
        return (O()+A()*(pc()*(1-6*p**2+4*p**3)+(1-pc())*(1-2*p))*(F0**2-(3*F0**2-F1**2-2*F0*F1)*(p-p**2)))

    return p0, fitfunc, fitfunc_str


def fit_QEC_process_curve_00(x,y, return_errorbar=False):

    guess_O = 0.5
    guess_A = 0.5
    guess_p = 1
    p0, fitfunc, fitfunc_str = fit_QEC_00(guess_A, guess_p,guess_O)

    fit_result = fit.fit1d(x, y, fit_QEC_00,
         guess_A, guess_p,guess_O,
            fixed=[],
            do_print=True, ret=True)

    p02, fitfunc2, fitfunc_str2 = fit_QEC_00( fit_result['params_dict']['A'], fit_result['params_dict']['pc'],g_O = fit_result['params_dict']['O'])

    x_temp      = np.linspace(min(x),max(x),200)
    y_temp      =  fitfunc2(x_temp)

    if return_errorbar == False:
        return x_temp, y_temp,fit_result['params_dict']['pc']
    else:
        return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['error_dict']['pc']

def fit_QEC_01(g_A, g_pc,g_O,F0=0.890,F1 = 0.988):
    '''Fit function for QEC process fidelity data
    g_O -  Offset, given by the fidelity of the state that is insensitive to errors
    g_A -  Amplitude, g_iven by the fidelity of the states that are sensitive
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''A*[pc(1-6p**2+4p**3)+(1-pc)(1-2p)]'''

    A   = fit.Parameter(g_A, 'A')
    O = fit.Parameter(g_O,'O')
    pc   = fit.Parameter(g_pc, 'pc')

    p0 = [A, pc,O]

    def fitfunc(p):
        '''test'''
        return (O()+A()*(pc()*(1-6*p**2+4*p**3)+(1-pc())*(1-2*p))*(F0*F1*(1-2*p+2*p**2)+(F0**2+F1**2)*(p-p**2)))

    return p0, fitfunc, fitfunc_str


def fit_QEC_process_curve_01(x,y, return_errorbar=False):

    guess_O = 0.5
    guess_A = 0.5
    guess_p = 1
    p0, fitfunc, fitfunc_str = fit_QEC_01(guess_A, guess_p,guess_O)

    fit_result = fit.fit1d(x, y, fit_QEC_01,
         guess_A, guess_p,guess_O,
            fixed=[],
            do_print=True, ret=True)

    p02, fitfunc2, fitfunc_str2 = fit_QEC_01( fit_result['params_dict']['A'], fit_result['params_dict']['pc'],g_O = fit_result['params_dict']['O'])

    x_temp      = np.linspace(min(x),max(x),200)
    y_temp      =  fitfunc2(x_temp)

    if return_errorbar == False:
        return x_temp, y_temp,fit_result['params_dict']['pc']
    else:
        return x_temp, y_temp,fit_result['params_dict']['pc'],fit_result['error_dict']['pc']


def fit_QEC_2_rounds2(g_A, g_P, g_P2):
    '''Fit function for QEC process fidelity data
    g_C -  Initial contrast
    g_A -  Extra uncorrected/introduced errors
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''test'''

    A   = fit.Parameter(g_A, 'A')
    P   = fit.Parameter(g_P, 'P')
    # P2   = fit.Parameter(g_P2, 'P2')

    p0 = [A,P]

    def fitfunc(p):
        '''test'''

        p2 = 1/2.*(1-np.sqrt(1-2*p))
        # return (A()*(P()*P()*(1-6*p2**2+4*p2**3)**2+
        #     P()*(1-P())*(1-6*p2**2+4*p2**3)*(1-2*p2)+
        #     (1-P())*P()*(1-6*p**2+4*p**3)
        #     +(1-P())*(1-P())*(1-2*p)))
        return (A()*(P()*(1-6*p2**2+4*p2**3)**2+
            (1-P())*(1-2*p)))

    return p0, fitfunc, fitfunc_str

def fit_QEC_2_rounds_curve2(x,y, return_errorbar=False,plot_guess = False):

    
    guess_A = 0.7
    guess_P = 1
    guess_P2 = 1

    p0, fitfunc, fitfunc_str = fit_QEC_2_rounds2(guess_A, guess_P, guess_P2)

    x_temp      = np.linspace(x[0],x[-1],200)

    if plot_guess == True:
        y_temp      =  fitfunc(x_temp)

        return x_temp, y_temp
    else:
        fit_result = fit.fit1d(x, y, fit_QEC_2_rounds2,
                guess_A, guess_P, guess_P2,
                fixed=[],
                do_print=True, ret=True)

        p02, fitfunc2, fitfunc_str2 = fit_QEC_2_rounds2( fit_result['params_dict']['A'],fit_result['params_dict']['P'], guess_P2)

        x_temp      = np.linspace(x[0],x[-1],200)
        y_temp      =  fitfunc2(x_temp)        
        if return_errorbar == False:
            return x_temp, y_temp,fit_result['params_dict']['p']
        else:
            return x_temp, y_temp,fit_result['params_dict']['P'],fit_result['error_dict']['P']


def fit_QEC_2_rounds2_11(g_A, g_P, g_P2,F0=0.890,F1 = 0.988):
    '''Fit function for QEC process fidelity data
    g_C -  Initial contrast
    g_A -  Extra uncorrected/introduced errors
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''test'''

    A   = fit.Parameter(g_A, 'A')
    P   = fit.Parameter(g_P, 'P')
    # P2   = fit.Parameter(g_P2, 'P2')

    p0 = [A,P]

    def fitfunc(p):
        '''test'''

        p2 = 1/2.*(1-np.sqrt(1-2*p))
        return (A()*(P()*(1-6*p2**2+4*p2**3)**2+
           +(1-P())*(1-2*p))*(F1**2-(3*F1**2-F0**2-2*F0*F1)*(p2-p2**2)))

    return p0, fitfunc, fitfunc_str

def fit_QEC_2_rounds_curve2_11(x,y, return_errorbar=False,plot_guess = False):

    
    guess_A = 0.7
    guess_P = 1
    guess_P2 = 1

    p0, fitfunc, fitfunc_str = fit_QEC_2_rounds2_11(guess_A, guess_P, guess_P2)

    x_temp      = np.linspace(x[0],x[-1],200)

    if plot_guess == True:
        y_temp      =  fitfunc(x_temp)

        return x_temp, y_temp
    else:
        fit_result = fit.fit1d(x, y, fit_QEC_2_rounds2_11,
                guess_A, guess_P, guess_P2,
                fixed=[],
                do_print=True, ret=True)

        p02, fitfunc2, fitfunc_str2 = fit_QEC_2_rounds2_11( fit_result['params_dict']['A'],fit_result['params_dict']['P'], guess_P2)

        x_temp      = np.linspace(x[0],x[-1],200)
        y_temp      =  fitfunc2(x_temp)        
        if return_errorbar == False:
            return x_temp, y_temp,fit_result['params_dict']['p']
        else:
            return x_temp, y_temp,fit_result['params_dict']['P'],fit_result['error_dict']['P']

def fit_QEC_3_rounds2(g_A,g_P):
    '''Fit function for QEC process fidelity data
    g_C -  Initial contrast
    g_A -  Extra uncorrected/introduced errors
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''test'''

    A   = fit.Parameter(g_A, 'A')
    P   = fit.Parameter(g_P, 'P')


    p0 = [A,P]

    def fitfunc(p):
        '''test'''
        p3 = 1/2.*(1-(1-2*p)**(1/3.))
        p2 = 1/2.*(1-(1-2*p)**(2/3.))
        # return (A()*(P()*P()*P()*(1-6*p3**2+4*p3**3)**3+
        #             P()*P()*(1-P())*(1-6*p3**2+4*p3**3)**2*(1-2*p3)+
        #             P()*(1-P())*P()*(1-6*p3**2+4*p3**3)*(1-6*p2**2+4*p2**3)+
        #              (1-P())*P()*P()*(1-6*p2**2+4*p2**3)*(1-6*p3**2+4*p3**3)+
        #              P()*(1-P())*(1-P())*(1-6*p3**2+4*p3**3)*(1-2*p2)+
        #              (1-P())*P()*(1-P())*(1-6*p2**2+4*p2**3)*(1-2*p3)+
        #                  (1-P())*(1-P())*P()*(1-6*p**2+4*p**3)+
        #                     (1-P())*(1-P())*(1-P())*(1-2*p)))

        return (A()*(P()*(1-6*p3**2+4*p3**3)**3+
                    (1-P())*(1-2*p)))
        
    return p0, fitfunc, fitfunc_str

def fit_QEC_3_rounds_curve2(x,y, return_errorbar=False,plot_guess = False):

    guess_A = 0.75
    guess_P = 1


  
    x_temp      = np.linspace(x[0],x[-1],200)
    
    p0, fitfunc, fitfunc_str = fit_QEC_3_rounds2( guess_A, guess_P)

    if plot_guess == True:
        y_temp      =  fitfunc(x_temp)

        return x_temp, y_temp
    else:
        
        fit_result = fit.fit1d(x, y, fit_QEC_3_rounds2,
                guess_A,  guess_P,
                fixed=[],
                do_print=True, ret=True)
        p02, fitfunc2, fitfunc_str2 = fit_QEC_3_rounds2( fit_result['params_dict']['A'], fit_result['params_dict']['P'])

        y_temp      =  fitfunc2(x_temp)

        if return_errorbar == False:
            return x_temp, y_temp,fit_result['params_dict']['P1']
        else:
            return x_temp, y_temp,fit_result['params_dict']['P'],fit_result['error_dict']['P']            


def fit_QEC_3_rounds2_11(g_A,g_P,F0=0.890,F1 = 0.988):
    '''Fit function for QEC process fidelity data
    g_C -  Initial contrast
    g_A -  Extra uncorrected/introduced errors
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''test'''

    A   = fit.Parameter(g_A, 'A')
    P   = fit.Parameter(g_P, 'P')


    p0 = [A,P]

    def fitfunc(p):
        '''test'''
        p3 = 1/2.*(1-(1-2*p)**(1/3.))
        p2 = 1/2.*(1-(1-2*p)**(2/3.))
        return (A()*(P()*(1-6*p3**2+4*p3**3)**3
                    +(1-P())*(1-2*p))*(F1**2-(3*F1**2-F0**2-2*F0*F1)*(p3-p3**2))**2)

    return p0, fitfunc, fitfunc_str

def fit_QEC_3_rounds_curve2_11(x,y, return_errorbar=False,plot_guess = False):

    guess_A = 0.75
    guess_P = 1


  
    x_temp      = np.linspace(x[0],x[-1],200)
    
    p0, fitfunc, fitfunc_str = fit_QEC_3_rounds2_11( guess_A, guess_P)

    if plot_guess == True:
        y_temp      =  fitfunc(x_temp)

        return x_temp, y_temp
    else:
        
        fit_result = fit.fit1d(x, y, fit_QEC_3_rounds2_11,
                guess_A,  guess_P,
                fixed=[],
                do_print=True, ret=True)
        p02, fitfunc2, fitfunc_str2 = fit_QEC_3_rounds2_11( fit_result['params_dict']['A'], fit_result['params_dict']['P'])

        y_temp      =  fitfunc2(x_temp)

        if return_errorbar == False:
            return x_temp, y_temp,fit_result['params_dict']['P1']
        else:
            return x_temp, y_temp,fit_result['params_dict']['P'],fit_result['error_dict']['P']                        

def fit_no_error_detection(g_Pin):

    '''Fit function for QEC error detection curve 'no error'
    g_pin - guess for initial error on encoded state
    '''

    fitfunc_str = '''test'''

    Pin   = fit.Parameter(g_Pin, 'Pin')


    p0 = [Pin]

    def fitfunc(p):
        '''test'''

        return 1 - 3 *(p + Pin() - 2 *p *Pin()) + 3* (p + Pin() - 2 *p *Pin())**2

    return p0, fitfunc, fitfunc_str    


def fit_error_detection(g_Pin):

    '''Fit function for QEC error detection curve 'no error'
    g_pin - guess for initial error on encoded state
    '''

    fitfunc_str = '''test'''

    Pin   = fit.Parameter(g_Pin, 'Pin')


    p0 = [Pin]

    def fitfunc(p):
        '''test'''

        return (Pin() + p - 2*Pin()*p) - (Pin() + p - 2*Pin()*p)**2

    return p0, fitfunc, fitfunc_str   


def fit_no_error_detection_curve(x,y,pin = 0.127, return_errorbar=False,plot_guess = False):

    
    guess_Pin = pin


    p0, fitfunc, fitfunc_str = fit_no_error_detection(guess_Pin)

    x_temp      = np.linspace(x[0],x[-1],200)

    if plot_guess == True:
        y_temp      =  fitfunc(x_temp)
        # print x_temp
        # print y_temp
        return x_temp, y_temp
    else:
        fit_result = fit.fit1d(x, y, fit_no_error_detection, guess_Pin,
                fixed=[0],
                do_print=True, ret=True)

        p02, fitfunc2, fitfunc_str2 = fit_no_error_detection( fit_result['params_dict']['Pin'])

        x_temp      = np.linspace(x[0],x[-1],200)
        y_temp      =  fitfunc2(x_temp) 

        if return_errorbar == False:
            return x_temp, y_temp,fit_result['params_dict']['Pin']
        else:
            return x_temp, y_temp,fit_result['params_dict']['Pin'],fit_result['error_dict']['Pin']

def fit_error_detection_curve(x,y, pin = 0.127 ,return_errorbar=False,plot_guess = False):

    
    guess_Pin = pin


    p0, fitfunc, fitfunc_str = fit_error_detection(guess_Pin)

    x_temp      = np.linspace(x[0],x[-1],200)

    if plot_guess == True:
        y_temp      =  fitfunc(x_temp)

        return x_temp, y_temp
    else:
        fit_result = fit.fit1d(x, y, fit_error_detection,
                guess_Pin,
                fixed=[],
                do_print=True, ret=True)

        p02, fitfunc2, fitfunc_str2 = fit_error_detection( fit_result['params_dict']['Pin'])

        x_temp      = np.linspace(x[0],x[-1],200)
        y_temp      =  fitfunc2(x_temp) 
               
        if return_errorbar == False:
            return x_temp, y_temp,fit_result['params_dict']['Pin']
        else:
            return x_temp, y_temp,fit_result['params_dict']['Pin'],fit_result['error_dict']['Pin']

def fit_timesweep_QEC_1round(g_A,g_T1,g_T2,g_T3,g_p1,g_p2,g_p3):

    '''Fit function for QEC process fidelity data
    g_A -  Amplitude
    g_T1 - Decay constant qubit 1
    g_T2 - Decay constant qubit 2
    g_T3 - Decay constant qubit 3
    '''

    fitfunc_str = ''''''

    A   = fit.Parameter(g_A, 'A')
    T1   = fit.Parameter(g_T1, 'T1')
    T2   = fit.Parameter(g_T2, 'T2')
    T3   = fit.Parameter(g_T3, 'T3')
    p1   = fit.Parameter(g_p1, 'p1')
    p2   = fit.Parameter(g_p2, 'p2')
    p3   = fit.Parameter(g_p3, 'p3')

    p0 = [A, T1,T2,T3,p1,p2,p3]

    print g_p1
    print g_p2
    print g_p3

    # print g_T1*1e3
    # print g_T2*1e3
    # print g_T3*1e3

    def fitfunc(t):
        '''test'''
        return A()*(np.exp(-t**2*(1/T1()**2+1/T2()**2+1/T3()**2))*(1-p1()/2.-p2()/2.-p3()/2.)
                    +p1()/2.*(-np.exp(-(t/T1())**2)+np.exp(-(t/T2())**2)+np.exp(-(t/T3())**2))
                    +p2()/2.*(+np.exp(-(t/T1())**2)-np.exp(-(t/T2())**2)+np.exp(-(t/T3())**2))
                    +p3()/2.*(+np.exp(-(t/T1())**2)+np.exp(-(t/T2())**2)-np.exp(-(t/T3())**2)))

    return p0, fitfunc, fitfunc_str


def fit_timesweep_QEC_1round_curve(x,y,A = 0.87, T1 = 9e-3,T2 = 10e-3, T3 = 19e-3,p1 = 0.95, p2 = 0.9, p3 = 0.85, return_errorbar=False,plot_guess = False):

    guess_A    = A
    guess_p1   = p1
    guess_p2   = p2
    guess_p3   = p3

    p0, fitfunc, fitfunc_str = fit_timesweep_QEC_1round(guess_A, T1, T2, T3, guess_p1, guess_p2, guess_p3)

    fit_result = fit.fit1d(x, y, fit_timesweep_QEC_1round,
            guess_A, T1,T2,T3,guess_p1, guess_p2, guess_p3,
            fixed=[1,2,3],
            do_print=True, ret=True)

    p02, fitfunc2, fitfunc_str2 = fit_timesweep_QEC_1round(fit_result['params_dict']['A'], T1,T2,T3,fit_result['params_dict']['p1'], fit_result['params_dict']['p2'], fit_result['params_dict']['p3'])

    x_temp      = np.linspace(min(x),max(x),200)
    y_temp      =  fitfunc2(x_temp)

    if plot_guess == True:
        y_temp      =  fitfunc(x_temp)
        print 'test'
        # print x_temp
        # print y_temp
        return x_temp, y_temp
    else:
        return x_temp, y_temp

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

def fit_timesweep_single(x,y,a = 0, A = 0.4, x0 =0, T = 21,n = 2, return_errorbar=False,plot_guess = False):

        g_a =a
        g_A =A
        g_x0 =x0
        g_T =T
        g_n =n

        p0, fitfunc, fitfunc_str = fit_general_exponential(g_a, g_A, g_x0, g_T, g_n)

        fit_result = fit.fit1d(x, y, fit_general_exponential,
                a, A, x0,T,n,
                fixed=[0,2],
                do_print=True, ret=True)

        p02, fitfunc2, fitfunc_str2 = fit_general_exponential(a,
                        fit_result['params_dict']['A'], x0, fit_result['params_dict']['T'],fit_result['params_dict']['n'])

        print (x[-1]-x[0])*1e3
        x_temp      = np.linspace(x[0],x[-1],(x[-1]-x[0])*1e3)
        y_temp      =  fitfunc2(x_temp)

        if plot_guess == True:
            y_temp      =  fitfunc(x_temp)
            return x_temp, y_temp
        else:
            if return_errorbar == False:
                return x_temp, y_temp,fit_result['params_dict']['T']
            else:
                return x_temp, y_temp,fit_result['params_dict']['T'],fit_result['error_dict']['T']


''' these functions are used to open/close save/load from and to HDF5 files '''

def openfile(name = ''):
   datafile = h5py.File(os.path.join(r'K:\ns\qt\Diamond\Projects\QEC LT\QEC data10', name))
   return datafile

def closefile(datafile):
    datafile.close()

def openfile_single_state_RO_run(sym = '00', RO = 1, state = 'Z', error_sign = 1,el_RO = 'positive', run = 1):
    name = 'run_'+str(run)+ '_error_sym_' + sym + '_state_'+state + '_RO_'+str(RO)+ '_sign'+str(error_sign)+'_el_RO_'+el_RO+'.hdf5'
    # print name
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

def save_QEC_dataset_single_sign_single_elRO(do_p = False,older_than = None,sym = '11', RO = 1, state = 'Z', error_sign = 1, el_RO = 'positive', run =1, sweep_time =False ):


    QEC_temp_dict, folder = QEC_create_data_dict_single_error_single_elRO(do_p = do_p,older_than = older_than, RO = RO, state = state,
                                                                        len_k = 6, sym = sym,error_sign = error_sign, el_RO = el_RO,sweep_time = sweep_time,run = run)

    if sweep_time == False:
        print 'yup'
        datafile = openfile_single_state_RO_run(sym = sym, RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)
    elif sweep_time ==True:
        datafile = openfile_single_state_RO_run(sym = sym  +'sweep_time_', RO = RO, state = state, error_sign = 0,el_RO = el_RO, run = run)

    save_data_single_state_RO_hdf5file(datafile, QEC_temp_dict)

    closefile(datafile)

def save_no_QEC_dataset_single_sign_single_elRO(run = 0,do_p = False,idle = False, sweep_time = False, older_than = None, RO = 1, state = 'Z', error_sign = 1, el_RO = 'positive'):


    QEC_temp_dict, folder = no_QEC_create_data_dict_single_error_single_elRO(do_p = do_p,older_than = older_than, RO = RO, idle = idle, sweep_time = sweep_time,
                             state = state, error_sign = error_sign, el_RO = el_RO, run = run)


    if idle == True:
        datafile = openfile_single_state_RO_run(sym = 'no_correction_idle', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)
    elif sweep_time == True:
        datafile = openfile_single_state_RO_run(sym = 'no_correction_sweep_time', RO = RO, state = state, error_sign = 0,el_RO = el_RO, run = 1)
    else:
        datafile = openfile_single_state_RO_run(sym = 'no_correction', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)
    save_data_single_state_RO_hdf5file(datafile, QEC_temp_dict)

    closefile(datafile)

def save_no_QEC_extra_dataset_single_sign_single_elRO(run = 4,idle = False, sweep_time = False, older_than = None, RO = 1, state = 'Z', error_sign = 1, el_RO = 'positive'):


    QEC_temp_dict = no_QEC_extra_data_single_state_RO(run = run,older_than = older_than,
                                                       state = state,RO = RO, load_set = False)
    datafile = openfile_single_state_RO_run(sym = 'no_correction', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)
    save_data_single_state_RO_hdf5file(datafile, QEC_temp_dict)

    closefile(datafile)

def save_single_Qubit_QEC_dataset_single_sign_single_elRO(run = 0,older_than = None, Qubit =1, sweep_time = False,state = 'Z', error_sign = 1, el_RO = 'positive'):

    if state == 'X' or state == 'mX':
        RO = 2+(Qubit-1)*3
    elif state == 'Y' or state == 'mY':
        RO = 1+(Qubit-1)*3
    elif state == 'Z' or state == 'mZ':
        RO = 0+(Qubit-1)*3
        if sweep_time == True:
            RO = Qubit-1

    QEC_temp_dict, folder = single_Qubit_QEC_create_data_dict_single_error_single_elRO(older_than = older_than, Qubit = Qubit,
                                                            state = state, error_sign = error_sign, el_RO = el_RO, sweep_time =sweep_time,run = run)
    if sweep_time == False:
        datafile = openfile_single_state_RO_run(sym = 'single_qubit', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)
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

def load_no_QEC_dataset_single_sign_single_elRO(idle = False, sweep_time = False,RO = 1, state = 'Z', error_sign = 1, el_RO = 'positive', run = 0):

    if idle == False and sweep_time == False:
        datafile = openfile_single_state_RO_run(sym = 'no_correction', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)
    elif sweep_time == True:

        datafile = openfile_single_state_RO_run(sym = 'no_correction_sweep_time', RO = RO, state = state, error_sign = 0,el_RO = el_RO, run = 1)
    else:
        # print 'oke'
        # print run
        datafile = openfile_single_state_RO_run(sym = 'no_correction_idle', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)

    QEC_temp_dict = load_single_data_hdf5file(datafile)

    closefile(datafile)

    return QEC_temp_dict

def load_single_Qubit_QEC_dataset_single_sign_single_elRO(Qubit =1, state = 'Z', error_sign = 1,sweep_time = False, el_RO = 'positive',run = 0):

    if state == 'X' or state == 'mX':
        RO = 2+(Qubit-1)*3
    elif state == 'Y' or state == 'mY':
        RO = 1+(Qubit-1)*3
    elif state == 'Z' or state == 'mZ':
        RO = 0+(Qubit-1)*3
        if sweep_time == True:
            RO = Qubit-1

    if sweep_time == False:
        datafile = openfile_single_state_RO_run(sym = 'single_qubit', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)
    elif sweep_time == True:
        datafile = openfile_single_state_RO_run(sym = 'single_qubit_sweep_time', RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)

    QEC_temp_dict = load_single_data_hdf5file(datafile)

    closefile(datafile)

    return QEC_temp_dict

''' from here you can plot data taken from an existing HDF5 file '''

def QEC_sum_data_single_state_RO(run = 1, no_error = '00',state = 'Z',RO = 0,load_set = True, older_than = '201501010101',do_p = False,sweep_time = False):

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
            if load_set == True:
                QEC_dict[str(error_sign)][el_RO] = load_QEC_dataset_single_sign_single_elRO(sym = no_error, RO = RO, state = state, error_sign = error_sign,el_RO = el_RO, run = run)
            else:
                QEC_dict[str(error_sign)][el_RO] , folder = QEC_create_data_dict_single_error_single_elRO(run = run,older_than = older_than, RO = RO, state = state,
                                                                                                len_k = 6, sym = no_error,error_sign = error_sign, el_RO = el_RO, sweep_time = sweep_time,do_p = do_p)
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
            QEC_dict[str(error_sign)][el_RO] = load_QEC_dataset_single_sign_single_elRO(sym = no_error, RO = RO, state = state,
                                                error_sign = error_sign,el_RO = el_RO, run = run,sweep_time = sweep_time)
        elif load_set == False:
            QEC_dict[str(error_sign)][el_RO] , folder = QEC_create_data_dict_single_error_single_elRO(older_than = older_than, RO = RO, state = state,
                                                                                                len_k = 6, sym = no_error,error_sign = error_sign, el_RO = el_RO, sweep_time = sweep_time)


    for v in range(len(y_list)):
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

def no_QEC_data_single_state_RO(run = 0,idle = False, sweep_time = False, older_than = None,state = 'Z',RO = 0, load_set = False,do_p = False):

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
                # print 'ok'
                QEC_dict[str(error_sign)][el_RO], folder = no_QEC_create_data_dict_single_error_single_elRO(idle = idle,do_p = do_p, sweep_time = sweep_time, older_than = older_than, RO = RO, state = state, error_sign = error_sign, el_RO = el_RO,run = run)
            else:
                QEC_dict[str(error_sign)][el_RO]= load_no_QEC_dataset_single_sign_single_elRO(run = run,idle = idle, sweep_time = sweep_time, RO = RO, state = state, error_sign = error_sign,el_RO = el_RO)
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

def no_QEC_extra_data_single_state_RO(run = 4,idle = False, sweep_time = False, older_than = None,state = 'Z',RO = 0, load_set = False):

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

        QEC_dict[str(error_sign)][el_RO] = {}

        timestamp, folder = toolbox.latest_data(contains = 'no_correct' +'_'+el_RO+'_RO'+str(RO)+'_k'+str(3)+'_sign'+ str(error_sign)
                                                +'_'+state, older_than = older_than,return_timestamp = True)
               # print folder
        SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)

        QEC_dict[str(error_sign)][el_RO] = load_QEC_data(folder, SSRO_folder, post_select = False)


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

def no_QEC_sum_data_single_state_RO(run_list = [0,1,2],add_4 =False,idle = False, sweep_time = False, older_than = None,state = 'Z',RO = 0, load_set = False,do_weighted = False):

    QEC_dict = {}
    y_list = ['y']
    y_err_list = ['y_err']
    QEC_temp_dict = {}

    if do_weighted == True: #then we assume we want to take ALL data into account and not specific runs, so run 1,2,3 and 4!
        for i, run in enumerate(run_list):
            QEC_temp_dict[run] = no_QEC_data_single_state_RO(run = run,idle = idle, sweep_time = sweep_time, older_than = older_than,state = state,RO = RO, load_set = load_set)
           
        # a bit stupid but the error signs and el_RO were already summed up (sloppy....)
        four_dict = load_QEC_dataset_single_sign_single_elRO(sym = 'no_correction', RO = RO, state = state, error_sign = 1,el_RO = 'positive', run = 4)
                                    
        QEC_dict['y'] = np.zeros(12)
        QEC_dict['y_err'] = np.zeros(12)

        for ii in range(6):
            QEC_dict['y'][ii] = np.average([QEC_temp_dict[0]['y'][ii],QEC_temp_dict[1]['y'][ii],QEC_temp_dict[2]['y'][ii],QEC_temp_dict[3]['y'][ii]],
                                    weights =[1./(QEC_temp_dict[0]['y_err'][ii])**2,1./(QEC_temp_dict[1]['y_err'][ii])**2,1./(QEC_temp_dict[2]['y_err'][ii])**2,1./(QEC_temp_dict[3]['y_err'][ii])**2])
            QEC_dict['y_err'][ii] = 1./(1./(QEC_temp_dict[0]['y_err'][ii])**2+1./(QEC_temp_dict[1]['y_err'][ii])**2+1./(QEC_temp_dict[2]['y_err'][ii])**2+1./(QEC_temp_dict[3]['y_err'][ii])**2)**0.5

        QEC_dict['y'][5:6] = np.average([four_dict['y'][0:1],QEC_temp_dict[3]['y'][5:6]],
                                weights = [1./(four_dict['y_err'][0:1])**2,1./(QEC_temp_dict[3]['y_err'][5:6])**2])
        QEC_dict['y_err'][5] = 1./(1./(four_dict['y_err'][0])**2+1./(QEC_temp_dict[3]['y_err'][5])**2)**0.5

        for ii in range(6):
            ii = ii+6
            QEC_dict['y'][ii] = np.average([QEC_temp_dict[0]['y'][ii-1],QEC_temp_dict[1]['y'][ii-1],QEC_temp_dict[2]['y'][ii-1],QEC_temp_dict[3]['y'][ii]],
                                    weights = [1./(QEC_temp_dict[0]['y_err'][ii-1])**2,1./(QEC_temp_dict[1]['y_err'][ii-1])**2,1./(QEC_temp_dict[2]['y_err'][ii-1])**2,1./(QEC_temp_dict[3]['y_err'][ii])**2])

            QEC_dict['y_err'][ii] = 1./(1./(QEC_temp_dict[0]['y_err'][ii-1])**2+1./(QEC_temp_dict[1]['y_err'][ii-1])**2+1./(QEC_temp_dict[2]['y_err'][ii-1])**2+1./(QEC_temp_dict[3]['y_err'][ii])**2)**0.5

        
        QEC_dict['x'] = QEC_temp_dict[3]['x']

    elif do_weighted == False:
        for i, run in enumerate(run_list):
            if run !=3:
                # print 'ok'
                QEC_temp_dict = no_QEC_data_single_state_RO(run = run,idle = idle, sweep_time = sweep_time, older_than = older_than,state = state,RO = RO, load_set = load_set)
                if i == 0:
                    QEC_dict = {}
                    QEC_dict['y'] = 1/float(len(run_list))*QEC_temp_dict['y']
                    QEC_dict['y_err_temp'] = QEC_temp_dict['y_err']**2
                else:
                    QEC_dict['y'] += 1/float(len(run_list))*QEC_temp_dict['y']
                    QEC_dict['y_err_temp'] += QEC_temp_dict['y_err']**2

           

        
        if 3 in run_list and len(run_list)>1:
            QEC_dict['y_err'] = 1/float(len(run_list)-1)*QEC_dict['y_err_temp']**0.5
        elif 3 not in run_list:
            QEC_dict['y_err'] = 1/float(len(run_list))*QEC_dict['y_err_temp']**0.5

        if add_4 == True:
            # a bit stupid but the error signs and el_RO were already summed up (sloppy....)
            four_dict = load_QEC_dataset_single_sign_single_elRO(sym = 'no_correction', RO = RO, state = state, error_sign = 1,el_RO = 'positive', run = 4)
                                        
            for item in ['y','y_err']:
                if 3 not in run_list:
                    QEC_dict[item] = np.concatenate((QEC_dict[item][0:5],four_dict[item],QEC_dict[item][5:]), axis=0)
                elif 3 in run_list:
                    QEC_dict[item] = np.concatenate((QEC_dict[item][0:5],1/2.*four_dict[item],QEC_dict[item][5:]), axis=0)
            
        if 3 in run_list and len(run_list)>1:
            QEC_temp_dict = no_QEC_data_single_state_RO(run = run,idle = idle, sweep_time = sweep_time, older_than = older_than,state = state,RO = RO, load_set = load_set)
            QEC_dict['y'] +=  np.concatenate((1/float(len(run_list))*QEC_temp_dict['y'][0:5],1/2.*QEC_temp_dict['y'][5:6],1/float(len(run_list))*QEC_temp_dict['y'][6:]), axis=0)
            QEC_dict['y_err'] =  1/float(len(run_list))* (float(len(run_list)-1)*QEC_dict['y_err']**2+ QEC_temp_dict['y_err']**2)**0.5
        elif 3 in run_list and len(run_list)==1:
            QEC_temp_dict = no_QEC_data_single_state_RO(run = run,idle = idle, sweep_time = sweep_time, older_than = older_than,state = state,RO = RO, load_set = load_set)
            QEC_dict['y'] = QEC_temp_dict['y']
            QEC_dict['y_err'] = QEC_temp_dict['y_err']

        
        QEC_dict['x'] = QEC_temp_dict['x']
    
    return QEC_dict


def no_QEC_data_single_state_RO_single_error_sign(idle = False, sweep_time = False, older_than = None,state = 'Z',RO = 0,error_sign = 0, load_set = False,run = 0):

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

            QEC_dict[str(error_sign)][el_RO]= load_no_QEC_dataset_single_sign_single_elRO(run = run,idle = idle, sweep_time = sweep_time, RO = RO, state = state, error_sign = error_sign,el_RO = el_RO)
    for v in range(1):
        QEC_data_dict[y_list[v]] = {}
        QEC_data_dict[y_err_list[v]] = {}
        QEC_data_dict[y_list[v]] = RO_C*(QEC_dict[str(error_sign)]['positive'][c_list[v]]-
                                                        QEC_dict[str(error_sign)]['negative'][c_list[v]])/2



        QEC_data_dict[y_err_list[v]] = RO_C*(QEC_dict[str(error_sign)]['positive'][u_list[v]]**2+
                                                        QEC_dict[str(error_sign)]['negative'][u_list[v]]**2)**0.5/2


    QEC_data_dict['x'] = QEC_dict[str(error_sign)]['positive']['x']
    return QEC_data_dict

def single_qubit_no_QEC_data_single_state_RO(run = 0,older_than = None,state = 'Z',Qubit = 1,sweep_time =False, load_set = False):
    
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
                QEC_dict[str(error_sign)][el_RO], folder = single_Qubit_QEC_create_data_dict_single_error_single_elRO(run = run,older_than = older_than, Qubit = Qubit,sweep_time=sweep_time, state = state, error_sign = error_sign, el_RO = el_RO)
            else:
                QEC_dict[str(error_sign)][el_RO]= load_single_Qubit_QEC_dataset_single_sign_single_elRO(run = run,Qubit = Qubit, state = state, error_sign = error_sign,el_RO = el_RO)
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


def single_qubit_avg_state(state = 'Z',run = 1):
    single_no_QEC_data_dict_Z = {}
    for Qubit in [1,2,3]:
        single_no_QEC_data_dict_Z[str(Qubit)]  =  single_qubit_no_QEC_data_single_state_RO(state = state,Qubit = Qubit, load_set = True,run = run)

    single_no_QEC_data_dict_Z['y'] = 1/3.*(single_no_QEC_data_dict_Z[str(1)]['y'] +single_no_QEC_data_dict_Z[str(2)]['y'] +single_no_QEC_data_dict_Z[str(3)]['y'] )
    single_no_QEC_data_dict_Z['y_err'] = 1/3.*(single_no_QEC_data_dict_Z[str(1)]['y_err']**2 +single_no_QEC_data_dict_Z[str(2)]['y_err']**2 +single_no_QEC_data_dict_Z[str(3)]['y_err']**2 )**0.5
    single_no_QEC_data_dict_Z['x'] = single_no_QEC_data_dict_Z[str(1)]['x']

    return single_no_QEC_data_dict_Z




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

def QEC_state_sum_all(state = 'Z', RO = 0,run_list_00 = [1,2,3],run_list_01 = [1,2,3],run_list_10 = [2],run_list_11 = [3]):

    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']
    p_list = ['p00','p01','p10','p11']
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
                    state_dict[run][yy] = 1/float(len(run_list))*state_dict_single[run][yy]
                    state_dict[run]['temp'+ y_err_list[k]] = state_dict_single[run][ y_err_list[k]]**2
                for jj,p in enumerate(p_list):
                    state_dict[run][p] =  1/float(len(run_list))*state_dict_single[run][p]

            else:
                for k, yy in enumerate(y_list):
                    state_dict[run][ yy] += 1/float(len(run_list))*state_dict_single[run][ yy]
                    state_dict[run]['temp'+  y_err_list[k]] += state_dict_single[run][ y_err_list[k]]**2

                for jj,p in enumerate(p_list):
                    state_dict[run][p] +=  1/float(len(run_list))*state_dict_single[run][p]
            for k, yy in enumerate(y_list):
                state_dict[run][ y_err_list[k]] = 1/float(len(run_list))*(state_dict[run]['temp'+ y_err_list[k]])**0.5



    state_dict['x'] = state_dict_single[run]['x']


        # sum over all
    for j, run in enumerate(['00','01','10','11']):
        state_dict_single = state_dict[run]

        if j ==0:
            summed_dict = {}

            for k, yy in enumerate(y_list):
                summed_dict[yy] = 1/4.*state_dict[run][yy]
                summed_dict['temp'+y_err_list[k]] = state_dict[run][y_err_list[k]]**2
            for jj,p in enumerate(p_list):
                    summed_dict[p] =  1/4.*state_dict[run][p]
        else:

            for k, yy in enumerate(y_list):
                summed_dict[yy] += 1/4.*state_dict[run][yy]
                summed_dict['temp'+ y_err_list[k]] += state_dict[run][y_err_list[k]]**2

            for jj,p in enumerate(p_list):
                    summed_dict[p] +=  1/4.*state_dict[run][p]
        for k, yy in enumerate(y_list):
            summed_dict[ y_err_list[k]] = 1/4.*(summed_dict['temp'+ y_err_list[k]])**0.5





    summed_dict['x'] = state_dict['x']

    return summed_dict


def QEC_state_sum_RO_ZmZ(state = 'Z'):

    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']
    p_list = ['p00','p01','p10','p11']
    state_dict_single = {}
    state_dict = {}

    for RO in [0,1,2]:
        
        state_dict_single = QEC_state_sum_all(state = state, RO = RO,run_list_00 = [1,2,3],run_list_01 = [1,2,3],run_list_10 = [2],run_list_11 = [3])
        if RO == 0:
            for k, yy in enumerate(y_list):
                state_dict[yy] = 1/float(3)*state_dict_single[yy]
                state_dict['temp'+ y_err_list[k]] = state_dict_single[ y_err_list[k]]**2
            for jj,p in enumerate(p_list):
                state_dict[p] =  1/float(3)*state_dict_single[p]

        else:
            for k, yy in enumerate(y_list):
                state_dict[ yy] += 1/float(3)*state_dict_single[ yy]
                state_dict['temp'+  y_err_list[k]] += state_dict_single[ y_err_list[k]]**2

            for jj,p in enumerate(p_list):
                state_dict[p] +=  1/float(3)*state_dict_single[p]

        for k, yy in enumerate(y_list):
            state_dict[ y_err_list[k]] = 1/float(3)*(state_dict['temp'+ y_err_list[k]])**0.5



    state_dict['x'] = state_dict_single['x']


    return state_dict

''' plot single QEC / no QEC lines '''

def QEC_plot_single_state_RO_saved_data(run = 1, no_error = '00',append_encode = False,append_single = False,state = 'Z',RO = 0, plot_separate = False,plot_guide = True,plot_no_correct = False):

    if no_error != '':
        dataset_dict_full = QEC_sum_data_single_state_RO(run = run, no_error = no_error,state = state,RO = RO)
    elif no_error == '':
        dataset_dict_full =  QEC_state_sum_all(state = state, RO = RO)
    QEC_data_dict  = dataset_dict_full

    folder  = r'D:\measuring\data\QEC_data\figs\Z_only'

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
    ax.errorbar(x,y,yerr=y_err,color = 'k',ls = '', marker = 'o', ms = 4)
    x_fit, y_fit, p_err = fit_QEC_curve(x,y)
    ax.plot(x_fit,y_fit, color = 'k', label = 'QEC, p_c = '+str(int(p_err*100)/100.) )
    if plot_guide == True:
        ax.plot(x_g,y_g,color = 'g' )

    if plot_no_correct == True:
        y_no_corr = undo_correction_single_state_RO(run = run, no_error = no_error,state = state,RO = RO)
        print y_no_corr
        ax.errorbar(x,y_no_corr,yerr = y_err, color = 'm')


    if append_encode == True:
        QEC_data_dict =  no_QEC_data_single_state_RO(state = state,RO = RO, load_set = True)
        QEC_idle_data_dict =  no_QEC_data_single_state_RO(idle = True,state = state,RO = RO, load_set = True)


        x = QEC_data_dict['x']
        y = QEC_data_dict['y']
        y_idle = QEC_idle_data_dict['y']

        y_err = QEC_data_dict['y_err']
        y_idle_err = QEC_idle_data_dict['y_err']

        x_fit, y_fit, p_err = fit_QEC_curve(x,y)
        ax.plot(x_fit,y_fit, color = 'b', label = 'Encode, p_c = '+str(int(p_err*100)/100.) )
        ax.errorbar(x,y,yerr=y_err,color = 'b' ,ls = '', marker = 'o', ms = 4)
        # ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'b',ls = '--', label = 'idling for QEC time' )

    if append_single == True:
        QEC_data_dict =  single_qubit_no_QEC_data_single_state_RO(state = state,Qubit = RO+1, load_set = True)

        x = QEC_data_dict['x']
        y = QEC_data_dict['y']


        y_err = QEC_data_dict['y_err']
        x_fit, y_fit, p_err = fit_QEC_curve(x,y)
        ax.plot(x_fit,y_fit ,color = 'g', label = 'single qubit, p_c = '+str(int(p_err*100)/100.))
        ax.errorbar(x,y,yerr=y_err,color = 'g',ls = '', marker = 'o', ms = 4)

    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_QEC')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')
    ax.legend()
    if plot_no_correct == True:
        try:
            fig.savefig(
                os.path.join(folder,'undo_correct_error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_all'+'.png'))
        except:
            print 'Figure has not been saved.'
    else:
        try:
            fig.savefig(
                os.path.join(folder,'error_syn_'+no_error+'_run_'+str(run)+'_state_'+state+'_RO_'+str(RO)+'_all'+'.pdf'))
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
    y_list = ['y','y_00','y_01','y_10','y_11','y_no_corr']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11','y_err']
    dataset_dict_full= {}

    for i, syndrome in enumerate(run_list):

        dataset_dict_full[syndrome] = QEC_sum_data_single_state_RO_single_error_sign(no_error = syndrome, run =1,state = state,RO = RO,load_set = True,sweep_time = True)
        dataset_dict_full[syndrome]['y_no_corr'] = undo_correction_single_state_RO_error_sign(run = 2, no_error = syndrome,state = state,RO = RO, error_sign = 1,sweep_time=True)

        if i ==0:
            summed_dict = {}
            for k, yy in enumerate(y_list):
                summed_dict[yy] = 1/float(len(run_list))*dataset_dict_full[syndrome][yy]
                summed_dict['temp'+ y_err_list[k]] = dataset_dict_full[syndrome][y_err_list[k]]**2
            for jj, p in enumerate(p_list):
                    summed_dict[p_list[jj]]=1/float(len(run_list))*dataset_dict_full[syndrome][p]

        else:
            for k, yy in enumerate(y_list):
                summed_dict[yy] += 1/float(len(run_list))*dataset_dict_full[syndrome][yy]
                summed_dict['temp'+ y_err_list[k]] += dataset_dict_full[syndrome][y_err_list[k]]**2
            for jj, p in enumerate(p_list):
                    summed_dict[p_list[jj]]+=1/float(len(run_list))*dataset_dict_full[syndrome][p]

    for k, yy in enumerate(y_list):
        summed_dict[y_err_list[k]] = 1/float(len(run_list))*summed_dict['temp'+ y_err_list[k]]**0.5

    summed_dict['x'] = dataset_dict_full[syndrome]['x']

    return summed_dict

def QEC_plot_single_state_sweep_time(older_than = '20150107_090000',run = 1, no_error = '',no_error_list = [],state = 'Z',add_encode = False, add_single =False, plot_no_correct = False,load_set = False):
    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()
    color = ['r','g','b']
    dataset_dict_full = {}
    no_QEC_data_dict = {}
    QEC_single_data_dict = {}

    for RO in [0]:
        if no_error_list != []:
            dataset_dict_full[RO] = QEC_sweep_time_sum_error_syns(state = state,RO = RO,run_list = no_error_list)
        else:
            dataset_dict_full[RO] = QEC_sum_data_single_state_RO_single_error_sign(no_error = no_error,state = state,RO = RO,load_set = load_set, older_than = older_than,sweep_time = True, run = run)
        if add_single == True:
            no_QEC_data_dict[RO] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = state,RO = RO, load_set = load_set,error_sign = 0)

    for RO in [0]:
        fig1,ax1 = plt.subplots()
        fig2,ax2 = plt.subplots()
        fig3,ax3 = plt.subplots()

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
            if no_error_list == []:
                y_no_corr = undo_correction_single_state_RO_error_sign(run = 1, no_error = no_error,state = state,RO = RO,error_sign = 1,sweep_time=True)
                dataset_dict_full[RO]['y_no_corr'] = y_no_corr
            else:
                y_no_corr = dataset_dict_full[RO]['y_no_corr']
                # ax1.errorbar(x,y_no_corr,yerr = y_err, color = color[RO],ls = '-.', label = 'undo correction, Q'+str(RO+1))



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
            # no_QEC_data_dict[RO] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = state,RO = RO, load_set = load_set,error_sign = 0)
            # if RO == 1:
            x = no_QEC_data_dict[RO]['x']
            y = no_QEC_data_dict[RO]['y']
            y_err = no_QEC_data_dict[RO]['y_err']
            # ax1.errorbar(x,y,yerr=y_err,color = color[RO],ls = '--', label = 'Encoding, sweep time, Q' + str(RO+1) )

        if add_single == True:
            QEC_single_data_dict[RO] =  single_qubit_no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,state = state,sweep_time = True, error_sign = -1, Qubit = RO+1, load_set = True)
            x_single = QEC_single_data_dict[RO]['x']
            y_single = QEC_single_data_dict[RO]['y']
            y_single_err = QEC_single_data_dict[RO]['y_err']

            # ax1.errorbar(x_single,y_single,yerr=y_single_err,color = color[RO],ls = ':', label = 'Single Qubit, sweep time, Q' + str(RO+1) )


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

        if add_toffoli == True:
            if no_error_list != []:
                dataset_dict_full[6] = QEC_sweep_time_sum_error_syns(state = state,RO = 6,run_list = no_error_list)
            else:
                dataset_dict_full[6] = QEC_sum_data_single_state_RO_single_error_sign(no_error = no_error,state = state,RO = 6,load_set = load_set, older_than = older_than,sweep_time = True)

            y_toff = 1/2.*(dataset_dict_full[0]['y']+dataset_dict_full[1]['y']+dataset_dict_full[2]['y']-dataset_dict_full[6]['y'])
            y_toff_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
            x = dataset_dict_full[6]['x']+ np.ones(len(dataset_dict_full[6]['x']))*parity_time
            ax1.errorbar(x,y_toff,yerr=y_toff_err,color = 'k', label = 'QEC+ majority vote' )
            if plot_no_correct == True:
                if no_error_list == []:
                    y_no_corr = undo_correction_single_state_RO_error_sign(run = 1, no_error = no_error,state = state,RO = 6,error_sign = 1,sweep_time=True)
                    dataset_dict_full[6]['y_no_corr'] = y_no_corr
                y_toff = 1/2.*(dataset_dict_full[0]['y_no_corr']+dataset_dict_full[1]['y_no_corr']+dataset_dict_full[2]['y_no_corr']-dataset_dict_full[6]['y_no_corr'])
                y_toff_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
                x = dataset_dict_full[6]['x']+ np.ones(len(dataset_dict_full[6]['x']))*parity_time
                ax1.errorbar(x,y_toff,yerr=y_toff_err,color = 'k',ls = '-.', label = 'undo QEC+ toffoli' )
            if add_single == True:
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
                    os.path.join(folder,'sweep_time_undo_correct_error_syn_'+no_error+'_state_'+state+'qubit'+str(RO+1)+'_all'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
            except:
                print 'Figure has not been saved.'
        else:
            try:
                fig1.savefig(
                    os.path.join(folder,'sweep_time_error_syn_'+no_error+'_state_'+state+'qubit'+str(RO+1)+'_all'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
            except:
                print 'Figure has not been saved.'

        try:
            fig2.savefig(
                os.path.join(folder,'sweep_time_error_syn_'+no_error+'_state_'+state+'qubit'+str(RO+1)+'_ps'+'.png'))
        except:
            print 'Figure has not been saved.'

        try:
            fig3.savefig(
                os.path.join(folder,'sweep_time_error_syn_'+no_error+'_state_'+state+'qubit'+str(RO+1)+'_probs'+'.png'))
        except:
            print 'Figure has not been saved.'

def QEC_plot_Z_mZ_sweep_time(run = 1,no_error_list = [],add_encode = False, add_single =False, plot_no_correct = False,add_toffoli = True,load_set = False):

    folder  = r'D:\measuring\data\QEC_data\figs\timesweep'
    parity_time = 2*(4.996e-6*34 +11.312e-6*48) +2*(13.616e-6*34+4.996e-6*34) + 2* 150e-6

    color = ['r','g','b']
    dataset_dict_full = {}
    no_QEC_data_dict = {}
    QEC_single_data_dict = {}

    fig1, ax1 = plt.subplots()

    for RO in [0,1,2,6]:
        print RO
        dataset_dict_full[RO] = {}
        no_QEC_data_dict[RO] = {}
        QEC_single_data_dict[RO] = {}
        for state in ['Z','mZ']:
                print RO
                print state
                dataset_dict_full[RO][state] = QEC_sweep_time_sum_error_syns(state = state,RO = RO,run_list = no_error_list)
                no_QEC_data_dict[RO][state] =  no_QEC_data_single_state_RO_single_error_sign(sweep_time = True,idle = False,state = state,RO = RO, load_set = True,error_sign = 0)
                if RO != 6:
                    QEC_single_data_dict[RO][state] =  single_qubit_no_QEC_data_single_state_RO_single_error_sign(state = state,sweep_time = True, error_sign = -1, Qubit = RO+1, load_set = True)
                    data_list = [dataset_dict_full,no_QEC_data_dict,QEC_single_data_dict]
        # average Z and mZ data
        dataset_dict_full[RO]['x'] = dataset_dict_full[RO]['Z']['x']
        dataset_dict_full[RO]['y'] = 1/2.*(dataset_dict_full[RO]['Z']['y']-dataset_dict_full[RO]['mZ']['y'])
        dataset_dict_full[RO]['y_no_corr'] = 1/2.*(dataset_dict_full[RO]['Z']['y_no_corr']-dataset_dict_full[RO]['mZ']['y_no_corr'])
        dataset_dict_full[RO]['y_err']= 1/2.*(dataset_dict_full[RO]['Z']['y_err']**2+dataset_dict_full[RO]['mZ']['y_err']**2)**0.5

        no_QEC_data_dict[RO]['x'] = no_QEC_data_dict[RO]['Z']['x']
        no_QEC_data_dict[RO]['y'] = 1/2.*(no_QEC_data_dict[RO]['Z']['y']-no_QEC_data_dict[RO]['mZ']['y'])
        no_QEC_data_dict[RO]['y_err']= 1/2.*(no_QEC_data_dict[RO]['Z']['y_err']**2+no_QEC_data_dict[RO]['mZ']['y_err']**2)**0.5

        if RO !=6:
            QEC_single_data_dict[RO]['x'] = QEC_single_data_dict[RO]['Z']['x']
            QEC_single_data_dict[RO]['y'] = 1/2.*(QEC_single_data_dict[RO]['Z']['y']-QEC_single_data_dict[RO]['mZ']['y'])
            QEC_single_data_dict[RO]['y_err']= 1/2.*(QEC_single_data_dict[RO]['Z']['y_err']**2+QEC_single_data_dict[RO]['mZ']['y_err']**2)**0.5

    x = dataset_dict_full[RO]['x']+ np.ones(len(dataset_dict_full[RO]['x']))*parity_time


    y_toff_QEC = 1/2.*(dataset_dict_full[0]['y']+dataset_dict_full[1]['y']+dataset_dict_full[2]['y']-dataset_dict_full[6]['y'])
    y_toff_QEC_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
    ax1.errorbar(x,y_toff_QEC,yerr=y_toff_QEC_err,color = 'r',lw = 2,marker = 'o', ms = 4, label = 'QEC + majority vote' )

    y_toff_parity = 1/2.*(dataset_dict_full[0]['y_no_corr']+dataset_dict_full[1]['y_no_corr']+dataset_dict_full[2]['y_no_corr']-dataset_dict_full[6]['y_no_corr'])
    y_toff_parity_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
    x = dataset_dict_full[6]['x']+ np.ones(len(dataset_dict_full[6]['x']))*parity_time
    ax1.errorbar(x,y_toff_parity,yerr=y_toff_parity_err,color = 'b',ls = '-',lw = 2,marker = 'o', ms = 4, label = 'undo QEC + majority vote' )

    y_toff_encode = 1/2.*(no_QEC_data_dict[0]['y']+no_QEC_data_dict[1]['y']+no_QEC_data_dict[2]['y']-no_QEC_data_dict[6]['y'])
    y_toff_encode_err = 1/2.*(no_QEC_data_dict[0]['y_err']**2+no_QEC_data_dict[1]['y_err']**2+no_QEC_data_dict[2]['y_err']**2+no_QEC_data_dict[6]['y_err']**2)**0.5
    x_enc = no_QEC_data_dict[0]['x']

    ax1.errorbar(x_enc,y_toff_encode,yerr=y_toff_encode_err,color = 'k',ls ='-',marker = 'o', ms = 4, label = 'Encoding + majority vote' )



    # add best single qubit
    x_single = QEC_single_data_dict[1]['x']
    y_single = QEC_single_data_dict[1]['y']
    y_single_err = QEC_single_data_dict[1]['y_err']
    ax1.errorbar(x_single,y_single,yerr=y_single_err,color = 'g',ls = ':',marker = 'o', ms = 4,label = 'Single Qubit  ' + str(2) )



    ax1.hlines([0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    # ax1.vlines([x_enc[4],x[7]],-0.1,1.1,color = '0.5')
    # plt.axvspan(x_enc[4],x[7], facecolor='k', alpha=0.1)
    ax1.set_ylim(-0.1,1.1)
    ax1.set_xlim(-1e-3,35e-3)
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('Contrast')
    lgd = ax1.legend()#loc = 2, bbox_to_anchor = (1,1))

    try:
        fig1.savefig(
            os.path.join(folder,'QEC_sweep_time_syndromes'+str(no_error_list)+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
        fig1.savefig(
            os.path.join(folder,'QEC_sweep_time_syndromes'+str(no_error_list)+'.pdf'),bbox_extra_artists = (lgd,),bbox_inches='tight')
    except:
        print 'Figure has not been saved.'

def QEC_plot_Z_mZ_sweep_time_single_RO(run = 1,no_error_list =['11'],load_set = True):

    folder  = r'D:\measuring\data\QEC_data\figs\timesweep'
    parity_time = 2*(4.996e-6*34 +11.312e-6*48) +2*(13.616e-6*34+4.996e-6*34) + 2* 150e-6

    color = ['r','g','b','r','g','b','k']
    dataset_dict_full = {}
    no_QEC_data_dict = {}
    QEC_single_data_dict = {}

    fig1, ax1 = plt.subplots()

    for RO in [0,1,2,6]:
        print RO
        dataset_dict_full[RO] = {}
        no_QEC_data_dict[RO] = {}
        QEC_single_data_dict[RO] = {}
        for state in ['Z','mZ']:
                print RO
                print state
                dataset_dict_full[RO][state] = QEC_sweep_time_sum_error_syns(state = state,RO = RO,run_list = no_error_list)

        # average Z and mZ data
        dataset_dict_full[RO]['x'] = dataset_dict_full[RO]['Z']['x']
        dataset_dict_full[RO]['y'] = 1/2.*(dataset_dict_full[RO]['Z']['y']-dataset_dict_full[RO]['mZ']['y'])
        dataset_dict_full[RO]['y_no_corr'] = 1/2.*(dataset_dict_full[RO]['Z']['y_no_corr']-dataset_dict_full[RO]['mZ']['y_no_corr'])
        dataset_dict_full[RO]['y_err']= 1/2.*(dataset_dict_full[RO]['Z']['y_err']**2+dataset_dict_full[RO]['mZ']['y_err']**2)**0.5


        x = dataset_dict_full[RO]['x']+ np.ones(len(dataset_dict_full[RO]['x']))*parity_time

        ax1.errorbar(x,dataset_dict_full[RO]['y_no_corr'],yerr=dataset_dict_full[RO]['y_err'],color = color[RO],ls ='-',marker = 'o', ms = 4, label = 'RO '+str(RO) )




    ax1.hlines([0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    # ax1.vlines([x_enc[4],x[7]],-0.1,1.1,color = '0.5')
    # plt.axvspan(x_enc[4],x[7], facecolor='k', alpha=0.1)
    ax1.set_ylim(-0.1,1.1)
    ax1.set_xlim(-1e-3,35e-3)
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('Contrast')
    lgd = ax1.legend()#loc = 2, bbox_to_anchor = (1,1))

    try:
        fig1.savefig(
            os.path.join(folder,'QEC_sweep_time_single_RO'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
        fig1.savefig(
            os.path.join(folder,'QEC_sweep_time_single_RO'+'.pdf'),bbox_extra_artists = (lgd,),bbox_inches='tight')
    except:
        print 'Figure has not been saved.'

def QEC_plot_Z_mZ_sweep_time_compare_syndromes(no_error_list = ['00','01','10','11'],add_encode = False,encode_run = 0, add_single =False, plot_no_correct = True,add_toffoli = True,load_set = True):

    folder  = r'D:\measuring\data\QEC_data\figs\timesweep'
    parity_time = 2*(4.996e-6*34 +11.312e-6*48) +2*(13.616e-6*34+4.996e-6*34) + 2* 150e-6

    color = ['c','k','m','b']
    dataset_dict_full = {}
    no_QEC_data_dict = {}
    QEC_single_data_dict = {}

    fig1, ax1 = plt.subplots(figsize = (10,8))

    for i, syndrome in enumerate(no_error_list):

        for RO in [0,1,2,6]:
            print RO
            dataset_dict_full[RO] = {}
            no_QEC_data_dict[RO] = {}
            QEC_single_data_dict[RO] = {}
            for state in ['Z','mZ']:
                    dataset_dict_full[RO][state] = QEC_sweep_time_sum_error_syns(state = state,RO = RO,run_list = [syndrome])
                    no_QEC_data_dict[RO][state] =  no_QEC_data_single_state_RO_single_error_sign(run = encode_run,sweep_time = True,idle = False,state = state,RO = RO, load_set = True,error_sign = 0)
                    if RO != 6:
                        QEC_single_data_dict[RO][state] =  single_qubit_no_QEC_data_single_state_RO_single_error_sign(state = state,sweep_time = True, error_sign = -1, Qubit = RO+1, load_set = True)
                        data_list = [dataset_dict_full,no_QEC_data_dict,QEC_single_data_dict]
            # average Z and mZ data
            dataset_dict_full[RO]['x'] = dataset_dict_full[RO]['Z']['x']
            dataset_dict_full[RO]['y'] = 1/2.*(dataset_dict_full[RO]['Z']['y']-dataset_dict_full[RO]['mZ']['y'])
            dataset_dict_full[RO]['y_no_corr'] = 1/2.*(dataset_dict_full[RO]['Z']['y_no_corr']-dataset_dict_full[RO]['mZ']['y_no_corr'])
            dataset_dict_full[RO]['y_err']= 1/2.*(dataset_dict_full[RO]['Z']['y_err']**2+dataset_dict_full[RO]['mZ']['y_err']**2)**0.5

            no_QEC_data_dict[RO]['x'] = no_QEC_data_dict[RO]['Z']['x']
            no_QEC_data_dict[RO]['y'] = 1/2.*(no_QEC_data_dict[RO]['Z']['y']-no_QEC_data_dict[RO]['mZ']['y'])
            no_QEC_data_dict[RO]['y_err']= 1/2.*(no_QEC_data_dict[RO]['Z']['y_err']**2+no_QEC_data_dict[RO]['mZ']['y_err']**2)**0.5

            if RO !=6:
                QEC_single_data_dict[RO]['x'] = QEC_single_data_dict[RO]['Z']['x']
                QEC_single_data_dict[RO]['y'] = 1/2.*(QEC_single_data_dict[RO]['Z']['y']-QEC_single_data_dict[RO]['mZ']['y'])
                QEC_single_data_dict[RO]['y_err']= 1/2.*(QEC_single_data_dict[RO]['Z']['y_err']**2+QEC_single_data_dict[RO]['mZ']['y_err']**2)**0.5

        x = dataset_dict_full[RO]['x']+ np.ones(len(dataset_dict_full[RO]['x']))*parity_time


        y_toff_QEC = 1/2.*(dataset_dict_full[0]['y']+dataset_dict_full[1]['y']+dataset_dict_full[2]['y']-dataset_dict_full[6]['y'])
        y_toff_QEC_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
        ax1.errorbar(x,y_toff_QEC,yerr=y_toff_QEC_err,color = color[i],lw = 2,marker = 'o', ms = 4, label = 'QEC + majority vote '+ syndrome )

        y_toff_parity = 1/2.*(dataset_dict_full[0]['y_no_corr']+dataset_dict_full[1]['y_no_corr']+dataset_dict_full[2]['y_no_corr']-dataset_dict_full[6]['y_no_corr'])
        y_toff_parity_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
        x = dataset_dict_full[6]['x']+ np.ones(len(dataset_dict_full[6]['x']))*parity_time
        ax1.errorbar(x,y_toff_parity,yerr=y_toff_parity_err,color =color[i],ls = ':',lw = 2,marker = 'o', ms = 4, label = 'undo QEC + majority vote '+ syndrome )

    # y_toff_encode = 1/2.*(no_QEC_data_dict[0]['y']+no_QEC_data_dict[1]['y']+no_QEC_data_dict[2]['y']-no_QEC_data_dict[6]['y'])
    # y_toff_encode_err = 1/2.*(no_QEC_data_dict[0]['y_err']**2+no_QEC_data_dict[1]['y_err']**2+no_QEC_data_dict[2]['y_err']**2+no_QEC_data_dict[6]['y_err']**2)**0.5
    # x_enc = no_QEC_data_dict[0]['x']
    # ax1.errorbar(x_enc,y_toff_encode,yerr=y_toff_encode_err,color = 'k',ls ='-',marker = 'o', ms = 4, label = 'Encoding + majority vote' )

    # # add best single qubit
    # x_single = QEC_single_data_dict[1]['x']
    # y_single = QEC_single_data_dict[1]['y']
    # y_single_err = QEC_single_data_dict[1]['y_err']
    # ax1.errorbar(x_single,y_single,yerr=y_single_err,color = 'g',ls = '-',marker = 'o', ms = 4, label = 'Single Qubit  ' + str(2) )
    # ax1.errorbar(x,dataset_dict_full[1]['y'],yerr=dataset_dict_full[1]['y_err'],color = 'g',ls ='-',marker = '*', ms = 4, label = 'QEC, decode to Qubit 2' )


    # y_no_corr = dataset_dict_full[1]['y_no_corr']
    # ax1.errorbar(x,y_no_corr,yerr = dataset_dict_full[1]['y_err'], color = 'g',ls = '-.', label = 'undo correction, Q'+str(RO+1))

    ax1.set_title('QEC sweep time compare syndromes')
    ax1.hlines([0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    # ax1.vlines([x_enc[4],x[7]],-0.1,1.1,color = '0.5')
    # plt.axvspan(x_enc[4],x[7], facecolor='k', alpha=0.1)
    ax1.set_ylim(-0.1,1.1)
    ax1.set_xlim(-1e-3,35e-3)
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('Contrast')
    lgd = ax1.legend()#loc = 2, bbox_to_anchor = (1,1))

    try:
        fig1.savefig(
            os.path.join(folder,'QEC_sweep_time_compare_syndromes.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
        fig1.savefig(
            os.path.join(folder,'QEC_sweep_time_compare_syndromes.pdf'),bbox_extra_artists = (lgd,),bbox_inches='tight')
    except:
        print 'Figure has not been saved.'

    # x = QEC_data_dict['x']+ np.ones(len(QEC_data_dict['x']))*parity_time

    # ax1.errorbar(x,y,yerr=y_err,color = color[RO], label = 'QEC, decode to Qubit '+str(RO+1))

    # if plot_no_correct == True:
    #     y_no_corr = dataset_dict_full[RO]['y_no_corr']
    #     ax1.errorbar(x,y_no_corr,yerr = y_err, color = color[RO],ls = '-.', label = 'undo correction, Q'+str(RO+1))

    # if add_encode == True:
    #     x = no_QEC_data_dict[RO]['x']
    #     y = no_QEC_data_dict[RO]['y']
    #     y_err = no_QEC_data_dict[RO]['y_err']
    #     ax1.errorbar(x,y,yerr=y_err,color = color[RO],ls = '--', label = 'Encoding, sweep time, Q' + str(RO+1) )

    # if add_single == True:
    #     x_single = QEC_single_data_dict[RO]['x']
    #     y_single = QEC_single_data_dict[RO]['y']
    #     y_single_err = QEC_single_data_dict[RO]['y_err']
    #     ax1.errorbar(x_single,y_single,yerr=y_single_err,color = color[RO],ls = ':', label = 'Single Qubit, sweep time, Q' + str(RO+1) )


    # if no_error_list != []:
    #     dataset_dict_full[6] = QEC_sweep_time_sum_error_syns(state = state,RO = 6,run_list = no_error_list)
    # else:
    #     dataset_dict_full[6] = QEC_sum_data_single_state_RO_single_error_sign(no_error = no_error,state = state,RO = 6,load_set = True, older_than = older_than,sweep_time = True)

    # if add_toffoli == True:
    #     y_toff = 1/2.*(dataset_dict_full[0]['y']+dataset_dict_full[1]['y']+dataset_dict_full[2]['y']-dataset_dict_full[6]['y'])
    #     y_toff_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
    #     x = dataset_dict_full[6]['x']+ np.ones(len(dataset_dict_full[6]['x']))*parity_time
    #     ax1.errorbar(x,y_toff,yerr=y_toff_err,color = 'k', label = 'QEC+ majority vote' )
    #     if plot_no_correct == True:
    #         if no_error_list == []:
    #             y_no_corr = undo_correction_single_state_RO_error_sign(run = 1, no_error = no_error,state = state,RO = 6,error_sign = 1,sweep_time=True)
    #             dataset_dict_full[6]['y_no_corr'] = y_no_corr
    #         y_toff = 1/2.*(dataset_dict_full[0]['y_no_corr']+dataset_dict_full[1]['y_no_corr']+dataset_dict_full[2]['y_no_corr']-dataset_dict_full[6]['y_no_corr'])
    #         y_toff_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
    #         x = dataset_dict_full[6]['x']+ np.ones(len(dataset_dict_full[6]['x']))*parity_time
    #         ax1.errorbar(x,y_toff,yerr=y_toff_err,color = 'k',ls = '-.', label = 'undo QEC+ toffoli' )
    #     if add_single == True:
    #         no_QEC_data_dict[6] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = state,RO = 6, load_set = load_set,error_sign = 0)

    #         y_toff = 1/2.*(no_QEC_data_dict[0]['y']+no_QEC_data_dict[1]['y']+no_QEC_data_dict[2]['y']-no_QEC_data_dict[6]['y'])
    #         y_toff_err = 1/2.*(no_QEC_data_dict[0]['y_err']**2+no_QEC_data_dict[1]['y_err']**2+no_QEC_data_dict[2]['y_err']**2+no_QEC_data_dict[6]['y_err']**2)**0.5
    #         x = no_QEC_data_dict[0]['x']
    #         ax1.errorbar(x,y_toff,yerr=y_toff_err,color = 'k',ls ='--', label = 'Encoding, sweep time, majority vote' )

# def no_QEC_plot_single_state_RO(state = 'Z',RO = 0, load_set = False, older_than = None):

#     QEC_data_dict =  no_QEC_data_single_state_RO(older_than = older_than,state = state,RO = RO, load_set = load_set)
#     QEC_idle_data_dict =  no_QEC_data_single_state_RO(idle = True,older_than = older_than,state = state,RO = RO, load_set = load_set)




    # if plot_no_correct == True:
    #     try:
    #         fig1.savefig(
    #             os.path.join(folder,'sweep_time_undo_correct_error_syn_'+no_error+'_state_'+state+'qubit'+str(RO+1)+'_all'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
    #     except:
    #         print 'Figure has not been saved.'
    # else:
    #     try:
    #         fig1.savefig(
    #             os.path.join(folder,'sweep_time_error_syn_'+no_error+'_state_'+state+'qubit'+str(RO+1)+'_all'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
    #     except:
    #         print 'Figure has not been saved.'

    # try:
    #     fig2.savefig(
    #         os.path.join(folder,'sweep_time_error_syn_'+no_error+'_state_'+state+'qubit'+str(RO+1)+'_ps'+'.png'))
    # except:
    #     print 'Figure has not been saved.'

    # try:
    #     fig3.savefig(
    #         os.path.join(folder,'sweep_time_error_syn_'+no_error+'_state_'+state+'qubit'+str(RO+1)+'_probs'+'.png'))
    # except:
    #     print 'Figure has not been saved.'

def no_QEC_plot_single_state_RO(state = 'Z',RO = 0, load_set = False, older_than = None, run = 0):

    print state
    print RO
    QEC_data_dict =  no_QEC_data_single_state_RO(older_than = older_than,state = state,RO = RO, load_set = load_set, run = run)
    # QEC_idle_data_dict =  no_QEC_data_single_state_RO(idle = True,older_than = older_than,state = state,RO = RO, load_set = load_set)


    folder  = r'D:\measuring\data\QEC_data\figs\Encoding'



    x = QEC_data_dict['x']
    y = QEC_data_dict['y']
    # y_idle = QEC_idle_data_dict['y']


    y_err = QEC_data_dict['y_err']
    # y_idle_err = QEC_idle_data_dict['y_err']

    fig,ax = plt.subplots()
    ax.errorbar(x,y,yerr=y_err,color = 'k', label = 'no idling' )
    # ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'b', label = 'idling for QEC time' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('state_'+state+'_RO_'+str(RO)+'_noQEC')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')
    ax.legend()
    try:
        fig.savefig(
            os.path.join(folder,'run'+str(run)+'no_QEC'+'_state_'+state+'_RO_'+str(RO)+'.png'))
    except:
        print 'Figure has not been saved.'

def no_QEC_sweep_time_plot_single_state(state = 'Z',load_set = False, older_than = None):
    fig,ax = plt.subplots()
    color = ['r','g','b','k','k','k','k']
    QEC_data_dict = {}
    QEC_single_data_dict = {}

    for  RO in [0,1,2,6]:
        print RO
        QEC_data_dict[RO] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = state,RO = RO, load_set = load_set,error_sign = 0)
        # QEC_single_data_dict[RO] =  single_qubit_no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,state = state,sweep_time = True, error_sign = -1, Qubit = RO+1, load_set = True)

        folder  = r'D:\measuring\data\QEC_data\figs\Encoding'


        x = QEC_data_dict[RO]['x']
        y = QEC_data_dict[RO]['y']
        y_err = QEC_data_dict[RO]['y_err']


        # x_single = QEC_single_data_dict[RO]['x']
        # y_single = QEC_single_data_dict[RO]['y']
        # y_single_err = QEC_single_data_dict[RO]['y_err']

        if RO == 2:
            extra_time = 2*(4.996e-6*34 +11.312e-6*48) # 2*(13.616e-6*34)+116e-6
        elif RO == 0:
             extra_time = 2*(13.616e-6*34+11.312e-6*48)# 2*(4.996e-6*34)+116e-6
        elif RO == 1:
            extra_time =  2*(4.996e-6*34 +13.616e-6*34)#2*(11.312e-6*48)+116e-6

        # x_single = x_single #+ np.ones(len(x_single))* extra_time

        ax.errorbar(x,y,yerr=y_err,color = color[RO], label = 'Encoding, sweep time, RO ' + str(RO) )
        # ax.errorbar(x_single,y_single,yerr=y_single_err,color = color[RO],ls = ':', label = 'Single Qubit, sweep time, Q' + str(RO+1) )

    # QEC_data_dict[6] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = state,RO = 6, load_set = load_set,error_sign = 0)

    # y_toff = 1/2.*(QEC_data_dict[0]['y']+QEC_data_dict[1]['y']+QEC_data_dict[2]['y']-QEC_data_dict[6]['y'])
    # y_toff_err = 1/2.*(QEC_data_dict[0]['y_err']**2+QEC_data_dict[1]['y_err']**2+QEC_data_dict[2]['y_err']**2+QEC_data_dict[6]['y_err']**2)**0.5

    # ax.errorbar(x,y_toff,yerr=y_toff_err,color = 'k', label = 'Encoding, sweep time, toff' )
    ax.set_ylim(-.1,1.1)
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

def no_QEC_sweep_time_plot_single_RO(load_set = False, older_than = None):
    fig, ax = plt.subplots()
    color = ['r','b','g','','','','k']
    QEC_data_dict = {}
    y_temp = {}
    for RO in [0,1,2,6]:
        QEC_data_dict['Z'] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = 'Z',RO = RO, load_set = load_set,error_sign = 0)
        QEC_data_dict['mZ'] =  no_QEC_data_single_state_RO_single_error_sign(older_than = older_than,sweep_time = True,idle = False,state = 'mZ',RO = RO, load_set = load_set,error_sign = 0)

    # y_toff = 1/2.*(QEC_data_dict[0]['y']+QEC_data_dict[1]['y']+QEC_data_dict[2]['y']-QEC_data_dict[6]['y'])
    # y_toff_err = 1/2.*(QEC_data_dict[0]['y_err']**2+QEC_data_dict[1]['y_err']**2+QEC_data_dict[2]['y_err']**2+QEC_data_dict[6]['y_err']**2)**0.5

        x = QEC_data_dict['Z']['x']
        y = 1/2.*(QEC_data_dict['Z']['y']-QEC_data_dict['mZ']['y'])
        y_err = 1/2.*(QEC_data_dict['Z']['y_err']**2+QEC_data_dict['mZ']['y_err']**2)**0.5

        folder  = r'D:\measuring\data\QEC_data\figs\Encoding'

      
        x = x*1e3
        ax.errorbar(x,y,yerr=y_err, label = 'Encoding, RO '+str(RO), color = color[RO],ls = '')
        x_temp, y_temp[str(RO)],T, T_err = fit_timesweep_single(x,y,return_errorbar = True)
        ax.plot(x_temp,y_temp[str(RO)],color = color[RO],ls = '-',lw = 1)

    ax.set_ylim(-0.1,1.1)
    # ax.set_xlim(-0.001,0.035)
    # ax.set_title('state_'+state+'_RO_'+str(RO)+'_noQEC_sweep_time')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('time (ms)')
    ax.set_ylabel('Contrast')
    # fitfunc_dict = {}
    # fitfunc_dict = pickle.load( open( "ramseys.p", "rb" ) )

    # ax.plot(x_temp,(y_temp['0']*y_temp['1']*y_temp['2'])*y_temp[str(6)][0]/(y_temp['0'][0]*y_temp['1'][0]*y_temp['2'][0]), color = 'k', lw = 2, ls = ':', label = 'multiplied, not shifted')
    # y_temp['0'] = y_temp['0'][0:-(400+1090)]
    # y_temp['1'] = y_temp['1'][(400):(-1090)]
    # y_temp['2'] = y_temp['2'][(400+1090):]
    # x_temp = x_temp[0:-(400+1090)]
    y_temp['0'] = y_temp['0'][(1090+980):-400]
    y_temp['1'] = y_temp['1'][(400+980):-1090]
    y_temp['2'] = y_temp['2'][(400+1090):-980]
    x_temp = x_temp[0:-(400+1090+980)]
    print size(y_temp['0'])
    print size(x_temp)
    
    ax.plot(x_temp,(y_temp['0']*y_temp['1']*y_temp['2'])*y_temp[str(6)][0]/(y_temp['0'][0]*y_temp['1'][0]*y_temp['2'][0]),
     color = 'k', lw = 2, ls = '--', label = 'multiplied, shifted')

    # x_temp, y_temp,T, T_err = fit_timesweep_single(x_temp,y_temp['0']*y_temp['1']*y_temp['2'],return_errorbar = True)
    # ax.plot(fitfunc_dict['x']*1e3,fitfunc_dict['y_decRamsey_C1_ms0'], color ='r')
    # ax.plot(fitfunc_dict['x']*1e3,fitfunc_dict['y_decRamsey_C2_ms0'], color ='b')
    # ax.plot(fitfunc_dict['x']*1e3,fitfunc_dict['y_decRamsey_C5_ms0'], color ='g')
    # ax.plot(fitfunc_dict['x']*1e3,fitfunc_dict['y_decRamsey_C1_ms0']*fitfunc_dict['y_decRamsey_C2_ms0']*fitfunc_dict['y_decRamsey_C5_ms0'], color ='k')

    # if state == 'Z':
    #     location = 4
    # elif state == 'mZ':
    location =1


    legend = ax.legend(loc = location)

    ax.set_title('Compare_RO')
    for label in legend.get_texts():
        label.set_fontsize('x-small')

    try:
        fig.savefig(
            os.path.join(folder,'Compare_RO_shifted_RO_at_end.png'))
    except:
        print 'Figure has not been saved.'

def single_Qubit_no_QEC_plot_single_state_RO(run_list = [0],state = 'Z',Qubit = 1, load_set = False, older_than = None):
    fig,ax = plt.subplots()

    for run in run_list:
        QEC_data_dict =  single_qubit_no_QEC_data_single_state_RO(run = run,older_than = older_than,state = state,Qubit = Qubit, load_set = load_set)


        folder  = r'D:\measuring\data\QEC_data\figs\Encoding'


        x = QEC_data_dict['x']
        y = QEC_data_dict['y']


        y_err = QEC_data_dict['y_err']

        ax.errorbar(x,y,yerr=y_err, label = 'run '+str(run) )
    ax.legend()
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('state_'+state+'Qubit_'+str(Qubit)+'_single_qubit')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')
    try:
        fig.savefig(
            os.path.join(folder,'single_qubit'+'_state_'+state+'_Qubit_'+str(Qubit)+'_runs_'+str(run_list)+'.png'))
    except:
        print 'Figure has not been saved.'

def single_Qubit_sweep_time_no_QEC_plot_single_state(state = 'Z', load_set = False, older_than = None):
    fig,ax = plt.subplots()
    color = ['r','g','b']
    for Qubit in [2]:#[1,2,3]:

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

    if state == 'X' or state == 'mX':
        RO_list = [6]
    if state == 'Y' or state == 'mY':
        RO_list = [4,5,6]
    if state == 'Z' or state == 'mZ':
        RO_list = [0,1,2]
    folder  = r'D:\measuring\data\QEC_data\figs'

    for ii,RO in enumerate(RO_list):
        Qubit = RO+1

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

        # no_QEC_idle_data_dict =  no_QEC_data_single_state_RO(idle = True,state = state,RO = RO, load_set = True)
        # y = no_QEC_idle_data_dict['y']
        # y_err = no_QEC_idle_data_dict['y_err']
        # x = no_QEC_idle_data_dict['x']
        # ax.errorbar(x,sign*y,yerr=y_err,color = 'k',ls =':', label =  'Encode+idle')

        # toff_process_dict_idle = no_QEC_toffoli_fids(idle = True)
        # y = toff_process_dict_idle['toff_'+state+'y']
        # y_err = toff_process_dict_idle['toff_'+state+'y_err']
        # x = toff_process_dict_idle['x']
        # ax.errorbar(x,y,yerr=y_err,color = 'k', ls = '--', label =  'Encode toffoli +idle')

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

def plot_single_state_QEC_vs_noQEC_11_vs_idle(older_than = '20150206070800',state = 'Z',RO = '0', load_set = True,run_QEC =4, run_enc = 1):
    
    dataset_11_QEC   = QEC_sum_data_single_state_RO(no_error = '11',state = state,RO = RO,load_set = load_set, older_than = older_than,do_p = True,run = run_QEC)
    dataset_idle_ENC = no_QEC_data_single_state_RO(older_than = older_than,state = state,RO = RO, idle = True,load_set = load_set,do_p = True,run = run_enc)

    folder  = r'D:\measuring\data\QEC_data\figs\11vsIdle'


    fig,ax = plt.subplots()

    x = dataset_11_QEC['x']
    y = dataset_11_QEC['y']

    y_err = dataset_11_QEC['y_err']
    
    ax.errorbar(x,y,yerr=y_err,color = 'r', label = '11 QEC' )
    x = dataset_idle_ENC['x']
    y = dataset_idle_ENC['y']

    y_err = dataset_idle_ENC['y_err']
    
    ax.errorbar(x,y,yerr=y_err,color = 'k', label = 'idle' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title('state_'+state+'_RO_'+str(RO)+ '11 vs Idle')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')
    ax.legend()
    try:
        fig.savefig(
            os.path.join(folder,'11_QECvsIdle_runs_'+str([run_QEC,run_enc])+'_state_'+state+'_RO_'+str(RO)+'.png'))
    except:
        print 'Figure has not been saved.'



def plot_single_state_QEC_vs_noQEC_11_vs_idle_compare_runs(state = 'Z',RO = '0', load_set = True,run_QEC_list =[3,5,6,7], run_enc_list = [0,2,3,4]):
    
    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    color_list = ['r','g','b','k']

    for ii in range(len(run_QEC_list)):
        run_QEC = run_QEC_list[ii]
        run_enc = run_enc_list[ii]
        dataset_11_QEC   = QEC_sum_data_single_state_RO(no_error = '11',state = state,RO = RO,load_set = load_set, do_p = True,run = run_QEC)
        dataset_idle_ENC = no_QEC_data_single_state_RO(state = state,RO = RO, idle = True,load_set = load_set,do_p = True,run = run_enc)

        folder  = r'D:\measuring\data\QEC_data\figs\11vsIdle'



        x = dataset_11_QEC['x']
        y = dataset_11_QEC['y']

        y_err = dataset_11_QEC['y_err']
        
        ax1.errorbar(x,y,yerr=y_err,color = color_list[ii], label = 'run ' +str(run_QEC) )
        x = dataset_idle_ENC['x']
        y = dataset_idle_ENC['y']

        y_err = dataset_idle_ENC['y_err']
        
        ax2.errorbar(x,y,yerr=y_err,color = color_list[ii], label = 'run ' +str(run_enc))
    
    ax1.set_ylim(-1.1,1.1)
    ax1.set_xlim(-0.1,1.1)
    ax1.set_title('state_'+state+'_RO_'+str(RO)+ '11 QEC runs')
    ax1.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax1.set_xlabel('error probability')
    ax1.set_ylabel('Contrast')
    ax1.legend()
    ax2.set_ylim(-1.1,1.1)
    ax2.set_xlim(-0.1,1.1)
    ax2.set_title('state_'+state+'_RO_'+str(RO)+ 'Idle runs')
    ax2.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax2.set_xlabel('error probability')
    ax2.set_ylabel('Contrast')
    ax2.legend()
    try:
        fig1.savefig(
            os.path.join(folder,'QEC_compare_runs_state_'+state+'_RO_'+str(RO)+'.png'))
    except:
        print 'Figure has not been saved.'

    try:
        fig2.savefig(
            os.path.join(folder,'ENC_compare_runs_state_'+state+'_RO_'+str(RO)+'.png'))
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


    p_list = ['p00','p01','p10','p11']
    p_dict = {}
    for p in range(4):
            p_dict[p_list[p]] = {}

            p_dict[p_list[p]] = dataset_dict[state]
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

def QEC_sum_probs(run_list = [1] ,no_error = '00'):
    p_list = ['p00','p01','p10','p11']
    p_dict = {}

    for mm, run in enumerate(run_list):
        dataset_dict = {}
        for ii, state in enumerate(['Z','mZ','Y','mY', 'X','mX']):
            dataset_dict[state] = {}
            if state == 'X' or state == 'mX':
                RO_list = [6]
            if state == 'Y' or state == 'mY':
                RO_list = [4,5,6]
            if state == 'Z' or state == 'mZ':
                RO_list = [0,1,2]

            for jj,RO in enumerate(RO_list):
                dataset_dict[state]['Tomo_'+str(RO)] = {}

                dataset_dict[state]['Tomo_'+str(RO)] = QEC_sum_data_single_state_RO(run = run, no_error = no_error,state = state,RO = RO)



                for p in range(4):
                    if mm == 0 and ii == 0 and jj == 0:
                        p_dict[p_list[p]] = {}
                        p_dict[p_list[p]] = 1/(len(run_list)*14.)*dataset_dict[state]['Tomo_'+str(RO)][p_list[p]]
                    else:
                        p_dict[p_list[p]] += 1/(len(run_list)*14.)*dataset_dict[state]['Tomo_'+str(RO)][p_list[p]]

    p_dict['x'] = dataset_dict[state]['Tomo_'+str(RO)]['x']
    return p_dict

def QEC_timesweep_sum_probs(no_error = '11'):
    p_list = ['p00','p01','p10','p11']
    p_dict = {}


    dataset_dict = {}
    for ii, state in enumerate(['Z','mZ']):
        dataset_dict[state] = {}
        RO_list = [0,1,2]

        for jj,RO in enumerate(RO_list):
            dataset_dict[state]['Tomo_'+str(RO)] = {}

            dataset_dict[state]['Tomo_'+str(RO)] = QEC_sum_data_single_state_RO_single_error_sign(no_error = no_error, run =1,state = state,RO = RO,load_set = True,sweep_time = True)



            for p in range(4):
                if  ii == 0 and jj == 0:
                    p_dict[p_list[p]] = {}
                    # print dataset_dict[state]['Tomo_'+str(RO)][p_list[p]]
                    p_dict[p_list[p]] = 1/(6.)*dataset_dict[state]['Tomo_'+str(RO)][p_list[p]]
                else:
                    p_dict[p_list[p]] += 1/(6.)*dataset_dict[state]['Tomo_'+str(RO)][p_list[p]]
    p_dict['x'] = dataset_dict[state]['Tomo_'+str(RO)]['x']
    return p_dict

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

def QEC_process_fids_sum_all(syndrome_list = ['00','01','10','11'],run_list_00 = [1,2,3],run_list_01 = [1,2,3],run_list_10 = [2],run_list_11 = [3]):
    dec_list = ['dec_1_','dec_2_','dec_3_','dec_avg_']
    y_list = ['y','y_00','y_01','y_10','y_11','y_new']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11','y_err']

    for i, run in enumerate(syndrome_list):
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
                        process_dict[item+ yy] = 1/float(len(syndrome_list))*process_dict_single[item+ yy]
                        process_dict['temp'+item+ y_err_list[k]] = process_dict_single[item+ y_err_list[k]]**2
            else:
                for ii, item in enumerate(dec_list):
                    for k, yy in enumerate(y_list):
                        process_dict[item+ yy] += 1/float(len(syndrome_list))*process_dict_single[item+ yy]
                        process_dict['temp'+item+ y_err_list[k]] += process_dict_single[item+ y_err_list[k]]**2
            for ii, item in enumerate(dec_list):
                for k, yy in enumerate(y_list):
                    process_dict[item+ y_err_list[k]] = 1/float(len(syndrome_list))*(process_dict['temp'+item+ y_err_list[k]])**0.5

    process_dict['x'] = process_dict_single['x']

    return process_dict

def no_QEC_process_fids(idle = False, run = 1):

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

            dataset_dict[state]['Tomo_'+str(RO)] = no_QEC_data_single_state_RO(idle = idle,run = run,state = state,RO = RO,load_set = True)
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

    process_dict['dec_avg_'+'y_err'] = 1/3.*1/8.*((dataset_dict['Z']['Tomo_'+str(0)]['y_err']**2 + dataset_dict['mZ']['Tomo_'+str(0)]['y_err']**2
                    + dataset_dict['Y']['Tomo_'+str(5)]['y_err']**2 +  dataset_dict['mY']['Tomo_'+str(5)]['y_err']**2
                    + dataset_dict['X']['Tomo_'+str(6)]['y_err']**2 + dataset_dict['mX']['Tomo_'+str(6)]['y_err']**2 ) +

                                (dataset_dict['Z']['Tomo_'+str(1)]['y_err']**2 + dataset_dict['mZ']['Tomo_'+str(1)]['y_err']**2
                                        + dataset_dict['Y']['Tomo_'+str(6)]['y_err']**2 +  dataset_dict['mY']['Tomo_'+str(6)]['y_err']**2
                                        + dataset_dict['X']['Tomo_'+str(6)]['y_err']**2 + dataset_dict['mX']['Tomo_'+str(6)]['y_err']**2 )+

                               (dataset_dict['Z']['Tomo_'+str(2)]['y_err']**2 + dataset_dict['mZ']['Tomo_'+str(2)]['y_err']**2
                                        + dataset_dict['Y']['Tomo_'+str(4)]['y_err']**2 +  dataset_dict['mY']['Tomo_'+str(4)]['y_err']**2
                                        + dataset_dict['X']['Tomo_'+str(6)]['y_err']**2 + dataset_dict['mX']['Tomo_'+str(6)]['y_err']**2 ))**0.5


    process_dict['x'] = dataset_dict['Z']['Tomo_'+str(0)]['x']
    return process_dict

def no_QEC_process_fids_sum_runs(idle = False, run_list = [2,3,4]):

    dec_list = ['dec_1_','dec_2_','dec_3_','dec_avg_']
    y_list = ['y']
    y_err_list = ['y_err']


    for i in range(len(run_list)):
        process_dict_single = {}
        process_dict_single = no_QEC_process_fids(run = run_list[i], idle = idle)
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

def no_QEC_toffoli_fids(idle = False,run_list = [0],state_list=['Z','mZ','Y','mY', 'X','mX'],add_4 = False,only_4=False,do_weighted = False):

    dataset_dict = {}
    toff_dict = {}

    for state in state_list:
        dataset_dict[state] = {}
        if state == 'X' or state == 'mX':
            RO_list = [6]
        if state == 'Y' or state == 'mY':
            RO_list = [3,4,5,6]
        if state == 'Z' or state == 'mZ':
            RO_list = [6,0,1,2]
        for RO in RO_list:
            dataset_dict[state]['Tomo_'+str(RO)] = {}

            if only_4 == False:
                dataset_dict[state]['Tomo_'+str(RO)] = no_QEC_sum_data_single_state_RO(run_list = run_list,idle = idle,state = state,RO = RO,load_set = True,add_4 = add_4,do_weighted = do_weighted)
            elif only_4 == True:
                dataset_dict[state]['Tomo_'+str(RO)] = load_QEC_dataset_single_sign_single_elRO(sym = 'no_correction', RO = RO, state = state, error_sign = 1,el_RO = 'positive', run = 4)



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

    if state_list==['Z','mZ','Y','mY', 'X','mX']:
        toff_dict['toff_process_y'] = 1/4.+1/8.*(toff_dict['toff_'+'Z'+'y']-toff_dict['toff_'+'mZ'+'y']+toff_dict['toff_'+'Y'+'y']-toff_dict['toff_'+'mY'+'y']+toff_dict['toff_'+'X'+'y']-toff_dict['toff_'+'mX'+'y'])
        toff_dict['toff_process_y_err'] = 1/8.*(toff_dict['toff_'+'Z'+'y_err']**2+toff_dict['toff_'+'mZ'+'y_err']**2+toff_dict['toff_'+'Y'+'y_err']**2+toff_dict['toff_'+'mY'+'y_err']**2+toff_dict['toff_'+'X'+'y_err']**2+toff_dict['toff_'+'mX'+'y_err']**2)**0.5
    # print toff_dict['toff_process_y'][0]

    toff_dict['x'] = dataset_dict[state]['Tomo_'+str(RO)]['x']
    return toff_dict

def no_QEC_toffoli_fids_one_state(state = 'Z',idle = False):

    dataset_dict = {}
    toff_dict = {}


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


    toff_dict['x'] = dataset_dict[state]['Tomo_'+str(RO)]['x']
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

def QEC_plot_process_fids_11_vs_idle(run_QEC = [5,6,7],run_enc = [2,3,4]):

    process_dict = QEC_process_fids_sum_runs(run_list = run_QEC,no_error = '11')
    fig,ax = plt.subplots()
    x = process_dict['x']
    folder  = r'D:\measuring\data\QEC_data\figs\11vsIdle'
    process_dict_idle = no_QEC_process_fids_sum_runs(idle = True, run_list = [2,3,4])

    y = process_dict['dec_'+'avg'+'_y']
    y_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
    ax.plot(x_fit, y_fit, 'r', lw=1, label =  'QEC, $p_c$='+str(int(p_c*1000)/1000.)+'('+str(int(p_c_err*1000))+')')
    ax.errorbar(x,y,yerr=y_err,color = 'r',ls = '',marker = 'o',ms = 2)#,label =  'QEC')

    y_undo = process_dict['dec_'+'avg'+'_y_new']
    y_undo_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit_undo, y_fit_undo, p_c, p_c_err = fit_QEC_process_curve(x,y_undo,return_errorbar = True)
    ax.plot(x_fit_undo, y_fit_undo, 'r',ls = ':', lw=1, label =  'Undo QEC, $p_c$='+str(int(p_c*1000)/1000.)+'('+str(int(p_c_err*1000))+')')
    ax.errorbar(x,y_undo,yerr=y_undo_err,color = 'r', ls = '',marker = 'o',ms = 2)#,label =  'undo QEC')


    y_idle = process_dict_idle['dec_'+'avg'+'_y']
    y_idle_err = process_dict_idle['dec_'+'avg'+'_y_err']
    x_fit_idle, y_fit_idle, p_err = fit_QEC_process_curve(x,y_idle)    
    ax.plot(x_fit_idle, y_fit_idle, 'k',ls = '-', lw=1, label =  'Idle, $p_c$='+str(int(p_err*100)/100.)) 
    ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'k', ls = '',marker = 'o',ms = 2)#,label =  'Idle')

    ax.set_ylim(0.3,0.7)
    ax.set_xlim(-0.01,0.51)
    ax.set_title('all_summed_process_fidelity dec to '+'avg_QEC_Run_'+str(run_QEC)+'_Enc_Run_'+str(run_enc)+'.png')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Process Fidelity')
    ax.legend()#loc = 2, bbox_to_anchor = (1,1))

    try:
        fig.savefig(
            os.path.join(folder,'11_vs_idle_QEC_Run_'+str(run_QEC)+'_Enc_Run_'+str(run_enc)+'.png'))
    except:
        print 'Figure has not been saved.'

    fig,ax = plt.subplots()
    ax.errorbar(x,y-y_idle,yerr=(y_err**2+y_idle_err**2)**0.5,color = 'b',ls = '',marker = 'o',ms = 8)
    ax.plot(x_fit,y_fit-y_fit_idle,color = 'b',ls = '-')
    ax.set_ylim(0.0,0.04)
    ax.set_xlim(-0.01,0.51)
    ax.set_title('11 vs idle difference.png')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Deviation')
    try:
        fig.savefig(
            os.path.join(folder,'11_vs_idle_QEC_difference.png'))
    except:
        print 'Figure has not been saved.'

def QEC_plot_process_fids_11_vs_idle_full_curve(run_QEC = [3],run_enc = [0]):
    fig,ax = plt.subplots()
    process_dict = QEC_process_fids_sum_runs(run_list = [3],no_error = '11')

    x = process_dict['x']
    folder  = r'D:\measuring\data\QEC_data\figs\11vsIdle'
    process_dict_idle = no_QEC_process_fids(idle = True, run = 0)

    y = process_dict['dec_'+'avg'+'_y']
    y_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
    ax.plot(x_fit, y_fit, 'r', lw=1, label =  'QEC, $p_c$='+str(int(p_c*1000)/1000.)+'('+str(int(p_c_err*1000))+')')
    ax.errorbar(x,y,yerr=y_err,color = 'r',ls = '',marker = 'o',ms = 2)#,label =  'QEC')

    y_undo = process_dict['dec_'+'avg'+'_y_new']
    y_undo_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit_undo, y_fit_undo, p_c, p_c_err = fit_QEC_process_curve(x,y_undo,return_errorbar = True)
    ax.plot(x_fit_undo, y_fit_undo, 'r',ls = ':', lw=1, label =  'Undo QEC, $p_c$='+str(int(p_c*1000)/1000.)+'('+str(int(p_c_err*1000))+')')
    ax.errorbar(x,y_undo,yerr=y_undo_err,color = 'r', ls = '',marker = 'o',ms = 2)#,label =  'undo QEC')


    y_idle = process_dict_idle['dec_'+'avg'+'_y']
    y_idle_err = process_dict_idle['dec_'+'avg'+'_y_err']
    x_fit_idle, y_fit_idle, p_err = fit_QEC_process_curve(x,y_idle)    
    ax.plot(x_fit_idle, y_fit_idle, 'k',ls = '-', lw=1, label =  'Idle, $p_c$='+str(int(p_err*100)/100.)) 
    ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'k', ls = '',marker = 'o',ms = 2)#,label =  'Idle')

    ax.set_ylim(0,1)
    ax.set_xlim(-0.02,1.02)
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = 'b')
    ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = 'b')
    ax.set_title('all_summed_process_fidelity dec to '+'avg_QEC_Run_'+str(run_QEC)+'_Enc_Run_'+str(run_enc)+'.png')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Process Fidelity')
    ax.legend()#loc = 2, bbox_to_anchor = (1,1))

    try:
        fig.savefig(
            os.path.join(folder,'11_vs_idle_QEC_Run_'+str(run_QEC)+'_Enc_Run_'+str(run_enc)+'.png'))
    except:
        print 'Figure has not been saved.'



def QEC_plot_process_fids_sum_new_fits(append_no_QEC =True, syndrome_list = ['00','01','10','11']):

    process_dict = QEC_process_fids_sum_all(syndrome_list = syndrome_list)
    toff_process_dict = no_QEC_toffoli_fids()
    x = process_dict['x']
    folder  = r'D:\measuring\data\QEC_data\figs\process fidelities2'

    t_list = ['1','2','3','avg']
    color_list = ['c','r','b','g']

    if append_no_QEC == True:
       no_process_dict = no_QEC_process_fids(run = 0)
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

            x_fit, y_fit, p_err = fit_QEC_process_curve(x,y)
            ax.errorbar(x,y,yerr=y_err,color = color_list[i], ls = '',marker = '*')
            ax.plot(x_fit, y_fit, color_list[i], lw=1, ls = ':',  label =  'no QEC, decode to '+ t_list[i]+ ', p_correction: '+str(int(p_err*100)/100.))

            # x_fit, y_idle_fit, p_err = fit_QEC_process_curve(x,y_idle)
            # ax.errorbar(x,y_idle,yerr=y_idle_err,color = color_list[i], ls = '', marker = '*')
            # ax.plot(x_fit, y_idle_fit, color_list[i], lw=1, ls = '--', label =  'with idling, no QEC, decode to '+ t_list[i]+ ', p_correction: '+str(int(p_err*100)/100.))


            if i ==3:
                y = toff_process_dict['toff_process_y']
                y_err = toff_process_dict['toff_process_y_err']
                x_fit, y_fit, p_err = fit_QEC_process_curve(x,y)
                ax.errorbar(x,y,yerr=y_err,color = 'k', ls = '',marker = 'o')
                ax.plot(x_fit, y_fit, 'k', lw=1, ls = ':',label =  'toffoli decoded'+ ', p_correction: '+str(int(p_err*100)/100.))

                y_idle = toff_dict_idle['toff_process_y']
                y_idle_err = toff_dict_idle['toff_process_y_err']
                x_fit, y_idle_fit, p_err = fit_QEC_process_curve(x,y_idle)
                ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'k',  ls = '',marker = '*')

                ax.plot(x_fit, y_idle_fit, 'k', lw=1, ls = '--', label =  'with idling, toffoli decoded'+ ', p_correction: '+str(int(p_err*100)/100.))

    ax.set_ylim(-0,1)
    ax.set_xlim(-0.02,1.02)
    ax.set_title('all_summed_process_fidelity dec to '+t_list[i]+'.png')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = 'b')
    ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = 'b')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Process Fidelity')
    lgd = ax.legend()#loc = 2, bbox_to_anchor = (1,1))

    try:
        fig.savefig(
            os.path.join(folder,'all_summed_process_fidelity_fit'+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
    except:
        print 'Figure has not been saved.'




    for i in range(4):
        fig,ax = plt.subplots(figsize = (10,8))
        y = process_dict['dec_'+t_list[i]+'_y']
        y_err = process_dict['dec_'+t_list[i]+'_y_err']
        x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
        ax.plot(x_fit, y_fit, 'r', lw=1, label =  'QEC, $p_c$='+str(int(p_c*1000)/1000.)+'('+str(int(p_c_err*1000))+')')
        ax.errorbar(x,y,yerr=y_err,color = 'r',ls = '',marker = 'o')
        # x_g = [x[0],x[-1]]
        # y_g = [y[0],y[-1]]
        # ax.plot(x_g,y_g,color = 'k', linestyle = ':' )

        y = process_dict['dec_'+t_list[i]+'_y_new']
        y_err = process_dict['dec_'+t_list[i]+'_y_err']
        x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
        ax.plot(x_fit, y_fit, 'r',ls = ':', lw=1, label =  'Undo QEC, $p_c$='+str(int(p_c*1000)/1000.)+'('+str(int(p_c_err*1000))+')')
        ax.errorbar(x,y,yerr=y_err,color = 'r', ls = '',marker = '.')

        y = no_process_dict['dec_'+t_list[i]+'_y']
        y_err = no_process_dict['dec_'+t_list[i]+'_y_err']

        x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
        ax.plot(x_fit, y_fit, 'b',ls = '-', lw=1, label =  'Encode, $p_c$='+str(int(p_c*100)/100.)+'('+str(int(p_c_err*100))+')')
        ax.errorbar(x,y,yerr=y_err,color = 'b', ls = '',marker ='.')

        y = toff_process_dict['toff_process_y']
        y_err = toff_process_dict['toff_process_y_err']
        x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
        ax.plot(x_fit, y_fit, 'b',ls = '-', lw=1,label =   'Encode + majority vote, $p_c$='+str(int(p_c*100)/100.)+'('+str(int(p_c_err*100))+')')
        ax.errorbar(x,y,yerr=y_err,color = 'b', ls = '',marker = 'o')

        # y_idle = process_dict_idle['dec_'+t_list[i]+'_y']
        # y_idle_err = process_dict_idle['dec_'+t_list[i]+'_y_err']
        # x_fit, y_fit, p_err = fit_QEC_process_curve(x,y_idle)
        # ax.plot(x_fit, y_fit, 'k',ls = '-', lw=1, label =  'with idling, no QEC, decode to '+ t_list[i]+ ', p_correction: '+str(int(p_err*100)/100.))
        # ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'k', ls = '',marker = '.')

        # y_idle = toff_dict_idle['toff_process_y']
        # y_idle_err = toff_dict_idle['toff_process_y_err']
        # x_fit, y_fit, p_err = fit_QEC_process_curve(x,y_idle)
        # ax.plot(x_fit, y_fit, 'k',ls = '-', lw=1, label =  'with idling, toffoli decoded'+ ', p_correction: '+str(int(p_err*100)/100.))
        # ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'k',  ls = '', marker = '*')

        y = single_process_dict['dec_'+t_list[i]+'_y']
        y_err = process_dict['dec_'+t_list[i]+'_y_err']
        x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
        ax.plot(x_fit, y_fit, 'g',ls = '-', lw=1, label =  'single Qubit, $p_c$='+str(int(p_c*1000)/1000.)+'('+str(int(p_c_err*1000))+')')
        ax.errorbar(x,y,yerr=y_err,color = 'g', ls = '',marker = 'o')

        ax.set_ylim(-0,1)
        ax.set_xlim(-0.02,1.02)
        ax.set_title('QEC process fidelities, dec to ' +t_list[i]+', syndrome_list '+str(syndrome_list))
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = 'b')
        ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = 'b')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Process Fidelity')
        lgd = ax.legend()#loc = 2, bbox_to_anchor = (1,1))

        try:
            fig.savefig(
                os.path.join(folder,'all_summed_process_fidelity_dec_to_' +t_list[i]+'syndrome_list'+str(syndrome_list)+'.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
        except:
            print 'Figure has not been saved.'            



# def plot_prob_timesweep_single_syndrome(syndrome = '11'):
#     folder  = r'D:\measuring\data\QEC_data\figs\Probabilities'

#     p_list = ['p11','p01','p00','p10']
#     p_dict = QEC_timesweep_sum_probs(no_error = syndrome)

#     fig,ax = plt.subplots()
#     color =  [c_green,c_orange,c_red,'r']
#     label_list = ['no error', 'Qubit 1','Qubit 2', 'Qubit 3']

#     for jj,p in enumerate(p_list):
#         ax.plot(p_dict['x'],p_dict[p],color = color[jj],label = label_list[jj])

#     ax.legend()
#     ax.set_ylim(-0.01,1.01)
#     ax.set_xlim(-1e-3,35e-3)
#     print 'ok'
#     ax.hlines([0.301401517586],x[0]-1,x[-1]+1,linestyles='dotted', color = c_green,lw = 0.5)
#     ax.hlines([0.203400853559],x[0]-1,x[-1]+1,linestyles='dotted', color = c_red,lw = 0.5)
#     ax.hlines([0.247598809228],x[0]-1,x[-1]+1,linestyles='dotted', color = c_orange,lw = 0.5)
#     ax.hlines([0.247598819488],x[0]-1,x[-1]+1,linestyles='dotted', color = 'r',lw = 0.5)
#     ax.set_xlabel('Time (s)')
#     ax.set_ylabel('Measured outcome probability')    

#     try:
#         fig.savefig(
#             os.path.join(folder,'timesweep_probability_plot_'+syndrome+'.png'))
#     except:
#         print 'Figure has not been saved.' 



def plot_prob_different_error():
    syndrome_list = ['00','01','10','11']
    p_dict = {}
    folder  = r'D:\measuring\data\QEC_data\figs\Probabilities'
    for syndrome in syndrome_list:

        if syndrome == '00' or syndrome == '01':
            run_list = [1,2,3]
        elif syndrome == '10':
            run_list = [2]
        elif syndrome == '11':
            run_list = [3]

        p_list = ['p00','p01','p10','p11']
        p_dict[syndrome] = QEC_sum_probs(run_list = run_list ,no_error = syndrome)

    p_dict['no qubit'] = 1/4.*(p_dict['00']['p00']+p_dict['01']['p01']+p_dict['10']['p10']+p_dict['11']['p11'])
    p_dict['Q2'] = 1/4.*(p_dict['00']['p01']+p_dict['01']['p00']+p_dict['10']['p11']+p_dict['11']['p10'])
    p_dict['Q3'] = 1/4.*(p_dict['00']['p10']+p_dict['01']['p11']+p_dict['10']['p00']+p_dict['11']['p01'])
    p_dict['Q1'] = 1/4.*(p_dict['00']['p11']+p_dict['01']['p10']+p_dict['10']['p01']+p_dict['11']['p00'])

    fig,ax = plt.subplots()
    for error in ['no qubit','Q1','Q2','Q3']:
        ax.plot(p_dict[syndrome]['x'],p_dict[error],label = 'Error on ' + error)
    ax.legend()
    ax.set_xlim([-0.01,1.01])
    ax.set_ylim([-0.01,1.01])
    ax.set_xlabel('Error probability')
    ax.set_ylabel('Measured error probability')    

    try:
        fig.savefig(
            os.path.join(folder,'error_probabilities.png'))
    except:
        print 'Figure has not been saved.'


def plot_timesweep_prob_different_error():
    syndrome_list = ['00','01','10','11']
    p_dict = {}
    folder  = r'D:\measuring\data\QEC_data\figs\Probabilities'
    for syndrome in syndrome_list:

        p_list = ['p00','p01','p10','p11']
        p_dict[syndrome] = QEC_timesweep_sum_probs(no_error = syndrome)

    p_dict['no qubit'] = 1/4.*(p_dict['00']['p00']+p_dict['01']['p01']+p_dict['10']['p10']+p_dict['11']['p11'])
    p_dict['Q2'] = 1/4.*(p_dict['00']['p01']+p_dict['01']['p00']+p_dict['10']['p11']+p_dict['11']['p10'])
    p_dict['Q3'] = 1/4.*(p_dict['00']['p10']+p_dict['01']['p11']+p_dict['10']['p00']+p_dict['11']['p01'])
    p_dict['Q1'] = 1/4.*(p_dict['00']['p11']+p_dict['01']['p10']+p_dict['10']['p01']+p_dict['11']['p00'])

    fig,ax = plt.subplots()
    for error in ['no qubit','Q1','Q2','Q3']:
        ax.plot(p_dict[syndrome]['x'],p_dict[error],label = 'Error on ' + error)
    ax.legend()
    ax.set_ylim(-0.01,1.01)
    ax.set_xlim(-1e-3,35e-3)
    ax.set_xlabel('Time (s)')
    ax.set_xlabel('Error probability')
    ax.set_ylabel('Measured error probability')    

    try:
        fig.savefig(
            os.path.join(folder,'timesweep_error_probabilities.png'))
    except:
        print 'Figure has not been saved.'

def plot_prob_fids_all():
    p_dict = {}
    p_list = ['p00','p01','p10','p11']
    syndrome_list = ['00','01','10','11']
    n = 0
    color = ['c','k','m','b']
    n_list = []

    folder  = r'D:\measuring\data\QEC_data\figs\Probabilities'
    fig,ax = plt.subplots(figsize=(10,8))
    for i,syndrome in enumerate(syndrome_list):
        pp = []
        if syndrome == '00' or syndrome == '01':
            run_list = [1,2,3]
        elif syndrome == '10':
            run_list = [2]
        elif syndrome == '11':
            run_list = [3]
        p_dict[syndrome] = QEC_sum_probs(run_list = run_list ,no_error = syndrome)

        for p in p_list:
            pp += [p_dict[syndrome][p][0]]
            n_list +=[n]
            n+=1

        ax.bar(n_list[n-4:n],pp,color = color[i],align ='center',label = 'no error assigned to '+syndrome)

    pp = []
    for p in p_list:
        pp +=[ 1/4.*(p_dict['00'][p][5]+p_dict['01'][p][5]+p_dict['10'][p][5]+p_dict['11'][p][5])]
        n_list +=[n]
        n+=1
    ax.bar(n_list[n-4:n],pp,color = '0.5',align ='center',label = 'average outcome at pe = 0.5')

    tick_list = 5*syndrome_list
    ax.set_xticks(n_list)
    ax.set_xticklabels(tick_list)

    ax.set_xlim([n_list[0]-0.6,n_list[-1]+0.6])
    ax.set_ylim([0,1])
    ax.set_xlabel('Parity measurement outcome')
    ax.set_ylabel('Probability')
    ax.legend(loc = 2)
    try:
        fig.savefig(
            os.path.join(folder,'probability_barplot.png'))
    except:
        print 'Figure has not been saved.'

def no_QEC_plot_toff_fids():

    toff_dict1 = no_QEC_toffoli_fids(state_list = ['Z','mZ'],run_list=[0,1,2],add_4 = False)
    toff_dict2 = no_QEC_toffoli_fids(state_list = ['Y','mY', 'X','mX'],run_list=[0],add_4 = False)
    # toff_dict_idle = no_QEC_toffoli_fids(idle = True)

    color =['b','g','r','c','m','k']

    folder  = r'D:\measuring\data\QEC_data\figs\Encoding'

    fig,ax = plt.subplots()
    for i, state in enumerate(['Z','mZ','Y','mY', 'X','mX']):
        if state in ['Z','mZ']:
            toff_dict = toff_dict1
        elif state in ['Y','mY', 'X','mX']:
            toff_dict = toff_dict2
        x = toff_dict['x']
        y = toff_dict['toff_'+state+'y']
        if state[0] == 'm':
            y = -1*y
        y_err = toff_dict['toff_'+state+'y_err']
        if state in ['Z','mZ','Y','mY']:
            x_fit, y_fit, p_err = fit_QEC_curve(x,y)
            ax.errorbar(x,y,yerr = y_err, ls = '', marker = 'o',color = color[i],label = 'state '+state +' p_c = '+str(int(p_err*100)/100.))
            ax.plot(x_fit,y_fit,color[i])
        else:ax.errorbar(x,y,yerr = y_err, marker = 'o',color = color[i],label = 'state '+state)


    ax.set_title('toffoli_decoded.png')
    ax.legend(loc = 3)
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')

    # ax.plot([0,1],[1,-1],'k')


    try:
        fig.savefig(
            os.path.join(folder,'toffoli_decoded_all.png'))
    except:
        print 'Figure has not been saved.'

def no_QEC_sum_Z_mZ():
    toff_process_dict = no_QEC_toffoli_fids(run_list = [0,1,2,3],state_list=['Z','mZ'],add_4 = True,do_weighted = True)
    y_Z1              = toff_process_dict['toff_Zy']
    y_err_Z1          = toff_process_dict['toff_Zy_err']
    y_mZ1             = toff_process_dict['toff_mZy']
    y_err_mZ1         = toff_process_dict['toff_mZy_err']
    x1                = toff_process_dict['x']

    fig,ax = plt.subplots()
    ax.errorbar(x1, (y_Z1-y_mZ1)/2., yerr=(y_err_Z1**2+y_err_mZ1**2)**0.5/2.,color = 'g', label = '1 round' )

    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.hlines([-1,0,1],x1[0]-1,x1[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')
    folder  = r'D:\measuring\data\QEC_data\figs\Encoding'
    # ax.plot([0,1],[1,-1],'k')
    try:
        fig.savefig(
            os.path.join(folder,'toffoli_decoded_state_Z_mZ_summed'+'.png'))
    except:
        print 'Figure has not been saved.'

def no_QEC_compare_Z_mZ(run_list = [0,1,2,3]):
    fig,ax = plt.subplots()
    for run in run_list:
        if run !=4:
            toff_process_dict = no_QEC_toffoli_fids(run_list = [run],state_list=['Z','mZ'],add_4 = False)
            y_Z1              = toff_process_dict['toff_Zy']
            y_err_Z1          = toff_process_dict['toff_Zy_err']
            y_mZ1             = toff_process_dict['toff_mZy']
            y_err_mZ1         = toff_process_dict['toff_mZy_err']
            x1                = toff_process_dict['x']


            ax.errorbar(x1, (y_Z1-y_mZ1)/2., yerr=(y_err_Z1**2+y_err_mZ1**2)**0.5/2., label = '1 round, run '+str(run) )
        elif run ==4:
            toff_process_dict = no_QEC_toffoli_fids(run_list = [run],state_list=['Z','mZ'],only_4 = True)
            y_Z1              = toff_process_dict['toff_Zy']
            y_err_Z1          = toff_process_dict['toff_Zy_err']
            y_mZ1             = toff_process_dict['toff_mZy']
            y_err_mZ1         = toff_process_dict['toff_mZy_err']
            x1                = toff_process_dict['x']

            ax.errorbar(x1, (y_Z1-y_mZ1)/2., yerr=(y_err_Z1**2+y_err_mZ1**2)**0.5/2., label = '1 round, run '+str(run) )


    toff_process_dict = no_QEC_toffoli_fids(run_list = [0,1,2,3],state_list=['Z','mZ'],add_4 = True,do_weighted = True)
    y_Z1              = toff_process_dict['toff_Zy']
    y_err_Z1          = toff_process_dict['toff_Zy_err']
    y_mZ1             = toff_process_dict['toff_mZy']
    y_err_mZ1         = toff_process_dict['toff_mZy_err']
    x1                = toff_process_dict['x']



    ### Averaging over Z and mZ
    y_tot1  = (y_Z1-y_mZ1)/2; y_err_tot1 = (y_err_Z1**2+y_err_mZ1**2)**0.5/2

    ### Fitting
    x_fit1, y_fit1, fit_result1,u_fit_result1 = fit_QEC_curve(x1,y_tot1, return_errorbar=True)

    ax.plot(x_fit1, (y_fit1),color = 'k')

    ax.errorbar(x1, (y_Z1-y_mZ1)/2., yerr=(y_err_Z1**2+y_err_mZ1**2)**0.5/2.,ls = '', label = 'weighted sum')
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    # ax.set_xlim(0.37,0.53)
    # ax.set_ylim(-0.1,0.25)
    ax.hlines([-1,0,1],x1[0]-1,x1[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.legend()
    ax.set_title('toffoli_decoded_state_Z_mZ'+'.png')
    ax.set_ylabel('Contrast')
    folder  = r'D:\measuring\data\QEC_data\figs\Encoding'
    try:
        fig.savefig(
            os.path.join(folder,'compare_toffoli_decoded_state_Z_mZ_summed'+'.png'))
    except:
        print 'Figure has not been saved.'

def no_QEC_plot_toff_fids_compare_runs(state_list = ['Z','mZ','Y','mY', 'X','mX'],run_list = [0]):

    folder  = r'D:\measuring\data\QEC_data\figs\Encoding'


    for state in state_list:
        fig,ax = plt.subplots()
        for run in run_list:

            if run !=4:
                toff_dict = no_QEC_toffoli_fids(run_list = [run],state_list=state_list)
                # toff_dict_idle = no_QEC_toffoli_fids(idle = True)

                x = toff_dict['x']
                y = toff_dict['toff_'+state+'y']
                y_err = toff_dict['toff_'+state+'y_err']
                ax.errorbar(x,y,yerr = y_err, label = 'toffoli state '+state + 'run '+str(run))

        if 4 in run_list:
            run =4
            toff_dict = no_QEC_toffoli_fids(run_list = [run],state_list=state_list,only_4 =True)
            # toff_dict_idle = no_QEC_toffoli_fids(idle = True)

            x = toff_dict['x']
            y = toff_dict['toff_'+state+'y']
            y_err = toff_dict['toff_'+state+'y_err']
            ax.errorbar(x,y,yerr = y_err, label = 'toffoli state '+state + 'run '+str(run))

        ax.set_title('toffoli_decoded_state_'+state+'.png')

        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        # ax.set_xlim(0.37,0.53)
        # ax.set_ylim(-0.25,0.1)
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')
        ax.legend()
        try:
            fig.savefig(
                os.path.join(folder,'compare_toffoli_decoded_state_'+state+'.png'))
        except:
            print 'Figure has not been saved.'

def no_QEC_plot_process_fids():

    process_dict = no_QEC_process_fids(run = 0)
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




########################################################################
''' analysis of 2Round experiments with singly applied error, by THT '''
########################################################################

''' This data set is found at older_than = '20150114_012244' '''

def QEC_2rounds0_load_data(older_than = None, load_from_data = False, len_k = 2):
    if load_from_data:
        data = {}
        for state in ['Z','mZ']:
            for RO in [0,1,2,6]:
                for syndrome in ['00','01','10','11']:
                    data[state + str(RO) + syndrome], folder = QEC_data_single_state_RO(older_than = older_than,state = state,
                                                                                        RO = RO, sym = syndrome, len_k=len_k)
        pickle.dump(data, open( "2rounds0.p", "wb" ) )
    else:
        data = pickle.load( open( "2rounds0.p", "rb" ) )

    return data

def QEC_2rounds0_combine_syndromes(data):

    for state in ['Z','mZ']:
        for RO in [0,1,2,6]:
            for i, syndrome in enumerate(['00','01','10','11']):

                if i == 0:
                    y       = data[state + str(RO) + syndrome]['y']
                    y_err   = data[state + str(RO) + syndrome]['y_err']**2

                    y_pos       = data[state + str(RO) + syndrome]['y_pos']
                    y_pos_err   = data[state + str(RO) + syndrome]['y_pos_err']**2

                    y_neg       = data[state + str(RO) + syndrome]['y_neg']
                    y_neg_err   = data[state + str(RO) + syndrome]['y_neg_err']**2

                else:
                    y      = y+data[state + str(RO) + syndrome]['y']
                    y_err  = y_err+data[state + str(RO) + syndrome]['y_err']**2

                    y_pos       = y_pos     + data[state + str(RO) + syndrome]['y_pos']
                    y_pos_err   = y_pos_err + data[state + str(RO) + syndrome]['y_pos_err']**2

                    y_neg       = y_neg     + data[state + str(RO) + syndrome]['y_neg']
                    y_neg_err   = y_neg_err + data[state + str(RO) + syndrome]['y_neg_err']**2

            y       = y/4.; y_pos       = y_pos/4.; y_neg       = y_neg/4.
            y_err   = (y_err**0.5)/4. ; y_pos_err   = (y_pos_err**0.5)/4. ; y_neg_err   = (y_neg_err**0.5)/4.

            data[state + str(RO)] = {}
            data[state + str(RO)]['y'] = y
            data[state + str(RO)]['y_err'] = y_err
            data[state + str(RO)]['y_pos'] = y_pos
            data[state + str(RO)]['y_pos_err'] = y_pos_err
            data[state + str(RO)]['y_neg'] = y_neg
            data[state + str(RO)]['y_neg_err'] = y_neg_err
            data[state + str(RO)]['x'] = data[state + str(RO) + syndrome]['x']

    return data

def QEC_2rounds0_apply_final_error(data):

    for state in ['Z','mZ']:
        for RO in [0,1,2,6]:
            error_list         = (1 - (1-2*data[state + str(RO)]['x'])**0.5 )/2.
            contrast_reduction =  1 - 2*error_list

            if RO == 6:
                power = 3
            else:
                power =1

            for syndrome in ['', '00','01','10','11']:
                data[state + str(RO) + syndrome]['y']          = data[state + str(RO) + syndrome]['y']     * contrast_reduction**power
                data[state + str(RO) + syndrome]['y_err']      = data[state + str(RO) + syndrome]['y_err'] * contrast_reduction**power
            for syndrome in ['']:
                data[state + str(RO) + syndrome]['y_pos']      = data[state + str(RO) + syndrome]['y_pos']     * contrast_reduction**power
                data[state + str(RO) + syndrome]['y_pos_err']  = data[state + str(RO) + syndrome]['y_pos_err'] * contrast_reduction**power
                data[state + str(RO) + syndrome]['y_neg']      = data[state + str(RO) + syndrome]['y_neg']     * contrast_reduction**power
                data[state + str(RO) + syndrome]['y_neg_err']  = data[state + str(RO) + syndrome]['y_neg_err'] * contrast_reduction**power

    return data

def QEC_2rounds0_apply_final_QEC(data):

    for state in ['Z','mZ']:

        for syndrome in ['', '00','01','10','11']:

            data[state + syndrome] = {}

            for error_sign in ['','_pos','_neg']:
                data[state + syndrome]['y' + error_sign]  =      (data[state + '0' + syndrome]['y' + error_sign] +
                                          data[state + '1' + syndrome]['y' + error_sign] +
                                          data[state + '2' + syndrome]['y' + error_sign] -
                                          data[state + '6' + syndrome]['y' + error_sign])/2.

                data[state + syndrome]['y'+error_sign + '_err'] = ((data[state + '0'+syndrome]['y'+error_sign + '_err']**2 +
                                         data[state + '1'+syndrome]['y'+error_sign + '_err']**2 +
                                         data[state + '2'+syndrome]['y'+error_sign + '_err']**2 +
                                         data[state + '6'+syndrome]['y'+error_sign + '_err']**2)**0.5)/2.
    return data

def QEC_2rounds0_analysis(older_than = '20150114_012244',load_from_data = False,len_k = 2):
    ''' this functions returns a dictionairy that contains all data from the QEC_2rounds0_experiments.
        TO DO: describe the dictionairy entries
        General_idea: anything not mentioned is averaged over'''

    data_dict  = QEC_2rounds0_load_data(older_than = older_than, load_from_data = load_from_data,len_k = len_k)
    data_dict  = QEC_2rounds0_combine_syndromes(data_dict)
    data_dict  = QEC_2rounds0_apply_final_error(data_dict)
    data_dict  = QEC_2rounds0_apply_final_QEC(data_dict)

    return data_dict

def QEC_2rounds0_plot(older_than = '20150114_012244',load_from_data = False, save_folder = r'D:\measuring\data\QEC_data\figs'):

    data_dict = QEC_2rounds0_analysis(older_than=older_than, load_from_data=load_from_data)
    x          = data_dict['Z0']['x']
    y_Z        = data_dict['Z']['y']
    y_err_Z    = data_dict['Z']['y_err']
    y_mZ       = data_dict['mZ']['y']
    y_err_mZ   = data_dict['mZ']['y_err']

    ### Full result
    if 1:
        fig1,ax = plt.subplots()
        ax.errorbar(x, y_Z, yerr=y_err_Z,color = 'b' )
        ax.plot([x[0],x[-1]], [y_Z[0],y_Z[-1]],'k:' )
        ax.plot([0,0.5], [1,0],'k:' )
        ax.set_ylim(-0.05,1.05)
        ax.set_xlim(-0.05,0.55)
        ax.set_title('QEC_data_2Rounds0 older_than = ' + older_than + '\n' + 'State = Z')
        ax.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Expectation value')

        try:
            fig1.savefig(
                os.path.join(save_folder,'StateZ.png'))
        except:
            print 'Figure has not been saved.'

        fig2,ax2 = plt.subplots()
        ax2.errorbar(x, y_mZ, yerr=y_err_mZ,color = 'b' )
        ax2.plot([x[0],x[-1]], [y_mZ[0],y_mZ[-1]],'k:' )
        ax2.plot([0,0.5], [-1,0],'k:' )
        ax2.set_ylim(-1.05,0.05)
        ax2.set_xlim(-0.05,0.55)
        ax2.set_title('QEC_data_2Rounds0 older_than = ' + older_than + '\n' + 'State = -Z')
        ax2.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax2.set_xlabel('error probability')
        ax2.set_ylabel('Expectation value')

        try:
            fig2.savefig(
                os.path.join(save_folder,'StatemZ.png'))
        except:
            print 'Figure has not been saved.'

    ### Comparing 4 error syndromes
    if 1:
        fig3,ax3 = plt.subplots()
        ax3.errorbar(x, data_dict['Z00']['y'], yerr=y_err_Z,color = 'b' )
        ax3.errorbar(x, data_dict['Z01']['y'], yerr=y_err_Z,color = 'm' )
        ax3.errorbar(x, data_dict['Z10']['y'], yerr=y_err_Z,color = 'g' )
        ax3.errorbar(x, data_dict['Z11']['y'], yerr=y_err_Z,color = 'y' )

        ax3.plot([x[0],x[-1]], [y_Z[0],y_Z[-1]],color = 'k' )
        ax3.plot([0,0.5], [1,0],color = 'k' )
        ax3.set_ylim(-0.05,1.05)
        ax3.set_xlim(-0.05,0.55)
        ax3.set_title('QEC_data_2Rounds0 Zstate older_than = ' + older_than + '\n' + 'State = Z, 4 error syndromes')
        ax3.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax3.set_xlabel('error probability')
        ax3.set_ylabel('Expectation value')

        try:
            fig2.savefig(
                os.path.join(save_folder,'4Syndromes.png'))
        except:
            print 'Figure has not been saved.'



    ### Positive vs negative errors
    if 1:
        fig5,ax5 = plt.subplots()
        ax5.errorbar(x, data_dict['Z']['y_pos'], yerr=data_dict['Z']['y_pos_err'],color = 'b' )
        ax5.errorbar(x, data_dict['Z']['y_neg'], yerr=data_dict['Z']['y_neg_err'],color = 'k' )
        # ax5.plot([x[0],x[-1]], [y_Z[0],y_Z[-1]],color = 'k' )
        # ax5.plot([0,0.5], [1,0],color = 'k' )
        ax5.set_ylim(-1.05,1.05)
        ax5.set_xlim(-0.05,0.55)
        ax5.set_title('QEC_data_2Rounds Zstate older_than = ' + older_than + '\n' + 'State = Z, pos vs neg errors')
        ax5.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax5.set_xlabel('error probability')
        ax5.set_ylabel('Expectation value')

        fig6,ax6 = plt.subplots()
        ax6.errorbar(x, data_dict['Z']['y_pos'], yerr=data_dict['Z']['y_pos_err'],color = 'b' )
        ax6.errorbar(x, data_dict['Z']['y_neg'], yerr=data_dict['Z']['y_neg_err'],color = 'k' )
        # ax6.plot([x[0],x[-1]], [y_Z[0],y_Z[-1]],color = 'k' )
        # ax6.plot([0,0.5], [1,0],color = 'k' )
        ax6.set_ylim(-1.05,1.05)
        ax6.set_xlim(-0.05,0.55)
        ax6.set_title('QEC_data_2Rounds Zstate older_than = ' + older_than + '\n' + 'State = Z, pos vs neg errors')
        ax6.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax6.set_xlabel('error probability')
        ax6.set_ylabel('Expectation value')

        try:
            fig2.savefig(
                os.path.join(save_folder,'PosVSNeg_errors.png'))
        except:
            print 'Figure has not been saved.'


    ### Averaging the states
    if 1:
        fig4,ax4 = plt.subplots()
        ax4.errorbar(x, (y_Z-y_mZ)/2., yerr=(y_err_Z**2+y_err_mZ**2)**0.5/2.,color = 'b' )

        ax4.plot([x[0],x[-1]], [y_Z[0],y_Z[-1]],color = 'k' )
        ax4.plot([0,0.5], [1,0],color = 'k' )
        ax4.set_ylim(-0.05,1.05)
        ax4.set_xlim(-0.05,0.55)
        ax4.set_title('QEC_data_2Rounds Zstate older_than = ' + older_than +  '\n' + 'State = (Z+mZ)/2')
        ax4.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax4.set_xlabel('error probability')
        ax4.set_ylabel('Expectation value')

        try:
            fig2.savefig(
                os.path.join(save_folder,'Analyze_script.png'))
        except:
            print 'Figure has not been saved.'

    plt.show()
    plt.close('all')


#########################################################################
''' Analysis of 2Round experiments with doubly applied errors, by THT '''
#########################################################################

''' This data set is found at:

    Set 1: 20150114_232135 to 20150116_085607. Contains pe = [0,0.1,0.2,0.3,0.4,0.5]
    Set 2: 20150125_183842 to 20150125_234653. Contains pe = [0.2, 0.3, 0.4, 0.45]
    Set 3: 20150201_011834 to 20150201_132340. Contains pe = [0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5]
    Currently only analyzing the '11' error syndrome
'''

def QEC_2rounds_get_data_dict(older_than = None, RO = 0, state = 'Z',  sym = '11', error_signs = '11',electron_RO = 'positive', len_k = range(6)):
    QEC_dict = {}
    k_dict = {}

    for k in len_k:

        timestamp, folder = toolbox.latest_data(contains = '2Rounds_syn_' + sym + '_' + electron_RO +'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_signs)
                                                +'_'+state, older_than = older_than,return_timestamp = True)
        print folder
        SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
        print SSRO_folder

        k_dict['k_'+str(k)] ={}
        k_dict['k_'+str(k)] = load_QEC_data(folder, SSRO_folder, post_select = True)

    for item in k_dict['k_' + str(len_k[0])]:

        if len(len_k) ==2:
            QEC_dict[item] = np.concatenate((k_dict['k_'+str(len_k[0])][item],k_dict['k_'+str(len_k[1])][item]), axis=0)
        elif len(len_k) ==3:
            QEC_dict[item] = np.concatenate((k_dict['k_'+str(len_k[0])][item],k_dict['k_'+str(len_k[1])][item],k_dict['k_'+str(len_k[2])][item]), axis=0)
        elif len(len_k) == 4:
            QEC_dict[item] = np.concatenate((k_dict['k_'+str(len_k[0])][item],k_dict['k_'+str(len_k[1])][item],k_dict['k_'+str(len_k[2])][item], k_dict['k_'+str(len_k[3])][item]), axis=0)
        elif len(len_k) == 6:
            QEC_dict[item] = np.concatenate((k_dict['k_'+str(len_k[0])][item],k_dict['k_'+str(len_k[1])][item],k_dict['k_'+str(len_k[2])][item], k_dict['k_'+str(len_k[3])][item], k_dict['k_'+str(len_k[4])][item], k_dict['k_'+str(len_k[5])][item]), axis=0)

    return QEC_dict,folder

def QEC_2rounds_load_data(run = 1, load_from_data = False):
    ''' example of a key for the dictioniary: 'ZRO6S11E-1-1eROneg'
        read as: State = Z, RO = 6, Syndrome = 11,  Error_signs = -1,-1, electron_RO = negative.
        Keys with less values generally are averaged over the missing degree of freedom'''

    if run ==1:
        older_than = '20150116_085607'
        len_k = range(2)
    elif run ==2:
        older_than = '20150125_234653'
        len_k = [4,5]
    elif run ==3:
        older_than = '20150201_132340'
        len_k = [6,7,8]


    if load_from_data:
        data = {}
        for state in ['Z','mZ']:
            for RO in [0,1,2,6]:
                for syndrome in ['11']:
                    for error_signs in ['11','1-1','-11','-1-1']:
                        for electron_RO in ['positive', 'negative']:
                            if electron_RO =='positive':
                                eRO = 'pos'
                            elif electron_RO == 'negative':
                                eRO = 'neg'

                            data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eRO'+eRO], folder = QEC_2rounds_get_data_dict(
                                                                                        older_than = older_than,
                                                                                        state = state,
                                                                                        RO = RO, sym = syndrome, error_signs = error_signs,
                                                                                        electron_RO = electron_RO,
                                                                                        len_k=len_k)
        if run == 1:
            pickle.dump(data, open( "2rounds_run1.p", "wb" ) )
        if run == 2:
            pickle.dump(data, open( "2rounds_run2.p", "wb" ) )
        if run == 3:
            pickle.dump(data, open( "2rounds_run3.p", "wb" ) )

    else:
        if run == 1:
            data = pickle.load( open( "2rounds_run1.p", "rb" ) )
        elif run == 2:
            data = pickle.load( open( "2rounds_run2.p", "rb" ) )
        elif run == 3:
            data = pickle.load( open( "2rounds_run3.p", "rb" ) )

    return data

def QEC_2rounds_combine_eRO(data):
    '''combines the positive and negative electron RO data'''

    for state in ['Z','mZ']:
            for RO in [0,1,2,6]:
                for syndrome in ['11']:
                    for error_signs in ['11','1-1','-11','-1-1']:

                            y     = (data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROpos']['c0'] -
                                    data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROneg']['c0'])/2.
                            y_err = (data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROpos']['c0_u']**2 +
                                    data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROneg']['c0_u']**2)**0.5/2.

                            data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs] = {}
                            data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y'] = y
                            data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y_err'] = y_err
                            data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['x'] = data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROpos']['x']

                            prop_list = ['p00','p01','p10','p11']

                            for ii in range(4):
                                data[state + 'RO'+str(RO) + 'S'+syndrome + 'E' + error_signs][prop_list[ii]] = (
                                         data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROpos'][prop_list[ii]] +
                                         data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROneg'][prop_list[ii]])/2.

    return data

def QEC_2rounds_combine_error_signs(data):
    '''combines the 4 different erro sign combinations'''

    prop_list = ['p00','p01','p10','p11']

    for state in ['Z','mZ']:
        for RO in [0,1,2,6]:
            for syndrome in ['11']:
                data[state + 'RO'+str(RO) + 'S'+syndrome] = {}
                for i,error_signs in enumerate(['11','1-1','-11','-1-1']):

                    if i == 0:
                        y     = data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y']
                        y_err = data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y_err']**2

                        for ii in range(len(prop_list)):
                            data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]] = data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs][prop_list[ii]]
                    else:
                        y     = y + data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y']
                        y_err = y_err + data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y_err']**2

                        for ii in range(len(prop_list)):
                            data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]] = (data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]] +
                                                                                    data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs][prop_list[ii]])
                y = y/4.; y_err = (y_err**0.5)/4.

                data[state + 'RO'+str(RO) + 'S'+syndrome]['y']     = y
                data[state + 'RO'+str(RO) + 'S'+syndrome]['y_err'] = y_err
                data[state + 'RO'+str(RO) + 'S'+syndrome]['x']     = data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['x']

                for ii in range(len(prop_list)):
                    data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]] = data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]]/4.

    return data

def QEC_2rounds_combine_syndromes(data):
    prop_list = ['p00','p01','p10','p11']

    for state in ['Z','mZ']:
        for RO in [0,1,2,6]:
            data[state + 'RO'+str(RO)] = {}
            for i, syndrome in enumerate(['11']):

                if i == 0:
                    y       = data[state + 'RO'+str(RO) + 'S'+syndrome]['y']
                    y_err   = data[state + 'RO'+str(RO) + 'S'+syndrome]['y_err']**2

                    for ii in range(len(prop_list)):
                        data[state + 'RO'+str(RO)][prop_list[ii]] = data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]]

                else:
                    y      = y+data[state + 'RO'+str(RO) + 'S'+syndrome]['y']
                    y_err  = y_err+data[state + 'RO'+str(RO) + 'S'+syndrome]['y_err']**2

                    for ii in range(len(prop_list)):
                        data[state + 'RO'+str(RO)][prop_list[ii]] = (data[state + 'RO'+str(RO)][prop_list[ii]] +
                                                                                    data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]])


            y = y; y_err   = (y_err**0.5)

            data[state + 'RO'+str(RO)]['y'] = y
            data[state + 'RO'+str(RO)]['y_err'] = y_err
            data[state + 'RO'+str(RO)]['x'] = data[state + 'RO'+str(RO) + 'S'+syndrome]['x']

            for ii in range(len(prop_list)):
                data[state + 'RO'+str(RO)][prop_list[ii]] = data[state + 'RO'+str(RO)][prop_list[ii]]

    return data

def QEC_2rounds_apply_final_QEC(data):

    prop_list = ['p00','p01','p10','p11']

    for state in ['Z','mZ']:
        for syndrome in ['']:#, 'S00','S01','S10','S11']:
            for error_sign in ['']:#,'E11','E-11','E1-1','E-1-1']:
                data[state + syndrome + error_sign] = {}

                data[state + syndrome + error_sign]['y']  =   (data[state + 'RO0' + syndrome + error_sign]['y'] +
                                          data[state + 'RO1' + syndrome + error_sign]['y'] +
                                          data[state + 'RO2' + syndrome + error_sign]['y'] -
                                          data[state + 'RO6' + syndrome + error_sign]['y'])/2.

                data[state + syndrome + error_sign]['y_err'] = ((data[state + 'RO0'+syndrome + error_sign]['y_err']**2 +
                                         data[state + 'RO1'+syndrome + error_sign]['y_err']**2 +
                                         data[state + 'RO2'+syndrome + error_sign]['y_err']**2 +
                                         data[state + 'RO6'+syndrome + error_sign]['y_err']**2)**0.5)/2.

                data[state + syndrome + error_sign]['x'] = data[state + 'RO0' + syndrome + error_sign]['x']

                for ii in range(len(prop_list)):
                    data[state + syndrome + error_sign][prop_list[ii]] = (data[state+'RO0'+syndrome][prop_list[ii]] +
                                                                          data[state+'RO1'+syndrome][prop_list[ii]] +
                                                                          data[state+'RO2'+syndrome][prop_list[ii]] +
                                                                          data[state+'RO6'+syndrome][prop_list[ii]])/4
    return data

def QEC_2rounds_analysis(run =1,load_from_data = False):
    ''' this functions returns a dictionairy that contains all data from the QEC_2rounds_experiments.
        TO DO: describe the dictionairy entries
        General_idea: anything not mentioned is averaged over'''

    data_dict  = QEC_2rounds_load_data(run = run, load_from_data = load_from_data)
    data_dict  = QEC_2rounds_combine_eRO(data_dict)
    data_dict  = QEC_2rounds_combine_error_signs(data_dict)
    data_dict  = QEC_2rounds_combine_syndromes(data_dict)
    data_dict  = QEC_2rounds_apply_final_QEC(data_dict)

    return data_dict

def QEC_2rounds_combined_runs(runs=[1,2]):

    x          = {}
    y_Z        = {}
    y_err_Z    = {}
    y_mZ       = {}
    y_err_mZ   = {}

    if len(runs) == 1:
        data_dict1 = QEC_2rounds_analysis(run=runs[0], load_from_data=False)
        x          = data_dict1['Z']['x']
        y_Z        = data_dict1['Z']['y']
        y_err_Z    = data_dict1['Z']['y_err']
        y_mZ       = data_dict1['mZ']['y']
        y_err_mZ   = data_dict1['mZ']['y_err']

    elif len(runs) == 2:
        data_dict1  = QEC_2rounds_analysis(run=runs[0], load_from_data=False)
        x1          = data_dict1['Z']['x']
        y_Z1        = data_dict1['Z']['y']
        y_err_Z1    = data_dict1['Z']['y_err']
        y_mZ1       = data_dict1['mZ']['y']
        y_err_mZ1   = data_dict1['mZ']['y_err']

        data_dict2  = QEC_2rounds_analysis(run=runs[1], load_from_data=False)
        x2          = data_dict2['Z']['x']
        y_Z2        = data_dict2['Z']['y']
        y_err_Z2    = data_dict2['Z']['y_err']
        y_mZ2       = data_dict2['mZ']['y']
        y_err_mZ2   = data_dict2['mZ']['y_err']

        weight_Z1    =  1/(y_err_Z1)**2
        weight_mZ1   =  1/(y_err_mZ1)**2
        
        weight_Z2    =  1/(y_err_Z2)**2
        weight_mZ2   =  1/(y_err_mZ2)**2


        if runs[0] == 1 and runs[1] == 2:
            x   = np.array([ 0.,   0.1,  0.2,  0.3,  0.4,  0.45, 0.5])

            ## Combine without weight
            y_Z     = np.array([y_Z1[0],    y_Z1[1],  (y_Z1[2]+y_Z2[0])/2.,   (y_Z1[3]+y_Z2[1])/2.,   (y_Z1[4]+y_Z2[2])/2.,   y_Z2[3],  y_Z1[5]])
            y_mZ    = np.array([y_mZ1[0],   y_mZ1[1], (y_mZ1[2]+y_mZ2[0])/2., (y_mZ1[3]+y_mZ2[1])/2., (y_mZ1[4]+y_mZ2[2])/2., y_mZ2[3], y_mZ1[5]])

            y_err_Z = np.array([y_err_Z1[0],   y_err_Z1[1],  (y_err_Z1[2]**2+y_err_Z2[0]**2)**0.5/2.,   (y_err_Z1[3]**2+y_err_Z2[1]**2)**0.5/2.,   (y_err_Z1[4]**2+y_err_Z2[2]**2)**0.5/2.,   y_err_Z2[3],  y_err_Z1[5]])
            y_err_mZ = np.array([y_err_mZ1[0],   y_err_mZ1[1],  (y_err_mZ1[2]**2+y_err_mZ2[0]**2)**0.5/2.,   (y_err_mZ1[3]**2+y_err_mZ2[1]**2)**0.5/2.,   (y_err_mZ1[4]**2+y_err_mZ2[2]**2)**0.5/2.,   y_err_mZ2[3],  y_err_mZ1[5]])

            ## Weighted average
            y_Z     = np.array([y_Z1[0],
                                y_Z1[1],  
                                np.average( [y_Z1[2],y_Z2[0]], weight = [weight_Z1[2],weight_Z2[0]]),   
                                np.average( [y_Z1[3],y_Z2[1]], weight = [weight_Z1[3],weight_Z2[1]]),   
                                np.average( [y_Z1[4],y_Z2[2]], weight = [weight_Z1[4],weight_Z2[2]]),   
                                y_Z2[3],  
                                y_Z1[5]])
            
            y_err_Z = np.array([y_err_Z1[0],   
                                y_err_Z1[1],  
                                1/(weight_Z1[2]+weight_Z2[0])**0.5,
                                1/(weight_Z1[3]+weight_Z2[1])**0.5,
                                1/(weight_Z1[4]+weight_Z2[2])**0.5,
                                y_err_Z2[3],  
                                y_err_Z1[5]])


            y_mZ     = np.array([y_mZ1[0],
                                y_mZ1[1],  
                                np.average( [y_mZ1[2],y_mZ2[0]], weight = [weight_mZ1[2],weight_mZ2[0]]),   
                                np.average( [y_mZ1[3],y_mZ2[1]], weight = [weight_mZ1[3],weight_mZ2[1]]),   
                                np.average( [y_mZ1[4],y_mZ2[2]], weight = [weight_mZ1[4],weight_mZ2[2]]),   
                                y_mZ2[3],  
                                y_mZ1[5]])

            y_err_mZ = np.array([y_err_mZ1[0],   
                                y_err_mZ1[1],  
                                1/(weight_mZ1[2]+weight_mZ2[0])**0.5,
                                1/(weight_mZ1[3]+weight_mZ2[1])**0.5,
                                1/(weight_mZ1[4]+weight_mZ2[2])**0.5,
                                y_err_mZ2[3],  
                                y_err_mZ1[5]])


        elif runs[0] == 1 and runs[1] == 3:

            x   = np.array([ 0.,   0.1,  0.2,  0.3, 0.35,  0.4,  0.45, 0.5])

            ## Combine without weight
            y_Z = np.array([   (y_Z1[0]+y_Z2[0])/2.,    (y_Z1[1]+y_Z2[1])/2.,    (y_Z1[2]+y_Z2[2])/2.,     (y_Z1[3]+y_Z2[3])/2.,    y_Z2[4],   (y_Z1[4]+y_Z2[5])/2.,   y_Z2[6],  (y_Z1[5]+y_Z2[7])/2.])
            y_mZ = np.array([(y_mZ1[0]+y_mZ2[0])/2.,  (y_mZ1[1]+y_mZ2[1])/2.,  (y_mZ1[2]+y_mZ2[2])/2.,   (y_mZ1[3]+y_mZ2[3])/2.,   y_mZ2[4], (y_mZ1[4]+y_mZ2[5])/2.,  y_mZ2[6],  (y_mZ1[5]+y_mZ2[7])/2.])

            y_err_Z  = np.array([ (y_err_Z1[0]**2+y_err_Z2[0]**2)**0.5/2.,    (y_err_Z1[1]**2+y_err_Z2[1]**2)**0.5/2.,    (y_err_Z1[2]**2+y_err_Z2[2]**2)**0.5/2.,     (y_err_Z1[3]**2+y_err_Z2[3]**2)**0.5/2.,    y_err_Z2[4],   (y_err_Z1[4]**2+y_err_Z2[5]**2)**0.5/2.,   y_err_Z2[6],  (y_err_Z1[5]**2+y_err_Z2[7]**2)**0.5/2.])
            y_err_mZ = np.array([ (y_err_mZ1[0]**2+y_err_mZ2[0]**2)**0.5/2.,  (y_err_mZ1[1]**2+y_err_mZ2[1]**2)**0.5/2.,  (y_err_mZ1[2]**2+y_err_mZ2[2]**2)**0.5/2.,   (y_err_mZ1[3]**2+y_err_mZ2[3]**2)**0.5/2.,   y_err_mZ2[4], (y_err_mZ1[4]**2+y_err_mZ2[5]**2)**0.5/2.,  y_err_mZ2[6],  (y_err_mZ1[5]**2+y_err_mZ2[7]**2)**0.5/2.])

            ## Weighted average
            y_Z = np.array([    np.average( [y_Z1[0],y_Z2[0]], weight = [weight_Z1[0],weight_Z2[0]]),    
                                np.average( [y_Z1[1],y_Z2[1]], weight = [weight_Z1[1],weight_Z2[1]]),    
                                np.average( [y_Z1[2],y_Z2[2]], weight = [weight_Z1[2],weight_Z2[2]]),     
                                np.average( [y_Z1[3],y_Z2[3]], weight = [weight_Z1[3],weight_Z2[3]]),    
                                y_Z2[4],   
                                np.average( [y_Z1[4],y_Z2[5]], weight = [weight_Z1[4],weight_Z2[5]]),   
                                y_Z2[6],  
                                np.average( [y_Z1[5],y_Z2[7]], weight = [weight_Z1[5],weight_Z2[7]])])

            y_err_Z  = np.array([   1/(weight_Z1[0]+weight_Z2[0])**0.5,    
                                    1/(weight_Z1[1]+weight_Z2[1])**0.5,    
                                    1/(weight_Z1[2]+weight_Z2[2])**0.5,     
                                    1/(weight_Z1[3]+weight_Z2[3])**0.5,    
                                    y_err_Z2[4],   
                                    1/(weight_Z1[4]+weight_Z2[5])**0.5,   
                                    y_err_Z2[6],  
                                    1/(weight_Z1[5]+weight_Z2[7])**0.5])

            y_mZ = np.array([   np.average( [y_mZ1[0],y_mZ2[0]], weight = [weight_mZ1[0],weight_mZ2[0]]),    
                                np.average( [y_mZ1[1],y_mZ2[1]], weight = [weight_mZ1[1],weight_mZ2[1]]),    
                                np.average( [y_mZ1[2],y_mZ2[2]], weight = [weight_mZ1[2],weight_mZ2[2]]),     
                                np.average( [y_mZ1[3],y_mZ2[3]], weight = [weight_mZ1[3],weight_mZ2[3]]),    
                                y_mZ2[4],   
                                np.average( [y_mZ1[4],y_mZ2[5]], weight = [weight_mZ1[4],weight_mZ2[5]]),   
                                y_mZ2[6],  
                                np.average( [y_mZ1[5],y_mZ2[7]], weight = [weight_mZ1[5],weight_mZ2[7]])])

            y_err_mZ  = np.array([  1/(weight_mZ1[0]+weight_mZ2[0])**0.5,    
                                    1/(weight_mZ1[1]+weight_mZ2[1])**0.5,    
                                    1/(weight_mZ1[2]+weight_mZ2[2])**0.5,     
                                    1/(weight_mZ1[3]+weight_mZ2[3])**0.5,    
                                    y_err_mZ2[4],   
                                    1/(weight_mZ1[4]+weight_mZ2[5])**0.5,   
                                    y_err_mZ2[6],  
                                    1/(weight_mZ1[5]+weight_mZ2[7])**0.5])

    elif len(runs) == 3:

        data_dict1  = QEC_2rounds_analysis(run=runs[0], load_from_data=False)
        x1          = data_dict1['Z']['x']
        y_Z1        = data_dict1['Z']['y']
        y_err_Z1    = data_dict1['Z']['y_err']
        y_mZ1       = data_dict1['mZ']['y']
        y_err_mZ1   = data_dict1['mZ']['y_err']

        data_dict2  = QEC_2rounds_analysis(run=runs[1], load_from_data=False)
        x2          = data_dict2['Z']['x']
        y_Z2        = data_dict2['Z']['y']
        y_err_Z2    = data_dict2['Z']['y_err']
        y_mZ2       = data_dict2['mZ']['y']
        y_err_mZ2   = data_dict2['mZ']['y_err']

        data_dict3  = QEC_2rounds_analysis(run=runs[2], load_from_data=False)
        x3          = data_dict3['Z']['x']
        y_Z3        = data_dict3['Z']['y']
        y_err_Z3    = data_dict3['Z']['y_err']
        y_mZ3       = data_dict3['mZ']['y']
        y_err_mZ3   = data_dict3['mZ']['y_err']

        weight_Z1    =  1/(y_err_Z1)**2
        weight_mZ1   =  1/(y_err_mZ1)**2
        
        weight_Z2    =  1/(y_err_Z2)**2
        weight_mZ2   =  1/(y_err_mZ2)**2

        weight_Z3    =  1/(y_err_Z3)**2
        weight_mZ3   =  1/(y_err_mZ3)**2

        if runs[0] == 1 and runs[1] == 2 and runs[2] == 3:
            x   = np.array([ 0.,   0.1,  0.2,  0.3, 0.35,  0.4,  0.45, 0.5])

            ### not weighted average
            y_Z = np.array([   (y_Z1[0]+y_Z3[0])/2.,
                                (y_Z1[1]+y_Z3[1])/2.,    (y_Z1[2]+y_Z2[0]+y_Z3[2])/3.,     (y_Z1[3]+y_Z2[1]+y_Z3[3])/3.,
                                y_Z3[4],   (y_Z1[4]+y_Z2[2]+y_Z3[5])/3.,   (y_Z2[3]+y_Z3[6])/2,  (y_Z1[5]+y_Z3[7])/2.])
            
            y_mZ = np.array([   (y_mZ1[0]+y_mZ3[0])/2.,    (y_mZ1[1]+y_mZ3[1])/2.,    (y_mZ1[2]+y_mZ2[0]+y_mZ3[2])/3.,     (y_mZ1[3]+y_mZ2[1]+y_mZ3[3])/3.,
                                y_mZ3[4],   (y_mZ1[4]+y_mZ2[2]+y_mZ3[5])/3.,   (y_mZ2[3]+y_mZ3[6])/2,  (y_mZ1[5]+y_mZ3[7])/2.])

            y_err_Z = np.array([   (y_err_Z1[0]**2+y_err_Z3[0]**2)**0.5/2.,    (y_err_Z1[1]**2+y_err_Z3[1]**2)**0.5/2.,    (y_err_Z1[2]**2+y_err_Z2[0]**2+y_err_Z3[2]**2)**0.5/3.,     (y_err_Z1[3]**2+y_err_Z2[1]**2+y_err_Z3[3]**2)**0.5/3.,
                                y_err_Z3[4],   (y_err_Z1[4]**2+y_err_Z2[2]**2+y_err_Z3[5]**2)**0.5/3.,   (y_err_Z2[3]**2+y_err_Z3[6]**2)**0.5/2,  (y_err_Z1[5]**2+y_err_Z3[7]**2)**0.5/2.])
            y_err_mZ = np.array([   (y_err_mZ1[0]**2+y_err_mZ3[0]**2)**0.5/2.,    (y_err_mZ1[1]**2+y_err_mZ3[1]**2)**0.5/2.,    (y_err_mZ1[2]**2+y_err_mZ2[0]**2+y_err_mZ3[2]**2)**0.5/3.,     (y_err_mZ1[3]**2+y_err_mZ2[1]**2+y_err_mZ3[3]**2)**0.5/3.,
                                y_err_mZ3[4],   (y_err_mZ1[4]**2+y_err_mZ2[2]**2+y_err_mZ3[5]**2)**0.5/3.,   (y_err_mZ2[3]**2+y_err_mZ3[6]**2)**0.5/2,  (y_err_mZ1[5]**2+y_err_mZ3[7]**2)**0.5/2.])

            ### weighted average
            y_Z = np.array([    np.average( [y_Z1[0],y_Z3[0]],          weights = [weight_Z1[0],weight_Z3[0]]),
                                np.average( [y_Z1[1],y_Z3[1]],          weights = [weight_Z1[1],weight_Z3[1]]),    
                                np.average( [y_Z1[2],y_Z2[0],y_Z3[2]],  weights = [weight_Z1[2],weight_Z2[0],weight_Z3[2]]),     
                                np.average( [y_Z1[3],y_Z2[1],y_Z3[3]],  weights = [weight_Z1[3],weight_Z2[1],weight_Z3[3]]),
                                y_Z3[4],   
                                np.average( [y_Z1[4],y_Z2[2],y_Z3[5]],  weights = [weight_Z1[4],weight_Z2[2],weight_Z3[5]]),   
                                np.average( [y_Z2[3],y_Z3[6]],          weights = [weight_Z2[3],weight_Z3[6]]),  
                                np.average( [y_Z1[5],y_Z3[7]],          weights = [weight_Z1[5],weight_Z3[7]])])

            y_err_Z = np.array([    1/(weight_Z1[0]+weight_Z3[0])**0.5,    
                                    1/(weight_Z1[1]+weight_Z3[1])**0.5,    
                                    1/(weight_Z1[2]+weight_Z2[0]+weight_Z3[2])**0.5,     
                                    1/(weight_Z1[3]+weight_Z2[1]+weight_Z3[3])**0.5,
                                    y_err_Z3[4],   
                                    1/(weight_Z1[4]+weight_Z2[2]+weight_Z3[5])**0.5,   
                                    1/(weight_Z2[3]+weight_Z3[6])**0.5,  
                                    1/(weight_Z1[5]+weight_Z3[7])**0.5])

            y_mZ = np.array([    np.average( [y_mZ1[0],y_mZ3[0]],          weights = [weight_mZ1[0],weight_mZ3[0]]),
                                np.average( [y_mZ1[1],y_mZ3[1]],          weights = [weight_mZ1[1],weight_mZ3[1]]),    
                                np.average( [y_mZ1[2],y_mZ2[0],y_mZ3[2]],  weights = [weight_mZ1[2],weight_mZ2[0],weight_mZ3[2]]),     
                                np.average( [y_mZ1[3],y_mZ2[1],y_mZ3[3]],  weights = [weight_mZ1[3],weight_mZ2[1],weight_mZ3[3]]),
                                y_mZ3[4],   
                                np.average( [y_mZ1[4],y_mZ2[2],y_mZ3[5]],  weights = [weight_mZ1[4],weight_mZ2[2],weight_mZ3[5]]),   
                                np.average( [y_mZ2[3],y_mZ3[6]],          weights = [weight_mZ2[3],weight_mZ3[6]]),  
                                np.average( [y_mZ1[5],y_mZ3[7]],          weights = [weight_mZ1[5],weight_mZ3[7]])])

            y_err_mZ = np.array([    1/(weight_mZ1[0]+weight_mZ3[0])**0.5,    
                                    1/(weight_mZ1[1]+weight_mZ3[1])**0.5,    
                                    1/(weight_mZ1[2]+weight_mZ2[0]+weight_mZ3[2])**0.5,     
                                    1/(weight_mZ1[3]+weight_mZ2[1]+weight_mZ3[3])**0.5,
                                    y_err_mZ3[4],   
                                    1/(weight_mZ1[4]+weight_mZ2[2]+weight_mZ3[5])**0.5,   
                                    1/(weight_mZ2[3]+weight_mZ3[6])**0.5,  
                                    1/(weight_mZ1[5]+weight_mZ3[7])**0.5])
            
    return x, y_Z, y_err_Z, y_mZ, y_err_mZ

def QEC_2rounds_plot_final_curves(runs = [1,3], load_from_data = False, save_folder = r'D:\measuring\data\QEC_data\figs\multiple_rounds'):


    x, y_Z, y_err_Z, y_mZ, y_err_mZ =  QEC_2rounds_combined_runs(runs=runs)

    ### Z and -Z states seperately
    if 1:
        fig1,ax = plt.subplots()
        ax.errorbar(x, y_Z, yerr=y_err_Z,color = 'b' )
        ax.plot([x[0],x[-1]], [y_Z[0],y_Z[-1]],'k:' )
        ax.plot([0,0.5], [1,0],'k:' )
        ax.set_ylim(-0.05,1.05)
        ax.set_xlim(-0.05,0.55)
        ax.set_title('QEC_data_2Rounds' + '\n' + 'State = Z')
        ax.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Expectation value')

        if save_folder != None:
            try:
                fig1.savefig(
                    os.path.join(save_folder,'StateZ.png'))
            except:
                print 'Figure has not been saved.'

        fig2,ax2 = plt.subplots()
        ax2.errorbar(x, y_mZ, yerr=y_err_mZ,color = 'b' )
        ax2.plot([x[0],x[-1]], [y_mZ[0],y_mZ[-1]],'k:' )
        ax2.plot([0,0.5], [-1,0],'k:' )
        ax2.set_ylim(-1.05,0.05)
        ax2.set_xlim(-0.05,0.55)
        ax2.set_title('QEC_data_2Rounds' + '\n' + 'State = -Z')
        ax2.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax2.set_xlabel('error probability')
        ax2.set_ylabel('Expectation value')

        if save_folder != None:
            try:
                fig2.savefig(
                    os.path.join(save_folder,'StatemZ.png'))
            except:
                print 'Figure has not been saved.'

    ### Averaging the +Z and -Z states
    if 1:
        fig4,ax4 = plt.subplots()
        ax4.errorbar(x, (y_Z-y_mZ)/2., yerr=(y_err_Z**2+y_err_mZ**2)**0.5/2.,color = 'b' )
        ax4.plot([x[0],x[-1]], [(y_Z[0]-y_mZ[0])/2,(y_Z[-1]-y_mZ[-1])/2],'k:' )
        ax4.plot([0,0.5], [1,0],'k:' )
        ax4.set_ylim(-0.05,1.05)
        ax4.set_xlim(-0.05,0.55)
        ax4.set_title('QEC_data_2Rounds Zstate' +  '\n' + 'State = (Z+mZ)/2')
        ax4.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax4.set_xlabel('error probability')
        ax4.set_ylabel('Expectation value')

        if save_folder != None:
            try:
                fig2.savefig(
                    os.path.join(save_folder,'Analyze_script.png'))
            except:
                print 'Figure has not been saved.'

    ### Compare runs
    if 0:
        x1, y_Z1, y_err_Z1, y_mZ1, y_err_mZ1 =  QEC_2rounds_combined_runs(runs=[1])
        x2, y_Z2, y_err_Z2, y_mZ2, y_err_mZ2 =  QEC_2rounds_combined_runs(runs=[2])

        fig4,ax4 = plt.subplots()
        ax4.errorbar(x1, (y_Z1-y_mZ1)/2., yerr=(y_err_Z1**2+y_err_mZ1**2)**0.5/2.,color = 'b',ls ='', label = 'run1')
        ax4.errorbar(x2, (y_Z2-y_mZ2)/2., yerr=(y_err_Z2**2+y_err_mZ2**2)**0.5/2.,color = 'r',ls ='', label = 'run2' )

        ax4.plot([x[0],x[-1]], [(y_Z[0]-y_mZ[0])/2,(y_Z[-1]-y_mZ[-1])/2],'k:' )
        ax4.plot([0,0.5], [1,0],'k:' )
        ax4.set_ylim(-0.05,1.05)
        ax4.set_xlim(-0.05,0.55)
        ax4.set_title('QEC_data_2Rounds Zstate' +  '\n' + 'State = (Z+mZ)/2')
        ax4.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax4.set_xlabel('error probability')
        ax4.set_ylabel('Expectation value')
        ax4.legend()

        if save_folder != None:
            try:
                fig4.savefig(
                    os.path.join(save_folder,'2Rounds_compare_runs.pdf'))
            except:
                print 'Figure has not been saved.'

    ### Compare with and without run 2
    if 0:
        x1, y_Z1, y_err_Z1, y_mZ1, y_err_mZ1 =  QEC_2rounds_combined_runs(runs=[1])
        x2, y_Z2, y_err_Z2, y_mZ2, y_err_mZ2 =  QEC_2rounds_combined_runs(runs=[1,2])

        fig4,ax4 = plt.subplots()
        ax4.errorbar(x1, (y_Z1-y_mZ1)/2., yerr=(y_err_Z1**2+y_err_mZ1**2)**0.5/2.,color = 'b',ls ='', label = 'only run 1')
        ax4.errorbar(x2, (y_Z2-y_mZ2)/2., yerr=(y_err_Z2**2+y_err_mZ2**2)**0.5/2.,color = 'r',ls ='', label = 'combined' )

        ax4.plot([x[0],x[-1]], [(y_Z[0]-y_mZ[0])/2,(y_Z[-1]-y_mZ[-1])/2],'k:' )
        ax4.plot([0,0.5], [1,0],'k:' )
        ax4.set_ylim(-0.05,1.05)
        ax4.set_xlim(-0.05,0.55)
        ax4.set_title('QEC_data_2Rounds Zstate' +  '\n' + 'State = (Z+mZ)/2')
        ax4.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax4.set_xlabel('error probability')
        ax4.set_ylabel('Expectation value')
        ax4.legend()

        if save_folder != None:
            try:
                fig4.savefig(
                    os.path.join(save_folder,'2Rounds_with_without_run2.pdf'))
            except:
                print 'Figure has not been saved.'

    ### Outcome probabilities
    if 0:
        data_dict = QEC_2rounds_analysis(run=1, load_from_data=False)
        x         = data_dict['Z']['x']

        data_dict2 = QEC_2rounds_analysis(run=2, load_from_data=False)
        x2         = data_dict2['Z']['x']


        prop_list = ['p00','p01','p10','p11']
        colors = cm.rainbow(np.linspace(0, 1, len(prop_list)))

        fig,ax = plt.subplots()
        for ii in range(len(prop_list)):
            p = data_dict['Z'][prop_list[ii]]
            ax.plot(x,p,'o',ls='',color = colors[ii], label = 'run1_' + prop_list[ii])

        colors = cm.rainbow(np.linspace(0, 1, len(prop_list)))
        for ii in range(len(prop_list)):
            p = data_dict2['Z'][prop_list[ii]]
            ax.plot(x2,p,'b', color = colors[ii], label = 'run_2' + prop_list[ii])
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.show()
    plt.close('all')

def QEC_2rounds_plot_other_curves(older_than = '20150116_085607',load_from_data = False, save_folder = r'D:\measuring\data\QEC_data\figs'):

    data_dict  = QEC_2rounds_analysis(older_than=older_than, load_from_data=load_from_data)
    x          = data_dict['Z']['x']
    y_Z        = data_dict['Z']['y']
    y_err_Z    = data_dict['Z']['y_err']
    y_mZ       = data_dict['mZ']['y']
    y_err_mZ   = data_dict['mZ']['y_err']

    # ### Comparing 4 error syndromes
    # if 1:
    #     fig3,ax3 = plt.subplots()
    #     ax3.errorbar(x, data_dict['Z00']['y'], yerr=y_err_Z,color = 'b' )
    #     ax3.errorbar(x, data_dict['Z01']['y'], yerr=y_err_Z,color = 'm' )
    #     ax3.errorbar(x, data_dict['Z10']['y'], yerr=y_err_Z,color = 'g' )
    #     ax3.errorbar(x, data_dict['Z11']['y'], yerr=y_err_Z,color = 'y' )

    #     ax3.plot([x[0],x[-1]], [y_Z[0],y_Z[-1]],color = 'k' )
    #     ax3.plot([0,0.5], [1,0],color = 'k' )
    #     ax3.set_ylim(-0.05,1.05)
    #     ax3.set_xlim(-0.05,0.55)
    #     ax3.set_title('QEC_data_2Rounds Zstate older_than = ' + older_than + '\n' + 'State = Z, 4 error syndromes')
    #     ax3.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
    #     ax3.set_xlabel('error probability')
    #     ax3.set_ylabel('Expectation value')

    #     try:
    #         fig2.savefig(
    #             os.path.join(save_folder,'4Syndromes.png'))
    #     except:
    #         print 'Figure has not been saved.'

    # ### Positive vs negative errors
    # if 1:
    #     fig5,ax5 = plt.subplots()
    #     ax5.errorbar(x, data_dict['Z']['y_pos'], yerr=data_dict['Z']['y_pos_err'],color = 'b' )
    #     ax5.errorbar(x, data_dict['Z']['y_neg'], yerr=data_dict['Z']['y_neg_err'],color = 'k' )
    #     # ax5.plot([x[0],x[-1]], [y_Z[0],y_Z[-1]],color = 'k' )
    #     # ax5.plot([0,0.5], [1,0],color = 'k' )
    #     ax5.set_ylim(-1.05,1.05)
    #     ax5.set_xlim(-0.05,0.55)
    #     ax5.set_title('QEC_data_2Rounds Zstate older_than = ' + older_than + '\n' + 'State = Z, pos vs neg errors')
    #     ax5.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
    #     ax5.set_xlabel('error probability')
    #     ax5.set_ylabel('Expectation value')

    #     fig6,ax6 = plt.subplots()
    #     ax6.errorbar(x, data_dict['Z']['y_pos'], yerr=data_dict['Z']['y_pos_err'],color = 'b' )
    #     ax6.errorbar(x, data_dict['Z']['y_neg'], yerr=data_dict['Z']['y_neg_err'],color = 'k' )
    #     # ax6.plot([x[0],x[-1]], [y_Z[0],y_Z[-1]],color = 'k' )
    #     # ax6.plot([0,0.5], [1,0],color = 'k' )
    #     ax6.set_ylim(-1.05,1.05)
    #     ax6.set_xlim(-0.05,0.55)
    #     ax6.set_title('QEC_data_2Rounds Zstate older_than = ' + older_than + '\n' + 'State = Z, pos vs neg errors')
    #     ax6.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
    #     ax6.set_xlabel('error probability')
    #     ax6.set_ylabel('Expectation value')

    #     try:
    #         fig2.savefig(
    #             os.path.join(save_folder,'PosVSNeg_errors.png'))
    #     except:
    #         print 'Figure has not been saved.'

###############################################################################
''' Analysis of 3Round experiments with triply applied errors, added by THT '''
###############################################################################

''' This data set is found at:
    Set 1: older_than = '20150122_043634', Contains pe = [0,0.1,0.2,0.3,0.4,0.5]
    Set 2: 20150125_001141 to 20150125_145530, Contains pe = [0.2, 0.3, 0.4, 0.45]
'''

def QEC_3rounds_get_data_dict(older_than = None, RO = 0, state = 'Z',  sym = '11', error_signs = '11',lectron_RO = 'positive', len_k = range(6)):
    QEC_dict = {}
    k_dict = {}

    for k in len_k:

        print older_than
        print ('3Rounds_syn_' + sym + '_' + electron_RO +'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_signs)
                                                +'_'+state)

        timestamp, folder = toolbox.latest_data(contains = '3Rounds_syn_' + sym + '_' + electron_RO +'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_signs)
                                                +'_'+state, older_than = older_than,return_timestamp = True)
        print folder
        SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
        print SSRO_folder

        k_dict['k_'+str(k)] ={}
        k_dict['k_'+str(k)] = load_QEC_data(folder, SSRO_folder, post_select = True, nr_of_parity_msmnts=4)

    for item in k_dict['k_' + str(len_k[0])]:

        if len(len_k) ==2:
            QEC_dict[item] = np.concatenate((k_dict['k_'+str(len_k[0])][item],k_dict['k_'+str(len_k[1])][item]), axis=0)
        elif len(len_k) == 4:
            QEC_dict[item] = np.concatenate((k_dict['k_'+str(len_k[0])][item],k_dict['k_'+str(len_k[1])][item],k_dict['k_'+str(len_k[2])][item], k_dict['k_'+str(len_k[3])][item]), axis=0)
        elif len(len_k) == 6:
            QEC_dict[item] = np.concatenate((k_dict['k_'+str(len_k[0])][item],k_dict['k_'+str(len_k[1])][item],k_dict['k_'+str(len_k[2])][item], k_dict['k_'+str(len_k[3])][item], k_dict['k_'+str(len_k[4])][item], k_dict['k_'+str(len_k[5])][item]), axis=0)

    return QEC_dict,folder

def QEC_3rounds_load_data(run = 1, load_from_data = False):
    ''' example of a key for the dictioniary: 'ZRO6S11E-1-1eROneg'
        read as: State = Z, RO = 6, Syndrome = 11,  Error_signs = -1,-1, electron_RO = negative.
        Keys with less values generally are averaged over the missing degree of freedom'''

    if run ==1:
        older_than = '20150123_031715'
        len_k = range(6)
    elif run ==2:
        older_than = '20150125_145530'
        len_k = [2,3,4,6]


    if load_from_data:
        data = {}
        for state in ['Z','mZ']:
            for RO in [0,1,2,6]:
                for syndrome in ['11']:#,'01','10','11']:
                    for error_signs in ['111','1-11','-111','-1-11', '11-1','1-1-1','-11-1','-1-1-1']:
                        for electron_RO in ['positive', 'negative']:
                            if electron_RO =='positive':
                                eRO = 'pos'
                            elif electron_RO == 'negative':
                                eRO = 'neg'

                            data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eRO'+eRO], folder = QEC_3rounds_get_data_dict(
                                                                                        older_than = older_than,
                                                                                        state = state,
                                                                                        RO = RO, sym = syndrome, error_signs = error_signs,
                                                                                        electron_RO = electron_RO,
                                                                                        len_k=len_k)
        if run == 1:
            pickle.dump(data, open( "3rounds_run1.p", "wb" ) )
        elif run == 2:
            pickle.dump(data, open( "3rounds_run2.p", "wb" ) )

    else:
        if run == 1:
            data = pickle.load( open( "3rounds_run1.p", "rb" ) )
        elif run == 2:
            data = pickle.load( open( "3rounds_run2.p", "rb" ) )


    return data

def QEC_3rounds_combine_eRO(data):
    '''combines the positive and negative electron RO data'''

    for state in ['Z','mZ']:
            for RO in [0,1,2,6]:
                for syndrome in ['11']:#,'11','01']:#,'01','10','11']:
                    for error_signs in ['111','1-11','-111','-1-11', '11-1','1-1-1','-11-1','-1-1-1']:

                            # RO data (not postselected)
                            y     = (data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROpos']['c0'] -
                                    data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROneg']['c0'])/2.
                            y_err = (data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROpos']['c0_u']**2 +
                                    data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROneg']['c0_u']**2)**0.5/2.

                            data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs] = {}
                            data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y'] = y
                            data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y_err'] = y_err
                            data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['x'] = data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROpos']['x']

                            prop_list = ['p0000','p0100','p1000','p1100','p0010','p0110','p1010','p1110',
                                         'p0001','p0101','p1001','p1101','p0011','p0111','p1011','p1111']
                            for ii in range(16):
                                data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs][prop_list[ii]] = (
                                         data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROpos'][prop_list[ii]] +
                                         data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs+ 'eROneg'][prop_list[ii]])/2.

    return data

def QEC_3rounds_combine_error_signs(data):
    '''combines the 4 different erro sign combinations'''

    prop_list = ['p0000','p0100','p1000','p1100','p0010','p0110','p1010','p1110',
                 'p0001','p0101','p1001','p1101','p0011','p0111','p1011','p1111']

    for state in ['Z','mZ']:
        for RO in [0,1,2,6]:
            for syndrome in ['11']:#,'11','01']:#,'01','10','11']:
                data[state + 'RO'+str(RO) + 'S'+syndrome] = {}

                for i,error_signs in enumerate(['111','1-11','-111','-1-11', '11-1','1-1-1','-11-1','-1-1-1']):

                    if i == 0:
                        y     = data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y']
                        y_err = data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y_err']**2

                        for ii in range(len(prop_list)):
                            data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]] = data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs][prop_list[ii]]

                    else:
                        y     = y + data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y']
                        y_err = y_err + data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['y_err']**2

                        for ii in range(len(prop_list)):
                            data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]] = (data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]] +
                                                                                        data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs][prop_list[ii]])

                y = y/8.; y_err = (y_err**0.5)/8.

                data[state + 'RO'+str(RO) + 'S'+syndrome]['y']     = y
                data[state + 'RO'+str(RO) + 'S'+syndrome]['y_err'] = y_err
                data[state + 'RO'+str(RO) + 'S'+syndrome]['x']     = data[state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs]['x']

                for ii in range(len(prop_list)):
                    data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]] = data[state + 'RO'+str(RO) + 'S'+syndrome][prop_list[ii]]/8.

    return data

def QEC_3rounds_combine_syndromes(data):

    for state in ['Z','mZ']:
        for RO in [0,1,2,6]:
            data[state + 'RO'+str(RO)] = {}
            for i, syndrome in enumerate(['11']):#['00','01','10','11']):

                if i == 0:
                    y       = data[state + 'RO'+str(RO) + 'S'+syndrome]['y']
                    y_err   = data[state + 'RO'+str(RO) + 'S'+syndrome]['y_err']**2
                else:
                    y      = y+data[state + 'RO'+str(RO) + 'S'+syndrome]['y']
                    y_err  = y_err+data[state + 'RO'+str(RO) + 'S'+syndrome]['y_err']**2

            # y       = y/2.; y_err   = (y_err**0.5)/2.
            y       = y; y_err   = (y_err**0.5)

            data[state + 'RO'+str(RO)]['y'] = y
            data[state + 'RO'+str(RO)]['y_err'] = y_err
            data[state + 'RO'+str(RO)]['x'] = data[state + 'RO'+str(RO) + 'S'+syndrome]['x']

    return data

def QEC_3rounds_apply_final_QEC(data):

    # state + 'RO'+str(RO) + 'S'+syndrome + 'E'+error_signs

    for state in ['Z','mZ']:
        for syndrome in ['']:#['', 'S00','S01','S10','S11']:
            for error_sign in ['']:#,'E11','E-11','E1-1','E-1-1']:

                data[state +  error_sign] = {}

                data[state + syndrome + error_sign]['y']  =   (data[state + 'RO0' + syndrome + error_sign]['y'] +
                                          data[state + 'RO1' + syndrome + error_sign]['y'] +
                                          data[state + 'RO2' + syndrome + error_sign]['y'] -
                                          data[state + 'RO6' + syndrome + error_sign]['y'])/2.

                data[state + syndrome + error_sign]['y_err'] = ((data[state + 'RO0'+syndrome + error_sign]['y_err']**2 +
                                         data[state + 'RO1' + syndrome + error_sign]['y_err']**2 +
                                         data[state + 'RO2' + syndrome + error_sign]['y_err']**2 +
                                         data[state + 'RO6' + syndrome + error_sign]['y_err']**2)**0.5)/2.

                data[state + syndrome + error_sign]['x'] = data[state + 'RO0' + syndrome + error_sign]['x']
    return data

def QEC_3rounds_analyze_outcome_probabilities(data):

    prop_list = ['p0000','p0100','p1000','p1100','p0010','p0110','p1010','p1110',
                 'p0001','p0101','p1001','p1101','p0011','p0111','p1011','p1111']

    for state in ['Z','mZ']:
        for syndrome in ['11']:
            data[state + syndrome + '_outcome_probs'] = {}
            for i,RO in enumerate([0,1,2,6]):

                if i == 0:
                    for ii in range(len(prop_list)):
                        data[state + syndrome +'_outcome_probs'][prop_list[ii]] = data[state + 'RO'+str(RO) + 'S' +syndrome][prop_list[ii]]
                        data[state + syndrome +'_outcome_probs'][prop_list[ii]] = data[state + 'RO'+str(RO) + 'S' +syndrome][prop_list[ii]]
                else:
                    for ii in range(len(prop_list)):
                        data[state + syndrome +'_outcome_probs'][prop_list[ii]] = (data[state + syndrome +'_outcome_probs'][prop_list[ii]] +
                                                                         data[state + 'RO'+str(RO)+ 'S' + syndrome][prop_list[ii]])
            for ii in range(len(prop_list)):
                data[state + syndrome +'_outcome_probs'][prop_list[ii]] = data[state + syndrome + '_outcome_probs'][prop_list[ii]]/4.
    
    
    return data

def QEC_3rounds_analysis(run =1,load_from_data = False):
    ''' this functions returns a dictionairy that contains all data from the QEC_3rounds_experiments.
        TO DO: describe the dictionairy entries
        General_idea: anything not mentioned is averaged over'''

    data_dict  = QEC_3rounds_load_data(run = run, load_from_data = load_from_data)
    data_dict  = QEC_3rounds_combine_eRO(data_dict)
    data_dict  = QEC_3rounds_combine_error_signs(data_dict)
    data_dict  = QEC_3rounds_combine_syndromes(data_dict)
    data_dict  = QEC_3rounds_apply_final_QEC(data_dict)
    data_dict  = QEC_3rounds_analyze_outcome_probabilities(data_dict)

    return data_dict

def QEC_3rounds_combined_runs(runs = [1,2]):

    prop_list  = ['p0000','p0100','p1000','p1100','p0010','p0110','p1010','p1110',
                  'p0001','p0101','p1001','p1101','p0011','p0111','p1011','p1111']

    error_probs   = {}
    error_probs1  = {}
    error_probs2  = {}

    if len(runs) == 1:

        data_dict1  = QEC_3rounds_analysis(run=runs[0], load_from_data=False)
        x          = data_dict1['Z']['x']
        y_Z        = data_dict1['Z']['y']
        y_err_Z    = data_dict1['Z']['y_err']
        y_mZ       = data_dict1['mZ']['y']
        y_err_mZ   = data_dict1['mZ']['y_err']

        for ii in range(len(prop_list)):
            error_probs[prop_list[ii]] = (data_dict1['Z11'+'_outcome_probs'][prop_list[ii]] + data_dict1['mZ11'+'_outcome_probs'][prop_list[ii]])/2.

    else:
        data_dict1  = QEC_3rounds_analysis(run=runs[0], load_from_data=False)
        x1          = data_dict1['Z']['x']
        y_Z1        = data_dict1['Z']['y']
        y_err_Z1    = data_dict1['Z']['y_err']
        y_mZ1       = data_dict1['mZ']['y']
        y_err_mZ1   = data_dict1['mZ']['y_err']

        for ii in range(len(prop_list)):
            error_probs1[prop_list[ii]] = (data_dict1['Z11'+'_outcome_probs'][prop_list[ii]] +data_dict1['mZ11'+'_outcome_probs'][prop_list[ii]])/2.

        weight_Z1    =  1/(y_err_Z1)**2
        weight_mZ1   =  1/(y_err_mZ1)**2
        
        data_dict2  = QEC_3rounds_analysis(run=runs[1], load_from_data=False)
        x2          = data_dict2['Z']['x']
        y_Z2        = data_dict2['Z']['y']
        y_err_Z2    = data_dict2['Z']['y_err']
        y_mZ2       = data_dict2['mZ']['y']
        y_err_mZ2   = data_dict2['mZ']['y_err']

        for ii in range(len(prop_list)):
            error_probs2[prop_list[ii]] = (data_dict2['Z11'+'_outcome_probs'][prop_list[ii]] +data_dict2['mZ11'+'_outcome_probs'][prop_list[ii]])/2.

        weight_Z2    =  1/(y_err_Z2)**2
        weight_mZ2   =  1/(y_err_mZ2)**2

        x   = np.array([ 0.,   0.1,  0.2,  0.3,  0.4,  0.45, 0.5])

        ### Combine error probabilities
        for ii in range(len(prop_list)):
            error_probs[prop_list[ii]] = np.array([error_probs1[prop_list[ii]][0],   
                                                    error_probs1[prop_list[ii]][1],  
                                                    (error_probs1[prop_list[ii]][2]+error_probs2[prop_list[ii]][0])/2.,   
                                                    (error_probs1[prop_list[ii]][3]+error_probs2[prop_list[ii]][1])/2.,   
                                                    (error_probs1[prop_list[ii]][4]+error_probs2[prop_list[ii]][2])/2.,   
                                                    error_probs2[prop_list[ii]][3],  
                                                    error_probs1[prop_list[ii]][5]])

        ### Combine without weight
        y_Z = np.array([y_Z1[0],   y_Z1[1],  (y_Z1[2]+y_Z2[0])/2.,   (y_Z1[3]+y_Z2[1])/2.,   (y_Z1[4]+y_Z2[2])/2.,   y_Z2[3],  y_Z1[5]])
        y_mZ = np.array([y_mZ1[0], y_mZ1[1], (y_mZ1[2]+y_mZ2[0])/2., (y_mZ1[3]+y_mZ2[1])/2., (y_mZ1[4]+y_mZ2[2])/2., y_mZ2[3], y_mZ1[5]])

        y_err_Z = np.array([y_err_Z1[0],   y_err_Z1[1],  (y_err_Z1[2]**2+y_err_Z2[0]**2)**0.5/2.,   (y_err_Z1[3]**2+y_err_Z2[1]**2)**0.5/2.,   (y_err_Z1[4]**2+y_err_Z2[2]**2)**0.5/2.,   y_err_Z2[3],  y_err_Z1[5]])
        y_err_mZ = np.array([y_err_mZ1[0],   y_err_mZ1[1],  (y_err_mZ1[2]**2+y_err_mZ2[0]**2)**0.5/2.,   (y_err_mZ1[3]**2+y_err_mZ2[1]**2)**0.5/2.,   (y_err_mZ1[4]**2+y_err_mZ2[2]**2)**0.5/2.,   y_err_mZ2[3],  y_err_mZ1[5]])

        ### Combine with weight
        y_Z = np.array([y_Z1[0],   
                        y_Z1[1],  
                        np.average([y_Z1[2],y_Z2[0]], weights = [y_Z1[2], y_Z2[0]]),    
                        np.average([y_Z1[3],y_Z2[1]], weights = [y_Z1[3], y_Z2[1]]),    
                        np.average([y_Z1[4],y_Z2[2]], weights = [y_Z1[4], y_Z2[2]]),    
                        y_Z2[3],  
                        y_Z1[5]])

        y_err_Z = np.array([y_err_Z1[0],   
                            y_err_Z1[1],  
                            1/(weight_Z1[2]+weight_Z2[0])**0.5,   
                            1/(weight_Z1[3]+weight_Z2[1])**0.5,   
                            1/(weight_Z1[4]+weight_Z2[2])**0.5,   
                            y_err_Z2[3],  
                            y_err_Z1[5]])


        y_mZ = np.array([y_mZ1[0],   
                        y_mZ1[1],   
                        np.average([y_mZ1[2],y_mZ2[0]], weights = [y_mZ1[2], y_mZ2[0]]),    
                        np.average([y_mZ1[3],y_mZ2[1]], weights = [y_mZ1[3], y_mZ2[1]]),    
                        np.average([y_mZ1[4],y_mZ2[2]], weights = [y_mZ1[4], y_mZ2[2]]),    
                        y_mZ2[3],   
                        y_mZ1[5]])
        y_err_mZ = np.array([y_err_mZ1[0],   
                            y_err_mZ1[1],  
                            1/(weight_mZ1[2]+weight_mZ2[0])**0.5,   
                            1/(weight_mZ1[3]+weight_mZ2[1])**0.5,   
                            1/(weight_mZ1[4]+weight_mZ2[2])**0.5,   
                            y_err_mZ2[3],  
                            y_err_mZ1[5]])

    return x, y_Z, y_err_Z, y_mZ, y_err_mZ, error_probs

def QEC_3rounds_plot_final_curves(runs=[1,2],load_from_data = False, save_folder = r'D:\measuring\data\QEC_data\figs'):

    x, y_Z, y_err_Z, y_mZ, y_err_mZ, error_probs = QEC_3rounds_combined_runs(runs = runs)

    ### Full result

    if 0:
        fig1,ax = plt.subplots()
        ax.errorbar(x, y_Z, yerr=y_err_Z,color = 'b' )
        ax.plot([x[0],x[-1]], [y_Z[0],y_Z[-1]],'k:' )
        ax.plot([0,0.5], [1,0],'k:' )
        ax.set_ylim(-0.05,1.05)
        ax.set_xlim(-0.05,0.55)
        ax.set_title('QEC_data_3Rounds' +'\n' + 'State = Z')
        ax.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Expectation value')

        if save_folder != None:
            try:
                fig1.savefig(
                    os.path.join(save_folder,'StateZ.png'))
            except:
                print 'Figure has not been saved.'

        fig2,ax2 = plt.subplots()
        ax2.errorbar(x, y_mZ, yerr=y_err_mZ,color = 'b' )
        ax2.plot([x[0],x[-1]], [y_mZ[0],y_mZ[-1]],'k:' )
        ax2.plot([0,0.5], [-1,0],'k:' )
        ax2.set_ylim(-1.05,0.05)
        ax2.set_xlim(-0.05,0.55)
        ax2.set_title('QEC_data_3Rounds' +'\n' + 'State = -Z')
        ax2.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax2.set_xlabel('error probability')
        ax2.set_ylabel('Expectation value')

        if save_folder != None:
            try:
                fig2.savefig(
                    os.path.join(save_folder,'StatemZ.png'))
            except:
                print 'Figure has not been saved.'


    ### Averaging the states
    if 0:
        fig4,ax4 = plt.subplots()
        ax4.errorbar(x, (y_Z-y_mZ)/2., yerr=(y_err_Z**2+y_err_mZ**2)**0.5/2.,color = 'b' )
        ax4.plot([x[0],x[-1]], [(y_Z[0]-y_mZ[0])/2,(y_Z[-1]-y_mZ[-1])/2],'k:' )
        ax4.plot([0,0.5], [1,0],'k:' )
        ax4.set_ylim(-0.05,1.05)
        ax4.set_xlim(-0.05,0.55)
        ax4.set_title('QEC_data_3Rounds Zstate' + '\n' + 'State = (Z+mZ)/2')
        ax4.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax4.set_xlabel('error probability')
        ax4.set_ylabel('Expectation value')

        if save_folder != None:
            try:
                fig2.savefig(
                    os.path.join(save_folder,'Analyze_script.png'))
            except:
                print 'Figure has not been saved.'

    if 1:
        x1, y_Z1, y_err_Z1, y_mZ1, y_err_mZ1, error_probs =  QEC_3rounds_combined_runs(runs=[1])
        x2, y_Z2, y_err_Z2, y_mZ2, y_err_mZ2, error_probs =  QEC_3rounds_combined_runs(runs=[2])

        fig4,ax4 = plt.subplots()
        ax4.errorbar(x1, (y_Z1-y_mZ1)/2., yerr=(y_err_Z1**2+y_err_mZ1**2)**0.5/2.,color = 'b', label = 'run1')
        ax4.errorbar(x2, (y_Z2-y_mZ2)/2., yerr=(y_err_Z2**2+y_err_mZ2**2)**0.5/2.,color = 'r', label = 'run2' )

        ax4.plot([x[0],x[-1]], [(y_Z[0]-y_mZ[0])/2,(y_Z[-1]-y_mZ[-1])/2],'k:' )
        ax4.plot([0,0.5], [1,0],'k:' )
        ax4.set_ylim(-0.05,1.05)
        ax4.set_xlim(-0.05,0.55)
        ax4.set_title('QEC_data_3Rounds Zstate' +  '\n' + 'State = (Z+mZ)/2')
        ax4.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
        ax4.set_xlabel('error probability')
        ax4.set_ylabel('Expectation value')
        ax4.legend()

        if save_folder != None:
            try:
                fig4.savefig(
                    os.path.join(save_folder,'3Rounds_compare_runs.pdf'))
            except:
                print 'Figure has not been saved.'

    plt.show()
    plt.close('all')

def QEC_3rounds_plot_outcome_probability_curves(runs = [1,2], load_from_data = False, save_folder = r'D:\measuring\data\QEC_data\figs'):

    x, y_Z, y_err_Z, y_mZ, y_err_mZ, error_probs = QEC_3rounds_combined_runs(runs = runs)

    prop_list  = ['p0000','p0100','p1000','p1100','p0010','p0110','p1010','p1110',
                  'p0001','p0101','p1001','p1101','p0011','p0111','p1011','p1111']


    #####################################################
    ### A. Probabilties for the 16 different outcomes ###
    #####################################################
    
    fig1,ax = plt.subplots()
    colors  = cm.rainbow(np.linspace(0, 1, len(prop_list)))

    for ii in range(len(prop_list)):
            p = error_probs[prop_list[ii]]
            ax.plot(x,p,'o',ls=':',color = colors[ii], label = prop_list[ii])

    ax.set_title('QEC_data_3Rounds outcome probabilities' + '\n' + 'syndrome 11')
    ax.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('outcome probability')     
    ax.set_ylim(-0.05,.7)
    ax.set_xlim(-0.05,0.55+0.2)
    ax.legend( bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.legend(loc=1)
   
    try:
        fig1.savefig(
            os.path.join(save_folder,'Error_probs_3rounds'+'.png'))
    except:
        print 'Figure has not been saved.'

    ###################################
    ### probabilities per QEC round ###
    ###################################

    ### Round 1 ###
    p_round1_00 = (error_probs['p0000'] +
                   error_probs['p0010'] +
                   error_probs['p0001'] +
                   error_probs['p0011'])
    p_round1_01 = (error_probs['p0100'] +
                   error_probs['p0110'] +
                   error_probs['p0101'] +
                   error_probs['p0111'])
    p_round1_10 = (error_probs['p1000'] +
                   error_probs['p1010'] +
                   error_probs['p1001'] +
                   error_probs['p1011'])
    p_round1_11 = (error_probs['p1100'] +
                   error_probs['p1110'] +
                   error_probs['p1101'] +
                   error_probs['p1111'])

    fig6,ax = plt.subplots()
    ax.plot(x,p_round1_00,'bo:', label = 'R1 00')
    ax.plot(x,p_round1_01,'ko:', label = 'R1 01')
    ax.plot(x,p_round1_10,'go:', label = 'R1 10')
    ax.plot(x,p_round1_11,'ro:', label = 'R1 11')

    ax.set_title('outcome probabilities round 1')
    ax.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('outcome probability')
    ax.legend( bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.legend(loc=1)
    ax.set_ylim(-0.05,.8)
    ax.set_xlim(-0.05,0.55)

    ### Round 2 ###
    p_round2_00 = (error_probs['p0000'] +
                   error_probs['p1000'] +
                   error_probs['p0100'] +
                   error_probs['p1100'])
    p_round2_01 = (error_probs['p0001'] +
                   error_probs['p1001'] +
                   error_probs['p0101'] +
                   error_probs['p1101'])
    p_round2_10 = (error_probs['p0010'] +
                   error_probs['p1010'] +
                   error_probs['p0110'] +
                   error_probs['p1110'])
    p_round2_11 = (error_probs['p0011'] +
                   error_probs['p1011'] +
                   error_probs['p0111'] +
                   error_probs['p1111'])

    # fig6,ax = plt.subplots()
    ax.plot(x,p_round2_00,'bo-', mfc='none', label = 'R2 00')
    ax.plot(x,p_round2_01,'ko-', mfc='none', label = 'R2 01')
    ax.plot(x,p_round2_10,'go-', mfc='none', label = 'R2 10')
    ax.plot(x,p_round2_11,'ro-', mfc='none', label = 'R2 11')

    ax.set_title('outcome probabilities for separate rounds')
    ax.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('outcome probability')
    ax.legend( bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.legend( loc=1)
    ax.set_ylim(-0.05,.8)
    ax.set_xlim(-0.05,0.55+0.1)

    try:
        fig6.savefig(
            os.path.join(save_folder,'Error_probs_per_round_3rounds'+'.png'))
    except:
        print 'Figure has not been saved.'

    ####################################
    ### post selected probabilities? ###
    ####################################

    plt.show()
    plt.close('all')

    return x, p_round1_11, p_round2_11

def QEC_3rounds_outcome_probability(runs = [1,2], load_from_data = False, save_folder = r'D:\measuring\data\QEC_data\figs'):

    x, y_Z, y_err_Z, y_mZ, y_err_mZ, error_probs = QEC_3rounds_combined_runs(runs = runs)

    prop_list  = ['p0000','p0100','p1000','p1100','p0010','p0110','p1010','p1110',
                  'p0001','p0101','p1001','p1101','p0011','p0111','p1011','p1111']



    ###################################
    ### probabilities per QEC round ###
    ###################################

    ### Round 1 ###
    p_round1_00 = (error_probs['p0000'] +
                   error_probs['p0010'] +
                   error_probs['p0001'] +
                   error_probs['p0011'])
    p_round1_01 = (error_probs['p0100'] +
                   error_probs['p0110'] +
                   error_probs['p0101'] +
                   error_probs['p0111'])
    p_round1_10 = (error_probs['p1000'] +
                   error_probs['p1010'] +
                   error_probs['p1001'] +
                   error_probs['p1011'])
    p_round1_11 = (error_probs['p1100'] +
                   error_probs['p1110'] +
                   error_probs['p1101'] +
                   error_probs['p1111'])


    ### Round 2 ###
    p_round2_00 = (error_probs['p0000'] +
                   error_probs['p1000'] +
                   error_probs['p0100'] +
                   error_probs['p1100'])
    p_round2_01 = (error_probs['p0001'] +
                   error_probs['p1001'] +
                   error_probs['p0101'] +
                   error_probs['p1101'])
    p_round2_10 = (error_probs['p0010'] +
                   error_probs['p1010'] +
                   error_probs['p0110'] +
                   error_probs['p1110'])
    p_round2_11 = (error_probs['p0011'] +
                   error_probs['p1011'] +
                   error_probs['p0111'] +
                   error_probs['p1111'])

    return x, p_round1_11, p_round2_11

def QEC_3rounds_plot_outcome_probability_curves_old(load_from_data = False, save_folder = r'D:\measuring\data\QEC_data\figs'):

    syndrome = '11'

    data_dict = QEC_3rounds_analysis(run=1, load_from_data=False)
    x         = data_dict['Z']['x']

    data_dict2 = QEC_3rounds_analysis(run=2, load_from_data=False)
    x2         = data_dict2['Z']['x']

    prop_list  = ['p0000','p0100','p1000','p1100','p0010','p0110','p1010','p1110',
                  'p0001','p0101','p1001','p1101','p0011','p0111','p1011','p1111']


    ### Probabilties for the 16 different outcomes
    fig1,ax = plt.subplots()
    colors  = cm.rainbow(np.linspace(0, 1, len(prop_list)))

    for ii in range(len(prop_list)):
            p = (data_dict['Z' +syndrome + '_outcome_probs'][prop_list[ii]] + data_dict['mZ' +syndrome + '_outcome_probs'][prop_list[ii]])/2
            ax.plot(x,p,'o',ls='',color = colors[ii], label = 'run1_' + prop_list[ii])

    for ii in range(len(prop_list)):
            p = (data_dict2['Z'+syndrome + '_outcome_probs'][prop_list[ii]]+data_dict2['Z'+syndrome + '_outcome_probs'][prop_list[ii]])/2
            ax.plot(x2,p,'b', color = colors[ii], label = 'run_2' + prop_list[ii])
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    ax.set_title('QEC_data_3Roundsoutcome probabilities' + '\n' + 'State = Z, syndrome ' + syndrome)
    ax.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('outcome probability')
    ax.legend( bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_ylim(-0.05,.7)
    ax.set_xlim(-0.05,0.55)

    try:
        fig1.savefig(
            os.path.join(save_folder,'Error_probs_Z_' + syndrome +'.pdf'))
    except:
        print 'Figure has not been saved.'

    fig1,ax = plt.subplots()

    colors = cm.rainbow(np.linspace(0, 1, len(prop_list)))

    for ii in range(len(prop_list)):
        p = data_dict['mZ' + syndrome+'_outcome_probs'][prop_list[ii]]
        ax.plot(x,p,'b', color = colors[ii], label = prop_list[ii])

    ax.set_title('QEC_data_3Rounds outcome probabilities' + '\n' + 'State = mZ, syndrome '+syndrome)
    ax.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('outcome probability')
    ax.legend( bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_ylim(-0.05,.7)
    ax.set_xlim(-0.05,0.55)

    try:
        fig1.savefig(
            os.path.join(save_folder,'Error_probs_Z_' + syndrome +'.pdf'))
    except:
        print 'Figure has not been saved.'

    ###############################
    ### probabilities per round ###
    ###############################

    state   = 'Z'
    syndrome = '11'

    ### Round 1 ###
    p_round1_00 = (data_dict[state + syndrome+'_outcome_probs']['p0000'] +
                   data_dict[state + syndrome+'_outcome_probs']['p0010'] +
                   data_dict[state + syndrome+'_outcome_probs']['p0001'] +
                   data_dict[state + syndrome+'_outcome_probs']['p0011'])
    p_round1_01 = (data_dict[state + syndrome+'_outcome_probs']['p0100'] +
                   data_dict[state + syndrome+'_outcome_probs']['p0110'] +
                   data_dict[state + syndrome+'_outcome_probs']['p0101'] +
                   data_dict[state + syndrome+'_outcome_probs']['p0111'])
    p_round1_10 = (data_dict[state + syndrome+'_outcome_probs']['p1000'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1010'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1001'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1011'])
    p_round1_11 = (data_dict[state + syndrome+'_outcome_probs']['p1100'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1110'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1101'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1111'])

    fig6,ax = plt.subplots()
    ax.plot(x,p_round1_00,'b', label = '00')
    ax.plot(x,p_round1_01,'k', label = '01')
    ax.plot(x,p_round1_10,'g', label = '10')
    ax.plot(x,p_round1_11,'r', label = '11')

    ax.set_title('outcome probabilities round 1' + '\n' + 'State = ' +state+ ', syndrome' +syndrome)
    ax.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('outcome probability')
    ax.legend( bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_ylim(-0.05,.8)
    ax.set_xlim(-0.05,0.55)

    ### Round 2 ###
    p_round2_00 = (data_dict[state + syndrome+'_outcome_probs']['p0000'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1000'] +
                   data_dict[state + syndrome+'_outcome_probs']['p0100'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1100'])
    p_round2_01 = (data_dict[state + syndrome+'_outcome_probs']['p0001'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1001'] +
                   data_dict[state + syndrome+'_outcome_probs']['p0101'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1101'])
    p_round2_10 = (data_dict[state + syndrome+'_outcome_probs']['p0010'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1010'] +
                   data_dict[state + syndrome+'_outcome_probs']['p0110'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1110'])
    p_round2_11 = (data_dict[state + syndrome+'_outcome_probs']['p0011'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1011'] +
                   data_dict[state + syndrome+'_outcome_probs']['p0111'] +
                   data_dict[state + syndrome+'_outcome_probs']['p1111'])

    fig6,ax = plt.subplots()
    ax.plot(x,p_round2_00,'b', label = '00')
    ax.plot(x,p_round2_01,'k', label = '01')
    ax.plot(x,p_round2_10,'g', label = '10')
    ax.plot(x,p_round2_11,'r', label = '11')

    ax.set_title('outcome probabilities round 2' + '\n' + 'State = ' +state+ ', syndrome' +syndrome)
    ax.hlines([-1,0,1],x[0]-0.05,x[-1]+0.05,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('outcome probability')
    ax.legend( bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_ylim(-0.05,.8)
    ax.set_xlim(-0.05,0.55)

    ###################################
    ### post selected probabilities ###
    ###################################


    plt.show()
    plt.close('all')


###############
''' Figures '''
###############

########################################################
''' Combined plots, multiple rounds (0,1,2,3 rounds) '''
########################################################

def QEC_multiple_rounds_plot_combined_curves(save_folder = r'D:\measuring\data\QEC_data\figs'):

    ### load 0 round data
    qubit = 2
    single_no_QEC_data_dict_Z_Q1  =  single_qubit_avg_state(state = 'Z',run = 1)
    single_no_QEC_data_dict_mZ_Q1 =  single_qubit_avg_state(state = 'mZ',run = 1)

    y_Z0         = single_no_QEC_data_dict_Z_Q1['y']
    y_err_Z0     = single_no_QEC_data_dict_Z_Q1['y_err']
    y_mZ0        = single_no_QEC_data_dict_mZ_Q1['y']
    y_err_mZ0    = single_no_QEC_data_dict_mZ_Q1['y_err']
    x0           = single_no_QEC_data_dict_mZ_Q1['x']

    ### load 1 round data (Majority voting/Toffoli)
    toff_process_dict = no_QEC_toffoli_fids(run_list = [0,1,2,3],state_list=['Z','mZ'],add_4 = True,do_weighted = True)
    y_Z1              = toff_process_dict['toff_Zy']
    y_err_Z1          = toff_process_dict['toff_Zy_err']
    y_mZ1             = toff_process_dict['toff_mZy']
    y_err_mZ1         = toff_process_dict['toff_mZy_err']
    x1                = toff_process_dict['x']

    # load 1 round data (QEC)
    data_dict_Z       = QEC_state_sum_all(state = 'Z', RO = qubit-1)
    data_dict_mZ      = QEC_state_sum_all(state = 'mZ', RO = qubit-1)

    y_Z1b         = data_dict_Z['y']
    y_err_Z1b     = data_dict_Z['y_err']
    y_mZ1b        = data_dict_mZ['y']
    y_err_mZ1b    = data_dict_mZ['y_err']
    x1b           = data_dict_mZ['x']

    ### load 2 rounds data

    x2, y_Z2, y_err_Z2, y_mZ2, y_err_mZ2, error_probs =  QEC_2rounds_combined_runs(runs=[1,2,3])

    ### load 3 rounds data
    x3, y_Z3, y_err_Z3, y_mZ3, y_err_mZ3, error_probs =  QEC_3rounds_combined_runs(runs=[1,2])

    ### Averaging over Z and mZ
    y_tot0  = (y_Z0-y_mZ0)/2; y_err_tot0 = (y_err_Z0**2+y_err_mZ0**2)**0.5/2
    y_tot1  = (y_Z1-y_mZ1)/2; y_err_tot1 = (y_err_Z1**2+y_err_mZ1**2)**0.5/2
    y_tot1b = (y_Z1b-y_mZ1b)/2; y_err_tot1b = (y_err_Z1b**2+y_err_mZ1b**2)**0.5/2
    y_tot2  = (y_Z2-y_mZ2)/2; y_err_tot2 = (y_err_Z2**2+y_err_mZ2**2)**0.5/2
    y_tot3  = (y_Z3-y_mZ3)/2; y_err_tot3 = (y_err_Z3**2+y_err_mZ3**2)**0.5/2

    print type(y_tot1)
    print type(y_tot2)

    ### Fitting
    x_fit0, y_fit0, fit_result0,u_fit_result0 = fit_QEC_curve(x0,y_tot0, return_errorbar=True)
    x_fit1, y_fit1, fit_result1,u_fit_result1 = fit_QEC_curve(x1,y_tot1, return_errorbar=True)
    x_fit2, y_fit2, fit_result2,u_fit_result2 = fit_QEC_2_rounds_curve2(x2,y_tot2, return_errorbar=True)
    x_fit3, y_fit3, fit_result3,u_fit_result3 = fit_QEC_3_rounds_curve2(x3,y_tot3, return_errorbar=True)
    # x_fit2, y_fit2 = fit_QEC_2_rounds_curve(x2,y_tot2, return_errorbar=False,plot_guess = True)
    # x_fit3, y_fit3 = fit_QEC_3_rounds_curve2(x3,y_tot3, return_errorbar=False,plot_guess = True)
    ### Complete result (avaraged over Z and mZ)
    fig4,ax = plt.subplots(figsize=(8,8))

    ax.errorbar(x0, (0.5*y_tot0+0.5), yerr= 0.5*y_err_tot0, color = 'm', ls = '', marker = 'o', ms = 4, label = '0 rounds, $p_c$='+str(int(fit_result0*1000)/1000.)+'('+str(int(u_fit_result0*1000))+')')
    ax.plot(x_fit0, (0.5*y_fit0+0.5), color = 'm')

    ax.errorbar(x1, (0.5*y_tot1+0.5), yerr= 0.5*y_err_tot1,color = 'g', ls = '', marker = 'o', ms = 4, label = '1 round, $p_c$='+str(int(fit_result1*1000)/1000.)+'('+str(int(u_fit_result1*1000))+')')
    ax.plot(x_fit1, (0.5*y_fit1+0.5), color = 'g')

    # ax.errorbar(x1b, (0.5*y_tot1n+0.5), yerr= 0.5*y_err_tot1b,color = 'm', label = '1 rounds QEC' )

    ax.errorbar(x2, (0.5*y_tot2+0.5), yerr=0.5*y_err_tot2,color = 'b', marker = 'o', ms = 4, ls = '',label = '2 rounds, $p_c$='+str(int(fit_result2*1000)/1000.)+'('+str(int(u_fit_result2*1000))+')')
    ax.plot(x_fit2, (0.5*y_fit2+0.5), color = 'b')

    ax.errorbar(x3, (0.5*y_tot3+0.5), yerr=0.5*y_err_tot3,color = 'r', marker = 'o', ms = 4,ls = '', label = '3 rounds, $p_c$='+str(int(fit_result3*1000)/1000.)+'('+str(int(u_fit_result3*1000))+')')
    ax.plot(x_fit3, (0.5*y_fit3+0.5), color = 'r')

    # ax.plot([x[0],x[-1]], [(y_Z[0]-y_mZ[0])/2,(y_Z[-1]-y_mZ[-1])/2],'k:' )

    ax.plot([0,0.5], [1,0.5],'k:' )
    ax.set_ylim(0.45,1.05)
    ax.set_xlim(-0.05,0.55)
    ax.set_title('QEC_data_multiple rounds, averaged over Z and -Z')
    ax.hlines([-1,0,1],x1[0]-0.05,x0[-1]+0.05,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Average fidelity')
    ax.legend()


    if save_folder != None:
        try:
            fig4.savefig(
                os.path.join(save_folder,'Multiple_rounds_Combined11.pdf'))
            fig4.savefig(
                os.path.join(save_folder,'Multiple_rounds_Combined11.png'))
        except:
            print 'Figure has not been saved.'

    plt.show()
    # plt.close('all')

#######################################
''' Error probability curves fitted '''
#######################################
folder = r'D:\measuring\data\QEC_data\figs\final figures'

def QEC_fit_error_probability_curves():

    fig,ax = plt.subplots(figsize = (10,10))
    syndrome_list = ['00','01','10','11']
    p_dict = {}
    for syndrome in syndrome_list:

        if syndrome == '00' or syndrome == '01':
            run_list = [1,2,3]
        elif syndrome == '10':
            run_list = [2]
        elif syndrome == '11':
            run_list = [3]

        p_list = ['p00','p01','p10','p11']
        p_dict[syndrome] = QEC_sum_probs(run_list = run_list ,no_error = syndrome)

    p_dict['no error'] = 1/4.*(p_dict['00']['p00']+p_dict['01']['p01']+p_dict['10']['p10']+p_dict['11']['p11'])
    p_dict['Q3'] = 1/4.*(p_dict['00']['p01']+p_dict['01']['p00']+p_dict['10']['p11']+p_dict['11']['p10'])  # Carbon 5
    p_dict['Q1'] = 1/4.*(p_dict['00']['p10']+p_dict['01']['p11']+p_dict['10']['p00']+p_dict['11']['p01'])  # Carbon 2
    p_dict['Q2'] = 1/4.*(p_dict['00']['p11']+p_dict['01']['p10']+p_dict['10']['p01']+p_dict['11']['p00'])  # Carbon 1


    color = [c_orange,c_green,c_red,'r']
    pin = [0.127,0.123,0.145,0.113]
    for i,error in enumerate(['no error','Q1','Q2','Q3']):
        if error == 'no error':
            # x_fit,y_fit, pin, u_pin = fit_no_error_detection_curve(p_dict[syndrome]['x'],p_dict[error], return_errorbar=True,plot_guess = False)
            x_fit,y_fit= fit_no_error_detection_curve(linspace(0,1,100),0.25*np.ones(100), return_errorbar=False,plot_guess = True)

        else:
            # x_fit,y_fit, pin, u_pin = fit_error_detection_curve(p_dict[syndrome]['x'],p_dict[error], return_errorbar=True,plot_guess = False)
            x_fit,y_fit= fit_error_detection_curve(linspace(0,1,100),0.25*np.ones(100),pin = pin[i], return_errorbar=False,plot_guess = True)
        ax.plot(p_dict[syndrome]['x'],p_dict[error],label = error, color = color[i], marker = 'o',markersize = 5,ls = '',markeredgecolor = color[i])

        ax.plot(x_fit,y_fit,color = color[i])
    lgd = ax.legend(loc = 9,frameon=False)
    for label in lgd.get_texts():
        label.set_fontsize(25)

    plt.xlim([-0.01,1.01])
    plt.xlabel('$p_e$',fontsize = 25)
    plt.ylim([-0.01,1.01])
    plt.xticks([0,0.5,1])
    plt.yticks([0,0.5,1])
    plt.ylabel('Normalized \n occurence',fontsize = 25)

    plt.xticks(np.arange(0,1.1,0.5))
    # plt.xticks(np.arange(0,1.1,0.1), minor = True)
    plt.yticks(np.arange(0,1.1,0.5))
    # plt.yticks(np.arange(0,1.1,0.1), minor = True)
    plt.tick_params(axis='x', which='major', labelsize=25)
    plt.tick_params(axis='y', which='major', labelsize=25)
    plt.tick_params('both', length=6, width=1, which='major')
    plt.tick_params('both', length=4, width=1, which='minor')

def plot_prob_single_syndrome(syndrome = '11',run_list = [],add_simulation=True):
    if syndrome == '00' or syndrome == '01':
        run_list = [1,2,3]
    elif syndrome == '10':
        run_list = [2]
    elif syndrome == '11':
        if run_list == []:
            run_list = [3]

    if syndrome == '00':
        p_list = ['p00','p10','p11','p01']

    if syndrome == '11':
        p_list = ['p11','p01','p00','p10']

    if syndrome == '01':
        p_list = ['p01','p11','p10','p00']

    if syndrome == '10':
        p_list = ['p10','p00','p01','p11']

    p_dict = QEC_sum_probs(run_list = run_list ,no_error = syndrome)

    if add_simulation == True:
        pin_no_alt= 0.093 # 0 #
        # pin_c1 =0.072 # 0 #
        # pin_c2 = 0.093 # 0 #
        # pin_c5 = 0.115 # 0 #

        # pin_no =0.035# 0.093 # 0 #
        # pin_c1 =0.035#0.072 # 0 #
        # pin_c2 =0.035# 0.093 # 0 #
        # pin_c5 =0.035# 0.115 # 0 #

        # # Obtained from encoding Corrected for el RO and basis rotations
        # pin_no =0.062# 0.093 # 0 #
        # pin_c1 =0.0766#0.072 # 0 #
        # pin_c2 =0.0483# 0.093 # 0 #
        # pin_c5 =0.0614# 0.115 # 0 #

        # # Obtained from encoding Corrected for el RO 
        # pin_no =0.0776# 0.093 # 0 #
        # pin_c1 =0.0917#0.072 # 0 #
        # pin_c2 =0.0645# 0.093 # 0 #
        # pin_c5 =0.0770# 0.115 # 0 #


        # pin_no      =0.0856612#0.062# 0.093 # 0 #
        # pin_no_alt  = 0.093
        # pin_c1      =0.0650188#0.0766#0.072 # 0 #
        # pin_c2      =0.0822674#0.0483# 0.093 # 0 #
        # pin_c5      =0.11562#0.0614# 0.115 # 0 #


        # obtained from corrected way 150428
        pin_no = 0.0775644
        pin_c1 = 0.0640236
        pin_c2 = 0.0912455
        pin_c5 = 0.077205

        F1 = 0.988
        F0 = 0.890

        p = linspace(0,1,1000)

        ptot_no = p + pin_no-2*p*pin_no
        # ptot_no_alt = p + pin_no_alt-2*p*pin_no_alt
        ptot_c1 = p + pin_c1-2*p*pin_c1
        ptot_c2 = p + pin_c2-2*p*pin_c2
        ptot_c5 = p + pin_c5-2*p*pin_c5

        # p_no_error = 1-3*ptot_no+3*ptot_no**2

        # p_c1       = ptot_c1-ptot_c1**2
        # p_c2       = ptot_c2-ptot_c2**2
        # p_c5       = ptot_c5-ptot_c5**2

        #correct way to do this 150428
        p_no_error = 1-ptot_c1-ptot_c2-ptot_c5+ptot_c1*ptot_c2+ptot_c1*ptot_c5+ptot_c2*ptot_c5
        p_c1 = ptot_c1-ptot_c1*ptot_c2-ptot_c1*ptot_c5+ptot_c2*ptot_c5
        p_c2 = ptot_c2-ptot_c2*ptot_c1-ptot_c2*ptot_c5+ptot_c1*ptot_c5
        p_c5 = ptot_c5-ptot_c5*ptot_c1-ptot_c5*ptot_c2+ptot_c1*ptot_c2

        #### For 11 assignment
        if syndrome == '11':
            P_D_no_error = p_no_error*F1**2 + p_c1*(1-F0)**2 + (p_c2+p_c5)*F1*(1-F0)
            # P_D_no_error_alt = p_no_error_alt*F1**2 + p_c1*(1-F0)**2 + (p_c2+p_c5)*F1*(1-F0)
            P_D_c1 = p_no_error*(1-F1)**2 + p_c1*F0**2 + (p_c2+p_c5)*F0*(1-F1)
            P_D_c2 = p_no_error*F1*(1-F1) + p_c1*(1-F0)*F0 + (p_c2)*F1*F0+ (p_c5)*(1-F1)*(1-F0)
            P_D_c5 = p_no_error*F1*(1-F1) + p_c1*(1-F0)*F0 + (p_c5)*F1*F0+ (p_c2)*(1-F1)*(1-F0)

        #### For 00 assignment
        if syndrome == '00':
            P_D_no_error = p_no_error*F0**2 + p_c1*(1-F1)**2 + (p_c2+p_c5)*F0*(1-F1)
            # P_D_no_error_alt = p_no_error_alt*F0**2 + p_c1*(1-F1)**2 + (p_c2+p_c5)*F0*(1-F1)
            P_D_c1 = p_no_error*(1-F0)**2 + p_c1*F1**2 + (p_c2+p_c5)*F1*(1-F0)
            P_D_c2 = p_no_error*F0*(1-F0) + p_c1*(1-F1)*F1 + (p_c2)*F0*F1+ (p_c5)*(1-F0)*(1-F1)
            P_D_c5 = p_no_error*F0*(1-F0) + p_c1*(1-F1)*F1 + (p_c5)*F0*F1+ (p_c2)*(1-F0)*(1-F1)

        #### For 01 assignment
        if syndrome == '01':
            P_D_no_error = p_no_error*F0*F1 + p_c1*(1-F0)*(1-F1) +  p_c2*F1*(1-F0)+ p_c2*F0*(1-F1)
            # P_D_no_error_alt = p_no_error_alt*F0*F1 + p_c1*(1-F0)*(1-F1) +  p_c2*F1*(1-F0)+ p_c2*F0*(1-F1)
            P_D_c1 = p_no_error*(1-F0)*(1-F1)+ p_c1*F1*F0 + p_c2*F0*(1-F1)+ p_c5*F1*(1-F0)
            P_D_c2 = p_no_error*F1*(1-F0) + p_c1*(1-F0)*F1 + (p_c2)*F1*F1+ (p_c5)*(1-F0)*(1-F1)
            P_D_c5 = p_no_error*F0*(1-F1) + p_c1*(1-F1)*F0 + (p_c5)*F0*F0+ (p_c2)*(1-F0)*(1-F1)

        #### For 10 assignment
        if syndrome == '10':
            P_D_no_error = p_no_error*F0*F1 + p_c1*(1-F0)*(1-F1) +  p_c2*F1*(1-F0)+ p_c2*F0*(1-F1)
            # P_D_no_error_alt = p_no_error_alt*F0*F1 + p_c1*(1-F0)*(1-F1) +  p_c2*F1*(1-F0)+ p_c2*F0*(1-F1)
            P_D_c1 = p_no_error*(1-F0)*(1-F1)+ p_c1*F1*F0 + p_c2*F0*(1-F1)+ p_c5*F1*(1-F0)
            P_D_c2 = p_no_error*F0*(1-F1) + p_c1*(1-F1)*F0 + (p_c2)*F0*F0+ (p_c5)*(1-F0)*(1-F1)
            P_D_c5 = p_no_error*F1*(1-F0) + p_c1*(1-F0)*F1 + (p_c5)*F1*F1+ (p_c2)*(1-F0)*(1-F1)
        

    fig,ax = plt.subplots(figsize = (5,5))
    color = [c_green,c_orange,c_red,'r']
    label = ['No error','Qubit 1','Qubit 2', 'Qubit 3']
    if add_simulation == True:
        plt.plot(p,P_D_no_error, color = color[0])
        # plt.plot(p,P_D_no_error_alt, color = color[0],ls = ':')
        plt.plot(p,P_D_c1, color = color[2])
        plt.plot(p,P_D_c2, color = color[1])
        plt.plot(p,P_D_c5, color = color[3])

        for p_D in [P_D_no_error,P_D_c1,P_D_c2,P_D_c5]:
            print p_D[len(p_D)/2.]

    for i,p in enumerate(p_list):
        ax.plot(p_dict['x'],p_dict[p],label = label[i], color = color[i], marker = 'o',
            markersize = 5,ls = '',markeredgecolor = color[i])


        # ax.plot(x_fit,y_fit,color = color[i])
    lgd = ax.legend(loc = 9,frameon=False)
    for label in lgd.get_texts():
        label.set_fontsize(25)

    plt.xlim([-0.01,1.01])
    plt.xlabel('$p_e$',fontsize = 25)
    plt.ylim([-0.0,1.0])
    plt.xticks([0,0.5,1])
    plt.yticks([0,0.5,1])
    plt.ylabel('Occurence',fontsize = 25)

    plt.xticks(np.arange(0,1.1,0.5))
    # plt.xticks(np.arange(0,1.1,0.1), minor = True)
    plt.yticks(np.arange(0,1.1,0.5))
    # plt.yticks(np.arange(0,1.1,0.1), minor = True)
    plt.tick_params(axis='x', which='major', labelsize=25)
    plt.tick_params(axis='y', which='major', labelsize=25)
    plt.tick_params('both', length=6, width=1, which='major')
    plt.tick_params('both', length=4, width=1, which='minor')

    if syndrome == '00': 
        color = c_green
    if syndrome == '01': 
        color = c_blue
    if syndrome == '10': 
        color = c_red
    if syndrome == '11': 
        color = c_orange

    # plt.setp(ax.spines.values(), color=color)
    # plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=color)

    try:
        fig.savefig(
            os.path.join(folder,'Probability_plot_syn'+syndrome+'.pdf'))
    except:
        print 'Figure has not been saved.' 

def QEC_compare_syndromes():
    fig,ax = plt.subplots(figsize= (10,10))
    mpl.rcParams['pdf.fonttype'] = 42
    
    color = [c_green, c_blue,c_red,c_orange]

    for i, no_error in enumerate(['11']):#,'01','10','11']):
        if no_error == '00':
            run_list = [1,2,3]
        elif no_error == '01':
            run_list = [1,2,3]
        elif no_error == '10':
            run_list = [2]
        elif no_error == '11':
            run_list = [3]
        process_dict = QEC_process_fids_sum_runs(run_list = run_list,no_error = no_error)


        x = process_dict['x']

        y = process_dict['dec_'+'avg'+'_y']
        y_new = process_dict['dec_'+'avg'+'_y_new']
        y_err = process_dict['dec_'+'avg'+'_y_err']
        
        if no_error == '11':
            # x_fit1, y_fit1, p_c, p_c_err= fit_QEC_process_curve_11(x,y,return_errorbar = True)
            # # x_fit1, y_fit1 = fit_QEC_process_curve_11(x,y,A=A,pc=p,O=O,return_guess=True)
            # # print p, A,O
            # ax.plot(x_fit1, y_fit1, color = color[i], lw=1, label =  'QEC, no error: '+ no_error)
            
            x_fit, y_fit, p_c, p_c_err= fit_QEC_process_curve_11(x,y_new,return_errorbar = True)
            ax.plot(x_fit, y_fit, color = color[i], lw=1, label =  'QEC, undo, no error: '+ no_error)
            print 'CORRECTION PROBABILITY UNDO'+no_error
            print p_c
            print p_c_err
            x_fit, y_fit, p_c, p_c_err= fit_QEC_process_curve_11(x,y,return_errorbar = True)

            # x_fit, y_fit = fit_QEC_process_curve(x,y,A=A,pc=p,O=O,return_guess=True)
            # print p, A,O
            ax.plot(x_fit, y_fit, color = color[i], lw=1, label =  'QEC, no error: '+ no_error)
            # ax.plot(x_fit,y_fit1-y_fit)
            # print y_fit1-y_fit
        if no_error == '00':
            # x_fit1, y_fit1, p_c, p_c_err = fit_QEC_process_curve_00(x,y,return_errorbar = True)
            # ax.plot(x_fit1, y_fit1, color = color[i], lw=1, label =  'QEC, no error: '+ no_error)
            x_fit, y_fit, p_c, p_c_err= fit_QEC_process_curve_00(x,y_new,return_errorbar = True)
            ax.plot(x_fit, y_fit, color = color[i], lw=1, label =  'QEC, undo, no error: '+ no_error)
            print 'CORRECTION PROBABILITY UNDO'+no_error
            print p_c
            print p_c_err
            x_fit, y_fit, p_c, p_c_err= fit_QEC_process_curve_00(x,y,return_errorbar = True)

            # x_fit, y_fit = fit_QEC_process_curve(x,y,A=A,pc=p,O=O,return_guess=True)
            # print p, A,O
            ax.plot(x_fit, y_fit, color = color[i], lw=1, label =  'QEC, no error: '+ no_error)
            # ax.plot(x_fit,y_fit1-y_fit)
        if no_error == '10' or no_error == '01':
            # x_fit1, y_fit1, p_c, p_c_err = fit_QEC_process_curve_01(x,y,return_errorbar = True)
            # ax.plot(x_fit1, y_fit1, color = color[i], lw=1, label =  'QEC, no error: '+ no_error)
            x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
            ax.plot(x_fit, y_fit, color = color[i], lw=1, label =  'QEC, no error: '+ no_error)
            # ax.plot(x_fit,y_fit1-y_fit)
            

        print 'CORRECTION PROBABILITY '+no_error
        print p_c
        print p_c_err

        # ax.plot(x_fit, y_fit, color = color[i], lw=1, label =  'QEC, no error: '+ no_error)#, $p_c$='+str(round(p_c*100)/100.)+'('+str(int(round(p_c_err*100)))+')')
        (_,caps,_)=ax.errorbar(x,y,yerr=y_err,color = color[i],markeredgecolor = color[i],ls = '',marker = 'o', ms = 7,capsize = 6)
        for cap in caps:
            cap.set_markeredgewidth(1)
        # y = process_dict['dec_'+'avg'+'_y_new']
        # y_err = process_dict['dec_'+'avg'+'_y_err']
        # x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
        # ax.plot(x_fit, y_fit, color = color[i],ls = '-.', lw=1)#, $p_c$='+str(round(p_c*100)/100.)+'('+str(int(round(p_c_err*100)))+')')
        # (_,caps,_)=ax.errorbar(x,y,yerr=y_err,color = color[i],markeredgecolor = color[i], ls = '',marker = 'o', ms = 7,capsize = 6)
        # for cap in caps:
        #     cap.set_markeredgewidth(1)
    

    ax.set_ylim(-0,1)
    ax.set_xlim(-0.01,1.01)

    ax.set_xticks(np.arange(0,1.1,0.5))
    ax.set_xticks(np.arange(0,1.1,0.1), minor = True)
    ax.set_yticks(np.arange(0,1.1,0.5))
    ax.set_yticks(np.arange(0,1.1,0.1), minor = True)
    ax.tick_params(axis='x', which='major', labelsize=25)
    ax.tick_params(axis='y', which='major', labelsize=25)
    ax.tick_params('both', length=6, width=1, which='major')
    ax.tick_params('both', length=4, width=1, which='minor')
    ax.set_title('QEC process fidelities')
    ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = '0.5',lw = 0.5)
    ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = '0.5',lw = 0.5)
    ax.set_xlabel('Error probability',fontsize=25)
    ax.set_ylabel('Process fidelity',fontsize=25)

    mpl.rcParams['axes.linewidth'] = 1

    lgd = ax.legend(loc = 1,frameon = False)
    for label in lgd.get_texts():
        label.set_fontsize(25)

    try:
        fig.savefig(
            os.path.join(folder,'Compare_syndromes.png'))
        fig.savefig(
            os.path.join(folder,'Compare_syndromes.pdf'))
    except:
        print 'Figure has not been saved.'

############################################
############################################
''' FIGURES THAT WE WANT TO PLOT FOR REAL'''
############################################
############################################


folder = r'D:\measuring\data\QEC_data\figs\final figures'


c_green = (9/255.,232/255.,94/255.)
c_grey = (64/255.,78/255.,77/255.)#(240/255.,242/255.,166/255.)
c_blue = (68/255.,204/255.,255/255.)
c_red = (150/255.,52/255.,132/255.)
c_orange = (242/255.,129/255.,35/255.)
c_orange_2 = (242/255.,129/255.,35/255.)

# c_green = (46/255.,198/255.,98/255.)
# c_grey = (64/255.,78/255.,77/255.)
# c_blue = (144/255.,170/255.,208/255.)
# c_red = (112/255.,22/255.,60/255.)
# c_orange = (242/255.,129/255.,35/255.)

def QEC_plot_process_fids_final():
    syndrome_list = ['00','01','10','11']

    process_dict = QEC_process_fids_sum_all(syndrome_list = syndrome_list)
    toff_process_dict = no_QEC_toffoli_fids()
    x = process_dict['x']

    t_list = ['1','2','3','avg']
    color_list = ['c','r','b','g']

    no_process_dict = no_QEC_process_fids(run = 0)
    single_process_dict = single_Qubit_no_QEC_process_fids()
    process_dict_idle = no_QEC_process_fids(idle = True)

    toff_dict_idle = no_QEC_toffoli_fids(idle = True)


    fig,ax = plt.subplots(figsize = (10,10))
    mpl.rcParams['pdf.fonttype'] = 42
    y = single_process_dict['dec_'+'avg'+'_y']
    y_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
    ax.plot(x_fit, y_fit, color = c_green,ls = '-', lw=1, label =  'Single qubit')#, $p_c$='+str(round(p_c*100)/100.)+'('+str(int(round(p_c_err*100)))+')')
    (_,caps,_)=ax.errorbar(x,y,yerr=y_err,color = c_green,markeredgecolor = c_green, ls = '',marker = 'o', ms = 7,capsize = 6)
    for cap in caps:
        cap.set_markeredgewidth(1)

    y = no_process_dict['dec_'+'avg'+'_y']
    y_err = no_process_dict['dec_'+'avg'+'_y_err']
    x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
    ax.plot(x_fit, y_fit, color = c_blue,ls = '-', lw=1, label =  'Encoded state')#, $p_c$='+str(round(p_c*100)/100.)+'('+str(int(round(p_c_err*100)))+')')
    (_,caps,_)=ax.errorbar(x,y,yerr=y_err,color = c_blue,markeredgecolor = c_blue, ls = '',marker = 'o', ms = 7,capsize = 6)
    for cap in caps:
        cap.set_markeredgewidth(1)

    y = process_dict['dec_'+'avg'+'_y']
    y_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
    ax.plot(x_fit, y_fit, color = c_red, lw=1, label =  'QEC, symmetrized read-out')#, $p_c$='+str(round(p_c*100)/100.)+'('+str(int(round(p_c_err*100)))+')')
    (_,caps,_)=ax.errorbar(x,y,yerr=y_err,color = c_red,markeredgecolor = c_red,ls = '',marker = 'o', ms = 7,capsize = 6)
    for cap in caps:
        cap.set_markeredgewidth(1)

    y = process_dict['dec_'+'avg'+'_y_new']
    y_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
    ax.plot(x_fit, y_fit, color = c_red,ls = '-.', lw=1, label =  'No QEC')#, $p_c$='+str(round(p_c*100)/100.)+'('+str(int(round(p_c_err*100)))+')')
    (_,caps,_)=ax.errorbar(x,y,yerr=y_err,color = c_red,markeredgecolor = c_red, ls = '',marker = 'o', ms = 7,capsize = 6)
    for cap in caps:
        cap.set_markeredgewidth(1)








    ax.set_ylim(-0,1)
    ax.set_xlim(-0.01,1.01)

    # ax.set_xticks([0,0.25,1])
    # ax.set_yticks([0,0.25,1])

    ax.set_xticks(np.arange(0,1.1,0.5))
    ax.set_xticks(np.arange(0,1.1,0.1), minor = True)
    ax.set_yticks(np.arange(0,1.1,0.5))
    ax.set_yticks(np.arange(0,1.1,0.1), minor = True)
    ax.tick_params(axis='x', which='major', labelsize=25)
    ax.tick_params(axis='y', which='major', labelsize=25)
    ax.tick_params('both', length=6, width=1, which='major')
    ax.tick_params('both', length=4, width=1, which='minor')
    # ax.set_title('QEC process fidelities')
    # ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = '0.5',lw = 0.5)
    # ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = '0.5',lw = 0.5)
    ax.set_xlabel('Error probability',fontsize=25)
    ax.set_ylabel('Process fidelity',fontsize=25)

    mpl.rcParams['axes.linewidth'] = 1.

    lgd = ax.legend(loc = 3,frameon=False)#loc = 2, bbox_to_anchor = (1,1))
    for label in lgd.get_texts():
        label.set_fontsize(25)

    #############################
    ###### INSET ################
    #############################

    a = axes([.525, .525, .35, .35])
    syndrome_list = ['00','01','10','11']
    p_dict = {}
    for syndrome in syndrome_list:

        if syndrome == '00' or syndrome == '01':
            run_list = [1,2,3]
        elif syndrome == '10':
            run_list = [2]
        elif syndrome == '11':
            run_list = [3]

        p_list = ['p00','p01','p10','p11']
        p_dict[syndrome] = QEC_sum_probs(run_list = run_list ,no_error = syndrome)

    p_dict['no error'] = 1/4.*(p_dict['00']['p00']+p_dict['01']['p01']+p_dict['10']['p10']+p_dict['11']['p11'])
    p_dict['Q3'] = 1/4.*(p_dict['00']['p01']+p_dict['01']['p00']+p_dict['10']['p11']+p_dict['11']['p10'])  # Carbon 5
    p_dict['Q1'] = 1/4.*(p_dict['00']['p10']+p_dict['01']['p11']+p_dict['10']['p00']+p_dict['11']['p01'])  # Carbon 2
    p_dict['Q2'] = 1/4.*(p_dict['00']['p11']+p_dict['01']['p10']+p_dict['10']['p01']+p_dict['11']['p00'])  # Carbon 1

    color = [c_green,c_orange,c_red,'r']
    pin = [0.127,0.113,0.123,0.145]

    # pin_no =0.0776# 0.093 # 0 #
    # pin_c1 =0.0917#0.072 # 0 #
    # pin_c2 =0.0645# 0.093 # 0 #
    # pin_c5 =0.0770# 0.115 # 0 #

    F = 0.939

    # p = linspace(0,1,1000)

    # ptot_no = p + pin_no-2*p*pin_no
    # ptot_c1 = p + pin_c1-2*p*pin_c1
    # ptot_c2 = p + pin_c2-2*p*pin_c2
    # ptot_c5 = p + pin_c5-2*p*pin_c5

    # p_no_error = 1-3*ptot_no+3*ptot_no**2
    # p_c1       = ptot_c1-ptot_c1**2
    # p_c2       = ptot_c2-ptot_c2**2
    # p_c5       = ptot_c5-ptot_c5**2

    # P_D_no_error = p_no_error*F**2 + p_c1*(1-F)**2 + (p_c2+p_c5)*F*(1-F)
    # P_D_c1 = p_no_error*(1-F)**2 + p_c1*F**2 + (p_c2+p_c5)*F*(1-F)
    # P_D_c2 = p_no_error*F*(1-F) + p_c1*(1-F)*F + (p_c2)*F*F+ (p_c5)*(1-F)*(1-F)
    # P_D_c5 = p_no_error*F*(1-F) + p_c1*(1-F)*F + (p_c5)*F*F+ (p_c2)*(1-F)*(1-F)


    # obtained from corrected way 150428
    pin_no = 0.0775644
    pin_c1 = 0.0640236
    pin_c2 = 0.0912455
    pin_c5 = 0.077205

    F1 = 0.988
    F0 = 0.890

    p = linspace(0,1,1000)

    ptot_no = p + pin_no-2*p*pin_no
    # ptot_no_alt = p + pin_no_alt-2*p*pin_no_alt
    ptot_c1 = p + pin_c1-2*p*pin_c1
    ptot_c2 = p + pin_c2-2*p*pin_c2
    ptot_c5 = p + pin_c5-2*p*pin_c5

    # p_no_error = 1-3*ptot_no+3*ptot_no**2

    # p_c1       = ptot_c1-ptot_c1**2
    # p_c2       = ptot_c2-ptot_c2**2
    # p_c5       = ptot_c5-ptot_c5**2

    #correct way to do this 150428
    p_no_error = 1-ptot_c1-ptot_c2-ptot_c5+ptot_c1*ptot_c2+ptot_c1*ptot_c5+ptot_c2*ptot_c5
    p_c1 = ptot_c1-ptot_c1*ptot_c2-ptot_c1*ptot_c5+ptot_c2*ptot_c5
    p_c2 = ptot_c2-ptot_c2*ptot_c1-ptot_c2*ptot_c5+ptot_c1*ptot_c5
    p_c5 = ptot_c5-ptot_c5*ptot_c1-ptot_c5*ptot_c2+ptot_c1*ptot_c2

    P_D_no_error = p_no_error*F**2 + p_c1*(1-F)**2 + (p_c2+p_c5)*F*(1-F)
    P_D_c1 = p_no_error*(1-F)**2 + p_c1*F**2 + (p_c2+p_c5)*F*(1-F)
    P_D_c2 = p_no_error*F*(1-F) + p_c1*(1-F)*F + (p_c2)*F*F+ (p_c5)*(1-F)*(1-F)
    P_D_c5 = p_no_error*F*(1-F) + p_c1*(1-F)*F + (p_c5)*F*F+ (p_c2)*(1-F)*(1-F)

    color = [c_green,c_orange,c_red,'r']
    label = ['No error','Qubit 1','Qubit 2', 'Qubit 3']

    plt.plot(p,P_D_no_error, color = color[0])
    # plt.plot(p,P_D_no_error_alt, color = color[0],ls = ':')
    plt.plot(p,P_D_c1, color = color[2])
    plt.plot(p,P_D_c2, color = color[1])
    plt.plot(p,P_D_c5, color = color[3])

    for i,error in enumerate(['no error','Q1','Q2','Q3']):
    #     if error == 'no error':
    #         x_fit,y_fit= fit_no_error_detection_curve(linspace(0,1,100),0.25*np.ones(100), return_errorbar=False,plot_guess = True)
    #     else:
    #         x_fit,y_fit= fit_error_detection_curve(linspace(0,1,100),0.25*np.ones(100),pin = pin[i], return_errorbar=False,plot_guess = True)
        
        plot(p_dict[syndrome]['x'],p_dict[error],label = error, color = color[i], marker = 'o',markersize = 5,ls = '',markeredgecolor = color[i])
        
    lgd = a.legend(loc = 9,frameon=False)
    for label in lgd.get_texts():
        label.set_fontsize(25)
    
    plt.xlim([-0.01,1.01])
    plt.xlabel('$p_e$',fontsize = 25)
    plt.ylim([0,1])
    plt.xticks([0,0.5,1])
    plt.yticks([0,0.5,1])
    plt.ylabel('Normalized \n occurence',fontsize = 25)

    plt.xticks(np.arange(0,1.1,0.5))
    # plt.xticks(np.arange(0,1.1,0.1), minor = True)
    plt.yticks(np.arange(0,1.1,0.5))
    # plt.yticks(np.arange(0,1.1,0.1), minor = True)
    plt.tick_params(axis='x', which='major', labelsize=25)
    plt.tick_params(axis='y', which='major', labelsize=25)
    plt.tick_params('both', length=6, width=1, which='major')
    plt.tick_params('both', length=4, width=1, which='minor')

    print p_dict['no error'][0]
    print p_dict['Q1'][0]
    print p_dict['Q2'][0]
    print p_dict['Q3'][0]
    print 1/3.*(p_dict['Q1'][0]+p_dict['Q2'][0]+p_dict['Q3'][0])


    try:
        fig.savefig(
            os.path.join(folder,'Process_fidelity_full_curve.png'),bbox_extra_artists = (lgd,),bbox_inches='tight')
        fig.savefig(
            os.path.join(folder,'Process_fidelity_full_curve.pdf'),bbox_extra_artists = (lgd,),bbox_inches='tight')
    except:
        print 'Figure has not been saved.'            

def QEC_plot_process_fids_11_vs_idle_full():
    fig,ax = plt.subplots(figsize= (10,10))
    mpl.rcParams['pdf.fonttype'] = 42
    process_dict = QEC_process_fids_sum_runs(run_list = [3],no_error = '11')

    x = process_dict['x']
    process_dict_idle = no_QEC_process_fids(idle = True, run = 0)

    single_process_dict = single_Qubit_no_QEC_process_fids()

    y_idle = process_dict_idle['dec_'+'avg'+'_y']
    y_idle_err = process_dict_idle['dec_'+'avg'+'_y_err']
    x_fit_idle, y_fit_idle, p_c, p_c_err= fit_QEC_process_curve(x,y_idle,return_errorbar = True)    
    ax.plot(x_fit_idle, y_fit_idle, color = c_grey,ls = '-', lw=2,label =  'Encoded state, idling')#, $p_c$='+str(round(p_c*100)/100.)+'('+str(int(round(p_c_err*100)))+')')
    (_,caps,_)=ax.errorbar(x,y_idle,yerr=y_idle_err,color = c_grey,markeredgecolor = c_grey, ls = '',marker = 'o', ms = 7,capsize = 6)
    for cap in caps:
        cap.set_markeredgewidth(1)

    y = process_dict['dec_'+'avg'+'_y']
    y_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve_11(x,y,return_errorbar = True)
    ax.plot(x_fit, y_fit, color = c_red, lw=2,label =  'QEC, optimized read-out')#, $p_c$='+str(round(p_c*100)/100.)+'('+str(int(round(p_c_err*100)))+')')
    (_,caps,_)=ax.errorbar(x,y,yerr=y_err,color = c_red,markeredgecolor = c_red,ls = '',marker = 'o', ms = 7,capsize = 6)
    for cap in caps:
        cap.set_markeredgewidth(1)
    y = process_dict['dec_'+'avg'+'_y_new']
    y_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve_11(x,y,return_errorbar = True)
    ax.plot(x_fit, y_fit, color = c_red,ls = '-.', lw=2,label =  'No QEC')#, $p_c$='+str(round(p_c*100)/100.)+'('+str(int(round(p_c_err*100)))+')')
    (_,caps,_)=ax.errorbar(x,y,yerr=y_err,color = c_red,markeredgecolor = c_red, ls = '',marker = 'o', ms = 7,capsize = 6)
    for cap in caps:
        cap.set_markeredgewidth(1)
    


    # y = single_process_dict['dec_'+'avg'+'_y']
    # y_err = process_dict['dec_'+'avg'+'_y_err']
    # x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
    # ax.plot(x_fit, y_fit, color = c_green,ls = '-', lw=2,label =  'Single qubit, $p_c$='+str(round(p_c*100)/100.)+'('+str(int(round(p_c_err*100)))+')')
    # (_,caps,_)=ax.errorbar(x,y,yerr=y_err,color = c_green,markeredgecolor = c_green, ls = '',marker = 'o', ms = 7,capsize = 6)
    # for cap in caps:
    #     cap.set_markeredgewidth(1)

    ax.set_ylim(-0,1)
    ax.set_xlim(-0.01,1.01)

    ax.set_xticks(np.arange(0,1.1,0.5))
    ax.set_xticks(np.arange(0,1.1,0.1), minor = True)
    ax.set_yticks(np.arange(0,1.1,0.5))
    ax.set_yticks(np.arange(0,1.1,0.1), minor = True)
    ax.tick_params(axis='x', which='major', labelsize=25)
    ax.tick_params(axis='y', which='major', labelsize=25)
    ax.tick_params('both', length=6, width=1, which='major')
    ax.tick_params('both', length=4, width=1, which='minor')
    # ax.set_title('QEC process fidelities')
    # ax.hlines([0.25,0.5],x[0]-1,x[-1]+1,linestyles='dotted', color = '0.5',lw 2 .5)
    # ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = '0.5',lw 2 .5)
    ax.set_xlabel('Error probability',fontsize=25)
    ax.set_ylabel('Process fidelity',fontsize=25)

    mpl.rcParams['axes.linewidth'] = 1

    lgd = ax.legend(loc = 3,frameon = False)
    for label in lgd.get_texts():
        label.set_fontsize(25)

    rectangle = plt.Rectangle((0.08, 0.5), 0.19, 0.15,edgecolor = '0.6', fill = None, lw = 2 )
    plt.gca().add_patch(rectangle)

    #############################
    ###### INSET ################
    #############################

    a = axes([.515, .520, .35, .35])
    process_dict = QEC_process_fids_sum_runs(run_list = [5,6,7],no_error = '11')
    
    x = process_dict['x']
    process_dict_idle = no_QEC_process_fids_sum_runs(idle = True, run_list = [2,3,4])

    y = process_dict['dec_'+'avg'+'_y']
    y_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve_11(x,y,return_errorbar = True)
    plot(x_fit, y_fit, color = c_red, lw=1.5, label =  'QEC, $p_c$='+str(int(p_c*1000)/1000.)+'('+str(int(p_c_err*1000))+')')
    (_,caps,_)=errorbar(x,y,yerr=y_err,color = c_red,markeredgecolor = c_red, ls = '',marker = 'o', ms = 5,capsize = 5, elinewidth = 2)
    for cap in caps:
        cap.set_markeredgewidth(1)
    y_idle = process_dict_idle['dec_'+'avg'+'_y']
    y_idle_err = process_dict_idle['dec_'+'avg'+'_y_err']
    x_fit_idle, y_fit_idle, p_err = fit_QEC_process_curve(x,y_idle)    
    plot(x_fit_idle, y_fit_idle, color = c_grey,ls = '-', lw=1.5, label =  'Idle, $p_c$='+str(int(p_err*100)/100.)) 
    (_,caps,_)=errorbar(x,y_idle,yerr=y_idle_err,color = c_grey,markeredgecolor = c_grey, ls = '',marker = 'o', ms = 5,capsize = 5, elinewidth = 2)
    for cap in caps:
        cap.set_markeredgewidth(1)
    # lgd = a.legend(loc = 9,frameon=False)
    # for label in lgd.get_texts():
    #     label.set_fontsize(18)
    
    plt.xlim([0.08,0.27])
    # plt.xlabel('$p_e$',fontsize = 25)
    plt.ylim([0.5,0.65])
    plt.xticks([0.1,0.2,0.3])
    plt.yticks([0.5,0.6])
    # plt.ylabel('Process fidelity',fontsize = 25)

    # plt.xticks(np.arange(0,1.1,0.5))
    # plt.xticks(np.arange(0,1.1,0.1), minor = True)
    # plt.yticks(np.arange(0,1.1,0.5))
    # plt.yticks(np.arange(0,1.1,0.1), minor = True)
    plt.tick_params(axis='x', which='major', labelsize=25)
    plt.tick_params(axis='y', which='major', labelsize=25)

    try:
        fig.savefig(
            os.path.join(folder,'11_vs_idle_full_curve_2.png'))
        fig.savefig(
            os.path.join(folder,'11_vs_idle_full_curve_2.pdf'))
    except:
        print 'Figure has not been saved.'

def QEC_plot_process_fids_11_vs_idle_zoom():
    fig,ax = plt.subplots()

    process_dict = QEC_process_fids_sum_runs(run_list = [5,6,7],no_error = '11')
    
    x = process_dict['x']
    process_dict_idle = no_QEC_process_fids_sum_runs(idle = True, run_list = [2,3,4])

    y = process_dict['dec_'+'avg'+'_y']
    y_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit, y_fit, p_c, p_c_err = fit_QEC_process_curve(x,y,return_errorbar = True)
    ax.plot(x_fit, y_fit, 'r', lw=1, label =  'QEC, $p_c$='+str(int(p_c*1000)/1000.)+'('+str(int(p_c_err*1000))+')')
    ax.errorbar(x,y,yerr=y_err,color = 'r',ls = '',marker = 'o',ms = 2)#,label =  'QEC')

    y_undo = process_dict['dec_'+'avg'+'_y_new']
    y_undo_err = process_dict['dec_'+'avg'+'_y_err']
    x_fit_undo, y_fit_undo, p_c, p_c_err = fit_QEC_process_curve(x,y_undo,return_errorbar = True)
    ax.plot(x_fit_undo, y_fit_undo, 'r',ls = ':', lw=1, label =  'Undo QEC, $p_c$='+str(int(p_c*1000)/1000.)+'('+str(int(p_c_err*1000))+')')
    ax.errorbar(x,y_undo,yerr=y_undo_err,color = 'r', ls = '',marker = 'o',ms = 2)#,label =  'undo QEC')

    y_idle = process_dict_idle['dec_'+'avg'+'_y']
    y_idle_err = process_dict_idle['dec_'+'avg'+'_y_err']
    x_fit_idle, y_fit_idle, p_err = fit_QEC_process_curve(x,y_idle)    
    ax.plot(x_fit_idle, y_fit_idle, 'k',ls = '-', lw=1, label =  'Idle, $p_c$='+str(int(p_err*100)/100.)) 
    ax.errorbar(x,y_idle,yerr=y_idle_err,color = 'k', ls = '',marker = 'o',ms = 2)#,label =  'Idle')

    ax.set_ylim(0.3,0.7)
    ax.set_xlim(-0.01,0.51)
    ax.set_title('11_vs_Idle_zoom.png')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Process Fidelity')
    ax.legend()

    try:
        fig.savefig(
            os.path.join(folder,'11_vs_idle_zoom.png'))
    except:
        print 'Figure has not been saved.'

    fig,ax = plt.subplots()
    ax.errorbar(x,y-y_idle,yerr=(y_err**2+y_idle_err**2)**0.5,color = 'b',ls = '',marker = 'o',ms = 8)
    ax.plot(x_fit,y_fit-y_fit_idle,color = 'b',ls = '-')
    ax.set_ylim(0.0,0.04)
    ax.set_xlim(-0.01,0.51)
    ax.set_title('11 vs idle difference.png')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Deviation')
    try:
        fig.savefig(
            os.path.join(folder,'11_vs_idle_QEC_difference.png'))
        fig.savefig(
            os.path.join(folder,'11_vs_idle_QEC_difference.pdf'))
    except:
        print 'Figure has not been saved.'

def QEC_plot_sweep_time():
    no_error_list = ['11']
    parity_time = 2*(4.996e-6*34 +11.312e-6*48) +2*(13.616e-6*34+4.996e-6*34) + 2* 150e-6

    color = ['r','g','b']
    dataset_dict_full = {}
    no_QEC_data_dict = {}
    QEC_single_data_dict = {}

    fig1, ax1 = plt.subplots(figsize=(10,10))

    for RO in [0,1,2,6]:
        print RO
        dataset_dict_full[RO] = {}
        no_QEC_data_dict[RO] = {}
        QEC_single_data_dict[RO] = {}
        for state in ['Z','mZ']:
                print RO
                print state
                dataset_dict_full[RO][state] = QEC_sweep_time_sum_error_syns(state = state,RO = RO,run_list = no_error_list)
                no_QEC_data_dict[RO][state] =  no_QEC_data_single_state_RO_single_error_sign(sweep_time = True,idle = False,state = state,RO = RO, load_set = True,error_sign = 0)
                if RO != 6:
                    QEC_single_data_dict[RO][state] =  single_qubit_no_QEC_data_single_state_RO_single_error_sign(state = state,sweep_time = True, error_sign = -1, Qubit = RO+1, load_set = True)
                    data_list = [dataset_dict_full,no_QEC_data_dict,QEC_single_data_dict]
        # average Z and mZ data
        dataset_dict_full[RO]['x'] = dataset_dict_full[RO]['Z']['x']
        dataset_dict_full[RO]['y'] = 1/2.*(dataset_dict_full[RO]['Z']['y']-dataset_dict_full[RO]['mZ']['y'])
        dataset_dict_full[RO]['y_no_corr'] = 1/2.*(dataset_dict_full[RO]['Z']['y_no_corr']-dataset_dict_full[RO]['mZ']['y_no_corr'])
        dataset_dict_full[RO]['y_err']= 1/2.*(dataset_dict_full[RO]['Z']['y_err']**2+dataset_dict_full[RO]['mZ']['y_err']**2)**0.5

        no_QEC_data_dict[RO]['x'] = no_QEC_data_dict[RO]['Z']['x']
        no_QEC_data_dict[RO]['y'] = 1/2.*(no_QEC_data_dict[RO]['Z']['y']-no_QEC_data_dict[RO]['mZ']['y'])
        no_QEC_data_dict[RO]['y_err']= 1/2.*(no_QEC_data_dict[RO]['Z']['y_err']**2+no_QEC_data_dict[RO]['mZ']['y_err']**2)**0.5

        if RO !=6:
            QEC_single_data_dict[RO]['x'] = QEC_single_data_dict[RO]['Z']['x']
            QEC_single_data_dict[RO]['y'] = 1/2.*(QEC_single_data_dict[RO]['Z']['y']-QEC_single_data_dict[RO]['mZ']['y'])
            QEC_single_data_dict[RO]['y_err']= 1/2.*(QEC_single_data_dict[RO]['Z']['y_err']**2+QEC_single_data_dict[RO]['mZ']['y_err']**2)**0.5

    # add best single qubit
    x_single = QEC_single_data_dict[1]['x']
    y_single = QEC_single_data_dict[1]['y']
    y_single_err = QEC_single_data_dict[1]['y_err']
    x_single = x_single*1000.
    x_temp, y_temp,T, T_err = fit_timesweep_single(x_single[0:-4],y_single[0:-4],return_errorbar = True)
    
    (_,caps,_) = ax1.errorbar(x_single[0:-4],1/2.*(y_single[0:-4]+1),yerr=1/2.*y_single_err[0:-4],
                color = c_green,markeredgecolor = c_green, ls = '',lw = 1,marker = 'o', ms = 7,capsize = 6, label = 'Un-encoded qubit')
    for cap in caps:
        cap.set_markeredgewidth(1)
    ax1.plot(x_temp,1/2.*(1+y_temp),color = c_green,ls = '-',lw = 1)

    y_toff_encode = 1/2.*(no_QEC_data_dict[0]['y']+no_QEC_data_dict[1]['y']+no_QEC_data_dict[2]['y']-no_QEC_data_dict[6]['y'])
    y_toff_encode_err = 1/2.*(no_QEC_data_dict[0]['y_err']**2+no_QEC_data_dict[1]['y_err']**2+no_QEC_data_dict[2]['y_err']**2+no_QEC_data_dict[6]['y_err']**2)**0.5
    x_enc = no_QEC_data_dict[0]['x']
    x_enc = x_enc*1000.
    x_temp, y_temp,T, T_err = fit_timesweep_single(x_enc[0:-1],y_toff_encode[0:-1],return_errorbar = True)
    
    (_,caps,_) = ax1.errorbar(x_enc[0:-1],1/2.*(y_toff_encode[0:-1]+1),yerr=1/2.*y_toff_encode_err[0:-1],
                color = c_blue,markeredgecolor = c_blue, ls = '',lw = 1,marker = 'o', ms = 7,capsize = 6, label = '1 round')
    for cap in caps:
        cap.set_markeredgewidth(1)
    ax1.plot(x_temp,1/2.*(1+y_temp),color = c_blue,ls = '-',lw = 1)

    x = dataset_dict_full[RO]['x']+ np.ones(len(dataset_dict_full[RO]['x']))*parity_time

    x = x*1000.

    fit_data_QEC = loadtxt('QEC.txt')
    fit_data_parity = loadtxt('parity.txt')

    y_toff_QEC = 1/2.*(dataset_dict_full[0]['y']+dataset_dict_full[1]['y']+dataset_dict_full[2]['y']-dataset_dict_full[6]['y'])
    y_toff_QEC_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
    (_,caps,_) = ax1.errorbar(x[0:-3],1/2.*(1+y_toff_QEC[0:-3]),yerr=1/2.*y_toff_QEC_err[0:-3],
                color = c_red,markeredgecolor = c_red, ls = '',lw = 1,marker = 'o', ms = 7,capsize = 6, label = '2 rounds')
    for cap in caps:
        cap.set_markeredgewidth(1)
    ax1.plot(fit_data_QEC[:,0][6:55],(fit_data_QEC[:,1][6:55]+1)/2.,color = c_red, ls = '-', lw = 1)

    y_toff_parity = 1/2.*(dataset_dict_full[0]['y_no_corr']+dataset_dict_full[1]['y_no_corr']+dataset_dict_full[2]['y_no_corr']-dataset_dict_full[6]['y_no_corr'])
    y_toff_parity_err = 1/2.*(dataset_dict_full[0]['y_err']**2+dataset_dict_full[1]['y_err']**2+dataset_dict_full[2]['y_err']**2+dataset_dict_full[6]['y_err']**2)**0.5
    x = dataset_dict_full[6]['x']+ np.ones(len(dataset_dict_full[6]['x']))*parity_time
    x = x*1000.
    (_,caps,_) = ax1.errorbar(x[0:-3],1/2.*(1+y_toff_parity[0:-3]),yerr=1/2.*y_toff_parity_err[0:-3],
                color = c_red,markeredgecolor = c_red, ls = '',lw = 1,marker = '*', ms = 9,capsize = 6, label = 'No feedback')
    for cap in caps:
        cap.set_markeredgewidth(1)
    ax1.plot(fit_data_parity[:,0][6:55],(fit_data_parity[:,1][6:55]+1)/2.,color = c_red, ls = '-.', lw = 1)



    # print T
    # print T_err



    print T
    print T_err

    ax1.set_xticks(np.arange(0,36,10))
    ax1.set_yticks(np.arange(0.5,1.1,0.25))

    ax1.tick_params(axis='x', which='major', labelsize=25)
    ax1.tick_params(axis='y', which='major', labelsize=25)

    ax1.hlines([0.5],x[0]-10,x[-1]+10,linestyles='dotted',color = '0.5', lw = 0.5)
    ax1.vlines([x[1],x[7]],-0.1,1.5,color = '0.5',lw = 1,linestyles = 'dashed')
    plt.axvspan(x[1],x[7], facecolor='y', alpha=0.1)
    # plt.axvspan(-1,x[1], facecolor='k', alpha=0.05)
    # plt.axvspan(x[7],35, facecolor='k', alpha=0.05)
    ax1.set_ylim(0.48,1.0)
    ax1.set_xlim(-1,30)
    ax1.set_xlabel('Time (ms)',fontsize = 25)
    ax1.set_ylabel('Average state fidelity',fontsize = 25)

    ax1.set_yticks(np.arange(0.5,1.05,0.05), minor = True)
    ax1.set_xticks(np.arange(0,31,2), minor = True)

    ax1.tick_params('both', length=4, width=1, which='minor')
    lgd = ax1.legend(loc = (0.60,0.65),frameon = False)
    for label in lgd.get_texts():
        label.set_fontsize(25)

    fig1.tight_layout()

    print x[1]
    print x[7]
    try:
        fig1.savefig(
            os.path.join(folder,'QEC_sweep_time.png'))
        fig1.savefig(
            os.path.join(folder,'QEC_sweep_time.pdf'))
    except:
        print 'Figure has not been saved.'

    # data = {}
    # data['x_QEC'] = x
    # data['y_toff_QEC'] = y_toff_QEC
    # data['y_toff_QEC_err'] = y_toff_QEC_err
    # data['x_parity'] = x
    # data['y_toff_parity'] = y_toff_parity
    # data['y_toff_parity_err'] = y_toff_parity_err
    # data['x_enc'] = x_enc
    # data['y_toff_encode'] = y_toff_encode
    # data['y_toff_encode_err'] = y_toff_encode_err
    # data['x_single'] = x_single
    # data['y_single'] = y_single
    # data['y_single_err'] = y_single_err
    # pickle.dump(data, open( "timesweep_data.p", "wb" ) )

def plot_prob_timesweep_11():
    syndrome = '11'
    
    parity_time = 2*(4.996e-6*34 +11.312e-6*48) +2*(13.616e-6*34+4.996e-6*34) + 2* 150e-6
    p_list = ['p11','p01','p00','p10']
    p_dict = QEC_timesweep_sum_probs(no_error = syndrome)

    fig,ax = plt.subplots(figsize = (7,7))
    color =  [c_green,c_orange,c_red,'r']
    label_list = ['no error', 'Qubit 1','Qubit 2', 'Qubit 3']

    x = (p_dict['x'][0:-3]+ np.ones(len(p_dict['x'][0:-3]))*parity_time)*1e3



    pin_no      =0.0796758#0.0856612
    # pin_no_alt  = 0.093
    pin_c1      =0.0492344#0.0650188
    pin_c2      =0.0803633#0.0822674
    pin_c5      =0.110111#0.11562

    pin_no      =0.0796758#0.0856612
    # pin_no_alt  = 0.093
    pin_c1      =0.0642991#0.0650188
    pin_c2      =0.0827476#0.0822674
    pin_c5      =0.11351#0.11562


    F1 = 0.988
    F0 = 0.890

    t  = np.linspace(x[0],x[-1],1000)
    
    T21 = 9.6
    T22 = 12.0
    T25 = 18.2


    p1h = 1/2.*(1-np.exp(-(0.5*t/T21)**2))
    p2h = 1/2.*(1-np.exp(-(0.5*t/T22)**2))
    p5h = 1/2.*(1-np.exp(-(0.5*t/T25)**2))

    # pavg = (p1h+p2h+p5h)/3.

    # ptot_no = pavg + pin_no-2*pavg*pin_no
    # ptot_no_alt = pavg + pin_no_alt-2*pavg*pin_no_alt
    ptot_c1 =p1h+ pin_c1-2*p1h*pin_c1
    ptot_c2 =p2h+ pin_c2-2*p2h*pin_c2
    ptot_c5 =p5h+ pin_c5-2*p5h*pin_c5

    # p_no_error = 1-3*ptot_no+3*ptot_no**2
    # p_no_error_alt = 1-3*ptot_no_alt+3*ptot_no_alt**2
    # p_c1       = ptot_c1-ptot_c1**2
    # p_c2       = ptot_c2-ptot_c2**2
    # p_c5       = ptot_c5-ptot_c5**2

    #correct way to do this 150428
    p_no_error = 1-ptot_c1-ptot_c2-ptot_c5+ptot_c1*ptot_c2+ptot_c1*ptot_c5+ptot_c2*ptot_c5
    p_c1 = ptot_c1-ptot_c1*ptot_c2-ptot_c1*ptot_c5+ptot_c2*ptot_c5
    p_c2 = ptot_c2-ptot_c2*ptot_c1-ptot_c2*ptot_c5+ptot_c1*ptot_c5
    p_c5 = ptot_c5-ptot_c5*ptot_c1-ptot_c5*ptot_c2+ptot_c1*ptot_c2

    #### For 11 assignment
    if syndrome == '11':
        P_D_no_error = p_no_error*F1**2 + p_c1*(1-F0)**2 + (p_c2+p_c5)*F1*(1-F0)
        # P_D_no_error_alt = p_no_error_alt*F1**2 + p_c1*(1-F0)**2 + (p_c2+p_c5)*F1*(1-F0)
        P_D_c1 = p_no_error*(1-F1)**2 + p_c1*F0**2 + (p_c2+p_c5)*F0*(1-F1)
        P_D_c2 = p_no_error*F1*(1-F1) + p_c1*(1-F0)*F0 + (p_c2)*F1*F0+ (p_c5)*(1-F1)*(1-F0)
        P_D_c5 = p_no_error*F1*(1-F1) + p_c1*(1-F0)*F0 + (p_c5)*F1*F0+ (p_c2)*(1-F1)*(1-F0)



    for jj,p in enumerate(p_list):
        ax.plot(x,p_dict[p][0:-3],color = color[jj],label = label_list[jj],
            marker = 'o',ms = 5,markeredgecolor = color[jj],ls = '')
        print label_list[jj]
        print x[0]
        print p_dict[p][0]
    plt.plot(t,P_D_no_error, color = color[0])
    # plt.plot(t,P_D_no_error_alt, color = color[0],ls = ':')
    plt.plot(t,P_D_c1, color = color[2])
    plt.plot(t,P_D_c2, color = color[1])
    plt.plot(t,P_D_c5, color = color[3])

    ax.set_xticks(np.arange(0,30,10))
    ax.set_xticks(np.arange(0,30,5), minor = True)
    ax.set_yticks(np.arange(0,1.1,0.5))
    ax.set_yticks(np.arange(0,1.1,0.25),minor = True)
    ax.set_ylim(0,1.0)
    ax.set_xlim(0,28)
    ax.set_xlabel('Time (ms)',fontsize = 25)
    ax.set_ylabel('Occurence',fontsize = 25)
    ax.tick_params(axis='x', which='major', labelsize=25)
    ax.tick_params(axis='y', which='major', labelsize=25)
    # ax.set_yticks(np.arange(0.5,1.05,0.05), minor = True)
    # ax.set_xticks(np.arange(0,29,2), minor = True)

    ax.tick_params('both', length=4, width=1, which='minor')
    lgd = ax.legend(loc = (0.60,0.65),frameon = False)
    for label in lgd.get_texts():
        label.set_fontsize(25)

    fig.tight_layout()  

    ax.hlines([0.301401466364],-1,29,linestyles='dotted', color = c_green)
    ax.hlines([0.203400846453],-1,29,linestyles='dotted', color = c_red)
    ax.hlines([0.247598833613],-1,29,linestyles='dotted', color = c_orange)
    ax.hlines([0.247598857149],-1,29,linestyles='dotted', color = 'r')

    try:
        fig.savefig(
            os.path.join(folder,'timesweep_probability_plot_'+syndrome+'.pdf'))
    except:
        print 'Figure has not been saved.' 

def QEC_multiple_rounds():
    save_folder = folder
    ### load 0 round data

    single_no_QEC_data_dict_Z_Q1  =  single_qubit_avg_state(state = 'Z',run = 1)
    single_no_QEC_data_dict_mZ_Q1 =  single_qubit_avg_state(state = 'mZ',run = 1)

    y_Z0         = single_no_QEC_data_dict_Z_Q1['y']
    y_err_Z0     = single_no_QEC_data_dict_Z_Q1['y_err']
    y_mZ0        = single_no_QEC_data_dict_mZ_Q1['y']
    y_err_mZ0    = single_no_QEC_data_dict_mZ_Q1['y_err']
    x0           = single_no_QEC_data_dict_mZ_Q1['x']

    ### load 1 round data (Majority voting/Toffoli)
    toff_process_dict = no_QEC_toffoli_fids(run_list = [0,1,2,3],state_list=['Z','mZ'],add_4 = True,do_weighted = True)
    y_Z1              = toff_process_dict['toff_Zy']
    y_err_Z1          = toff_process_dict['toff_Zy_err']
    y_mZ1             = toff_process_dict['toff_mZy']
    y_err_mZ1         = toff_process_dict['toff_mZy_err']
    x1                = toff_process_dict['x']

    # load 1 round data (QEC)
    data_dict_Z       = QEC_state_sum_RO_ZmZ(state = 'Z')
    data_dict_mZ      = QEC_state_sum_RO_ZmZ(state = 'mZ')

    y_Z1b         = data_dict_Z['y']
    y_err_Z1b     = data_dict_Z['y_err']
    y_mZ1b        = data_dict_mZ['y']
    y_err_mZ1b    = data_dict_mZ['y_err']
    x1b           = data_dict_mZ['x']

    ### load 2 rounds data

    x2, y_Z2, y_err_Z2, y_mZ2, y_err_mZ2 =  QEC_2rounds_combined_runs(runs=[1,2,3])

    ### load 3 rounds data
    x3, y_Z3, y_err_Z3, y_mZ3, y_err_mZ3, error_probs =  QEC_3rounds_combined_runs(runs=[1,2])

    ### Averaging over Z and mZ
    y_tot0  = (y_Z0-y_mZ0)/2; y_err_tot0 = (y_err_Z0**2+y_err_mZ0**2)**0.5/2
    y_tot1  = (y_Z1-y_mZ1)/2; y_err_tot1 = (y_err_Z1**2+y_err_mZ1**2)**0.5/2
    y_tot1b = (y_Z1b-y_mZ1b)/2; y_err_tot1b = (y_err_Z1b**2+y_err_mZ1b**2)**0.5/2
    y_tot2  = (y_Z2-y_mZ2)/2; y_err_tot2 = (y_err_Z2**2+y_err_mZ2**2)**0.5/2
    y_tot3  = (y_Z3-y_mZ3)/2; y_err_tot3 = (y_err_Z3**2+y_err_mZ3**2)**0.5/2

    print type(y_tot1)
    print type(y_tot2)

    ### Fitting
    x_fit0, y_fit0, fit_result0,u_fit_result0 = fit_QEC_curve(x0,y_tot0, return_errorbar=True)
    x_fit1, y_fit1, fit_result1,u_fit_result1 = fit_QEC_curve(x1,y_tot1, return_errorbar=True)
    x_fit2, y_fit2, fit_result2,u_fit_result2 = fit_QEC_2_rounds_curve2(x2,y_tot2, return_errorbar=True)
    x_fit3, y_fit3, fit_result3,u_fit_result3 = fit_QEC_3_rounds_curve2(x3,y_tot3, return_errorbar=True)
    
    x_fit1, y_fit1, fit_result1,u_fit_result1 = fit_QEC_curve(x1[0:7],y_tot1[0:7], return_errorbar=True)
    x_fit2, y_fit2, fit_result2,u_fit_result2 = fit_QEC_2_rounds_curve2_11(x2,y_tot2, return_errorbar=True)
    x_fit3, y_fit3, fit_result3,u_fit_result3 = fit_QEC_3_rounds_curve2_11(x3,y_tot3, return_errorbar=True)
    

    print '1 round'
    print fit_result1,u_fit_result1
    print '2 round'
    print fit_result2,u_fit_result2
    print '3 round'
    print fit_result3,u_fit_result3

    # x_fit1, y_fit1= fit_QEC_curve(x1,y_tot1, return_guess=True)

    fig4,ax = plt.subplots(figsize=(10,10))


    (_,caps,_)=ax.errorbar(x0, (0.5*y_tot0+0.5), yerr= 0.5*y_err_tot0, color = c_green,markeredgecolor = c_green, ls = '',marker = 'o', ms = 7,capsize = 6) 
    ax.plot(x_fit0, (0.5*y_fit0+0.5), color = c_green, lw = 1,label = 'Single qubit')#), $p_c$='+str(round(fit_result0*100)/100.)+'('+str(int(round(u_fit_result0*100)))+')')
    for cap in caps:
        cap.set_markeredgewidth(1)

    (_,caps,_)=ax.errorbar(x1, (0.5*y_tot1+0.5), yerr= 0.5*y_err_tot1,color = c_blue,markeredgecolor = c_blue, ls = '',marker = 'o', ms = 7,capsize = 6)
    ax.plot(x_fit1, (0.5*y_fit1+0.5), color = c_blue, lw = 1,label = 'Round C')#), $p_c$='+str(round(fit_result1*100)/100.)+'('+str(int(round(u_fit_result1*100)))+')')
    for cap in caps:
        cap.set_markeredgewidth(1)

    (_,caps,_)=ax.errorbar(x2, (0.5*y_tot2+0.5), yerr=0.5*y_err_tot2,color = c_red,markeredgecolor = c_red, ls = '',marker = 'o', ms = 7,capsize = 6)
    ax.plot(x_fit2, (0.5*y_fit2+0.5), color = c_red, lw = 1,label = 'Round B + C')#), $p_c$='+str(round(fit_result2*100)/100.)+'('+str(int(round(u_fit_result2*100)))+')')
    for cap in caps:
        cap.set_markeredgewidth(1)

    (_,caps,_)=ax.errorbar(x3, (0.5*y_tot3+0.5), yerr=0.5*y_err_tot3,color = c_orange_2,markeredgecolor = c_orange_2, ls = '',marker = 'o', ms = 7,capsize = 6)
    ax.plot(x_fit3, (0.5*y_fit3+0.5), color = c_orange_2, lw = 1,label = 'Round A + B + C')#), $p_c$='+str(round(fit_result3*100)/100.)+'('+str(int(round(u_fit_result3*100)))+')')
    for cap in caps:
        cap.set_markeredgewidth(1)

    ax.plot([0,0.5,1],[1,0.5,0],color = '0.5',ls = 'dotted')
    
    ax.set_ylim(0.48,1.0)
    ax.set_xlim(-0.01,0.51)


    ax.set_xticks(np.arange(0,0.51,0.25))
    # ax.set_xticks(np.arange(0,0.51,0.1), minor = True)
    ax.set_yticks(np.arange(0.5,1.01,0.25))
    # ax.set_yticks(np.arange(0.5,1.01,0.1), minor = True)

    ax.tick_params(axis='x', which='major', labelsize=25)
    ax.tick_params(axis='y', which='major', labelsize=25)
    ax.tick_params('both', length=6, width=1, which='major')
    # ax.tick_params('both', length=4, width=1, which='minor')
    # ax.set_title('QEC process fidelities')
    ax.hlines([0.5,1],x0[0]-1,x0[-1]+1,linestyles='dotted', color = '0.5',lw = 0.5)
    # ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = '0.5',lw = 0.5)
    ax.set_xlabel('Error probability $p_e$',fontsize=25)
    ax.set_ylabel('Average state fidelity',fontsize=25)

    mpl.rcParams['axes.linewidth'] = 1.

    fig4.tight_layout()

    lgd = ax.legend(loc = [0.01,0.03],frameon=False)#loc = 2, bbox_to_anchor = (1,1))
    for label in lgd.get_texts():
        label.set_fontsize(20)

    mpl.rcParams['pdf.fonttype'] = 42


    x, p_R1_11, p_R2_11 = QEC_3rounds_outcome_probability()
    a = axes([.62, .62, .3, .3])

    pin_list      =[0.0924663, 0.0856549]
    color = [c_red,c_orange_2]

    F1 = 0.988
    F0 = 0.890

    p_plot = linspace(0,1,1000)
    p = 1/2.*(1-(1-2*p_plot)**(1/3.))

    for j,pin in enumerate(pin_list):
        ptot = p + pin-2*p*pin


        p_no_error = 1-3*ptot+3*ptot**2
        p_c1       = ptot-ptot**2
        p_c2       = ptot-ptot**2
        p_c5       = ptot-ptot**2

        P_D_no_error = p_no_error*F1**2 + p_c1*(1-F0)**2 + (p_c2+p_c5)*F1*(1-F0)
        # P_D_c1 = p_no_error*(1-F1)**2 + p_c1*F0**2 + (p_c2+p_c5)*F0*(1-F1)
        # P_D_c2 = p_no_error*F1*(1-F1) + p_c1*(1-F0)*F0 + (p_c2)*F1*F0+ (p_c5)*(1-F1)*(1-F0)
        # P_D_c5 = p_no_error*F1*(1-F1) + p_c1*(1-F0)*F0 + (p_c5)*F1*F0+ (p_c2)*(1-F1)*(1-F0)
        plot(p_plot,P_D_no_error,color = color[j])


    plot(x, p_R1_11,label = 'In round 1', color = c_red, marker = 'o',markersize = 5,ls = '',markeredgecolor = c_red)
    plot(x, p_R2_11,label = 'In round 2', color = c_orange_2 , marker = 'o',markersize = 5,ls = '',markeredgecolor = c_orange_2 )
    
    print p_R1_11[0]
    print p_R2_11[0]

    lgd = a.legend(loc = 9,frameon=False)
    for label in lgd.get_texts():
        label.set_fontsize(25)
    
    plt.xlim([-0.02,0.52])
    plt.xlabel('$p_e$',fontsize = 25)
    plt.ylim([0.25,1.0])
    plt.ylabel('Normalized \n occurence',fontsize = 25)

    plt.xticks(np.arange(0.0,0.52,0.25))
    # plt.xticks(np.arange(0,1.1,0.1), minor = True)
    plt.yticks(np.arange(0.25,1.1,0.25))
    # plt.yticks(np.arange(0,1.1,0.1), minor = True)
    plt.tick_params(axis='x', which='major', labelsize=25)
    plt.tick_params(axis='y', which='major', labelsize=25)
    plt.tick_params('both', length=6, width=1, which='major')
    plt.tick_params('both', length=4, width=1, which='minor')

    ax.set_yticks(np.arange(0.5,1.05,0.05), minor = True)
    ax.set_xticks(np.arange(0,0.55,0.05), minor = True)

    ax.tick_params('both', length=4, width=1, which='minor')

    if save_folder != None:
        try:
            fig4.savefig(
                os.path.join(save_folder,'Multiple_rounds_Combined11.pdf'))
            fig4.savefig(
                os.path.join(save_folder,'Multiple_rounds_Combined11.png'))
        except:
            print 'Figure has not been saved.'

    plt.show()
    # plt.close('all')


def QEC_multiple_rounds_alt_fit():
    save_folder = folder
    ### load 0 round data

    single_no_QEC_data_dict_Z_Q1  =  single_qubit_avg_state(state = 'Z',run = 1)
    single_no_QEC_data_dict_mZ_Q1 =  single_qubit_avg_state(state = 'mZ',run = 1)

    y_Z0         = single_no_QEC_data_dict_Z_Q1['y']
    y_err_Z0     = single_no_QEC_data_dict_Z_Q1['y_err']
    y_mZ0        = single_no_QEC_data_dict_mZ_Q1['y']
    y_err_mZ0    = single_no_QEC_data_dict_mZ_Q1['y_err']
    x0           = single_no_QEC_data_dict_mZ_Q1['x']

    ### load 1 round data (Majority voting/Toffoli)
    toff_process_dict = no_QEC_toffoli_fids(run_list = [0,1,2,3],state_list=['Z','mZ'],add_4 = True,do_weighted = True)
    y_Z1              = toff_process_dict['toff_Zy']
    y_err_Z1          = toff_process_dict['toff_Zy_err']
    y_mZ1             = toff_process_dict['toff_mZy']
    y_err_mZ1         = toff_process_dict['toff_mZy_err']
    x1                = toff_process_dict['x']

    # load 1 round data (QEC)
    data_dict_Z       = QEC_state_sum_RO_ZmZ(state = 'Z')
    data_dict_mZ      = QEC_state_sum_RO_ZmZ(state = 'mZ')

    y_Z1b         = data_dict_Z['y']
    y_err_Z1b     = data_dict_Z['y_err']
    y_mZ1b        = data_dict_mZ['y']
    y_err_mZ1b    = data_dict_mZ['y_err']
    x1b           = data_dict_mZ['x']

    ### load 2 rounds data

    x2, y_Z2, y_err_Z2, y_mZ2, y_err_mZ2 =  QEC_2rounds_combined_runs(runs=[1,2,3])

    ### load 3 rounds data
    x3, y_Z3, y_err_Z3, y_mZ3, y_err_mZ3, error_probs =  QEC_3rounds_combined_runs(runs=[1,2])

    ### Averaging over Z and mZ
    y_tot0  = (y_Z0-y_mZ0)/2; y_err_tot0 = (y_err_Z0**2+y_err_mZ0**2)**0.5/2
    y_tot1  = (y_Z1-y_mZ1)/2; y_err_tot1 = (y_err_Z1**2+y_err_mZ1**2)**0.5/2
    y_tot1b = (y_Z1b-y_mZ1b)/2; y_err_tot1b = (y_err_Z1b**2+y_err_mZ1b**2)**0.5/2
    y_tot2  = (y_Z2-y_mZ2)/2; y_err_tot2 = (y_err_Z2**2+y_err_mZ2**2)**0.5/2
    y_tot3  = (y_Z3-y_mZ3)/2; y_err_tot3 = (y_err_Z3**2+y_err_mZ3**2)**0.5/2

    print type(y_tot1)
    print type(y_tot2)

    ### Fitting
    # x_fit0, y_fit0, fit_result0,u_fit_result0 = fit_QEC_curve(x0,y_tot0, return_errorbar=True)
    # x_fit1, y_fit1, fit_result1,u_fit_result1 = fit_QEC_curve(x1[0:7],y_tot1[0:7], return_errorbar=True)
    # x_fit1, y_fit1, fit_result1,u_fit_result1 = fit_QEC_curve(x1,y_tot1, return_errorbar=True)
    # x_fit2, y_fit2, fit_result2,u_fit_result2 = fit_QEC_2_rounds_curve2(x2,y_tot2, return_errorbar=True)
    # x_fit3, y_fit3, fit_result3,u_fit_result3 = fit_QEC_3_rounds_curve2(x3,y_tot3, return_errorbar=True)
    
    x_fit1, y_fit1, fit_result1,u_fit_result1 = fit_QEC_curve_11(x1,y_tot1, return_errorbar=True)
    x_fit2, y_fit2, fit_result2,u_fit_result2 = fit_QEC_curve_11_2(x2,y_tot2, return_errorbar=True)
    x_fit3, y_fit3, fit_result3,u_fit_result3 = fit_QEC_curve_11_3(x3,y_tot3, return_errorbar=True)
    

    print '1 round'
    print fit_result1,u_fit_result1
    print '2 round'
    print fit_result2,u_fit_result2
    print '3 round'
    print fit_result3,u_fit_result3

    # x_fit1, y_fit1= fit_QEC_curve(x1,y_tot1, return_guess=True)

    fig4,ax = plt.subplots(figsize=(10,10))


    # (_,caps,_)=ax.errorbar(x0, (0.5*y_tot0+0.5), yerr= 0.5*y_err_tot0, color = c_green,markeredgecolor = c_green, ls = '',marker = 'o', ms = 7,capsize = 6) 
    # ax.plot(x_fit0, (0.5*y_fit0+0.5), color = c_green, lw = 1,label = 'Single qubit')#), $p_c$='+str(round(fit_result0*100)/100.)+'('+str(int(round(u_fit_result0*100)))+')')
    # for cap in caps:
    #     cap.set_markeredgewidth(1)

    # (_,caps,_)=ax.errorbar(x1, (0.5*y_tot1+0.5), yerr= 0.5*y_err_tot1,color = c_blue,markeredgecolor = c_blue, ls = '',marker = 'o', ms = 7,capsize = 6)
    # ax.plot(x_fit1, (0.5*y_fit1+0.5), color = c_blue, lw = 1,label = 'Round C')#), $p_c$='+str(round(fit_result1*100)/100.)+'('+str(int(round(u_fit_result1*100)))+')')
    # for cap in caps:
    #     cap.set_markeredgewidth(1)

    # (_,caps,_)=ax.errorbar(x1, y_tot1, yerr= y_err_tot1,color = c_blue,markeredgecolor = c_blue, ls = '',marker = 'o', ms = 7,capsize = 6)
    # ax.plot(x_fit1, y_fit1, color = c_blue, lw = 1,label = 'Round C')#), $p_c$='+str(round(fit_result1*100)/100.)+'('+str(int(round(u_fit_result1*100)))+')')
    # for cap in caps:
    #     cap.set_markeredgewidth(1)


    (_,caps,_)=ax.errorbar(x2, (0.5*y_tot2+0.5), yerr=0.5*y_err_tot2,color = c_red,markeredgecolor = c_red, ls = '',marker = 'o', ms = 7,capsize = 6)
    ax.plot(x_fit2, (0.5*y_fit2+0.5), color = c_red, lw = 1,label = 'Round B + C')#), $p_c$='+str(round(fit_result2*100)/100.)+'('+str(int(round(u_fit_result2*100)))+')')
    for cap in caps:
        cap.set_markeredgewidth(1)

    (_,caps,_)=ax.errorbar(x3, (0.5*y_tot3+0.5), yerr=0.5*y_err_tot3,color = c_orange_2,markeredgecolor = c_orange_2, ls = '',marker = 'o', ms = 7,capsize = 6)
    ax.plot(x_fit3, (0.5*y_fit3+0.5), color = c_orange_2, lw = 1,label = 'Round A + B + C')#), $p_c$='+str(round(fit_result3*100)/100.)+'('+str(int(round(u_fit_result3*100)))+')')
    for cap in caps:
        cap.set_markeredgewidth(1)

    ax.plot([0,0.5,1],[1,0.5,0],color = '0.5',ls = 'dotted')
    
    ax.set_ylim(0.48,1.0)
    ax.set_xlim(-0.01,0.51)


    ax.set_xticks(np.arange(0,0.51,0.25))
    # ax.set_xticks(np.arange(0,0.51,0.1), minor = True)
    ax.set_yticks(np.arange(0.5,1.01,0.25))
    # ax.set_yticks(np.arange(0.5,1.01,0.1), minor = True)

    ax.tick_params(axis='x', which='major', labelsize=25)
    ax.tick_params(axis='y', which='major', labelsize=25)
    ax.tick_params('both', length=6, width=1, which='major')
    # ax.tick_params('both', length=4, width=1, which='minor')
    # ax.set_title('QEC process fidelities')
    ax.hlines([0.5,1],x0[0]-1,x0[-1]+1,linestyles='dotted', color = '0.5',lw = 0.5)
    # ax.vlines([0.5],-0.1,1.1,linestyles='dotted', color = '0.5',lw = 0.5)
    ax.set_xlabel('Error probability $p_e$',fontsize=25)
    ax.set_ylabel('Average state fidelity',fontsize=25)

    mpl.rcParams['axes.linewidth'] = 1.

    fig4.tight_layout()

    lgd = ax.legend(loc = [0.01,0.03],frameon=False)#loc = 2, bbox_to_anchor = (1,1))
    for label in lgd.get_texts():
        label.set_fontsize(20)

    mpl.rcParams['pdf.fonttype'] = 42


    # x, p_R1_11, p_R2_11 = QEC_3rounds_outcome_probability()
    # a = axes([.62, .62, .3, .3])

    # pin_list      =[0.0924663, 0.0856549]
    # color = [c_red,c_orange_2]

    # F1 = 0.988
    # F0 = 0.890

    # p_plot = linspace(0,1,1000)
    # p = 1/2.*(1-(1-2*p_plot)**(1/3.))

    # for j,pin in enumerate(pin_list):
    #     ptot = p + pin-2*p*pin


    #     p_no_error = 1-3*ptot+3*ptot**2
    #     p_c1       = ptot-ptot**2
    #     p_c2       = ptot-ptot**2
    #     p_c5       = ptot-ptot**2

    #     P_D_no_error = p_no_error*F1**2 + p_c1*(1-F0)**2 + (p_c2+p_c5)*F1*(1-F0)
    #     # P_D_c1 = p_no_error*(1-F1)**2 + p_c1*F0**2 + (p_c2+p_c5)*F0*(1-F1)
    #     # P_D_c2 = p_no_error*F1*(1-F1) + p_c1*(1-F0)*F0 + (p_c2)*F1*F0+ (p_c5)*(1-F1)*(1-F0)
    #     # P_D_c5 = p_no_error*F1*(1-F1) + p_c1*(1-F0)*F0 + (p_c5)*F1*F0+ (p_c2)*(1-F1)*(1-F0)
    #     plot(p_plot,P_D_no_error,color = color[j])


    # plot(x, p_R1_11,label = 'In round 1', color = c_red, marker = 'o',markersize = 5,ls = '',markeredgecolor = c_red)
    # plot(x, p_R2_11,label = 'In round 2', color = c_orange_2 , marker = 'o',markersize = 5,ls = '',markeredgecolor = c_orange_2 )
    
    # print p_R1_11[0]
    # print p_R2_11[0]

    # lgd = a.legend(loc = 9,frameon=False)
    # for label in lgd.get_texts():
    #     label.set_fontsize(25)
    
    # plt.xlim([-0.02,0.52])
    # plt.xlabel('$p_e$',fontsize = 25)
    # plt.ylim([0.24,1.1])
    # plt.ylabel('Normalized \n occurence',fontsize = 25)

    # plt.xticks(np.arange(0.0,0.52,0.25))
    # # plt.xticks(np.arange(0,1.1,0.1), minor = True)
    # plt.yticks(np.arange(0.25,1.1,0.25))
    # # plt.yticks(np.arange(0,1.1,0.1), minor = True)
    # plt.tick_params(axis='x', which='major', labelsize=25)
    # plt.tick_params(axis='y', which='major', labelsize=25)
    # plt.tick_params('both', length=6, width=1, which='major')
    # plt.tick_params('both', length=4, width=1, which='minor')

    # ax.set_yticks(np.arange(0.5,1.05,0.05), minor = True)
    # ax.set_xticks(np.arange(0,0.55,0.05), minor = True)

    ax.tick_params('both', length=4, width=1, which='minor')

    if save_folder != None:
        try:
            fig4.savefig(
                os.path.join(save_folder,'Multiple_rounds_Combined11_fitted_single_roundvs2.pdf'))
            fig4.savefig(
                os.path.join(save_folder,'Multiple_rounds_Combined11_fitted_single_roundvs2.png'))
        except:
            print 'Figure has not been saved.'

    plt.show()
    # plt.close('all')

