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

def Contrast_Plot_QEC(timestamps=[None, None], measurement_name = ['adwindata'],folder_name ='QEC',
        post_select_QEC = False, ssro_calib_timestamp =None, do_plot = False, return_data = False) :

    '''
    Function that makes a plot with errorbars of MBI type data that has been measured with a positive
    and negative RO.
    '''       

    a_p, x, c0_p, u_c0_p, c0_00_p, u_c0_00_p, c0_01_p, u_c0_01_p, c0_10_p, u_c0_10_p, c0_11_p, u_c0_11_p, x_labels, folder_p = Plot_QEC(timestamp = timestamps[0], 
            measurement_name = measurement_name, folder_name = 'positive',
            ssro_calib_timestamp = ssro_calib_timestamp) 

    a_n, x, c0_n, u_c0_n, c0_00_n, u_c0_00_n, c0_01_n, u_c0_01_n, c0_10_n, u_c0_10_n, c0_11_n, u_c0_11_n, x_labels, folder_n = Plot_QEC(timestamp = timestamps[1], 
            measurement_name = measurement_name, folder_name = 'negative',
            ssro_calib_timestamp =ssro_calib_timestamp) 
        
    ### Combine data

        ## all data
    y = (c0_p - c0_n)/2.
    y_err =  1./2*(u_c0_p**2 + u_c0_n**2)**0.5      
    if do_plot == True:
        fig,ax = plt.subplots() 
        rects = ax.errorbar(x,y,yerr=y_err,color = 'k' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title(str(folder_p)+'/'+'\n' + script_name)
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')

        try:
            fig.savefig(
                os.path.join(folder_p,'QEC.png'))
        except:
            print 'Figure has not been saved.'

        ## ms=0 data
    y_00 = (c0_00_p - c0_00_n)/2.
    y_00_err =  1./2*(u_c0_00_p**2 + u_c0_00_n**2)**0.5 

    y_01 = (c0_01_p - c0_01_n)/2.
    y_01_err =  1./2*(u_c0_01_p**2 + u_c0_01_n**2)**0.5 

    y_10 = (c0_10_p - c0_10_n)/2.
    y_10_err =  1./2*(u_c0_10_p**2 + u_c0_10_n**2)**0.5 

    y_11 = (c0_11_p - c0_11_n)/2.
    y_11_err =  1./2*(u_c0_11_p**2 + u_c0_11_n**2)**0.5 
    if do_plot == True:
        fig,ax = plt.subplots() 
        ax.errorbar(x,y_00,yerr=y_00_err,color = 'c', label = 'y_00' )
        ax.errorbar(x,y_01,yerr=y_01_err,color = 'k', label = 'y_01' )
        ax.errorbar(x,y_10,yerr=y_10_err,color = 'm', label = 'y_10' )
        ax.errorbar(x,y_11,yerr=y_11_err,color = 'b', label = 'y_11' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        plt.legend()
        ax.set_title(str(folder_p)+'/'+'\n postselect')
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')


        try:
            fig.savefig(
                os.path.join(folder_p,'QEC_ps.png'))
        except:
            print 'Figure has not been saved.'


    # process the probabilities

    p00_avg = (a_p.p00+a_n.p00)/2
    p01_avg = (a_p.p01+a_n.p01)/2
    p10_avg = (a_p.p10+a_n.p10)/2
    p11_avg = (a_p.p11+a_n.p11)/2
    if do_plot == True:
        fig,ax = plt.subplots()
        ax.set_title(str(folder_p)+'/'+ '\n probabilities')
        ax.plot(x,p00_avg, 'c', label = 'p00')
        ax.plot(x,p01_avg, 'k', label = 'p01')
        ax.plot(x,p10_avg, 'm', label = 'p10')
        ax.plot(x,p11_avg, 'b', label = 'p11')
        plt.legend()
        ax.set_xlabel('error probability')
        ax.set_ylabel('outcome probability')

        try:
            fig.savefig(
                os.path.join(folder_p,'QEC_probabilities.png'))
        except:
            print 'Figure has not been saved.'

    if return_data == True:
        return x, y, y_err, y_00, y_00_err, p00_avg, y_01, y_01_err, p01_avg, y_10, y_10_err, p10_avg, y_11, y_11_err, p11_avg

def Contrast_Plot_Encoding(timestamps=[None, None], measurement_name = ['adwindata'],folder_name ='QEC',
        post_select_QEC = False, ssro_calib_timestamp =None, do_plot = False, return_data = False) :

    '''
    Function that makes a plot with errorbars of MBI type data that has been measured with a positive
    and negative RO.
    '''       

    a_p, x, c0_p, u_c0_p, x_labels, folder_p = Plot_Encoding(timestamp = timestamps[0], 
            measurement_name = measurement_name, folder_name = 'positive',
            ssro_calib_timestamp = ssro_calib_timestamp, return_raw = True) 

    a_n, x, c0_n, u_c0_n, x_labels, folder_n = Plot_Encoding(timestamp = timestamps[1], 
            measurement_name = measurement_name, folder_name = 'negative',
            ssro_calib_timestamp =ssro_calib_timestamp, return_raw = True) 
        
    ### Combine data

        ## all data
    y = (c0_p - c0_n)/2.
    y_err =  1./2*(u_c0_p**2 + u_c0_n**2)**0.5      
    if do_plot == True:
        fig,ax = plt.subplots() 
        ax.errorbar(x,y,yerr=y_err,color = 'k' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title(str(folder_p)+'/'+'\n' + script_name)
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')

        try:
            fig.savefig(
                os.path.join(folder_p,'QEC.png'))
        except:
            print 'Figure has not been saved.'
    
    if return_data == True:
        return x, y, y_err, y_00, y_00_err, p00_avg, y_01, y_01_err, p01_avg, y_10, y_10_err, p10_avg, y_11, y_11_err, p11_avg


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
    print folder_a
    print folder_b
    
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



''' These functions both load and plot individual measurements TODO THT, seperate loading from plotting'''

def Plot_QEC(timestamp = None, measurement_name = ['adwindata'],folder_name ='QEC',
        plot_post_select = False,
        ssro_calib_timestamp = None, save = True,
        do_plots = False, title =None ,fontsize = 10, post_select_QEC = True, return_raw = True, return_dict = False) :
    '''
    Function that makes a bar plot with errorbars of MBI type data
    '''
    plt.rc('font', size=fontsize)
    if timestamp == None:
        timestamp, folder   = toolbox.latest_data(folder_name,return_timestamp =True)
    else:
        folder = toolbox.data_from_time(timestamp)

    if ssro_calib_timestamp == None:
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
        # print ssro_calib_folder
    a = CP.ConditionalParityAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata',post_select_QEC = False)
    # print ssro_calib_folder
    a.get_electron_ROC(ssro_calib_folder)

    x_labels = a.sweep_pts.reshape(-1)

    ''' all data '''

    c0,u_c0 = a.convert_fidelity_to_contrast(a.p0,a.u_p0)
    x = x_labels

    if do_plots ==True:
        fig,ax = plt.subplots()
        ax.errorbar(x,c0,yerr=u_c0,color = 'k' )
        # ax.set_xticks(x)
        if title == None:
            ax.set_title(str(folder)+'/'+str(timestamp))
        else:
            ax.set_title(title)
        # print x_labels
        # ax.set_xticklabels(x_labels.tolist())
        ax.set_ylim(-1,1)
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

        try:
            fig.savefig(os.path.join(folder,'QEC.png'))
            fig.savefig(os.path.join(folder, title+'.pdf'),
                    format='pdf',bbox_inches='tight')
        except:
            print 'Figure A has not been saved.'
    
    if post_select_QEC == False and return_dict == True:
            data_dict = {}
            # data_dict['a'] = a
            data_dict['x'] =x 
            data_dict['c0'] =c0 
            data_dict['u_c0'] =u_c0 
            data_dict['x_labels'] =x_labels 
            # data_dict['folder'] =folder 

            return data_dict, folder 

    ''' postselected data '''
    if post_select_QEC == True:
        a = CP.ConditionalParityAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata',post_select_QEC = post_select_QEC)
        a.get_electron_ROC(ssro_calib_folder, post_select_QEC = post_select_QEC)
                
        c0_00,u_c0_00 =  a.convert_fidelity_to_contrast(a.p0_00,a.u_p0_00)
        c0_01,u_c0_01 =  a.convert_fidelity_to_contrast(a.p0_01,a.u_p0_01)
        c0_10,u_c0_10 =  a.convert_fidelity_to_contrast(a.p0_10,a.u_p0_10)
        c0_11,u_c0_11 =  a.convert_fidelity_to_contrast(a.p0_11,a.u_p0_11)
        
        if plot_post_select ==True:
            fig,ax = plt.subplots()
            ax.errorbar(x,c0_00,yerr=u_c0_00,color = 'k' )
            ax.errorbar(x,c0_01,yerr=u_c0_01,color = 'c' )
            ax.errorbar(x,c0_10,yerr=u_c0_10,color = 'g' )
            ax.errorbar(x,c0_11,yerr=u_c0_11,color = 'r' )
            # ax.set_xticks(x)
            if title == None:
                ax.set_title(str(folder)+'/'+str(timestamp))
            else:
                ax.set_title(title)
            # print x_labels
            # ax.set_xticklabels(x_labels.tolist())
            # ax.set_ylim(-1,1)
            ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

            try:
                fig.savefig(os.path.join(folder,'QEC.png'))
                fig.savefig(os.path.join(folder, title+'.pdf'),
                        format='pdf',bbox_inches='tight')
            except:
                print 'Figure A has not been saved.'
        if return_raw == True:
            return a, x, c0, u_c0, c0_00, u_c0_00, c0_01, u_c0_01, c0_10, u_c0_10, c0_11, u_c0_11, x_labels, folder
        
        elif return_dict:
            data_dict = {}
            # data_dict['a'] = a
            data_dict['x'] =x 
            data_dict['c0'] =c0 
            data_dict['u_c0'] =u_c0 
            data_dict['c0_00'] =c0_00 
            data_dict['u_c0_00'] =u_c0_00 
            data_dict['c0_01'] =c0_01 
            data_dict['u_c0_01'] =u_c0_01 
            data_dict['c0_10'] =c0_10 
            data_dict['u_c0_10'] =u_c0_10 
            data_dict['c0_11'] =c0_11 
            data_dict['u_c0_11'] =u_c0_11 
            data_dict['x_labels'] =x_labels 
            # data_dict['folder'] =folder 
            data_dict['p00'] = a.p00
            data_dict['p01'] = a.p01
            data_dict['p10'] = a.p10
            data_dict['p11'] = a.p11

            return data_dict, folder

def Plot_Encoding(timestamp = None, measurement_name = ['adwindata'], folder_name ='QEC',
        plot_post_select = False,
        ssro_calib_timestamp = None, save = True,
        do_plots = False, title =None ,fontsize = 10, return_raw = True, return_dict = False) :
    '''
    Load and plot encoding data (no postselection, because no parity measurements)
    '''

    ### Find data ###
    if timestamp == None:
        timestamp, folder   = toolbox.latest_data(folder_name,return_timestamp =True)
    else:
        folder = toolbox.data_from_time(timestamp)

    if ssro_calib_timestamp == None:
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
    
    ### Load data ###
    a = CP.ConditionalParityAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata',post_select_QEC = False)
    a.get_electron_ROC(ssro_calib_folder)

    x_labels = a.sweep_pts.reshape(-1)

    c0,u_c0 = a.convert_fidelity_to_contrast(a.p0,a.u_p0)
    x = x_labels


    ### Plotting (optional) ###
    if do_plots ==True:
        plt.rc('font', size=fontsize)
        fig,ax = plt.subplots()
        ax.errorbar(x,c0,yerr=u_c0,color = 'k' )
        # ax.set_xticks(x)
        if title == None:
            ax.set_title(str(folder)+'/'+str(timestamp))
        else:
            ax.set_title(title)
        # ax.set_xticklabels(x_labels.tolist())
        ax.set_ylim(-1,1)
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

        try:
            fig.savefig(os.path.join(folder,'QEC.png'))
            fig.savefig(os.path.join(folder, title+'.pdf'),
                    format='pdf',bbox_inches='tight')
        except:
            print 'Figure A has not been saved.'
    
    ### Return dictionairy (optional) ###
    if return_raw == True:
        return a, x, c0, u_c0, x_labels, folder
    
    elif return_dict == True:
            data_dict = {}
            # data_dict['a'] = a
            data_dict['x'] =x 
            data_dict['c0'] =c0 
            data_dict['u_c0'] =u_c0 
            data_dict['x_labels'] =x_labels 
            # data_dict['folder'] =folder 

            return data_dict, folder 
   

def QEC_create_data_dict(older_than = None, RO = 0, state = 'Z'):
    QEC_dict = {}
    k_dict = {}
    
    for error_sign in [1,-1]:
        # print 'sign_'+str(error_sign)
        QEC_dict[str(error_sign)] ={}
        for direction in ['positive','negative']:
            QEC_dict[str(error_sign)][direction] = {}

            for k in range(4):
                # print 'k_'+str(k)
                

                timestamp, folder = toolbox.latest_data(contains = direction+'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_sign)+'_'+state, older_than = older_than,return_timestamp = True)

                SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)
                # print SSRO_timestamp
                k_dict['k_'+str(k)] ={}
                k_dict['k_'+str(k)], folder = Plot_QEC(timestamp = timestamp, folder_name = folder,
                    ssro_calib_timestamp = SSRO_timestamp, return_raw = False, return_dict = True, post_select_QEC = True) 
                
                
            for item in k_dict['k_0']:
                QEC_dict[str(error_sign)][direction][item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item], k_dict['k_3'][item]), axis=0)

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

def openfile(name = 'QEC_141121_final.hdf5'):
   datafile = h5py.File(os.path.join(r'K:\ns\qt\Diamond\Projects\QEC LT\QEC data', name)) 
   return datafile

def closefile(datafile):
    datafile.close()


def save_data_hdf5file(datafile, data_dict, state, RO):
    f = datafile
    f_grp = f.create_group('state_'+state+'_RO_'+str(RO))

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



''' here you save new data '''

def QEC_data_single_state_RO(older_than = None,state = 'Z',RO = 0):

    QEC_data_dict = {}
    u_list = ['u_c0', 'u_c0_00','u_c0_01','u_c0_10','u_c0_11']
    c_list = ['c0', 'c0_00','c0_01','c0_10','c0_11']
    p_list = ['p00','p01','p10','p11']
    y_list = ['y','y_00','y_01','y_10','y_11']
    y_err_list = ['y_err','y_err_00','y_err_01','y_err_10','y_err_11']


    QEC_dict, folder = QEC_create_data_dict(older_than = older_than, RO = RO, state = state)
    for v in range(5):
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
    for p in range(4):
            QEC_data_dict[p_list[p]] = {}
            
            QEC_data_dict[p_list[p]] = (QEC_dict[str(-1)]['positive'][p_list[p]]+
                                                            QEC_dict[str(1)]['positive'][p_list[p]]+
                                                            QEC_dict[str(-1)]['negative'][p_list[p]]+
                                                            QEC_dict[str(1)]['negative'][p_list[p]])/4

    QEC_data_dict['x'] = QEC_dict[str(1)]['positive']['x']
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

def save_QEC_dataset(older_than = None, no_error = '11_1'):
    
    datafile = openfile(name = 'QEC_'+older_than[0:8]+'_error_syn_'+no_error+'.hdf5')
    QEC_temp_dict = {}

    for state in ['mY','Y','Z','mZ','X','mX']:
        QEC_temp_dict[state] = {}
        for RO  in range(7):
            QEC_temp_dict[state]['Tomo_'+str(RO)] = {}
            QEC_temp_dict[state]['Tomo_'+str(RO)], folder = QEC_data_single_state_RO(older_than = older_than, RO = RO, state = state)

            save_data_hdf5file(datafile,QEC_temp_dict[state]['Tomo_'+str(RO)], state, RO)
       
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



''' from here you can plot data taken from an existing HDF5 file '''

def QEC_plot_single_state_RO(date = '20141120', no_error = '00',state = 'Z',RO = 0, load_set = True, older_than = None):        
    
    
    if load_set == True:
        QEC_data_dict = {}   
        dataset_dict = load_QEC_dataset(date = date, no_error = no_error)
        QEC_data_dict  = dataset_dict[state]['Tomo_'+str(RO)]
        
    else:
        QEC_data_dict, folder =  QEC_data_single_state_RO(older_than = older_than,state = state,RO = RO)
    
    folder  = r'K:\ns\qt\Diamond\Projects\QEC LT\QEC data'   

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

    fig,ax = plt.subplots() 
    ax.errorbar(x,y,yerr=y_err,color = 'k' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title(date+'_errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_QEC')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')
    try:
        fig.savefig(
            os.path.join(folder,date+'_errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_all'+'.png'))
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
    ax.set_title(date+'_errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_QEC_PS')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')

    try:
        fig.savefig(
            os.path.join(folder,date+'_errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_ps'+'.png'))
    except:
        print 'Figure has not been saved.'




    fig,ax = plt.subplots()
    ax.set_title(str(folder)+'/'+ '\n probabilities')
    ax.plot(x,p_00, 'c', label = 'p00')
    ax.plot(x,p_01, 'k', label = 'p01')
    ax.plot(x,p_10, 'm', label = 'p10')
    ax.plot(x,p_11, 'b', label = 'p11')
    plt.legend()
    ax.set_xlabel('error probability')
    ax.set_ylabel('outcome probability')  
    ax.set_title(date+'_errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_QEC_probs')                
    
    try:
        fig.savefig(
            os.path.join(folder,date+'_errorsyn_'+no_error+'_state_'+state+'_RO_'+str(RO)+'_probs'+'.png'))
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
    ax.set_title(date+'_state_'+state+'_RO_'+str(RO)+'_noQEC')
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
        ax.set_title(date+'_errorsyn_'+no_error+'_state_'+state+'QEC_'+sum_type)
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


        # process_dict['dec_toff_'+y_list[v]] = 1/2.+1/16.*( dataset_dict['X']['Tomo_'+str(6)][y_list[v]] - dataset_dict['mX']['Tomo_'+str(6)][y_list[v]]
        #                 + dataset_dict['Z']['Tomo_'+str(0)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(0)][y_list[v]]
        #                 - dataset_dict['Y']['Tomo_'+str(5)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(5)][y_list[v]]
        #                 +dataset_dict['Z']['Tomo_'+str(1)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(1)][y_list[v]]
        #                 - dataset_dict['Y']['Tomo_'+str(6)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(6)][y_list[v]]
        #                 + dataset_dict['Z']['Tomo_'+str(3)][y_list[v]] - dataset_dict['mZ']['Tomo_'+str(3)][y_list[v]]
        #                 - dataset_dict['Y']['Tomo_'+str(4)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(4)][y_list[v]]
        #                 -dataset_dict['Z']['Tomo_'+str(6)][y_list[v]] + dataset_dict['mZ']['Tomo_'+str(6)][y_list[v]]
        #                 - dataset_dict['Y']['Tomo_'+str(3)][y_list[v]] + dataset_dict['mY']['Tomo_'+str(3)][y_list[v]]
        #                 )

        # process_dict['dec_toff_'+y_err_list[v]] = 1/16.*(dataset_dict['X']['Tomo_'+str(6)][y_list[v]]**2 + dataset_dict['mX']['Tomo_'+str(6)][y_list[v]]**2
        #                 + dataset_dict['Z']['Tomo_'+str(0)][y_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(0)][y_list[v]]**2
        #                 + dataset_dict['Y']['Tomo_'+str(5)][y_list[v]]**2 + dataset_dict['mY']['Tomo_'+str(5)][y_list[v]]**2
        #                 +dataset_dict['Z']['Tomo_'+str(1)][y_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(1)][y_list[v]]**2
        #                 + dataset_dict['Y']['Tomo_'+str(6)][y_list[v]]**2 + dataset_dict['mY']['Tomo_'+str(6)][y_list[v]]**2
        #                 + dataset_dict['Z']['Tomo_'+str(3)][y_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(3)][y_list[v]]**2
        #                 + dataset_dict['Y']['Tomo_'+str(4)][y_list[v]]**2 + dataset_dict['mY']['Tomo_'+str(4)][y_list[v]]**2
        #                 +dataset_dict['Z']['Tomo_'+str(6)][y_list[v]]**2 + dataset_dict['mZ']['Tomo_'+str(6)][y_list[v]]**2
        #                 + dataset_dict['Y']['Tomo_'+str(3)][y_list[v]]**2 + dataset_dict['mY']['Tomo_'+str(3)][y_list[v]]**2)**0.5


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