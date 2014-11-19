import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from analysis.lib.QEC import ConditionalParity as CP
reload (CP)
from matplotlib import pyplot as plt
script_name = 'three_qubit_QEC_analysis.py'


def Plot_QEC(timestamp = None, measurement_name = ['adwindata'],folder_name ='QEC',
        plot_post_select = False,
        ssro_calib_timestamp = None, save = True,
        do_plots = False, title =None ,fontsize = 10, post_select_QEC = False, return_raw = True, return_dict = False) :
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
        ax.errorbar(x,c0,yerr=u_c0,ecolor = 'k' )
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

    ''' postselected data '''

    a = CP.ConditionalParityAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata',post_select_QEC = True)
    a.get_electron_ROC(ssro_calib_folder, post_select_QEC = True)
            
    c0_00,u_c0_00 =  a.convert_fidelity_to_contrast(a.p0_00,a.u_p0_00)
    c0_01,u_c0_01 =  a.convert_fidelity_to_contrast(a.p0_01,a.u_p0_01)
    c0_10,u_c0_10 =  a.convert_fidelity_to_contrast(a.p0_10,a.u_p0_10)
    c0_11,u_c0_11 =  a.convert_fidelity_to_contrast(a.p0_11,a.u_p0_11)

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
        rects = ax.errorbar(x,y,yerr=y_err,ecolor = 'k' )
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
        ax.errorbar(x,y_00,yerr=y_00_err,ecolor = 'c', label = 'y_00' )
        ax.errorbar(x,y_01,yerr=y_01_err,ecolor = 'k', label = 'y_01' )
        ax.errorbar(x,y_10,yerr=y_10_err,ecolor = 'm', label = 'y_10' )
        ax.errorbar(x,y_11,yerr=y_11_err,ecolor = 'b', label = 'y_11' )
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


def Contrast_Plot_QEC_yesno(timestamp = None , measurement_name = ['adwindata'],folder_name ='QEC',
        post_select_QEC = False, ssro_calib_timestamp =None, do_plot = False, return_data = False,return_timestamp = False) :

    '''
    Function that makes a plot with errorbars of MBI type data that has been measured with a positive
    and negative RO.
    '''       
    timestamp_pos_n, folder_a = toolbox.latest_data(contains = 'positive_', older_than = timestamp,return_timestamp = True)
    timestamp_neg_n, folder_b = toolbox.latest_data(contains = 'negative_', older_than = timestamp,return_timestamp = True)
    # print folder_a
    # print folder_b
    timestamp_pos_y, folder_a = toolbox.latest_data(contains = 'positive_', older_than = timestamp_pos_n,return_timestamp = True)
    timestamp_neg_y, folder_b = toolbox.latest_data(contains = 'negative_', older_than = timestamp_pos_n,return_timestamp = True)
    # print folder_a
    # print folder_b
    a_p_n, x, c0_p_n, u_c0_p_n, c0_00_p_n, u_c0_00_p_n, c0_01_p_n, u_c0_01_p_n, c0_10_p_n, u_c0_10_p_n, c0_11_p_n, u_c0_11_p_n, x_labels, folder_p_n = Plot_QEC(timestamp = timestamp_pos_n, 
            measurement_name = measurement_name, folder_name = 'positive',
            ssro_calib_timestamp = ssro_calib_timestamp) 

    a_n_n, x, c0_n_n, u_c0_n_n, c0_00_n_n, u_c0_00_n_n, c0_01_n_n, u_c0_01_n_n, c0_10_n_n, u_c0_10_n_n, c0_11_n_n, u_c0_11_n_n, x_labels, folder_n_n = Plot_QEC(timestamp = timestamp_neg_n, 
            measurement_name = measurement_name, folder_name = 'negative',
            ssro_calib_timestamp =ssro_calib_timestamp) 
    
    a_p_y, x, c0_p_y, u_c0_p_y, c0_00_p_y, u_c0_00_p_y, c0_01_p_y, u_c0_01_p_y, c0_10_p_y, u_c0_10_p_y, c0_11_p_y, u_c0_11_p_y, x_labels, folder_p_y = Plot_QEC(timestamp =timestamp_pos_y, 
            measurement_name = measurement_name, folder_name = 'positive',
            ssro_calib_timestamp = ssro_calib_timestamp) 

    a_n_y, x, c0_n_y, u_c0_n_y, c0_00_n_y, u_c0_00_n_y, c0_01_n_y, u_c0_01_n_y, c0_10_n_y, u_c0_10_n_y, c0_11_n_y, u_c0_11_n_y, x_labels, folder_n_y = Plot_QEC(timestamp = timestamp_neg_y, 
            measurement_name = measurement_name, folder_name = 'negative',
            ssro_calib_timestamp =ssro_calib_timestamp)     
    ### Combine data

        ## all data
    y_y = (c0_p_y - c0_n_y)/2.
    y_err_y =  1./2*(u_c0_p_y**2 + u_c0_n_y**2)**0.5    
    y_n = (c0_p_n - c0_n_n)/2.
    y_err_n =  1./2*(u_c0_p_n**2 + u_c0_n_n**2)**0.5

    if do_plot == True:
        fig,ax = plt.subplots() 
        ax.errorbar(x,y_y,yerr=y_err_y,ecolor = 'k' )
        ax.errorbar(x,y_n,yerr=y_err_n,ecolor = 'r' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title(str(folder_a)+'/'+'\n' + script_name)
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')

        try:
            fig.savefig(
                os.path.join(folder_a,'QEC.png'))
        except:
            print 'Figure has not been saved.'

        ## ms=0 data
    y_00_y = (c0_00_p_y - c0_00_n_y)/2.
    y_err_00_y =  1./2*(u_c0_00_p_y**2 + u_c0_00_n_y**2)**0.5    
    y_00_n = (c0_00_p_n - c0_00_n_n)/2.
    y_err_00_n =  1./2*(u_c0_00_p_n**2 + u_c0_00_n_n**2)**0.5

    if do_plot == True:
        fig,ax = plt.subplots() 
        ax.errorbar(x,y_00_y,yerr=y_err_00_y,ecolor = 'k' )
        ax.errorbar(x,y_00_n,yerr=y_err_00_n,ecolor = 'r' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title(str(folder_a)+'/'+'\n postselect_00')
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')


      

        try:
            fig.savefig(
                os.path.join(folder_a,'QEC_ps_00.png'))
        except:
            print 'Figure has not been saved.'


        ## ms=1 data
    y_01_y = (c0_01_p_y - c0_01_n_y)/2.
    y_err_01_y =  1./2*(u_c0_01_p_y**2 + u_c0_01_n_y**2)**0.5    
    y_01_n = (c0_01_p_n - c0_01_n_n)/2.
    y_err_01_n =  1./2*(u_c0_01_p_n**2 + u_c0_01_n_n**2)**0.5

    if do_plot == True:
        fig,ax = plt.subplots() 
        ax.errorbar(x,y_01_y,yerr=y_err_01_y,ecolor = 'k' )
        ax.errorbar(x,y_01_n,yerr=y_err_01_n,ecolor = 'r' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title(str(folder_a)+'/'+ '\n postselect_01')
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')

        try:
            fig.savefig(
                os.path.join(folder_a,'QEC_ps_01.png'))
        except:
            print 'Figure has not been saved.'


        ## ms=0 data
    y_10_y = (c0_10_p_y - c0_10_n_y)/2.
    y_err_10_y =  1./2*(u_c0_10_p_y**2 + u_c0_10_n_y**2)**0.5    
    y_10_n = (c0_10_p_n - c0_10_n_n)/2.
    y_err_10_n =  1./2*(u_c0_10_p_n**2 + u_c0_10_n_n**2)**0.5

    if do_plot == True:
        fig,ax = plt.subplots() 
        ax.errorbar(x,y_10_y,yerr=y_err_10_y,ecolor = 'k' )
        ax.errorbar(x,y_10_n,yerr=y_err_10_n,ecolor = 'r' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)

        ax.set_title(str(folder_a)+'/'+'\n postselect_10')
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')

        try:
            fig.savefig(
                os.path.join(folder_a,'QEC_ps_10.png'))
        except:
            print 'Figure has not been saved.'


        ## ms=1 data
    y_11_y = (c0_11_p_y - c0_11_n_y)/2.
    y_err_11_y =  1./2*(u_c0_11_p_y**2 + u_c0_11_n_y**2)**0.5    
    y_11_n = (c0_11_p_n - c0_11_n_n)/2.
    y_err_11_n =  1./2*(u_c0_11_p_n**2 + u_c0_11_n_n**2)**0.5

    if do_plot == True:
        fig,ax = plt.subplots() 
        ax.errorbar(x,y_11_y,yerr=y_err_11_y,ecolor = 'k' )
        ax.errorbar(x,y_11_n,yerr=y_err_11_n,ecolor = 'r' )
        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-0.1,1.1)
        ax.set_title(str(folder_a)+'/'+ '\n postselect_11')
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        ax.set_xlabel('error probability')
        ax.set_ylabel('Contrast')

        try:
            fig.savefig(
                os.path.join(folder_a,'QEC_ps_11.png'))
        except:
            print 'Figure has not been saved.'

    # process the probabilities

    # p00_avg = (a_p.p00+a_n.p00)/2
    # p01_avg = (a_p.p01+a_n.p01)/2
    # p10_avg = (a_p.p10+a_n.p10)/2
    # p11_avg = (a_p.p11+a_n.p11)/2
    # if do_plot == True:
    #     fig,ax = plt.subplots()
    #     ax.set_title(str(folder_a)+'/'+ '\n probabilities')
    #     ax.plot(x,p00_avg, 'c', label = 'p00')
    #     ax.plot(x,p01_avg, 'k', label = 'p01')
    #     ax.plot(x,p10_avg, 'm', label = 'p10')
    #     ax.plot(x,p11_avg, 'b', label = 'p11')
    #     plt.legend()
    #     ax.set_xlabel('error probability')
    #     ax.set_ylabel('outcome probability')

        # try:
        #     fig.savefig(
        #         os.path.join(folder_a,'QEC_probabilities.png'))
        # except:
        #     print 'Figure has not been saved.'

    # if return_data == True:
    #     return x, y, y_err, y_00, y_00_err, p00_avg, y_01, y_01_err, p01_avg, y_10, y_10_err, p10_avg, y_11, y_11_err, p11_avg
    if return_timestamp == True:
        return timestamp_pos_y





def Contrast_Plot_QEC_loop(older_than = None, newer_than = None, ssro_calib_timestamp =None) :

    timestamp = older_than

    while timestamp != newer_than:
        timestamp = Contrast_Plot_QEC_yesno(timestamp = timestamp , measurement_name = ['adwindata'],folder_name ='QEC',
            post_select_QEC = False, ssro_calib_timestamp =ssro_calib_timestamp, do_plot = True, return_data = False,return_timestamp = True)




def Contrast_Plot_QEC_full(timestamp = None, measurement_name = ['adwindata'],folder_name ='RO1',
        ssro_calib_timestamp =None, save = True,
        do_plot  = True):
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
        rects = ax.errorbar(x,y,yerr=y_err,ecolor = 'k' )
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
        rects = ax.errorbar(x,y_00,yerr=y_00_err,ecolor = 'k' )
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
        rects = ax.errorbar(x,y_01,yerr=y_01_err,ecolor = 'k' )
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
        rects = ax.errorbar(x,y_10,yerr=y_10_err,ecolor = 'k' )
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
        rects = ax.errorbar(x,y_11,yerr=y_11_err,ecolor = 'k' )
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



''' analysis for no QEC'''


def Plot_errorcurve_no_QEC(timestamp = None, measurement_name = ['adwindata'],folder_name ='QEC',
        ssro_calib_timestamp =None, save = True,
        plot_fit = True, return_data = False) :
    '''
    Function that makes a bar plot with errorbars of MBI type data that has been measured with a positive
    and negative RO.
    '''

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
        ax.errorbar(x,y,yerr=y_err,ecolor = 'k' )
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
        return x_labels, x, y, y_err


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

                k_dict['k_'+str(k)] ={}
                k_dict['k_'+str(k)], folder = Plot_QEC(timestamp = timestamp, folder_name = folder,
                    ssro_calib_timestamp = SSRO_timestamp, return_raw = False, return_dict = True) 
                
                
            for item in k_dict['k_0']:
                QEC_dict[str(error_sign)][direction][item] = np.concatenate((k_dict['k_0'][item],k_dict['k_1'][item],k_dict['k_2'][item], k_dict['k_3'][item]), axis=0)

    return QEC_dict,folder

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
                                                        QEC_dict[str(1)]['negative'][u_list[v]]**2)/4
    for p in range(4):
            QEC_data_dict[p_list[p]] = {}
            
            QEC_data_dict[p_list[p]] = (QEC_dict[str(-1)]['positive'][p_list[p]]+
                                                            QEC_dict[str(1)]['positive'][p_list[p]]+
                                                            QEC_dict[str(-1)]['negative'][p_list[p]]+
                                                            QEC_dict[str(1)]['negative'][p_list[p]])/4

    QEC_data_dict['x'] = QEC_dict[str(1)]['positive']['x']

    return QEC_data_dict, folder

def QEC_plot_single_state_RO(older_than = None,state = 'Z',RO = 0):        
   

    QEC_data_dict, folder = QEC_data_single_state_RO(older_than = older_than, RO = RO, state = state)

    print folder
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
    y_err_11 = QEC_data_dict['y_11']

    p_00 = QEC_data_dict['p00']
    p_01 = QEC_data_dict['p01']
    p_10 = QEC_data_dict['p10']
    p_11 = QEC_data_dict['p11']

    fig,ax = plt.subplots() 
    ax.errorbar(x,y,yerr=y_err,ecolor = 'k' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.set_title(str(folder)+'/'+'\n' + script_name)
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')
    try:
        fig.savefig(
            os.path.join(folder,'QEC.png'))
    except:
        print 'Figure has not been saved.'

    fig,ax = plt.subplots() 
    ax.errorbar(x,y_00,yerr=y_err_00,ecolor = 'c', label = 'y_00' )
    ax.errorbar(x,y_01,yerr=y_err_01,ecolor = 'k', label = 'y_01' )
    ax.errorbar(x,y_10,yerr=y_err_10,ecolor = 'm', label = 'y_10' )
    ax.errorbar(x,y_11,yerr=y_err_11,ecolor = 'b', label = 'y_11' )
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-0.1,1.1)
    plt.legend()
    ax.set_title(str(folder)+'/'+'\n postselect')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_xlabel('error probability')
    ax.set_ylabel('Contrast')

    try:
        fig.savefig(
            os.path.join(folder,'QEC_ps_full.png'))
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
    
    try:
        fig.savefig(
            os.path.join(folder,'QEC_probs.png'))
    except:
        print 'Figure has not been saved.'


def plot_all_state_Z(older_than = None):        
   
    for RO  in range(7):
        QEC_state_dict['Tomo_'+RO] = {}
        QEC_state_dict['Tomo_'+RO], folder = QEC_data_single_state_RO(older_than = older_than, RO = RO, state = state)

    # three_qubit_fid = 1/8.*(np.ones(len(QEC_state_dict['Tomo_0']))+QEC_state_dict['Tomo_'+RO]+QEC_state_dict['Tomo_'+RO]+QEC_state_dict['Tomo_'+RO]