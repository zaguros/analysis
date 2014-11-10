import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from analysis.lib.QEC import ConditionalParity as CP
reload (CP)
from matplotlib import pyplot as plt

def Plot_QEC(timestamp = None, measurement_name = ['adwindata'],folder_name ='Tomo',
        plot_post_select = False,
        ssro_calib_timestamp = None, save = True,
        do_plots = False, title =None ,fontsize = 10, post_select_QEC = False) :
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
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'

    a = CP.ConditionalParityAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata',post_select_QEC = False)
    a.get_electron_ROC(ssro_calib_folder)

    x_labels = a.sweep_pts.reshape(-1)

    ''' all data '''

    c0,u_c0 = a.convert_fidelity_to_contrast(a.p0,a.u_p0)
    x = range(len(c0))

    if do_plots ==True:
        fig,ax = plt.subplots()
        ax.plot(x,c0,yerr=u_c0,align ='center',ecolor = 'k' )
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

    x = range(len(c0_00))

    # if do_plots == True and plot_post_select == True:
        


    #     fig_00,ax_00 = plt.subplots()
    #     rects = ax_00.plot(x,c0_00,yerr=u_c0_00,align ='center',ecolor = 'k' )
    #     # ax_00.bar(x,a.p0_00,yerr = u_c0_00,align = 'center',ecolor = 'k')
    #     ax_00.set_xticks(x)
    #     if title == None:
    #         ax_00.set_title(str(folder)+'/'+str(timestamp)+'0')
    #     else:
    #         ax_00.set_title(str(title)+'0')

    #     ax_00.set_xticklabels(x_labels.tolist())
    #     # ax_00.set_ylim(-1,1)
    #     ax_00.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')


    #     fig_01,ax_01 = plt.subplots()
    #     rects = ax_01.bar(x,c0_01,yerr=u_c0_01,align ='center',ecolor = 'k' )
    #     ax_01.set_xticks(x)
    #     if title == None:
    #         ax_01.set_title(str(folder)+'/'+str(timestamp)+'1')
    #     else:
    #         ax_01.set_title(str(title)+'1')
    #     ax_01.set_xticklabels(x_labels.tolist())
    #     # ax_01.set_ylim(-1,1)
    #     ax_01.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')


    #         # print values on bar plot
    #     def autolabel(rects):
    #         for ii,rect in enumerate(rects):
    #             height = rect.get_height()
    #             plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(c0_01[ii],2)) +'('+ str(int(round(u_c0_01[ii]*100))) +')',
    #                 ha='center', va='bottom')
    #     autolabel(rects)

    #     fig_10,ax_10 = plt.subplots()
    #     rects = ax_10.bar(x,c0_10,yerr=u_c0_10,align ='center',ecolor = 'k' )
    #     # ax_10.bar(x,a.p0_10,yerr = u_c0_10,align = 'center',ecolor = 'k')
    #     ax_10.set_xticks(x)
    #     if title == None:
    #         ax_10.set_title(str(folder)+'/'+str(timestamp)+'0')
    #     else:
    #         ax_10.set_title(str(title)+'0')

    #     ax_10.set_xticklabels(x_labels.tolist())
    #     # ax_10.set_ylim(-1,1)
    #     ax_10.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')


    #         # print values on bar plot
    #     def autolabel(rects):
    #         for ii,rect in enumerate(rects):
    #             height = rect.get_height()
    #             plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(c0_10[ii],2)) +'('+ str(int(round(u_c0_10[ii]*100))) +')',
    #                 ha='center', va='bottom')
    #     autolabel(rects)

    #     fig_11,ax_11 = plt.subplots()
    #     rects = ax_11.bar(x,c0_11,yerr=u_c0_11,align ='center',ecolor = 'k' )
    #     ax_11.set_xticks(x)
    #     if title == None:
    #         ax_11.set_title(str(folder)+'/'+str(timestamp)+'1')
    #     else:
    #         ax_11.set_title(str(title)+'1')
    #     ax_11.set_xticklabels(x_labels.tolist())
    #     # ax_11.set_ylim(-1,1)
    #     ax_11.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')


    #         # print values on bar plot
    #     def autolabel(rects):
    #         for ii,rect in enumerate(rects):
    #             height = rect.get_height()
    #             plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(c0_11[ii],2)) +'('+ str(int(round(u_c0_11[ii]*100))) +')',
    #                 ha='center', va='bottom')
    #     autolabel(rects)



        # try:
        #     fig_00.savefig(os.path.join(folder, str(title)+'0.pdf'),
        #             format='pdf',bbox_inches='tight')

        #     fig_01.savefig(os.path.join(folder, str(title)+'1.pdf'),
        #             format='pdf',bbox_inches='tight')
        #     fig_10.savefig(os.path.join(folder, str(title)+'0.pdf'),
        #             format='pdf',bbox_inches='tight')

        #     fig_11.savefig(os.path.join(folder, str(title)+'1.pdf'),
        #             format='pdf',bbox_inches='tight')
        # except:
        #     print 'Figure B has not been saved.'

    return x, c0, u_c0, c0_00, u_c0_00, c0_01, u_c0_01, c0_10, u_c0_10, c0_11, u_c0_11, x_labels, folder



def Contrast_Plot_QEC(timestamps=[None, None], measurement_name = ['adwindata'],folder_name ='Tomo',
        post_select_QEC = False, ssro_calib_timestamp =None) :

    '''
    Function that makes a plot with errorbars of MBI type data that has been measured with a positive
    and negative RO.
    '''
    x, c0_p, u_c0_p, c0_00_p, u_c0_00_p, c0_01_p, u_c0_01_p, c0_10_p, u_c0_10_p, c0_11_p, u_c0_11_p, x_labels, folder_p = Plot_QEC(timestamp = timestamps[0], 
            measurement_name = measurement_name, folder_name = 'positive',
            ssro_calib_timestamp = ssro_calib_timestamp) 

    x, c0_n, u_c0_n, c0_00_n, u_c0_00_n, c0_01_n, u_c0_01_n, c0_10_n, u_c0_10_n, c0_11_n, u_c0_11_n, x_labels, folder_n = Plot_QEC(timestamp = timestamps[1], 
            measurement_name = measurement_name, folder_name = 'negative',
            ssro_calib_timestamp =ssro_calib_timestamp) 
        
    ### Combine data

        ## all data
    y = (c0_p - c0_n)/2.
    y_err =  1./2*(u_c0_p**2 + u_c0_n**2)**0.5      

    fig,ax = plt.subplots() 
    rects = ax.plot(x,y,yerr=y_err,align ='center',ecolor = 'k' )
    ax.set_ylim(-1.1,1.1)
    ax.set_title(str(folder_p)+'/'+'\n' + script_name)
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
   

    try:
        fig.savefig(
            os.path.join(folder_p,'QEC.png'))
    except:
        print 'Figure has not been saved.'

        ## ms=0 data
    y_00 = (c0_00_p - c0_00_n)/2.
    y_00_err =  1./2*(u_c0_00_p**2 + u_c0_00_n**2)**0.5 

    fig,ax = plt.subplots() 
    rects = ax.plot(x,y_00,yerr=y_00_err,align ='center',ecolor = 'k' )
    ax.set_ylim(-1.1,1.1)
    ax.set_title(str(folder_p)+'/'+'\n postselect_00')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

  

    try:
        fig.savefig(
            os.path.join(folder_p,'QEC_ps_00.png'))
    except:
        print 'Figure has not been saved.'


        ## ms=1 data
    y_01 = (c0_01_p - c0_01_n)/2.
    y_01_err =  1./2*(u_c0_01_p**2 + u_c0_01_n**2)**0.5 

    fig,ax = plt.subplots() 
    rects = ax.plot(x,y_01,yerr=y_01_err,align ='center',ecolor = 'k' )
    ax.set_ylim(-1.1,1.1)
    ax.set_title(str(folder_p)+'/'+ '\n postselect_01')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

    try:
        fig.savefig(
            os.path.join(folder_p,'QEC_ps_01.png'))
    except:
        print 'Figure has not been saved.'


        ## ms=0 data
    y_10 = (c0_10_p - c0_10_n)/2.
    y_10_err =  1./2*(u_c0_10_p**2 + u_c0_10_n**2)**0.5 

    fig,ax = plt.subplots() 
    rects = ax.plot(x,y_10,yerr=y_10_err,align ='center',ecolor = 'k' )
    ax.set_ylim(-1.1,1.1)
    ax.set_title(str(folder_p)+'/'+'\n postselect_10')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

    try:
        fig.savefig(
            os.path.join(folder_p,'QEC_ps_10.png'))
    except:
        print 'Figure has not been saved.'


        ## ms=1 data
    y_11 = (c0_11_p - c0_11_n)/2.
    y_11_err =  1./2*(u_c0_11_p**2 + u_c0_11_n**2)**0.5 

    fig,ax = plt.subplots() 
    rects = ax.plot(x,y_11,yerr=y_11_err,align ='center',ecolor = 'k' )
    ax.set_ylim(-1.1,1.1)
    ax.set_title(str(folder_p)+'/'+ '\n postselect_11')
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

    try:
        fig.savefig(
            os.path.join(folder_p,'QEC_ps_11.png'))
    except:
        print 'Figure has not been saved.'