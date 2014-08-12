import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from analysis.lib.QEC import conditionalparity as CP
reload (CP)
from matplotlib import pyplot as plt

def BarPlotTomo(timestamp = None, measurement_name = ['adwindata'],folder_name ='Tomo',
        post_select = False,
        ssro_calib_timestamp =None, save = True,
        plot_fit = True, title ='' ) :
    '''
    Function that makes a bar plot with errorbars of MBI type data
    '''
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
    a.get_readout_results(name='adwindata',post_select = post_select)
    a.get_electron_ROC(ssro_calib_folder)

    x_labels = a.sweep_pts.reshape(-1)

    if post_select == False:

        y,y_err = a.convert_fidelity_to_contrast(a.p0,a.u_p0)
        x = range(len(y))

        if plot_fit ==True:
            fig,ax = plt.subplots()
            ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
            ax.set_xticks(x)
            if title == '':
                ax.set_title(str(folder)+'/'+str(timestamp))
            else:
                ax.set_title(title)
            print x_labels
            ax.set_xticklabels(x_labels.tolist())
            ax.set_ylim(-1,1)
            ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

        if save and ax != None:
            try:
                fig.savefig(
                    os.path.join(folder,'tomo.png'))
            except:
                print 'Figure has not been saved.'
    else: #if post_select == True:
        c0_0,u_c0_0 =  a.convert_fidelity_to_contrast(a.p0_0,a.u_p0_0)
        c0_1,u_c0_1 =  a.convert_fidelity_to_contrast(a.p0_1,a.u_p0_1)

        x = range(len(c_0))

        if plot_fit == True:
            fig_0,ax_0 = plt.subplots
            ax_0.bar(x,c0_0,yerr=u_c0_0,align ='center',ecolor = 'k' )
            ax_0.set_xticks(x)
            if title == '':
                ax_0.set_title(str(folder)+'/'+str(timestamp))
            else:
                ax_0.set_title(title)

            ax_0.set_xticklabels(x_labels.tolist())
            ax_0.set_ylim(-1,1)
            ax_0.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

            fig_1,ax_1 = plt.subplots
            ax_1.bar(x,c0_1,yerr=u_c0_1,align ='center',ecolor = 'k' )
            ax_1.set_xticks(x)
            if title == '':
                ax_1.set_title(str(folder)+'/'+str(timestamp))
            else:
                ax_1.set_title(title)
            ax_1.set_xticklabels(x_labels.tolist())
            ax_1.set_ylim(-1,1)
            ax_1.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')




def BarPlotTomoContrast(timestamps=None, measurement_name = ['adwindata'],folder_name ='Tomo',
        ssro_calib_timestamp =None, save = True,
        plot_fit = True) :
    '''
    Function that makes a bar plot with errorbars of MBI type data that has been measured with a positive
    and negative RO.


    !!!! PROBABLY DEPRECIATED, Remove if no longer used

    '''

    ### SSRO calibration
    if ssro_calib_timestamp == None:
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'

    ### Obtain and analyze data
    folder_a = toolbox.data_from_time(timestamps[0])
    a = mbi.MBIAnalysis(folder_a)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder)
    y_a= ((a.p0.reshape(-1)[:])-0.5)*2
    y_err_a = 2*a.u_p0.reshape(-1)[:]


    folder_b = toolbox.data_from_time(timestamps[1])
    b = mbi.MBIAnalysis(folder_b)
    b.get_sweep_pts()
    b.get_readout_results(name='adwindata')
    b.get_electron_ROC(ssro_calib_folder)
    y_b= ((b.p0.reshape(-1)[:])-0.5)*2
    y_err_b = 2*b.u_p0.reshape(-1)[:]

    x_labels = a.sweep_pts.reshape(-1)[:]
    x = range(len(y_a))

    ### Combine data
    y = (y_a - y_b)/2.
    y_err =  1./2*(y_err_a**2 + y_err_b**2)**0.5
    print y

    ### Fidelities
    # F_ZZ  = (1 + y[2] + y[5] + y[14])/4
    # F_ZmZ     = (1 + y[2] - y[5] - y[14])/4
    # F_ent     = (1 + y[0] -y[4] -y[8])/4
    F_ent   = (1 + y[0] +y[1] +y[2])/4
    # print 'Fidelity with ZZ  = ' + str(F_ZZ)
    # print 'Fidelity with ZmZ  = ' + str(F_ZmZ)
    print 'Fidelity with ent = ' + str(F_ent)



    if plot_fit ==True:
        fig,ax = plt.subplots()
        ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels.tolist())
        ax.set_ylim(-1.1,1.1)
        ax.set_title(str(folder_a)+'/'+str(timestamps[0]))
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

    if save and ax != None:
        try:
            fig.savefig(
                os.path.join(folder_a,'tomo.png'))
        except:
            print 'Figure has not been saved.'
