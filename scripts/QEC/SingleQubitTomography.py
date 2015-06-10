import numpy as np
import os
from analysis.lib.tools import toolbox
reload(toolbox)
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(fit)


def OneQubitTomo(timestamp=None, measurement_name = ['adwindata'], ssro_calib_timestamp =None,
        frequency = [1,1,1], offset =[ 0.5,0.5,0.5], amplitude =[ 0.5,0.5,0.5],  phase =[0.5,0.5,0.5],
        fixed = [],
        plot_fit = False, do_print = False, show_guess = True):
    '''
    Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed')
    [freq, offset, Amplitude, phase]
    '''
    # timestampZ = '20141023_174418'
    # folderZ = toolbox.data_from_time(timestampZ)
    # timestampX = '20141023_173547'
    # folderX = toolbox.data_from_time(timestampX)
    # timestampY = '20141023_173918'
    # folderY = toolbox.data_from_time(timestampY)
    
    if timestamp != None:
        timestampZ = timestamp
        folderZ = toolbox.data_from_time(timestamp)
    else:
        timestampZ, folderZ   = toolbox.latest_data('CarbonR',return_timestamp =True)
    timestampY, folderY = toolbox.latest_data('CarbonR',older_than = timestampZ, return_timestamp =True)
    timestampX, folderX = toolbox.latest_data('CarbonR',older_than = timestampY, return_timestamp =True)
    folders = [folderX,folderY, folderZ]

    if ssro_calib_timestamp == None:
        ssro_calib_folder = toolbox.latest_data('SSRO',older_than = timestampX)
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'

    order_of_bases = ['X', 'Y','Z']
    x = [None]*len(folders)
    y = [None]*len(folders)
    fit_result = [None]*len(folders)
    for i,folder in enumerate(folders):
        fit_results = []

        print '*'*60
        print order_of_bases[i] + ' Tomography'
        print 'folder %s' %folder


        for k in range(0,len(measurement_name)):
            a = mbi.MBIAnalysis(folder)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder)
            if i ==0:
                ax = a.plot_results_vs_sweepparam(ret='ax',labels = order_of_bases[i])
                # a.plot_results_vs_sweepparam(ax=ax)
            else:
                a.plot_results_vs_sweepparam(ax=ax,labels = order_of_bases[i])
            x[i] = a.sweep_pts.reshape(-1)[:]
            y[i]= a.p0.reshape(-1)[:]


            print frequency[i]
            p0, fitfunc, fitfunc_str = common.fit_cos(frequency[i], offset[i], amplitude [i],phase[i] )
            try:
                fit_result[i] = fit.fit1d(x[i],y[i], None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=fixed)
                if plot_fit == True:
                    plot.plot_fit1d(fit_result[i], np.linspace(x[i][0],x[i][-1],201), ax=ax,
                            plot_data=False,print_info = False)
                fit.write_to_file(fit_result[i],folder,fitname = str(order_of_bases[i])+'-tomography')
            except:
                pass
            if show_guess:
                ax.plot(np.linspace(x[i][0],x[i][-1],201), fitfunc(np.linspace(x[i][0],x[i][-1],201)), ':', lw=2)


    print 'fitfunction: '+fitfunc_str
    if plot_fit ==True:
        ax.legend(('X data','X-fit','Y data','Y-fit','Z data','Z-fit'),fontsize='x-small')
    elif plot_fit == False:
        ax.legend(('X data','Y data','Z data'),fontsize='small')

    ## plot data and fit as function of total time

    fit_results.append(fit_result[i])

    plt.savefig(os.path.join(folder, 'analyzed_tomography_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_tomography_result.png'),
    format='png')

    return fit_results

def BarPlotTomo(timestamp = None, measurement_name = ['adwindata'],folder_name ='Tomo',
        ssro_calib_timestamp =None, save = True,
        plot_fit = True) :
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
        ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil8'
        print ssro_calib_folder
    # <<<<<<< HEAD
        # ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
    # =======
        


    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder)

    x_labels = a.sweep_pts.reshape(-1)
    y= ((a.p0.reshape(-1))-0.5)*2
    x = range(len(y)) 
    y_err = 2*a.u_p0.reshape(-1)
    print 'y', y
    print 'err', y_err
    if plot_fit ==True: 
        fig,ax = plt.subplots() 
        rects = ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
        ax.set_xticks(x)
        # ax.title = timestamp
        print x_labels
        ax.set_xticklabels(x_labels.tolist())
        ax.set_ylim(-1.1,1.1)
        ax.set_title(str(folder)+'/'+str(timestamp))
        # ax.grid()
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')


            # print values on bar plot
    def autolabel(rects):
        for ii,rect in enumerate(rects):
            height = rect.get_height()
            plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2)) +'('+ str(int(round(y_err[ii]*100))) +')',
                ha='center', va='bottom')
    autolabel(rects)

    if save and ax != None:
        try:
            fig.savefig(
                os.path.join(folder,'tomo.png'))
        except:
            print 'Figure has not been saved.'

def BarPlotTomoContrast(timestamps = [None,None], tag = '', measurement_name = ['adwindata'],folder_name ='Tomo',
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
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'


    if timestamps[0] == None: 
        folder_a = toolbox.latest_data(contains='positive' + tag)
        folder_b = toolbox.latest_data(contains='negative' + tag)
    elif len(timestamps)==1:        
        folder_b = toolbox.data_from_time(timestamps[0])      
        print folder_b
        folder_a = toolbox.latest_data(contains = 'pos', older_than = timestamps[0])   
        print folder_a
    else:
        folder_a = toolbox.data_from_time(timestamps[0])      
        folder_b = toolbox.data_from_time(timestamps[1])           
    
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

    x_labels = a.sweep_pts.reshape(-1)[:]
    x = range(len(y_a)) 


    
    ### Combine data
    y = (y_a - y_b)/2.
    y_err =  1./2*(y_err_a**2 + y_err_b**2)**0.5 
    
    print tag
    print y
    print y_err
    # print folder_a
    # print folder_b

    # print y_a
    # print y_b
    # print y


    ### Fidelities
    # F_ZZ  = (1 + y[2] + y[5] + y[14])/4
    # F_ZmZ     = (1 + y[2] - y[5] - y[14])/4
    # F_ent     = (1 + y[0] -y[4] -y[8])/4
    # F_ent     = (1 + y[0] +y[1] +y[2])/4
    # print 'Fidelity with ZZ  = ' + str(F_ZZ)
    # print 'Fidelity with ZmZ  = ' + str(F_ZmZ)
    # print 'Fidelity with ent = ' + str(F_ent)

    # print 'XY = ' +str( (y[0]**2 + y[1]**2)**0.5)

    if plot_fit ==True: 
        fig,ax = plt.subplots() 
        rects = ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels.tolist())
        ax.set_ylim(-1.1,1.1)
        print 'test'
        print folder_a
        ax.set_title(str(folder_a)+'/'+str(timestamps[0]))
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

            # print values on bar plot
        def autolabel(rects):
            for ii,rect in enumerate(rects):
                height = rect.get_height()
                plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2)) +'('+ str(int(round(y_err[ii]*100))) +')',
                    ha='center', va='bottom')
        autolabel(rects)

    if save and ax != None:
        try:
            fig.savefig(
                os.path.join(folder_a,'tomo.png'))
        except:
            print 'Figure has not been saved.'

    if return_data == True:
        return x_labels, x, y, y_err


BarPlotTomoContrast(timestamps = ['20141219_123403','20141219_123526'], tag = '1_MBI', measurement_name = ['adwindata'],folder_name ='Tomo',
        ssro_calib_timestamp ='20141219_122724', save = True,
        plot_fit = True, return_data = False)        

BarPlotTomoContrast(timestamps = ['20141219_123647','20141219_123811'], tag = '2_MBI', measurement_name = ['adwindata'],folder_name ='Tomo',
        ssro_calib_timestamp ='20141219_122724', save = True,
        plot_fit = True, return_data = False)        

BarPlotTomoContrast(timestamps = ['20141219_123945','20141219_124114'], tag = '5_MBI', measurement_name = ['adwindata'],folder_name ='Tomo',
        ssro_calib_timestamp ='20141219_122724', save = True,
        plot_fit = True, return_data = False)        

plt.show()