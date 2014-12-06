import numpy as np
import os,sys

sys.path.append("/measuring/")
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import matplotlib.cm as cm
reload(common)
reload(plot)


def Carbon_gate_optimization(timestamp2=None, measurement_name = ['adwindata'], 
            A = [0.5, 0.5],
            fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = True):
    ''' Function to analyze data for optimization of the number of pulses for a controlled C13 gate. 
        The timestamp given should be for the no_pulse measurement
    '''
       
    ################################
    ### Measurements on Carbon 1 ###
    ################################
        #set1 (complete)
    t_stamp =  '20141128_222540' #older than time stamp to seperate different data sets
    tau_list = ['4.984', '4.986', '4.988', '4.99e-06', '4.992', '4.994', '4.996', '4.998', '5e-06', '5.002e-06', '5.004e-06' ]
    ssro_calib_folder = 'D:\\measuring\\data\\20141128\\082058_AdwinSSRO_SSROCalibration_111_1_sil18'
    pulse_list =['x','-x','y','-y']      
        #set2 (complete)
    t_stamp =  '20141130_041027' #older than time stamp to seperate different data sets
    tau_list = ['4.984', '4.986', '4.988', '4.99e-06', '4.992', '4.994', '4.996', '4.998', '5e-06', '5.002e-06', '5.004e-06' ]
    ssro_calib_folder = 'D:\\measuring\\data\\20141129\\184353_AdwinSSRO_SSROCalibration_111_1_sil18'
    pulse_list =['x','-x','y','-y']  
    ################################
    ### Measurements on Carbon 5 ###
    ################################


        # set 1 (Aborted due to crasd)
    # t_stamp =  '20141129_080859'
    # tau_list = ['8.916e-06', '8.918e-06', '8.92e-06', '8.922e-06', '8.924e-06', '8.926e-06', '8.928e-06']#, '8.93e-06', '8.932e-06', '8.934e-06', '8.934e-06' ]
    # ssro_calib_folder = 'D:\\measuring\\data\\20141128\\082058_AdwinSSRO_SSROCalibration_111_1_sil18'
    # pulse_list =['x','-x','y','-y']  

        # set 2 (Aborted due to corrupted data)
    # t_stamp = '20141129_182419'
    # tau_list = ['8.916e-06', '8.918e-06', '8.92e-06', '8.922e-06', '8.924e-06']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141129\\141457_AdwinSSRO_SSROCalibration_111_1_sil18'
    # pulse_list =['x','-x','y','-y']  

        # set 3 (Full)
    # t_stamp = '20141129_234902'
    # tau_list = ['8.916e-06', '8.918e-06', '8.92e-06', '8.922e-06', '8.924e-06', '8.926e-06', '8.928e-06', '8.93e-06', '8.932e-06', '8.934e-06', '8.936e-06' ]
    # # tau_list = ['8.924e-06', '8.926e-06', '8.928e-06', '8.93e-06']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141129\\184353_AdwinSSRO_SSROCalibration_111_1_sil18'
    # pulse_list =['x','-x','y','-y'] 

    ### Measurements on Carbon 2
        #set1 (complete)
    # t_stamp =  '20141130_112114' #older than time stamp to seperate different data sets
    # tau_list = ['1.0038e-05', '1.004e-05', '1.0042e-05', '1.0044e-05', '1.0046e-05', '1.0048e-05', '1.005e-05', 
    # '1.0052e-05', '1.0054e-05', '1.0056e-05', '1.0058e-05',' 1.006e-05', '1.0062e-05','1.0064e-05','1.0066e-05','1.0068e-05','1.007e-05','1.0072e-05',
    # '1.0074e-05','1.0076e-05','1.0078e-05']
    # tau_list = ['1.0058e-05', '1.005e-05']
    # ssro_calib_folder = 'D:\\measuring\data\\20141130\\114938_AdwinSSRO_SSROCalibration_111_1_sil18'
    # pulse_list =['x','-x','y','-y']     
        
        #set2 (smaller, for optimization)
    t_stamp =  '20141203_092830'
    tau_list = ['1.0048e-05', '1.005e-05', 
    '1.0052e-05', '1.0054e-05', '1.0056e-05', '1.0058e-05','1.006e-05', '1.0062e-05','1.0064e-05','1.0066e-05','1.0068e-05']
    tau_list = ['1.0048e-05', '1.005e-05', 
    '1.0052e-05', '1.0054e-05', '1.0056e-05', '1.0058e-05','1.006e-05']
    ssro_calib_folder = 'D:\\measuring\data\\20141203\\041410_AdwinSSRO_SSROCalibration_111_1_sil18'
    pulse_list =['x','-x','y','-y','no_pulse']    



    for ii, tau in enumerate(tau_list):

        if timestamp2 != None:
             folder1 = toolbox.data_from_time(timestamp2)
        else:
            timestamp, folder1 = toolbox.latest_data(tau, return_timestamp =True, older_than=t_stamp)

        folder_list=[folder1]

        print folder1

        for k in pulse_list:
            folder1 = toolbox.latest_data(contains ='sweep_N_'+ k +'_tau', older_than=timestamp)
            folder_list=folder_list+[folder1]

        y_all = []
        y_all_u = []
        for i,folder in enumerate(folder_list):
            for k in range(0,len(measurement_name)):
                a = mbi.MBIAnalysis(folder)
                a.get_sweep_pts()
                a.get_readout_results(name='adwindata')
                a.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)
                
                if i==0:
                    ax = a.plot_results_vs_sweepparam(ret='ax', markersize = 4, save=False, fmt = 'o-')
                if i!= 0: 
                    a.plot_results_vs_sweepparam(ax= ax,  markersize = 4, save=False, fmt = 'o-')
                
                ax.set_ylim(0.2,0.8)
                ax.set_xlim(a.sweep_pts[1],a.sweep_pts[-1])
                ax.axhspan(0,0.5,fill=False,ls='dotted')
                    

                x = a.sweep_pts.reshape(-1)[:]
                y = a.p0.reshape(-1)[:]
                y_u = a.u_p0.reshape(-1)[:]

                y_all.append(2*y-1)
                y_all_u.append(2*y_u)
               
        plt.savefig(os.path.join(folder_list[0], 'analyzed_result.pdf'),
        format='pdf')
        plt.savefig(os.path.join(folder_list[0], 'analyzed_result.png'),
        format='png')

        ### Analyze stokes vector length
        x_vec = (-y_all[1]+y_all[2])/2
        x_vec_u = ((y_all_u[1]**2+y_all_u[2]**2)**0.5)/2
        y_vec = (y_all[3]-y_all[4])/2
        y_vec_u = ((y_all_u[3]**2+y_all_u[4]**2)**0.5)/2
        z_vec = (y_all[0]-y_all[5])/2
        z_vec_u = ((y_all_u[0]**2+y_all_u[5]**2)**0.5)/2

        amplitude = (x_vec**2 + y_vec**2 + z_vec**2)**0.5
        amplitude_u =   (1./(x_vec**2 + y_vec**2 + z_vec**2)*(
                        x_vec**2 * x_vec_u**2 + 
                        y_vec**2 * y_vec_u**2 + 
                        z_vec**2 * z_vec_u**2))**0.5

        if ii == 0:
            fig = a.default_fig(figsize=(7.5,5))
            ax2 = a.default_ax(fig)
            start, end = ax2.get_xlim()
            ax2.axhspan(0,1,fill=False,ls='dotted')
            ax2.set_ylim(-0.05,0.5)
            ax2.set_xlim(x[1]-1,x[-1]+1)
            colors = cm.rainbow(np.linspace(0, 1, len(tau_list)))

        ax2.errorbar(x[1::], amplitude[1::], amplitude_u[1::], 0*np.ones(len(amplitude[1::])), '.-', lw=1, color = colors[ii],
                label = 'tau = ' +str(tau) + '_min = ' + str(np.around(amplitude[np.argmin(amplitude)],3)) + ' +/- ' + str(np.around(amplitude_u[np.argmin(amplitude)],3)) + ' at N = ' + str(x[np.argmin(amplitude)])) #N = 16
        ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
      
        # print 'minimum = ' + str(amplitude[np.argmin(amplitude)]) + 'at N = ' + str(x[np.argmin(amplitude)])

        x_av    = (y_all[1]+y_all[2])/2
        x_av_u  = ((y_all_u[1]**2+y_all_u[2]**2)**0.5)/2
        y_av    = (y_all[3]+y_all[4])/2
        y_av_u  = ((y_all_u[3]**2+y_all_u[4]**2)**0.5)/2
        z_av    = y_all[0]
        z_av_u  = y_all_u[0]


        fig3 = a.default_fig(figsize=(5,5))
        ax3 = a.default_ax(fig3)
        start, end = ax3.get_xlim()
        ax3.set_xlim(x[0]-1,x[-1]+1)
        ax3.errorbar(x, x_av, x_av_u,0*np.ones(len(amplitude)),  '.-r', lw=1,label = 'x') #N = 16
        ax3.errorbar(x, y_av, x_av_u,0*np.ones(len(amplitude)), '.-b', lw=1,label = 'y') #N = 16
        ax3.errorbar(x, z_av, x_av_u,0*np.ones(len(amplitude)), '.-g', lw=1,label = 'z') #N = 16
           
        ax3.legend()




    # diff = np.abs(y - 0.5)
    # print diff
    # print 'Optimum number of pulses N = ' + str(x[np.argmin(diff)])
    # print 'with y-0.5 = ' + str(y[np.argmin(diff)]-0.5) + ' +/- ' + str(y_u[np.argmin(diff)])

            # freq = fit_results[0]['params_dict']['f1']
            # period = 1/freq 
            # print 'Period is %s pulses ' %(period)
            # # N_pi = round(period*.5/2)*2.0 
            # N_pi2 = round(period/2*.25)*2.0
            # # print 'Pi pulse: %s pulses' %N_pi
            # print 'Pi2 pulse: %s pulses' %N_pi2
    # return fit_results


