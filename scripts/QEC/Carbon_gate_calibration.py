import numpy as np
import os,sys
# import qt

sys.path.append("/measuring/")
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import matplotlib.cm as cm
reload(common)
reload(plot)


# SAMPLE = qt.exp_params['samples']['current']
# SAMPLE_CFG = qt.exp_params['protocols']['current']

def Carbon_gate_optimization_XY(measurement_name = ['adwindata'], 
            A = [0.5, 0.5],
            fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = True):
    ''' Function to analyze data for optimization of the number of pulses for a controlled C13 gate. 
        The timestamp given should be for the no_pulse measurement
    '''
       
    ################################
    ### Measurements on Carbon 1 ###
    ################################

    # t_stamp =  '20141208_040429'    #older than time stamp to seperate different data sets
    # tau_list = ['4.988e-06', '4.99e-06', '4.992e-06', '4.994e-06', '4.996e-06', '4.998e-06', '5e-06', '5.002e-06', '5.004e-06' ]
    # tau_list = ['4.992e-06', '4.994e-06', '4.996e-06', '4.998e-06']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141207\\223337_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(32,51,2)


    ### other resonances for C1??
    
    t_stamp =  '20150329_120000'
    ssro_calib_folder = 'D:\\measuring\\data\\20150329\\100452_AdwinSSRO_SSROCalibration_111_1_sil18'
    

    # tau0 = 7.22e-6
    # tau_range = 8e-9
    # tau_list = np.arange(tau0 - tau_range,tau0 + tau_range,2e-9)
    # N_list   = np.arange(34,60,2)

    
    tau0 = 9.436e-6
    tau_range = 8e-9
    tau_list = np.arange(tau0 - tau_range,tau0 + tau_range,2e-9)
    N_list   = np.arange(34,60,2)
    

    ################################
    ### Measurements on Carbon 5 ###
    ################################
        ### Measurment set 1 = Bugged!
    # t_stamp =  '20141208_102528'    #older than time stamp to seperate different data sets
    # tau_list = ['8.918e-06', '8.92e-06', '8.922e-06', '8.924e-06', '8.926e-06', '8.928e-06', '8.93e-06', '8.932e-06', '8.934e-06']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141207\\223337_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(34,55,2)

        ### Measurement set 2
    # t_stamp =  '20141210_080551'    #older than time stamp to seperate different data sets
    # tau_list = ['8.916e-06', '8.918e-06', '8.92e-06', '8.922e-06', '8.924e-06', '8.926e-06', '8.928e-06', 
    #             '8.93e-06', '8.932e-06', '8.934e-06', '8.936e-06']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141210\\071632_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(34,55,2)

    # t_stamp =  '20141214_103403'    #older than time stamp to seperate different data sets
    # tau_list = ['6.524e-06', '6.526e-06', '6.528e-06', '6.53e-06', '6.532e-06', '6.534e-06', '6.536e-06', '6.538e-06',
    # '6.54e-06', '6.542e-06', '6.544e-06','6.548e-06','6.55e-06','6.552e-06','6.554e-06',
    # '6.556e-06', '6.558e-06', '6.56e-06','6.562e-06','6.564e-06']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141210\\071632_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(28,49,2)

    # t_stamp =  '20141215_083533'    #older than time stamp to seperate different data sets
    # tau_list = ['6.53e-06', '6.532e-06', '6.534e-06', '6.536e-06', '6.538e-06',
    # '6.54e-06', '6.542e-06', '6.544e-06','6.548e-06','6.55e-06','6.552e-06','6.554e-06',
    # '6.556e-06', '6.558e-06']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141210\\071632_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(20,35,2)

    


    # t_stamp =  '20150226_223255'    #older than time stamp to seperate different data sets
    # ssro_calib_folder = 'D:\\measuring\\data\\20150225\\231447_AdwinSSRO_SSROCalibration_111_1_sil18'
    # tau_list = np.arange(9.652e-6, 9.658e-6, 2e-9)
    
    # tau_list = np.r_[tau_list, np.arange(9.678e-6, 9.684e-6, 2e-9), np.arange(9.671e-6,9.6735e-6,0.5e-9)]
    # N_list   = np.arange(9, 16, 1)

    # t_stamp =  '20141217_065954'    #older than time stamp to seperate different data sets
    # tau_list = ['1.1288e-05', '1.129e-05', '1.1292e-05', '1.1294e-05',
    #      '1.1296e-05',   '1.1298e-05',   '1.13e-05',
    #      '1.1302e-05',   '1.1304e-05',   '1.1306e-05',
    #      '1.1308e-05',   '1.131e-05',   '1.1312e-05',
    #      '1.1314e-05',   '1.1316e-05']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141216\\040057_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(32,53,2)



    

    ################################
    ### Measurements on Carbon 2 ###
    ###############################
        ### Measurment set 1: Bugged!
    # t_stamp =  '20141208_140824'    #older than time stamp to seperate different data sets
    # tau_list = ['1.0048e-05', '1.005e-05', '1.0052e-05', '1.0054e-05', '1.0056e-05', '1.0058e-05','1.006e-05', '1.0062e-05','1.0064e-05','1.0066e-05']
    # tau_list = ['1.0054e-05', '1.0056e-05', '1.0058e-05','1.006e-05', '1.0062e-05']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141207\\223337_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(16,25,2)

        ### Measurement set 2
    # t_stamp =  '20141209_203041'    #older than time stamp to seperate different data sets
    # tau_list = ['1.0044e-05', '1.0046e-05','1.0048e-05', '1.005e-05', '1.0052e-05', '1.0054e-05', '1.0056e-05', 
    #         '1.0058e-05','1.006e-05', '1.0062e-05','1.0064e-05','1.0066e-05', '1.0068e-05', '1.007e-05', '1.0072e-05']
    # tau_list = ['1.0056e-05', '1.0058e-05','1.006e-05', '1.0066e-05']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141209\\145307_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(14,27,2)
    # N_list = np.arange(16,25,2)

         ### Measurement set 
    # t_stamp =  '20141209_232400'    #older than time stamp to seperate different data sets
    # tau_list = ['1.0056e-05', '1.0058e-05','1.006e-05', '1.0066e-05']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141209\\145307_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(16,25,2)

    # t_stamp =  '20141211_153429'    #older than time stamp to seperate different data sets
    # tau_list = ['1.3586e-05','1.3588e-05','1.359e-05','1.3592e-05','1.3594e-05','1.3596e-05','1.3598e-05','1.36e-05','1.3602e-05','1.3604e-05'
    #             ,'1.3608e-05','1.361e-05','1.3612e-05','1.3614e-05','1.3616e-05','1.3618e-05']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141209\\145307_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(14,33,2)

    # t_stamp =  '20141211_130244'    #older than time stamp to seperate different data sets
    # tau_list = ['1.122e-05','1.1222e-05','1.1224e-05','1.1226e-05','1.1228e-05', '1.123e-05','1.1232e-05','1.1234e-05','1.1236e-05','1.1238e-05','1.124e-05','1.1242e-05', '1.1244e-05', '1.1246e-05', '1.1248e-05']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141209\\145307_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(16,33,2)

    # t_stamp =  '20141211_195841'    #older than time stamp to seperate different data sets
    # tau_list = ['1.1228e-05', '1.123e-05','1.1232e-05','1.1234e-05','1.1236e-05']
    # ssro_calib_folder = 'D:\\measuring\\data\\20141209\\145307_AdwinSSRO_SSROCalibration_111_1_sil18'
    # N_list = np.arange(10,17,2)



    ################################
    ### Measurements on Carbon 3 ###
    ################################

    ###
    ###Set 1
    ###

    # t_stamp =  '20150320_120000'    #older than time stamp to seperate different data sets
    # ssro_calib_folder = 'D:\\measuring\\data\\20150320\\022926_AdwinSSRO_SSROCalibration_111_1_sil18'
    # tau_list = np.arange(9.668e-6, 9.668e-6+16e-9, 2e-9)
    # N_list   = np.arange(6, 16, 2)
    
    # tau_list = np.arange(11.944e-6-10e-9, 11.944e-6+10e-9, 2e-9)
    # N_list   = np.arange(10,24, 2)

    # tau_list = np.arange(13.082e-6-10e-9, 13.082e-6+10e-9, 2e-9)
    # N_list   = np.arange(10, 24, 2)



    ###
    ### Set 2
    ### --> FAILED!!!

    # t_stamp =  '20150322_120000'
    # ssro_calib_folder = 'D:\\measuring\\data\\20150320\\163938_AdwinSSRO_SSROCalibration_111_1_sil18'
    
    # tau_list = np.arange(9.668e-6-4e-9, 9.668e-6+16e-9, 2e-9)
    # N_list   = np.arange(6, 18, 2)

    # tau0=10.806
    # tau_list = np.arange(tau0*1e-6-14e-9, tau0*1e-6+14e-9, 2e-9)
    # N_list   = np.arange(6, 22, 2)


    # tau0=15.362
    # tau_list = np.arange(tau0*1e-6-14e-9, tau0*1e-6+14e-9, 2e-9)
    # N_list   = np.arange(6, 22, 2)

    # tau0=16.5
    # tau_list = np.arange(tau0*1e-6-14e-9, tau0*1e-6+14e-9, 2e-9)
    # N_list   = np.arange(6, 22, 2)


    ###
    ### Set 3
    ###
    # t_stamp =  '20150327_120000'    #older than time stamp to seperate different data sets
    # ssro_calib_folder = 'D:\\measuring\\data\\20150327\\024343_AdwinSSRO_SSROCalibration_111_1_sil18'

    # tau_list = np.arange(9.668e-6-4e-9, 9.668e-6+16e-9, 2e-9)
    # N_list   = np.arange(6, 18, 2)

    # tau0=15.362
    # tau_list = np.arange(tau0*1e-6-14e-9, tau0*1e-6+14e-9, 2e-9)
    # N_list   = np.arange(12, 24, 2)

    # tau0=16.5
    # tau_list = np.arange(tau0*1e-6-14e-9, tau0*1e-6+14e-9, 2e-9)
    # N_list   = np.arange(12, 24, 2)


#######################################################
#### Evaluation                                    ####
#######################################################
    for ii, tau in enumerate(tau_list):

        folder_list_positive=[]
        folder_list_negative=[]


        # if tau in ['1.3586e-05','1.3588e-05','1.359e-05','1.3592e-05','1.3594e-05','1.3596e-05','1.3598e-05','1.36e-05','1.3602e-05','1.3604e-05'
        #         ,'1.3608e-05']:
        #     t_stamp =  '20141211_130244'  
        # else:
        #     t_stamp =  '20141211_161313'

        # if tau in ['1.361e-05','1.3612e-05','1.3614e-05']:
        #     # print 'yes'
        #     N_list = np.arange(14,39,2)
        # elif tau in ['1.3616e-05','1.3618e-05']:
        #     # print 'yes'
        #     N_list = np.arange(30,39,2)

        # if tau in ['1.1228e-05', '1.123e-05','1.1232e-05','1.1234e-05','1.1236e-05']:
        #     N_list = np.arange(10,33,2)
        #     t_stamp =  '20141211_195841' 
        # else:
        #     t_stamp =  '20141211_130244'
        #     N_list = np.arange(16,33,2)




        for k in N_list:
            # print tau
            # print k 
            
            folder_positive = toolbox.latest_data(contains = 'positive_tau'+ str(tau) + '_N'+ str(k), older_than=t_stamp)
            #print folder_positive
            folder_negative = toolbox.latest_data(contains = 'negative_tau'+ str(tau) + '_N'+ str(k) , older_than=t_stamp)
            #print folder_negative
            ### The folders below can be used to analyze the results for carbon 3 (faulty folder naming fo these measurements)
            # folder_negative = toolbox.latest_data(contains = '_tau'+ str(tau) + '_N'+ str(k)+'negative' , older_than=t_stamp)
            # folder_positive = toolbox.latest_data(contains = '_tau'+ str(tau) + '_N'+ str(k)+'positive', older_than=t_stamp)
            folder_list_positive = folder_list_positive+[folder_positive]
            folder_list_negative = folder_list_negative+[folder_negative]

        y_all = []
        y_all_u = []
        for i,folder in enumerate(folder_list_positive):
                                   
            # print folder_list_positive[i]     
            # print folder_list_negative[i]           
            
            a = mbi.MBIAnalysis(folder_list_positive[i])
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)
            
            b = mbi.MBIAnalysis(folder_list_negative[i])
            b.get_sweep_pts()
            b.get_readout_results(name='adwindata')
            b.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)

            a.p0 = 2*a.p0-1; a.u_p0 = 2*a.u_p0
            b.p0 = 2*b.p0-1; b.u_p0 = 2*b.u_p0

            a.p0 = (a.p0 - b.p0)/2
            a.u_p0 = ((a.u_p0**2 + b.u_p0**2)**0.5)/2

            a.u_p0 =   (1./(a.p0[0]**2 + a.p0[1]**2) * (a.p0[0]**2 * a.u_p0[0]**2 + a.p0[1]**2 * a.u_p0[1]**2))**0.5
            a.p0 = (a.p0[0]**2 + a.p0[1]**2)**0.5
             
            y_all.append(a.p0[0])
            y_all_u.append(a.u_p0[0])
       
        if ii == 0:
            fig = a.default_fig(figsize=(7.5,5))
            ax2 = a.default_ax(fig)
            start, end = ax2.get_xlim()
            ax2.axhspan(0,1,fill=False,ls='dotted')
            ax2.set_ylim(0.4,1)
            ax2.set_xlim(N_list[0],N_list[-1]+1)
            ax2.set_xlim(N_list[0]-2,N_list[-1]+1+2)
            colors = cm.rainbow(np.linspace(0, 1, len(tau_list)))

        # ax2.errorbar(N_list, y_all, y_all_u, 0*np.ones(len(y_all_u)), '.-', lw=1, color = colors[ii],
        #         label = 'tau = ')
        ax2.errorbar(N_list, y_all, y_all_u, y_all_u,'.-', lw=1, color = colors[ii],
                label = 'tau = ' +str(tau) + '_max = ' + str(np.around(y_all[np.argmax(y_all)],3)) + ' +/- ' + str(np.around(y_all_u[np.argmax(y_all)],3)) + ' at N = ' + str(N_list[np.argmax(y_all)])) #N = 16
        ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

                #N = 16
        # ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
      
        # print 'minimum = ' + str(amplitude[np.argmin(amplitude)]) + 'at N = ' + str(x[np.argmin(amplitude)])

        # x_av    = (y_all[1]+y_all[2])/2
        # x_av_u  = ((y_all_u[1]**2+y_all_u[2]**2)**0.5)/2
        # y_av    = (y_all[3]+y_all[4])/2
        # y_av_u  = ((y_all_u[3]**2+y_all_u[4]**2)**0.5)/2
        # z_av    = y_all[0]
        # z_av_u  = y_all_u[0]


        # fig3 = a.default_fig(figsize=(5,5))
        # ax3 = a.default_ax(fig3)
        # start, end = ax3.get_xlim()
        # ax3.set_xlim(x[0]-1,x[-1]+1)
        # ax3.errorbar(x, x_av, x_av_u,0*np.ones(len(amplitude)),  '.-r', lw=1,label = 'x') #N = 16
        # ax3.errorbar(x, y_av, x_av_u,0*np.ones(len(amplitude)), '.-b', lw=1,label = 'y') #N = 16
        # ax3.errorbar(x, z_av, x_av_u,0*np.ones(len(amplitude)), '.-g', lw=1,label = 'z') #N = 16
           
        # ax3.legend()




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


def Carbon_gate_optimization_routine(older_than=None, measurement_name = ['adwindata'], carbon = 1,ssro_calib_folder = None):
    ''' Function to analyze data for optimization of the number of pulses for a controlled C13 gate. 
        The timestamp given should be for the no_pulse measurement
    '''
    folder_list_positive= []
    folder_list_negative = []

    if carbon == 1:
        tau_list = [7.218e-6,4.994e-6,4.994e-6,4.996e-6,4.996e-6,
                               4.996e-6,4.998e-6,4.998e-6,7.214e-6]
        N_list  = [40,34,36,32,34,36,34,36,42]
        # tau_list = [4.998e-6,4.998e-6,5.000e-6]
        # N_list =[38,40,40]
    if carbon == 2:
        tau_list = [13.612e-6,13.612e-6,13.612e-6,13.614e-6,13.614e-6,13.614e-6,13.616e-6
                                ,13.616e-6,13.616e-6]
        N_list  = [26,28,30,30,32,34,32,34,36]

    if carbon == 5:

        tau_list = [6.536e-6,6.536e-6,6.536e-6,6.538e-6,6.538e-6, 6.538e-6,6.540e-6,6.540e-6,6.540e-6]
        N_list   = [30,32,34,30,32,34,30,32,34]

        tau_list =  [11.308e-6, 11.308e-6, 11.308e-6, 11.310e-6, 11.310e-6, 11.310e-6, 11.312e-6, 11.312e-6, 11.312e-6]
        N_list   = [44,46,48,46,48,50,46,48,50]

    if carbon == 3:
        tau_list = np.arange(9.668e-6-10e-9, 9.668e-6+10e-9, 2e-9)
        N_list=[]
        for i,t in enumerate(tau_list):
            N_list.extend(np.arange(12, 24, 2))


    x_tick_labels = []
    for i in range(len(tau_list)):
        if carbon == 3:
            for j in range(len(N_list)):
                x_tick_labels = x_tick_labels+ ['tau: '+str(tau_list[i]*1e6)+' us, N: '+str(N_list[j])]
        else:
            x_tick_labels = x_tick_labels+ ['tau: '+str(tau_list[i]*1e6)+' us, N: '+str(N_list[i])]
    # ssro_calib_folder = 'D:\\measuring\data\\20141203\\041410_AdwinSSRO_SSROCalibration_111_1_sil18'

    for ii in range(len(tau_list)):
        if carbon == 3:
            for jj,Nj in enumerate(np.arange(12, 24, 2)):
                folder_positive = toolbox.latest_data(contains = '_positive_tau'+str(tau_list[ii])+'_N'+str(Nj), older_than=older_than)
                folder_negative = toolbox.latest_data(contains = '_negative_tau'+str(tau_list[ii])+'_N'+str(Nj), older_than=older_than)
                
                folder_list_positive = folder_list_positive+[folder_positive]
                folder_list_negative = folder_list_negative+[folder_negative]
        else:
            folder_positive = toolbox.latest_data(contains = 'Gate_calibration_Sil18_C'+str(carbon)+'_ii_'+ str(ii) +'_positive', older_than=older_than)
            folder_negative = toolbox.latest_data(contains = 'Gate_calibration_Sil18_C'+str(carbon)+'_ii_'+ str(ii) +'_negative', older_than=older_than)
            folder_list_positive = folder_list_positive+[folder_positive]
            folder_list_negative = folder_list_negative+[folder_negative]
        # folder_positive = toolbox.latest_data(contains = '_ii_'+ str(ii) +'_positive', older_than=older_than)
        # folder_negative = toolbox.latest_data(contains = '_ii_'+ str(ii) +'_negative', older_than=older_than)
         

        

        y_all = []
        y_all_u = []
    print len(folder_list_positive)

    for i,folder in enumerate(folder_list_positive):

        a = mbi.MBIAnalysis(folder_list_positive[i])
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC()
        
        b = mbi.MBIAnalysis(folder_list_negative[i])
        b.get_sweep_pts()
        b.get_readout_results(name='adwindata')
        b.get_electron_ROC()

        a.p0 = 2*a.p0-1; a.u_p0 = 2*a.u_p0
        b.p0 = 2*b.p0-1; b.u_p0 = 2*b.u_p0

        a.p0 = (a.p0 - b.p0)/2
        a.u_p0 = ((a.u_p0**2 + b.u_p0**2)**0.5)/2

        a.u_p0 =   (1./(a.p0[0]**2 + a.p0[1]**2)*(
        a.p0[0]**2 * a.u_p0[0]**2 + 
        a.p0[1]**2 * a.u_p0[1]**2))**0.5
        a.p0 = (a.p0[0]**2 + a.p0[1]**2)**0.5


        y_all.append(a.p0[0])
        y_all_u.append(a.u_p0[0])

    # print len(np.arange(12, 24, 2))
    # print y_all

    if carbon == 3:
        fig = a.default_fig(figsize=(7.5,5))
        ax2 = a.default_ax(fig)
        ax2.axhspan(0,1,fill=False,ls='dotted')
        ax2.set_ylim(0.4,1)
        colors = cm.rainbow(np.linspace(0, 1, len(tau_list)))

        for i,t in enumerate(tau_list):
            x = np.arange(12, 24, 2)
            y = y_all[i*len(np.arange(12, 24, 2)):(i+1)*len(np.arange(12, 24, 2))]
            y_u = y_all_u[i*len(np.arange(12, 24, 2)):(i+1)*len(np.arange(12, 24, 2))]
            label = 'tau = ' +str(t) + '_min = ' + str(np.around(y[np.argmax(y)],3)) + ' +/- ' + str(np.around(y_u[np.argmax(y)],3)) + ' at N = ' + str(np.arange(12,24,2)[np.argmax(y)])
            ax2.errorbar(x, y, yerr=y_u,fmt='o',label=label, color = colors[i])

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        # start, end = ax2.get_xlim()

        # x_tick_labels = []
        # for i in range(len(tau_list)):
        #     x_tick_labels = x_tick_labels+ ['tau: '+str(tau_list[i]*1e6)+' us, N: '+str(N_list[i])+ '\n height: ' + str(np.around(y_all[i],3))]

        # ax2.bar(range(len(tau_list)), y_all, width=0.8,yerr = y_all_u)
        # ax2.set_xticks(np.arange(9) + 0.8/2)
        # ax2.set_xticklabels(x_tick_labels, rotation=90)
        # ax2.get_xaxis().set_tick_params(direction='in',pad = -120, labelcolor = 'white')
        # ax2.text(0.8,0.9,'max at '+ x_tick_labels[np.argmax(y_all)])

    else:
        fig = a.default_fig(figsize=(7.5,5))
        ax2 = a.default_ax(fig)
        start, end = ax2.get_xlim()
        ax2.axhspan(0,1,fill=False,ls='dotted')
        # ax2.set_ylim(0.4,1)

        x_tick_labels = []
        for i in range(len(tau_list)):
            x_tick_labels = x_tick_labels+ ['tau: '+str(tau_list[i]*1e6)+' us, N: '+str(N_list[i])+ '\n height: ' + str(np.around(y_all[i],3))]

        ax2.bar(range(len(tau_list)), y_all, width=0.8,yerr = y_all_u)
        ax2.set_xticks(np.arange(9) + 0.8/2)
        ax2.set_xticklabels(x_tick_labels, rotation=90)
        ax2.get_xaxis().set_tick_params(direction='in',pad = -120, labelcolor = 'white')
        ax2.text(0.8,0.9,'max at '+ x_tick_labels[np.argmax(y_all)])

    plt.savefig(os.path.join(folder_list_positive[0], 'analyzed_result.pdf'),
    format='pdf')