'''
This scripts analysis multiple DESR measurements to determine the long term magnetic field stability.
It is initially made to work for an NV center with a strongly couple C13 spin (it fits with a double gaussian).
by THT, adapted from JC
'''

import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import analysis.lib.QEC.nuclear_spin_characterisation as SC 
import matplotlib.cm as cm
import os
from matplotlib import pyplot as plt
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.tools import analysis_magnet_tools as amt
from analysis.lib.m2.ssro import sequence
from analysis.lib.fitting import dark_esr_auto_analysis; 
import matplotlib.mlab as mlab
from tempfile import TemporaryFile
reload(dark_esr_auto_analysis)
# reload(toolbox)
# reload(amt)
# reload(fit)

# Start of data 20160712 183728
# end of data

def fit_B_msmt_loop(older_than = None, newer_than = None, filename = 'magnet_Zpos_optimize_init', contains= 'magnet_'):
    ''' Just loads the data'''
    f0m = []; u_f0m = []; #f0p = [] ;u_f0p = []
    # Bx_field_measured = []
    # Bz_field_measured = []
    # f_centre_list = []; f_diff_list=[]
    timestamp_list = []
    # f_centre_error_list= []
    it_list = []
    # f_diff_error_list = []
    absolute_time_list = []
    temperature_list = []
    magnet_position_list = []

    #msm
    print 'start'
    older_than_SSRO = older_than
    print older_than, newer_than
    iteration = 0
    while toolbox.latest_data(contains=contains, older_than=older_than, newer_than=newer_than,raise_exc = False)!=False:
        print 'ITERATION '+str(iteration)
        ## Find the data folder
        
            ### This was needed for some of the meaurements where both a fine and a coarse measurement was done, THT 160725
        # timestamp,folder = toolbox.latest_data(contains='magnet_msm1_coarse', older_than=older_than, newer_than=newer_than,return_timestamp = True)
        # older_than = str(int(timestamp)-1)

        timestamp,folder = toolbox.latest_data(contains=contains, older_than=older_than, newer_than=newer_than,return_timestamp = True)
        print 'm folder '+folder
        
        ## Load the data ##
        ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration', older_than=older_than_SSRO)
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')

        print ssro_calib_folder
        a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

        print ' check point      '

        f0m_temp, u_f0m_temp = dark_esr_auto_analysis.analyze_dark_esr_double(do_plot=False, add_folder = folder)

        # #msp
        # folder = toolbox.latest_data(contains='msmt_msp_', older_than=older_than, newer_than=newer_than,)
        # print folder
        # a = sequence.SequenceAnalysis(folder)
        # a.get_sweep_pts()
        # a.get_readout_results('ssro')
        # a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

        # f0p_temp,u_f0p_temp = dark_esr_auto_analysis.analyze_dark_esr(None, 2.196*1e-3, add_folder = folder )
        
        # print f0p_temp >0 and f0m_temp >0
        # if f0p_temp >0 and f0m_temp >0:
            # Bz_measured, Bx_measured = amt.get_B_field(msm1_freq=f0m_temp*1e9, msp1_freq=f0p_temp*1e9, u_msm1_freq =u_f0m_temp ,u_msp1_freq=u_f0p_temp)
            # f_centre = (f0m_temp+f0p_temp)/2
            # f_centre_error = np.sqrt(u_f0m_temp**2+u_f0p_temp**2)/2
            # f_diff = (f_centre-ZFS*1e-9)*1e6
            # f_diff_error = f_centre_error*1e6
        f0m.append(f0m_temp)
        u_f0m.append(u_f0m_temp)
            # f0p.append(f0p_temp)
            # u_f0p.append(u_f0p_temp)
            # f_centre_list.append(f_centre)
            # f_centre_error_list.append(f_centre_error)
            # f_diff_list.append(f_diff)
            # f_diff_error_list.append(f_diff_error)
            # Bx_field_measured.append(Bx_measured)
            # Bz_field_measured.append(Bz_measured)
        timestamp_list.append(timestamp)
        it_list.append(iteration)


        ### Get the temperature data
        temperature = (a.g.attrs['temp']-100)/0.385
        temperature_list.append(temperature)

        magnet_position = (a.g.attrs['magnet_position'])
        magnet_position_list.append(magnet_position)

        absolute_time = int(timestamp[-2:]) + int(timestamp[-4:-2])*60 + int(timestamp[-6:-4])*60**2 + int(timestamp[6:8])*(60**2)*24
        absolute_time_list.append(absolute_time)

        older_than = str(int(timestamp)-1)
        print 'older than: ' + older_than
        iteration = iteration+1
    
    it_list = np.linspace(0,len(it_list)-1,len(it_list))
    outfile = TemporaryFile()

    f0m = list(reversed(f0m))
    u_f0m = list(reversed(u_f0m))
    temperature_list = list(reversed(temperature_list))
    timestamp_list = list(reversed(timestamp_list))
    absolute_time_list = list(reversed(absolute_time_list))
    magnet_position_list = list(reversed(magnet_position_list))
         
    
    np.savez(filename,f0m=f0m,
            u_f0m=u_f0m,
            # f0p=f0p,
            # u_f0p=u_f0p,
            # f_centre_list=f_centre_list,
            # f_centre_error_list=f_centre_error_list,
            # f_diff_list=f_diff_list,
            # f_diff_error_list=f_diff_error_list,
            # Bx_field_measured=Bx_field_measured,
            # Bz_field_measured=Bz_field_measured,
            timestamp_list=timestamp_list,
            absolute_time_list=absolute_time_list,
            it_list=it_list, 
            temperature_list = temperature_list,
            magnet_position_list = magnet_position_list)




def plot_meas_B_loop(filename = 'backlash_calib_data',filename2 = 'backlash_calib_data', optimization_target = 10, plot_out_of_range = True):

    d = np.load(filename+'.npz')
    f0m = d['f0m']; u_f0m = d['u_f0m']; #f0p = d['f0p'] ;u_f0p = d['u_f0p']
    # Bx_field_measured = d['Bx_field_measured']
    # Bz_field_measured = d['Bz_field_measured']
    # f_centre_list = d['f_centre_list']; f_diff_list=d['f_diff_list']
    timestamp_list = d['timestamp_list']
    # f_centre_error_list= d['f_centre_error_list']
    it_list = d['it_list']
    # f_diff_error_list = d['f_diff_error_list']
    abs_time = d['absolute_time_list']
    temperature_list = d['temperature_list']
    magnet_position_list = d['magnet_position_list']
    # print 'magnet position is   ' +str(magnet_position_list)

    # f0m_c = f0m - (temperature_list-22.24)*.001133


    d2 = np.load(filename2+'.npz')
    f0p = d2['f0m']; u_f0m = d['u_f0m'];
    dZFS= ((f0p+f0m)/2 - 2.877623)

    print 'dZFS is' + str(1e6*dZFS)


    ### Deleting all the failed fits 

    indices_failed_fits         = [i for i,j in enumerate(u_f0m*1e6) if j>5]
    f0m_filtered                = np.delete(f0m,indices_failed_fits)
    u_f0m_filtered              = np.delete(u_f0m,indices_failed_fits)
    timestamp_list_filtered     = np.delete(timestamp_list,indices_failed_fits)
    it_list_filtered            = np.delete(it_list,indices_failed_fits)
    abs_time_filtered           = np.delete(abs_time,indices_failed_fits)
    temperature_list_filtered   = np.delete(temperature_list,indices_failed_fits)
    magnet_position_list        = np.delete(magnet_position_list,indices_failed_fits)
   

    print 'Number of failed fits = ' + str(len(indices_failed_fits))


    mean_fms0         = np.mean(f0m_filtered)   
    stdev_fms0        = np.std(f0m_filtered)

    # stdev_temperature_list = np.std(temperature_list)

    ##folder to save the figures
    folder = toolbox.data_from_time(timestamp=timestamp_list[0])
    print 'The graphs are saved in:'
    print folder

    print 'Mean frequency = ' + str(mean_fms0)
    print 'Std = ' + str(stdev_fms0)


    # f0m_corrected = f0m_filtered - (temperature_list-22.24)*.001133
    relative_time   = (abs_time_filtered-abs_time[0])/60.
    relative_freq   = (f0m_filtered-1.746)*1e6
    relative_freq_u = u_f0m_filtered*1e6
    # # relative_freq_c = (f0m_corrected-1.746)*1e6

    # f0m_corrected = f0m_filtered - (temperature_list-22.24)*.001133
    # relative_freq_corrected   = (f0m_corrected-1.746)*1e6
    # relative_freq_unfiltered  =  (f0m-1.746)*1e6
    # iteration_index = list(reversed(it_list_filtered))

    # print 'iteration index   '+str(iteration_index)
    # print len(iteration_index)
    # print len(relative_freq_corrected)
    
    ### Figure freq vs time filtered
    fig = plt.figure(1,figsize=(18,5))
    ax = fig.add_subplot(111)
    ax.errorbar(relative_time,relative_freq,u_f0m_filtered*1e6)
    ax.set_xlabel('minutes from ' + timestamp_list[-1])
    ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
    # splt.ylim(48s,500)
    # plt.ylim(np.max(relative_freq)-400,np.max(relative_freq)+50)
    # plt.xlim(-300,-100)
    plt.title('Magnetic field vs time')
    plt.savefig(folder+'/freq_vs_time.png',format='png')




    # plot ZFS distribution

    mean_dZFS         = np.mean(dZFS)   
    stdev_dZFS        = np.std(dZFS)

    plt.figure(2)
    n, bins, patches = plt.hist((dZFS)*1e6,50,normed = 1)
    bincenters = 0.5*(bins[1:]+bins[:-1])
    # y = mlab.normpdf( bincenters, 1e6*(mean_fms0-1.746), 1e6*stdev_fms0)
    # plt.plot(bincenters, y, 'r--', linewidth=1)
    plt.xlabel('binned dZFS (kHz)')
    plt.title('Mean '+str(round(mean_dZFS*1e6,6))+' KHz, stdev '+str(round(stdev_dZFS*1e6,2))+' kHz')
    plt.savefig(folder+'/binned_dZFS.png',format='png')



   



    # plot magnetic field distribution 

    plt.figure(0)
    n, bins, patches = plt.hist((f0m_filtered-1.746)*1e6,50,normed = 1)
    bincenters = 0.5*(bins[1:]+bins[:-1])
    # y = mlab.normpdf( bincenters, 1e6*(mean_fms0-1.746), 1e6*stdev_fms0)
    # plt.plot(bincenters, y, 'r--', linewidth=1)
    plt.xlabel('binned relative center freq (kHz)')
    plt.title('Mean '+str(round(mean_fms0,6))+' GHz, stdev '+str(round(stdev_fms0*1e6,2))+' kHz')
    plt.savefig(folder+'/binned_freq.png',format='png')


    # Plot magnet position vs time
    plt.figure(3,figsize=(18,5))
    plt.plot(relative_time,magnet_position_list)
    plt.xlabel('minutes from ' + timestamp_list[-1])
    plt.ylabel('freq (kHz) (offset 1.746 GHz)')
    plt.title('Magnet position vs time')
    plt.savefig(folder+'/freq_vs_time.png',format='png')
    


    ## Figure temperature vs time filtered
    plt.figure(4,figsize=(18,5))
    plt.plot(relative_time,temperature_list_filtered)
    plt.xlabel('minutes from ' + timestamp_list[-1])
    plt.ylabel('Temperature')
    plt.title('Stability measurement M1 - temperature')
    plt.savefig(folder+'/temperature_vs_time.png',format='png')


    #### trial
    plt.figure(5,figsize=(25,10))
    plt.plot(relative_time,relative_freq, 'r')
    plt.plot(relative_time,700*temperature_list_filtered+(relative_freq[0]-700*temperature_list_filtered[0]),'b')
    plt.plot(relative_time,70000*magnet_position_list+(relative_freq[0]-70000*magnet_position_list[0]),'g')
    plt.grid(True)
    plt.show


    if plot_out_of_range == True:

        meas_num = range(0,len(relative_freq))
        out_of_range_flag = len(relative_freq)*[0]
        out_of_range_time = 0

        for i in range(len(relative_freq)-1):

            if abs(relative_freq[i]-666) > optimization_target:
                out_of_range_flag[i] = 1

            else: 
                out_of_range_flag[i] = 0


            out_of_range_time = out_of_range_time + (relative_time[i+1]-relative_time[i])*out_of_range_flag[i]

        out_of_range_perc = 100*out_of_range_time/(relative_time[len(relative_time)-1]-relative_time[0])

        # print relative_time

        print '.......                        ............                 ..............'
        print ' out of range time in % = ' + str(out_of_range_perc)

        print ''
        print '.............................................................................'



        b = 0
        s = 0
        delta_t = 0 
        delta_time = []

        # for i in range(0,len(out_of_range_flag)-1):

        #     if out_of_range_flag[i]-b > 0:
        #         b = 1
        #         delta_t = (relative_time[i]-relative_time[s])
        #         delta_time.append(delta_t)
        #         s = i 
        #     elif out_of_range_flag[i]-b < 0:
        #         b = 0



        for i in range(0,len(out_of_range_flag)-1):

            if out_of_range_flag[i]-b > 0:
                b = 1
                delta_t = (relative_time[i]-relative_time[s])
                delta_time.append(delta_t)
                 
            elif out_of_range_flag[i]-b < 0:
                b = 0
                s = i

        average_delta_time = np.average(delta_time)

        print 'average delta time is ' + str(average_delta_time)

        # print delta_time
            # plotting the in range flag vs time
        plt.figure(10)
        plt.plot(relative_time,out_of_range_flag)
        plt.xlabel('time')
        plt.ylabel('out of range flag')
        # plt.xlim(-300,-100)
        plt.ylim(0,1.1)
        plt.title('out of range flag vs time')
        plt.savefig(folder+'/out of range flag vs time.png',format='png')


        plt.figure(11)
        plt.plot(meas_num,out_of_range_flag)
        plt.xlabel('time')
        plt.ylabel('In range flag')
        # plt.xlim(0,200)
        plt.ylim(0,1.1)
        plt.title('out of  range flag vs measurement number')
        plt.savefig(folder+'/out of range flag vs measurement number.png',format='png')

        # plot magnetic field distribution 

        mean_delta_time   = np.mean(delta_time)   
        stdev_delta_time  = np.std(delta_time)

        print 'time between two consecutive optimizations' + str(delta_time)

        print 'mean time between two consecutive optimizations is ' + str(mean_delta_time)
        print 'standard deviation is ' + str(stdev_delta_time)

        plt.figure(12)
        n, bins, patches = plt.hist(delta_time,50,normed = 1, facecolor= 'green', alpha=0.75)
        #bincenters = 0.5*(bins[1:]+bins[:-1])
        # y = mlab.normpdf( bincenters, 1e6*(mean_fms0-1.746), 1e6*stdev_fms0)
        # plt.plot(bincenters, y, 'r--', linewidth=1)

        # y = mlab.normpdf( bins, mean_delta_time, stdev_delta_time)
        # l = plt.plot(bins, y, 'r--', linewidth=1)


        plt.xlabel('binned delta_time in minutes')
        plt.title('Mean '+str(mean_delta_time)+' minutes, stdev '+str(stdev_delta_time)+'minutes')
        plt.grid(True)
        # plt.show()
        
        plt.savefig(folder+'/binned_delta_time.png',format='png')








def plot_steps(filename = 'magnet_z_loop_gr100_steps_tempcorr', optimization_target = 10):

    d = np.load(filename+'.npz')
    f0m = d['f0m']; u_f0m = d['u_f0m']; #f0p = d['f0p'] ;u_f0p = d['u_f0p']
    # Bx_field_measured = d['Bx_field_measured']
    # Bz_field_measured = d['Bz_field_measured']
    # f_centre_list = d['f_centre_list']; f_diff_list=d['f_diff_list']
    timestamp_list = d['timestamp_list']
    # f_centre_error_list= d['f_centre_error_list']
    it_list = d['it_list']
    # f_diff_error_list = d['f_diff_error_list']
    abs_time = d['absolute_time_list']
    # temperature_list = d['temperature_list']
    magnet_position_list = d['magnet_position_list']
    # print 'magnet position is   ' +str(magnet_position_list)

    # f0m_c = f0m - (temperature_list-22.24)*.001133


    ### Deleting all the failed fits 

    indices_failed_fits         = [i for i,j in enumerate(u_f0m*1e6) if j>5]
    f0m_filtered                = np.delete(f0m,indices_failed_fits)
    u_f0m_filtered              = np.delete(u_f0m,indices_failed_fits)
    timestamp_list_filtered     = np.delete(timestamp_list,indices_failed_fits)
    it_list_filtered            = np.delete(it_list,indices_failed_fits)
    abs_time_filtered           = np.delete(abs_time,indices_failed_fits)
    # temperature_list            = np.delete(temperature_list,indices_failed_fits)
    magnet_position_list        = np.delete(magnet_position_list,indices_failed_fits)



    relative_time   = (abs_time_filtered-abs_time[0])/60.
    relative_freq   = (f0m_filtered-1.746)*1e6
    relative_freq_u = u_f0m_filtered*1e6
   

    print 'Number of failed fits = ' + str(len(indices_failed_fits))


    mean_fms0         = np.mean(f0m_filtered)   
    stdev_fms0        = np.std(f0m_filtered)

    # stdev_temperature_list = np.std(temperature_list)

    ##folder to save the figures
    folder = toolbox.data_from_time(timestamp=timestamp_list[0])
    print 'The graphs are saved in:'
    print folder



         ### Figure freq vs time filtered
    fig = plt.figure(1,figsize=(18,5))
    ax = fig.add_subplot(111)
    ax.errorbar(relative_time,relative_freq,u_f0m_filtered*1e6)
    ax.set_xlabel('minutes from ' + timestamp_list[-1])
    ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
    # splt.ylim(48s,500)
    # plt.ylim(np.max(relative_freq)-400,np.max(relative_freq)+50)
    # plt.xlim(-300,-100)
    plt.title('Magnetic field vs time')
    plt.savefig(folder+'/freq_vs_time.png',format='png')

    delta_time_steps =[]

    for i in range(0,len(relative_time)-1):

                delta_t_steps = (relative_time[i+1]-relative_time[i])
                delta_time_steps.append(delta_t_steps)
                 
            
    

    average_delta_time_steps = np.average(delta_time_steps)

    print 'average delta time is ' + str(average_delta_time_steps)


    mean_delta_time   = np.mean(delta_time_steps)   
    stdev_delta_time  = np.std(delta_time_steps)

    print 'time between two consecutive optimizations' + str(delta_time_steps)

    print 'mean time between two consecutive optimizations is ' + str(mean_delta_time)
    print 'standard deviation is ' + str(stdev_delta_time)




    plt.figure(10)
    n, bins, patches = plt.hist(delta_time_steps,50,normed = 1, facecolor= 'green', alpha=0.75)
    plt.xlabel('binned delta_time in minutes')
    plt.title('Mean '+str(mean_delta_time)+' minutes, stdev '+str(stdev_delta_time)+'minutes')
    plt.grid(True)
    plt.show()
    
    plt.savefig(folder+'/binned_delta_time.png',format='png')










     # .............................................................................................. #







    # fig = plt.figure(2,figsize=(18,5))
    # ax = fig.add_subplot(111)
    # ax.errorbar(relative_time,relative_freq_c,u_f0m_filtered*1e6)
    # ax.set_xlabel('minutes from ' + timestamp_list[-1])
    # ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
    # # splt.ylim(48s,500)
    # plt.ylim(np.max(relative_freq_c)-200,np.max(relative_freq_c)+60)
    # # plt.xlim(0,220)
    # plt.title('Temperature corrected Magnetic field vs time')
    # plt.savefig(folder+'/Temperature corrected Magnetic field vs time.png',format='png')


    # ### Figure freq vs time filtered
    # fig = plt.figure(0,figsize=(18,5))
    # ax = fig.add_subplot(111)
    # ax.errorbar(relative_time,relative_freq,u_f0m_filtered*1e6)
    # ax.set_xlabel('minutes from ' + timestamp_list[-1])
    # ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
    # # splt.ylim(440s,500)
    # # plt.ylim(np.min(relative_freq)-80,np.min(relative_freq)+200)
    # # plt.xlim(0,220)
    # plt.title('Magnetic field vs time')
    # plt.savefig(folder+'/freq_vs_time.png',format='png')





  

    #############################################

# ## plot magnetic field vs measurement number

    # # print 'iteration_index =' +str(iteration_index)

    # fig = plt.figure(3,figsize=(18,5))
    # ax = fig.add_subplot(111)
    # ax.errorbar(iteration_index,relative_freq_unfiltered,u_f0m_filtered*1e6*10)
    # ax.set_xlabel('measurement number')
    # ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
    # # splt.ylim(440s,500)
    # # plt.ylim(300,450)
    # # plt.xlim(0,200)
    # plt.title('temperature corrected magnetic field vs measurement number')
    # plt.savefig(folder+'/temperature corrected magnetic field vs measurement number.png',format='png')
    




    ### Figure temperature vs time filtered
    # fig = plt.figure(4,figsize=(18,5))
    # ax = fig.add_subplot(111)
    # ax.errorbar(relative_time,temperature_list,[0]*len(temperature_list))
    # ax.set_xlabel('minutes from ' + timestamp_list[-1])
    # ax.set_ylabel('Temperature')
    # plt.title('Stability measurement M1 - temperature')
    # plt.savefig(folder+'/temperature_vs_time.png',format='png')


    
    # plt.figure(5)
    # n, bins, patches = plt.hist(temperature_list,normed = 1)
    # bincenters = 0.5*(bins[1:]+bins[:-1])
    # # y = mlab.normpdf( bincenters, 1e6*(mean_fms0-1.746), 1e6*stdev_fms0)
    # # plt.plot(bincenters, y, 'r--', linewidth=1)
    # plt.xlabel('binned relative center freq (kHz)')
    # plt.title('stdev '+str(round(stdev_temperature_list*1e3))+' mK')
    # plt.savefig(folder+'/binned_temperature.png',format='png')



    ## plot magnetic field vs position

    # fig = plt.figure(6,figsize=(7,7))
    # ax = fig.add_subplot(111)
    # ax.plot(magnet_position_list,(f0m_filtered-1.746)*1e6,'ro')
    # plt.title('magnetic field vs position')
    # plt.savefig(folder+'/magnetic field vs position.png',format='png')



     ## plot temperature corrected magnetic field vs position

    # fig = plt.figure(7,figsize=(7,7))
    # ax = fig.add_subplot(111)
    # f0m_corrected = f0m_filtered - (temperature_list-22.24)*.001133
    # ax.plot(magnet_position_list,(f0m_corrected-1.746)*1e6,'bo')
    # plt.title('temperature corrected magnetic field vs position')
    # plt.savefig(folder+'/temperature corrected magnetic field vs position.png',format='png')







    # Plot and fit temperature corrected magnetic field vs position
  
    # fig = plt.figure(10,figsize=(7,7))
    # ax = fig.add_subplot(111)
    # ax.plot(magnet_position_list,(f0m_corrected-1.746)*1e6,'ro')

        # fitting correlation to a line
    # A_guess = 200       # kHz
    # B_guess = 100         # kHz/K

    # A = fit.Parameter(A_guess, 'A')
    # B = fit.Parameter(B_guess, 'B')

    # def fitfunc(x):
    #             return A() + B()*x

    # print 'running fit'
    # fit_result = fit.fit1d(magnet_position_list, (f0m_corrected-1.746)*1e6, None, p0 = [A, B],
    #                 fitfunc = fitfunc, do_print=True, ret=True, fixed=[])

    # A0 = fit_result['params_dict']['A']
    # B0 = fit_result['params_dict']['B']

    #     # plotting the fitted function
    # plot.plot_fit1d(fit_result, np.linspace(min(magnet_position_list), max(magnet_position_list), 1000), ax=ax, plot_data=True)
    # ax.set_xlabel('magnet position')
    # ax.set_ylabel('ESR Frequency')
    # plt.title('fitting of temperature corrected magnetic field vs position')
    # plt.savefig(folder+'/fitting of temperature corrected magnetic field vs position.png',format='png')






    # fig = plt.figure(12,figsize=(18,5))
    # ax = fig.add_subplot(111)
    # ax.plot(point_number,f_average_b,'bo-')
    # ax.set_xlabel('point number')
    # ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
    # # splt.ylim(440s,500)
    # # plt.ylim(300,450)
    # # plt.xlim(0,200)
    # plt.title('averaged magnetic field vs point number - bottom data')
    # plt.savefig(folder+'/averaged magnetic field vs point number - bottom data.png',format='png')







    # print f_bl
    # print f0m

    # print '# Back lash calibration with temperature correction'
    # # Back lash calibration with temperature correction
    # f_bl_c = [] 
    # j = 6
    # # f0m_c = f0m - (temperature_list-22.24)*.001133
    # f0m_c_reversed = list(reversed(f0m_c))
    # Back_lash_c = []

    # #print f0m

    # for i in range (num_turning_points):
    #     f_bl_c.append(f0m_c_reversed[j-1])
    #     f_bl_c.append(f0m_c_reversed[j])
    #     f_bl_c.append(f0m_c_reversed[j+1])
    #     delta_f = (f0m_c_reversed[j]-f0m_c_reversed[j+1])/(f0m_c_reversed[j]-f0m_c_reversed[j-1])
    #     Back_lash_c.append((1-delta_f)*step)
    #     j = j+12
    

    # Average_Back_lash_c =  np.average(Back_lash_c)     
    # Back_lash_c_st_dev = np.std(Back_lash_c)

    # print '............................................'

    # print '............................................'
    # print 'Average_Back_lash_corrected is    '+ str(Average_Back_lash_c) 

    # print '............................................'
    # print 'Standard_deviation_of_Back_lash_corrected is    '+ str(Back_lash_c_st_dev)
    # print '............................................'
    # print '............................................'
    # print '............................................'
    # print 'Back_lash_corrected is    '+ str(Back_lash_c) 
    # print '............................................'
    

    #print f_bl
    # print f0m


     ## plot averaged magnetic field 

    


    # print 'f0m_filtered is   '+ str(f0m_filtered)
    # print 'f0m_corrected is  ' + str(f0m_corrected)

        # fitting correlation to a line
    # A_guess = -22000        # kHz
    # B_guess = 1000          # kHz/K

    # A = fit.Parameter(A_guess, 'A')
    # B = fit.Parameter(B_guess, 'B')

    # def fitfunc(x):
    #             return A() + B()*x

    # print 'running fit'
    # fit_result = fit.fit1d(magnet_position_list, (f0m_filtered-1.746)*1e6, None, p0 = [A, B],
    #                 fitfunc = fitfunc, do_print=True, ret=True, fixed=[])

    # A0 = fit_result['params_dict']['A']
    # B0 = fit_result['params_dict']['B']

    #     # plotting the fitted function
    # plot.plot_fit1d(fit_result, np.linspace(min(magnet_position_list), max(magnet_position_list), 1000), ax=ax, plot_data=True)
    # ax.set_xlabel('Position')
    # ax.set_ylabel('ESR Frequency')
    # plt.title('Correlation_fit: Temperature and Frequency')

    # plt.savefig(folder+'/Correlation_fit.png',format='png')

    
    

    ## Plot and fit correlation function
  
    # fig = plt.figure(6,figsize=(7,7))
    # ax = fig.add_subplot(111)
    # ax.plot(magnet_positions_list,(f0m_filtered-1.746)*1e6,'ro')

    #     # fitting correlation to a line
    # A_guess = -22000        # kHz
    # B_guess = 1000          # kHz/K

    # A = fit.Parameter(A_guess, 'A')
    # B = fit.Parameter(B_guess, 'B')

    # def fitfunc(x):
    #             return A() + B()*x

    # print 'running fit'
    # fit_result = fit.fit1d(magnet_positions_list, (f0m_filtered-1.746)*1e6, None, p0 = [A, B],
    #                 fitfunc = fitfunc, do_print=True, ret=True, fixed=[])

    # A0 = fit_result['params_dict']['A']
    # B0 = fit_result['params_dict']['B']

    #     # plotting the fitted function
    # plot.plot_fit1d(fit_result, np.linspace(min(magnet_positions_list), max(magnet_positions_list), 1000), ax=ax, plot_data=True)
    # ax.set_xlabel('Temperature')
    # ax.set_ylabel('ESR Frequency')
    # plt.title('Correlation_fit: Temperature and Frequency')

    # plt.savefig(folder+'/Correlation_fit.png',format='png')


    # ### Comparison of field and temperature
    # fig = plt.figure(20,figsize=(18,5))
    # ax = fig.add_subplot(111)
    # ax.errorbar(relative_time,relative_freq,relative_freq_u)
    # ax.set_xlabel('minutes from ' + timestamp_list[-1])
    # ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
    # # splt.ylim(440s,500)
    # plt.ylim(300,450)

    ### Figure vs time filtered
    # ax.plot(relative_time,B0*temperature_list+A0,'r-')
    # # plt.xlim(500,700)
    # plt.title('Stability measurement M1 - temperature')
    # plt.savefig(folder+'/temperature_vs_time.png',format='png')


## A+B*temp(t-t0) to fre(t)

## freq(t) = B*(temp(t-t0)-temp0) + C*(magnet_pos(t)-magnet_pos0) + A
## freq(t) = B*temp(t-t0) + C*magnet_pos(t) + A


    ### Some other plots
    # ### Figure vs time filtered
    # ax2 = ax.twinx()
    # ax2.plot(relative_time,temperature_list+0.02,'r-')
    # ax2.set_xlabel('minutes from ' + timestamp_list[-1])
    # ax2.set_ylabel('Temperature')
    # # plt.xlim(500,700)
    # plt.title('Stability measurement M1 - temperature')
    # plt.savefig(folder+'/temperature_vs_time.png',format='png')

        ### Figure vs time filtered
    # fig = plt.figure(3,figsize=(18,5))
    # ax = fig.add_subplot(111)
    # ax.errorbar(relative_time_averaged,relative_freq_averaged,relative_freq_averaged_u)
    # ax.set_xlabel('minutes from ' + timestamp_list[-1])
    # ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
    # # splt.ylim(440s,500)
    # # plt.ylim(300,450)
    # # plt.xlim(0,200)
    # plt.title('Stability meaurment M1 averaged')
    # plt.savefig(folder+'/freq_vs_time.png',format='png')


  
  # # ## plot temperature corrected magnetic field vs time

  #   fig = plt.figure(2,figsize=(18,5))
  #   ax = fig.add_subplot(111)
  #   ax.errorbar(relative_time,relative_freq_corrected,u_f0m_filtered*1e6)
  #   ax.set_xlabel('minutes from ' + timestamp_list[-1])
  #   ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
  #   # splt.ylim(440s,500)
  #   plt.ylim(2000,2.3e3)
  #   plt.xlim(0,220)
  #   plt.title('Temperature corrected Magnetic field vs time')
  #   plt.savefig(folder+'/temp_corrected_freq_vs_time.png',format='png')



    