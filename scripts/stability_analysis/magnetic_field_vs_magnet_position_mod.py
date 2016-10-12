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

def fit_B_msmt_loop(older_than = None, newer_than = None, filename = 'backlash_calib_data'):
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
    while toolbox.latest_data(contains='magnet_scanner_calib_fine', older_than=older_than, newer_than=newer_than,raise_exc = False)!=False:
        print 'ITERATION '+str(iteration)
        ## Find the data folder
        
            ### This was needed for some of the meaurements where both a fine and a coarse measurement was done, THT 160725
        # timestamp,folder = toolbox.latest_data(contains='magnet_msm1_coarse', older_than=older_than, newer_than=newer_than,return_timestamp = True)
        # older_than = str(int(timestamp)-1)

        timestamp,folder = toolbox.latest_data(contains='magnet_scanner_calib_fine', older_than=older_than, newer_than=newer_than,return_timestamp = True)
        print 'm folder '+folder
        
        ## Load the data ##
        ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration', older_than=older_than_SSRO)
        a = sequence.SequenceAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results('ssro')
        a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

        print ' check point      '

        f0m_temp, u_f0m_temp = dark_esr_auto_analysis.analyze_dark_esr_double(do_plot=True, add_folder = folder)

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

def plot_meas_B_loop(filename = 'backlash_calib_data'):

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

    f0m_c = f0m - (temperature_list-22.24)*.001133


    ### Deleting all the failed fits 

    indices_failed_fits         = [i for i,j in enumerate(u_f0m*1e6) if j>5]
    f0m_filtered                = np.delete(f0m,indices_failed_fits)
    u_f0m_filtered              = np.delete(u_f0m,indices_failed_fits)
    timestamp_list_filtered     = np.delete(timestamp_list,indices_failed_fits)
    it_list_filtered            = np.delete(it_list,indices_failed_fits)
    abs_time_filtered           = np.delete(abs_time,indices_failed_fits)
    temperature_list            = np.delete(temperature_list,indices_failed_fits)
    magnet_position_list        = np.delete(magnet_position_list,indices_failed_fits)
   

    print 'Number of failed fits = ' + str(len(indices_failed_fits))

    # Bz_diff_list = [j-304.21 for j in Bz_field_measured]
    # Bz_error = [1/(4.*ZFS*g_factor*Bz_field_measured[j])*
    #         (f0p[j]**2*u_f0p[j]**2+f0m[j]**2*u_f0m[j]**2)**(1/2.) for j in range(len(Bz_field_measured))]
    # Bx_error = [1/Bx_field_measured[j]*(f0m[j]**2*u_f0m[j]**2/g_factor**2
    #         +(ZFS-g_factor*Bz_field_measured[j])**2*Bz_error[j]**2) for j in range(len(Bx_field_measured))]

    # print len(Bz_field_measured)
    # total_B = [(Bx_field_measured[j]**2+Bz_field_measured[j]**2)**(1/2.) for j in range(len(Bz_field_measured))]
    # total_B_error = [1/total_B[j]*(Bx_field_measured[j]**2*Bx_error[j]**2+Bz_field_measured[j]**2*Bz_error[j]**2)**(1/2.) for j in range(len(Bz_field_measured))]

    mean_fms0         = np.mean(f0m_filtered)   
    stdev_fms0        = np.std(f0m_filtered)

    stdev_temperature_list = np.std(temperature_list)

    ##folder to save the figures
    folder = toolbox.data_from_time(timestamp=timestamp_list[0])
    print 'The graphs are saved in:'
    print folder

    print 'Mean frequency = ' + str(mean_fms0)
    print 'Std = ' + str(stdev_fms0)


    f0m_corrected = f0m_filtered - (temperature_list-22.24)*.001133
    relative_time   = (abs_time_filtered-abs_time[-1])/60.
    relative_freq   = (f0m_filtered-1.746)*1e6
    relative_freq_u = u_f0m_filtered*1e6
    relative_freq_c = (f0m_corrected-1.746)*1e6


    relative_freq_averaged          = []
    relative_freq_averaged_u        = []
    relative_time_averaged          = []
    nr_to_average_over = 3
    for i in range(len(relative_freq)/nr_to_average_over):
        relative_freq_averaged.append( (relative_freq[3*i] + relative_freq[3*i+1] + relative_freq[3*i+2])/3.)
        relative_freq_averaged_u.append( (relative_freq_u[3*i]**2 + relative_freq_u[3*i+1]**2 + relative_freq_u[3*i+2]**2)**0.5/3.)
        relative_time_averaged.append(relative_time[3*i+1])



    f0m_corrected = f0m_filtered - (temperature_list-22.24)*.001133
    relative_freq_corrected   = (f0m_corrected-1.746)*1e6
    relative_freq_unfiltered  =  (f0m-1.746)*1e6
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
    plt.ylim(np.max(relative_freq)-200,np.max(relative_freq)+60)
    # plt.xlim(0,220)
    plt.title('Magnetic field vs time')
    plt.savefig(folder+'/freq_vs_time.png',format='png')



    fig = plt.figure(2,figsize=(18,5))
    ax = fig.add_subplot(111)
    ax.errorbar(relative_time,relative_freq_c,u_f0m_filtered*1e6)
    ax.set_xlabel('minutes from ' + timestamp_list[-1])
    ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
    # splt.ylim(48s,500)
    plt.ylim(np.max(relative_freq_c)-200,np.max(relative_freq_c)+60)
    # plt.xlim(0,220)
    plt.title('Temperature corrected Magnetic field vs time')
    plt.savefig(folder+'/Temperature corrected Magnetic field vs time.png',format='png')


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
    fig = plt.figure(4,figsize=(18,5))
    ax = fig.add_subplot(111)
    ax.errorbar(relative_time,temperature_list,[0]*len(temperature_list))
    ax.set_xlabel('minutes from ' + timestamp_list[-1])
    ax.set_ylabel('Temperature')
    plt.title('Stability measurement M1 - temperature')
    plt.savefig(folder+'/temperature_vs_time.png',format='png')


    
    # plt.figure(5)
    # n, bins, patches = plt.hist(temperature_list,normed = 1)
    # bincenters = 0.5*(bins[1:]+bins[:-1])
    # # y = mlab.normpdf( bincenters, 1e6*(mean_fms0-1.746), 1e6*stdev_fms0)
    # # plt.plot(bincenters, y, 'r--', linewidth=1)
    # plt.xlabel('binned relative center freq (kHz)')
    # plt.title('stdev '+str(round(stdev_temperature_list*1e3))+' mK')
    # plt.savefig(folder+'/binned_temperature.png',format='png')



    ## plot magnetic field vs position

    fig = plt.figure(6,figsize=(7,7))
    ax = fig.add_subplot(111)
    ax.plot(magnet_position_list,(f0m_filtered-1.746)*1e6,'ro')
    plt.title('magnetic field vs position')
    plt.savefig(folder+'/magnetic field vs position.png',format='png')



     ## plot temperature corrected magnetic field vs position

    fig = plt.figure(7,figsize=(7,7))
    ax = fig.add_subplot(111)
    f0m_corrected = f0m_filtered - (temperature_list-22.24)*.001133
    ax.plot(magnet_position_list,(f0m_corrected-1.746)*1e6,'bo')
    plt.title('temperature corrected magnetic field vs position')
    plt.savefig(folder+'/temperature corrected magnetic field vs position.png',format='png')







    # Plot and fit temperature corrected magnetic field vs position
  
    # fig = plt.figure(10,figsize=(7,7))
    # ax = fig.add_subplot(111)
    # ax.plot(magnet_position_list,(f0m_corrected-1.746)*1e6,'ro')

        # fitting correlation to a line
    A_guess = 200       # kHz
    B_guess = 100         # kHz/K

    A = fit.Parameter(A_guess, 'A')
    B = fit.Parameter(B_guess, 'B')

    def fitfunc(x):
                return A() + B()*x

    print 'running fit'
    fit_result = fit.fit1d(magnet_position_list, (f0m_corrected-1.746)*1e6, None, p0 = [A, B],
                    fitfunc = fitfunc, do_print=True, ret=True, fixed=[])

    A0 = fit_result['params_dict']['A']
    B0 = fit_result['params_dict']['B']

        # plotting the fitted function
    plot.plot_fit1d(fit_result, np.linspace(min(magnet_position_list), max(magnet_position_list), 1000), ax=ax, plot_data=True)
    ax.set_xlabel('magnet position')
    ax.set_ylabel('ESR Frequency')
    plt.title('fitting of temperature corrected magnetic field vs position')
    plt.savefig(folder+'/fitting of temperature corrected magnetic field vs position.png',format='png')



def backlash_calib(step=0.2, num_turning_points=10,starting_point=6,filename = 'backlash_calib_data'):


    d = np.load(filename+'.npz')
    f0m = d['f0m']; u_f0m = d['u_f0m']; #f0p = d['f0p'] ;u_f0p = d['u_f0p']
    timestamp_list = d['timestamp_list']
    it_list = d['it_list']
    abs_time = d['absolute_time_list']
    temperature_list = d['temperature_list']
    magnet_position_list = d['magnet_position_list']
    f0m_c = f0m - (temperature_list-22.24)*.001133


    folder = toolbox.data_from_time(timestamp=timestamp_list[0])
    print 'The graphs are saved in:'
    print folder

    # Back lash calibration
    j = starting_point
    Back_lash = []
    f_average = 10*[0]
    f_average_b = 10*[0]
    indices_failed_fits= [i for i,j in enumerate(u_f0m*1e6) if j>5]
    average_over = 10*[0]
    average_over_b = 10*[0]

    ###################################################################################################

    for i in range (num_turning_points):

        delta_f = (f0m[j]-f0m[j+1])/(f0m[j]-f0m[j-1])
        Back_lash.append((1-delta_f)*step)

        # j = j+12
        j = j+5
         
    Average_Back_lash =  np.average(Back_lash)   
    Back_lash_st_dev  =  np.std(Back_lash)

    #averaging the data points we have
    
    # j2= starting_point
    # for i in range (num_turning_points/2):
    #     for k in range (11):
    #         #think of adding condition over failed fits
    #         f_average[k] = f_average[k] + ((f0m[k+j2-5])/num_turning_points)
    #         f_average_b[k] = f_average_b[k] + ((f0m[k+j2+12-5])/num_turning_points)

    #     j2= j2+24

    #####################################################################################################
    
    
    j2= starting_point

    print 'indices of failed fits   '+ str(indices_failed_fits)

    for i in range (num_turning_points/2):
        for k in range (10):
            #think of adding condition over failed fits
            if (k+j2-5) not in indices_failed_fits:
                average_over[k]=average_over[k] + 1    
                f_average[k] = f_average[k] + ((f0m[k+j2-5]))

            #commented because not necessary for the current measurement
            # if (k+j2+12-5) not in indices_failed_fits:
            #     average_over_b[k] = average_over_b [k]+1
            #     f_average_b[k] = f_average_b[k] + ((f0m[k+j2+12-5]))


        # j2= j2+24
        j2 = j2+10

    for k in range (10):   
        f_average[k] = f_average[k]/average_over[k]
        # f_average_b[k] = f_average_b[k]/average_over_b[k]


    print 'step is   ' +str(step)
    print '............................................'
    print 'Average_Back_lash is    '+ str(Average_Back_lash) 
    print '............................................'
    print 'Standard_deviation_of_Back_lash is    '+ str(Back_lash_st_dev)
    print '............................................'
    print '............................................'

    print 'Back_lash is    '+ str(Back_lash) 
    print '............................................'

    print '............................................'
    print 'f_average is    '+ str(f_average) 
    print '............................................'


    point_number = list(range(1,11))
    print'..................trial...........'
    print point_number
    print f_average
    print'..................trial...........'
    fig = plt.figure(10,figsize=(18,5))
    ax = fig.add_subplot(111)
    ax.plot(point_number,f_average,'bo-')
    ax.set_xlabel('point number')
    ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
    # splt.ylim(440s,500)
    # plt.ylim(300,450)
    # plt.xlim(0,200)
    plt.title('averaged magnetic field vs point number')
    plt.savefig(folder+'/averaged magnetic field vs point number.png',format='png')




     #####################################################################################################

     # for temperature corrected
def backlash_calib_temp_c(step=0.2, num_turning_points=10,starting_point=6,filename = 'backlash_calib_data'):


    d = np.load(filename+'.npz')
    f0m = d['f0m']; u_f0m = d['u_f0m']; #f0p = d['f0p'] ;u_f0p = d['u_f0p']
    timestamp_list = d['timestamp_list']
    it_list = d['it_list']
    abs_time = d['absolute_time_list']
    temperature_list = d['temperature_list']
    magnet_position_list = d['magnet_position_list']
    f0m_c = f0m - (temperature_list-22.24)*.001133


    print type(temperature_list)

    folder = toolbox.data_from_time(timestamp=timestamp_list[0])
    print 'The graphs are saved in:'
    print folder

    # Back lash calibration
    j = starting_point
    Back_lash = []
    f_average = 10*[0]
    f_average_b = 10*[0]
    indices_failed_fits= [i for i,j in enumerate(u_f0m*1e6) if j>5]
    average_over = 10*[0]
    average_over_b = 10*[0]

    j = starting_point
    Back_lash = []
    f_average = 10*[0]
    f_average_b = 10*[0]
    indices_failed_fits= [i for i,j in enumerate(u_f0m*1e6) if j>5]
    average_over = 10*[0]
    average_over_b = 10*[0]
    j2= starting_point
    f0m_corr = f0m - (temperature_list-22.24)*.001133
    print 'indices of failed fits   '+ str(indices_failed_fits)

    for i in range (num_turning_points/2):
        for k in range (10):
            #think of adding condition over failed fits
            if (k+j2-5) not in indices_failed_fits:
                average_over[k]=average_over[k] + 1    
                f_average[k] = f_average[k] + ((f0m_corr[k+j2-5]))

            #commented because not necessary for the current measurement
            # if (k+j2+12-5) not in indices_failed_fits:
            #     average_over_b[k] = average_over_b [k]+1
            #     f_average_b[k] = f_average_b[k] + ((f0m[k+j2+12-5]))


        # j2= j2+24
        j2 = j2+10

    for k in range (10):   
        f_average[k] = f_average[k]/average_over[k]
        # f_average_b[k] = f_average_b[k]/average_over_b[k]


    print 'step is   ' +str(step)
    # print '............................................'
    # print 'Average_Back_lash is    '+ str(Average_Back_lash) 
    # print '............................................'
    # print 'Standard_deviation_of_Back_lash is    '+ str(Back_lash_st_dev)
    # print '............................................'
    # print '............................................'

  


    point_number = list(range(1,11))
    print'..................trial...........'
    print point_number
    print f_average
    print'..................trial...........'
    fig = plt.figure(15,figsize=(18,5))
    ax = fig.add_subplot(111)
    ax.plot(point_number,f_average,'bo-')
    ax.set_xlabel('point number')
    ax.set_ylabel('freq (kHz) (offset 1.746 GHz)')
    # splt.ylim(440s,500)
    # plt.ylim(300,450)
    # plt.xlim(0,200)
    plt.title('averaged magnetic field corrected vs point number')
    plt.savefig(folder+'/averaged magnetic field corrected vs point number.png',format='png')


    # Plot and fit correlation function
  
    fig = plt.figure(16,figsize=(7,7))
    ax = fig.add_subplot(111)
    # ax.plot(point_number[6:11],f_average[6:11],'ro')

        # fitting correlation to a line
    A_guess = 1.7        # kHz
    B_guess = -160e-6          # kHz/K

    A = fit.Parameter(A_guess, 'A')
    B = fit.Parameter(B_guess, 'B')

    def fitfunc(x):
                return A() + B()*x

    
    x1 = 0.2*np.array(point_number[1:6])
    y1 = (np.array(f_average[1:6])-1.76)*1e6
    print 'x1   ' +str(x1)
    print 'y1   ' +str(y1)
    print fitfunc(x1)

    print 'running fit'
    fit_result = fit.fit1d(x1, y1, None, p0 = [A, B],
                    fitfunc = fitfunc, do_print=True, ret=True, fixed=[])

    A0 = fit_result['params_dict']['A']
    B0 = fit_result['params_dict']['B']

        # plotting the fitted function
    plot.plot_fit1d(fit_result, np.linspace(min(x1),max(x1) , 100), ax=ax, plot_data=True)
    ax.set_xlabel('position in micro meter')
    ax.set_ylabel('Average ESR Frequency')
    plt.title('fit: average frequence and position')

    plt.savefig(folder+'/average_frequence_and_position.png',format='png')


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



    

