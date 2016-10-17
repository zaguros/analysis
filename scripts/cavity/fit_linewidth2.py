"""
Script to analyze the data obtained in a sweep of the montana sync delays measurement
it is based on the general fit formulas in analysis_linewidth_fsr 
- SvD 28-1-2016
"""
import os
import h5py
import numpy as np
from analysis.lib.tools import toolbox as tb
from matplotlib import pyplot as plt

from analysis.scripts.cavity import peakdetect as pd; reload(pd)

from analysis.scripts.cavity import cavity_general_analysis as cga; reload(cga)

from analysis.scripts.cavity import cavity_fitting as cf; reload(cf)

from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot

#def fit_linewidth(folder = tb.latest_data('piezo_scan'),threshold=0.01,delta=0.03):
def fit_linewidth(folder = None,threshold=0.01,delta=0.03, **kw):

    print folder
    show_plots = kw.pop('show_plots',True) #if the keyword show_plots exists, set to the given value, otherwise default is True
    a = cga.cavity_analysis(folder)
    print a.f
    a.get_x_pts() #the x points of the laser of voltage scan
    a.get_data()

    #x_datas contains nr_repetitions copies of a.x_pts
    x_datas = np.full((a.nr_scans,len(a.x_pts)),a.x_pts)  

    single_scan_lw = np.zeros(a.nr_scans)
    u_single_scan_lw = np.zeros(a.nr_scans)
 
    #new_delays = np.zeros(len(a.sweep_pts)*len(a.nr_scans))
    name = ''#reps_0_1_10_11_20_21_30_31_40_41'

    print a.nr_scans
    print x_datas
    print a.data


    for i in np.arange(a.nr_scans):
        single_scan_lw[i],u_single_scan_lw[i]=cf.fit_piezo_plots(folder,x_datas[i],a.data[i],
        show_plots = True, threshold=threshold, averaging = False)

    return single_scan_lw, u_single_scan_lw #SvD: you want to return not the index i, but the whole array

#-------------------------------------------------------------------

n_xticks = 8
n_yticks = 8
EOM_freq = 6  #GHz
EOM_ON = True
filename='LaserScanData'

#linewidth fitting for 3 lorentzians for a single scan
def fit_3_lorentzians(indexmax,date,g_a1,g_A1,g_x01,g_gamma1,g_dx,g_A2,folder=None,fixed=[],threshold=0.01,delta=0.03, **kw):
    print folder
    show_plots = kw.pop('show_plots',True) #if the keyword show_plots exists, set to the given value, otherwise default is True
    a = cga.cavity_analysis(folder)
    print a.f
    a.get_x_pts() #the x points of the laser of voltage scan
    a.get_data()

    #x_datas contains nr_repetitions copies of a.x_pts
    #x_datas = np.full((a.nr_scans,len(a.x_pts)),a.x_pts)  

    single_scan_lw = np.zeros(a.nr_scans)
    u_single_scan_lw = np.zeros(a.nr_scans)
 
    #new_delays = np.zeros(len(a.sweep_pts)*len(a.nr_scans))
    name = ''#reps_0_1_10_11_20_21_30_31_40_41'

    #Now we check to see that we have the correct voltage and amplitude points. 
#    print a.x_pts
#    print a.data

    #We want to analyze only one of the peaks that pops up for the linewidth. So we select a portion of the data and do the analysis on this
    #First we have to find the length of the data
    print ('Length of Data')
    print len(a.x_pts)

    #The entire range of the data
    x=a.x_pts
    y=a.data[0]
    print len(x)
    print len(y)

    #Plot the entire range of the data
    fig,ax = plt.subplots(figsize=(14,8))
    ax.plot(x,y)
    ax.set_xlabel("Voltage (V)]", fontsize = 14)
    ax.set_ylabel("Intensity (a.u.)", fontsize = 14)

    ax.set_title('Raw'+filename+date)
    plt.savefig(os.path.join(folder,filename+'_'+date+'.png'))

    #Selecting one peak of the data for the proper fit 
    #indexmax=y.argmax(axis=0)
    #print indexmax
    #indexmax=30
    dataleft=indexmax-(0.05*len(x))
    dataright=indexmax+(0.05*len(x))
    x=x[dataleft:dataright]
    y=y[dataleft:dataright]
    print len(x)
    print len(y)


    #x = 1.e3*np.a.x_pts[0]#[3000:] # multiplied by factor 1.e3  to get ms  
    #y = np.a#[3000:]

    #Plot just one of the peaks 
    fig,ax = plt.subplots(figsize=(8,4))
    ax.plot(x,y)
    ax.set_xlabel("Voltage (V)]", fontsize = 14)
    ax.set_ylabel("Intensity (a.u.)", fontsize = 14)
    #ax.set_xlim(x[0],x[-1])
    #X_min_freq = g_x01-ax.get_xlim()[0]
    #X_max_freq = g_x01+ax.get_xlim()[-1]
    #xticklabels = np.linspace(X_min_freq,X_max_freq,n_xticks)
    #ax.set_xticks(xticks)

    ax.set_title('Select Data'+filename+date)
    plt.savefig(os.path.join(folder,filename+'Select Data'+'_'+date+'.png'))

    print 'fitting data to 3 lorentzians'
    p0, fitfunc, fitfunc_str = common.fit_3lorentz_symmetric(g_a1, g_A1, g_x01, g_gamma1, g_dx, g_A2)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)

#    # x01 = fit_result['params_dict']['x01']
    dx = fit_result['params_dict']['dx']
    gamma1 = fit_result['params_dict']['gamma1']
#    # gamma2 = fit_result['params_dict']['gamma2']
    u_gamma1 = fit_result['error_dict']['gamma1']
#    # u_gamma2 = fit_result['error_dict']['gamma2']

    scaling = EOM_freq/dx #scale the known EOM freq with the separation here.
    linewidth = gamma1*scaling #scale the linewidth to get linewidht in frequency
    u_linewidth = u_gamma1*scaling
    linewidth_string = 'gamma = '+str(round(linewidth,2))+'+-'+str(round(u_linewidth,3))+'GHz'
    print linewidth_string

 #Plotting

    fig,ax = plt.subplots(figsize=(8,4))
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],len(x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)
    ax.set_xlabel("Frequency [GHz]", fontsize = 14)
    ax.set_ylabel("Intensity (a.u.)", fontsize = 14)
    ax.set_xlim(x[0],x[-1])
    print ax.set_xlim(x[0],x[-1])
    xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],n_xticks)
    #rescaling for x-axis in GHz
    X_min_freq = ax.get_xlim()[0]*scaling
    X_max_freq = ax.get_xlim()[-1]*scaling
    print g_x01*scaling
    print X_min_freq
    print X_max_freq
    print scaling

    xticklabels = np.linspace(X_min_freq,X_max_freq,n_xticks)

    xticklabels_round=[]
    for j in xticklabels:
      round_ = round(j,3)
      xticklabels_round = np.append(xticklabels_round,round_)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels_round)

    ax.set_title('Fitted'+filename+date+'\n'+linewidth_string)
    plt.savefig(os.path.join(folder,filename+'_'+date+'_fit.png'))

##---------------------------------------------------------------------------------------

n_xticks = 8
n_yticks = 8
EOM_freq = 6  #GHz
EOM_ON = True
filename='LaserScanData'

#def fit_3_lorentzians_params(date):

#linewidth fitting for 3 lorentzians for a single scan
def fit_3_lorentzians(indexmax,date,g_a1,g_A1,g_x01,g_gamma1,g_dx,g_A2,folder=None,fixed=[],threshold=0.01,delta=0.03, **kw):
    print folder
    show_plots = kw.pop('show_plots',True) #if the keyword show_plots exists, set to the given value, otherwise default is True
    a = cga.cavity_analysis(folder)
    print a.f
    a.get_x_pts() #the x points of the laser of voltage scan
    a.get_data()

    #x_datas contains nr_repetitions copies of a.x_pts
    #x_datas = np.full((a.nr_scans,len(a.x_pts)),a.x_pts)  

    single_scan_lw = np.zeros(a.nr_scans)
    u_single_scan_lw = np.zeros(a.nr_scans)
 
    #new_delays = np.zeros(len(a.sweep_pts)*len(a.nr_scans))
    name = ''#reps_0_1_10_11_20_21_30_31_40_41'

    #Now we check to see that we have the correct voltage and amplitude points. 
#    print a.x_pts
#    print a.data

    #We want to analyze only one of the peaks that pops up for the linewidth. So we select a portion of the data and do the analysis on this
    #First we have to find the length of the data
    print ('Length of Data')
    print len(a.x_pts)

    #The entire range of the data
    x=a.x_pts
    y=a.data[0]
    print len(x)
    print len(y)

    #Plot the entire range of the data
    fig,ax = plt.subplots(figsize=(14,8))
    ax.plot(x,y)
    ax.set_xlabel("Voltage (V)]", fontsize = 14)
    ax.set_ylabel("Intensity (a.u.)", fontsize = 14)

    ax.set_title('Raw'+filename+date)
    plt.savefig(os.path.join(folder,filename+'_'+date+'.png'))

    #Selecting one peak of the data for the proper fit 
    #indexmax=y.argmax(axis=0)
    #print indexmax
    #indexmax=30
    dataleft=indexmax-(0.05*len(x))
    dataright=indexmax+(0.05*len(x))
    x=x[dataleft:dataright]
    y=y[dataleft:dataright]
    print len(x)
    print len(y)


    #x = 1.e3*np.a.x_pts[0]#[3000:] # multiplied by factor 1.e3  to get ms  
    #y = np.a#[3000:]

    #Plot just one of the peaks 
    fig,ax = plt.subplots(figsize=(8,4))
    ax.plot(x,y)
    ax.set_xlabel("Voltage (V)]", fontsize = 14)
    ax.set_ylabel("Intensity (a.u.)", fontsize = 14)
    #ax.set_xlim(x[0],x[-1])
    #X_min_freq = g_x01-ax.get_xlim()[0]
    #X_max_freq = g_x01+ax.get_xlim()[-1]
    #xticklabels = np.linspace(X_min_freq,X_max_freq,n_xticks)
    #ax.set_xticks(xticks)

    ax.set_title('Select Data'+filename+date)
    plt.savefig(os.path.join(folder,filename+'Select Data'+'_'+date+'.png'))

    print 'fitting data to 3 lorentzians'
    p0, fitfunc, fitfunc_str = common.fit_3lorentz_symmetric(g_a1, g_A1, g_x01, g_gamma1, g_dx, g_A2)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)

#    # x01 = fit_result['params_dict']['x01']
    dx = fit_result['params_dict']['dx']
    gamma1 = fit_result['params_dict']['gamma1']
#    # gamma2 = fit_result['params_dict']['gamma2']
    u_gamma1 = fit_result['error_dict']['gamma1']
#    # u_gamma2 = fit_result['error_dict']['gamma2']

    scaling = EOM_freq/dx #scale the known EOM freq with the separation here.
    linewidth = gamma1*scaling #scale the linewidth to get linewidht in frequency
    u_linewidth = u_gamma1*scaling
    linewidth_string = 'gamma = '+str(round(linewidth,2))+'+-'+str(round(u_linewidth,3))+'GHz'
    print linewidth_string

 #Plotting

    fig,ax = plt.subplots(figsize=(8,4))
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],len(x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)
    ax.set_xlabel("Frequency [GHz]", fontsize = 14)
    ax.set_ylabel("Intensity (a.u.)", fontsize = 14)
    ax.set_xlim(x[0],x[-1])
    print ax.set_xlim(x[0],x[-1])
    xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],n_xticks)
    #rescaling for x-axis in GHz
    X_min_freq = ax.get_xlim()[0]*scaling
    X_max_freq = ax.get_xlim()[-1]*scaling
    print g_x01*scaling
    print X_min_freq
    print X_max_freq
    print scaling

    xticklabels = np.linspace(X_min_freq,X_max_freq,n_xticks)

    xticklabels_round=[]
    for j in xticklabels:
      round_ = round(j,3)
      xticklabels_round = np.append(xticklabels_round,round_)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels_round)

    ax.set_title('Fitted'+filename+date+'\n'+linewidth_string)
    plt.savefig(os.path.join(folder,filename+'_'+date+'_fit.png'))

def find_average_linewidth(from_folder='20160204145330', until_folder='20160204145336'):

    avg_lw = np.array([])
    u_avg_lw = 0

    i = int(until_folder)+1
    print i

    while i > int(from_folder):
        print i
        lws, u_lws = fit_linewidth(folder = tb.latest_data('piezo_scan', older_than = str(i)),threshold=0.01,delta=0.03)
        avg_lw = np.append(avg_lw,lws)
        print avg_lw
        i = i - 1

    print 'Voila! Hvala!'

#-------------------------------------------------------------------

    # fig,ax = plt.subplots(figsize=(6,4.7))
    # ax.errorbar(np.arange(a.nr_scans), avg_lws, yerr=u_avg_lws)
    # ax.legend()
    # ax.set_xlabel('scan number',fontsize=14)
    # ax.set_ylabel('average linewidth (nm)',fontsize=14)
    # ax.tick_params(axis = 'both', which = 'major',labelsize = 14)
    # ax.set_xlim(1,a.nr_scans)
    # ax.set_ylim(0,max(avg_lws))

    # ax.set_title(folder )

    # fig.savefig(folder +'/'+"lw_vs_scan_number .png")
    # plt.show()
    # plt.close()
    

# if __name__ == '__main__':
#     find_average_linewidth()