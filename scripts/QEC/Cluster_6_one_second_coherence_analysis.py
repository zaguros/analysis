import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
from matplotlib import gridspec
reload(common)
reload(fit) 
reload(toolbox)
reload(mbi)



def CosineSum_MBI_data(folder=None,timestamp=[], measurement_name = ['adwindata'],ssro_calib_folder=None, ssro_calib_timestamp =None,
        new_tsmp = '20170612_230600',old_tsmp = '20170614_071700',two_cos=True,title='Cluster',
        x_ticks=np.arange(0,100,5), y_ticks=np.arange(0.2,0.8,0.1),color='b',
        c=1,t=1,frequency = [1,1], offset =0.5, amplitude =[ 0.5,0.5],  phase =[0,0], 
        fixed = [], 
        plot_fit = False, do_print = False, show_guess = True, xlim1=None, xlim2=None):
    ''' 
    Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    folder: allows to specify specific folder for the data(this overwrites the timestamp input)
    ssro_folder: allows to specify specific folder for ssro calib(this overwrites the ssro_timestamp input)
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    [freq, offset, Amplitude, phase] 
    '''
    timestamp=[]
    plot_fit    = True
    show_guess  = False
    exponent=2
    #two_cos=True
    fit_results = []

    cum_pts = 0
    cum_sweep_pts       = np.empty(0)
    cum_p0              = np.empty(0)
    cu                  =np.empty(0) 
    cu_u                =np.empty(0) 
    cum_u_p0            = np.empty(0)
    cum_tau_list            = np.empty(0)

    # new_tsmp = '20170611_165140' ## newer than
    # old_tsmp = '20170611_235700' ## older than

    # new_tsmp = '20170612_005040' ## newer than
    # old_tsmp = '20170612_015800' ## older than

    # new_tsmp = '20170612_230600' ## newer than
    # old_tsmp = '20170614_071700' ## older than

    # new_tsmp = '20170616_190600' ## newer than
    # old_tsmp = '20170617_053700' ## older than

    # new_tsmp = '20170621_183900' ## newer than
    # old_tsmp = '20170621_185800' ## older than

    # new_tsmp = '20170719_161700' ## newer than
    # old_tsmp = '20170719_165500' ## older than

    # new_tsmp = '20170720_171400' ## newer than
    # old_tsmp = '20170720_181700' ## older than

    # new_tsmp = '20170720_182900' ## newer than
    # old_tsmp = '20170720_195700' ## older than

    search_string = 'Carbon'
    while toolbox.latest_data(contains=search_string,
                                            return_timestamp =True,
                                            older_than=old_tsmp,
                                            newer_than=new_tsmp,
                                            raise_exc=False) != False:
        old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                            return_timestamp =True,
                                            older_than=old_tsmp,
                                            newer_than=new_tsmp,
                                            raise_exc=False)

        timestamp.append(old_tsmp)


    timestamp = timestamp[::-1]

    
    if ssro_calib_folder==None:
        if ssro_calib_timestamp == None: 
            ssro_calib_folder = toolbox.latest_data('SSRO')
        else:
            ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
            ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
    print ssro_calib_folder
    
 
    for kk in range(len(timestamp)):
        folder = toolbox.data_from_time(timestamp[kk])

        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC(ssro_calib_folder)
        cum_pts += a.pts
        # temperature = (a.g.attrs['temp']-100)/0.385
        # temperature_list.append(temperature)


        if kk == 0:
            cum_sweep_pts = a.sweep_pts
            cum_p0 = a.p0
            cum_u_p0 = a.u_p0
            #cum_tau_list = a.tau_list
        #elif kk in [5,10,15,21,26,31]: 
        else:
            cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
            cum_p0 = np.concatenate((cum_p0, a.p0))
            cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))
            #cum_tau_list = np.concatenate((cum_tau_list, a.tau_list))



    a.pts   = cum_pts 
    a.sweep_pts = 1.0e3*cum_sweep_pts 
    a.p0    = cum_p0
    a.u_p0  = cum_u_p0

    
    #ax = a.plot_results_vs_sweepparam(ret='ax',fmt='o',color='brown',figsize=(46,12))
    e= a.u_p0
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
        
    # fig = plt.figure(1,figsize=(10,2))
    # ax2 = fig.add_subplot(111)


    # ax2.set_xlabel('Free evolution time (ms)')
    # ax2.set_ylabel('State Fidelity')

    #f = plt.figure(1,figsize=(60,2))
    #f,(ax1,ax2,ax3,ax4) = plt.subplots(1,4,sharey=True, facecolor='w',figsize=(10,3.5))

    fig = plt.figure(figsize=(9, 3.3)) 
    gs = gridspec.GridSpec(1, 4, width_ratios=[3, 3,2,1]) 
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
    ax4 = plt.subplot(gs[3])

    ax1.errorbar(x.flatten(),y.flatten(),yerr=e,fmt='o',label='',color=color,lw=1,markersize=3)
    ax2.errorbar(x.flatten(),y.flatten(),yerr=e,fmt='o',label='',color=color,lw=1,markersize=3)
    ax3.errorbar(x.flatten(),y.flatten(),yerr=e,fmt='o',label='',color=color,lw=1,markersize=3)
    ax4.errorbar(x.flatten(),y.flatten(),yerr=e,fmt='o',label='',color=color,lw=1,markersize=3)



    fit_result = [None]

    #p0, fitfunc, fitfunc_str = common.fit_sum_2cos(offset,amplitude[0],frequency[0],phase[0],amplitude[1],frequency[1],phase[1]) 
    if two_cos==True:
        p0, fitfunc, fitfunc_str = common.fit_gaussian_decaying_2cos(offset,c,t,amplitude[0],
            frequency[0],phase[0],amplitude[1],frequency[1],phase[1]) 
    else:
        p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(offset,amplitude[0], 0,t,
            exponent,frequency[0] ,phase[0] )
    if show_guess:
        ax1.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=fixed)
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],8001), ax=ax1, color=color,
                plot_data=False,print_info = False,lw=1.5)

        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],8001), ax=ax2, color=color,
                plot_data=False,print_info = False,lw=1.5)

        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],8001), ax=ax3, color=color,
                plot_data=False,print_info = False,lw=1.5)
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],8001), ax=ax4, color=color,
                plot_data=False,print_info = False,lw=1.5)



    fit.write_to_file(fit_result,folder,fitname = 'Sum of cosine fit') 


    ## plot data and fit as function of total time
    #plt.xticks(x_ticks)
    # plt.yticks(np.arange(min(y),(max(y)),0.1))
    #plt.yticks(y_ticks)
    ax1.set_xlim(0,31)
    ax2.set_xlim(xlim1,xlim2)
    ax3.set_xlim(609,631)
    ax4.set_xlim(910,920)

    #ax2.set_ylim(0,1)

    #ax2.set_xlabel('Free evolution time (ms)')
    #ax2.set_ylabel('State Fidelity')
    d = .045 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    ax1.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax1.yaxis.tick_left()
    ax1.tick_params(labelright='off')
    ax2.tick_params(labelleft='off')
    ax2.tick_params(labelright='off')
    ax3.tick_params(labelleft='off')
    ax3.tick_params(labelright='off')
    ax4.tick_params(labelleft='off')
    ax4.tick_params(labelright='off')
    #ax4.tick_params(labelbottom='off')
    #ax1.tick_params(labelbottom='off')
    #ax2.tick_params(labelbottom='off')
    #ax3.tick_params(labelbottom='off')

    ax4.spines['left'].set_visible(False)
    plt.sca(ax1)
    plt.xticks([10,20,30])
    plt.yticks([0.4,0.5,0.6,0.7])
    plt.sca(ax2)
    plt.xticks([200,210,220,230])
    plt.yticks([])
    plt.sca(ax3)
    plt.xticks([610,620,630])
    plt.yticks([])
    plt.sca(ax4)
    plt.xticks([910,920])
    plt.yticks([])


    #plt.step(ax1,ax2,xticks=[10,20,30],yticks=[0.4,0.5,0.6,0.7])
    #ax4.yaxis.tick_right()
    #ax4.xaxis.tick_bottom()

    # kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    # ax1.plot((1-d,1+d), (-d,+d), **kwargs)
    # ax1.plot((1-d,1+d),(1-d,1+d), **kwargs)

    # kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    # ax2.plot((-d,+d), (1-d,1+d), **kwargs)
    # ax2.plot((-d,+d), (-d,+d), **kwargs)
    # ax2.plot((1-d,1+d), (-d,+d), **kwargs)
    # ax2.plot((1-d,1+d),(1-d,1+d), **kwargs)

    # kwargs = dict(transform=ax3.transAxes, color='k', clip_on=False)
    # ax3.plot((-d,+d), (1-d,1+d), **kwargs)
    # ax3.plot((-d,+d), (-d,+d), **kwargs)
    # ax3.plot((1-d,1+d), (-d,+d), **kwargs)
    # ax3.plot((1-d,1+d),(1-d,1+d), **kwargs)

    # kwargs.update(transform=ax4.transAxes)  # switch to the bottom axes
    # ax4.plot((-d,+d), (1-d,1+d), **kwargs)
    # ax4.plot((-d,+d), (-d,+d), **kwargs)
        
    folder='C:\Users\TUD277931\Dropbox\TaminiauLab\Projects\Coherence in multi-qubit systems\Paper\Figures'
    plt.savefig(os.path.join(folder, title + '.pdf'),
    format='pdf',bbox_inches='tight',pad_inches=0.2,transparent=True)

    # return fit_result
