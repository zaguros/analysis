"""
This file provides tools to analyze Spin photon correlations as a function of the superposition angle 'theta'
Functions should be executed on the computer that stores the PQ data (othewise use pq_folder = 'xxxx' when calling instances of purify_pq)
Based on the analysis class purify_pq

NK 2017
"""
import numpy as np
from analysis.lib.purification import purify_pq as ppq; reload(ppq)
import Analysis_params_SCE as analysis_params; reload(analysis_params)
from analysis.lib.pq import pq_tools,pq_plots; reload(pq_plots)
from analysis.lib.tools import toolbox as tb
from matplotlib import pyplot as plt


def analyze_spcorrs(contains,is_remote_lt3_measurement = False,**kw):
    

    #### kws
    plot_filter = kw.pop('plot_filter',False)
    do_plot = kw.pop('do_plot',False)

    analysis_file,filtering_file,ssro_calib_folder,trans = get_data_objects(contains,is_remote_lt3_measurement,**kw)
    sn_lt,st_fltr_c0,st_fltr_c1 = temporal_filtering(filtering_file,plot_filter = plot_filter)

    if analysis_params.SPCorr_settings['photon_channel'] == 0:
        analysis_file.sn_filtered = sn_lt[st_fltr_c0]
    elif analysis_params.SPCorr_settings['photon_channel'] == 1:
        analysis_file.sn_filtered = sn_lt[st_fltr_c1]
    else:
        analysis_file.sn_filtered = sn_lt[np.logical_or(st_fltr_c0,st_fltr_c1)]
    
    #################################
    ##### electron RO filtering #####
    #################################

    adwin_filter,adwin_syncs = analysis_file.filter_adwin_data_from_pq_syncs(analysis_file.sn_filtered) ## adwin syncs within window
    results = analysis_file.agrp['ssro_results'].value
    adwin_filt_bool = np.in1d(range(len(results)),adwin_filter) ### convert the adwin filter to boolean
    analysis_file.get_time_filtered_RO_results(adwin_filt_bool)
    analysis_file.get_e_ROC(ssro_calib_folder,trans = trans)
    analysis_file.get_sweep_pts()

    ###### plotting ######
    if do_plot:
        analysis_file.plot_results_vs_sweepparam()

    else: 
        return analysis_file


def sweep_analysis_parameters(contains,parameter_name, parameter_range,is_remote_lt3_measurement = False,**kw):
    
    plot_tail = kw.pop('plot_tail', False)
    
    analysis_file,filtering_file,ssro_calib_folder,trans = get_data_objects(contains,is_remote_lt3_measurement,**kw)
    x,y,z_list,u_z_list,tail_list,u_tail_list= get_sweep_analysis_results(analysis_file,filtering_file,parameter_name,parameter_range,ssro_calib_folder,trans)
    xlabel = analysis_file.sweep_name

    ylim = (-0.05,1.05)#kw.get('ylim', (-0.05, 1.05))
    fig,ax = plt.subplots(1,1)

    ax = plt.subplot(111)
    ax.set_title(analysis_file.timestamp+'\n'+analysis_file.measurementstring + '\n' + 'parameter_name: ' + parameter_name)
    

    for i,z,u_z,yparam in zip(range(len(z_list)),z_list,u_z_list,y):
        ax.errorbar(x,z,u_z,label = str(yparam),fmt = '-o',color = (float(i)/len(z_list),0.1,0.1))
        

    ## formatting
    ax.axhspan(0,1,fill=False,ls='dotted')
    ax.set_ylabel(r'$F(|0\rangle)$')
    ax.set_xlabel(analysis_file.sweep_name)
    ax.set_ylim(ylim)

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    fig2,ax2 = plt.subplots(1,1)
    
    cols = ax2.imshow(z_list,interpolation = 'none',extent=[np.amin(x),np.amax(x),np.amin(y),np.amax(y)],
               aspect = 'equal',origin='lower')
    cax = fig2.add_axes([0.82, 0.15, 0.052, 0.7])
    ax2.set_aspect((np.amax(x)-np.amin(x))/(np.amax(y)-np.amin(y)))
    ax2.set_xlabel(analysis_file.sweep_name)
    ax2.set_ylabel('window start (ps)')
    fig2.colorbar(cols,cax=cax,label=r'$F(|0\rangle)$',orientation='vertical')
    
    # plt.show()


    if plot_tail:
        fig3,ax3 = plt.subplots(1,1)

        for i,t,u_t,yparam in zip(range(len(z_list)),tail_list,u_tail_list,y):
            ax3.errorbar(x,t,u_t,label = str(yparam),fmt = '-o',color = (float(i)/len(z_list),0.1,0.1))

        ## formatting
        ax3.set_ylabel(r'Tail counts *10^4')
        ax3.set_xlabel(analysis_file.sweep_name)
        ax3.set_title(analysis_file.timestamp+'\n'+analysis_file.measurementstring + '\n' + 'parameter_name: ' + parameter_name)



        ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        ### density plot
        fig4,ax4 = plt.subplots(1,1)
        cols = ax4.imshow(tail_list,interpolation = 'none',extent=[np.amin(x),np.amax(x),np.amin(y),np.amax(y)],
                   aspect = 'equal',origin='lower')
        cax = fig4.add_axes([0.82, 0.15, 0.052, 0.7])
        ax4.set_aspect((np.amax(x)-np.amin(x))/(np.amax(y)-np.amin(y)))
        ax4.set_xlabel(analysis_file.sweep_name)
        ax4.set_ylabel('window start (ps)')
        fig4.colorbar(cols,cax = cax,label = 'tail_cts * 10^4',orientation='vertical')


        return fig,fig2,fig3,fig4
    return fig,fig2

def get_sweep_analysis_results(analysis_file,filtering_file,parameter_name,parameter_range,ssro_calib_folder,trans):
    
    

    analysis_file.get_sweep_pts()
    x = analysis_file.sweep_pts
    y = parameter_range
    z = []
    u_z = []
    tail = []
    u_tail = []


    for p in parameter_range:

        analysis_params.SPCorr_settings[parameter_name] = p
        ### do one analysis cycle for the chosen parameter setting

        #  temporal filtering
        sn_lt,st_fltr_c0,st_fltr_c1 = temporal_filtering(filtering_file)

        if analysis_params.SPCorr_settings['photon_channel'] == 0:
            analysis_file.sn_filtered = sn_lt[st_fltr_c0]
        elif analysis_params.SPCorr_settings['photon_channel'] == 1:
            analysis_file.sn_filtered = sn_lt[st_fltr_c1]
        else:
            analysis_file.sn_filtered = sn_lt[np.logical_or(st_fltr_c0,st_fltr_c1)]
        

        #RO filtering and read-out correction
        adwin_filter,adwin_syncs = analysis_file.filter_adwin_data_from_pq_syncs(analysis_file.sn_filtered) ## adwin syncs within window
        results = analysis_file.agrp['ssro_results'].value
        adwin_filt_bool = np.in1d(range(len(results)),adwin_filter) ### convert the adwin filter to boolean

        analysis_file.get_time_filtered_RO_results(adwin_filt_bool)
        analysis_file.get_e_ROC(ssro_calib_folder,trans = trans)
        analysis_file.get_sweep_pts()

        ### store sweep results
        z.append(analysis_file.p0);u_z.append(analysis_file.u_p0)

        analysis_file = get_tail_counts(analysis_file)
        tail.append(analysis_file.tail_cts); u_tail.append(analysis_file.u_tail_cts)

    #restore analysis parameters
    reload(analysis_params)
    return x,y,z,u_z,tail,u_tail

### helper functions

def get_data_objects(contains,is_remote_lt3_measurement,**kw):


    analysis_computer = kw.pop('analysis_computer', 'lt4')
    ssro_calib_contains = kw.pop('ssro_calib_contains', 'SSROCalib')

    if analysis_computer == 'lt4':

        folder=tb.latest_data(contains,**kw)
        a = ppq.purifyPQAnalysis(folder, hdf5_mode='r')

        if is_remote_lt3_measurement:
            lt3_folder = tb.latest_data(contains, folder ='Z:\data')
            b = ppq.purifyPQAnalysis(lt3_folder, hdf5_mode='r')
            analysis_file = b
            filtering_file = a
            ssro_calib_folder  = tb.latest_data(ssro_calib_contains, folder ='Z:\data')
            if 'p' in a.g.attrs['electron_transition']:
                trans = 'msp1'
            else:
                trans = 'msm1'

        else:
            analysis_file = a
            filtering_file = a
            ssro_calib_folder = tb.latest_data(ssro_calib_contains)
            trans = None  

    if analysis_computer == 'lt3_analysis':

        folder=tb.latest_data(contains,folder= r'X:\data',**kw)
        a = ppq.purifyPQAnalysis(folder, hdf5_mode='r')

        if is_remote_lt3_measurement:
            lt3_folder = tb.latest_data(contains, folder ='Y:\data')
            b = ppq.purifyPQAnalysis(lt3_folder, hdf5_mode='r')
            analysis_file = b
            filtering_file = a
            ssro_calib_folder  = tb.latest_data(ssro_calib_contains, folder ='Y:\data')
            if 'p' in a.g.attrs['electron_transition']:
                trans = 'msp1'
            else:
                trans = 'msm1'

        else:
            analysis_file = a
            filtering_file = a
            ssro_calib_folder = tb.latest_data(ssro_calib_contains,folder= r'X:\data')
            trans = None  

    return analysis_file,filtering_file,ssro_calib_folder,trans


def temporal_filtering(filtering_file, **kw):
    """ 
    Apply temporal filtering
    Can manually override the dict by passing the arguement temporal_filter_params in the form {'param (e.g. st_start)' : value}
    """
    plot_filter = kw.pop('plot_filter',False)
    filter_params = kw.pop('temporal_filter_params',{})
    if len(filter_params):
        old_params = analysis_params.SPCorr_settings
    for key,value in filter_params.iteritems():
        analysis_params.SPCorr_settings[key] = value

    if plot_filter:
        plot_ph_hist_and_fltr(filtering_file)

    photon_channel =  analysis_params.SPCorr_settings['photon_channel'] # 0 or 1; 2 means both HH detectors
    st_start = analysis_params.SPCorr_settings['st_start']
    st_len       =  analysis_params.SPCorr_settings['st_len']
    ch1_offset = analysis_params.SPCorr_settings['ch1_offset']

    if len(filter_params):
        analysis_params.SPCorr_settings = old_params

    sp_lt = filtering_file.pqf['/PQ_special-1'].value # Special: marker channel
    ch_lt = filtering_file.pqf['/PQ_channel-1'].value # Channel: photon channel
    sn_lt = filtering_file.pqf['/PQ_sync_number-1'].value # Sync number: the number of the last sync signal
    st_lt = filtering_file.pqf['/PQ_sync_time-1'].value # Sync time: time that has passed since the last sync signal
    
    uniqueVals, uniqueCounts = np.unique(sn_lt[(sp_lt == 0)], return_counts  = True) # Remove syncs with two clicks
    dupVals = uniqueVals[uniqueCounts>1]
    uniqueSync = np.logical_not(np.in1d(sn_lt,dupVals))
    st_fltr_c0 = uniqueSync & (ch_lt == 0)&(st_lt > st_start)  & (st_lt < (st_start  + st_len)) & (sp_lt == 0)  
    st_fltr_c1 = uniqueSync & (ch_lt == 1)&(st_lt > st_start+ch1_offset)  & (st_lt < (st_start  + st_len+ch1_offset)) & (sp_lt == 0) 

    return sn_lt,st_fltr_c0,st_fltr_c1


def get_tail_counts(a):
    """ 
    Added by PH. overrides the code from purify_pq, as in this case we keep trying until we get a success - different number of AWG reps per sweep pt. Need to calculate all this.
    """
    counted_awg_reps = a.agrp['counted_awg_reps'].value

    awg_reps_per_attempt = counted_awg_reps - np.append([0],counted_awg_reps[:-1])
    sweep_length = a.g.attrs['sweep_length']
    repetitions = a.g.attrs['repetitions']
    ### PH - added to make sure that has a full set of repetitions
    awg_reps_per_attempt = awg_reps_per_attempt[:(sweep_length*(len(awg_reps_per_attempt)/sweep_length))]
    repsPerClick = np.sum(np.reshape(awg_reps_per_attempt,[-1,sweep_length]),axis=0)

    counts_per_pt = a.no_of_ro
    a.tail_cts = 10**4 * counts_per_pt/repsPerClick
    a.tail_cts_u = 10**4 * np.sqrt(counts_per_pt)/repsPerClick
    return a

def plot_ph_hist_and_fltr(a):
    st_start = analysis_params.SPCorr_settings['st_start']
    st_len = analysis_params.SPCorr_settings['st_len']
    ch1_offset = analysis_params.SPCorr_settings['ch1_offset']
    f,(ax0,ax1) = plot_marker_filter_comparison(a,
                          mrkr_chan = 1,
                          start = st_start-20e3,
                          length= st_len+100e3,
                          hist_binsize = 1e2,save = False,log=True,ret=True)
    ax0.vlines(np.array([st_start,st_start+st_len])/1e3,0,1000,color='r',lw=2)
    ax1.vlines(np.array([st_start+ch1_offset,st_start+st_len+ch1_offset])/1e3,0,1000,color= 'r',lw=2)
    # ax1.set_xlim([(st_start-20e3+ch1_offset)*1e-3,(st_start+st_len+40e3+ch1_offset)*1e-3])
    



def plot_marker_filter_comparison(a,mrkr_chan = 2,ret=False,**kw):
    pqf = a.pqf

    # get the PLU marked photons first
    is_ph_ch0, is_ph_ch1 = pq_tools.get_photons(pqf)
    is_ph = is_ph_ch0 | is_ph_ch1
    is_ph_with_PLU_mrkr = is_ph & pq_tools.filter_marker(pqf, mrkr_chan)

    # sp_lt = pqf['/PQ_special-1'].value # Special: marker channel
    # sn_lt = pqf['/PQ_sync_number-1'].value # Sync number: the number of the last sync signal
   
    # uniqueVals, uniqueCounts = np.unique(sn_lt[(sp_lt == 0)], return_counts  = True)
    # dupVals = uniqueVals[uniqueCounts>1]
    # uniqueSync = np.logical_not(np.in1d(sn_lt,dupVals))
    # is_ph_with_PLU_mrkr = is_ph_with_PLU_mrkr & uniqueSync

    # adwin_filter,adwin_syncs = a.filter_adwin_data_from_pq_syncs(sn_lt)
    # reps = (a.agrp['counted_awg_reps'].value)
    # LDE_attempts = (((reps - np.append([0],reps[:-1]))-1) % (a.g.attrs['LDE_attempts'])) + 1
    # LDE_attempts = LDE_attempts[adwin_filter]
    
    # is_ph_with_PLU_mrkr = is_ph_with_PLU_mrkr & (LDE_attempts < a.g.attrs['LDE_attempts'])

    if ret:
        return pq_plots.plot_photon_hist_filter_comparison(pqf,fltr =is_ph_with_PLU_mrkr,**kw)
    else:
        pq_plots.plot_photon_hist_filter_comparison(pqf,fltr =is_ph_with_PLU_mrkr,**kw)