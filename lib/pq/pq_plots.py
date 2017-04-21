### This is mainly an adaption from Wolfgang's teleportation.plots functions. They all expect a pqf file

import os

from matplotlib import pyplot as plt
import numpy as np

from analysis.lib.pq import pq_tools
from analysis.lib.tools import toolbox as tb

##############################################################################
### Plotting photon histograms
##############################################################################

def _plot_photon_hist(ax, h, b, log=True, **kw):
    label = kw.pop('label', '')
    x_label = kw.pop('x_label','time (ns)')
    _h = h.astype(float)
    _h[_h<=1e-1] = 1e-1
    _h = np.append(_h, _h[-1])
           
    ax.plot(b/1000., _h, drawstyle='steps-post', label=label)
    if log:
        ax.set_yscale('log')
    ax.set_xlabel(x_label)
    ax.set_ylabel('events')
    ax.set_ylim(bottom=0.1)
    ax.set_xlim(min(b/1000.), max(b/1000.))
    ax.grid(True)

def plot_photon_hist(pqf, **kw):    


    ret = kw.pop('ret', 'subplots')

    (h0, b0), (h1, b1) = pq_tools.get_photon_hist(pqf, **kw)
   
    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,8))
    _plot_photon_hist(ax0, h0, b0, **kw)
    _plot_photon_hist(ax1, h1, b1, **kw)
    print np.sum(h0)+np.sum(h1)
    ax0.set_title('photons channel 0')
    ax1.set_title('photons channel 1')

    fp = pq_tools.fp_from_pqf(pqf)
    
    fig.suptitle(tb.get_msmt_header(fp) + ' -- Photon histogram')
    
    if ret == 'subplots':
        return fig, (ax0, ax1)

    
def plot_photon_hist_filter_comparison(pqf, fltr, **kw):
    ret = kw.pop('ret', 'subplots')
    
    (h0, b0), (h1, b1) = pq_tools.get_photon_hist(pqf, **kw)
    (h0f, b0f), (h1f, b1f) = pq_tools.get_photon_hist(pqf, fltr=fltr, **kw)
    
    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,8))
    _plot_photon_hist(ax0, h0, b0, label='unfiltered', **kw)
    _plot_photon_hist(ax0, h0f, b0f, label='filtered', **kw)

    _plot_photon_hist(ax1, h1, b1, label='unfiltered', **kw)
    _plot_photon_hist(ax1, h1f, b1f, label='filtered', **kw)
    

    if kw.get('plot_threshold_ch0', 0) !=0:
        ax0.vlines([kw.get('plot_threshold_ch0', 0)],0,1000,color = 'r')
    if kw.get('plot_threshold_ch1', 0) !=0:
        ax1.vlines([kw.get('plot_threshold_ch1', 0)],0,1000,color='r')

    # if kw.get('window_ch0',[0,0]) != [0,0]:
    #     ax0.

    ax0.set_title('photons channel 0')
    ax1.set_title('photons channel 1')

    ax0.legend()
    ax1.legend()

    fp = pq_tools.fp_from_pqf(pqf)

    fig.suptitle(tb.get_msmt_header(fp) + ' -- Photon histogram, comparison w/ filter')
        
    if ret=='subplots':
        return fig, (ax0, ax1)


def plot_marker_filter_comparison(pqf,mrkr_chan = 2,ret=False,**kw):

    # get the PLU marked photons first
    is_ph_ch0, is_ph_ch1 = pq_tools.get_photons(pqf)
    is_ph = is_ph_ch0 | is_ph_ch1
    is_ph_with_PLU_mrkr = is_ph & pq_tools.filter_marker(pqf, mrkr_chan)
    if ret:
        return plot_photon_hist_filter_comparison(pqf,fltr =is_ph_with_PLU_mrkr,**kw)
    else:
        plot_photon_hist_filter_comparison(pqf,fltr =is_ph_with_PLU_mrkr,**kw)
    


def plot_autocorrelation_histogram(pqf,start = 0,length = 2000,index = 1, pq_device = '', binsize = 1e3,**kw):
    """
    does not exlcude no photon (i.e. markers) events at the moment.
    """
    sync_time_name = pq_device + '/PQ_sync_time-' + str(index)
    tot_time_name =  pq_device + '/PQ_time-' + str(index)
    sync_num_name = pq_device + '/PQ_sync_number-' + str(index)

    tot_time_0 = np.append(np.array([0]),pqf[tot_time_name])
    tot_time_1 = np.append(pqf[tot_time_name],np.array([0]))

    time_diff = (tot_time_1-tot_time_0)/1e6

    binedges = np.arange(start,start+length, binsize)
    h0, b0 = np.histogram(time_diff, bins=binedges)

    fig, ax0 = plt.subplots(1,1, figsize=(12,8))
    _plot_photon_hist(ax0, h0, b0,x_label = 'time (s)',**kw)

    fp = pq_tools.fp_from_pqf(pqf)

    fig.suptitle(tb.get_msmt_header(fp) + ' -- autocorrelation hist in total time')

    

#def plot_tail(fp, **kw):
#    ret = kw.pop('ret', 'subplots')
#    
#    (h0, b0), (h1, b1) = events.get_photon_hist(fp, **kw)
#    
#    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,8))
#    _plot_photon_hist(ax0, h0, b0)
#    _plot_photon_hist(ax1, h1, b1)
#    
#    ax0.set_xlim(settings.CH0_START-5, settings.CH0_START + 100)
#    ax1.set_xlim(settings.CH1_START-5, settings.CH1_START + 100)
#    ax0.set_title('tail channel 0')
#    ax1.set_title('tail channel 1')
#    
#    fitres0, fitres1 = events.fit_tail(fp)
#    
#    ax0.plot(fitres0['fit_x'], fitres0['fit_y'], lw=2)
#    ax1.plot(fitres1['fit_x'], fitres1['fit_y'], lw=2)
#    
#    ax0.text(settings.CH0_START, 1, fitres0['info'], va='bottom', ha='left', color='k', size='x-small',
#             bbox = {'facecolor' : 'white', 'alpha' : 0.5})
#    ax1.text(settings.CH1_START, 1, fitres1['info'], va='bottom', ha='left', color='k', size='x-small',
#             bbox = {'facecolor' : 'white', 'alpha' : 0.5})
#    
#    starts = stats.get_sequence_starts(fp)
#    ax0.text(settings.CH0_START+5, max(h0)/5., 
#             'tailcounts / shot: {:.2e} ({} x {} starts)'.format(fitres0['tcpsh'], starts, settings.SEQREPS))
#    ax1.text(settings.CH1_START+5, max(h1)/5., 
#             'tailcounts / shot: {:.2e} ({} x {} starts)'.format(fitres1['tcpsh'], starts, settings.SEQREPS))
#
#    fig.suptitle(files.get_msmt_header(fp) + ' -- Tail')
#
#    if ret=='subplots':
#        return fig, (ax0, ax1)

def plot_PLU_filter(pqf,chan = 2):

    # get the PLU marked photons first
    is_ph_ch0, is_ph_ch1 = pq_tools.get_photons(pqf)
    is_ph = is_ph_ch0 | is_ph_ch1
    is_ph_with_PLU_mrkr = is_ph & pq_tools.filter_marker(pqf, chan)

    # first window
    fig, (ax0, ax1) = plot_photon_hist_filter_comparison(
        pqf, save=False, fltr=is_ph_with_PLU_mrkr, log=True,
        binedges=settings.PHOTONHIST_BINEDGES)

    ax0.axvline(settings.CH0_START, color='k', linestyle='--')
    ax0.text(settings.CH0_START+1, 1, '{} ns'.format(settings.CH0_START), color='k')
    ax0.legend()
    ax0.set_title('PLU (hardware) filter first pi-pulse ch0')

    ax1.axvline(settings.CH1_START, color='k', linestyle='--')
    ax1.text(settings.CH1_START+1, 1, '{} ns'.format(settings.CH1_START), color='k')
    ax1.legend()
    ax1.set_title('PLU (hardware) filter first pi-pulse ch1')

    fig.savefig(os.path.join(folder, 'PLU_photons_window1.png'))

    # second window
    fig, (ax0, ax1) = plot_photon_hist_filter_comparison(
        pqf, save=False, fltr=is_ph_with_PLU_mrkr, log=True,
        binedges=settings.PHOTONHIST_BINEDGES + settings.PIPULSESEP)

    ax0.axvline(settings.CH0_START+settings.PIPULSESEP, color='k', linestyle='--')
    ax0.text(settings.CH0_START+settings.PIPULSESEP+1, 
        1, '{} ns'.format(settings.CH0_START+settings.PIPULSESEP), color='k')
    ax0.legend()
    ax0.set_title('PLU (hardware) filter second pi-pulse ch0')

    ax1.axvline(settings.CH1_START+settings.PIPULSESEP, color='k', linestyle='--')
    ax1.text(settings.CH1_START+settings.PIPULSESEP+1, 
        1, '{} ns'.format(settings.CH1_START+settings.PIPULSESEP), color='k')
    ax1.legend()
    ax1.set_title('PLU (hardware) filter second pi-pulse ch1')
     
    fig.savefig(os.path.join(folder, 'PLU_photons_window2.png'))


