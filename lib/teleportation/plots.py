import os

from matplotlib import pyplot as plt
import numpy as np

import files, events, settings, stats

##############################################################################
### Photon histograms in time
##############################################################################

def _plot_photon_hist(ax, h, b, log=True, **kw):
    label = kw.pop('label', '')

    _h = h.astype(float)
    _h[_h<=1e-1] = 1e-1
    _h = np.append(_h, _h[-1])
           
    ax.plot(b, _h, drawstyle='steps-post', label=label)
    if log:
        ax.set_yscale('log')
    ax.set_xlabel('time (ns)')
    ax.set_ylabel('events')
    ax.set_ylim(bottom=0.1)
    ax.set_xlim(min(b), max(b))

def plot_photon_hist(fp, **kw):    
    ret = kw.pop('ret', 'subplots')

    (h0, b0), (h1, b1) = events.get_photon_hist(fp, **kw)
   
    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,8))
    _plot_photon_hist(ax0, h0, b0)
    _plot_photon_hist(ax1, h1, b1)

    ax0.set_title('photons channel 0')
    ax1.set_title('photons channel 1')
    
    fig.suptitle(files.get_msmt_header(fp) + ' -- Photon histogram')
    
    if ret == 'subplots':
        return fig, (ax0, ax1)
    
def plot_photon_hist_filter_comparison(fp, fltr, **kw):
    ret = kw.pop('ret', 'subplots')
    
    (h0, b0), (h1, b1) = events.get_photon_hist(fp, **kw)
    (h0f, b0f), (h1f, b1f) = events.get_photon_hist(fp, fltr=fltr, **kw)
    
    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,8))
    _plot_photon_hist(ax0, h0, b0, label='unfiltered')
    _plot_photon_hist(ax0, h0f, b0f, label='filtered')

    _plot_photon_hist(ax1, h1, b1, label='unfiltered')
    _plot_photon_hist(ax1, h1f, b1f, label='filtered')
    
    ax0.set_title('photons channel 0')
    ax1.set_title('photons channel 1')

    fig.suptitle(files.get_msmt_header(fp) + ' -- Photon histogram, comparison w/ filter')
        
    if ret=='subplots':
        return fig, (ax0, ax1)


def plot_tail(fp, **kw):
    ret = kw.pop('ret', 'subplots')
    
    (h0, b0), (h1, b1) = events.get_photon_hist(fp, **kw)
    
    fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,8))
    _plot_photon_hist(ax0, h0, b0)
    _plot_photon_hist(ax1, h1, b1)
    
    ax0.set_xlim(settings.CH0_START-5, settings.CH0_START + 100)
    ax1.set_xlim(settings.CH1_START-5, settings.CH1_START + 100)
    ax0.set_title('tail channel 0')
    ax1.set_title('tail channel 1')
    
    fitres0, fitres1 = events.fit_tail(fp)
    
    ax0.plot(fitres0['fit_x'], fitres0['fit_y'], lw=2)
    ax1.plot(fitres1['fit_x'], fitres1['fit_y'], lw=2)
    
    ax0.text(settings.CH0_START, 1, fitres0['info'], va='bottom', ha='left', color='k', size='x-small',
             bbox = {'facecolor' : 'white', 'alpha' : 0.5})
    ax1.text(settings.CH1_START, 1, fitres1['info'], va='bottom', ha='left', color='k', size='x-small',
             bbox = {'facecolor' : 'white', 'alpha' : 0.5})
    
    starts = stats.get_sequence_starts(fp)
    ax0.text(settings.CH0_START+5, max(h0)/5., 
             'tailcounts / shot: {:.2e} ({} x {} starts)'.format(fitres0['tcpsh'], starts, settings.SEQREPS))
    ax1.text(settings.CH1_START+5, max(h1)/5., 
             'tailcounts / shot: {:.2e} ({} x {} starts)'.format(fitres1['tcpsh'], starts, settings.SEQREPS))

    fig.suptitle(files.get_msmt_header(fp) + ' -- Tail')

    if ret=='subplots':
        return fig, (ax0, ax1)


##############################################################################
### CR stat plotting
##############################################################################

def _plot_CR_hist(ax, heights, **kw):
    ax.bar(np.arange(len(heights)), heights, align='center', 
           width=1, color=settings.COLORS[0], lw=0)
    ax.set_xlabel('CR counts')
    ax.set_ylabel('Events')
    ax.set_title(kw.get('title', ''))
    ax.set_xlim(-0.5, len(heights)-0.5 if len(heights) > 0 else 0.5)
    
### plotting of histograms
def plot_CR_hist_total(fp):
    h1, b1 = stats.get_CR_hist_total(fp, 'lt1')
    h2, b2 = stats.get_CR_hist_total(fp, 'lt2')
    
    fig, (ax1, ax2) = files.subplots(1,2, figsize=(12,5))
    _plot_CR_hist(ax1, h1, title='LT1')
    _plot_CR_hist(ax2, h2, title='LT2')
    
    fig.suptitle(files.get_msmt_header(fp) + ' -- CR histogram total')
    folder, _fn = os.path.split(fp)
    fig.savefig(os.path.join(folder, 'CR_histogram_total.png'))
    
def plot_CR_hist_sequence_timeout(fp):
    h1, b1 = stats.get_CR_hist_sequence_timeout(fp, 'lt1')
    h2, b2 = stats.get_CR_hist_sequence_timeout(fp, 'lt2')
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))
    _plot_CR_hist(ax1, h1, title='LT1')
    _plot_CR_hist(ax2, h2, title='LT2')
    
    fig.suptitle(files.get_msmt_header(fp) + ' -- CR histogram sequence timeout')
    folder, _fn = os.path.split(fp)
    fig.savefig(os.path.join(folder, 'CR_histogram_sequence_timeout.png'))


