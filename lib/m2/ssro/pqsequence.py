import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, sequence
from analysis.lib.math import error
from analysis.lib.m2 import m2
from analysis.lib.tools import toolbox

class PQSequenceAnalysis(sequence.SequenceAnalysis):

    def get_photons(self):
        """
        returns two filters (1d-arrays): whether events are ch0-photons/ch1-photons
        """
        channel = self.f['/PQ_channel-1'].value
        special = self.f['/PQ_special-1'].value

        is_not_special = special==0
        is_channel_0 = channel==0
        is_channel_1 = channel==1

        is_photon_0 = np.logical_and(is_not_special, is_channel_0)
        is_photon_1 = np.logical_and(is_not_special, is_channel_1)

        return is_photon_0, is_photon_1

    def get_tail_vs_sweep(self, syncs_per_sweep, channel, start, tail_length, pq_binsize=1e-3, hist_binsize=1.0, verbose= False):#ns
        is_ph = self.get_photons()[channel]
        sync_time_ns = self.f['/PQ_sync_time-1'].value * pq_binsize
        sync_nrs=self.f['/PQ_sync_number-1'].value  

        if not hasattr(self, 'sweep_pts'):
            self.sweep_pts = np.arange(len(self.ssro_results)) + 1
            self.sweep_name = 'sweep parameter'
        self.sweep_length = len(self.sweep_pts)
        self.syncs_per_sweep = syncs_per_sweep
        sweep_idxs=np.mod(np.floor(sync_nrs/syncs_per_sweep),self.sweep_length)
        
        hist_bins = np.arange(start-hist_binsize*.5,start+1*tail_length+hist_binsize,hist_binsize)
        
        self.tail_hist_h=np.zeros((self.sweep_length,len(hist_bins)-1))
        
        st_fltr = (start  <= sync_time_ns) &  (sync_time_ns< (start + tail_length))
        if verbose:
            print 'total_photons in channel', channel, ':', len(sync_time_ns[np.where(is_ph)])  
            print 'total_photons in window:', len(sync_time_ns[np.where(is_ph & st_fltr)]) 
        valid_tail_idxs = sweep_idxs[np.where(is_ph & st_fltr)]
        if verbose:

            print 'total_sweeps in window:', len(valid_tail_idxs) 
            print 'total ph in window with sweep element 0:', len(np.where(valid_tail_idxs==0)[0])
            print 'div factor:', (self.reps*syncs_per_sweep/self.sweep_length)
        self.tail_cts_per_sweep_idx=np.zeros(self.sweep_length)
        
        for sweep_idx in range(self.sweep_length):
            self.tail_cts_per_sweep_idx[sweep_idx]= \
                                float(len(np.where(valid_tail_idxs==sweep_idx)[0])) / (self.reps*syncs_per_sweep/self.sweep_length)
            self.tail_hist_h[sweep_idx], self.tail_hist_b = \
                                np.histogram(sync_time_ns[np.where(is_ph & (sweep_idxs == sweep_idx))], bins=hist_bins)

    
    def plot_tail_cts_vs_sweep(self, name='', save=True, **kw):
        ret = kw.get('ret', None)
        ax = kw.get('ax', None)

        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            save = False
        y_err=np.sqrt(self.tail_cts_per_sweep_idx)/np.sqrt((self.reps*self.syncs_per_sweep/self.sweep_length))
        ax.errorbar(self.sweep_pts, self.tail_cts_per_sweep_idx*1e4, fmt='o',yerr=y_err*1e4)
    
        ax.set_xlabel(self.sweep_name)
        ax.set_ylabel('Tail counts per shot * 10^-4')

        if save:
            self.save_fig_incremental_filename(fig,'Tail_counts_vs_sweep')
        
        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig

    def plot_tail_hists_vs_sweep(self, name='', save=True, **kw):
        ret = kw.get('ret', None)
        ax = kw.get('ax', None)

        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            save = False

        ax.imshow(np.log10(self.tail_hist_h), interpolation='none', cmap='gray', origin='lower', aspect='auto',
                extent=[self.tail_hist_b[0],self.tail_hist_b[-1],self.sweep_pts[0],self.sweep_pts[-1]])
        #ax.colorbar()
        ax.set_xlabel('Time after sync [ns]')
        ax.set_ylabel(self.sweep_name)

        if save:
            self.save_fig_incremental_filename(fig,'plot_tail_hists_vs_sweep')
        
        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig

    def plot_tail_hist_integated(self, name='', save=True, log_plot=True, **kw):  
        ret = kw.get('ret', None)
        ax = kw.get('ax', None)

        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            save = False

        xx=self.tail_hist_b[:-1]
        yy=np.sum(self.tail_hist_h, axis=0)
        if log_plot:
            ax.semilogy(xx,yy,'-', color = 'k')
        else:
            ax.plot(xx,yy)
        #ax.colorbar()
        ax.set_xlabel('Time after sync [ns]')
        ax.set_ylabel('Counts')

        if save:
            self.save_fig_incremental_filename(fig,'plot_tail_hist_integated')
        
        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig

    def plot_tail_hist_all(self, name='', save=True, log_plot=True, offset=0, **kw):  
        ret = kw.get('ret', None)
        ax = kw.get('ax', None)

        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            save = False
        xx=self.tail_hist_b[:-1]
        indices = kw.pop('indices',range(self.sweep_length) )
        for i in indices:
            
            yy=self.tail_hist_h[i]+offset*i
            if log_plot:
                ax.semilogy(xx,yy,'-')
            else:
                ax.plot(xx,yy)
        #ax.colorbar()
        ax.set_xlabel('Time after sync [ns]')
        ax.set_ylabel('Counts')

        if save:
            self.save_fig_incremental_filename(fig,'plot_tail_hist_all')
        
        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig


def analyze_sweep(folder, name='', cr=False, roc=True):
    a = PQSequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name)
    if cr:
        a.get_cr_results(name)
    #get_tail_vs_sweep(self, syncs_per_sweep, channel, start, tail_length, pq_binsize=1e-3, hist_binsize=0.2):#ns
    a.get_tail_vs_sweep(250, 0, 0, 1000,pq_binsize=1)

    a.plot_tail_cts_vs_sweep()
    a.plot_tail_hists_vs_sweep()
    a.finish()