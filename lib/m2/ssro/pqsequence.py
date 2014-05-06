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

    def __init__(self, folder, **kw):
        sequence.SequenceAnalysis.__init__(self,folder, **kw)
        pq_folder=kw.pop('pq_folder', None)
        if pq_folder != None:
            h5filepath = toolbox.measurement_filename(pq_folder)
            h5mode=kw.get('hdf5_mode', 'r')
            self.pqf = h5py.File(h5filepath,h5mode)
        else:
            self.pqf=self.f

    def get_photons(self):
        """
        returns two filters (1d-arrays): whether events are ch0-photons/ch1-photons
        """
        channel = self.pqf['/PQ_channel-1'].value
        special = self.pqf['/PQ_special-1'].value

        is_not_special = special==0
        is_channel_0 = channel==0
        is_channel_1 = channel==1

        is_photon_0 = np.logical_and(is_not_special, is_channel_0)
        is_photon_1 = np.logical_and(is_not_special, is_channel_1)

        return is_photon_0, is_photon_1

    def get_sweep_idxs(self,  noof_syncs_per_sweep_pt=1):
        """
        Calculate the sweep-index for each PQ event, 
        based on the PQ sync number, sweep length and number of 
        sync pulses supplied per sweep point (i.e. per sequence start).
        """
        if not hasattr(self, 'sweep_pts'):
            self.sweep_pts = np.arange(len(self.ssro_results)) + 1
            self.sweep_name = 'sweep parameter'
        self.sync_nrs=self.pqf['/PQ_sync_number-1'].value  
        self.sweep_length = len(self.sweep_pts)
        self.syncs_per_sweep = noof_syncs_per_sweep_pt
        self.sweep_idxs=np.mod(np.floor((self.sync_nrs-1)/self.syncs_per_sweep),self.sweep_length)

class TailAnalysis(PQSequenceAnalysis):

    def get_tail_vs_sweep(self, channel, start_ns, tail_length_ns, pq_binsize_ns=1e-3, hist_binsize_ns=1.0, verbose= False):
        """
        Integrate & make a time histogram of all photon clicks on given channel,
        that arrive between start_ns and start_ns+tail_length_ns.
        Normalizes the integrated result to 1/(self.reps*syncs_per_sweep/self.sweep_length)
        """
        if not hasattr(self,'reps'):
            self.get_readout_results('ssro')
        if not hasattr(self, 'sweep_pts'):
            print 'get_sweep_pts first'
            return
        if not hasattr(self, 'sweep_idxs'):
            print 'get_sweep_idxs first'
            return

        is_ph = self.get_photons()[channel]
        sync_time_ns = self.pqf['/PQ_sync_time-1'].value * pq_binsize_ns

        hist_bins = np.arange(start_ns-hist_binsize_ns*.5,start_ns+1*tail_length_ns+hist_binsize_ns,hist_binsize_ns)
        
        self.tail_hist_h=np.zeros((self.sweep_length,len(hist_bins)-1))
        
        st_fltr = (start_ns  <= sync_time_ns) &  (sync_time_ns< (start_ns + tail_length_ns))
        if verbose:
            print 'total_photons in channel', channel, ':', len(sync_time_ns[np.where(is_ph)])  
            print 'total_photons in window:', len(sync_time_ns[np.where(is_ph & st_fltr)]) 
        valid_tail_idxs = self.sweep_idxs[np.where(is_ph & st_fltr)]
        if verbose:

            print 'total_sweeps in window:', len(valid_tail_idxs) 
            print 'total ph in window with sweep element 0:', len(np.where(valid_tail_idxs==0)[0])
            print 'div factor:', (self.reps*self.syncs_per_sweep/self.sweep_length)
        self.tail_cts_per_sweep_idx=np.zeros(self.sweep_length)
        
        for sweep_idx in range(self.sweep_length):
            self.tail_cts_per_sweep_idx[sweep_idx]= \
                                float(len(np.where(valid_tail_idxs==sweep_idx)[0])) / (self.reps*self.syncs_per_sweep/self.sweep_length)
            self.tail_hist_h[sweep_idx], self.tail_hist_b = \
                                np.histogram(sync_time_ns[np.where(is_ph & (self.sweep_idxs == sweep_idx))], bins=hist_bins)

    
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
        indices = kw.pop('indices',range(self.sweep_length))

        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            save = False
        xx=self.tail_hist_b[:-1]
        
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

class TailAnalysisIntegrated(TailAnalysis):  

    def get_sweep_idxs(self,  noof_syncs_per_sweep_pt=1):
        self.syncs_per_sweep = noof_syncs_per_sweep_pt

    def get_tail_vs_sweep(self, channel, start_ns, tail_length_ns, pq_binsize_ns=1e-3, verbose= False):
        """
        Integrate & make a time histogram of all photon clicks on given channel,
        that arrive between start_ns and start_ns+tail_length_ns.
        Normalizes the integrated result to 1/(self.reps*syncs_per_sweep/self.sweep_length)
        """
        if not hasattr(self,'reps'):
            self.get_readout_results('ssro')
        if not hasattr(self, 'sweep_pts'):
            print 'get_sweep_pts first'
            return
        if not hasattr(self, 'sweep_idxs'):
            print 'get_sweep_idxs first'
            return

        hist_bins = np.arange(start_ns-pq_binsize_ns*.5,start_ns+1*tail_length_ns+pq_binsize_ns,pq_binsize_ns)
        start=np.floor(start_ns/pq_binsize_ns)
        length=np.floor(tail_length_ns/pq_binsize_ns)

        self.tail_hist_h = (self.pqf['PQ_hist{}'.format(channel)].value)[start:start+length,:]
        self.tail_hist_b = np.linspace(start_ns, start_ns+tail_length_ns, len(self.tail_hist_h[:,0]))
        self.tail_cts_per_sweep_idx = np.sum(self.tail_hist_h, axis=0) / float(self.reps*self.syncs_per_sweep/self.sweep_length)
       

class FastSSROAnalysis(PQSequenceAnalysis):

    def get_sweep_pts(self):
        PQSequenceAnalysis.get_sweep_pts(self)
        self.sweep_pts=np.hstack((e,e) for e in self.sweep_pts)

    def get_fastssro_results(self,channel,pq_binsize_ns, hist_binsize_ns):
        if not hasattr(self,'reps'):
            self.get_readout_results('ssro')
        if not hasattr(self, 'sweep_pts'):
            print 'get_sweep_pts first'
            return
        if not hasattr(self, 'sweep_idxs'):
            print 'get_sweep_idxs first'
            return
        self.reps_per_sweep = float(self.reps)/len(self.sweep_pts)
        self.is_ph = self.get_photons()[channel]
        self.hist_binsize_ns = hist_binsize_ns
        self.sync_time_ns = self.pqf['/PQ_sync_time-1'].value * pq_binsize_ns
        self.extra_time_ns = (self.g.attrs['pq_sync_length']+self.g.attrs['wait_length'])*1e9

    def _get_relaxation(self, ms, sweep_index, start, length):
        bins = np.arange(start-self.hist_binsize_ns*.5,start+length,self.hist_binsize_ns)
        return np.histogram(self.sync_time_ns[np.where(self.is_ph & (self.sweep_idxs == 2*sweep_index+ms))], bins=bins)

    def _get_fidelity_and_mean_cpsh(self,ms,sweep_index, start, length):
        st_fltr = (start  <= self.sync_time_ns) &  (self.sync_time_ns< (start + length))
        vsync=self.sync_nrs[np.where(self.is_ph & st_fltr & (self.sweep_idxs == 2*sweep_index+ms))]
        cpsh = float(len(vsync))/self.reps_per_sweep
        if ms == 0:
            #print len(np.unique(vsync))
            return float(len(np.unique(vsync)))/self.reps_per_sweep,cpsh
        elif ms == 1:
            return float(self.reps_per_sweep-len(np.unique(vsync)))/self.reps_per_sweep,cpsh
        
    def _get_RO_window(self, ms, sweep_index):
        if ms == 0:
            start = start= self.g.attrs['A_SP_durations_AWG'][sweep_index]*1e9 + self.extra_time_ns
        elif ms ==1:
            start= self.g.attrs['E_SP_durations_AWG'][sweep_index]*1e9 + self.extra_time_ns
        length = self.g.attrs['E_RO_durations_AWG'][sweep_index]*1e9
        return start, length

    def _get_SP_window(self, ms, sweep_index):

        start = self.extra_time_ns
        if ms == 0:
            length = self.g.attrs['A_SP_durations_AWG'][sweep_index]*1e9
        elif ms ==1:
            length = self.g.attrs['E_SP_durations_AWG'][sweep_index]*1e9 
        return start, length

    def plot_relaxation(self,sweep_index,ms,st, save=True, **kw):
        ret = kw.get('ret', None)
        ax = kw.get('ax', None)
        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            save = False
        name= '_ms'+ str(ms) + '_sweep_' + str(sweep_index)
        start,length=getattr(self,'_get_'+st+'_window')(ms,sweep_index)
        y,x=self._get_relaxation(ms, sweep_index, start, length)
        x=x[:-1]
        y=y/float(self.reps_per_sweep)/(self.hist_binsize_ns*1e-9)
        print len(x), len(y)


        ax.plot(x,y)
        #ax.colorbar()
        ax.set_xlabel('Time [ns]')
        ax.set_ylabel('Countrate during {}, ms={}, [Hz]'.format(st,ms))

        if save:
            self.save_fig_incremental_filename(fig,st+'_relaxation_'+name)
        
        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig

    def plot_mean_fidelity(self, sweep_index, save=True, plot_points=100, **kw):
        ro_length = kw.pop('RO_length_ns', None)
        ret = kw.get('ret', None)
        ax = kw.get('ax', None)
        name= 'sweep_' + str(sweep_index)
        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            save = False
        start0, len0 = self._get_RO_window(0,sweep_index)
        if ro_length != None:
                len0=ro_length
        start1, len1 = self._get_RO_window(1,sweep_index)
        x=np.zeros(plot_points)
        f0=np.zeros(plot_points)
        f1=np.zeros(plot_points)
        mf=np.zeros(plot_points)

        for i,l in enumerate(np.linspace(0,len0,plot_points)):
            x[i]=l
            f0[i],_tmp = self._get_fidelity_and_mean_cpsh(0,sweep_index, start0, l)
            f1[i],_tmp = self._get_fidelity_and_mean_cpsh(1,sweep_index, start1, l)
            mf[i] = (f0[i]+f1[i])/2.
        ax.plot(x,f0, 'b-')
        ax.plot(x,f1, 'r-')
        ax.plot(x,mf, 'g-')
        #ax.colorbar()
        ax.set_xlabel('Time after RO start [ns]')
        ax.set_ylabel('Fidelity')

        if save:
            self.save_fig_incremental_filename(fig,'mean fidelity_'+name)
        
        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig
        
    def plot_fidelity_cpsh_vs_sweep(self, save=True, **kw):
        ret = kw.get('ret', None)
        ax = kw.get('ax', None)
        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            save = False
        ro_length = kw.pop('RO_length_ns', None)
        print self.sweep_pts
        f0=np.zeros(self.sweep_length/2)
        f1=np.zeros(self.sweep_length/2)
        mf=np.zeros(self.sweep_length/2)
        mcpsh0 = np.zeros(self.sweep_length/2)
        mcpsh1 = np.zeros(self.sweep_length/2)
        for i in range(self.sweep_length/2):
            #print i
            start0, length = self._get_RO_window(0,i)
            start1, _tmp = self._get_RO_window(1,i)
            if ro_length != None:
                length=ro_length
            f0[i],mcpsh0[i] = self._get_fidelity_and_mean_cpsh(0,i, start0, length)
            f1[i],mcpsh1[i] = self._get_fidelity_and_mean_cpsh(1,i, start1, length)
            mf[i] = (f0[i]+f1[i])/2.
            #etc
        x=self.sweep_pts[::2]
        print len(x), len(f0)
        #print x, f0
        ax.plot(x,f0, 'bo')
        ax.plot(x,f1, 'ro')
        ax.plot(x,mf, 'go')
        ax2=ax.twinx()
        ax2.plot(x, mcpsh0, 'b-')
        ax2.plot(x, mcpsh1, 'r-')
        ax.set_xlabel(self.sweep_name)
        ax.set_ylabel('Fidelity')
        ax2.set_ylabel('Ms=0 mean counts per shot')
        if save:
            self.save_fig_incremental_filename(fig,'fidelity_mcpsh_vs_sweep')
        
        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig
        
    def ssro_plots(self, sweep_index):
        self.plot_relaxation(sweep_index,ms=0, st='RO')
        self.plot_relaxation(sweep_index,ms=0, st='SP')
        self.plot_relaxation(sweep_index,ms=1, st='RO')
        self.plot_relaxation(sweep_index,ms=1, st='SP')
        self.plot_mean_fidelity(sweep_index)

class FastSSROAnalysisIntegrated(FastSSROAnalysis):

    def get_fastssro_results(self,channel,pq_binsize_ns):
        self.hist_binsize_ns = pq_binsize_ns
        self.hist = self.pqf['PQ_hist{}'.format(channel)].value
    
    def _get_relaxation(self, ms, sweep_index, start, length):
        start = np.floor(start / self.hist_binsize_ns)
        length = np.floor(length / self.hist_binsize_ns)
        return self.hist[start:start+length,2*sweep_index+ms]

    def _get_fidelity_and_mean_cpsh(self,ms,sweep_index, start, length):
        start = np.floor(start / self.hist_binsize_ns)
        length = np.floor(length / self.hist_binsize_ns)
        cpsh=np.sum(self.hist[start:start+length,2*sweep_index+ms])
        #fidelity not possible!

def analyze_tail(folder, name='ssro', cr=False, roc=True):
    a = TailAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name)
    if cr:
        a.get_cr_results(name)

    a.get_sweep_idxs(noof_syncs_per_sweep_pt=250)
    a.get_tail_vs_sweep(channel=0, start_ns=0, tail_length_ns=1000, pq_binsize_ns=1e-3, hist_binsize_ns=1.0, verbose= False)

    a.plot_tail_cts_vs_sweep()
    a.plot_tail_hists_vs_sweep()
    a.finish()



def fast_ssro_calib(folder, name='ssro',cr=True):
    a = FastSSROAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name)
    if cr:
        a.get_cr_results(name)

    a.get_sweep_idxs(noof_syncs_per_sweep_pt=1)
    a.get_fastssro_results(channel=0,pq_binsize_ns=1e-3, hist_binsize_ns=1.0)
    a.plot_fidelity_cpsh_vs_sweep(RO_length=None)
    a.ssro_plots(sweep_index = 0)
    a.finish()