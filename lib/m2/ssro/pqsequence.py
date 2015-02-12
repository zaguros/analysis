import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, sequence
from analysis.lib.m2 import m2
from analysis.lib.pq import pq_tools
from analysis.lib.tools import toolbox

from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot

class PQSequenceAnalysis(sequence.SequenceAnalysis):

    def __init__(self, folder, **kw):
        sequence.SequenceAnalysis.__init__(self,folder, **kw)
        pq_folder=kw.pop('pq_folder', None)
        if pq_folder != None:
            if pq_folder =='bs_remote':
                h5filepath = 'X:' + self.g.attrs['bs_data_path'][12:]
                print h5filepath
            else:
                h5filepath = toolbox.measurement_filename(pq_folder)
            h5mode=kw.get('hdf5_mode', 'r')
            self.pqf = h5py.File(h5filepath,h5mode)
        else:
            self.pqf=self.f

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

    def plot_histogram(self,channel,start=None,length=None,fltr=None,hist_binsize=1,save=True, **kw):
        ret = kw.get('ret', None)
        ax = kw.get('ax', None)
        log_plot=kw.get('log_plot',True)
        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            save = False

        if start==None:
            start=self.g.attrs['MIN_SYNC_BIN']
        if length == None:
            stop = self.g.attrs['MAX_SYNC_BIN']
        else:
            stop = start + length

        is_ph = pq_tools.get_photons(self.pqf)[channel]
        if fltr == None:
            fltr=is_ph
        else:
            fltr=fltr & is_ph

        bins = np.arange(start-.5,stop,hist_binsize)
        y,x=np.histogram(self.pqf['/PQ_sync_time-1'].value[np.where(fltr)], bins=bins)
        x=x[:-1]
        print 'Total clicks:', np.sum(y)
        y=y/float(self.reps)
        
        if log_plot:
            ax.semilogy(x,y)
        else:
            ax.plot(x,y)
        #ax.colorbar()
        ax.set_xlabel('Time [bins]')
        ax.set_ylabel('Counts per rep per bin')

        if save:
            self.save_fig_incremental_filename(fig,'histogram_chan_{}'.format(channel))
        
        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig


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

        self.start_ns = start_ns

        is_ph = pq_tools.get_photons(self.pqf)[channel]
        sync_time_ns = self.pqf['/PQ_sync_time-1'].value * pq_binsize_ns

        hist_bins = np.arange(self.start_ns-hist_binsize_ns*.5,self.start_ns+1*tail_length_ns+hist_binsize_ns,hist_binsize_ns)
        
        self.tail_hist_h=np.zeros((self.sweep_length,len(hist_bins)-1))
        
        st_fltr = (self.start_ns  <= sync_time_ns) &  (sync_time_ns< (self.start_ns + tail_length_ns))
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
        ax.text(self.sweep_pts[len(self.sweep_pts)/2],np.min(self.tail_cts_per_sweep_idx)*1e4,'Tail start at {:.1f} ns'.format(self.start_ns))

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
        sync_nrs=self.pqf['/PQ_sync_number-1'].value
        if sync_nrs[-1] != self.reps:
            print 'WARNING last sync number ({}) != noof reps ({})! Sync error?'.format(sync_nrs[-1],self.reps)
        self.reps_per_sweep = float(self.reps)/len(self.sweep_pts)
        self.is_ph = pq_tools.get_photons(self.pqf)[channel]
        self.hist_binsize_ns = hist_binsize_ns
        self.sync_time_ns = self.pqf['/PQ_sync_time-1'].value * pq_binsize_ns
        self.extra_time_ns = (self.g.attrs['wait_length'])*1e9-100 #self.g.attrs['pq_sync_length']+

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
            #print 'len vsync',len(vsync)
            #print 'unique vsync',len(np.unique(vsync))
            return float(self.reps_per_sweep-len(np.unique(vsync)))/self.reps_per_sweep,cpsh
        
    def _get_RO_window(self, ms, sweep_index):
        if ms == 0:
            start = start= self.g.attrs['A_SP_durations_AWG'][sweep_index]*1e9 + 2*self.extra_time_ns
        elif ms ==1:
            start= self.g.attrs['E_SP_durations_AWG'][sweep_index]*1e9 + 2*self.extra_time_ns
        length = self.g.attrs['E_RO_durations_AWG'][sweep_index]*1e9
        return start, length

    def _get_SP_window(self, ms, sweep_index):

        start = self.extra_time_ns
        if ms == 0:
            length = self.g.attrs['A_SP_durations_AWG'][sweep_index]*1e9
        elif ms ==1:
            length = self.g.attrs['E_SP_durations_AWG'][sweep_index]*1e9 
        return start, length

    def plot_relaxation(self, sweep_index, ms, st, save=True, **kw):
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
        #print len(x), len(y)


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

    def plot_relaxation_vs_sweep(self, ms, st, save=True, **kw):
        ret = kw.get('ret', None)
        ax = kw.get('ax', None)
        ax2 = kw.get('ax2', None)
        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
            fig2 = self.default_fig(figsize=(10,6))
            #fig2.title('Relaxation fit vs' + self.sweep_name)
        else:
            save = False
        name= '_ms'+ str(ms) + '_' + str(st)
        sweep=self.sweep_pts[::2]
        taus=np.zeros(len(sweep))
        for i,sw in enumerate(sweep):
            #print i
            start,length=getattr(self,'_get_'+st+'_window')(ms,i)
            y,x=self._get_relaxation(ms, i, start, length)
            x=x[:-1]
            res = fit.fit1d(x,y, common.fit_exp_decay_shifted_with_offset, 0, y[0], 1000, x[0], 
            ret=True, do_print=False, fixed=[3])
    
            if res != False:
                ax2=plt.subplot(3,np.ceil(float(len(sweep)+1)/3),i+1)
                plot.plot_fit1d(res, x, ax=ax2, plot_data=True, print_info=False)
                ax2.text(ax2.get_xlim()[-1]*0.8,ax2.get_ylim()[-1]*0.8, '{:.2f}'.format(sw),horizontalalignment='right')
                if (res['params_dict']['A'] > 10*res['params_dict']['a']):
                    taus[i]=res['params_dict']['tau']/1e3
                else:
                    print 'Sweep point {} ignored'.format(i)
            else:
                print 'Could not fit sweep point {}'.format(i)
        skip_points = kw.get('skip_points', 0)
        xx=sweep[taus>1.][skip_points:]
        yy=taus[taus>1.][skip_points:]
        
        #y(x) = A * exp(-(x-x0)/tau) + a
        #g_a : offset
        #g_A : initial Amplitude
        #g_tau : decay constant
        #g_x0 : x offset
        res2 = fit.fit1d(xx,yy, common.fit_exp_decay_with_offset, 10, yy[0], xx[len(xx)/2]/2.,  
                ret=True, do_print=True, fixed=[])
        ax.plot(xx,yy, 'o')
        #ax.plot(sweep,10+yy[0]*np.exp(-(sweep- xx[0])/(xx[len(xx)/2]/2.)))
        
        if res2 != False:
            plot.plot_fit1d(res2, xx, ax=ax, plot_data=False, print_info=True)
        #ax.colorbar()
        ax.set_xlabel(self.sweep_name)
        ax.set_ylabel('ms={} {}-relaxation time [us]'.format(ms,st))
        return res

        if save:
            self.save_fig_incremental_filename(fig,name+'_relaxation_vs_sweep')
        
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
        u_f0=np.zeros(plot_points)
        u_f1=np.zeros(plot_points)
        u_mf=np.zeros(plot_points)
        #np.sqrt(pzero*(1-pzero)/reps)
        for i,l in enumerate(np.linspace(0,len0,plot_points)):
            x[i]=l
            f0[i],_tmp = self._get_fidelity_and_mean_cpsh(0,sweep_index, start0, l)
            f1[i],_tmp = self._get_fidelity_and_mean_cpsh(1,sweep_index, start1, l)
            u_f0[i] = np.sqrt(f0[i]*(1-f0[i])/self.reps_per_sweep)
            u_f1[i] = np.sqrt(f1[i]*(1-f1[i])/self.reps_per_sweep)
            mf[i] = (f0[i]+f1[i])/2.
            u_mf[i] = np.sqrt( (0.5*u_f0[i])**2 + (0.5*u_f1[i])**2 )

        ax.errorbar(x,f0, fmt='b-', yerr=u_f0)
        ax.errorbar(x,f1, fmt='g-', yerr=u_f1)
        ax.errorbar(x,mf, fmt='r-', yerr=u_mf)
        #ax.colorbar()
        ax.set_ylim(0.5,1.01)
        ax.set_xlabel('Time after RO start [ns]')
        ax.set_ylabel('Fidelity')

        ax.text(len0,0.86,'Fidelity at {} = {:.2f}, \n \
                    and RO time = {} ns: \n \
                    ms1 {:.2f} $\pm$ {:.2f} \n \
                    ms0 {:.2f} $\pm$ {:.2f} \n \
                    mean {:.2f} $\pm$ {:.2f} '.format(self.sweep_name,self.sweep_pts[sweep_index*2],
                                                      len0,
                                                      f1[-1]*100.,u_f1[-1]*100.,
                                                      f0[-1]*100.,u_f0[-1]*100.,
                                                      mf[-1]*100.,u_mf[-1]*100.),
                        horizontalalignment='right',verticalalignment='top')
        ii = np.argmax(mf)
        #print 'Max mean fidelity at RO time = {:.1f}: ms0 {:.2f} $\pm$ {:.2f}, \
        #                                              ms1 {:.2f} $\pm$ {:.2f}, \
        #                                              mean {:.2f} $\pm$ {:.2f}'.format(x[ii],
        #                                                                               f0[ii]*100.,u_f0[ii]*100.,
        #                                                                               f1[ii]*100.,u_f1[ii]*100.,
        #                                                                               mf[ii]*100.,u_mf[ii]*100.)

        if save:
            dat0 = np.vstack((x, f0, u_f0)).T
            dat1 = np.vstack((x, f1, u_f1)).T
            datm = np.vstack((x, mf, u_mf)).T
            f = self.analysis_h5data()
            dat_name='fidelity-{}'.format(sweep_index)
            if not dat_name in f:
                f.create_group(dat_name)
            g = f['/'+dat_name]
            if 'ms0' in g:
                del  g['ms0']
            if 'ms1' in g:
                del  g['ms1']
            if 'mean' in g:
                del  g['mean']
            g['ms0'] = dat0
            g['ms1'] = dat1
            g['mean'] = datm
            g.attrs['max RO time']= len0
            g.attrs['sweep_value']= self.sweep_pts[sweep_index*2]
            g.attrs['sweep_name'] = self.sweep_name
            f.close()

            self.save_fig_incremental_filename(fig,'mean fidelity-{}'.format(sweep_index)+name)

        
        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig
        if ret == 'fid':
            return x,f0,f1,mf
        
    def plot_fidelity_cpsh_vs_sweep(self, save=True, **kw):
        ret = kw.get('ret', None)
        ax = kw.get('ax', None)
        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            save = False
        ro_length = kw.pop('RO_length_ns', None)
        #print self.sweep_pts
        f0=np.zeros(self.sweep_length/2)
        f1=np.zeros(self.sweep_length/2)
        mf=np.zeros(self.sweep_length/2)
        u_f0=np.zeros(self.sweep_length/2)
        u_f1=np.zeros(self.sweep_length/2)
        u_mf=np.zeros(self.sweep_length/2)
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
            u_f0[i] = np.sqrt(f0[i]*(1-f0[i])/self.reps_per_sweep)
            u_f1[i] = np.sqrt(f1[i]*(1-f1[i])/self.reps_per_sweep)
            mf[i] = (f0[i]+f1[i])/2.
            u_mf[i] = np.sqrt( (0.5*u_f0[i])**2 + (0.5*u_f1[i])**2 )
            #etc
        x=self.sweep_pts[::2]
        #print len(x), len(f0)
        #print x, f0 ax.errorbar(time, fid0, fmt='.', yerr=fid0_err, label='ms=0')
        ax.errorbar(x,f0, fmt='bo', yerr=u_f0)
        ax.errorbar(x,f1, fmt='go', yerr = u_f1)
        ax.errorbar(x,mf, fmt='ro',yerr = u_mf)
        ax.set_ylim(0.5,1.01)
        ii = np.argmax(mf)
        ax.text(x[-1], 0.86, 'Max fid. at sweep pt. = {:.1f}: \n \
                                    ms1 {:.2f} $\pm$ {:.2f},\n \
                                    ms0 {:.2f} $\pm$ {:.2f},\n \
                                    mean {:.2f} $\pm$ {:.2f}'.format(x[ii],
                                                                     f1[ii]*100.,u_f1[ii]*100.,
                                                                     f0[ii]*100.,u_f0[ii]*100.,
                                                                     mf[ii]*100.,u_mf[ii]*100.),
                horizontalalignment='right',verticalalignment='top')
        ax2=ax.twinx()
        ax2.plot(x, mcpsh0, 'b-')
        ax2.plot(x, mcpsh1, 'g-')
        ax.set_xlabel(self.sweep_name)
        ax.set_ylabel('Fidelity')
        ax2.set_ylabel('Ms=0 mean counts per shot')
        if save:
            self.save_fig_incremental_filename(fig,'fidelity_mcpsh_vs_sweep')
        
        if ret == 'ax':
            return ax
        elif ret == 'fig':
            return fig
        elif ret =='max_sweep_index':
            return ii
        
    def ssro_plots(self, sweep_index, **kw):
        self.plot_relaxation(sweep_index,ms=0, st='RO', ret='fig', **kw)
        plt.close()
        self.plot_relaxation(sweep_index,ms=0, st='SP', ret='fig', **kw)
        plt.close()
        self.plot_relaxation(sweep_index,ms=1, st='RO', ret='fig', **kw)
        plt.close()
        self.plot_relaxation(sweep_index,ms=1, st='SP', ret='fig', **kw)
        plt.close()
        self.plot_mean_fidelity(sweep_index, **kw)

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
        plt.close('all')

    a.get_sweep_idxs(noof_syncs_per_sweep_pt=250)
    a.get_tail_vs_sweep(channel=0, start_ns=0, tail_length_ns=1000, pq_binsize_ns=1e-3, hist_binsize_ns=1.0, verbose= False)

    a.plot_tail_cts_vs_sweep()
    a.plot_tail_hists_vs_sweep()
    a.finish()



def fast_ssro_calib(folder='', name='ssro',cr=True, RO_length_ns=3500, plot_sweep_index='max'):
    if folder=='':
        folder=toolbox.latest_data('FastSSRO')
    a = FastSSROAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name)
    if cr:
        a.get_cr_results(name)
        plt.close('all')

    a.get_sweep_idxs(noof_syncs_per_sweep_pt=1)
    a.get_fastssro_results(channel=0,pq_binsize_ns=1.0, hist_binsize_ns=100)
    swp_idx=a.plot_fidelity_cpsh_vs_sweep(RO_length_ns=RO_length_ns, ret = 'max_sweep_index')
    if plot_sweep_index != 'max':
        swp_idx = plot_sweep_index
    a.ssro_plots(sweep_index = swp_idx, RO_length_ns=RO_length_ns)
    a.finish()

def get_analysed_fast_ssro_calibration(folder, readout_time=None, sweep_index=None):
    fp = os.path.join(folder, 'analysis.hdf5')
    f = h5py.File(fp, 'r')

    key=None
    for k in f.keys():
        if 'fidelity' in k:
            if sweep_index==None:
                key=k
            elif k=='fidelity-{}'.format(sweep_index):
                key=k
    if key == None:
        print 'No analysis found, correct sweep-index specified?'
    g=f[key]

    times = g['mean'].value[:,0]
    fids0 = g['ms0'].value
    fids1 = g['ms1'].value

    if readout_time==None:
        tidx=len(times)-1
    else:
        tidx = np.argmin(abs(times-readout_time))

    f0 = fids0[tidx,1]
    u_f0 = fids0[tidx,2]
    f1 = fids1[tidx,1]
    u_f1 = fids1[tidx,2]

    return f0, u_f0, f1, u_f1