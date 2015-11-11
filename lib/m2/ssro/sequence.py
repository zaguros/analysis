import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro
from analysis.lib.m2.ssro import ssro as ssro_calib
from analysis.lib.math import error
from analysis.lib.m2 import m2
from analysis.lib.tools import toolbox

class SequenceAnalysis(m2.M2Analysis):
    
    def get_readout_results(self, name=''):
        self.result_corrected = False

        adwingrp = self.adwingrp(name)        
        self.reps = adwingrp['completed_reps'].value
        self.ssro_results = adwingrp['RO_data'].value
        self.normalized_ssro = self.ssro_results/(float(self.reps)/len(self.sweep_pts))
        self.u_normalized_ssro = \
            (self.normalized_ssro*(1.-self.normalized_ssro)/(float(self.reps)/len(self.sweep_pts)))**0.5  #this is quite ugly, maybe replace?
        
        return self.normalized_ssro
      

    def get_cr_results(self, name='', plot=True):
        adwingrp = self.adwingrp(name)   
        self.cr_before=adwingrp['CR_before'].value
        self.cr_after=adwingrp['CR_after'].value

        if plot:
            self._plot_cr('CR before '+name, self.cr_before)
            self._plot_cr('CR after '+name, self.cr_after)
            
    def _plot_cr(self,name,cr_counts):
        title_suffix = ': '+name if name != '' else ''
        fn_suffix = '_'+name if name != '' else ''

        fig = self.default_fig(figsize=(6,4))
        ax = self.default_ax(fig)
        ax.hist(cr_counts, abs(max(cr_counts)-min(cr_counts)+1), 
                normed=True)
        ax.set_xlabel('counts during CR check sequence')
        ax.set_ylabel('probability')
        ax.set_title(self.default_plot_title + title_suffix)
        
        fig.savefig(os.path.join(self.folder, 
            'cr_check_before'+fn_suffix+'.'+self.plot_format), 
            format=self.plot_format)

    def fix_sweep_pts(self, new_sweep_name, new_sweep_pts_attr):
        self.g.attrs['sweep_name']=new_sweep_name
        self.g.attrs['sweep_pts']=self.g.attrs[new_sweep_pts_attr]
        
    def get_sweep_pts(self):
        self.sweep_name = self.g.attrs['sweep_name']
        self.sweep_pts = self.g.attrs['sweep_pts']
    
    def get_mean_cr_cts(self, save=True, pts=1, max_cr=-1):
        ### plot the mean of the CR counts --- without the zero --- vs the sweep-param
       
        if max_cr < 0:
            max_cr = np.max(self.cr_after)

        self.sweep_CR_hist = np.zeros((pts, max_cr+1))
        self.sweep_CR_sum = np.zeros(pts)
        self.sweep_CR_variance = np.zeros(pts)
        
        for i in range(pts):
            cr = self.cr_after[i::pts]
            self.sweep_CR_hist[i,:], binedges = np.histogram(cr, 
                bins=np.arange(max_cr+2)-0.5,
                normed=True)
            self.sweep_CR_sum[i] = float(np.sum(cr))/len(np.where(cr>0)[0])
            self.sweep_CR_variance[i] = np.sqrt(np.sum((cr[np.where(cr>0)[0]] - self.sweep_CR_sum[i])**2)/len(np.where(cr>0)[0]))
       
        return (self.sweep_CR_hist, self.sweep_CR_sum, self.sweep_CR_variance)

    def plot_cr_vs_sweep(self, save=True, max_cr=-1, ionization_crit=1):
        
        pts=len(self.sweep_pts)
        
        if max_cr < 0:
            max_cr = np.max(self.cr_after)

        self.sweep_CR_hist, sweep_CR_sum, self.sweep_CR_variance = self.get_mean_cr_cts(save, pts, max_cr)

        ### plot the mean of the CR counts --- without the zero --- vs the sweep-param
       
        fig = self.default_fig(figsize=(6,4))
        ax = self.default_ax(fig)
        ax.plot(self.sweep_pts, self.sweep_CR_sum,'o')
        ax.set_xlabel(self.sweep_name)
        ax.set_ylabel('mean CR counts after RO')
        if save:
            fig.savefig(
                os.path.join(self.folder, 'post-CR_sum_vs_sweepparam.png'),
                format='png')
        
        ### plot the height of the zero CR counts bar vs the sweep-param

        zero_bar = self.sweep_CR_hist[:,0:ionization_crit].sum(axis=-1).reshape(-1)
        
        fig = self.default_fig(figsize=(6,4))
        ax = self.default_ax(fig)
        ax.plot(self.sweep_pts,zero_bar,'o')
        ax.set_xlabel(self.sweep_name)
        ax.set_ylabel('Ionization probability')
        if save:
            fig.savefig(
                os.path.join(self.folder, 'post-CR_zero_vs_sweepparam.png'),
                format='png')

        ### color plot of sweep-resolved CR counts

        fig = self.default_fig(figsize=(6,4))
        ax = self.default_ax(fig)

        # xaxis = 
        im = ax.imshow(self.sweep_CR_hist.transpose(), 
            origin = 'lower',
            aspect = 'auto',
            interpolation='nearest',
            cmap = 'gist_heat')

        xt = ax.get_xticks().astype(int)
        yt = ax.get_yticks()

        xt = xt[xt>=0]
        xt = xt[xt<pts]
        ax.set_xticks(xt)
        ax.set_xticklabels(self.sweep_pts[xt])

        ax.set_ylabel('CR counts')
        ax.set_xlabel(self.sweep_name)

        cbar = fig.colorbar(im)
        cbar.set_label('Occurrence')

        # fig = self.default_fig(figsize=(6,4))
        # ax = self.default_ax(fig)
        # off=0
        # offstep=1.2*max([np.max(h[0]) for h in sweep_CR_hist])
        # for i,x in enumerate(self.sweep_pts):
        #     xx=sweep_CR_hist[i][1][:-1]
        #     yy=sweep_CR_hist[i][0]
        #     ax.bar(xx,yy, bottom=off)
        #     ax.text(max_cr,off,self.sweep_name+':'+str(x),horizontalalignment='right')
        #     off=off+offstep
        
        if save:
            fig.savefig(
                os.path.join(self.folder, 'post-CR_vs_sweepparam.png'),
                format='png')
        
    def get_electron_ROC(self, **kw):
        ssro_calib_folder = kw.pop('ssro_calib_folder', toolbox.latest_data('SSROCalibration'))
        print ssro_calib_folder
        if ssro_calib_folder == '':
                ssro_calib_folder = toolbox.latest_data('SSROCalibration')
        self.p0 = np.zeros(self.normalized_ssro.shape)
        self.u_p0 = np.zeros(self.normalized_ssro.shape)
        
        ro_duration = self.g.attrs['SSRO_duration']
        roc = error.SingleQubitROC()
        roc.F0, roc.u_F0, roc.F1, roc.u_F1 = \
            ssro.get_SSRO_calibration(ssro_calib_folder, 
                    ro_duration)
        p0, u_p0 = roc.num_eval(self.normalized_ssro,
                self.u_normalized_ssro)
            
        self.p0 = p0
        self.u_p0 = u_p0
        
        self.result_corrected = True
    
    def plot_result_vs_sweepparam(self, name='', save=True, **kw):
        ret = kw.get('ret', None)
        ax = kw.get('ax', None)
        ylim = kw.get('ylim', (-0.05, 1.05))
       
        if not hasattr(self, 'sweep_pts'):
            self.sweep_pts = np.arange(len(self.ssro_results)) + 1
            self.sweep_name = 'sweep parameter'

        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            save = False
        
        if not self.result_corrected:
            ax.errorbar(self.sweep_pts, self.normalized_ssro, fmt='o-',
                yerr=self.u_normalized_ssro)
        else:
            ax.errorbar(self.sweep_pts, self.p0, fmt='o',
                yerr=self.u_p0)
    
        ax.set_xlabel(self.sweep_name)

        if not self.result_corrected:
            ax.set_ylabel('avg (uncorrected) outcome')
        else:
            ax.set_ylabel(r'$F(|0\rangle)$')

        #if self.result_corrected:
        #ax.set_ylim(-0.05, 1.05)
        if save:
            fig.savefig(
                os.path.join(self.folder, 'ssro_result_vs_sweepparam.png'),
                format='png')
        
        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig

    def finish(self):
        self.f.close()

class MagnetometrySequenceAnalysis(SequenceAnalysis):

    def get_magnetometry_phase_calibration(self, name='',ssro_calib_folder=''):

        self.result_corrected = False
        adwingrp = self.adwingrp(name)        
        self.reps = adwingrp['completed_reps'].value-1
        RO_clicks = adwingrp['RO_data'].value
        phase = adwingrp['set_phase'].value
        times=self.g.attrs['ramsey_time']

        phase_values = np.unique(phase)
        self.ssro_results = np.zeros(len(phase_values))

        
        for j in phase_values:
            self.ssro_results [np.min(np.where(phase==j))] = np.sum(RO_clicks[np.where(phase==j)])
            
        if len(self.ssro_results)==1:
            self.sweep_pts = self.g.attrs['sweep_pts']
            phase=np.array(list(np.arange(len(self.sweep_pts)))*int((self.reps+1)/float(len(self.sweep_pts))))
       
            phase_values = np.unique(phase)
            self.ssro_results = np.zeros(len(phase_values))
            
            for j in phase_values:

                self.ssro_results [np.min(np.where(phase==j))] = np.sum(RO_clicks[np.where(phase==j)])

        #self.sweep_pts = times
        #print 'len phases',len(phases)
        #print phases
        #print phase_values

        self.normalized_ssro = self.ssro_results*len(self.ssro_results)/(float(self.reps))
        self.u_normalized_ssro = (self.normalized_ssro*(1.-self.normalized_ssro)/(float(self.reps)))**0.5  #this is quite ugly, maybe replace?

        
        if ssro_calib_folder == '':
                ssro_calib_folder = toolbox.latest_data('SSROCalibration')
                print ssro_calib_folder
        self.p0 = self.normalized_ssro
        self.u_p0 = self.u_normalized_ssro
        
        ro_duration = self.g.attrs['SSRO_duration']
        roc = error.SingleQubitROC()
        roc.F0, roc.u_F0, roc.F1, roc.u_F1 = \
            ssro.get_SSRO_calibration(ssro_calib_folder, 
                    ro_duration)
        p0, u_p0 = roc.num_eval(self.normalized_ssro,
                self.u_normalized_ssro)
            
        self.p0 = p0
        self.u_p0 = u_p0
        
        self.result_corrected = True
        return self.p0,self.u_p0

    def get_magnetometry_data (self, name='', ssro = True):
        self.result_corrected = False

        adwingrp = self.adwingrp(name)        
        self.reps = adwingrp['completed_reps'].value
        RO_clicks = np.array(adwingrp['RO_data'].value)
        CR_after = np.array(adwingrp['CR_after'].value)
        set_phase = adwingrp['set_phase'].value
        theta = adwingrp['theta'].value
        #print theta
        self.phase_calibration_msmnt = False
        self.debug_pk = True

        try:
            self.theta_opt = np.array(adwingrp['theta'].value)
            self.p_tn = np.array(adwingrp['real_p_tn'].value) +1j*np.array(adwingrp['imag_p_tn'].value)
            self.p_2tn = np.array(adwingrp['real_p_2tn'].value) +1j*np.array(adwingrp['imag_p_2tn'].value)
        except:
            self.p_tn  = None
            self.p_2tn = None
            self.debug_pk = False
        self.repetitions=self.g.attrs['repetitions']
        self.sweep_pts = self.g.attrs['sweep_pts']
        self.ramsey_time = self.g.attrs['ramsey_time'] 
        self.N = self.g.attrs['adptv_steps']
        try:
            self.do_add_phase=self.g.attrs['do_add_phase']
        except:
            self.do_add_phase=0   

        if self.do_add_phase==0:
            self.phase_update=False
        else:
            self.phase_update=True     
        try:
            self.K = self.g.attrs['adptv_steps']-1
        except:
            self.K = self.N-1
        try:
            self.phase_update=self.g.attrs['phase_update']
        except:    
            self.phase_update=False
        try:
            self.phases_detuning = self.g.attrs['phases_detuning'] 
            self.set_detuning = self.g.attrs['set_detuning_value']
            self.M = self.g.attrs['M']  
            self.t0 = self.g.attrs['tau0']
            self.maj_reps = self.g.attrs['reps_majority_vote']
            self.maj_thr = self.g.attrs['threshold_majority_vote']
        except:
            self.phase_calibration_msmnt = True        

        try:
            self.test_adwin = self.g.attrs['do_ext_test'] 
            self.G = self.g.attrs['G'] 
            self.F = self.g.attrs['F'] 
            self.renorm_ssro = self.g.attrs['renorm_ssro_M'] 
            self.exp_fid0 = self.g.attrs['fid0'] 
            self.exp_fid1 = self.g.attrs['fid1']
            self.msmnt_type = 'realtime' 
        except:
            self.msmnt_type = 'table_based'
        try:
            self.T2_mult_t0 = self.g.attrs['T2']
        except:
            self.T2_mult_t0 = 46    
        #print self.T2_mult_t0

        try:
            timer_data = np.array(adwingrp['timer'].value)

            timer_data = timer_data[timer_data != 0]
            self.timer = np.mean (timer_data)*3.33333*1e-9*1000
        except:
            self.timer = None
        try:
            self.save_pk_n = self.g.attrs['save_pk_n']
            self.save_pk_m = self.g.attrs['save_pk_m']
        except:    

            self.save_pk_n = 0
            self.save_pk_m = 0

        try:
            self.swarm_opt = self.g.attrs['swarm_opt']
        except:
            self.swarm_opt= 0

        if self.save_pk_n>0:
            self.real_pk_adwin = np.array(adwingrp['real_p_k'].value)
            self.imag_pk_adwin = np.array(adwingrp['imag_p_k'].value)
        else:
            self.real_pk_adwin = None
            self.imag_pk_adwin = None

        n_points = len(self.sweep_pts)
        rows = int(self.reps/float(len(self.sweep_pts)))
        cols = len(self.sweep_pts)
        #print 'self reps',self.reps

        if self.do_add_phase==1:
            m_sweeps=self.N*self.G + self.N*self.F*(self.N-1)/2
            rows_addphase = int(len(RO_clicks)/float(m_sweeps))
            self.reps=rows_addphase
            cols_addphase = m_sweeps
            #print 'rows: ',rows_addphase, 'cols: ', cols_addphase, 'm_sweeps: ', m_sweeps, 'reps:', self.reps

            RO_clicks = RO_clicks[:rows_addphase*cols_addphase]
            self.clicks = np.squeeze(np.reshape(RO_clicks, (rows_addphase, cols_addphase)))
            theta = theta[:rows*cols]
            self.theta = np.squeeze(np.reshape(theta, (rows, cols)))

        elif self.swarm_opt:
            m_sweeps=self.N*self.G + self.N*self.F*(self.N-1)/2
            rows_addphase = int(len(RO_clicks)/float(m_sweeps))
            self.reps=rows_addphase
            cols_addphase = m_sweeps
            RO_clicks = RO_clicks[:rows_addphase*cols_addphase]
            self.clicks = np.squeeze(np.reshape(RO_clicks, (rows_addphase, cols_addphase)))
            theta = theta[:rows_addphase*cols_addphase]
            self.theta = np.squeeze(np.reshape(theta, (rows_addphase, cols_addphase)))
        else:
            RO_clicks = RO_clicks[:rows*cols]
            self.clicks = np.squeeze(np.reshape(RO_clicks, (rows, cols)))
            self.reps=rows
            theta = theta[:rows*cols]
            self.theta = np.squeeze(np.reshape(theta, (rows, cols)))

        set_phase = set_phase[:rows*cols]
        #'14-12-18 Only for calibration. put back again if adaptive magnetometry
        CR_after = CR_after[:rows*cols]

        
        self.set_phase = np.squeeze(np.reshape(set_phase, (rows, cols)))
        #'14-12-18 Only for calibration. put back again if adaptive magnetometry
        self.CR_after =  np.squeeze(np.reshape(CR_after, (rows, cols)))

        if ssro:
            self.normalized_ssro = np.sum(self.clicks, axis=0)/(float(self.reps/n_points))
            self.u_normalized_ssro = (self.normalized_ssro*(1.-self.normalized_ssro)/(float(self.reps/n_points)))**0.5  






def analyze_sweep(folder, name='', cr=False, roc=True):
    a = SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name)
    if cr:
        a.get_cr_results(name)
    if roc:
        a.get_electron_ROC()
    a.plot_result_vs_sweepparam()
    a.finish()
    
def analyze_cr_only(folder,name=''):
    a = SequenceAnalysis(folder)
    a.get_cr_results(name)
    a.finish()
   
def analyze_cr_sweep(folder,name=''):
    a = SequenceAnalysis(folder)
    a.get_cr_results(name)
    a.get_sweep_pts()
    a.plot_cr_vs_sweep()
    a.finish()