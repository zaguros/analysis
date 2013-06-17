import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro
from analysis.lib.math import error
from analysis.lib.m2 import m2
from measurement.lib.tools import toolbox

class SequenceAnalysis(m2.M2Analysis):
    
    def get_readout_results(self, name=''):
        self.result_corrected = False

        adwingrp = self.adwingrp(name)        
        self.reps = adwingrp['completed_reps'].value
        self.ssro_results = adwingrp['RO_data'].value
        self.normalized_ssro = self.ssro_results/(float(self.reps)/len(self.sweep_pts))
        self.u_normalized_ssro = \
            (self.normalized_ssro*(1.-self.normalized_ssro)/self.reps)**0.5
        
        return self.normalized_ssro
    
    def get_cr_resutls(self, name='', plot=True):
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
        
    def get_sweep_pts(self):
        self.sweep_name = self.g.attrs['sweep_name']
        self.sweep_pts = self.g.attrs['sweep_pts']
    
    def plot_cr_vs_sweep(self, save=True):
        pts=len(self.sweep_pts)
        sweep_CR_hist=[]
        max_cr=0
        self.sweep_CR_sum=np.zeros(pts)
        
        for i in range(pts):
            cr=self.cr_after[i::pts]
            max_cr=max(np.max(cr),max_cr)
            sweep_CR_hist.append(np.histogram(cr,bins=np.max(cr),normed=True))
            self.sweep_CR_sum[i]=float(np.sum(cr))/len(np.where(cr>0)[0])
        fig = self.default_fig(figsize=(6,4))
        ax = self.default_ax(fig)
        ax.plot(self.sweep_pts,self.sweep_CR_sum,'o')
        ax.set_xlabel(self.sweep_name)
        if save:
            fig.savefig(
                os.path.join(self.folder, 'cr_after_sum_vs_sweepparam.png'),
                format='png')
        
        zero_bar=np.array([h[0][0] for h in sweep_CR_hist])
        fig = self.default_fig(figsize=(6,4))
        ax = self.default_ax(fig)
        ax.plot(self.sweep_pts,zero_bar,'o')
        ax.set_xlabel(self.sweep_name)
        if save:
            fig.savefig(
                os.path.join(self.folder, 'cr_after_0_vs_sweepparam.png'),
                format='png')
                
        fig = self.default_fig(figsize=(6,4))
        ax = self.default_ax(fig)
        off=0
        offstep=1.2*max([np.max(h[0]) for h in sweep_CR_hist])
        for i,x in enumerate(self.sweep_pts):
            xx=sweep_CR_hist[i][1][:-1]
            yy=sweep_CR_hist[i][0]
            ax.bar(xx,yy, bottom=off)
            ax.text(max_cr,off,self.sweep_name+':'+str(x),horizontalalignment='right')
            off=off+offstep
        if save:
            fig.savefig(
                os.path.join(self.folder, 'cr_after_vs_sweepparam.png'),
                format='png')
        
    def get_electron_ROC(self, **kw):
        ssro_calib_folder = kw.pop('ssro_calib_folder', toolbox.latest_data('SSROCalibration'))
                
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

        ax.set_ylim(-0.05, 1.05)

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

def analyze_sweep(folder, name=''):
    a = SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name)
    a.get_cr_resutls(name)
    a.get_electron_ROC()
    a.plot_result_vs_sweepparam()
    a.finish()
    
def analyze_cr_only(folder,name=''):
    a = SequenceAnalysis(folder)
    a.get_cr_resutls(name)
    a.finish()
    
def analyze_cr_sweep(folder,name=''):
    a = SequenceAnalysis(folder)
    a.get_cr_resutls(name)
    a.get_sweep_pts()
    a.plot_cr_vs_sweep()
    a.finish()