import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, readout_correction
from analysis.lib.math import error
from analysis.lib.m2 import m2
from measurement.lib.tools import toolbox

class SequenceSSROAnalysis(m2.M2Analysis):
    def get_readout_results(self, name=''):
        self.result_corrected = False

        adwingrp = self.adwingrp(name)
        
        self.reps = adwingrp['completed_reps'].value
        idx0 = adwingrp.attrs['RO_duration'] - 2
        self.ssro_results = adwingrp['ssro_results']\
                [idx0::adwingrp.attrs['RO_duration']]
        self.normalized_ssro = self.ssro_results/(float(self.reps)/len(self.sweep_pts))
        self.u_normalized_ssro = \
            (self.normalized_ssro*(1.-self.normalized_ssro)/self.reps)**0.5
        
        return self.normalized_ssro

    def get_sweep_pts(self):
        self.sweep_name = self.g.attrs['sweep_name']
        self.sweep_pts = self.g.attrs['sweep_pts']

    def get_electron_ROC(self, ssro_calib_folder=toolbox.latest_data('SSROCalibration')):
        self.p0 = np.zeros(self.normalized_ssro.shape)
        self.u_p0 = np.zeros(self.normalized_ssro.shape)
        
        ro_duration = self.g.attrs['RO_duration']
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

        if not hasattr(self, 'sweep_pts'):
            self.sweep_pts = np.arange(len(self.ssro_results)) + 1
            self.sweep_name = 'sweep parameter'

        fig = self.default_fig(figsize=(6,4))
        ax = self.default_ax(fig)
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
                os.path.join(self.folder, 'ssro_result_vs_sweepparam.pdf'),
                format='pdf')
        
        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig

    def finish(self):
        self.f.close()

def analyze_sweep(folder, name=''):
    a = SequenceSSROAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name)
    a.get_electron_ROC()
    a.plot_result_vs_sweepparam()
    a.finish()



class DarkESRAnalysis(SequenceSSROAnalysis):

    sweep_name = 'MW frq. (Hz)'
    
    def get_sweep_pts(self, name=''):
        self.sweep_pts = np.linspace(self.g.attrs['ssbmod_frq_start'], 
            self.g.attrs['ssbmod_frq_stop'], 
            self.g.attrs['pts'] ) + self.g.attrs['mw_frq']

        return self.sweep_pts

class ElectronRabiAnalysis(SequenceSSROAnalysis):

    sweep_name = 'MW pulse length (ns)'
    
    def get_sweep_pts(self, name=''):
        self.sweep_pts = self.g.attrs['AWG_RO_MW_pulse_durations']

        return self.sweep_pts


def darkesr(folder=os.getcwd()):
    a = DarkESRAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results()
    a.plot_result_vs_sweepparam()

def rabi(folder=os.getcwd()):
    a = ElectronRabiAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results()
    a.plot_result_vs_sweepparam()
