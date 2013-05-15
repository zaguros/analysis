import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.spin import Nspin_twolevel_correction
from analysis.lib.m2 import m2
from analysis.lib.m2.ssro import ssro
from sequence_ssro import SequenceSSROAnalysis
from analysis.lib.math import error
from measurement.lib.tools import toolbox

class MBIAnalysis(m2.M2Analysis):
    def get_readout_results(self, name=''):
        self.result_corrected = False

        adwingrp = self.adwingrp(name)
        
        self.reps = self.g.attrs['reps_per_ROsequence']
        self.pts = adwingrp.attrs['sweep_length']
        self.readouts = adwingrp.attrs['nr_of_ROsequences']
        
        discards = len(adwingrp['ssro_results'].value) % (self.pts*self.readouts)
        self.reps = int(len(adwingrp['ssro_results'].value) / (self.pts*self.readouts))
        
        if discards > 0:
            results = adwingrp['ssro_results'].value.reshape(-1)[:-discards]
        else:
            results = adwingrp['ssro_results'].value

        self.ssro_results = results.reshape((-1,self.pts,self.readouts)).sum(axis=0)
        self.normalized_ssro = self.ssro_results/float(self.reps)
        self.u_normalized_ssro = \
            (self.normalized_ssro*(1.-self.normalized_ssro)/self.reps)**0.5

    def get_sweep_pts(self):
        self.sweep_name = self.g.attrs['sweep_name']
        self.sweep_pts = self.g.attrs['sweep_pts']

    def get_electron_ROC(self, ssro_calib_folder=toolbox.latest_data('SSRO')):
        self.p0 = np.zeros(self.normalized_ssro.shape)
        self.u_p0 = np.zeros(self.normalized_ssro.shape)
        
        ro_durations = self.g.attrs['E_RO_durations']
        roc = error.SingleQubitROC()

        for i in range(len(self.normalized_ssro[0])):
            roc.F0, roc.u_F0, roc.F1, roc.u_F1 = \
                ssro.get_SSRO_calibration(ssro_calib_folder, 
                        ro_durations[i])
            p0, u_p0 = roc.num_eval(self.normalized_ssro[:,i],
                    self.u_normalized_ssro[:,i])
            
            self.p0[:,i] = p0
            self.u_p0[:,i] = u_p0      
        
        self.result_corrected = True

    def get_N_ROC(self, F_init=1, u_F_init=0, F0_RO_pulse=1, u_F0_RO_pulse=0,
            F1_RO_pulse=1, u_F1_RO_pulse=0,
            ssro_calib_folder=toolbox.latest_data('SSRO')):

        self.p0 = np.zeros(self.normalized_ssro.shape)
        self.u_p0 = np.zeros(self.normalized_ssro.shape)
        
        ro_durations = self.g.attrs['E_RO_durations']
        roc = Nspin_twolevel_correction.NuclearSpinROC()
        roc.F_init = F_init
        roc.u_F_init = u_F_init
        roc.F0_RO_pulse = F0_RO_pulse
        roc.u_F0_RO_pulse = u_F0_RO_pulse
        roc.F1_RO_pulse = F1_RO_pulse
        roc.u_F1_RO_pulse = u_F1_RO_pulse
        
        for i in range(len(self.normalized_ssro[0])):
            roc.F0_ssro, roc.u_F0_ssro, roc.F1_ssro, roc.u_F1_ssro = \
                ssro.get_SSRO_calibration(ssro_calib_folder,
                        ro_durations[i])
            
            p0, u_p0 = roc.num_eval(self.normalized_ssro[:,i],
                    self.u_normalized_ssro[:,i])
          
            self.p0[:,i] = p0
            self.u_p0[:,i] = u_p0
        
        self.result_corrected = True

    def save(self, fname='analysis.hdf5'):
        f = h5py.File(os.path.join(self.folder, fname), 'w')
        g = f.create_group('readout_result')
        
        g.attrs['sweep_name'] = self.sweep_name

        f['/readout_result/sweep_pts'] = self.sweep_pts
        f['/readout_result/normalized_ssro'] = self.normalized_ssro
        f['/readout_result/u_normalized_ssro'] = self.u_normalized_ssro
        
        if self.result_corrected:
            f['readout_result/p0'] = self.p0
            f['readout_result/u_p0'] = self.u_p0

        f.close()
        

    def plot_results_vs_sweepparam(self, name='', save=True, **kw):
        labels = kw.get('labels', None)
        ret = kw.get('ret', None)
        ylim = kw.get('ylim', (-0.05, 1.05))
        
        if not hasattr(self, 'sweep_pts'):
            self.sweep_pts = np.arange(self.pts) + 1
            self.sweep_name = 'sweep parameter'
        
        fig = self.default_fig(figsize=(6,4))
        ax = self.default_ax(fig)
        
        if labels == None:
            labels = ['RO #%d' % i for i in range(self.readouts)]
        
        for i in range(self.readouts):
            if not self.result_corrected:
                ax.errorbar(self.sweep_pts, self.normalized_ssro[:,i], fmt='o-',
                    yerr=self.u_normalized_ssro[:,i], label=labels[i])
            else:
                ax.errorbar(self.sweep_pts, self.p0[:,i], fmt='o',
                    yerr=self.u_p0[:,i], label=labels[i])

        
        ax.set_xlabel(self.sweep_name)
        if not self.result_corrected:
            ax.set_ylabel('avg (uncorrected) outcome')
        else:
            ax.set_ylabel(r'$F(|0\rangle)$')
        
        ax.set_ylim(ylim)
        ax.legend(loc='best')

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

def analyze_single_sweep(folder, name='', correction='electron', **kw):
    F_init = kw.pop('F_init', 1.)
    u_F_init = kw.pop('u_F_init', 0.)
    F0_RO_pulse = kw.pop('F0_RO_pulse', 1.)
    u_F0_RO_pulse = kw.pop('u_F0_RO_pulse', 0.)
    F1_RO_pulse = kw.pop('F1_RO_pulse', 1.)
    u_F1_RO_pulse = kw.pop('u_F1_RO_pulse', 0.)

    a = MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name)
    
    if correction=='electron':
        a.get_electron_ROC()
    elif correction == 'N':
        a.get_N_ROC(F_init, u_F_init, F0_RO_pulse, u_F0_RO_pulse, F1_RO_pulse,
                u_F1_RO_pulse)

    a.save()    
    a.plot_results_vs_sweepparam(**kw)
    a.finish()


#### DEPRECATED
class PostInitDarkESRAnalysis(MBIAnalysis):
    
    sweep_name = 'MW frq. (MHz)'

    def get_sweep_pts(self, name=''):
        self.sweep_name = 'MW frq. (MHz) - %.3f GHz' % (self.g.attrs['mw_frq']*1e-9)
        self.sweep_pts = (self.g.attrs['RO_MW_pulse_ssbmod_frqs'])*1e-6

        return self.sweep_pts
    

def esr(folder=os.getcwd()):
    a = PostInitDarkESRAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results()
    a.plot_results_vs_sweepparam()

### ESR

class ConditionalPrecessionAnalysis(MBIAnalysis):

    sweep_name = 'time delay (us)'
    
    def get_sweep_pts(self, name=''):
        self.sweep_pts = self.g.attrs['AWG_wait_before_RO_durations'] / 1000.
        return self.sweep_pts


def condprec(folder=os.getcwd()):
    a = ConditionalPrecessionAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results()
    a.plot_results_vs_sweepparam()
        
### conditional Precession

class ElectronRabiAnalysis(MBIAnalysis):

    sweep_name = 'MW pulse duration (us)'

    def get_sweep_pts(self, name=''):
        self.sweep_pts = self.g.attrs['AWG_RO_MW_pulse_durations'] / 1000.
        return self.sweep_pts

def erabi(folder=os.getcwd()):
    a = ElectronRabiAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results()
    a.plot_results_vs_sweepparam()

### Electron Rabi
