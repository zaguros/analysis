import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting

from analysis.lib.m2 import m2
from sequence_ssro import SequenceSSROAnalysis


class MBIAnalysis(m2.M2Analysis):
    
    def get_readout_results(self, name=''):
        adwingrp = self.adwingrp(name)
        
        self.reps = self.g.attrs['reps_per_ROsequence']
        self.pts = adwingrp.attrs['sweep_length']
        self.readouts = adwingrp.attrs['nr_of_ROsequences']
        self.ssro_results = adwingrp['ssro_results'].value.reshape((-1,self.pts,self.readouts)).sum(axis=0)
        self.normalized_ssro = self.ssro_results/float(self.reps)
        self.u_normalized_ssro = \
            (self.normalized_ssro*(1.-self.normalized_ssro)/self.reps)**0.5

        return  self.normalized_ssro

    def plot_results_vs_sweepparam(self, name='', save=True, **kw):
        labels = kw.get('labels', None)
        
        if not hasattr(self, 'sweep_pts'):
            self.sweep_pts = np.arange(self.pts) + 1
            self.sweep_name = 'sweep parameter'
        
        fig = self.default_fig()
        ax = self.default_ax(fig)
        
        if labels == None:
            labels = ['RO #%d' % i for i in range(self.readouts)]
        
        for i in range(self.readouts):
            ax.errorbar(self.sweep_pts, self.normalized_ssro[:,i], fmt='o-',
                    yerr=self.u_normalized_ssro[:,i], label=labels[i])

        ax.set_xlabel(self.sweep_name)
        ax.set_ylabel('avg (uncorrected) outcome')
        ax.legend()

        if save:
            fig.savefig(
                os.path.join(self.folder, 'ssro_result_vs_sweepparam.pdf'),
                format='pdf')

    def finish(self):
        self.f.close()


class PostInitDarkESRAnalysis(MBIAnalysis):
    
    sweep_name = 'MW frq. (Hz)'

    def get_sweep_pts(self, name=''):
        self.sweep_pts = self.g.attrs['RO_MW_pulse_ssbmod_frqs'] + self.g.attrs['mw_frq']

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
