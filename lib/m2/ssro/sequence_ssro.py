import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting

from analysis.lib.m2 import m2
from measurement.lib.tools import toolbox

class SequenceSSROAnalysis(m2.M2Analysis):
    
    def get_readout_results(self, name=''):
        adwingrp = self.adwingrp(name)
        
        self.reps = adwingrp['completed_reps'].value
        idx0 = adwingrp.attrs['RO_duration'] - 2
        self.ssro_results = adwingrp['ssro_results']\
                [idx0::adwingrp.attrs['RO_duration']]
        self.normalized_ssro = self.ssro_results/(float(self.reps)/len(self.sweep_pts))
        
        return self.normalized_ssro
    
    
    def plot_result_vs_sweepparam(self, name='', save=True):

        if not hasattr(self, 'sweep_pts'):
            self.sweep_pts = np.arange(len(self.ssro_results)) + 1
            self.sweep_name = 'sweep parameter'

        fig = self.default_fig()
        ax = self.default_ax(fig)
        ax.plot(self.sweep_pts, self.normalized_ssro, 'o-')
        ax.set_xlabel(self.sweep_name)
        ax.set_ylabel('avg (uncorrected) outcome')

        if save:
            fig.savefig(
                os.path.join(self.folder, 'ssro_result_vs_sweepparam.pdf'),
                format='pdf')

    def finish(self):
        self.f.close()

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
