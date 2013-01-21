import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting

from measurement.lib.tools import toolbox

# TODO need also a function that simply does all
class SequenceSSROAnalysis:

    def __init__(self, folder):
        self.folder = folder
        self.h5filepath = toolbox.measurement_filename(folder)
        self.f = h5py.File(self.h5filepath, 'r')
        self.name = self.f.keys()[0]
        self.g = self.f[self.name]

        self.measurementstring = os.path.split(self.folder)[1]
        self.timestamp = os.path.split(os.path.split(self.folder)[0])[1] \
                + '/' + self.measurementstring[:6]
        self.measurementstring = self.measurementstring[7:]
        

    def adwingrp(self, name=''):
        
        if name != '':
            adwingrpname = name
        else:
            if len(self.g.keys()) == 1:
                adwingrpname = self.g.keys()[0]
            else:
                logging.error("More than one measurement. Please give a name")
                return False

        return self.g[adwingrpname]
    
    def get_ssro_results(self, name=''):
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

        self.sweepresult_fig = plt.figure()
        ax = self.sweepresult_fig.add_subplot(111)
        ax.plot(self.sweep_pts, self.normalized_ssro, 'ko-')
        ax.set_xlabel(self.sweep_name)
        ax.set_ylabel('avg (uncorrected) outcome')
        ax.set_title(self.timestamp+'\n'+self.measurementstring)

        if save:
            self.sweepresult_fig.savefig(
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
        self.sweep_pts = np.linspace(self.g.attrs['pulse_length_start'], 
            self.g.attrs['pulse_length_stop'], 
            self.g.attrs['pts'] )

        return self.sweep_pts


def darkesr(folder=os.getcwd()):
    a = DarkESRAnalysis(folder)
    a.get_sweep_pts()
    a.get_ssro_results()
    a.plot_result_vs_sweepparam()

def rabi(folder=os.getcwd()):
    a = ElectronRabiAnalysis(folder)
    a.get_sweep_pts()
    a.get_ssro_results()
    a.plot_result_vs_sweepparam()
