####################################################################################################################################################################
# A.M.J. Zwerver
#a.m.j.zwerver@student.tudelft.nl
# October 2015

import os
import h5py
import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
import pylab as plt
from analysis.lib import fitting
from analysis.lib.fitting import esr, common
from analysis.lib.fitting import fit,esr, common
from matplotlib import rc, cm
from scipy import interpolate
from scipy.optimize import curve_fit
import operator

from analysis.lib.m2 import m2

class cavity_analysis(m2.M2Analysis):

    def get_x_pts(self):
        self.nr_x_pts = self.f.attrs['nr_steps']
        self.x_pts = np.linspace(self.f.attrs['start_voltage'],self.f.attrs['end_voltage'],self.f.attrs['nr_steps'])

    def get_sweep_pts(self):
        self.sweep_name = self.f.attrs['sweep_name']
        self.sweep_pts = self.f.attrs['sweep_pts']
        self.nr_syncs_per_pt = self.f.attrs['nr_syncs_per_pt']
        print self.sweep_name, self.sweep_pts

    def get_sweep_data (self):
        """
        function that loads the data from a sweep measurement
        """
        self.nr_scans = self.f.attrs['nr_scans']
        self.nr_repetitions = self.f.attrs['nr_repetitions']

        self.sweep_data=np.zeros([len(self.sweep_pts),self.nr_repetitions,self.nr_x_pts])

        for i in np.arange(self.nr_repetitions):
            for j in np.arange(len(self.sweep_pts)):
                first_rep = str(int(i*self.nr_scans+1))
                last_rep = str(int(i*self.nr_scans+self.nr_scans))
                data_index_name = 'sweep_pt_'+str(j)+'_reps_'+first_rep+'-'+last_rep  
                grp = self.f['raw_data_'+data_index_name]
                for k in np.arange(self.nr_scans):
                    self.single_scan_data = grp['scannr_'+str(k+1)].value
                    #fill the sweep_data array for j = sweep pt, i*self.nr_scans+k = repetition nr
                    self.sweep_data[j,int(i*self.nr_scans+k)] = self.single_scan_data

        self.avg_sweep_data = np.average(self.sweep_data,axis=1)

        return self.sweep_data

    def get_lengthscan_data(self):
        """
        function that loads the data from a single length scan measurement
        """
        self.x_data = self.g['piezo_voltage']
        self.y_data = self.g['PD_signal']
        return self.x_data,self.y_data

    def get_laserscan_data(self):
        """
        function that loads the data from a single laserscan measurement
        """
        self.frq_data = self.g['laser_frequency']
        self.x_data = self.g['laser_tuning_voltage']
        self.y_data = np.array(self.g['PD_signal'])
        self.y_data_per_ms = self.g['PD_signal_per_ms'][:]
        return self.x_data,self.y_data,self.frq_data,self.y_data_per_ms

    def bin_data_per_ms(self,binsize=100):
        self.binsize = binsize
        if np.size(self.y_data_per_ms,1)%binsize !=0:
            del_points = np.remainder(np.size(self.y_data_per_ms,1),binsize)
            print 'deleting %d points to enable reshaping'%del_points
            self.y_data_per_ms = self.y_data_per_ms[:,:-del_points]
        self.nr_bins =(np.size(self.y_data_per_ms,1)/binsize)
        #print 'using %d bins of size %d'%(self.nr_bins,self.binsize)
        # try:
        y_data_per_ms_reshaped = self.y_data_per_ms.reshape(np.size(self.y_data_per_ms,0),self.nr_bins,binsize)
        # except ValueError:
        #     print 'valueerror'

        self.y_data_per_ms_binned = np.average(y_data_per_ms_reshaped,axis=2)
        return self.y_data_per_ms_binned

    def avg_bins_per_s(self,**kw):
        max_bins = kw.pop('max_bins',1000/self.binsize)
        if self.nr_bins%max_bins==0:
            bin_repetitions = self.nr_bins/max_bins
            self.nr_bins = max_bins
            y_data_per_ms_reshaped = self.y_data_per_ms_binned.reshape(np.size(self.y_data_per_ms,0),bin_repetitions,self.nr_bins)
            self.y_data_per_ms_binned = np.average(y_data_per_ms_reshaped,axis=1)
        return self.y_data_per_ms_binned

    def calc_sweeptime(self):
        self.total_sweep_time = self.f.attrs['nr_steps']*(self.f.attrs['ADC_averaging_cycles']+self.f.attrs['wait_cycles'])/self.f.attrs['cycle_duration']*3/100 # in seconds
        return self.total_sweep_time