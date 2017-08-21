####################################################################################################################################################################
# general analysis script for cavity measurements. 
# SvD, 2016

import numpy as np

from analysis.lib.m2 import m2;reload(m2)

class cavity_analysis(m2.M2Analysis):

    def get_x_pts(self):
        self.nr_x_pts = self.g.attrs['nr_steps']
        self.x_pts = np.linspace(self.g.attrs['start_voltage'],self.g.attrs['end_voltage'],self.g.attrs['nr_steps'])

    def get_sweep_pts(self):
        self.sweep_name = self.g.attrs['sweep_name']
        self.sweep_pts = self.g.attrs['sweep_pts']
        self.nr_syncs_per_pt = self.g.attrs['nr_syncs_per_pt']
        print self.sweep_name, self.sweep_pts

    def get_sweep_data (self):
        """
        function that loads the data from a sweep measurement
        """
        self.nr_scans = self.g.attrs['nr_scans']
        self.nr_repetitions = self.g.attrs['nr_repetitions']

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

    def get_scan_data(self, name):
        grp = self.adwingrp(name) 
        self.frq_data = grp['laser_frequency']
        self.x_data = np.linspace(self.g.attrs['start_voltage'], self.g.attrs['end_voltage'], self.g.attrs['nr_steps'])
        self.y_data = np.reshape(np.array(grp['photodiode_voltage']),(self.g.attrs['nr_scans'], self.g.attrs['nr_steps']))
        self.y_data_per_ms = np.reshape(grp['photodiode_voltage_ms'], (self.g.attrs['nr_scans'], self.g.attrs['nr_steps'], self.g.attrs['nr_ms_per_point']))
        self.average_data() #necessary since nr_scans can be >1
        return self.x_data,self.y_data,self.frq_data,self.y_data_per_ms

    def get_lengthscan_data(self):
        """
        function that loads the data from a single length scan measurement
        """
        self.get_scan_data('lengthscan')
        return self.x_data,self.y_data

    def get_laserscan_data(self):
        """
        function that loads the data from a single laserscan measurement
        """
        self.get_scan_data('laserscan')
        return self.x_data,self.y_data,self.frq_data,self.y_data_per_ms

    def average_data(self):
        self.y_data = np.mean(self.y_data,axis=0)
        self.y_data_per_ms = np.mean(self.y_data_per_ms,axis=0)

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
        self.total_sweep_time = self.g.attrs['nr_steps']*(self.g.attrs['ADC_averaging_cycles']+self.g.attrs['wait_cycles'])/self.g.attrs['cycle_duration']*3/100 # in seconds
        return self.total_sweep_time