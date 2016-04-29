
'''
Analysis program for the raw 2d data from the p7889
'''
import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
#from analysis.lib import fitting
#from analysis.lib.m2.ssro import ssro
#from analysis.lib.math import error
from analysis.lib.m2 import m2
from analysis.lib.tools import toolbox


class DarkESRanalysis(m2.M2Analysis):

    def get_data(self, name):     
        grp = self.adwingrp(name)
        self.sweep_name = grp.attrs['sweep_name']
        self.sweep_pts = grp.attrs['sweep_pts']
        self.pts = grp.attrs['pts']
        self.mmt_type = grp.attrs['measurement_type']
        self.counts = grp['Counts'].value
        self.frq = grp[self.sweep_name].value
        self.t = grp['t (ns)'].value


    def eval_data(self,Eval_ROI_start,Eval_ROI_end):

        #do post-processing
        #grp1=h5.DataGroup('analyzed data',self.h5data,base=self.h5base)
        counts = self.counts
        frq = self.frq
        time= self.t

        if Eval_ROI_start<np.amin(time):
            startind=[0]
        else:
            startind=np.where(time==Eval_ROI_start)
        if Eval_ROI_end>np.amax(time):
            endind=[time.size-1]
        else:
            endind=np.where(time==Eval_ROI_end)

        #strip the count array of points in time that lie outside of the specified ROIs
        print startind[0], endind[0]
        sliced=counts[:,startind[0]:endind[0]]
        self.summation=np.array(range(self.pts))

        for i,row in enumerate(sliced):
            self.summation[i]=np.sum(row)

        #grp1.add_coordinate(name=self.params['sweep_name'],data=y)
        #grp1.add_value(name='Counts in ROI',data=summation)

    def get_raw_counts(self):
        return self.counts    
    def plot_data(self,Eval_ROI_start,Eval_ROI_end,ax=None):
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.errorbar(self.sweep_pts, self.summation, yerr=np.sqrt(self.summation),fmt='.-')
        ax.set_xlabel(self.sweep_name)
        ax.set_ylabel('counts')
        ax.set_title(self.mmt_type+'ROI_'+str(Eval_ROI_start)+'-'+str(Eval_ROI_end))

        fig.savefig(os.path.join(self.folder,
            self.mmt_type+'_'+self.name+'_ROI_'+str(Eval_ROI_start)+'-'+str(Eval_ROI_end)+'.png'),
            format='png')
        return ax
    def save_data(self):
        f = self.analysis_h5data()
        if not 'processed_data' in f:
            f.create_group('processed_data')

        g = f['/processed_data']
        if 'counts_in_ROI' in g:
            del  g['counts_in_ROI']
        if 'frequencies' in g:
            del  g['frequencies']

        g['counts_in_ROI'] = self.summation
        g['frequencies'] = self.frq
        f.close()


def Analyse_DarkESR(ROI_start, ROI_end, guess_ctr, 
    folder='',
    guess_splitN = 2e-3,
    guess_offset = 1,
    guess_width = 0.2e-3,
    guess_amplitude = 0.3,
    ret = 'f0'):
    if folder=='':
        folder = toolbox.latest_data('Rabi')
    print folder
    a = DarkESRanalysis(folder)
    a.get_data('raw data')
    a.eval_data(ROI_start,ROI_end)
    a.plot_data(ROI_start,ROI_end)
    a.save_data()
    # plt.close('all')
    a.finish()
    return a.get_raw_counts()

def Fit_DarkESR(folder, 
    guess_ctr,
    ROI_start = 1670,
    ROI_end = 2200,
    guess_splitN = 2e6,
    guess_offset = 1,
    guess_width = 0.2e-3,
    guess_amplitude = 0.3,
    ret = 'f0'
    ):
    """
    Fit Gaussians to Dark ESR result.
    """

    # Fetch & postprocess data
    # a = dea.Analyse_DarkESR(folder)
    a = DarkESRanalysis(folder)
    a.get_data('raw_data')
    a.eval_data(ROI_start, ROI_end)
    x = a.sweep_pts
    y = a.summation
    y_normalized = y / np.amax(y)

    # Fit Gaussians
    j = 0
    while y_normalized[j]>0.97 and j < len(y)-2: # find indices for 'non-dips'
        k = j
        j = j+1

    if k > len(y)-5:
        print 'Could not find dip'
        return
    else:
        guess_ctr = x[k]+ guess_splitN #convert to GHz and go to middle dip
        print 'guess_ctr = '+str(guess_ctr)

    fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
            guess_amplitude, guess_width, guess_ctr,
            (3, guess_splitN),
            do_print=True, ret=True, fixed=[4])


    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']

        return f0, u_f0
