
'''
Analysis program for the raw 2d data from the p7889
'''
import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib.m2 import m2
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot


class Rabi_DESR_analysis(m2.M2Analysis):
    adwin_dataname = 'raw data'
    ROI_start = 1670
    ROI_end = 2200

    # def __init__(self, name, folder, ROI_start, ROI_end):
    #     m2.M2Analysis.__init__(self, folder)
    #     self.name = name
    #     self.ROI_start = ROI_start
    #     self.ROI_end = ROI_end
        
    def get_data(self):     
        grp = self.adwingrp(self.adwin_dataname)
        self.sweep_name = grp.attrs['sweep_name']
        self.sweep_pts = grp.attrs['sweep_pts']
        self.pts = grp.attrs['pts']
        self.mmt_type = grp.attrs['measurement_type']
        self.counts = grp['Counts'].value
        self.sweep = grp[self.sweep_name].value
        self.t = grp['t (ns)'].value

    def eval_data(self):#,Eval_ROI_start,Eval_ROI_end):

        #do post-processing
        counts = self.counts
        sweep = self.sweep
        time= self.t

        if self.ROI_start<np.amin(time):
            startind=[0]
        else:
            startind=np.where(time==self.ROI_start)
        if self.ROI_end>np.amax(time):
            endind=[time.size-1]
        else:
            endind=np.where(time==self.ROI_end)

        #strip the count array of points in time that lie outside of the specified ROIs
        print startind[0], endind[0]
        sliced=counts[:,startind[0]:endind[0]]
        self.summation=np.array(range(self.pts))

        for i,row in enumerate(sliced):
            self.summation[i]=np.sum(row)

    def get_raw_counts(self):
        return self.counts    

    def plot_data(self):#,Eval_ROI_start,Eval_ROI_end,ax=None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.errorbar(self.sweep_pts, self.summation, yerr=np.sqrt(self.summation),fmt='.-')
        ax.set_xlabel(self.sweep_name)
        ax.set_ylabel('counts')
        ax.set_title(self.mmt_type+'ROI_'+str(self.ROI_start)+'-'+str(self.ROI_end))

        fig.savefig(os.path.join(self.folder,
            self.mmt_type+'_'+self.name+'_ROI_'+str(self.ROI_start)+'-'+str(self.ROI_end)+'.png'),
            format='png')
        return fig, ax

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
        g['frequencies'] = self.sweep
        f.close()

    def prepare_data(self):
        self.get_data()
        self.eval_data()
        self.fig, self.ax = self.plot_data()
        self.x = self.sweep_pts
        self.y = self.summation

def Fit_Rabi(folder,
    ROI_start = 1670,
    ROI_end = 2200,
    guess_frq = 1./8,
    guess_amp = 0,
    guess_of = 1,
    guess_phi = 0,
    guess_k = 0.2,
    fixed = []
    ):
    """
    Fit damped oscillation to Rabi signal and plots results.

    Arguments:
    folder      --  data directory
    ROI_start   --  lower integration bound for photon counts
    ROI_end     --  upper integration bound 
    guess_frq   --  guessed Rabi frequency
    guess_amp   --  guessed amplitude. If 0, take max(photon_counts) as guess.
    guess_of    --  guessed offset 
    guess_phi   --  guessed phase
    guess_k     --  guessed decay constant
    fixed       --  fixes fitting parameters via [frq, amp, phi, o, k]   
    """

    a = Rabi_DESR_analysis(folder)
    a.ROI_start = ROI_start
    a.ROI_end = ROI_end

    a.prepare_data()

    # Fit & plot
    o = fit.Parameter(guess_of, 'o')
    f = fit.Parameter(guess_frq, 'f')
    if guess_amp == 0:
        A = fit.Parameter(np.amax(a.summation), 'A')
    else:
        A = fit.Parameter(guess_amp, 'A')
    phi = fit.Parameter(guess_phi, 'phi')
    k = fit.Parameter(guess_k, 'k')
    p0 = [f, A, phi, o, k]
    fitfunc_str = 'o - A + A*e^(-kx)*cos(2pi (fx-phi))'

    def fitfunc(x):
        return (o()-A()) + A() * np.exp(-k()*x) * np.cos(2*np.pi*(f()*x - phi()))

    fit_result = fit.fit1d(a.x, a.y, None, p0=p0, fitfunc=fitfunc, fixed=fixed,
            do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(0,a.x[-1],201), ax=a.ax,
        plot_data=False)

    a.fig.savefig(os.path.join(a.folder,
            a.mmt_type+'_'+ str(a.name) +'_ROI_'+str(a.ROI_start)+'-'+str(a.ROI_end)+'.png'),
            format='png')

    print "Saved figure to %s" % (os.path.join(a.folder,
            a.mmt_type+'_'+ str(a.name) +'_ROI_'+str(a.ROI_start)+'-'+str(a.ROI_end)+'.png'))

    a.finish()

def Fit_DESR(folder,
    guess_ctr,
    ROI_start = 1670,
    ROI_end = 2200,
    splitN = False,
    splitC = True,
    guess_splitN = 2e-3,
    guess_splitC = 10e-6,
    guess_offset = 0,
    guess_width = 50e-6,
    guess_amplitude = 0,
    fixed = [],
    ret = 'f0'
    ):
    """
    Fits Gaussians to Dark ESR measurement and plots results.

    Arguments (see also 'Fit_Rabi' function):
    guess_ctr       -- guessed center frequency of DESR dip
    splitN          -- nitrogen splitting present (True/False)
    splitC          -- carbon splitting present (True/False)
    guess_splitN    -- guessed nitrogen splitting
    guess_splitC    -- guessed carbon splitting
    guess_offset    -- guessed offset. If 0, take max(photon counts) as guess
    guess_width     -- guess width of resonance dip. Default is 50 KHz (carbon)
    guess_amplitude -- guess dip amplitude. If 0, take 3% of offset
    fixed           -- fixes parameters according to [offset, ampl, width, ctr, splitN, splitC]
    ret             -- If 'f0', output fitted results
    """

    a = Rabi_DESR_analysis(folder)
    a.ROI_start = ROI_start
    a.ROI_end = ROI_end

    a.prepare_data() 

    # Set default guesses
    if guess_offset == 0:
        guess_offset = np.amax(a.summation)
    if guess_amplitude == 0:
        guess_amplitude = guess_offset * 0.03

    # Fit & plot
    if splitN == True:
        fit_result = fit.fit1d(a.x, a.y, esr.fit_ESR_gauss, guess_offset,
        guess_amplitude, guess_width, guess_ctr,
        (3, guess_splitN),
        do_print=True, ret=True, fixed=fixed)
    elif splitC == True:
        fit_result = fit.fit1d(a.x, a.y, esr.fit_ESR_gauss, guess_offset,
        guess_amplitude, guess_width, guess_ctr,
        (2, guess_splitC),
        do_print=True, ret=True, fixed=fixed)
    else: # Fit without splitting
        fit_result = fit.fit1d(a.x, a.y, esr.fit_ESR_gauss, guess_offset,
        guess_amplitude, guess_width, guess_ctr,
        do_print=True, ret=True, fixed=fixed)


    plot.plot_fit1d(fit_result, np.linspace(np.amin(a.x), np.amax(a.x), 1000), ax=a.ax, plot_data=False)

    a.fig.savefig(os.path.join(a.folder,
            a.mmt_type+'_'+ str(a.name) +'_ROI_'+str(a.ROI_start)+'-'+str(a.ROI_end)+'.png'),
            format='png')

    print "Saved figure to %s" % (os.path.join(a.folder,
            a.mmt_type+'_'+ str(a.name) +'_ROI_'+str(a.ROI_start)+'-'+str(a.ROI_end)+'.png'))

    if ret == 'f0':
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']
    print fit_result['params_dict']['A']/np.sqrt(fit_result['params_dict']['A'])
    print fit_result['params_dict']['A']/fit_result['params_dict']['a']
    a.finish()
    print f0, u_f0

if __name__ == '__main__':
    Fit_Rabi(folder, fixed = [])    