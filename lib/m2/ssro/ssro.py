"""
module to postprocess and analyze the resulting data
from an adwin SSRO measurement (measurement2 generation).

"""

### imports
import sys, os, time
import numpy as np
import h5py
from matplotlib import pyplot as plt
from analysis.lib.tools import toolbox
from analysis.lib.m2 import m2


# some general tools
def get_SSRO_calibration(folder, readout_time):
    fp = os.path.join(folder, 'analysis.hdf5')
    f = h5py.File(fp, 'r')
    
    times = f['fidelity/time'].value
    fids0 = f['fidelity/ms0'].value
    fids1 = f['fidelity/ms1'].value

    tidx = np.argmin(abs(times-readout_time))
    f0 = fids0[tidx,1]
    u_f0 = fids0[tidx,2]
    f1 = fids1[tidx,1]
    u_f1 = fids1[tidx,2]

    return f0, u_f0, f1, u_f1


# analysis classes and shortcut functions
class SSROAnalysis(m2.M2Analysis):
    
    def get_run(self, name):
        grp = self.adwingrp(name)
        
        self.reps = grp['completed_reps'].value
        self.cycle_duration = grp.attrs['cycle_duration']
        
        self.repetitions = np.arange(self.reps)
        self.ro_time = self.ssro_time = np.arange(grp.attrs['SSRO_duration']) * \
                self.cycle_duration * 10./3.
        self.sp_time = np.arange(grp.attrs['SP_duration']) * \
                self.cycle_duration * 10./3.

        self.binsize = grp.attrs['cycle_duration']/300.
        self.ro_counts = grp['RO_data'].value.reshape(
                self.reps, grp.attrs['SSRO_duration'])
        self.cr_counts = grp['CR_after'].value
        self.sp_counts = grp['SP_hist'].value

    def cpsh_hist(self, ro_counts, reps, firstbin=0, lastbin=-1, plot=True, **kw):
        name = kw.pop('name', '')
 
        title_suffix = ': '+name if name != '' else ''
        fn_suffix = '_'+name if name != '' else ''

        cpsh = np.sum(ro_counts[:,firstbin:lastbin], axis=1)

        annotation = 'repetitions: %d' % reps
        mean_cpsh = sum(cpsh)/float(reps)
        annotation += "\nmean c.p.sh. = %.2f" % (mean_cpsh)
        pzero = len(np.where(cpsh==0)[0])/float(reps)
        annotation += "\np(0 counts) = %.2f" % (pzero)

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(cpsh, max(cpsh)+1, align='mid', label='counts', 
                    normed=True)
            ax.set_xlabel('cts/shot')
            ax.set_ylabel('probability')
            ax.set_title(self.default_plot_title + title_suffix)

            plt.figtext(0.85, 0.85, annotation, horizontalalignment='right',
                verticalalignment='top')

            fig.savefig(os.path.join(self.folder, 
                'cpsh'+fn_suffix+'.'+self.plot_format), 
                format=self.plot_format)

    def readout_relaxation(self, ro_time, ro_counts, reps, binsize, 
            lastbin=-1, plot=True, **kw):
        
        name = kw.pop('name', '')
        ret = kw.pop('ret', False)
 
        title_suffix = ': '+name if name != '' else ''
        fn_suffix = '_'+name if name != '' else ''

        ro_countrate = np.sum(ro_counts[:,:lastbin], axis=0) / \
                (binsize*1e-6*reps)
        ro_time = ro_time[:lastbin]

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(ro_time, ro_countrate, 'o')
            ax.set_xlabel('RO time [ns]')
            ax.set_ylabel('counts [Hz]')
            ax.set_title(self.default_plot_title + title_suffix)

            fig.savefig(os.path.join(self.folder, 
                'readout_relaxation'+fn_suffix+'.'+self.plot_format), 
                format=self.plot_format)

        if ret:
            return ro_time, ro_countrate

    
    def spinpumping(self, sp_time, sp_counts, reps, binsize,
            plot=True, **kw):

        name = kw.pop('name', '')
 
        title_suffix = ': '+name if name != '' else ''
        fn_suffix = '_'+name if name != '' else ''
        
        sp = np.array([j/(binsize*1e-6*reps) for j in sp_counts])
        
        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(sp_time+1E3, sp, 'o')
            ax.set_xlabel('spin pumping time [ns]')
            ax.set_ylabel('counts [Hz]')
            ax.set_title(self.default_plot_title + title_suffix)
            
            fig.savefig(os.path.join(self.folder, 
                'spinpumping'+fn_suffix+'.'+self.plot_format), 
                format=self.plot_format)

    def charge_hist(self, cr_counts, plot=True, **kw):
        
        name = kw.pop('name', '')
 
        title_suffix = ': '+name if name != '' else ''
        fn_suffix = '_'+name if name != '' else ''

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(cr_counts, abs(max(cr_counts)-min(cr_counts)+1), 
                    normed=True)
            ax.set_xlabel('counts during CR check')
            ax.set_ylabel('probability')
            ax.set_title(self.default_plot_title + title_suffix)
            
            fig.savefig(os.path.join(self.folder, 
                'cr_check'+fn_suffix+'.'+self.plot_format), 
                format=self.plot_format)

    def fidelity(self, ro_counts, reps, binsize, ms,
            lastbin=None, plot=True, **kw):
        
        name = kw.pop('name', '')
 
        title_suffix = ': '+name if name != '' else ''
        fn_suffix = '_'+name if name != '' else ''
        dataset_name = name if name != '' else 'fidelity'

        if lastbin == None:
            lastbin = len(ro_counts[0,:])

        fid_dat = np.zeros((0,3))

        for i in range(1,lastbin):
            t = i*binsize
            
            # d: hist of counts, c: counts per shot
            d = np.sum(ro_counts[:,:i], axis=1)
            c = np.sum(d)/float(reps)

            # we get the fidelity from the probability to get zero counts in a
            # shot 
            pzero = len(np.where(d==0)[0])/float(reps)
            pzero_err = np.sqrt(pzero*(1-pzero)/reps)
            fid = 1-pzero if ms == 0 else pzero # fidelity calc. depends on ms
            fid_dat = np.vstack((fid_dat, np.array([[t, fid, pzero_err]])))

        f = self.analysis_h5data()
        if not 'fidelity' in f:
            f.create_group('fidelity')

        g = f['/fidelity']
        if dataset_name in g:
            del  g[dataset_name]
        
        g[dataset_name] = fid_dat
        f.close()

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.errorbar(fid_dat[:,0], fid_dat[:,1], fmt='o', yerr=fid_dat[:,2])
            ax.set_xlabel('RO time [us]')
            ax.set_ylabel('ms = %d RO fidelity' % ms)
            ax.set_title(self.default_plot_title + title_suffix)

            fig.savefig(os.path.join(self.folder, 
                'fidelity'+fn_suffix+'.'+self.plot_format), 
                format=self.plot_format)

    def mean_fidelity(self, plot=True, **kw):
        
        f = self.analysis_h5data()
        g = f['/fidelity']
        _fid0 = g['ms0']
        _fid1 = g['ms1']

        time = _fid0[:,0]
        fid0 = _fid0[:,1]
        fid0_err = _fid0[:,2]
        fid1 = _fid1[:,1]
        fid1_err = _fid1[:,2]

        F = (fid0 + fid1)/2.
        F_err = np.sqrt(fid0_err**2 + fid1_err**2)
        F_max = max(F)
        t_max = time[F.argmax()]

        meanfid = (_fid0[:,1]+_fid1[:,1])*0.5
        meanfiderr = np.sqrt( (0.5*_fid0[:,2])**2 + (0.5*_fid1[:,2])**2 )

        try:
            del g['time']
            del g['mean_fidelity']
            del g['mean_fidelity_err']
        except:
            pass

        g['time'] = time
        g['mean_fidelity'] = meanfid
        g['mean_fidelity_err'] = meanfiderr
        f.close()

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.errorbar(time, fid0, fmt='.', yerr=fid0_err, label='ms=0')
            ax.errorbar(time, fid1, fmt='.', yerr=fid1_err, label='ms=1')
            ax.errorbar(time, F, fmt='.', yerr=F_err, label='mean')
            ax.set_xlabel('RO time (us)')
            ax.set_ylabel('RO fidelity')
            ax.set_ylim((0.5,1))
            ax.legend(loc=4)
            plt.figtext(0.8, 0.5, "max. F=%.2f at t=%.2f us" % (F_max, t_max),
                    horizontalalignment='right')

            ax.set_title(self.default_plot_title + ': mean RO fidelity')

            fig.savefig(os.path.join(self.folder, 
                'mean_fidelity.'+self.plot_format), 
                format=self.plot_format)


def ssrocalib(folder=''):
    if folder=='':
        folder=get_latest_data()
    a = SSROAnalysis(folder)
    
    for n,ms in zip(['ms0', 'ms1'], [0,1]): #zip((['ms0'], [0]):#
        a.get_run(n)
        a.cpsh_hist(a.ro_counts, a.reps, name=n)
        a.readout_relaxation(a.ro_time, a.ro_counts, a.reps, a.binsize, name=n)
        a.spinpumping(a.sp_time, a.sp_counts, a.reps, a.binsize, name=n)
        a.charge_hist(a.cr_counts, name=n)
        a.fidelity(a.ro_counts, a.reps, a.binsize, ms, name=n)
    
    plt.close('all')
    a.mean_fidelity()
    a.finish()
    
def get_latest_data(string = 'AdwinSSRO', datapath = ''):
    meas_folder = r'D:\measuring\data'
    currdate = time.strftime('%Y%m%d')
    
    if datapath == '':
        df = os.path.join(meas_folder, currdate)
    else:
        df = datapath
    
    right_dirs = list()

    if os.path.isdir(df):
        for k in os.listdir(df):
            if string in k:
                right_dirs.append(k)
        
        if len(right_dirs) > 0:
            latest_dir = os.path.join(df,right_dirs[len(right_dirs)-1])
        else:
            print 'No measurements containing %s in %s'%(string, df)
        
        print '\nAnalyzing data in %s'%latest_dir

    else:
        print 'Folder %s does not exist'%df
        latest_dir = False

    return latest_dir
    
class AWGSSROAnalysis(m2.M2Analysis):

    def get_count_probability(self, name):
        grp = self.adwingrp(name)
        
        self.reps = grp['completed_reps'].value
        self.pts = grp.attrs['pts']
        self.times = self.g.attrs['AWG_SSRO_durations']
        self.counts = grp['RO_data'].value.reshape((-1,self.pts))

        self.count_probability = np.array([float(len(np.where(self.counts[:,i]>0)[0])) \
                for i in range(self.pts)]) / grp.attrs['SSRO_repetitions']

        self.u_count_probability = np.sqrt(
                self.count_probability*(1.-self.count_probability)/\
                        grp.attrs['SSRO_repetitions'])

def awgssro(folder, ymax=1.):
    a = AWGSSROAnalysis(folder)

    fig = a.default_fig()
    ax = a.default_ax(fig)

    for n,ms in zip(['ms0', 'ms1'], [0,1]):
        a.get_count_probability(n)
        ax.errorbar(a.times, a.count_probability, yerr=a.u_count_probability,
                fmt='o', label=n)

    ax.set_xlim(0,max(a.times))
    ax.set_ylim(0,ymax)
    ax.set_xlabel('RO time (ns)')
    ax.set_ylabel('Cnt. prob.')
    ax.legend()

def awgssro_prjprob(folder, pop0=1./6):
    a = AWGSSROAnalysis(folder)

    fig = a.default_fig()
    ax = a.default_ax(fig)

    a.get_count_probability('ms0')
    cpr0 = a.count_probability

    a.get_count_probability('ms1')
    cpr1 = a.count_probability

    prjprob = cpr0*pop0 / (cpr0*pop0 + cpr1*(1.-pop0))
    ax.plot(a.times, prjprob, 'o')

    ax.set_xlabel('RO time (ns)')
    ax.set_ylabel('Prob. for projection into 0')

        
