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

# Some general tools

def get_SSRO_MWInit_calibration(folder, readout_time,el_state):
    if 'analysis.hdf5' in str(folder):
        fp = folder
    else:
        fp = os.path.join(folder, 'analysis.hdf5')
    f = h5py.File(fp, 'r')

    times = f['fidelity/time'].value
    fids0 = f['fidelity/ms0'].value
    fids1 = f['fidelity/ms' + str(el_state[-2:])].value
    
    tidx = np.argmin(abs(times-readout_time))
    f0 = fids0[tidx,1]
    u_f0 = fids0[tidx,2]
    f1 = fids1[tidx,1]
    u_f1 = fids1[tidx,2]

    return f0, u_f0, f1, u_f1


def get_SSRO_calibration(folder, readout_time):
    if 'analysis.hdf5' in str(folder):
        fp = folder
    else:
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
        save = kw.pop('save', True)

        title_suffix = ': '+name if name != '' else ''
        fn_suffix = '_'+name if name != '' else ''

        cpsh = np.sum(ro_counts[:,firstbin:lastbin], axis=1)

        annotation = 'repetitions: %d' % reps
        mean_cpsh = sum(cpsh)/float(reps)
        annotation += "\nmean c.p.sh. = %.2f" % (mean_cpsh)
        pzero = len(np.where(cpsh==0)[0])/float(reps)
        annotation += "\np(0 counts) = %.2f" % (pzero)

        if save:
            f = self.analysis_h5data()
            if not 'cpsh' in f:
                    f.create_group('cpsh')
            g = f['/cpsh']
            if name in g:
                del g[name]
            g[name]=  cpsh
            f.close()
        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(cpsh, np.arange(max(cpsh)+2)-0.5, align='mid', label='counts',
                    normed=True) # , stacked=True)
            ax.set_xlabel('cts/shot')
            ax.set_ylabel('probability')
            ax.set_title(self.default_plot_title + title_suffix)
            ax.set_xlim(-0.5, max(cpsh)+0.5)
            if name == 'ms1':
                ax.set_yscale('log')

            plt.figtext(0.85, 0.85, annotation, horizontalalignment='right',
                verticalalignment='top')
            if save:
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
        save = kw.pop('save', True)
        ret = kw.pop('ret', False)

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

        if save:
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

        if ret:
            return fid_dat


    def mean_fidelity(self, plot=True, plot_photon_ms0=True, **kw):

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
        maxidx = F.argmax()
        F_max = F[maxidx]
        F_max_err = F_err[maxidx]
        t_max = time[F.argmax()]

        meanfid = (_fid0[:,1]+_fid1[:,1])*0.5
        meanfiderr = np.sqrt( (0.5*_fid0[:,2])**2 + (0.5*_fid1[:,2])**2 )

        print 'SSRO calibration : ', self.timestamp
        print 'max. F = ({:.2f} +/- {:.2f})% at t={:.0f} us'.format(F_max*100., F_max_err*100., t_max)
        print '\tms_0 = ({:.2f} +/- {:.2f})%'.format(fid0[maxidx]*100, fid0_err[maxidx]*100)
        print '\tms_1 = ({:.2f} +/- {:.2f})%'.format(fid1[maxidx]*100, fid1_err[maxidx]*100)

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
            plt.figtext(0.8, 0.5, "max. F=({:.2f} +/- {:.2f})% at t={:.0f} us".format(F_max*100., F_max_err*100., t_max),
                    horizontalalignment='right')

            ax.set_title(self.default_plot_title + ': mean RO fidelity')

            fig.savefig(os.path.join(self.folder,
                'mean_fidelity.'+self.plot_format),
                format=self.plot_format)
            # fig.savefig(os.path.join(self.folder,
            #     'mean_fidelity.pdf'),
            #     format='pdf')

        if plot_photon_ms0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            Prob_ms0 = (fid0[1:])/(fid0[1:]+(1-fid1[1:]))
            max_Prob_ms0 = Prob_ms0.max()
            time_max_Prob_ms0 = time[Prob_ms0.argmax()+1]
            ax.errorbar(time[1:], Prob_ms0, fmt='.', yerr=0*Prob_ms0)
            ax.set_xlabel('RO time (us)')
            ax.set_ylabel('Prob. photon came from ms=0')
            ax.set_ylim((0.9,1))
            plt.figtext(0.8, 0.5, "max. {:.2f} at t={:.0f} us".format(max_Prob_ms0*100., time_max_Prob_ms0),
                    horizontalalignment='right')

            ax.set_title(self.default_plot_title + ': Probability_photon_from_ms0')

            fig.savefig(os.path.join(self.folder,
                'projectivity.'+self.plot_format),
                format=self.plot_format) 

    def mean_fidelity_MWInit(self, plot=True, plot_photon_ms0=True, **kw):


        # it's not beautiful but works atm. Should be made more versatile
        f = self.analysis_h5data()
        # print f
        g = f['/fidelity']

        # print [int(g) for g in str.split() if g.isdigit()]
        _fid0 = g['ms0']
        _fidp1 = g['msp1']
        _fidm1 = g['msm1']

        time = _fid0[:,0]
        fid0 = _fid0[:,1]
        fid0_err = _fid0[:,2]
        fidp1 = _fidp1[:,1]
        fidp1_err = _fidp1[:,2]
        fidm1 = _fidm1[:,1]
        fidm1_err = _fidm1[:,2]

        Fm = (fid0 + fidm1)/2.
        Fm_err = np.sqrt(fid0_err**2 + fidm1_err**2)
        Fp = (fid0 + fidp1)/2.
        Fp_err = np.sqrt(fid0_err**2 + fidp1_err**2)

        maxidxm = Fm.argmax()
        Fm_max = Fm[maxidxm]
        Fm_max_err = Fm_err[maxidxm]
        tm_max = time[Fm.argmax()]

        maxidxp = Fp.argmax()
        Fp_max = Fp[maxidxp]
        Fp_max_err = Fp_err[maxidxp]
        tp_max = time[Fp.argmax()] 

        meanfidm = (_fid0[:,1]+_fidm1[:,1])*0.5
        meanfidmerr = np.sqrt( (0.5*_fid0[:,2])**2 + (0.5*_fidm1[:,2])**2 )

        meanfidp = (_fid0[:,1]+_fidp1[:,1])*0.5
        meanfidperr = np.sqrt( (0.5*_fid0[:,2])**2 + (0.5*_fidp1[:,2])**2 )

        print 'SSRO calibration : ', self.timestamp
        print 'max. Fm = ({:.2f} +/- {:.2f})% at t={:.0f} us'.format(Fm_max*100., Fm_max_err*100., tm_max)
        print '\tms_0 = ({:.2f} +/- {:.2f})%'.format(fid0[maxidxm]*100, fid0_err[maxidxm]*100)
        print '\tms_1 = ({:.2f} +/- {:.2f})%'.format(fidm1[maxidxm]*100, fidm1_err[maxidxm]*100)

        print 'max. Fp = ({:.2f} +/- {:.2f})% at t={:.0f} us'.format(Fp_max*100., Fp_max_err*100., tp_max)
        print '\tms_0 = ({:.2f} +/- {:.2f})%'.format(fid0[maxidxp]*100, fid0_err[maxidxp]*100)
        print '\tms_1 = ({:.2f} +/- {:.2f})%'.format(fidp1[maxidxp]*100, fidp1_err[maxidxp]*100)

        try:
            del g['time']
            del g['mean_fidelity']
            del g['mean_fidelity_err']
        except:
            pass

        g['time'] = time
        g['mean_fidelity'] = meanfidm
        g['mean_fidelity_err'] = meanfidmerr
        f.close()

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.errorbar(time, fid0, fmt='.', yerr=fid0_err, label='ms=0')
            ax.errorbar(time, fidm1, fmt='.', yerr=fidm1_err, label='ms=-1')
            ax.errorbar(time, Fm, fmt='.', yerr=Fm_err, label='mean')
            ax.set_xlabel('RO time (us)')
            ax.set_ylabel('RO fidelity')
            ax.set_ylim((0,1))
            ax.legend(loc=4)
            plt.figtext(0.8, 0.5, "max. Fm1=({:.2f} +/- {:.2f})% at t={:.0f} us".format(Fm_max*100., Fm_max_err*100., tm_max),
                    horizontalalignment='right')

            ax.set_title(self.default_plot_title + ': mean RO fidelity')

            fig.savefig(os.path.join(self.folder,
                'mean_fidelity_msm1.'+self.plot_format),
                format=self.plot_format)
            fig.savefig(os.path.join(self.folder,
                'mean_fidelity_msm1.pdf'),
                format='pdf')

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.errorbar(time, fid0, fmt='.', yerr=fid0_err, label='ms=0')
            ax.errorbar(time, fidp1, fmt='.', yerr=fidp1_err, label='ms=+1')
            ax.errorbar(time, Fp, fmt='.', yerr=Fp_err, label='mean')
            ax.set_xlabel('RO time (us)')
            ax.set_ylabel('RO fidelity')
            ax.set_ylim((0,1))
            ax.legend(loc=4)
            plt.figtext(0.8, 0.5, "max. Fp1=({:.2f} +/- {:.2f})% at t={:.0f} us".format(Fp_max*100., Fp_max_err*100., tp_max),
                    horizontalalignment='right')

            ax.set_title(self.default_plot_title + ': mean RO fidelity')

            fig.savefig(os.path.join(self.folder,
                'mean_fidelity_msp1.'+self.plot_format),
                format=self.plot_format)
            fig.savefig(os.path.join(self.folder,
                'mean_fidelity_msp1.pdf'),
                format='pdf')


        if plot_photon_ms0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            Prob_ms0 = (fid0[1:])/(fid0[1:]+(1-fid1[1:]))
            max_Prob_ms0 = Prob_ms0.max()
            time_max_Prob_ms0 = time[Prob_ms0.argmax()+1]
            ax.errorbar(time[1:], Prob_ms0, fmt='.', yerr=0*Prob_ms0)
            ax.set_xlabel('RO time (us)')
            ax.set_ylabel('Prob. photon came from ms=0')
            ax.set_ylim((0.9,1))
            plt.figtext(0.8, 0.5, "max. {:.2f} at t={:.0f} us".format(max_Prob_ms0*100., time_max_Prob_ms0),
                    horizontalalignment='right')

            ax.set_title(self.default_plot_title + ': Probability_photon_from_ms0')

            fig.savefig(os.path.join(self.folder,
                'projectivity.'+self.plot_format),
                format=self.plot_format)    


def ssrocalib(contains = '',folder='', plot = True, plot_photon_ms0 = True):
    if folder=='' and contains == '':
        folder=toolbox.latest_data('AdwinSSRO')
    elif contains != '':
        folder = toolbox.latest_data(contains)

    a = SSROAnalysis(folder)

    for n,ms in zip(['ms0', 'ms1'], [0,1]): #zip((['ms0'], [0]):#
        a.get_run(n)
        a.cpsh_hist(a.ro_counts, a.reps, name=n, plot = plot)
        a.readout_relaxation(a.ro_time, a.ro_counts, a.reps, a.binsize, name=n, plot = plot)
        a.spinpumping(a.sp_time, a.sp_counts, a.reps, a.binsize, name=n, plot = plot)
        a.charge_hist(a.cr_counts, name=n, plot = plot)
        a.fidelity(a.ro_counts, a.reps, a.binsize, ms, name=n, plot = plot)
    #f = self.analysis_h5data()
    plt.close('all')
    a.mean_fidelity(plot,plot_photon_ms0)
    a.finish()


def ssrocalib_MWInit(folder='', plot = True, plot_photon_ms0 = True):
    if folder=='':
        folder=toolbox.latest_data('_SSRO_calib_MWInit_')
    a = SSROAnalysis(folder)

    for n,ms in zip(['ms0','msp1','msm1'], [0,1,-1]): #zip((['ms0'], [0]):#
        print n, ms
        a.get_run(n)
        a.cpsh_hist(a.ro_counts, a.reps, name=n, plot = plot)
        a.readout_relaxation(a.ro_time, a.ro_counts, a.reps, a.binsize, name=n, plot = plot)
        a.spinpumping(a.sp_time, a.sp_counts, a.reps, a.binsize, name=n, plot = plot)
        a.charge_hist(a.cr_counts, name=n, plot = plot)
        a.fidelity(a.ro_counts, a.reps, a.binsize, ms, name=n, plot = plot)
    #f = self.analysis_h5data()
    plt.close('all')
    a.mean_fidelity_MWInit(plot,plot_photon_ms0)
    a.finish()
    print 'Job\'s done!'

def thcalib(folder='', analyze_probe = False):
    if folder=='':
        folder=toolbox.latest_data('AdwinSSRO')
    a = SSROAnalysis(folder)

    if analyze_probe:
        ths = [1,2,3,4,6,8,10,15,20,25,30]
    else:
        ths = np.linspace(5,45,9)

    for i,th in enumerate(ths):

        if analyze_probe:
            name = 'th_pres_{}_probe_{}'.format(pres[0],th)
        else:
            name = 'th_pres_{}_probe_{}'.format(th,th)

        a = sequence.SequenceAnalysis(folder)
        a.get_cr_results(name, plot=True)
        plt.close('all')

        mean = a.get_mean_cr_cts()
        means.append(mean)

        stats = a.adwingrp(name)['statistics'].value
        fail = stats[2]
        percentage_pass = 5000./(5000. + fail) * 100 #5000 is the number of succesfull measurements.
        percentage_passes.append(percentage_pass)

        a.finish()

    fig = a.default_fig(figsize=(6,4))
    ax = a.default_ax(fig)
    ax.plot(ths, means,'o')
    ax.set_xlabel('sweep threshold')
    ax.set_ylabel('mean CR counts after sequence')
    if save:
        fig.savefig(
            os.path.join(folder, 'post-CR_sum_vs_sweepparam.png'),
            format='png')

    fig = a.default_fig(figsize=(6,4))
    ax = a.default_ax(fig)
    ax.plot(ths[:], percentage_passes[:],'o')
    ax.set_xlabel('sweep threshold')
    ax.set_ylabel('percentage CR passes')
    if save:
        fig.savefig(
            os.path.join(folder, 'percentage_CR_pass_vs_sweepparam.png'),
            format='png')








    for n,th in zip(['th_pres_30.0_probe_30.0', 'th_pres_40.0_probe_40.0',
            'th_pres_50.0_probe_50.0','th_pres_60.0_probe_60.0',
            'th_pres_70.0_probe_70.0','th_pres_80.0_probe_80.0',
            'th_pres_90.0_probe_90.0','th_pres_100.0_probe_100.0',
            'th_pres_110.0_probe_110.0','th_pres_120.0_probe_120.0'], [30,40,50,60,70,80,90,100,110,120]): #zip((['ms0'], [0]):#
        a.get_run(n)
        a.cpsh_hist(a.ro_counts, a.reps, name=n)
        a.readout_relaxation(a.ro_time, a.ro_counts, a.reps, a.binsize, name=n)
        a.spinpumping(a.sp_time, a.sp_counts, a.reps, a.binsize, name=n)
        a.charge_hist(a.cr_counts, name=n)
        a.fidelity(a.ro_counts, a.reps, a.binsize, th, name=n)

    plt.close('all')
    a.finish()


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


def sync_num_fast_SSRO_ph_events(fp_LT, RO_start, VERBOSE = True):
    """
    Returns a list with the sync numbers of events that have at least one photon in the time from the readout start 
    til the readout end. The length of the readout is taken from the data. The function als return the RO counts per shot.
    """

    f = h5py.File(fp_LT, 'r')
    sync_num_RO = f['/PQ_sync_number-1'].value
    special_RO = f['/PQ_special-1'].value
    sync_time_RO = f['/PQ_sync_time-1'].value

    # Get name of the group to find read out length
    group = toolbox.get_msmt_name(fp_LT)
    total_string_name = '/' + group + '/joint_params'
    RO_length = f[total_string_name].attrs['LDE_RO_duration']  * 1e9
    f.close()

    is_ph_RO = special_RO == 0   

    is_in_window = (RO_start  <= sync_time_RO) & (sync_time_RO < (RO_start + RO_length))
    is_ph_RO_in_ro_window = is_in_window & is_ph_RO

    sync_num_ph_events = np.unique(sync_num_RO[is_ph_RO_in_ro_window])

    if VERBOSE:
        print "The total number of events for which a photon was read out is:", len(sync_num_ph_events)
        print "The total number of photons that are readout for these events:", sum(is_ph_RO_in_ro_window)
        print "The average amount of photons detected for each sync number is:", float(sum(is_ph_RO_in_ro_window))\
                                                                                        /len(sync_num_ph_events)


    return sync_num_ph_events, is_ph_RO_in_ro_window, sync_num_RO