"""
module to postprocess and analyze the resulting data
from an adwin SSRO measurement (measurement2 generation).

"""

### imports
import sys, os, time
import numpy as np
import h5py
from matplotlib import pyplot as plt
from measurement.lib.tools import toolbox

FILEBASE = 'AdwinSSRO'

class SSROAnalysis:
    
    def __init__(self, folder):
        self.folder = folder
        self.h5filepath = toolbox.measurement_filename(folder)
        self.f = h5py.File(self.h5filepath, 'r')
        self.name = self.f.keys()[0]
        self.g = self.f[self.name]
        
        self.savebase = os.path.join(self.folder, 'analysis')
        self.analysis_file = h5py.File(self.savebase + '.hdf5', 'w')

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
    
    def analyze_run(self, ms):
        savebase_ms = self.savebase + '_ms' + str(ms)
        
        ### load data
        grp = self.adwingrp('ms'+str(ms))
        reps = grp['completed_reps'].value
        cycle_duration = grp.attrs['cycle_duration']
    
        repetitions = np.arange(reps)
        ssro_time = np.arange(grp.attrs['SSRO_duration']) * cycle_duration * 10./3.
        sp_time = np.arange(grp.attrs['SP_duration']) * cycle_duration * 10./3.

        binsize = grp.attrs['cycle_duration']/300.
        rocounts = grp['RO_data'].value.reshape(reps,grp.attrs['SSRO_duration'])
        crcounts = grp['CR_after'].value
        spcounts = grp['SP_hist'].value

        ### pack; historical reasons
        spdata = {'counts' : spcounts, 'time' : sp_time }
        rodata = {'counts' : rocounts, 'time' : ssro_time, 
                'repetitions' : repetitions }
        crdata = {'repetitions' : repetitions, 'counts' : crcounts }

        ### start analysis of the run
        firstbin = 0
        stop = ssro_time[-1]/1e3
        lastbin = int(stop/binsize)
        fig1 = plt.figure(figsize=(8,8))

        ### histogram for photon counts per shot (cpsh)
        annotation = 'repetitions: %d' % reps
        cpsh = np.sum(rocounts[:,firstbin:lastbin], axis=1) # sum over all columns
        
        ax = plt.subplot(221)
        plt.hist(cpsh, max(cpsh)+1, align='mid', label='counts')
        plt.xlabel('counts/shot')
        plt.ylabel('occurrences')
        plt.title('counts per shot, m_s = %d' % ms)

        mean_cpsh = sum(cpsh)/float(reps)
        annotation += "\nmean c.p.sh. = %.2f" % (mean_cpsh)

        pzero = len(np.where(cpsh==0)[0])/float(reps)
        annotation += "\np(0 counts) = %.2f" % (pzero)

        plt.figtext(0.45, 0.85, annotation, horizontalalignment='right',
                verticalalignment='top')

        ### spin relaxation during readout
        ro_time = ssro_time[:lastbin] # create time axis [us]
        ro_countrate = np.sum(rocounts[:,:lastbin], axis=0) / (binsize*1e-6*reps) # [Hz]

        ax = plt.subplot(222)
        plt.plot(ro_time, ro_countrate, 'o')
        plt.xlabel('RO time [ns]')
        plt.ylabel('counts [Hz]')
        plt.title('spin relaxation during readout')
        
        ### spin pumping data
        # convert counts to Hz
        sp = np.array([j/(binsize*1e-6*reps) for j in spcounts])
        ax = plt.subplot(223)
        plt.plot(spdata['time']+1E3, sp, 'o')
        plt.xlabel('spin pumping time [ns]')
        plt.ylabel('counts [Hz]')
        plt.title('spin relaxation during pumping')

        ### charge readout
        ax = plt.subplot(224)
        cr = crdata['counts']
        plt.hist(cr, abs(max(cr)-min(cr)+1), label='cr')
        #### plt.hist(cr[:,2], abs(max(cr[:,2])-min(cr[:,2]))+1, label='cr2')
        plt.xlabel('counts')
        plt.ylabel('occurences')
        plt.title('charge readout statistics')
        # plt.legend()

        ### fidelity analysis
        fid_dat = np.zeros((0,3))

        for i in range(1,lastbin):
            t = i*binsize
            
            # d: hist of counts, c: counts per shot
            d = np.sum(rocounts[:,:i], axis=1)
            c = np.sum(d)/float(reps)

            # we get the fidelity from the probability to get zero counts in a
            # shot 
            pzero = len(np.where(d==0)[0])/float(reps)
            pzero_err = np.sqrt(pzero*(1-pzero)/reps)
            fid = 1-pzero if ms == 0 else pzero # fidelity calc. depends on ms
            fid_dat = np.vstack((fid_dat, np.array([[t, fid, pzero_err]])))

        fig2 = plt.figure()
        plt.errorbar(fid_dat[:,0], fid_dat[:,1], fmt='o', yerr=fid_dat[:,2])
        plt.xlabel('RO time [us]')
        plt.ylabel('ms = %d RO fidelity' % ms)
        plt.title('SSRO fidelity')

        fig1.savefig(savebase_ms+'_autoanalysis.pdf', format='pdf')
        fig2.savefig(savebase_ms+'_fid_vs_ROtime.pdf', format='pdf')
        
        self.analysis_file['/fidelity/ms'+str(ms)] = fid_dat
        self.analysis_file.flush()

    def mean_fidelity(self):
        
        _fid0 = self.analysis_file['/fidelity/ms0']
        _fid1 = self.analysis_file['/fidelity/ms1']

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

        self.analysis_file['/fidelity/time'] = time
        self.analysis_file['/fidelity/mean_fidelity'] = meanfid
        self.analysis_file['/fidelity/mean_fidelity_err'] = meanfiderr
        self.analysis_file.flush()

        fig = plt.figure()
        plt.errorbar(time, fid0, fmt='.', yerr=fid0_err, label='ms=0')
        plt.errorbar(time, fid1, fmt='.', yerr=fid1_err, label='ms=1')
        plt.errorbar(time, F, fmt='.', yerr=F_err, label='mean')
        plt.xlabel('RO time (us)')
        plt.ylabel('RO fidelity')
        plt.ylim((0.5,1))
        plt.title('SSRO fidelity')
        plt.legend(loc=4)
        plt.figtext(0.8, 0.5, "max. F=%.2f at t=%.2f us" % (F_max, t_max),
                horizontalalignment='right')

        fig.savefig(self.savebase+'_fidelity.pdf', format='pdf')

    def analyze_calibration(self):
        self.analyze_run(0)
        self.analyze_run(1)
        self.mean_fidelity()
        

    def finish(self):
        self.analysis_file.close()

        
def ssro(folder=os.getcwd()):
    a = SSROAnalysis(folder)
    a.analyze_calibration()
    a.finish()


