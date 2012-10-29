import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
import matplotlib.cm as cm
from analysis.lib.fitting import fit
from analysis.lib.lde import sscorr, fidelities
from analysis import config

config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'
config.datadir = r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde'

rcParams['mathtext.default'] = 'regular'
rcParams['legend.numpoints'] = 1


hhpfilebase = 'hhp_data'
chmaxtime = 2300

# columns of the data field
T1 = 0
T2 = 1
CH1 = 2
CH2 = 3
RO1 = 4
RO2 = 5
CR1 = 6
CR2 = 7
GP = 8

class EntanglementEventAnalysis:

    def __init__(self, name=''):
        self.events = np.zeros((0,9))
        self.savedir = os.path.join(config.outputdir, time.strftime('%Y%m%d')+'-lde')
        self.name = name

    def add_events(self, fp, start=600, stop=1000):
        idx = os.path.splitext(fp)[0][-3:]
        hhpf = np.load(fp)
        hhp = hhpf['hhp']
        hhpf.close()

        folder,fn = os.path.split(fp)
        ssro1, ssro2, cr1, cr2, gate = sscorr.load_adwin_data(folder, DATIDX=int(idx))

        photons = np.logical_and(np.logical_and(hhp[:,3] == 0, hhp[:,1] >= start), hhp[:,1] <= stop)
        mrkrs = hhp[:,3] == 1
        for _i, nsync in np.ndenumerate(hhp[hhp[:,3]==1,0]):
            i = _i[0]
            w1 = hhp[np.logical_and(photons, hhp[:,0] == nsync-1)]
            w2 = hhp[np.logical_and(photons, hhp[:,0] == nsync)]           
            
            if len(w1) == 1 and len(w2) == 1:
                self.events = np.vstack((self.events, 
                    np.array([w1[0,1], w2[0,1], w1[0,2], w2[0,2], 
                        ssro1[i], ssro2[i], cr1[i], cr2[i], gate[i]])))

    def save_events(self, fn):
        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)
        np.savez(os.path.join(self.savedir,fn), events=self.events)

    def load_events(self, folder, fn):
        f = np.load(os.path.join(config.outputdir, folder, fn))
        self.events = f['events']
        f.close()

    def get_entanglement_events(self, folder):
        for (path, dirs, files) in os.walk(os.path.join(config.datadir,folder)):
            for fn in files:
                if hhpfilebase in fn and fn[0] != '.':
                    self.add_events(os.path.join(path, fn))

        print 'found %d events in total' % len(self.events)

    def correlations(self, plot=True, save=False):
        """
        gives correlation counts in the form [ms=00, 01, 10, 11]
        """
        self.psi1= np.zeros(4, dtype=int)
        self.psi2= np.zeros(4, dtype=int)
        for e in self.events:
            idx = int(e[RO1] == 0)*2**1 + int(e[RO2] == 0)*2**0
            if e[CH1] == e[CH2]:
                self.psi1[idx] += 1
            else:
                self.psi2[idx] += 1

        self.readout_correction()

        if save:
            if not os.path.exists(self.savedir):
                os.makedirs(self.savedir)
            np.savez(os.path.join(self.savedir, self.name+'_correlations'), 
                    psi1=self.psi1, psi2=self.psi2,
                    psi1_corrected=self.psi1_corrected, u_psi1_corrected=self.u_psi1_corrected,
                    psi2_corrected=self.psi2_corrected, u_psi2_corrected=self.u_psi2_corrected)
        
        if plot:
            self.plot_correlations()
        
    def plot_correlations(self, save=False):
        fig,[[ax1,ax2],[ax3,ax4]] = plt.subplots(2,2,figsize=(8,8))

        for ax,state,name in zip([ax1,ax2], [self.psi1, self.psi2], 
                ['Psi1', 'Psi2']):
            ax.bar(np.arange(4)-0.4, state, color='w', ec='k', width=0.8)
            ax.set_xticks(np.arange(4))
            ax.set_xticklabels(['00', '01', '10', '11'])
            ax.set_xlabel('State')
            ax.set_ylabel('Occurences')
            ax.set_xlim(-0.5,3.5)
            ax.set_title(name)

        for ax,state,err,name in zip([ax3,ax4], 
                [self.psi1_corrected, self.psi2_corrected],
                [self.u_psi1_corrected, self.u_psi2_corrected], 
                ['Psi1 corrected', 'Psi2 corrected']):
           
            ax.bar(np.arange(4)-0.4, state, color='w', ec='k', width=0.8,
                    yerr=err, ecolor='k')
            ax.set_xticks(np.arange(4))
            ax.set_xticklabels(['00', '01', '10', '11'])
            ax.set_xlabel('State')
            ax.set_ylabel('Probability')
            ax.set_xlim(-0.5,3.5)
            ax.set_title(name)
        
        plt.suptitle(self.name)
        plt.tight_layout()
        
        if save:
            if not os.path.exists(self.savedir):
                os.makedirs(self.savedir)
            fig.savefig(os.path.join(self.savedir, self.name+'_correlations.png'))

    def filter_times(self, mintimes=(637,666), window=150, dtmin=0, dtmax=50, 
            offset=(0,19)):
        badidxs = []
        for i,e in enumerate(self.events):
            ch1 = int(e[CH1])
            ch2 = int(e[CH2])
            if e[T1] < mintimes[ch1] or e[T1] > mintimes[ch1]+window:
                badidxs.append(i)
                continue

            if e[T2] < mintimes[ch2] or e[T2] > mintimes[ch2]+window:
                badidxs.append(i)
                continue
            
            dt = abs(e[T2] - offset[ch2] - e[T1] + offset[ch1])
            if dt < dtmin or dt > dtmax:
                badidxs.append(i)
                continue

        self.events = np.delete(self.events, badidxs, 0)

    def filter_gate_phase(self):
        self.events = self.events[self.events[:,8]>0]

    def readout_correction(self, F0_1=0.917, uF0_1=0.003, F1_1=0.993, uF1_1=0.001,
            F0_2=0.808, uF0_2=0.007, F1_2=0.989, uF1_2=0.005):
        
        p1, up1 = sscorr.ssro_correct_twoqubit_state_photon_numbers(
                self.psi1[::-1], F0_1, F0_2, F1_1, F1_2, dF0a=uF0_1, dF1a=uF1_1,
                dF0b=uF0_2, dF1b=uF1_2, return_error_bars=True, verbose=False)
        self.psi1_corrected = p1[::-1].reshape(-1)
        self.u_psi1_corrected = up1[::-1]

        p2,up2 = sscorr.ssro_correct_twoqubit_state_photon_numbers(
                self.psi2[::-1], F0_1, F0_2, F1_1, F1_2, dF0a=uF0_1, dF1a=uF1_1,
                dF0b=uF0_2, dF1b=uF1_2, return_error_bars=True, verbose=False)
        self.psi2_corrected = p2[::-1].reshape(-1)
        self.u_psi2_corrected = up2[::-1]


class FidelityAnalysis:
 
    def __init__(self, name=''):
        self.savedir = os.path.join(config.outputdir, time.strftime('%Y%m%d')+'-lde')
        self.name = name
       
        self.psi1_ZZ = None
        self.psi1_XX = None
        self.psi1_XmX = None
        self.psi2_ZZ = None
        self.psi2_XX = None
        self.psi2_XmX = None

        self.dtvals = np.arange(21)*5 + 10
        self.winvals = np.arange(21)*10 + 20
        

    def fidelity(self, F0_1=0.917, uF0_1=0.003, F1_1=0.993, uF1_1=0.001,
            F0_2=0.808, uF0_2=0.007, F1_2=0.989, uF1_2=0.005):
        self.F_psi1, self.u_F_psi1, _dict = fidelities.get_fidelity(
                self.psi1_ZZ, self.psi1_XX, self.psi1_XmX,
                psi1=True, F0a=F0_1, F0b=F0_2, F1a=F1_1, F1b=F1_2, dF0a=uF0_1, dF1a=uF1_1,
                dF0b=uF0_2, dF1b=uF1_2)
        self.F_psi2, self.u_F_psi2, _dict = fidelities.get_fidelity(
                self.psi2_ZZ, self.psi2_XX, self.psi2_XmX,
                psi1=False, F0a=F0_1, F0b=F0_2, F1a=F1_1, F1b=F1_2, dF0a=uF0_1, dF1a=uF1_1,
                dF0b=uF0_2, dF1b=uF1_2)

    def optimize_fidelity(self, folder='20121028-lde', ZZfn='ZZ.npz', XXfn='XX.npz',
            XmXfn='X-X.npz'):

        eZZ = EntanglementEventAnalysis('ZZ')
        eZZ.load_events(folder, ZZfn)
        eZZ.filter_gate_phase()
        ZZ0 = eZZ.events.copy()

        eXX = EntanglementEventAnalysis('XX')
        eXX.load_events(folder, XXfn)
        eXX.filter_gate_phase()
        XX0 = eXX.events.copy()
        
        eXmX = EntanglementEventAnalysis('X-X')
        eXmX.load_events(folder, XmXfn)
        eXmX.filter_gate_phase()
        XmX0 = eXmX.events.copy()

        self.psi1fids = np.zeros((len(self.dtvals), len(self.winvals)))
        self.u_psi1fids = np.zeros((len(self.dtvals), len(self.winvals)))
        self.psi2fids = np.zeros((len(self.dtvals), len(self.winvals)))
        self.u_psi2fids = np.zeros((len(self.dtvals), len(self.winvals)))

        for i,dt in enumerate(self.dtvals):
            for j,window in enumerate(self.winvals):

                eZZ.events = ZZ0.copy()
                eZZ.filter_times(dtmax=dt, window=window)
                eZZ.correlations(plot=False, save=False)
                eXX.events = XX0.copy()
                eXX.filter_times(dtmax=dt, window=window)
                eXX.correlations(plot=False, save=False)
                eXmX.events = XmX0.copy()
                eXmX.filter_times(dtmax=dt, window=window)
                eXmX.correlations(plot=False, save=False)
                
                self.psi1_ZZ = eZZ.psi1
                self.psi1_XX = eXX.psi1
                self.psi1_XmX = eXmX.psi1

                self.psi2_ZZ = eZZ.psi2
                self.psi2_XX = eXX.psi2
                self.psi2_XmX = eXmX.psi2

                self.fidelity()                
                self.psi1fids[i,j] = self.F_psi1
                self.u_psi1fids[i,j] = self.u_F_psi1
                self.psi2fids[i,j] = self.F_psi2
                self.u_psi2fids[i,j] = self.u_F_psi2

    def plot(self, save=True):
        if save:
            if not os.path.exists(self.savedir):
                os.makedirs(self.savedir)

        for fids,errs,name in zip([self.psi1fids, self.psi2fids], 
                [self.u_psi1fids, self.u_psi2fids], ['psi1', 'psi2']):
            fig = plt.figure()
            plt.imshow(fids[::-1,:], cmap=cm.bone, interpolation='nearest')
            plt.colorbar()
            plt.xlabel('Window length (bins)')
            plt.ylabel('dt max (bins)')
            plt.title(name+' fidelity')
            plt.xticks(arange(len(self.winvals))[::2], self.winvals[::2])
            plt.yticks(arange(len(self.dtvals))[::2], self.dtvals[::2][::-1])
            
            if save:
                fig.savefig(os.path.join(self.savedir, name+'_fidelites.png'))

            fig = plt.figure()
            plt.imshow((fids[::-1,:]-0.5)/errs[::-1,:], cmap=cm.bone, interpolation='nearest')
            plt.colorbar()
            plt.xlabel('Window length (bins)')
            plt.ylabel('dt max (bins)')
            plt.title(name+' sigma above 0.5')
            plt.xticks(arange(len(self.winvals))[::2], self.winvals[::2])
            plt.yticks(arange(len(self.dtvals))[::2], self.dtvals[::2][::-1])

            if save:
                fig.savefig(os.path.join(self.savedir, name+'_sigmas.png'))
        

    def save(self):
        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)
        
        np.savez(os.path.join(self.savedir, 'fidelities'), psi1fids=self.psi1fids,
                u_psi1fids=self.u_psi1fids, psi2fids=self.psi2fids,
                u_psi2fids=self.u_psi2fids, dtvals=self.dtvals, winvals=self.winvals)



if __name__ == '__main__':
    fid = FidelityAnalysis('Fidelity')
    fid.optimize_fidelity()
    fid.plot()
    fid.save()
