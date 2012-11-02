import os, sys, time
import pickle
import pprint
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import rcParams
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import ImageGrid

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

# indices used for bases
ZZidx = 0
XXidx = 1
XmXidx = 2

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

    def correlations(self, plot=True, save=False, **kw):
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

        self.readout_correction(**kw)

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

    def filter_times(self, mintimes=(639,668), window=150, dtmin=0, dtmax=50, 
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
            F0_2=0.809, uF0_2=0.007, F1_2=0.985, uF1_2=0.011):
        
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
        self.savedir = os.path.join(config.outputdir, time.strftime('%Y%m%d')+'-ldefidelity')
        self.name = name
       
        self.psi1_ZZ = None
        self.psi1_XX = None
        self.psi1_XmX = None
        self.psi2_ZZ = None
        self.psi2_XX = None
        self.psi2_XmX = None
        
        #self.psi1_XmX_pt2 = None
        #self.psi2_XmX_pt2 = None

        self.F0_1 = 0.917; self.u_F0_1 = 0.003
        self.F1_1 = 0.993; self.u_F1_1 = 0.001
        self.F0_2 = 0.808; self.u_F0_2 = 0.007
        self.F1_2 = 0.985; self.u_F1_2 = 0.011
        
        self.rofids_regular = [self.F0_1, self.F0_2, self.F1_1, self.F1_2, 
                self.u_F0_1, self.u_F0_2, self.u_F1_1, self.u_F1_2]

        self.dtvals = np.arange(21).astype(int)*5 + 10
        self.winvals = np.arange(21).astype(int)*10 + 20
        self.ch0starts = np.arange(6).astype(int)*2 + 631
        self.ch1starts = np.arange(6).astype(int)*2 + 660
        
    def _get_fidelity(self, zzhist, xxhist, xmxhist, state='psi1'):
        rofids_regular = []
        
        ZZ = fidelities.ro_c_F_odd(zzhist[::-1], *self.rofids_regular)
        u_ZZ = fidelities.ro_c_dF(zzhist[::-1], *self.rofids_regular)
        ZZS = fidelities.ro_c_F_S(zzhist[::-1], *self.rofids_regular)
        u_ZZS = fidelities.ro_c_dF_S(zzhist[::-1], *self.rofids_regular)
        
        xxfunc = fidelities.ro_c_F_even if state=='psi1' else fidelities.ro_c_F_odd
        xmxfunc = fidelities.ro_c_F_odd if state=='psi1' else fidelities.ro_c_F_even

        XX = xxfunc(xxhist[::-1], *self.rofids_regular)
        u_XX = fidelities.ro_c_dF(xxhist[::-1], *self.rofids_regular)
        XmX = xmxfunc(xmxhist[::-1], *self.rofids_regular)
        u_XmX = fidelities.ro_c_dF(xmxhist[::-1], *self.rofids_regular)

        nxx = float(sum(xxhist))/sum(xxhist+xmxhist)
        nxmx = float(sum(xmxhist))/sum(xxhist+xmxhist)

        XXavg = nxx*XX + nxmx*XmX
        u_XXavg = np.sqrt((nxx*u_XX)**2 + (nxmx*u_XmX)**2)

        F = ZZ/2. - ZZS + (XXavg - .5)
        u_F = np.sqrt(.25*u_ZZ**2 + u_ZZS**2 + u_XXavg**2)

        return F, u_F

    def _get_sqrtterm(self, zzhist):
        return fidelities.ro_c_F_S(zzhist[::-1], *self.rofids_regular)

    
    def fidelity(self):       
        self.F_psi1, self.u_F_psi1 = self._get_fidelity(
                self.psi1_ZZ, self.psi1_XX, self.psi1_XmX, state='psi1')
        
        self.F_psi2, self.u_F_psi2 = self._get_fidelity(
                self.psi2_ZZ, self.psi2_XX, self.psi2_XmX, state='psi2')


    def get_fidelities(self, folder='20121030-lde', ZZfn='ZZ.npz', XXfn='XX.npz',
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
        
        self.psi1fids = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))
        self.u_psi1fids = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))
        self.psi2fids = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))
        self.u_psi2fids = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))

        self.psi1sqrtterm = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))
        self.psi2sqrtterm = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))

        self.psi1_ZZ_npq = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))
        self.psi1_XX_npq = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))
        self.psi1_XmX_npq = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))

        self.psi2_ZZ_npq = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))
        self.psi2_XX_npq = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))
        self.psi2_XmX_npq = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))
        
        self.psi1_N = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))
        self.psi2_N = np.zeros((len(self.dtvals), len(self.winvals), 
            len(self.ch0starts), len(self.ch1starts)))

        self.rawpsi1correlations = np.zeros((len(self.dtvals), 
            len(self.winvals), len(self.ch0starts), len(self.ch1starts), 3, 4))
        self.rawpsi2correlations = np.zeros((len(self.dtvals), 
            len(self.winvals), len(self.ch0starts), len(self.ch1starts), 3, 4))
        self.correctedpsi1correlations = np.zeros((len(self.dtvals), 
            len(self.winvals), len(self.ch0starts), len(self.ch1starts), 3, 4))
        self.correctedpsi2correlations = np.zeros((len(self.dtvals), 
            len(self.winvals), len(self.ch0starts), len(self.ch1starts), 3, 4))
        self.u_correctedpsi1correlations = np.zeros((len(self.dtvals), 
            len(self.winvals), len(self.ch0starts), len(self.ch1starts), 3, 4))
        self.u_correctedpsi2correlations = np.zeros((len(self.dtvals), 
            len(self.winvals), len(self.ch0starts), len(self.ch1starts), 3, 4))

        cnt = 0
        for i,dt in enumerate(self.dtvals):
            for j,window in enumerate(self.winvals):
                for k,ch0start in enumerate(self.ch0starts):
                    cnt += 1
                    print 'chunk', cnt, '/', \
                            len(self.dtvals)*len(self.winvals)*len(self.ch0starts)

                    for l,ch1start in enumerate(self.ch1starts):
                        
                        
                        eZZ.events = ZZ0.copy()
                        eZZ.filter_times(dtmax=dt, window=window, 
                                mintimes=(ch0start,ch1start))
                        eZZ.correlations(plot=False, save=False)
                        eXX.events = XX0.copy()
                        eXX.filter_times(dtmax=dt, window=window, 
                                mintimes=(ch0start,ch1start))
                        eXX.correlations(plot=False, save=False)
                        eXmX.events = XmX0.copy()
                        eXmX.filter_times(dtmax=dt, window=window, 
                                mintimes=(ch0start,ch1start))
                        eXmX.correlations(plot=False, save=False)
                        
                        self.psi1_ZZ = eZZ.psi1
                        self.psi1_XX = eXX.psi1
                        self.psi1_XmX = eXmX.psi1

                        self.psi2_ZZ = eZZ.psi2
                        self.psi2_XX = eXX.psi2
                        self.psi2_XmX = eXmX.psi2

                        self.rawpsi1correlations[i,j,k,l,ZZidx] = self.psi1_ZZ
                        self.rawpsi1correlations[i,j,k,l,XXidx] = self.psi1_XX
                        self.rawpsi1correlations[i,j,k,l,XmXidx] = self.psi1_XmX
                        self.rawpsi2correlations[i,j,k,l,ZZidx] = self.psi2_ZZ
                        self.rawpsi2correlations[i,j,k,l,XXidx] = self.psi2_XX
                        self.rawpsi2correlations[i,j,k,l,XmXidx] = self.psi2_XmX

                        self.correctedpsi1correlations[i,j,k,l,ZZidx] = eZZ.psi1_corrected
                        self.correctedpsi1correlations[i,j,k,l,XXidx] = eXX.psi1_corrected
                        self.correctedpsi1correlations[i,j,k,l,XmXidx] = eXmX.psi1_corrected 
                        self.correctedpsi2correlations[i,j,k,l,ZZidx] = eZZ.psi2_corrected
                        self.correctedpsi2correlations[i,j,k,l,XXidx] = eXX.psi2_corrected
                        self.correctedpsi2correlations[i,j,k,l,XmXidx] = eXmX.psi2_corrected

                        self.u_correctedpsi1correlations[i,j,k,l,ZZidx] = eZZ.u_psi1_corrected
                        self.u_correctedpsi1correlations[i,j,k,l,XXidx] = eXX.u_psi1_corrected
                        self.u_correctedpsi1correlations[i,j,k,l,XmXidx] = eXmX.u_psi1_corrected 
                        self.u_correctedpsi2correlations[i,j,k,l,ZZidx] = eZZ.u_psi2_corrected
                        self.u_correctedpsi2correlations[i,j,k,l,XXidx] = eXX.u_psi2_corrected
                        self.u_correctedpsi2correlations[i,j,k,l,XmXidx] = eXmX.u_psi2_corrected                    

                        self.fidelity()                
                        self.psi1fids[i,j,k,l] = self.F_psi1
                        self.u_psi1fids[i,j,k,l] = self.u_F_psi1
                        self.psi2fids[i,j,k,l] = self.F_psi2
                        self.u_psi2fids[i,j,k,l] = self.u_F_psi2

                        self.psi1sqrtterm[i,j,k,l] = self._get_sqrtterm(self.psi1_ZZ)
                        self.psi2sqrtterm[i,j,k,l] = self._get_sqrtterm(self.psi2_ZZ)

                        self.psi1_ZZ_npq[i,j,k,l] = \
                                (eZZ.psi1[0]+eZZ.psi1[3])*(eZZ.psi1[1]+eZZ.psi1[2])/float(sum(eZZ.psi1))
                        self.psi1_XX_npq[i,j,k,l] = \
                                (eXX.psi1[0]+eXX.psi1[3])*(eXX.psi1[1]+eXX.psi1[2])/float(sum(eXX.psi1))
                        self.psi1_XmX_npq[i,j,k,l] = \
                                (eXmX.psi1[0]+eXmX.psi1[3])*(eXmX.psi1[1]+eXmX.psi1[2])/float(sum(eXmX.psi1))

                        self.psi2_ZZ_npq[i,j,k,l] = \
                                (eZZ.psi2[0]+eZZ.psi2[3])*(eZZ.psi2[1]+eZZ.psi2[2])/float(sum(eZZ.psi2))
                        self.psi2_XX_npq[i,j,k,l] = \
                                (eXX.psi2[0]+eXX.psi2[3])*(eXX.psi2[1]+eXX.psi2[2])/float(sum(eXX.psi2))
                        self.psi2_XmX_npq[i,j,k,l] = \
                                (eXmX.psi2[0]+eXmX.psi2[3])*(eXmX.psi2[1]+eXmX.psi2[2])//float(sum(eXmX.psi2))

                        self.psi1_N = sum(eZZ.psi1) + sum(eXX.psi1) + sum(eXmX.psi1)
                        self.psi2_N = sum(eZZ.psi2) + sum(eXX.psi2) + sum(eXmX.psi2)
                    
    def save_fidelities(self):
        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)
        
        np.savez(os.path.join(self.savedir, 'fidelities'), psi1fids=self.psi1fids,
                u_psi1fids=self.u_psi1fids, psi2fids=self.psi2fids,
                u_psi2fids=self.u_psi2fids, dtvals=self.dtvals, winvals=self.winvals,
                ch0starts=self.ch0starts, ch1starts=self.ch1starts, psi1_N=self.psi1_N, psi2_N=self.psi2_N,
                psi1_ZZ_npq=self.psi1_ZZ_npq, psi1_XX_npq=self.psi1_XX_npq, psi1_XmX_npq=self.psi1_XmX_npq,
                psi2_ZZ_npq=self.psi2_ZZ_npq, psi2_XX_npq=self.psi2_XX_npq, psi2_XmX_npq=self.psi2_XmX_npq )

        np.savez(os.path.join(self.savedir, 'correlations'), 
                rawpsi1correlations = self.rawpsi1correlations,
                rawpsi2correlations = self.rawpsi2correlations,
                correctedpsi1correlations = self.correctedpsi1correlations,
                correctedpsi2correlations = self.correctedpsi2correlations,
                u_correctedpsi1correlations = self.u_correctedpsi1correlations,
                u_correctedpsi2correlations = self.u_correctedpsi2correlations)

    def load_fidelities(self, folder, fns=['fidelities.npz', 'correlations.npz']):
        for fn in fns:
            f = np.load(os.path.join(config.outputdir, folder, fn))
            for k in f.keys():
                setattr(self, k, f[k])
            f.close()

    
    def plot_fidelity_map(self, ims, xticks, yticks, xlabel='x', ylabel='y'):
        
        fig = plt.figure(figsize=(8,8))
        grid = ImageGrid(fig, 111, # similar to subplot(111)
                nrows_ncols = (2,2),
                axes_pad = 0.75,
                # add_all=True,
                label_mode = "1",
                cbar_mode = 'each',
                cbar_size="5%",
                cbar_pad="1%",
                )
        
        vmins = [0.5, 3., 0.5, 3.]
        titles = ['F (psi1)', 'sigmas (psi1)', 'F (psi2)', 'sigmas (psi2)']

        for i,im in enumerate(ims):
            img = grid[i].imshow(im, cmap=cm.gist_earth, origin='lower', 
                    interpolation='nearest', vmin=vmins[i])
            fig.colorbar(img, cax=grid.cbar_axes[i])
            grid[i].set_title(titles[i])

        grid[2].set_xlabel(xlabel)
        grid[2].set_ylabel(ylabel)
        xt = grid[2].get_xticks().astype(int)[:-1]
        grid[2].set_xticklabels(xticks[xt])
        yt = grid[2].get_yticks().astype(int)[:-1]
        grid[2].set_yticklabels(yticks[yt])

        return fig

    def plot_map_starts(self, save=True, psi1dt=20, psi1win=70, psi2dt=40, psi2win=200):       
        psi1dtidx = argmin(abs(self.dtvals - psi1dt))
        psi1winidx = argmin(abs(self.winvals - psi1win))
        psi2dtidx = argmin(abs(self.dtvals - psi2dt))
        psi2winidx = argmin(abs(self.winvals - psi2win))

        im0 = self.psi1fids[psi1dtidx,psi1winidx,:,:]
        im1 = (self.psi1fids[psi1dtidx,psi1winidx,:,:]-0.5)/self.u_psi1fids[psi1dtidx,psi1winidx,:,:]
        im2 = self.psi2fids[psi2dtidx,psi2winidx,:,:]
        im3 = (self.psi2fids[psi2dtidx,psi2winidx,:,:]-0.5)/self.u_psi2fids[psi2dtidx,psi2winidx,:,:]

        xticks = self.ch1starts
        yticks = self.ch0starts

        fig = self.plot_fidelity_map([im0,im1,im2,im3], xticks, yticks, 'ch1 start', 'ch0 start')
        if save:
            if not os.path.exists(self.savedir):
                os.makedirs(self.savedir)
            fig.savefig(os.path.join(self.savedir, 'fidelities_vs_ch0_vs_ch1.png'))

    def plot_map_window(self, save=True, ch0start=637, ch1start=666):
        ch0idx = argmin(abs(self.ch0starts - ch0start))
        ch1idx = argmin(abs(self.ch1starts - ch1start))

        im0 = self.psi1fids[:,:,ch0idx,ch1idx]
        im1 = (self.psi1fids[:,:,ch0idx,ch1idx]-0.5)/self.u_psi1fids[:,:,ch0idx,ch1idx]
        im2 = self.psi2fids[:,:,ch0idx,ch1idx]
        im3 = (self.psi2fids[:,:,ch0idx,ch1idx]-0.5)/self.u_psi2fids[:,:,ch0idx,ch1idx]

        xticks = self.winvals
        yticks = self.dtvals

        fig = self.plot_fidelity_map([im0,im1,im2,im3], xticks, yticks, 'window', 'dt')
        if save:
            if not os.path.exists(self.savedir):
                os.makedirs(self.savedir)
            fig.savefig(os.path.join(self.savedir, 'fidelities_vs_dt_vs_win.png'))

    



if __name__ == '__main__':
    #### get all fidelities
    fid = FidelityAnalysis('Fidelity')
    fid.get_fidelities()
    fid.save_fidelities()
    #fid.load_fidelities('20121031-ldefidelity')
    #fid.plot_map_starts()
    #fid.plot_map_window()


    #### use this way to extract (and filter) entanglement events from the hhp-data
    # e = EntanglementEventAnalysis('X-X')
    # e.get_entanglement_events('X-X')
    # e.save_events('X-X')

    # e = EntanglementEventAnalysis('ZZ')
    # e.get_entanglement_events('ZZ')
    # e.save_events('ZZ')

    # e = EntanglementEventAnalysis('XX')
    # e.get_entanglement_events('XX')
    # e.save_events('XX')

    #### filtering and correlations
    #e.load_events('20121029-lde', 'X-X_highdc.npz')
    #e.filter_gate_phase()
    #e.filter_times()
    #e.correlations(save=False, F1_2=0.959, uF1_2=0.004)
