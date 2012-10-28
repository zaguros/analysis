import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from analysis.lib.fitting import fit
from analysis import config

config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'
config.datadir = r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/tpqi'

rcParams['mathtext.default'] = 'regular'
rcParams['legend.numpoints'] = 1

class TPQIAnalysis:

    def __init__(self):
        self.savedir = os.path.join(config.outputdir, time.strftime('%Y%m%d')+'-tpqi')
        self.binwidth = 0.256
        self.pulsesep = 600.
        self.fn = 'coincidences_win150.pkl'
        self.deltasinfile = range(-2,3)
        self.deltas = [-1,0,1]
        self.offset = 666.-637.
        self.central = 2
        self.tau = 12.5 / self.binwidth
        self.centernormfactor = 0.5
    
    def get_g2_hist(self, range=(-256,256), bins=64):
        # NOTE bins should be an even number!
        
        f = open(os.path.join(config.datadir,self.fn),'rb')
        self.coincidences = pickle.load(f)
        f.close()

        self.peaks = {}
        self.amplitudes = []
        self.normpeaks = {}

        for i,d in enumerate(self.deltasinfile):
            if d in self.deltas:
                hy,hx = np.histogram(np.array(self.coincidences[i])+self.offset, 
                        bins=bins, range=range)
                self.peaks[d] = (hy,hx)

                if d != 0:
                    A = fit.Parameter(max(hy))

                    def fitfunc(x):
                        return A()*exp(-abs(x)/self.tau)

                    fit.fit1d(hx[:-1]+(hx[1]-hx[0])/2., hy, None, fitfunc=fitfunc, p0=[A], do_print=True, 
                            fixed=[])

                    self.amplitudes.append(A())
                    
        # normalize peaks
        self.meanamp = np.mean(self.amplitudes)
        for d in self.peaks:
            if d != 0:
                self.normpeaks[d] = (self.peaks[d][0]/self.meanamp, self.peaks[d][1])
            else:
                # because we only have two coincidence windows; see labbook
                self.normpeaks[d] = (self.peaks[d][0]/self.meanamp/2., self.peaks[d][1])
                       
        return True   

    
    # NOTE here the deltas are actually hardcoded; fine for now
    def plot_g2(self, fits=True):
        
        self.g2fig,(ax1,ax2,ax3) = plt.subplots(1,3,sharey=True,figsize=(12,8))
              
        for d,ax in zip([-1,0,1], (ax1,ax2,ax3)):
            ax.plot(self.normpeaks[d][1][1:],self.normpeaks[d][0],
                    'k', drawstyle='steps')
            if fits:
                x = np.linspace(-300,300,101)
                if d == 0:
                    ax.plot(x, 0.5*np.exp(-abs(x)/self.tau), 'k', ls='--')
                else:
                    ax.plot(x, np.exp(-abs(x)/self.tau), 'k', ls='--')


        ax1.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        ax1.yaxis.tick_left()
        ax3.yaxis.tick_right()
        ax2.yaxis.set_ticks_position('none')
        ax1.set_ylabel(r'$g^2(\tau)$')
        ax2.set_xlabel(r'$\tau$ (bins)')
        plt.subplots_adjust(wspace=0.15)

    def projected_fidelity(self, crratio=1.1, corr=0.95):
        
        self.fids = np.array([])
        self.u_fids = np.array([])
        
        self.fids_crratio = np.array([])
        self.u_fids_crratio = np.array([])

        self.fids_corr = np.array([])
        self.u_fids_corr = np.array([])
        
        self.dts = np.array([])

        zp = self.peaks[0]
        p1 = self.peaks[1]
        pm1 = self.peaks[-1]
        dt = zp[1][1]-zp[1][0]
        length = len(zp[0])/2
        
        # NOTE formulas see labbook
        for i in range(1,length+1):
            self.dts = np.append(self.dts, i*dt)
            
            z = np.sum(zp[0][length-i:length+i])
            s = np.sum(p1[0][length-i:length+i])+np.sum(pm1[0][length-i:length+i])
            v = (s/4. - z/2.) / (s/4.)
            self.fids = np.append(self.fids, 0.5 + 0.5*v)
            
            dip = 2*crratio/(1+crratio)**2
            v_crratio = (dip*s/2. - z/2.) / (dip*s/2.)
            self.fids_crratio = np.append(self.fids_crratio, 0.5 + 0.5*v_crratio)
            
            u_z = np.sqrt(z)
            u_s = np.sqrt(s)
            u_v = np.sqrt((2./s * u_z)**2 + ((1./s - 4.*(s/4.-z/2.)/s**2) * u_s)**2)
            self.u_fids = np.append(self.u_fids, 0.5 * u_v)

            u_v_crratio = np.sqrt((1./s/dip * u_z)**2 + ((1./s - 2./dip*(s/2.*dip-z/2.)/s**2) * u_s)**2)
            self.u_fids_crratio = np.append(self.u_fids_crratio, 0.5*u_v_crratio)

        self.fids_corr = self.fids_crratio * corr
        self.u_fids_corr = self.u_fids_crratio * corr
        
        self.fidfig = plt.figure()
        ax = plt.subplot(111)
        #ax.errorbar(self.dts, self.fids, yerr=self.u_fids, fmt='o',
        #        label='identical emitters', mfc='w', mec='k', ecolor='k')
        ax.errorbar(self.dts, self.fids_crratio, yerr=self.u_fids_crratio, fmt='o',
                label='real emitters', mfc='w', mec='k', ecolor='k')
        ax.errorbar(self.dts, self.fids_corr, yerr=self.u_fids_corr, fmt='o', 
                label='real emitters, universal corr.', mfc='w', mec='r', ecolor='r')
        
        ax.set_xlabel(r'$| \tau |$ (bins)')
        ax.set_ylabel('projected fidelity')
        plt.legend()


    def save(self):
        if not os.path.isdir(self.savedir):
            os.makedirs(self.savedir)

        try:       
            self.g2fig.savefig(os.path.join(self.savedir,'g2.png'))
        except AttributeError:
            print 'no g2 plot, skip...'

        try:       
            self.fidfig.savefig(os.path.join(self.savedir,'fidelity.png'))
        except AttributeError:
            print 'no fidelity plot, skip...'


        f = open(os.path.join(self.savedir, 'peak_histograms.pkl'), 'wb')
        pickle.dump(self.peaks, f)
        f.close()

        f = open(os.path.join(self.savedir, 'peak_histograms_normalized.pkl'), 'wb')
        pickle.dump(self.normpeaks, f)
        f.close()

        np.savez(os.path.join(self.savedir, 'fidelities'), fid=self.fids_crratio, u_fid=self.u_fids_crratio,
                fid_corr=self.fids_corr, u_fid_corr=self.u_fids_corr, dts=self.dts)

    
       
if __name__ == '__main__':
    a = TPQIAnalysis()
    a.get_g2_hist()
    a.plot_g2()
    a.projected_fidelity()
    a.save()
