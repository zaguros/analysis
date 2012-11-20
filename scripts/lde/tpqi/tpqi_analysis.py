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
        self.fn = 'coincidences_winfull.pkl'
        self.deltasinfile = range(-2,3)
        self.deltas = [-1,0,1]
        self.offset = 666.-637.
        self.central = 2
        self.tau = 12.1 / self.binwidth
        self.centernormfactor = 0.5
    
    def get_g2_hist(self, range=(-260,260), bins=26):
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

    def visibility(self, crratio=1.1):
        
        self.visibility = np.array([])
        self.u_visibility = np.array([])
        self.fidelity = np.array([])
        self.u_fidelity = np.array([])
        self.dts = np.array([])

        zp = self.peaks[0]
        amp = 2*self.meanamp # this is also the expected height of the peak
                             # in the middle for coherent light; (only two
                             # windows)
        length = len(zp[0])/2
        dip = 2*crratio/(1+crratio)**2
        dt = zp[1][1]-zp[1][0]

        # NOTE formulas see labbook
        for i in range(0,length):
            self.dts = np.append(self.dts, i*dt)
            tpqicts = np.sum(zp[0][length+i:length+i+1]) + np.sum(zp[0][length-i-1:length-i])
            
            def notpqicts(x0,x1):
                return -2*dip*amp/dt*self.tau*(np.exp(-abs(x1)/self.tau) - np.exp(-abs(x0)/self.tau))
            
            n = notpqicts(i*dt, (i+1)*dt)
            v = (n-tpqicts)/n
            u_v = np.sqrt(tpqicts)/n
            F = 0.5*(1+v)
            u_F = 0.5*u_v

            self.visibility = np.append(self.visibility, v)
            self.u_visibility = np.append(self.u_visibility, u_v)
            self.fidelity = np.append(self.fidelity, F)
            self.u_fidelity = np.append(self.u_fidelity, u_F)
        
        self.fidfig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(np.append(0, self.dts+dt), np.append(self.visibility[0],self.visibility), 
                'k', drawstyle='steps')
        ax.errorbar(self.dts+dt/2., self.visibility, fmt='o', mec='k', mfc='w',
                yerr=self.u_visibility, ecolor='k')

        ax.set_xlabel(r'$| \tau |$ (bins)')
        ax.set_ylabel('visibility')
        ax.set_xlim(0,100)
        ax.set_ylim(-0.5,1)

        ax2 = ax.twinx()
        ax2.plot(np.append(0,self.dts+dt), np.append(self.fidelity[0], self.fidelity),
                'r', drawstyle='steps')
        ax2.errorbar(self.dts+dt/2., self.fidelity, fmt='o', mec='r', mfc='w',
                yerr=self.u_fidelity, ecolor='r')
        ax2.set_ylabel('fidelity', color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')

        ax2.set_xlim(0,100)
        ax2.set_ylim(0,1)

    def save(self):
        if not os.path.isdir(self.savedir):
            os.makedirs(self.savedir)

        try:       
            self.g2fig.savefig(os.path.join(self.savedir,'g2.png'))
        except AttributeError:
            print 'no g2 plot, skip...'

        try:       
            self.fidfig.savefig(os.path.join(self.savedir,'visibility.png'))
        except AttributeError:
            print 'no fidelity plot, skip...'


        f = open(os.path.join(self.savedir, 'peak_histograms.pkl'), 'wb')
        pickle.dump(self.peaks, f)
        f.close()

        f = open(os.path.join(self.savedir, 'peak_histograms_normalized.pkl'), 'wb')
        pickle.dump(self.normpeaks, f)
        f.close()

        np.savez(os.path.join(self.savedir, 'visibility'), 
                visibility=self.visibility, 
                u_visibility=self.u_visibility, 
                dts=self.dts,
                fidelity=self.fidelity,
                u_fidelity=self.u_fidelity)

    
       
if __name__ == '__main__':
    a = TPQIAnalysis()
    a.get_g2_hist()
    a.plot_g2()
    a.visibility()
    a.save()
