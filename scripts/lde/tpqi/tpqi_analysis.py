import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from analysis import fit

rcParams['mathtext.default'] = 'regular'

class TPQIAnalysis:

    def __init__(self):
        self.binwidth = 0.256
        self.pulsesep = 600.
        self.files = ['coincidences_all.pkl', 'coincidences_winfull.pkl', 'coincidences_win150.pkl']
        self.names = ['full window (incl laser)', 'tail (full)', 'tail (150 bins)']
        self.deltasinfile = range(-2,3)
        self.offset = (666.-637.)*0.256
        self.central = 2


    def all_g2s(self):
        
        for i,fn,n in zip(range(len(self.files)), self.files, self.names):

            f = open(fn,'rb')
            coincidences = pickle.load(f)
            f.close()

            fig = plt.figure()
            ax = subplot(111)
            ymax = 0

            # histograms
            for i,d in enumerate(self.deltasinfile):
                h = [0] 
                
                if len(coincidences[i] > 0):
                    h,bins,patches = plt.hist(
                            np.array(coincidences[i])*self.binwidth+d*600+self.offset,
                            histtype='step',
                            color='k', hatch='//', 
                            range=(-50+d*600,50+d*600), bins=101)                

                if max(h) > ymax:
                    ymax = max(h)

            # labels
            plt.xlabel('delay (ns)')
            plt.ylabel('events')
            plt.title(n)
            plt.ylim((0,ymax*1.1))

            for i,d in enumerate(self.deltasinfile):
                if len(coincidences[i] > 0):
                    y = ymax*0.95
                    plt.text(d*600, y, str(len(coincidences[i])), ha='center', va='top')

            fig.savefig(n+'.png')

            #plt.xlim((-50,50))
            #fig.savefig(n+' central peak.png')

    def g2(self, fidx=2, normcenter=0.5, tau=12.5):
        
        fn = self.files[fidx]
        f = open(fn,'rb')
        coincidences = pickle.load(f)
        f.close()

        fig,(ax1,ax2,ax3) = plt.subplots(1,3,sharey=True,figsize=(12,8))

        peaks = {
                -1: {
                    'hy' : np.array([]),
                    'hx' : np.array([]),
                    'A' : 0.,
                    'x' : np.array([]),
                    'ax' : ax1,
                    'norm' : 1.,
                    'fit' : np.array([]),
                    },
                0: {
                    'hy' : np.array([]),
                    'hx' : np.array([]),
                    'x' : np.array([]),
                    'A' : 0.,
                    'dw' : 0.,
                    'ax' : ax2,
                    'norm' : normcenter,
                    'fit' : np.array([]),
                    },
                1: {
                    'hy' : np.array([]),
                    'hx' : np.array([]),
                    'x' : np.array([]),
                    'A' : 0.,
                    'ax' : ax3,
                    'norm' : 1.,
                    'fit' : np.array([]),
                    },
                }
        
        for i,d in enumerate(self.deltasinfile):
            if d in peaks:
                hy,hx = np.histogram(np.array(coincidences[i])*self.binwidth+self.offset, 
                        bins=391, range=(-50,50))

                peaks[d]['hy'] = hy
                peaks[d]['hx'] = hx+d*600
                peaks[d]['x'] = hx[:-1]+d*600+(hx[1]-hx[0])*0.5


            if d in [-1,1]:
                A = fit.Parameter(max(hy))
                x0 = fit.Parameter(-0.3)

                def fitfunc(x):
                    return A()*exp(-abs(x-x0())/tau)

                fit.fit1d(peaks[d]['x']-d*600, hy, None, fitfunc=fitfunc, p0=[A,x0], do_print=True, 
                        fixed=[1])
                peaks[d]['A'] = A()
                peaks[d]['fit'] = fitfunc(peaks[d]['x']-d*600)

            if d == 0:
                A = fit.Parameter(5.6*10./9., 'A')
                dw = fit.Parameter(5e6/1e9, 'dw')
                x0 = fit.Parameter(-0.3, 'x0')
                B = fit.Parameter(1., 'B')
                def fitfunc(x):
                    return A() * exp(-abs(x-x0())/tau) * (1-B()*np.exp(-0.25 * dw()**2 * (x-x0())**2))

                ret = fit.fit1d(peaks[d]['x'], hy, None, fitfunc=fitfunc, p0=[A,dw,x0,B], 
                        do_print=True, fixed=[2], ret=True)

                peaks[d]['A'] = A()
                peaks[d]['dw'] = dw()
                peaks[d]['fit'] = fitfunc(peaks[d]['x'])
                peaks[d]['dw_err'] = ret['error_dict']['dw']

        A_mean = 0.5*(peaks[-1]['A']+peaks[1]['A'])
        print 'mean amplitude', A_mean
        for p in peaks:
            pk = peaks[p]
            ax = peaks[p]['ax']
            ax.plot(pk['hx'][:-1], pk['hy']/pk['A']*pk['norm'], 'k', drawstyle='steps')
            if p == 0:
                ax.plot(pk['x'], 0.5*np.exp(-abs(pk['x']/tau)), 'r--', lw=3)
                ax.plot(pk['x'], pk['fit']/pk['A']*pk['norm'], 'r', lw=3)
            else:
                ax.plot(pk['x'], pk['fit']/pk['A']*pk['norm'], 'r', lw=3)
            ax.set_xlim(p*600-50, p*600+50)
            pk['ax'] = ''


        ax1.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        ax1.yaxis.tick_left()
        ax3.yaxis.tick_right()
        ax2.yaxis.set_ticks_position('none')
        ax1.set_ylabel(r'$g^2(\tau)$')
        ax2.set_xlabel(r'$\tau$ (ns)')
        plt.subplots_adjust(wspace=0.15)

        plt.suptitle('TPQI (150 bins)')
        ax2.text(0,1.5, r'$\delta\omega = 2\pi \times$ (%.1f $\pm$ %.1f) MHz' % \
                (peaks[0]['dw']*1e3/2./np.pi, peaks[0]['dw_err']*1e3/2./np.pi), va='top', ha='center')
        # ax2.text(0,1.4, r'($2\pi \times 38$ MHz for fixed Amplitude)', va='top', ha='center')

        f = open('peaks.pkl', 'wb')
        pickle.dump(peaks, f)
        f.close()
        
if __name__ == '__main__':
    a = TPQIAnalysis()
    a.g2()
    # a.all_g2s()
