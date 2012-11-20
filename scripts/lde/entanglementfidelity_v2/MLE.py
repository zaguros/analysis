import os, sys, time
import pickle
import pprint
import numpy as np
from scipy.misc import factorial as fac
from scipy import integrate

from matplotlib import pyplot as plt
from matplotlib import rcParams

from analysis import config
from analysis.lib.fitting import fit

config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'

rcParams['mathtext.default'] = 'regular'
rcParams['legend.numpoints'] = 1


class Emily:

    def __init__(self):
        self.savedir = os.path.join(config.outputdir, 
                time.strftime('%Y%m%d')+'-ldefidelity') 

        self.N_ZZ = np.array([8,72,66,16])
        self.N_XX = np.array([9,34,37,21])
        self.N_XmX = np.array([30,12,15,34])
        self.state = 'psi2'

        # readout fidelities
        self.F0_1 = 0.923
        self.F1_1 = 0.996
        self.F0_2 = 0.818
        self.F1_2 = 0.996

        self.prob_pts = 101
        self.ab_pts = 101
        self.blocks = 1
        self.result_bins = 50
        self.F_bins = 50


    def diag_likelihood(self, N, p00,p01,p10):
        """
        returns the value of the likelihood function, which is given by
        a multinomial PDF (results being now parameters and the probabilities
        the variables) multiplied by the readout error. 
        p11 is given by 1-p00-p01-p10.
        """
        
        perr00,perr01,perr10,perr11 = self._apply_readout_error(
                p00,p01,p10,1.-p00-p01-p10)

        lk = fac(N.sum()) / \
                (fac(N[0])*fac(N[1])*fac(N[2])*fac(N[3])) * \
                perr00**N[0] * perr01**N[1] * perr10**N[2] * \
                (1.-perr00-perr01-perr10)**N[3]

        if type(p00) == np.ndarray:
            lk[p00+p01+p10>1]=0.; lk[p00<0]=0; lk[p01<0]=0; lk[p10<0]=0
        elif p00+p01+p10 > 1. or p00<0 or p01<0 or p10<0:
            return 0
        return lk

    def offdiag_likelihood(self, N, a, b):
        """
        the probabilities for the XX measurement have some constraints, and are
        in the end described by the parameters a and b (see lab book)
        """
        perr00,perr01,perr10,perr11 = self._apply_readout_error(
                .25+a-b,
                .25-a-b,
                .25-a+b,
                .25+a+b)

        lk = fac(N.sum()) / \
                (fac(N[0])*fac(N[1])*fac(N[2])*fac(N[3])) * \
                perr00**N[0] * perr01**N[1] * perr10**N[2] * \
                (1.-perr00-perr01-perr10)**N[3]
        if type(a) == np.ndarray:
            lk[abs(b) > .25-abs(a)] = 0.
        elif abs(b) > .25-abs(a):
            return 0
        return lk

    def _apply_readout_error(self, p00, p01, p10, p11):
        """
        Uses the readout fidelities given as member variables.
        """
        Err1 = np.array([[self.F0_1, 1.-self.F1_1],
            [1.-self.F0_1, self.F1_1]])
        Err2 = np.array([[self.F0_2, 1.-self.F1_2],
            [1.-self.F0_2, self.F1_2]])
        Err = np.kron(Err1, Err2)

        new = []
        for i in range(4):
            new.append(Err[i,0]*p00 + Err[i,1]*p01 + Err[i,2]*p10 + Err[i,3]*p11)
        return new  
    
    def _ZZlk(self):
        hist = np.zeros(self.result_bins)
        histLB = np.zeros(self.result_bins)
        for i in range(self.blocks):
            pp00,pp01,pp10 = \
                    np.mgrid[i*self.prob_pts/self.blocks:(i+1)*self.prob_pts/self.blocks,\
                    0:self.prob_pts,0:self.prob_pts]/float(self.prob_pts)
            plk = self.diag_likelihood(self.N_ZZ, pp00,pp01,pp10)
            hist += np.histogram(pp01+pp10, bins=self.result_bins,
                    range=(0,1), density=False, weights=plk)[0]
            histLB += np.histogram(pp01+pp10-2*np.sqrt(pp00*(1.-pp00-pp01-pp10)), 
                    bins=self.result_bins, range=(0,1), density=False, weights=plk)[0]
        I = integrate.trapz(hist, dx=self.ZZ_dx)
        hist /= I
        ILB = integrate.trapz(histLB, dx=self.ZZ_dx)
        histLB /= ILB
        return hist, histLB

    def _XXlk(self, XmX=False):
        """
        see labbook; a goes from -1/4 to 1/4. abs(b) <= 1/4-abs(a).
        """

        hist_XX = np.zeros(self.result_bins)
        if XmX: 
            hist_XmX = np.zeros(self.result_bins)
        for i in range(self.blocks):
            aa,bb = np.mgrid[i*self.ab_pts/self.blocks:(i+1)*self.ab_pts/self.blocks,\
                    0:self.ab_pts]/float(2*self.ab_pts)-0.25
            plk_XX = self.offdiag_likelihood(self.N_XX, aa, bb)
            hist_XX += np.histogram(4*aa, bins=self.result_bins,
                range=(-1,1), density=False, weights=plk_XX)[0]
            if XmX:
                plk_XmX = self.offdiag_likelihood(self.N_XmX, aa, bb)
                hist_XmX += np.histogram(-4*aa, bins=self.result_bins,
                        range=(-1,1), density=False, weights=plk_XmX)[0]
        
        I = integrate.trapz(hist_XX,dx=self.XX_dx)
        hist_XX /= I
        
        # we actually need the 
        if XmX:
            I = integrate.trapz(hist_XmX, dx=self.XX_dx)
            hist_XmX /= I

            XXvals, XmXvals = np.meshgrid(self.XX_x, self.XX_x)
            hhXX, hhXmX = np.meshgrid(hist_XX, hist_XmX)
            hist = np.histogram(0.5*(XXvals+XmXvals), bins=self.result_bins,
                    range=(-1,1), density=True, weights=hhXmX*hhXX)[0]
        else:
            hist = hist_XX

        return hist
    
    def likelihood(self):
        t0 = time.time()
        
        self.pvals = np.linspace(0.,1.,self.prob_pts)
        self.ZZ_bin_edges = np.linspace(0,1,self.result_bins+1)
        self.ZZ_dx = self.ZZ_bin_edges[1] - self.ZZ_bin_edges[0]
        self.ZZ_x = self.ZZ_bin_edges[:-1] + self.ZZ_dx/2.
        self.XX_bin_edges = np.linspace(-1,1,self.result_bins+1)
        self.XX_dx = self.XX_bin_edges[1] - self.XX_bin_edges[0]
        self.XX_x = self.XX_bin_edges[:-1] + self.XX_dx/2.
        sign = -1 if self.state == 'psi2' else 1

        self.F_bin_edges = np.linspace(0, 1, self.F_bins+1)
        self.F_dx = self.F_bin_edges[1] - self.F_bin_edges[0]
        self.F_x = self.F_bin_edges[:-1]+self.F_dx/2.
        
        self.ZZparity_lk, self.ZZparityLB_lk = self._ZZlk()
        self.XXterm_lk = self._XXlk(XmX=True)

        # calculate the fidelity likelihoods
        ZZvals, XXvals = np.meshgrid(self.ZZ_x,self.XX_x)
        ZZlks, XXlks = np.meshgrid(self.ZZparity_lk, self.XXterm_lk)
        self.F_lk, _x = np.histogram(0.5*(ZZvals+XXvals*sign), bins=self.F_bins,
                range=(0,1), density=True, weights=ZZlks*XXlks)
        
        ZZLBlks, XXlks = np.meshgrid(self.ZZparityLB_lk, self.XXterm_lk)
        self.FLB_lk, _x = np.histogram(0.5*(ZZvals+XXvals*sign), bins=self.F_bins,
                range=(0,1), density=True, weights=ZZLBlks*XXlks)

        self.duration = time.time() - t0
        print 'calculation took %.1f secs' % self.duration
    
    def analysis(self):
        ZZ_fitx = F_fitx = np.linspace(0,1,501)
        XX_fitx = np.linspace(-1,1,501)
        sign = -1 if self.state == 'psi2' else 1

        # fit ZZ parity
        sigma = fit.Parameter(0.03)
        x0 = fit.Parameter(0.9)
        def ff(x):
            return 1./(sigma()*np.sqrt(2*np.pi)) * np.exp(-0.5*(x-x0())**2/sigma()**2)
        fit.fit1d(self.ZZ_x, self.ZZparity_lk, None, fitfunc=ff,
                p0=[sigma,x0], do_print=True)

        zz = self.ZZ_x[argmax(self.ZZparity_lk)]
        u_zz = sigma()                
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.ZZ_x+self.ZZ_dx/2., self.ZZparity_lk, 'k', drawstyle='steps')
        # ax.plot(self.ZZ_x, self.ZZparity_lk, 'ko')
        ax.plot(ZZ_fitx, ff(ZZ_fitx), 'k-', lw='2')
        ax.set_xlabel('ZZ parity')
        ax.set_ylabel('Likelihood')
        ax.text(0.1, 1, '%.3f $\pm$ %.3f' % (zz, u_zz))

        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)
        fig.savefig(os.path.join(self.savedir, 'MLE_ZZparity.png'))


        ### same for the XX term
        sigma = fit.Parameter(0.03)
        x0 = fit.Parameter(0.5)
        def ff(x):
            return 1./(sigma()*np.sqrt(2*np.pi)) * np.exp(-0.5*(x-x0())**2/sigma()**2)
        fit.fit1d(self.XX_x, self.XXterm_lk, None, fitfunc=ff,
                p0=[sigma,x0], do_print=True)

        xx = x0() # x[argmax(self.XXterm_lk)]
        u_xx = sigma()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.XX_x+self.XX_dx/2., self.XXterm_lk, 'k', drawstyle='steps')
        # ax.plot(self.XX_x, self.XXterm_lk, 'ko')
        ax.plot(XX_fitx, ff(XX_fitx), 'k-', lw='2')
        ax.set_xlabel('XX term')
        ax.set_ylabel('Likelihood')
        ax.text(0.1, 1, '%.3f $\pm$ %.3f' % (xx, u_xx))

        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)
        fig.savefig(os.path.join(self.savedir, 'MLE_XXterm.png'))
        
        ### Fidelity
        sigma = fit.Parameter(0.05)
        x0 = fit.Parameter(0.75)
        def ff(x):
            return 1./(sigma()*np.sqrt(2*np.pi)) * np.exp(-0.5*(x-x0())**2/sigma()**2)
        fit.fit1d(self.F_x, self.F_lk, None, fitfunc=ff,
                p0=[sigma,x0], do_print=True)
        F = x0()
        u_F = sigma()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.F_x+self.F_dx/2., self.F_lk, 'k', drawstyle='steps')
        ax.plot(F_fitx, ff(F_fitx), 'k-', lw='2')
        
        ### Fidelity lower bound
        sigma = fit.Parameter(0.05)
        x0 = fit.Parameter(0.75)
        fit.fit1d(self.F_x, self.FLB_lk, None, fitfunc=ff,
                p0=[sigma,x0], do_print=True)
        FLB = x0()
        u_FLB = sigma()

        ax.plot(self.F_x+self.F_dx/2., self.FLB_lk, 'r', drawstyle='steps')
        ax.plot(F_fitx, ff(F_fitx), 'r-', lw='2')
        
        ax.set_xlabel('Fidelity')
        ax.set_ylabel('Likelihood')
        ax.text(0.1, 1, '$F = %.3f \pm %.3f$' % (F, u_F))
        ax.text(0.1, 0.5, '$F_{LB} = %.3f \pm %.3f$' % (FLB, u_FLB), color='r')


    def save(self):
        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)
        
        np.savez(os.path.join(self.savedir, 'MLEresults_%s' % self.state),
                pvals = self.pvals,
                result_bin_edges = self.result_bin_edges,
                ZZparity_lk = self.ZZparity_lk,
                XXterm_lk = self.XXterm_lk,
                )

    def load(self, folder):
        f = np.load(os.path.join(config.outputdir, folder, 'MLEresults_%s.npz' % self.state))
        for k in f.keys():
            setattr(self, k, f[k])
        f.close()
        
if __name__ == '__main__':
    emily = Emily()
    emily.likelihood()
    # emily.save()
    # emily.load('20121111-ldefidelity')
    emily.analysis()




        
