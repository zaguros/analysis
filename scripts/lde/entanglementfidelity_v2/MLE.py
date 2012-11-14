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
        self.state = 'psi2'

        # readout fidelities
        self.F0_1 = 0.923
        self.F1_1 = 0.996
        self.F0_2 = 0.818
        self.F1_2 = 0.996

        self.prob_pts = 200
        self.result_bins = 100


    def probability_likelihood(self, N, p00,p01,p10):
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
            lk[p00+p01+p10>1] = 0.  
        elif p00+p01+p10 > 1.:
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

    def _ZZlk(self,pp00,pp01,pp10):
        plk = self.probability_likelihood(self.N_ZZ, pp00, pp01, pp10)
        return np.histogram(pp01+pp10, bins=self.result_bins, range=(0,1),
                density=True, weights=plk)[0]

    def _XXlk_psi2(self,pp00,pp01,pp10):
        plk = self.probability_likelihood(self.N_XX, pp00, pp01, pp10)
        return np.histogram(pp01+pp10-pp00-(1-pp00-pp01-pp10), bins=self.result_bins,
                range=(0,1), density=True, weights=plk)[0]    
    
    def likelihood(self):
        t0 = time.time()

        self.pvals = np.linspace(0.,1.,self.prob_pts)
        self.result_bin_edges = np.linspace(0, 1, self.result_bins+1)

        # construct a grid for the probability values
        # TODO need to be able to chop this up for large arrays; alternative:
        # sparse arrays?
        pp00, pp01, pp10 = \
                np.mgrid[0:self.prob_pts,0:self.prob_pts,0:self.prob_pts]/\
                float((self.prob_pts-1))
        
        self.ZZparity_lk = self._ZZlk(pp00,pp01,pp10)
        self.XXterm_lk = self._XXlk_psi2(pp00,pp01,pp10)

        self.duration = time.time() - t0
        print 'calculation took %.1f secs' % self.duration
    
    def analysis(self, xpts=25):
        fitx = np.linspace(0,1,501)
        dx = self.result_bin_edges[1] - self.result_bin_edges[0]
        x = self.result_bin_edges[:-1]+dx/2.

        # fit ZZ parity
        sigma = fit.Parameter(0.03)
        x0 = fit.Parameter(0.9)
        def ff(x):
            return 1./(sigma()*np.sqrt(2*np.pi)) * np.exp(-0.5*(x-x0())**2/sigma()**2)
        fit.fit1d(x, self.ZZparity_lk, None, fitfunc=ff,
                p0=[sigma,x0], do_print=True)

        zz = x[argmax(self.ZZparity_lk)]
        u_zz = sigma()                
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x+dx/2., self.ZZparity_lk, 'k', drawstyle='steps')
        ax.plot(x, self.ZZparity_lk, 'ko')
        ax.plot(fitx, ff(fitx), 'r-')
        ax.set_xlabel('ZZ parity')
        ax.set_ylabel('Likelihood')
        ax.text(0.1, 1, '%.2f $\pm$ %.2f' % (zz, u_zz))

        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)
        fig.savefig(os.path.join(self.savedir, 'MLE_ZZparity.png'))


        ### same for the XX term       
        sigma = fit.Parameter(0.03)
        x0 = fit.Parameter(0.5)
        def ff(x):
            return 1./(sigma()*np.sqrt(2*np.pi)) * np.exp(-0.5*(x-x0())**2/sigma()**2)
        fit.fit1d(x, self.XXterm_lk, None, fitfunc=ff,
                p0=[sigma,x0], do_print=True)

        xx = x[argmax(self.XXterm_lk)]
        u_xx = sigma()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x+dx/2., self.XXterm_lk, 'k', drawstyle='steps')
        ax.plot(x, self.XXterm_lk, 'ko')
        ax.plot(fitx, ff(fitx), 'b-')
        ax.set_xlabel('XX term')
        ax.set_ylabel('Likelihood')
        ax.text(0.1, 1, '%.2f $\pm$ %.2f' % (xx, u_xx))

        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)
        fig.savefig(os.path.join(self.savedir, 'MLE_XXterm.png'))
        
        F = 0.5 * (zz + xx)
        u_F = np.sqrt(0.25*u_zz**2 + 0.25*u_xx**2)

        print 'Fidelity: F = %.2f +/- %.2f' % (F, u_F)
        

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
    emily.save()
    emily.load('20121111-ldefidelity')
    emily.analysis()




        
