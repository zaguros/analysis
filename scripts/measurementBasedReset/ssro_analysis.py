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

# analysis classes and shortcut functions
class SSROAnalysis(m2.M2Analysis):

    def fidelity(self, plot=True, plot_photon_ms0=True, **kw):

        counts_ms0 = np.array(self.g['ms0']['RO_data'])
        counts_ms1 = np.array(self.g['ms1']['RO_data'])
        time = np.array(self.g.attrs['sweep_pts'])
        reps = self.g.attrs['repetitions']
        readout_power = self.g.attrs['E_RO_amplitudes_AWG']*1e9

        if np.size(readout_power) > 1:
            print 'this script isnt designed to sweep power!'
            return

        pzero_ms0 = counts_ms0 / float(reps)
        pzero_ms1 = counts_ms1 / float(reps)

        fid0 = pzero_ms0 
        fid1 = 1-pzero_ms1 

        fid0_err = np.sqrt(pzero_ms0*(1-pzero_ms0)/reps)
        fid1_err = np.sqrt(pzero_ms1*(1-pzero_ms1)/reps)
        
        F = (fid0 + fid1)/2.
        F_err = np.sqrt(fid0_err**2 + fid1_err**2)
        maxidx = F.argmax()
        F_max = F[maxidx]
        F_max_err = F_err[maxidx]
        t_max = time[F.argmax()]

        print 'max. F = ({:.2f} +/- {:.2f})% at t={:.0f} us'.format(F_max*100., F_max_err*100., t_max)
        print '\tms_0 = ({:.2f} +/- {:.2f})%'.format(fid0[maxidx]*100, fid0_err[maxidx]*100)
        print '\tms_1 = ({:.2f} +/- {:.2f})%'.format(fid1[maxidx]*100, fid1_err[maxidx]*100)

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

            ax.set_title(self.default_plot_title + '\n mean RO fidelity at '+ str(readout_power) + 'nW')

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


def ssrocalib(contains = '',folder='', plot = True, plot_photon_ms0 = False):
    if folder=='' and contains == '':
        folder=toolbox.latest_data('FastSSROCalib_')
    elif contains != '':
        folder = toolbox.latest_data(contains)

    a = SSROAnalysis(folder)
    a.fidelity(plot,plot_photon_ms0)
    a.finish()

