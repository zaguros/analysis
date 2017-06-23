"""
Analysis of time-tagged measurements of the purification class.
Specifically this class allows filtering of adwin RO data based on individual time traces.

NK 2016
"""

from analysis.lib.m2.ssro import pqsequence
from analysis.lib.math import readout_correction as roc
import numpy as np
import os,h5py

class purifyPQAnalysis(pqsequence.PQSequenceAnalysis,):
    """
    Combines pq measurements with adwin ssro.
    """
    def __init__(self,folder,**kw):
        # print folder
        pqsequence.PQSequenceAnalysis.__init__(self,folder,**kw)

        if  'adwindata' in self.g.keys():
            self.agrp=self.adwingrp('adwindata')
        elif 'adwinadata' in self.g.keys(): # At some point this got misspelled! Annoying (PH)
            self.agrp=self.adwingrp('adwinadata')
        else:
            print self.g.keys()
            print folder
            raise Exception('Cant find adwin data!')


        if  'joint_params' in self.g.keys():
            self.joint_grp = self.adwingrp('joint_params')


        self.result_corrected = False


    def filter_pq_data_from_adwin_syncs(self, adwin_syncs = None, pq_syncs = None):

        """
        returns a boolean numpy array.
        True if the recored sync is in the adwin adwin_syncs
        False if the recored sync is not in the adwin_syncs (agrp['counted_awg_reps'])

        TODO: Needs to be generalized for longer PQ meausrements with more data sets.
        """

        if adwin_syncs == None:
            adwin_syncs = self.agrp['counted_awg_reps'].value

        if pq_syncs == None:
            pq_syncs = self.pqf['/PQ_sync_number-1'].value

        return np.in1d(pq_syncs,adwin_syncs) ## CHANGE back

    def filter_adwin_data_from_pq_syncs(self,filtered_sn):
        """
        takes the filtered pq syncs as input and returns a boolean array.
        This array serves as filter for the adwin RO results
        """

        adwin_syncs = self.agrp['counted_awg_reps'].value
        # print 'elen', len(filtered_sn)
        insert_pos = np.searchsorted(adwin_syncs,filtered_sn)
        #insert_pos = np.searchsorted(filtered_sn,adwin_syncs)
        return insert_pos, adwin_syncs[insert_pos] 

    def get_adwin_data_from_pq_syncs(self,filtered_syncs):
        """
        takes the filtered pq syncs as input and returns a boolean array.
        This array serves as filter for the adwin RO results
        """

        if np.sum(np.logical_not(np.in1d(filtered_sn,counted_awg_reps))) != 0:
            print 'Connecting pq syncs to adwin data seems to be going wrong!'

        adwin_syncs = self.agrp['counted_awg_reps'].value

        return np.searchsorted(adwin_syncs,sn_lt[filtered_syncs]) #does not return what we want.

    
    def get_sweep_pts(self):

        self.sweep_name = self.g.attrs['sweep_name']
        self.sweep_pts = self.g.attrs['sweep_pts']

    def get_time_filtered_RO_results(self,adwin_filt_bool):
        """
        also gathers the effective amount of tail counts
        """

        self.pts = self.g.attrs['sweep_length']
        results = self.agrp['ssro_results'].value
        self.readouts = self.g.attrs['nr_of_ROsequences']
        results = results[:(self.pts*(len(results)/self.pts))]
        adwin_filt_bool = adwin_filt_bool[:(self.pts*(len(results)/self.pts))]
        
        sorted_results = results.reshape((-1,self.pts,self.readouts))
        sorted_adwin_fltr = adwin_filt_bool.reshape((-1,self.pts,self.readouts))
        
        
        ### can't do any numpy magic anymore. from here it is brutforce 'for-looping'
        ### (all conceived arrays will have different lengths due to temporal filtering. this break most np methods)
        ### although vstack and hstack would probably work...
        
        self.ro0 = np.array([]) ## how many times did we ro ms=0
        self.no_of_ro = np.array([]) ### total number of ROs
        for i in range(self.pts): 
            self.ro0 = np.append(self.ro0,np.sum((sorted_results[:,i,:])[sorted_adwin_fltr[:,i,:]])) 
            self.no_of_ro = np.append(self.no_of_ro,np.sum(sorted_adwin_fltr[:,i,:]))

        self.normalized_ro = np.multiply(self.ro0,1./self.no_of_ro)
        self.u_normalized_ro = np.sqrt(np.multiply(self.normalized_ro*(1.-self.normalized_ro),1./self.no_of_ro))
        
            
            ### this is how to analyze the tail
        self.tail_cts  = self.pts*(self.no_of_ro * 1e4)/(self.agrp['counted_awg_reps'].value[-1])
        self.u_tail_cts = self.pts*np.sqrt(self.no_of_ro)*1e4/(self.agrp['counted_awg_reps'].value[-1])

    def get_e_ROC(self,ssro_calib_folder,trans = None):
        """
        Applies electron RO correction to a filtered data set.
        Should really make the class above also a child of 'mbi' or 'ssro'
        So far implemented in ipython notebooks
        """
        ssro_fids = get_analysed_ssro_calibration(ssro_calib_folder,
            e_transition = trans, 
            readout_time = self.agrp.attrs['E_RO_durations'][0])

        self.p0,self.u_p0 = roc.single_qubit_outcome_with_ROC_from_fraction(self.normalized_ro, self.u_normalized_ro,*ssro_fids)
        self.result_corrected = True


    def plot_results_vs_sweepparam(self, name='', save=True, **kw):
        mode = kw.get('mode', 'ssro_results')
        labels = kw.get('labels', None)
        ret = kw.get('ret', None)
        ylim = kw.get('ylim', (-0.05, 1.05))
        ax = kw.get('ax', None)
        fmt = kw.get('fmt', 'o')
        figsize=kw.get('figsize',(6,4))
        markersize = kw.get('markersize',6)
        capsize = kw.get('capsize',3)
        if labels != None:
            print 'labels', labels
        if not hasattr(self, 'sweep_pts'):
            self.sweep_pts = np.arange(self.pts) + 1
            self.sweep_name = 'sweep parameter'

        if ax == None:
            fig = self.default_fig(figsize=figsize)
            ax = self.default_ax(fig)
        else:
            ax.set_title(self.timestamp+'\n'+self.measurementstring)

        if labels == None:
            if mode != 'correlations':
                labels = ['RO #%d' % i for i in range(self.readouts)]
            else:
                labels = self.correlation_names

        if mode != 'correlations':
            
            if not self.result_corrected:
                print "Loud noises"
                ax.errorbar(self.sweep_pts, self.normalized_ssro, fmt='o-',
                    yerr=self.u_normalized_ssro, label=labels,markersize=markersize,capsize=capsize)
            else:
                ax.errorbar(self.sweep_pts, self.p0, fmt=fmt,
                    yerr=self.u_p0, label=labels,markersize=markersize,capsize=capsize)
                ax.axhspan(0,1,fill=False,ls='dotted')

            if not self.result_corrected:
                ax.set_ylabel('avg (uncorrected) outcome')
            else:
                ax.set_ylabel(r'$F(|0\rangle)$')

        else:
            for i in range(len(self.correlation_names)):
                if not self.result_correlation_corrected:
                    ax.errorbar(self.sweep_pts, self.normalized_correlations[:,i],
                        fmt=fmt, yerr=self.u_normalized_correlations[:,i],markersize=markersize,capsize=capsize,
                        label=labels[i])
                else:
                    ax.errorbar(self.sweep_pts, self.p_correlations[:,i], fmt=fmt,
                        yerr=self.u_p_correlations[:,i], label=labels[i],markersize=markersize,capsize=capsize)
                    ax.axhspan(0,1,fill=False,ls='dotted')

            if not self.result_correlation_corrected:
                ax.set_ylabel('avg (uncorrected) outcome')
            else:
                ax.set_ylabel('Probability')


        ax.set_xlabel(self.sweep_name)
        ax.set_ylim(ylim)

        if self.readouts > 1:
            ax.legend(loc='best')

        if save and ax != None:
            try:
                fig.savefig(
                    os.path.join(self.folder, '{}_vs_sweepparam.'.format(mode) + self.plot_format),
                    format=self.plot_format)
            except:
                print 'Figure has not been saved.'

        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig



def get_analysed_ssro_calibration(folder, readout_time=None, e_transition = None, sweep_index=None):
    """
    returns calibrated fidelities and uncertainties for ms = 0 and the dark state.
    """

    fp = os.path.join(folder, 'analysis.hdf5')
    f = h5py.File(fp, 'r')

    key=None
    for k in f.keys():
        if 'fidelity' in k:
            if sweep_index==None:
                key=k
            elif k=='fidelity-{}'.format(sweep_index):
                key=k
    if key == None:
        print 'No analysis found, correct sweep-index specified?'
    g=f[key]

    if e_transition == None:
        trans = 'ms1'
    else:
        trans = e_transition
        
    times = g['ms0'].value[:,0]
    fids0 = g['ms0'].value
    fids1 = g[trans].value

    if readout_time==None:
        tidx=len(times)-1
    else:
        tidx = np.argmin(abs(times-readout_time))

    f0 = fids0[tidx,1]
    u_f0 = fids0[tidx,2]
    f1 = fids1[tidx,1]
    u_f1 = fids1[tidx,2]

    return f0, u_f0, f1, u_f1