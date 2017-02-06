"""
Analysis of time-tagged measurements of the purification class.
Specifically this class allows filtering of adwin RO data based on individual time traces.

NK 2016
"""

from analysis.lib.m2.ssro import pqsequence
import numpy as np
import os,h5py

class purifyPQAnalysis(pqsequence.PQSequenceAnalysis):
	"""
	Combines pq measurements with adwin ssro.
	"""
	def __init__(self,folder,**kw):
		# print folder
		pqsequence.PQSequenceAnalysis.__init__(self,folder,**kw)

		self.agrp=self.adwingrp('adwindata')
		self.joint_grp = self.adwingrp('joint_params')

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

		TODO: generalize for arbitrary PQ data size
		"""
		adwin_syncs = self.agrp['counted_awg_reps'].value
		# print 'elen', len(filtered_sn)
		insert_pos = np.searchsorted(adwin_syncs,filtered_sn)
		#insert_pos = np.searchsorted(filtered_sn,adwin_syncs)

		return insert_pos, adwin_syncs[insert_pos] #does not return what we want.

	def get_adwin_data_from_pq_syncs(self,filtered_syncs):
		"""
		takes the filtered pq syncs as input and returns a boolean array.
		This array serves as filter for the adwin RO results

		TODO: generalize for arbitrary PQ data size
		"""
		sn_lt = self.pqf['/PQ_sync_number-1'].value
		adwin_syncs = self.agrp['counted_awg_reps'].value

		return np.searchsorted(adwin_syncs,sn_lt[filtered_syncs]) #does not return what we want.

	def get_filtered_ROC(normalized_ssro,u_normalized_ssro):
		"""
		Applies electron RO correction to a filtered data set.
		Should really make the class above also a child of 'mbi' or 'ssro'
		So far implemented in ipython notebooks
		"""
		pass

	def get_analysed_ssro_calibration(self,ssro_calib_folder):
		"""
		returns calibrated fidelities and uncertainties for ms = 0 and the dark state.
		"""

def get_analysed_fast_ssro_calibration(folder, readout_time=None, e_transition = None, sweep_index=None):

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