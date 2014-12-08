import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)

"""
Analyze the zeno data
NK 2014
"""

def Zeno_get_2Q_values(	timestamp=None, folder=None,
						measurement_name = ['adwindata'], 
						ssro_calib_timestamp =None):
	"""
	Returns the relevant 2qubit values for a given timestamp.
	"""

	### SSRO calibration
	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSRO')
	else:
	    ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
	    ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'

	### Obtain and analyze data
		### postive RO data
	if timestamps[0] == None: 
		folder_a = toolbox.latest_data(contains='positive')
	else:		
		folder_a = toolbox.data_from_time(timestamps[0])
		
def Zeno_2Q_state_fidelity():
	"""
	Plots the state fidelity for a 2-qubit state as a function of time
	"""

def Zeno_2Q_proc_fidelity():
	"""
	Plots the process fidelity for a 2-qubit state as a function of time
	"""

def Zeno_1Q_state_fidelity():
	"""
	Plots the state fidelity for a decoded qubit as a function of time (one parity expectation value)
	"""

def Zeno_1Q_proc_fidelity():
	"""
	Plots the process fidelity for a decoded qubit as a function of time
	"""