import os, sys
import numpy as np
import h5py
import logging

import matplotlib.cm as cm
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot
from analysis.lib.math import error
from analysis.lib.m2.ssro import  sequence, mbi 

import analysis.lib.qec.nuclear_spin_characterisation as SC #used for simulating FP response of spins
import analysis.lib.qec.hyperfine_params as hyperfine_params; reload(hyperfine_params)

reload(toolbox)


def get_data(name_contains, input_timestamp = None, input_id = None, ssro_calib_folder = ''):
	'''
	Retrieves msmt data for a SINGLE loop in SimpleDecoupling, i.e. 'tot' = 1.
	Inputs: 
	- String containing an element of the file name
	- (Optional): timestamp. Default is latest data 

	Returns:
	- 'x': array with sweeped parameter values
	- 'y': array with msmt outcomes (corrected for SSRO fidelity)
	- 'yerr': array with msmt errors
	- 'folder': path of data file
	- (optional)'timestamp': timestamp corresponding to data. Returned by default unless a timestamp is specified
	- (optional) 'ssro_calib_folder': path of the corresponding ssro calibration file. If '', latest calibration is used 
	'''

	if input_timestamp != None:
		folder = toolbox.data_from_time(input_timestamp)
	else:
		timestamp, folder = toolbox.latest_data(name_contains, return_timestamp = True)
	a = mbi.MBIAnalysis(folder)

	print 'DATA FOLDER & SSRO FOLDER:'
	print folder

	#ax = a.plot_results_vs_sweepparam(ret='ax', )


	# Retrieve data & reshape
	if input_id != None:
		name = 'measurement'+ str(input_id)
		adwingrp = a.adwingrp(name = name)
		a.adgrp = adwingrp
		sweep = a.adgrp.attrs['sweep_pts']
		x = sweep.reshape(-1)[:]
		a.get_readout_results(name = name)
		a.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)
		y = a.p0.reshape(-1)[:]
		yerr = a.u_p0.reshape(-1)[:]
	else:
		a.get_sweep_pts()
		a.get_readout_results(name='measurement0')
		a.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)
		x = a.sweep_pts.reshape(-1)[:]
		y = a.p0.reshape(-1)[:]
		yerr = a.u_p0.reshape(-1)[:]


	if input_timestamp != None: 
		return x, y, yerr, folder
	else:
		return x, y, yerr, folder, timestamp

def get_data_multiple_msmts(name_contains, nr_ids, input_timestamp = None, ssro_calib_folder = ''):
	"""
	Retrieves data for MULTIPLE loops in SimpleDecoupling, i.e. 'tot' > 1.
	For description of inputs & outputs: refer to 'get_data' function above.
	Additional input:
	- nr_ids: number of loops performed in SimpleDecoupling  
	"""

	if input_timestamp != None:
		folder = toolbox.data_from_time(input_timestamp)
		a = mbi.MBIAnalysis(folder)
	else:
		timestamp, folder = toolbox.latest_data(contains = name_contains, return_timestamp = True)
		a = mbi.MBIAnalysis(folder)

	print 'DATA FOLDER & SSRO FOLDER:'
	print folder

	# Retrieve data & reshape
	# x, y & yerr data are collected in lists & concatenated into 1 array afterwards
	x = []
	y = []
	yerr = []
	stitching_points = []

	# Retrieve sweep data 
	for ID in range(nr_ids):
		name = 'measurement'+ str(ID)
		print name
		adwingrp = a.adwingrp(name = name)
		a.adgrp = adwingrp
		sweep = a.adgrp.attrs['sweep_pts']
		x.append(sweep.reshape(-1)[:])
		a.get_readout_results(name = name)
		a.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)
		data = a.p0.reshape(-1)[:]
		y.append(data)
		data_err = a.u_p0.reshape(-1)[:]
		yerr.append(data_err)

	# Find stitching points
	for k in range(len(x)-1):
		a = set(x[k]) # List of x values for measurement k
		b = set(x[k+1])
		stitch = list(a.intersection(b))
		stitching_points.append(stitch)

	# print x[6]
	# print 'Minimum signal =', np.amin(y[6])
	# print 'Corresponding tau =', x[6][np.argmin(y[6])]

	# Concatenate
	x = np.hstack(x)
	y = np.hstack(y)
	yerr = np.hstack(yerr)
	stitching_points = np.hstack(stitching_points)

	if input_timestamp != None: 
		return x, y, yerr, folder, stitching_points
	else:
		return x, y, yerr, folder, timestamp, stitching_points

def get_data_multiple_files(data_stamps):
	"""
	Retrieves data from multiple HDF5 files, each of which can consist of multiple loops in SimpleDecoupling.
	Data is combined into one x, y and yerr array.
	Folder of the last data set is returned
	Input:
	- data_stamps: dictionary with timestamps as keys and ['nr_ids', 'ssro_calib_folder'] as entries. 
	NOTE: if ssro_calib_folder = '', the latest ssro calibration file is used
	Output:
	- x, y, yerr as described in "get_data_multiple_msmts" above 
	- folder of last timestamp in 'data_stamps'
	- points where multiple files/msmts within one file are 'stitched', i.e. doubly measured points
	"""

	timestamps = data_stamps.keys()

	x = []
	y = []
	yerr = []
	folders = [] # Currently unused
	all_stitching_points = []

	for nr, time in enumerate(timestamps):
		nr_id = data_stamps[time][0] # Number of measurements contained in one hdf5 file
		ssro_calib_folder = data_stamps[time][1]
		if nr_id == 1:
			tau, signal, d_signal, path = get_data(name_contains = '', input_timestamp = time, ssro_calib_folder = ssro_calib_folder)
		else:
			tau, signal, d_signal, path, stitching_points = get_data_multiple_msmts(name_contains = '', nr_ids = nr_id, input_timestamp = time, ssro_calib_folder = ssro_calib_folder)
		x.append(tau)
		y.append(signal)
		yerr.append(d_signal)
		folders.append(path) # UNUSED --> To do: check for most recent folder
		# NOTE: WILL FAIL IF NR_IDS = 1!
		all_stitching_points.append(stitching_points)

	for k in range(len(x)-1):
		a = set(x[k]) # List of x values for measurement k
		b = set(x[k+1])
		stitch = list(a.intersection(b))
		all_stitching_points.append(stitch)

	# Concatenate data
	x = np.hstack(x)
	y = np.hstack(y)
	yerr = np.hstack(yerr)
	if len(all_stitching_points) != 0:
		all_stitching_points = np.hstack(all_stitching_points)

	return x, y, yerr, path, all_stitching_points

def fit_fingerprints(spin_dict, Bx_list, Bz, N, sweep_parameter, name_contains, nr_ids = 1, timestamp = None, data_stamps = None, ssro_calib_folder = '', **kw ):
	"""
	Overlaps fingerprint DD data with simulations so that spin bath characteristics can be manually fitted.
	Inputs:
	- spin_dict*: dictionary containing dictionaries representing the coupled C-spins.
	Each spin dictionary should contain another dictionary with entries 'par' and 'perp',
	representing the parallel and perpendicular hyperfine component (in Hz) 
	- Bx_list*: list containing values for perpendicular component of static field.
	- Bz: parallel component of static field (in G)
	- N: number of pulses in XY decoupling scheme
	- sweep_parameter: string indicating which simulation should be multiply plotted. 
	Must be identical to the name of one of the asterisk-marked input parameters.
	If 'None': plots all coupled spins for a single Bx value (STILL NEEDS TO BE IMPLEMENTED)
	EXAMPLE 1: if sweep_parameter = 'Bx_list', all separate values in Bx_list are simulated & plotted for all spins in spin_dict
	EXAMPLE 2: if sweep_parameter = 'spin_dict', each spin in spin_dict is separately simulated & plotted for the given Bx
	- name_contains: string containing an element of the file name
	- (Optional) Timestamp: to import data from a specific time. Overrules 'name_contains' 
	- (Optional) data_stamps: dictionary containing info of multiple hdf5 files with data. See function 'get_data_multiple_files' above.
	If provided, then data is imported from corresponding timestamps (i.e. 'name_contains' and 'timestamp' are overruled) 
	- (Optional) ssro_calib_folder: for data from a single file, specify the corresponding ssro_calib_folder. If '', latest calibration is used
	- (Optional) kw: plot formatting & saving parameters. Can be one of the following:
		-- xlim
		-- ylim
		-- add_stitching_pts: if True, adds stitching points to fingerprint  
	"""

	# Check errors
	if sweep_parameter != 'spin_dict' and sweep_parameter != 'Bx_list' and sweep_parameter != None:
		raise Exception('No sweeping parameter selected!')

	# Retrieve measurement data
	if data_stamps != None:
		x, y, yerr, folder, all_stitching_points = get_data_multiple_files(data_stamps)
	elif timestamp != None:
		if nr_ids == 1:
			x, y, yerr, folder = get_data(name_contains, input_timestamp = timestamp, ssro_calib_folder = ssro_calib_folder)
		else: # multiple loops in SimpleDecoupling
			x, y, yerr, folder, all_stitching_points = get_data_multiple_msmts(name_contains, nr_ids, input_timestamp = timestamp, ssro_calib_folder = ssro_calib_folder)
	else:
		if nr_ids == 1: 
			x, y, yerr, folder, timestamp = get_data(name_contains, input_timestamp = timestamp, ssro_calib_folder = ssro_calib_folder)
		else:
			x, y, yerr, folder, timestamp, all_stitching_points = get_data_multiple_msmts(name_contains, nr_ids, input_timestamp = timestamp, ssro_calib_folder = ssro_calib_folder)

	# Plot measurement data

	if 'figsize' in kw:
		figsize = kw.pop('figsize')
		fig = plt.figure(figsize = figsize)
	else:
		fig = plt.figure(figsize=(25,6))
	ax = fig.add_subplot(111)

	plt.errorbar(x,y,yerr = yerr, fmt = 'o',color='RoyalBlue')

	ymin, ymax = ax.get_ylim()
	# Add stitching points as vertical lines (if they exist)
	add_stitching_pts = kw.pop('add_stitching_pts', True)
	if len(all_stitching_points) != 0 and add_stitching_pts:
		plt.vlines(all_stitching_points, ymin, ymax, linestyles = 'dashed', label = 'stitching_points')
		for val in all_stitching_points:
			plt.text(val, 1., str(val))

	# Set general parameters
	spin_names = spin_dict.keys()

	##################
	## SWEEP NR OF SPINS ##
	##################

	if sweep_parameter == 'spin_dict':

		# Check if only 1 value of Bx exists
		if len(Bx_list) > 1:
			raise Exception("You're sweeping spin_dict has more than 1 value for Bx. Cannot do it.")
		else:
			Bx = Bx_list[0]

		for s in range(len(spin_dict)):
			name = spin_names[s]
			label = '%s: par = %s KHz, perp = %s KHz' % (spin_names[s], spin_dict[name]['par'] *1e-3, spin_dict[name]['perp']*1e-3 )

			# Compute hyperfine components
			HF_par = [ spin_dict[name]['par']  - spin_dict[name]['perp'] * Bx/Bz]
			HF_perp = [ spin_dict[name]['perp'] + spin_dict[name]['par'] * Bx/Bz]

			# Perform simulation
			tau_list = x * 1e-6
			Mt = SC.dyn_dec_signal(HF_par, HF_perp, Bz, N, tau_list)
			FP_signal = ((Mt + 1)/2)

			# Plot simulation result
			ax.plot(x, FP_signal[0,:], '.-', lw=.8, label = label)


	##################
	## SWEEP Bx ##
	##################
	elif sweep_parameter == 'Bx_list':

		for b in range(len(Bx_list)):
			Bx = Bx_list[b]
			label = 'Bx = %s G with %s spins' % ( Bx, len(spin_dict) )

			for s in range(len(spin_dict)):
				name = spin_names[s]
				label = '%s: par = %s KHz, perp = %s KHz' % (spin_names[s], spin_dict[name]['par'] *1e-3, spin_dict[name]['perp']*1e-3 )

				# Compute hyperfine components
				HF_par = [ spin_dict[name]['par']  - spin_dict[name]['perp'] * Bx/Bz]
				HF_perp = [ spin_dict[name]['perp'] + spin_dict[name]['par'] * Bx/Bz]

				# Perform simulation
				tau_list = x * 1e-6
				if s == 0:
					Mt = SC.dyn_dec_signal(HF_par, HF_perp, Bz, N, tau_list)
				else:
					Mt = Mt * SC.dyn_dec_signal(HF_par, HF_perp, Bz, N, tau_list)
			
			FP_signal = ((Mt + 1)/2)

			ax.plot(x, FP_signal[0,:], '.-', lw=.8, label = label)

	##################
	## PLOT CURRENT VALUES ##
	##################
	# NOTE: DOESN'T WORK YET! (30-04-2015)

	# elif sweep_parameter == None:
	# 	label = '%s spin(s) '

	# Plot formatting
	if 'xlim' in kw:
		xlim = kw.pop('xlim')
		plt.xlim( [xlim[0], xlim[1]] )
	if 'ylim' in kw:
		ylim = kw.pop('ylim')
		plt.ylim( [ylim[0], ylim[1]] )

	plt.ylabel(r'Signal normalized to 1',fontsize=25)
	plt.xlabel('tau (us)',fontsize=25)
	plt.tick_params(axis='x', labelsize=25)
	plt.tick_params(axis='y', labelsize=25)
	plt.title(timestamp)
	plt.legend(loc = 4)

	# Save fit
	plt.savefig(os.path.join(folder, 'SimpleDecoupling_Fingerprint.png'))

	return 