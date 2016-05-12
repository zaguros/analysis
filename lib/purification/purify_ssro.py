'''
Analyzes all non PQ data of the purification project.
All data are analyzed according to the MBI analysis class with varying input names
2016 NK
'''

import numpy as np
import os
from analysis.lib.tools import toolbox, plot; 
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
# import matplotlib as mpl
from analysis.lib.fitting import fit, common

reload(fit);reload(mbi);reload(common);reload(toolbox)

CR_after_check = True # global variable that let's us post select whether or not the NV was ionized

def get_tstamp_from_folder(folder):
    return folder[18:18+15]

def create_plot(f,**kw):
	ylabel = kw.pop('ylabel',None)
	xlabel = kw.pop('xlabel',None)
	title = kw.pop('title',None)

	fig = plt.figure()
	ax = plt.subplot()

	if xlabel != None:
		plt.xlabel(xlabel)
	if ylabel != None:
		plt.ylabel(ylabel)
	else:
		plt.ylabel('Contrast')
	if title != None:
		plt.title(title + ' '+get_tstamp_from_folder(f))
	else:
		plt.title(get_tstamp_from_folder(f))
    
    
def save_and_close_plot(f):
    plt.savefig(os.path.join(f,'Results.pdf'),format='pdf')
    plt.savefig(os.path.join(f,'Results.png'),format='png')
    plt.show()
    plt.close('all')

def plot_data(x,y,**kw):
    label = kw.get('label',None)
    y_u = kw.pop('y_u',None)
    if y_u != None:
        plt.errorbar(x,y,y_u,fmt = 'o',label = label,**kw)
    else: plt.plot(x,y)


def quadratic_addition(X,Y,X_u,Y_u):
	'''
	takes RO results and averages them according to theiry geometric average.
	inputs are numpy arrays and are assumed to have the same shape
	'''

	res = np.sqrt(Y**2+X**2)
	res_u = np.sqrt((X*X_u)**2+(Y*Y_u)**2)/np.sqrt((X**2+Y**2))

	return res, res_u

def get_pos_neg_data(a,adwindata_str = '',ro_array = ['positive','negative'],**kw):
	'''
	Input: a : a data object of the class MBIAnalysis
	averages positive and negative data
	returns the sweep points, measured contrast and uncertainty 
	'''

	### get SSRO
	ssro_calib_folder=kw.pop('ssro_calib_folder',None)
	if ssro_calib_folder == None:
	    ssro_calib_timestamp = kw.pop('ssro_calib_timestamp',None)

	    if ssro_calib_timestamp == None: 
	        ssro_calib_folder = toolbox.latest_data('SSRO',**kw)#, older_than = older_than)
	    else:
	        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
	        ssro_calib_folder = toolbox.data_from_time(ssro_calib_timestamp)

	if adwindata_str == '':
		return

	##acquire pos_neg data
	for i,ro in enumerate(ro_array):
		a.get_sweep_pts()
		a.get_readout_results(name = adwindata_str+ro,CR_after_check=CR_after_check)
		a.get_electron_ROC(ssro_calib_folder)


		x_labels = a.sweep_pts.reshape(-1)
		if i == 0:
			res = ((a.p0.reshape(-1))-0.5)*2
			res_u = 2*a.u_p0.reshape(-1)

		else:
			y = ((a.p0.reshape(-1))-0.5)*2 # Contrast
			y_u = 2*a.u_p0.reshape(-1) # contrast
			res = [y0/2-y[ii]/2 for ii,y0 in enumerate(res)]
			res_u = [np.sqrt(y0**2+y_u[ii]**2)/2 for ii,y0 in enumerate(res_u)]

	return np.array(x_labels),np.array(res),np.array(res_u)

def average_repump_time(contains = '',do_fit = False,**kw):
	'''
	gets data from a folder whose name contains the contains variable.
	Does or does not fit the data with a gaussian function
	'''

	### folder choice
	if contains == '':
		contains = 'Sweep_Repump_time'
	elif len(contains) == 2:
		contains = 'Sweep_Repump_time'+contains
	elif len(contains) == 1:
		contains = 'Sweep_Repump_time_'+contains


	# older_than = kw.get('older_than',None) automatically handled by kws
	### acquire data
	f = toolbox.latest_data(contains,**kw)
	a = mbi.MBIAnalysis(f)
	
	if '_Z' in f:
		x,y,y_u = get_pos_neg_data(a,adwindata_str = 'Z_',**kw)
		ylabel = 'Z'
	else:
		x,y1,y1_u  = get_pos_neg_data(a,adwindata_str = 'X_',**kw)
		x2,y2,y2_u = get_pos_neg_data(a,adwindata_str = 'Y_',**kw)
		y,y_u = quadratic_addition(y1,y2,y1_u,y2_u)
		# y=y1
		# y_u = y1_u
		ylabel = 'Bloch vector length'


	### create a plot
	xlabel = a.g.attrs['sweep_name']
	x = a.g.attrs['sweep_pts'] # could potentially be commented out?
	create_plot(f,xlabel = xlabel,ylabel =ylabel,title = 'avg repump time')

	### fitting if you feel like it / still needs implementation
	if do_fit:
		pass

	## plot data
	plot_data(x,y,y_u=y_u)

	## save and close plot. We are done.
	save_and_close_plot(f)


def number_of_repetitions(contains = '', do_fit = False, **kw):
	'''
	gets data from a folder whose name contains the contains variable.
	Does or does not fit the data with a gaussian function
	'''

	### folder choice
	if contains == '':
		contains = '_sweep_number_of_reps'
	elif len(contains) == 2:
		contains = '_sweep_number_of_reps'+contains
	elif len(contains) == 1:
		contains = '_sweep_number_of_reps_'+contains


	# older_than = kw.get('older_than',None) automatically handled by kws
	### acquire data
	f = toolbox.latest_data(contains,**kw)
	a = mbi.MBIAnalysis(f)
	
	if '_Z' in f:
		x,y,y_u = get_pos_neg_data(a,adwindata_str = 'Z_',**kw)
		ylabel = 'Z'
	else:
		x,y1,y1_u  = get_pos_neg_data(a,adwindata_str = 'X_',**kw)
		x2,y2,y2_u = get_pos_neg_data(a,adwindata_str = 'Y_',**kw)
		y,y_u = quadratic_addition(y1,y2,y1_u,y2_u)
		# y=y1
		# y_u = y1_u
		ylabel = 'Bloch vector length'


	### create a plot
	xlabel = a.g.attrs['sweep_name']
	x = a.g.attrs['sweep_pts'] # could potentially be commented out?
	create_plot(f,xlabel = xlabel,ylabel =ylabel,title = 'avg repump time')

	### fitting if you feel like it / still needs implementation
	if do_fit:
		pass

	## plot data
	plot_data(x,y,y_u=y_u)

	## save and close plot. We are done.
	save_and_close_plot(f)


def el_to_c_swap(contains = '',input_el=['Z'], do_fit = False, **kw):
	'''
	gets data from a folder whose name contains the contains variable.
	Does or does not fit the data with a gaussian function
	'''

	### folder choice
	if contains == '':
		contains = 'Swap_el_to_C'

	# older_than = kw.get('older_than',None) automatically handled by kws
	### acquire data
	f = toolbox.latest_data(contains,**kw)
	a = mbi.MBIAnalysis(f)
	
	for el in input_el:
		ro_str = 'el_state_'+el+'_'

		ro_array = ['positive','negative']
		x,y,y_u = get_pos_neg_data(a,adwindata_str = ro_str,ro_array=ro_array,**kw)
		print get_tstamp_from_folder(f),el,x,y,y_u
	# ### create a plot XXX TODO
	# xlabel = a.g.attrs['sweep_name']
	# x = a.g.attrs['sweep_pts'] # could potentially be commented out?
	# create_plot(f,xlabel = xlabel,ylabel =ylabel,title = 'avg repump time')

	# ### fitting if you feel like it / still needs implementation
	# if do_fit:
	# 	pass

	# ## plot data
	# plot_data(x,y,y_u=y_u)

	# ## save and close plot. We are done.
	# save_and_close_plot(f)
