import numpy as np
import os,sys
# import qt

sys.path.append("/measuring/")
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import matplotlib.cm as cm
reload(common)
reload(plot)


def get_bloch_length(x,y,x_u,y_u):

	### init
	b_arr,b_u_arr = np.array([]),np.array([])

	##calc
	b_arr = np.sqrt(x**2+y**2)
	b_u_arr =  np.sqrt((x*x_u)**2+(y*y_u)**2)/np.sqrt((x**2+y**2))
	return b_arr,b_u_arr


def get_raw_data(carbon,**kw):
	"""
	extracts the data for one carbon from a positive and negative file.
	returns arrays for the contrast along x, y and the respective uncertainties.
	"""
	older_than = kw.pop('older_than',None)
	ssro_tstamp = kw.pop('ssro_tstamp',None)

	if ssro_tstamp == None:
		ssro_calib_folder = toolbox.latest_data(contains = 'SSRO')
		print ssro_calib_folder

	search_string_pos = 'Sweep_carbon_Gate_positive_C'+str(carbon)
	search_string_neg = 'Sweep_carbon_Gate_negative_C'+str(carbon)

	folder_pos = toolbox.latest_data(contains = search_string_pos, older_than=older_than)
	folder_neg = toolbox.latest_data(contains = search_string_neg, older_than=older_than)

	a = mbi.MBIAnalysis(folder_pos)
	a.get_sweep_pts()
	a.get_readout_results(name='adwindata')
	a.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)

	b = mbi.MBIAnalysis(folder_neg)
	b.get_sweep_pts()
	b.get_readout_results(name='adwindata')
	b.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)

	a.p0 = 2*a.p0-1; a.u_p0 = 2*a.u_p0
	b.p0 = 2*b.p0-1; b.u_p0 = 2*b.u_p0

	###combine positive & negative:

	a.p0 = (a.p0 - b.p0)/2
	a.u_p0 = ((a.u_p0**2 + b.u_p0**2)**0.5)/2


	x_arr,x_u_arr = np.array([]),np.array([])
	y_arr,y_u_arr = np.array([]),np.array([])
	gates = np.array([]) ### used parameters

	### sort into X and y lists.
	for (pt,val,val_u) in zip(a.sweep_pts.reshape(-1),a.p0.reshape(-1),a.u_p0.reshape(-1)):
		if 'X' in pt:
			x_arr = np.append(x_arr,val)
			x_u_arr = np.append(x_u_arr,val_u)
			gates = np.append(gates,'N = '+str(pt[2:4])+',\ntau = '+str(pt[5:]))
		elif 'Y' in pt:
			y_arr = np.append(y_arr,val)
			y_u_arr = np.append(y_u_arr,val_u)


	return gates,x_arr,y_arr,x_u_arr,y_u_arr,folder_pos

def get_gate_fidelity(carbon,**kw):
	"""
	gets data, plots it and prints the gate parameters for maximum bloch vector length.
	"""

	older_than = kw.pop('older_than',None)
	ssro_tstamp = kw.pop('ssro_tstamp',None)

	gates,x,y,x_u,y_u,folder_pos = get_raw_data(carbon,older_than = older_than,ssro_tstamp = ssro_tstamp)
	print folder_pos
	b,b_u = get_bloch_length(x,y,x_u,y_u)

	best_b = np.amax(b)
	best_b_ind = np.argmax(b)

	print 'best gate configuration at: ', gates[best_b_ind]
	print 'bloch vector length: ', best_b

	fig = plt.figure()
	ax = plt.subplot()

	rects= ax.bar(np.arange(len(gates)),b,yerr=b_u,align='center')
	ax.set_xticks(np.arange(len(gates)))
	ax.set_xticklabels(gates, rotation=90)
	plt.xlabel('Gate configuration')
	plt.ylabel('Bloch vector length')

	for ii,rect in enumerate(rects):
		plt.text(rect.get_x()+rect.get_width()/2., 0.4*rect.get_height(), 'F='+str(round(b[ii],3))+' $\pm$ '+str(round(b_u[ii],3)),
			ha='center', va='bottom',rotation='vertical',color='white')

	plt.savefig(os.path.join(folder_pos, 'Sweep_gates.png'), format='png')
	plt.show()
	plt.close('all')




