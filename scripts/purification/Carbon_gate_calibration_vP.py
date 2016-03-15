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


'''
### Instructions ~ SK 2016
### Reads out data from a gate N and tau sweep, where a number of central taus may have been specified
### 1. This script works in two ways. You can either seperately return the data 
### 	and then use the methods defined to plot it. Or use the function directly to plot all the stuff, just set kw's to true
### 2. plot_gate_time and line_plot_fidelity can both take lists to plot multiple datasets
### bar plot can not, since you will never want to overlap bar plots anyway (will be a mess)
'''

def get_bloch_length(x,y,x_u,y_u):

	### init
	b_arr,b_u_arr = np.array([]),np.array([])

	##calc
	b_arr = np.sqrt(x**2+y**2)
	b_u_arr =  np.sqrt((x*x_u)**2+(y*y_u)**2)/np.sqrt((x**2+y**2))
	return b_arr,b_u_arr


def get_raw_data_all_parts(carbon,**kw):
	"""
	extracts the data for one carbon from a positive and negative file.
	returns arrays for the contrast along x, y and the respective uncertainties.
	"""
	older_than = kw.pop('older_than',None)
	newer_than = kw.pop('newer_than',None)
	ssro_tstamp = kw.pop('ssro_tstamp',None)
	tau_nr		= kw.pop('tau_nr',None)

	if ssro_tstamp == None:
		ssro_calib_folder = toolbox.latest_data(contains = 'SSRO')
		print ssro_calib_folder
	else:
		ssro_calib_folder = toolbox.latest_data(contains = ssro_tstamp + '_AdwinSSRO')
		print ssro_calib_folder

	if tau_nr == None:
		search_string_pos = 'Sweep_carbon_Gate_positive_C'+str(carbon)
		search_string_neg = 'Sweep_carbon_Gate_negative_C'+str(carbon)
	else: 
		search_string_pos = 'Sweep_carbon_Gate__C'+str(carbon)+ '_positive_tau' + str(tau_nr) + '_'
		search_string_neg = 'Sweep_carbon_Gate__C'+str(carbon)+ '_negative_tau' + str(tau_nr) + '_'


	#list of all folders meeting the requirements
	#this was needed because I couldn't get the toolbox to work for fragmented folder names build of multiple parts SK
	folder_pos_list_reversed = toolbox.latest_data(contains = search_string_pos, older_than=older_than, newer_than=newer_than, return_all = True)
	folder_neg_list_reversed = toolbox.latest_data(contains = search_string_neg, older_than=older_than, newer_than=newer_than, return_all = True)
	print 'Number of parts: ' + str(len(folder_pos_list_reversed))

	#flip the folder list or we could loop from the back
	folder_pos_list = list(reversed(folder_pos_list_reversed))
	folder_neg_list = list(reversed(folder_neg_list_reversed))
	# print folder_pos_list
	# print folder_neg_list

	# Initialization of arrays
	x_arr,x_u_arr = np.array([]),np.array([])
	y_arr,y_u_arr = np.array([]),np.array([])
	gates = np.array([]) ### used parameters
	gate_values = [] ### used parameters but not in string form

	for f_index in range(len(folder_pos_list)):
		# print folder_pos_list[f_index]
		# print folder_neg_list[f_index]
		folder_pos = toolbox.latest_data(contains = folder_pos_list[f_index][-40:], older_than=older_than, newer_than=newer_than)
		folder_neg = toolbox.latest_data(contains = folder_neg_list[f_index][-40:], older_than=older_than, newer_than=newer_than)

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
	
		### sort into X and y lists.
		for (pt,val,val_u) in zip(a.sweep_pts.reshape(-1),a.p0.reshape(-1),a.u_p0.reshape(-1)):
			if 'X' in pt:
				x_arr = np.append(x_arr,val)
				x_u_arr = np.append(x_u_arr,val_u)
				#could be made more efficient
				gate_values.append([pt[2:4] ,pt[7:]])
				gates = np.append(gates,'N = '+str(pt[2:4])+',\ntau = '+str(pt[7:]))
			elif 'Y' in pt:
				y_arr = np.append(y_arr,val)
				y_u_arr = np.append(y_u_arr,val_u)

	return gates,x_arr,y_arr,x_u_arr,y_u_arr,folder_pos, gate_values

# def get_raw_data_original(carbon,**kw):
# 	"""
# 	extracts the data for one carbon from a positive and negative file.
# 	returns arrays for the contrast along x, y and the respective uncertainties.
# 	"""
# 	older_than = kw.pop('older_than',None)
# 	newer_than = kw.pop('newer_than',None)
# 	ssro_tstamp = kw.pop('ssro_tstamp',None)

# 	if ssro_tstamp == None:
# 		ssro_calib_folder = toolbox.latest_data(contains = 'SSRO')
# 		print ssro_calib_folder

# 	search_string_pos = 'Sweep_carbon_Gate__C'+str(carbon)+ '_positive'
# 	search_string_neg = 'Sweep_carbon_Gate__C'+str(carbon)+ '_negative'

# 	folder_pos = toolbox.latest_data(contains = search_string_pos, older_than=older_than)
# 	folder_neg = toolbox.latest_data(contains = search_string_neg, older_than=older_than)

# 	a = mbi.MBIAnalysis(folder_pos)
# 	a.get_sweep_pts()
# 	a.get_readout_results(name='adwindata')
# 	a.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)

# 	b = mbi.MBIAnalysis(folder_neg)
# 	b.get_sweep_pts()
# 	b.get_readout_results(name='adwindata')
# 	b.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)

# 	a.p0 = 2*a.p0-1; a.u_p0 = 2*a.u_p0
# 	b.p0 = 2*b.p0-1; b.u_p0 = 2*b.u_p0

# 	###combine positive & negative:

# 	a.p0 = (a.p0 - b.p0)/2
# 	a.u_p0 = ((a.u_p0**2 + b.u_p0**2)**0.5)/2


# 	x_arr,x_u_arr = np.array([]),np.array([])
# 	y_arr,y_u_arr = np.array([]),np.array([])
# 	gates = np.array([]) ### used parameters

# 	### sort into X and y lists.
# 	for (pt,val,val_u) in zip(a.sweep_pts.reshape(-1),a.p0.reshape(-1),a.u_p0.reshape(-1)):
# 		if 'X' in pt:
# 			x_arr = np.append(x_arr,val)
# 			x_u_arr = np.append(x_u_arr,val_u)
# 			gates = np.append(gates,'N = '+str(pt[2:4])+',\ntau = '+str(pt[5:]))
# 		elif 'Y' in pt:
# 			y_arr = np.append(y_arr,val)
# 			y_u_arr = np.append(y_u_arr,val_u)


# 	return gates,x_arr,y_arr,x_u_arr,y_u_arr,folder_pos

def get_gate_fidelity(carbon, **kw):
	"""
	gets data, plots it and prints the gate parameters for maximum bloch vector length.
	"""

	older_than = kw.pop('older_than',None)
	newer_than = kw.pop('newer_than',None)
	ssro_tstamp = kw.pop('ssro_tstamp',None)
	tau_nrs 	= kw.pop('tau_nrs', None)
	return_data = kw.pop('return_data',False)
	plot_fidelity = kw.pop('plot_fidelity',False)
	gate_time = kw.pop('gate_time',False)
	line_fidelity = kw.pop('line_fidelity',False)
	# loop over each tau

	if tau_nrs != None:
		for t in tau_nrs:
			
			gates,x,y,x_u,y_u,folder_pos, gate_values = get_raw_data_all_parts(carbon,older_than = older_than,newer_than = newer_than,
				ssro_tstamp = ssro_tstamp, tau_nr = t)
			
			b,b_u = get_bloch_length(x,y,x_u,y_u)

			best_b = np.amax(b)
			best_b_ind = np.argmax(b)

			print 'best gate configuration at: ', gates[best_b_ind]
			print 'bloch vector length: ', best_b

			if plot_fidelity:
				bar_plot_fidelity(gates,gate_values,b,b_u)

			if gate_time:
				plot_gate_time([gates],[gate_values],[b],[b_u])

			if line_fidelity:
				line_plot_fidelity([gates],[gate_values],[b],[b_u])

			if return_data:
				return gates,gate_values,b,b_u
				print 'data returned'
	else:
		gates,x,y,x_u,y_u,folder_pos, gate_values = get_raw_data_all_parts(carbon,older_than = older_than,newer_than = newer_than,
				ssro_tstamp = ssro_tstamp)
			
		b,b_u = get_bloch_length(x,y,x_u,y_u)

		best_b = np.amax(b)
		best_b_ind = np.argmax(b)

		print 'best gate configuration at: ', gates[best_b_ind]
		print 'bloch vector length: ', best_b
		
		if plot_fidelity:
			bar_plot_fidelity(gates,gate_values,b,b_u)

		if gate_time:
			plot_gate_time([gates],[gate_values],[b],[b_u])

		if line_fidelity:
			line_plot_fidelity([gates],[gate_values],[b],[b_u])

		if return_data:
			return gates,gate_values,b,b_u
			print 'data returned'
	


	print 'mission complete'


def bar_plot_fidelity(gates,gate_values,b,b_u):
	'''Bar plot fidelity'''
	# gates,gate_values,b,b_u = self.get_gate_fidelity(carbon=carbon)
	fig, ax = plt.subplots(figsize=(len(gates)/2.,5))
	# ax = plt.subplot()
	print len(b)
	print len(gates)
	rects= ax.bar(np.arange(len(gates)),b,yerr=b_u, width = 0.8)#,align='center')
	ax.set_xticks(np.arange(len(gates)))
	ax.set_xticks(np.array(range(len(gates))))
	ax.set_xticklabels(gates, rotation=90)
	plt.xlabel('Gate configuration')
	plt.ylabel('Bloch vector length')

	for ii,rect in enumerate(rects):
		plt.text(rect.get_x()+rect.get_width()/2., 0.4*rect.get_height(), 'F='+str(round(b[ii],3))+' $\pm$ '+str(round(b_u[ii],3)),
			ha='center', va='bottom',rotation='vertical',color='white')

def line_plot_fidelity(gates,gate_values,b,b_u):
	plt.figure()
	for i in range(len(gates)):
		plt.errorbar(range(len(gates[i])),b[i],b_u[i],label = 'gates' + str(i+1))
		# plt.savefig(os.path.join(folder_pos, 'Sweep_gates.png'), format='png')
	plt.xlabel('Gates, same order as barplot')
	plt.ylabel('Bloch vector length')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.show()

def plot_gate_time(gates,gate_values,b,b_u):
		##### plot gate time vs fidelity
	plt.figure()
	for i in range(len(gate_values)):
		gate_time = []
		for g in gate_values[i]:
			# print 'statements'
			# print int(g[0])
			# print float(g[1])
			# print int(g[0])*float(g[1])
			gate_time.append(int(g[0])*float(g[1]))
		# print gate_time
		plt.plot(gate_time, b[i], '.',label = 'gates' + str(i+1))
	plt.xlabel('Gate time (N*tau)')
	plt.ylabel('Bloch vector length')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


# # Should actually generalize this into all
# def compare_msmts(carbons = [5,5],
# 					newer_than = [None,None],
# 					older_than = [None,None],
# 					tau_nrs = [0,0]):
	
# 	plot_fidelity = kw.pop('plot_fidelity',False)
# 	plot_gate_time = kw.pop('plot_gate_time',False)
# 	plot_graph = kw.pop('plot_graph',False)

# 	for i in len(carbons):
# 		self.get_gate_fidelity(carbon = carbons[i], newer_than=newer_than[i], older_than= None, tau_nrs = [0],
#                      plot_fidelity = plot_fidelity, plot_gate_time = plot_gate_time, plot_graph=plot_graph)


print 'Mission complete'



