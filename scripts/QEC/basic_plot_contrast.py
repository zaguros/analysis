import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
from analysis.lib.fitting import fit, common

def simple_plot(timestamp = None, measurement_name = ['adwindata'],folder_name ='CarbonPiCal',
		ssro_calib_timestamp =None, save = True,
		plot_fit = True) :
	'''
	Function that makes a bar plot with errorbars of MBI type data 
	'''
	if timestamp == None:
	    timestamp, folder   = toolbox.latest_data(folder_name,return_timestamp =True)
	else: 
	    folder = toolbox.data_from_time(timestamp) 

	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSRO')
	else:
		ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
		ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil8'
		print ssro_calib_folder


	a = mbi.MBIAnalysis(folder)
	a.get_sweep_pts()
	a.get_readout_results(name='adwindata')
	a.get_electron_ROC(ssro_calib_folder)

	x_labels = a.sweep_pts.reshape(-1)
	y= ((a.p0.reshape(-1))-0.5)*2
	x = range(len(y)) 
	y_err = 2*a.u_p0.reshape(-1)

	fig,ax = plt.subplots() 
	rects = ax.errorbar(x,y,yerr=y_err)
	ax.set_xticks(x)
	ax.set_xticklabels(x_labels.tolist(),rotation=90)
	ax.set_ylim(-1.1,1.1)
	ax.set_title(str(folder)+'/'+str(timestamp))
	ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

	if save and ax != None:
		try:
		    fig.savefig(
		        os.path.join(folder,'simple_plot.png'))
		except:
		    print 'Figure has not been saved.'

def simple_plot_contrast(timestamps = [None,None], tag = '', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp =None, save = True,
		do_plot = True, return_data = False) :
	'''
	Function that makes a bar plot with errorbars of MBI type data that has been measured with a positive
	and negative RO.
	'''

	### SSRO calibration
	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSRO')
	else:
	    ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
	    ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil8'

	if timestamps[0] == None: 
		folder_a = toolbox.latest_data(contains='positive' + tag)
		folder_b = toolbox.latest_data(contains='negative' + tag)
	elif len(timestamps)==1:		
		folder_b = toolbox.data_from_time(timestamps[0])      
		print folder_b
		folder_a = toolbox.latest_data(contains = 'pos', older_than = timestamps[0])   
		print folder_a
	else:
		folder_a = toolbox.data_from_time(timestamps[0])      
		folder_b = toolbox.data_from_time(timestamps[1])     	   
 	
	a = mbi.MBIAnalysis(folder_a)
	a.get_sweep_pts()
	a.get_readout_results(name='adwindata')
	a.get_electron_ROC(ssro_calib_folder)
	y_a= ((a.p0.reshape(-1)[:])-0.5)*2
	y_err_a = 2*a.u_p0.reshape(-1)[:] 


	b = mbi.MBIAnalysis(folder_b)
	b.get_sweep_pts()
	b.get_readout_results(name='adwindata')
	b.get_electron_ROC(ssro_calib_folder)
	y_b= ((b.p0.reshape(-1)[:])-0.5)*2
	y_err_b = 2*b.u_p0.reshape(-1)[:] 

	x = a.sweep_pts.reshape(-1)[:]
	# x = range(len(y_a)) 


	
	### Combine data
	y = (y_a - y_b)/2.
	y_err =  1./2*(y_err_a**2 + y_err_b**2)**0.5 
	if do_plot == True:
		fig,ax = plt.subplots() 
		rects = ax.errorbar(x,y,yerr=y_err)
		ax.set_xticks(x)
		# ax.set_xticklabels(x_labels.tolist())
		ax.set_ylim(-1.1,1.1)
		print folder_a
		ax.set_title(str(folder_a)+'/'+str(timestamps[0]))
		ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')



	if save and ax != None:
		try:
		    fig.savefig(
		        os.path.join(folder_a,'simple_plot_contrast.png'))
		except:
		    print 'Figure has not been saved.'

	if return_data == True:
		return  x, y, y_err

def simple_plot_contrast_concatenate(older_than = None,tag_list = [], do_fit = False,measurement_name = ['adwindata'],
		ssro_calib_timestamp =None, save = True,
		do_plot = True,absolute_y = False,
		xlabel = 'time', ylabel = 'Contrast',
		**kw):
		### SSRO calibration

	for k in range(len(tag_list)):

		timestamp_pos, folder_a = toolbox.latest_data(contains = 'positive'+ tag_list[k], older_than = older_than,return_timestamp = True)
		timestamp_neg, folder_b = toolbox.latest_data(contains = 'negative'+ tag_list[k], older_than = older_than,return_timestamp = True)
	
		print folder_a
		print folder_b

		x_t, y_t, y_err_t  = simple_plot_contrast(timestamps = [timestamp_pos,timestamp_neg], measurement_name = ['adwindata'],
								ssro_calib_timestamp =None, save = False,
								do_plot = False, return_data = True)
		
		if k == 0:
			# print type(x_labels_t)
			# print list(x_labels_t)
			# x_labels = list(x_labels_t)
			x = list(x_t)
			y = list(y_t)
			y_err = list(y_err_t)
		else:
			# x_labels.extend(list(x_labels_t))
			x.extend(list(x_t))
			y.extend(list(y_t))
			y_err.extend(list(y_err_t))
	# x = range(len(y)) 
		# print x_labels.tolist
		# print y.tolist
	if absolute_y == True:
		y = np.absolute(y)

	if do_plot ==True: 
		fig,ax = plt.subplots(figsize=(4,4)) 
		if do_fit == False:
			rects = ax.errorbar(x,y,yerr=y_err )
		else: 
			ax.errorbar(x,y,yerr = y_err, marker = 'o',ls = '')
		ax.set_xticks(x)
		# ax.set_xticklabels(x_labels)
		ax.set_ylim(-1.1,1.1)
		if absolute_y == True:
			ax.set_ylim(-0.1,1.1)
		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)
		ax.set_title(str(folder_a)[18:])
		ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

		if do_fit == True:
		    o = fit.Parameter(kw.pop('guess_o',0), 'o')
		    f = fit.Parameter(kw.pop('guess_freq',1/10.), 'f')
		    A = fit.Parameter(kw.pop('guess_amp',1), 'A')
		    phi = fit.Parameter(kw.pop('guess_phi',0), 'phi')
		    k = fit.Parameter(kw.pop('guess_k',10.), 'k')
		    p0 = [A,o,f,phi,k]
		    fitfunc_str = ''
		    

		    def fitfunc(x):
		        return (o()-A()) + A() * np.exp(-(k()*x)**2) * np.cos(2*np.pi*(f()*x - phi()))

		    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=[1,3],
		            do_print=True, ret=True)

	if save and ax != None:
		try:
		    fig.savefig(
		        os.path.join(folder_a,'fullplot.png'))
		except:
		    print 'Figure has not been saved.'
