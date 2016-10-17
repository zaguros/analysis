import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi; reload(mbi)
from matplotlib import pyplot as plt

import h5py

def BarPlotTomo(timestamp = None, measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp =None, save = True,tag = '',
		plot_fit = True) :
	'''
	Function that makes a bar plot with errorbars of MBI type data 
	'''
	if timestamp == None:
		if tag != '':
			timestamp,folder = toolbox.latest_data(contains=tag,return_timestamp=True)
		else:
			timestamp, folder   = toolbox.latest_data(folder_name,return_timestamp =True)
	else: 
	    folder = toolbox.data_from_time(timestamp) 

	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSROCalibration')
	else:
		ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
		ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Pippin_SIL2'
		print ssro_calib_folder

	a = mbi.MBIAnalysis(folder)

	a.get_sweep_pts()
	a.get_readout_results(name=measurement_name[0])
	a.get_electron_ROC(ssro_calib_folder)
	# a.p0 = a.normalized_ssro
	# a.u_p0 = a.u_normalized_ssro
	x_labels = a.sweep_pts.reshape(-1)
	y= ((a.p0.reshape(-1))-0.5)*2
	x = range(len(y)) 
	y_err = 2*a.u_p0.reshape(-1)
	# print 'y', y
	# print 'err', y_err
	if plot_fit ==True: 
		fig,ax = plt.subplots() 
		rects = ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
		ax.set_xticks(x)
		# ax.title = timestamp
		print x_labels
		ax.set_xticklabels(x_labels.tolist())
		ax.set_ylim(-1.1,1.1)
		ax.set_title(str(folder)+'/'+str(timestamp))
		# ax.grid()
		ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')


			# print values on bar plot
	def autolabel(rects):
	    for ii,rect in enumerate(rects):
	        height = rect.get_height()
	        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2)) +'('+ str(int(round(y_err[ii]*100))) +')',
	            ha='center', va='bottom')
	autolabel(rects)

	if save and ax != None:
		try:
		    fig.savefig(
		        os.path.join(folder,'tomo.png'))
		    # fig.savefig(os.path.join(folder,'tomo.eps'))
		except:
		    print 'Figure has not been saved.'

def BarPlotTomoContrast(timestamps = [None,None], tag = '', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp =None, save = True,
		plot_fit = True, return_data = False,older_than = None) :
	'''
	Function that makes a bar plot with errorbars of MBI type data that has been measured with a positive
	and negative RO.
	'''

	### SSRO calibration
	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSROCalibration')
	else:
	    ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
	    ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'

	### Obtain and analyze data
		### postive RO data

	# if timestamps[0][0] == 'C': 
	# 	folder_a = toolbox.latest_data(contains=timestamps[0]+'_positive')
	# 	folder_b = toolbox.latest_data(contains=timestamps[0]+'_negative')	
	print 'older_than', older_than
	if timestamps[0] == None: 
		folder_a = toolbox.latest_data(contains='positive' + tag,older_than = older_than)
		folder_b = toolbox.latest_data(contains='negative' + tag,older_than = older_than)
		print folder_a
		print folder_b

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

	x_labels = a.sweep_pts.reshape(-1)[:]
	x = range(len(y_a)) 

	
	### Combine data
	y = (y_a - y_b)/2.
	y_err =  1./2*(y_err_a**2 + y_err_b**2)**0.5 
	
	# print folder_a
	# print folder_b

	# print y_a
	# print y_b
	# print y


	### Fidelities
	# F_ZZ 	= (1 + y[2] + y[5] + y[14])/4
	# F_ZmZ 	= (1 + y[2] - y[5] - y[14])/4
	# F_ent 	= (1 + y[0] -y[4] -y[8])/4
	# F_ent 	= (1 + y[0] +y[1] +y[2])/4
	# print 'Fidelity with ZZ  = ' + str(F_ZZ)
	# print 'Fidelity with ZmZ  = ' + str(F_ZmZ)
	# print 'Fidelity with ent = ' + str(F_ent)

	# print 'XY = ' +str( (y[0]**2 + y[1]**2)**0.5)

	if plot_fit ==True: 
		fig,ax = plt.subplots() 
		rects = ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
		ax.set_xticks(x)
		ax.set_xticklabels(x_labels.tolist()) ## if you want rotated labels: ,rotation=90
		ax.set_ylim(-1.1,1.1)
		print 'test'
		print folder_a
		ax.set_title(str(folder_a)+'/'+str(timestamps[0]))
		ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

			# print values on bar plot
		def autolabel(rects):
		    for ii,rect in enumerate(rects):
		        height = rect.get_height()
		        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2)) +'('+ str(int(round(y_err[ii]*100))) +')',
		            ha='center', va='bottom')
		autolabel(rects)

	if save and ax != None:
		try:
		    fig.savefig(
		        os.path.join(folder_a,'tomo.png'))
		except:
		    print 'Figure has not been saved.'

	if return_data == True:
		return x_labels, x, y, y_err

def CompleteTomo(timestamp = None, 
		state = 'Z', 
		measurement_name = ['adwindata'],
		folder_tag ='Tomo',
		ssro_calib_timestamp =None, 
		save = True,
		older_than = None,
		newer_than = None,
		return_num = "All",
		return_newest_or_oldest = 'Newest',
		plot_fit = True):

	### Initialize with empty lists
	x_labels =[]
	x = []
	y = []
	y_err = []

	# find the folders you'd like to use, contains both negative and postive
	folder_list = toolbox.latest_data(contains = folder_tag, 
			older_than = older_than, 
			newer_than = newer_than, 
			return_all = True)[::-1]
	print 'number of folders: '  + str(len(folder_list))

	folder_pos_list = [s for s in folder_list if "positive" in s]
	folder_neg_list = [s for s in folder_list if "negative" in s]

	if return_num != "All": # Return only the latest files.
		if return_newest_or_oldest == 'Newest':
			folder_pos_list = folder_pos_list[-return_num:]
			folder_neg_list = folder_neg_list[-return_num:]
		elif return_newest_or_oldest == 'Oldest':
			folder_pos_list = folder_pos_list[:return_num]
			folder_neg_list = folder_neg_list[:return_num]
		else:
			print "Unknown option for return_newest_or_oldest"


	# TODO, improve by doing correct axis etc. and maybe max only
	if len(folder_neg_list) == 0:
		pass
		# write so it is compatible with positive data only
		# BarPlotTomo
	else:
		for p,n in zip(folder_pos_list,folder_neg_list):
			
			# gets time stamps so it doesn't conflict with older tomo version
			tp = p.split('\\')[3] + '_' + p.split('\\')[4][:6]
			tn = n.split('\\')[3] + '_' + n.split('\\')[4][:6]

			# try:
			# 	a = ('mX'|'mY'|'mZ') in p.split('\\')[4]
			# 	print a
			# except:
			# 	pass
			# get tomo data from timestamp
			x_labels_t, x_t, y_t, y_err_t  = BarPlotTomoContrast(
				timestamps = [tp,tn],
				ssro_calib_timestamp =None, 
				save = False,
				plot_fit = False, 
				return_data = True)

			x_labels.extend(list(x_labels_t))
			# x.extend(list(x_t))
			print 10*'*'
			print x_labels_t
			print y_t
			print 10*'*'

			y.extend(list(y_t))
			y_err.extend(list(y_err_t))

	x = range(len(y))

	if plot_fit ==True: 
		fig,ax = plt.subplots(figsize=(16,4)) 
		rects = ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
		ax.set_xticks(x)
		ax.set_xticklabels(x_labels)
		ax.set_ylim(-1.1,1.1)
		ax.set_title(str(folder_pos_list[0])+'...\n'+str(folder_neg_list[-1]))
		ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

			# print values on bar plot
		def autolabel(rects):
		    for ii,rect in enumerate(rects):
		        height = rect.get_height()
		        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2)) +'('+ str(int(round(y_err[ii]*100))) +')',
		            ha='center', va='bottom')
		autolabel(rects)

	if save:# and ax != None:
		try:
		    fig.savefig(
		        os.path.join(folder_pos_list[0],'tomo.png'))
		except:
		    print 'Figure has not been saved.'

def BarPlotTomoContrastFull(timestamp = None, state = 'Z', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp =None, save = True,
		plot_fit = True):
		### SSRO calibration

	for k in range(21):
		
		print k 
		print 'positive'+ str(k)
		if k < 9:
			print 'no'

			timestamp_9, folder_9 = toolbox.latest_data(contains = 'state_' + state +'_positive_'+ str(9), older_than = timestamp,return_timestamp = True)
			timestamp_pos, folder_a = toolbox.latest_data(contains = 'state_' + state +'_positive_'+ str(k), older_than = timestamp_9,return_timestamp = True)
			timestamp_neg, folder_b = toolbox.latest_data(contains = 'state_' + state +'_negative_'+ str(k), older_than = timestamp_9,return_timestamp = True)
		else:
			print 'yes'
			timestamp_pos, folder_a = toolbox.latest_data(contains = 'state_' + state +'_positive_'+ str(k), older_than = timestamp,return_timestamp = True)
			timestamp_neg, folder_b = toolbox.latest_data(contains = 'state_' + state +'_negative_'+ str(k), older_than = timestamp,return_timestamp = True)
		
		print folder_a
		print folder_b

		x_labels_t, x_t, y_t, y_err_t  = BarPlotTomoContrast(timestamps = [timestamp_pos,timestamp_neg], measurement_name = ['adwindata'],folder_name ='Tomo',
								ssro_calib_timestamp =None, save = False,
								plot_fit = False, return_data = True)
		if k == 0:
			print type(x_labels_t)
			print list(x_labels_t)
			x_labels = list(x_labels_t)
			x = list(x_t)
			y = list(y_t)
			y_err = list(y_err_t)
		else:
			x_labels.extend(list(x_labels_t))
			x.extend(list(x_t))
			y.extend(list(y_t))
			y_err.extend(list(y_err_t))
	x = range(len(y)) 
		# print x_labels.tolist
		# print y.tolist

	if plot_fit ==True: 
		fig,ax = plt.subplots(figsize=(35,5)) 
		rects = ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
		ax.set_xticks(x)
		ax.set_xticklabels(x_labels)
		ax.set_ylim(-1.1,1.1)
		ax.set_title(str(folder_a)+'/'+str(timestamp_neg))
		ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

			# print values on bar plot
		def autolabel(rects):
		    for ii,rect in enumerate(rects):
		        height = rect.get_height()
		        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2)) +'('+ str(int(round(y_err[ii]*100))) +')',
		            ha='center', va='bottom')
		autolabel(rects)

	if save and ax != None:
		try:
		    fig.savefig(
		        os.path.join(folder_a,'tomo.png'))
		except:
		    print 'Figure has not been saved.'


def BarPlotTomoContrastFull_mult_msmts(timestamp = None, measurement_name = ['adwindata'],folder_name ='Tomo', older_than = None, newer_than = None, ssro_calib_timestamp =None, save = True, plot_fit = True):
	### SSRO calibration

	Tomo_str_list = []
	Tomo_basis_list = ([
	            ['X','I'],['Y','I'],['Z','I'],
	            ['I','X'],['I','Y'],['I','Z'],
	            ['X','X'],['X','Y'],['X','Z'],
	            ['Y','X'],['Y','Y'],['Y','Z'],
	            ['Z','X'],['Z','Y'],['Z','Z']])

	for BP in Tomo_basis_list:
		Tomo_str_list.append(''.join([str(s) for s in BP]))


	for k, BP in enumerate(Tomo_str_list):
		print 'Tomo' + BP +'_ROpositive_'
		timestamp_pos, folder_a = toolbox.latest_data(contains = 'sweep_phase_FET0.305s_auto_C1&2_Tomo' + BP +'_ROpositive_', older_than = older_than, newer_than = newer_than, return_timestamp = True)
		timestamp_neg, folder_b = toolbox.latest_data(contains = 'sweep_phase_FET0.305s_auto_C1&2_Tomo' + BP +'_ROnegative_', older_than = older_than, newer_than = newer_than, return_timestamp = True)
		print folder_a
		print folder_b
		x_labels_t, x_t, y_t, y_err_t  = BarPlotTomoContrast(timestamps = [timestamp_pos,timestamp_neg], measurement_name = ['adwindata'],
								ssro_calib_timestamp =None, save = False,
								plot_fit = False, return_data = True)
		if k == 0:
			print type(x_labels_t)
			print list(x_labels_t)
			x_labels = list(x_labels_t)
			x = list(x_t)
			y = list(y_t)
			y_err = list(y_err_t)
		else:
			x_labels.extend(list(x_labels_t))
			x.extend(list(x_t))
			y.extend(list(y_t))
			y_err.extend(list(y_err_t))
	
	x = range(len(y)) 
	# print x_labels.tolist
	# print y.tolist

	squared = 0.
	fidel = 1.
	for ii in x:
		if y[ii]>0.2:
			print x_labels[ii]
			print y[ii]
			fidel += y[ii]
			squared += y_err[ii]**2
			print y_err[ii]

	print 'rooterror', squared**0.5
	print 'fidelity' , fidel / 4.
	if plot_fit ==True: 
		fig,ax = plt.subplots(figsize=(7,5)) 
		rects = ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
		print y
		y2 = [0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,1.,0.]

		rects2 = ax.bar(x,y2,align ='center', alpha=0.2 )
		ax.set_xticks(x)
		ax.set_xticklabels(x_labels,fontsize=15,rotation='vertical')
		ax.tick_params(axis='y', which='major', labelsize=15)
		ax.set_ylim(-1.1,1.1)
		# ax.set_title(str(folder_a)+'/'+str(timestamp_neg))
		ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted', linewidth=2)
		ax.set_ylabel('Expectation value', fontsize=15)
		# plt.savefig(r'D:\Dropbox\QEC LT\Decoupling memory\00_Thesis_plots\Tomo12tight.pdf',bbox_inches='tight')
			# print values on bar plot
		def autolabel(rects):
		    for ii,rect in enumerate(rects):
		        height = rect.get_height()
		        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2)) +'('+ str(int(round(y_err[ii]*100))) +')',
		            ha='center', va='bottom', fontsize=12)
		autolabel(rects)

	if save and ax != None:
		try:
		    fig.savefig(
		        os.path.join(folder_a,'tomo.png'))
		except:
		    print 'Figure has not been saved.'