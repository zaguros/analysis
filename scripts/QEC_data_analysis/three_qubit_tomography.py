import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import matplotlib as mpl
from analysis.scripts.QEC import Two_Qubit_Tomography_Postselection as tomo_ps
reload(tomo_ps)
from analysis.lib.QEC import ConditionalParity as CP
import matplotlib as mpl; reload(mpl)


folder = r'D:\measuring\data\Analyzed figures\Three qubit states'
folder = r'D:\measuring\data\QEC_data\figs\final figures'

def BarPlotTomo(timestamp = None, measurement_name = ['adwindata'],folder_name ='Tomo',
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
	    ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'

	a = mbi.MBIAnalysis(folder)
	a.get_sweep_pts()
	a.get_readout_results(name='adwindata')
	a.get_electron_ROC(ssro_calib_folder)

	x_labels = a.sweep_pts.reshape(-1)
	y= ((a.p0.reshape(-1))-0.5)*2
	x = range(len(y)) 
	y_err = 2*a.u_p0.reshape(-1)

	if plot_fit ==True: 
		fig,ax = plt.subplots() 
		rects = ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
		ax.set_xticks(x)
		# ax.title = timestamp
		# print x_labels
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
		except:
		    print 'Figure has not been saved.'

def BarPlotTomoContrast(timestamps = [None,None], tag = '', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp =None, save = True,
		plot_fit = True, return_data = False) :
	'''
	Function that makes a bar plot with errorbars of MBI type data that has been measured with a positive
	and negative RO.
	'''

	### SSRO calibration
	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSRO')
	else:
	    ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
	    ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
	    # print ssro_calib_folder
	### Obtain and analyze data
		### postive RO data

	# if timestamps[0][0] == 'C': 
	# 	folder_a = toolbox.latest_data(contains=timestamps[0]+'_positive')
	# 	folder_b = toolbox.latest_data(contains=timestamps[0]+'_negative')	

	if timestamps[0] == None: 
		folder_a = toolbox.latest_data(contains='positive' + tag)
		folder_b = toolbox.latest_data(contains='negative' + tag)
	elif len(timestamps)==1:		
		folder_b = toolbox.data_from_time(timestamps[0])      
		folder_a = toolbox.latest_data(contains = 'positive', older_than = timestamps[0])   
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
		ax.set_xticklabels(x_labels.tolist())
		ax.set_ylim(-1.1,1.1)
		ax.set_title(str(folder_a)+'/'+str(timestamps[0]))
		ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')

			# print values on bar plot
		def autolabel(rects):
		    for ii,rect in enumerate(rects):
		        height = rect.get_height()
		        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2)) +'('+ str(int(round(y_err[ii]*100))) +')',
		            ha='center', va='bottom')
		# autolabel(rects)

	if save and ax != None:
		try:
		    fig.savefig(
		        os.path.join(folder_a,'tomo.png'))
		    fig.savefig(
		        os.path.join(folder_a,'tomo.pdf'))
		except:
		    print 'Figure has not been saved.'

	if return_data == True:
		return x_labels, x, y, y_err

def BarPlotTomoContrastFull(timestamp = None, state = 'Z', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp =None, save = True,
		plot_fit = True, color = 'r', plot_type = ''):
		### SSRO calibration
	if state != '000init':
		k_range = 21
	else:
		k_range = 9

	for k in range(k_range):
		
		# print k 
		# print 'positive'+ str(k)
		if k < 9:
			if k_range >9:
				timestamp_9, folder_9 = toolbox.latest_data(contains = 'state_' + state +'_positive_'+ str(9), older_than = timestamp,return_timestamp = True)
				timestamp_pos, folder_a = toolbox.latest_data(contains = 'state_' + state +'_positive_'+ str(k), older_than = timestamp_9,return_timestamp = True)
				timestamp_neg, folder_b = toolbox.latest_data(contains = 'state_' + state +'_negative_'+ str(k), older_than = timestamp_9,return_timestamp = True)
			else:
				timestamp_pos, folder_a = toolbox.latest_data(contains =  '_positive_'+ state + str(k), older_than = timestamp,return_timestamp = True)
				timestamp_neg, folder_b = toolbox.latest_data(contains =  '_negative_'+ state + str(k), older_than = timestamp,return_timestamp = True)
		else:
			timestamp_pos, folder_a = toolbox.latest_data(contains = 'state_' + state +'_positive_'+ str(k), older_than = timestamp,return_timestamp = True)
			timestamp_neg, folder_b = toolbox.latest_data(contains = 'state_' + state +'_negative_'+ str(k), older_than = timestamp,return_timestamp = True)

		
		# print folder_a
		# print folder_b


		x_labels_t, x_t, y_t, y_err_t  = BarPlotTomoContrast(timestamps = [timestamp_pos,timestamp_neg], measurement_name = ['adwindata'],folder_name ='Tomo',
								ssro_calib_timestamp =ssro_calib_timestamp, save = False,
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
	# x_labels = x_labels_t
		# print x_labels.tolist
		# print y.tolist

	Tomo = {}
	Tomo['X_list'] = [9,18,27,50,52,58,62]
	Tomo['Y_list'] = [9,18,27,49,53,59,61]
	Tomo['Z_list'] = [0,3,6,9,18,27,36]
	Tomo['mX_list'] = [9,18,27,50,52,58,62]
	Tomo['mY_list'] = [9,18,27,49,53,59,61]
	Tomo['mZ_list'] = [0,3,6,9,18,27,36]
	Tomo['000init_list'] = [2,5,8,17,26,35,62]

	Tomo['X_decoded'] = [62,62,62]
	Tomo['Y_decoded'] = [53,59,61]
	Tomo['Z_decoded'] = [0,3,6]
	Tomo['mX_decoded'] = [62,62,62]
	Tomo['mY_decoded'] = [53,59,61]
	Tomo['mZ_decoded'] = [0,3,6]

	state_tick_list = x_labels

	if plot_fit ==True :
		if plot_type != 'compressed': 
			fig,ax = plt.subplots(figsize=(32,5)) 
			x1 = [x[i] for i in np.linspace(0,8,9).astype(int)]
			y1 = [y[i] for i in np.linspace(0,8,9).astype(int)]
			y_err_1 = [y_err[i] for i in np.linspace(0,8,9).astype(int)]
			x2 = [x[i] for i in np.linspace(9,35,27).astype(int)]
			y2 = [y[i] for i in np.linspace(9,35,27).astype(int)]
			y_err_2 = [y_err[i] for i in np.linspace(9,35,27).astype(int)]
			x3 = [x[i] for i in np.linspace(36,62,27).astype(int)]
			y3 = [y[i] for i in np.linspace(36,62,27).astype(int)]
			y_err_3 = [y_err[i] for i in np.linspace(36,62,27).astype(int)]

			y_fid = [abs(y[i]) for i in Tomo[state+'_list']]

			y_fid_err = [y_err[i]**2 for i in Tomo[state+'_list']]

			y_id = np.zeros(len(y))
			for ii in Tomo[state+'_list']:
				y_id[ii] = np.sign(y[ii])
			print y_id
			y_fidelity = 1/8.*(1+np.sum(y_fid))
			y_fidelity_err = 1/8.*(np.sum(y_fid_err))**0.5

			print state
			print 'FIDELITY'
			print y_fidelity
			print y_fidelity_err

			y_fid = [abs(y[i])/0.93 for i in [9,18,27]]

			y_fid_err = [y_err[i] for i in [9,18,27]]

			y_id = np.zeros(len(y))
			for ii in Tomo[state+'_list']:
				y_id[ii] = np.sign(y[ii])
			print y_id

			y_fidelity_cs1 = 1/4.*(1+y_fid[0]+y_fid[1]+y_fid[2])
			y_fidelity_cs2 = 1/4.*(1-y_fid[0]-y_fid[1]+y_fid[2]) ## Carbon 1
			y_fidelity_cs3 = 1/4.*(1-y_fid[0]+y_fid[1]-y_fid[2]) ## Carbon 2
			y_fidelity_cs4 = 1/4.*(1+y_fid[0]-y_fid[1]-y_fid[2]) ## Carbon 5
			y_fidelity_cs_err = 1/4.*(y_fid_err[0]**2+y_fid_err[1]**2+y_fid_err[2]**2)**0.5

			print state
			print 'CODE SPACE FIDELITY'
			print y_fidelity_cs1, y_fidelity_cs2, y_fidelity_cs3, y_fidelity_cs4
			print y_fidelity_cs_err

			

			# y_fid = [(abs(y[i])+1)/2. for i in Tomo[state+'_decoded']]

			# y_fid_err = [y_err[i]**2/2. for i in Tomo[state+'_decoded']]

			# y_fidelity = 1/3.*(np.sum(y_fid))
			# y_fidelity_err = 1/3.*(np.sum(y_fid_err))**0.5
			
			# if state == 'X' or state =='mX':
			# 	y_fidelity_err = y_err[62]/2.

			# print state
			# print 'DECODED FIDELITY'
			# print y_fidelity
			# print y_fidelity_err

			ax.bar(x,y_id,align ='center',color = color, alpha = 0.2,ecolor = color, linewidth = 0)
			# ax.bar(x1,y1,yerr=y_err_1,align ='center',ecolor = color ,color = color,alpha = 1,linewidth = 1)
			# ax.bar(x2,y2,yerr=y_err_2,align ='center',ecolor = color ,color = color,alpha = 1,linewidth = 1)
			# ax.bar(x3,y3,yerr=y_err_3,align ='center',ecolor = color ,color = color, alpha = 1,linewidth = 1)

			rects = ax.bar(x,y,yerr=y_err,align ='center',ecolor = color ,color = color, alpha = 1,linewidth = 1)
			ax.set_xticks(np.arange(0,63,1))
			# ax.bar(x,y_id,align ='center',color = '0.9')
			ax.set_xticklabels(state_tick_list, rotation = 'vertical')
			ax.set_ylim(-1.1,1.1)
			ax.set_ylim(-1,1)
			ax.set_xlim(-0.5,62.5)
			# ax.set_xticks(x, minor = True)
			yticks = np.linspace(-1,1,5)
			plt.yticks(yticks)
			# rects.set_linewidth(1.5)
			# ax.hlines([0],x[0]-1,x[-1]+1,linestyles='dotted',linewidth = 3)
			ax.tick_params(axis='x', which='major', labelsize=20)
			ax.tick_params(axis='y', which='major', labelsize=30)
			# autolabel(rec4s)
			ax.set_ylabel('Expectation value',fontsize = 30)
			mpl.rcParams['axes.linewidth'] = 1
			mpl.rcParams['pdf.fonttype'] = 42

			def autolabel(rects):
			    for ii,rect in enumerate(rects):
			        height = rect.get_height()
			        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2)) +'('+ str(int(round(y_err[ii]*100))) +')',
			            ha='center', va='bottom')
			# autolabel(rects)
			if save and ax != None:
				print 'SAVED '+ folder
				try:
				    fig.savefig(
				        os.path.join(folder,state+'tomo_full.png'),bbox_inches='tight')
				    fig.savefig(
				        os.path.join(folder,state+'tomo_full.pdf'),bbox_inches='tight')
				except:
				    print 'Figure has not been saved.'
			# return y_fidelity, y_fidelity_err
			return  y_fidelity_cs1, y_fidelity_cs2, y_fidelity_cs3, y_fidelity_cs4, y_fidelity_cs_err

		elif plot_type == 'compressed':
			fig,ax = plt.subplots(figsize=(7,5)) 
			if state == '000init':
				fig,ax = plt.subplots(figsize=(7,5)) 

			plot_list = [0,3,6,9,18,27,53,59,61,62] 
			# [9,18,27,49,53,59,61]
			if state == '000init':
				plot_list = [2,5,8,17,26,35,62]			
			y_id = np.zeros(len(y))
			for ii in Tomo[state+'_list']:
				y_id[ii] = np.sign(y[ii])
			print y_id

			x1 = np.arange(0,len(plot_list),1)
			y1 = [y[i] for i in plot_list]
			y_err_1 = [y_err[i] for i in plot_list]

			y_id_1 = [y_id[i] for i in plot_list]

			state_tick_list_1 = [state_tick_list[i] for i in plot_list]

			y_fid = [abs(y[i]) for i in Tomo[state+'_list']]

			y_fid_err = [y_err[i] for i in Tomo[state+'_list']]



			y_fidelity = 1/8.*(1+np.sum(y_fid))
			y_fidelity_err = 1/8.*(np.sum(y_fid_err))**0.5

			print state
			print 'FIDELITY'
			print y_fidelity
			print y_fidelity_err
			ax.bar(x1,y_id_1,align ='center',color = color, alpha = 0.2,ecolor = color, linewidth = 1,edgecolor = color)
			rects = ax.bar(x1,y1,yerr=y_err_1,align ='center',ecolor = 'k' ,color = color,alpha = 1,linewidth = 1)

			ax.set_xticks(np.arange(0,len(plot_list),1))
			# ax.bar(x,y_id,align ='center',color = '0.9')
			ax.set_xticklabels(state_tick_list_1, rotation = 'vertical')
			

			ax.set_ylim(-1,1)
			if state == '000init':
				ax.set_ylim(0,1)
			ax.set_xlim(-0.5,len(plot_list)-0.5)
			# ax.set_xticks(x, minor = True)
			if state != '000init':
				yticks = np.linspace(-1,1,5)
			else:
				yticks = np.linspace(0,1,2)
			plt.yticks(yticks)
			# rects.set_linewidth(1.5)
			# ax.hlines([0],x[0]-1,x[-1]+1,linestyles='dotted',linewidth = 3)
			ax.tick_params(axis='x', which='major', labelsize=30)
			ax.tick_params(axis='y', which='major', labelsize=30)
			# autolabel(rec4s)
			ax.set_ylabel('Expectation value',fontsize = 30)
			mpl.rcParams['axes.linewidth'] = 1
			mpl.rcParams['pdf.fonttype'] = 42

			
		# ax.tick_params('x', length=3, width=1.5, which='minor')
			# print values on bar plot
			def autolabel(rects):
			    for ii,rect in enumerate(rects):
			        height = rect.get_height()
			        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y1[ii],2)) +'('+ str(int(round(y_err_1[ii]*100))) +')',
			            ha='center', va='bottom')
			# autolabel(rects)
	ax.tick_params('both', length=3, width=1, which='major')
	if save and ax != None:
		print 'SAVED '+ folder
		try:
		    fig.savefig(
		        os.path.join(folder,state+'tomo_full.png'),bbox_inches='tight')
		    fig.savefig(
		        os.path.join(folder,state+'tomo_full.pdf'),bbox_inches='tight')
		except:
		    print 'Figure has not been saved.'

c_grey = (240/255.,242/255.,166/255.)
c_green = (9/255.,232/255.,94/255.)
c_grey = (64/255.,78/255.,77/255.)#(240/255.,242/255.,166/255.)
c_blue = (68/255.,204/255.,255/255.)
c_red = (150/255.,52/255.,132/255.)
c_orange = (240/255.,242/255.,166/255.)
c_orange_2 = (242/255.,129/255.,35/255.)

color_list = ['b',c_green,c_blue,c_red,c_orange_2,'r']

y_fidelity_cs1 = np.zeros(6)
y_fidelity_cs2 = np.zeros(6)
y_fidelity_cs3 = np.zeros(6)
y_fidelity_cs4 = np.zeros(6)
y_fidelity_err = np.zeros(6)

y_fidelity_cs1[0],y_fidelity_cs2[0],y_fidelity_cs3[0], y_fidelity_cs4[0],y_fidelity_err[0] = BarPlotTomoContrastFull(timestamp = '20141225_224000', state = 'Z', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp ='20141225_150151', save = True,
		plot_fit = True,color = c_blue)

y_fidelity_cs1[1],y_fidelity_cs2[1],y_fidelity_cs3[1], y_fidelity_cs4[1],y_fidelity_err[1] = BarPlotTomoContrastFull(timestamp = '20141226_015500', state = 'mZ', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp ='20141225_150151', save = True,
		plot_fit = True,color = c_blue)

y_fidelity_cs1[2],y_fidelity_cs2[2],y_fidelity_cs3[2], y_fidelity_cs4[2],y_fidelity_err[2] = BarPlotTomoContrastFull(timestamp = '20141230_103000', state = 'Y', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp ='20141230_064816', save = True,
		plot_fit = True,color = c_blue)

y_fidelity_cs1[3],y_fidelity_cs2[3],y_fidelity_cs3[3], y_fidelity_cs4[3],y_fidelity_err[3] = BarPlotTomoContrastFull(timestamp = '20141230_144000', state = 'mY', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp ='20141230_064816', save = True,
		plot_fit = True,color = c_blue)

y_fidelity_cs1[4],y_fidelity_cs2[4],y_fidelity_cs3[4], y_fidelity_cs4[4],y_fidelity_err[4] = BarPlotTomoContrastFull(timestamp = '20141226_050300', state = 'X', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp ='20141225_150151', save = True,
		plot_fit = True,color = c_blue)

y_fidelity_cs1[5],y_fidelity_cs2[5],y_fidelity_cs3[5], y_fidelity_cs4[5],y_fidelity_err[5] = BarPlotTomoContrastFull(timestamp = '20141226_090000', state = 'mX', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp ='20141225_150151', save = True,
		plot_fit = True,color = c_blue)

# print sum(y_fidelity_cs1)/6.
# print sum(y_fidelity_cs2)/6.
# print sum(y_fidelity_cs3)/6.
# print sum(y_fidelity_cs4)/6.
# print sum(y_fidelity_err**2)**0.5/6.

BarPlotTomoContrastFull(timestamp = '20141225_190000', state = '000init', measurement_name = ['adwindata'],folder_name ='Tomo',
		ssro_calib_timestamp ='20141225_150151', save = True,
		plot_fit = True,color = c_blue)


# BarPlotTomoContrastFull(timestamp = '20141225_224000', state = 'Z', measurement_name = ['adwindata'],folder_name ='Tomo',
# 		ssro_calib_timestamp ='20141225_150151', save = True,
# 		plot_fit = True,color = c_blue, plot_type = 'compressed')

# # BarPlotTomoContrastFull(timestamp = '20141226_015500', state = 'mZ', measurement_name = ['adwindata'],folder_name ='Tomo',
# # 		ssro_calib_timestamp ='20141225_150151', save = True,
# # 		plot_fit = True,color = c_blue, plot_type = 'compressed')

# BarPlotTomoContrastFull(timestamp = '20141230_103000', state = 'Y', measurement_name = ['adwindata'],folder_name ='Tomo',
# 		ssro_calib_timestamp ='20141230_064816', save = True,
# 		plot_fit = True,color = c_blue, plot_type = 'compressed')

# BarPlotTomoContrastFull(timestamp = '20141230_144000', state = 'mY', measurement_name = ['adwindata'],folder_name ='Tomo',
# 		ssro_calib_timestamp ='20141230_064816', save = True,
# 		plot_fit = True,color = c_blue, plot_type = 'compressed')

# BarPlotTomoContrastFull(timestamp = '20141226_050300', state = 'X', measurement_name = ['adwindata'],folder_name ='Tomo',
# 		ssro_calib_timestamp ='20141225_150151', save = True,
# 		plot_fit = True,color = c_blue, plot_type = 'compressed')

# # BarPlotTomoContrastFull(timestamp = '20141226_090000', state = 'mX', measurement_name = ['adwindata'],folder_name ='Tomo',
# 		ssro_calib_timestamp ='20141225_150151', save = True,
# 		plot_fit = True,color = c_blue, plot_type = 'compressed')

plt.show()