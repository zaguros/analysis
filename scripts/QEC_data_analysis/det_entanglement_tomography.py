import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import matplotlib as mpl
from analysis.scripts.QEC import Two_Qubit_Tomography_Postselection as tomo_ps
reload(tomo_ps)
from analysis.lib.QEC import ConditionalParity as CP

def append_data(timestamps_1 = [None,None],timestamps_2 = [],timestamps_3 = [],ssro_folder = None, det = False, ms = None):
	
	y_total = ([])
	y_err_total = ([])
	x_total = ([])
	x_labels_total = ([])
	p0 = ([])
	p1 = ([])
	time_list = [timestamps_1]
	if timestamps_2 != []: 
		time_list = [timestamps_1,timestamps_2]
	if timestamps_3 !=[]:
		time_list = [timestamps_1,timestamps_2,timestamps_3]

	for timestamps in time_list:
		folder_a = toolbox.data_from_time(timestamps[0])      
		folder_b = toolbox.data_from_time(timestamps[1])
		
		if det == False:
			a = mbi.MBIAnalysis(folder_a)
		else:
			a = CP.ConditionalParityAnalysis(folder_a)

		a.get_sweep_pts()
		a.get_readout_results(name='adwindata')
		a.get_electron_ROC(ssro_folder)
		y_a= ((a.p0.reshape(-1)[:])-0.5)*2
		y_err_a = 2*a.u_p0.reshape(-1)[:] 




		if det == False:
			b = mbi.MBIAnalysis(folder_b)
		else:
			b = CP.ConditionalParityAnalysis(folder_b)
		b.get_sweep_pts()
		b.get_readout_results(name='adwindata')
		b.get_electron_ROC(ssro_folder)
		
		y_b= ((b.p0.reshape(-1)[:])-0.5)*2
		y_err_b = 2*b.u_p0.reshape(-1)[:] 


		x_labels = a.sweep_pts.reshape(-1)[:]
		# x = range(len(y_a)) 
		
		y = (y_a - y_b)/2.
		y_err =  1./2*(y_err_a**2 + y_err_b**2)**0.5 




		if ms == 0:
		    a = CP.ConditionalParityAnalysis(folder_a)

		    a.get_sweep_pts()
		    a.get_readout_results(name='adwindata',post_select = True)
		    a.get_electron_ROC(ssro_folder)
		            
		    c0_0_p,u_c0_0_p =  a.convert_fidelity_to_contrast(a.p0_0,a.u_p0_0)
		    c0_1_p,u_c0_1_p =  a.convert_fidelity_to_contrast(a.p0_1,a.u_p0_1)

		    b = CP.ConditionalParityAnalysis(folder_b)
		    b.get_sweep_pts()
		    b.get_readout_results(name='adwindata',post_select = True)
		    b.get_electron_ROC(ssro_folder)
	            
		    c0_0_n,u_c0_0_n =  b.convert_fidelity_to_contrast(b.p0_0,b.u_p0_0)
		    c0_1_n,u_c0_1_n =  b.convert_fidelity_to_contrast(b.p0_1,b.u_p0_1)	

	        ## ms=0 data
		    y = (c0_0_p - c0_0_n)/2.
		    y_err =  1./2*(u_c0_0_p**2 + u_c0_0_n**2)**0.5 
		        
		if ms == 1:
		    a = CP.ConditionalParityAnalysis(folder_a)
		    
		    a.get_sweep_pts()
		    a.get_readout_results(name='adwindata',post_select = True)
		    a.get_electron_ROC(ssro_folder)
		    c0_0_p,u_c0_0_p =  a.convert_fidelity_to_contrast(a.p0_0,a.u_p0_0)
		    c0_1_p,u_c0_1_p =  a.convert_fidelity_to_contrast(a.p0_1,a.u_p0_1)

		    b = CP.ConditionalParityAnalysis(folder_b)
		    b.get_sweep_pts()
		    b.get_readout_results(name='adwindata',post_select = True)

		    b.get_electron_ROC(ssro_folder)
		            
		    c0_0_n,u_c0_0_n =  b.convert_fidelity_to_contrast(b.p0_0,b.u_p0_0)
		    c0_1_n,u_c0_1_n =  b.convert_fidelity_to_contrast(b.p0_1,b.u_p0_1)	
		        ## ms=1 data
		    y = (c0_1_p - c0_1_n)/2.
		    y_err =  1./2*(u_c0_1_p**2 + u_c0_1_n**2)**0.5 

		y_total = np.append(y_total,y)

		if det == True and ms != None:
			print 'PROBABILITIES'
			print 1/2.*(a.p0+b.p0)
			print 1/2.*(a.p1+b.p1)
			p0 = np.append(p0, 1/2.*(a.p0+b.p0))
			p1 = np.append(p1, 1/2.*(a.p1+b.p1))
		
		y_err_total = np.append(y_err_total,y_err)
		# x_total = np.append(x_total,x)
		x_labels_total = np.append(x_labels_total,x_labels)
	if det == True and ms != None:
		print 'PROBABILITIES'
		print p0
		print p1
		print np.average(p0)
		print np.average(p1)
	# print 1/2.*(a.p0+b.p0)

	# if ms ==1:
	# 	sign = [1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,-1]
	# 	print len(sign)
	# 	print len(y_total)
	# 	y_total = [y_total[i]*sign[i] for i in range(len(y_total))]	

	x_total= range(len(y_total))
	return x_total, x_labels_total, y_total,y_err_total

###########################################################



folder 							= r'D:\measuring\data\Analyzed figures\Deterministic Entanglement'
SSRO_timestamp, ssro_folder     = toolbox.latest_data(contains = 'AdwinSSRO', older_than = '20141230_134013',return_timestamp = True)
figure_name_list 				= ['00_init','deterministic_entanglement','probabilistic_entanglement']
folder = r'D:\measuring\data\QEC_data\figs\final figures'

###########################
''' Paper figures'''
###########################
figure_name_list = ['00_init_12','determinstic_entanglement_12','probabilistic_entanglement_12','00_init_15','determinstic_entanglement_15','probabilistic_entanglement_15']
# figure_name_list = ['00_init']
# figure_name_list = ['probabilistic_entanglement']
# figure_name_list = ['determinstic_entanglement_15']

### Carbons 5 and 1
c_grey = (240/255.,242/255.,166/255.)
c_green = (9/255.,232/255.,94/255.)
c_grey = (64/255.,78/255.,77/255.)#(240/255.,242/255.,166/255.)
c_blue = (68/255.,204/255.,255/255.)
c_red = (150/255.,52/255.,132/255.)
c_orange = (240/255.,242/255.,166/255.)
c_orange_2 = (242/255.,129/255.,35/255.)

# color_list = ['b',c_green,c_blue,c_red,c_orange_2,'r']
color_list = [c_green]+2*[c_red]+[c_green]+2*[c_red]
color_list = 6*[c_red]

for ii, figure_name in enumerate(figure_name_list):
	# print figure_name
	c_color = color_list[ii]
	if figure_name == '00_init_12':
		timestamps_1 = ['20141230_154629','20141230_155355']
		timestamps_2 = []
		timestamps_3 = []
		det = False
		ms_list = [None]

	elif figure_name == 'probabilistic_entanglement_12':
		timestamps_1 = ['20150102_162737','20150102_163437']
		timestamps_2 = ['20141230_162443','20141230_163020']
		timestamps_3 = []	
		det = False
		ms_list = [None]
	elif figure_name == 'determinstic_entanglement_12':
		timestamps_1 = ['20141230_161531','20141230_161851']
		timestamps_2 = ['20141230_160309','20141230_160842']	
		timestamps_3 = []	
		det = True
		ms_list = [None,0,1]
	if figure_name == '00_init_15':
		timestamps_1 = ['20141230_180459','20141230_181258']
		timestamps_2 = []
		timestamps_3 = []
		det = False
		ms_list = [None]
	elif figure_name == 'probabilistic_entanglement_15':
		timestamps_1 = ['20150102_161207','20150102_161932']
		timestamps_2 = ['20150102_154713','20150102_155914']
		timestamps_3 = []	
		det = False
		ms_list = [None]
	elif figure_name == 'determinstic_entanglement_15':
		timestamps_1 = ['20141230_183900','20141230_184242']
		timestamps_2 = ['20141230_182409','20141230_182754']	
		timestamps_3 = ['20141230_183209','20141230_183420']		
		det = True
		ms_list = [None,0,1]


	print 'test'
	print figure_name
	print ms_list
	for jj,ms in enumerate(ms_list):
		hatch_list = ['','/','\\']
		print 'jj is ' + str(jj)
		if jj >0:
			hatch = hatch_list[jj]
		else:
			hatch = ''
					# print 'test'
		if 'probabilistic_entanglement' in figure_name:
			hatch = '//'

		print ms
		print hatch
		x, x_labels, y, y_err = append_data(ms = ms, 
								timestamps_1 = timestamps_1,timestamps_2 = timestamps_2,timestamps_3 =timestamps_3, 
								ssro_folder = ssro_folder, det = det)
		fig,ax = plt.subplots() 
		state_tick_list = x_labels

		x1 = [x[i] for i in np.linspace(0,5,6).astype(int)]
		y1 = [y[i] for i in np.linspace(0,5,6).astype(int)]
		y_err_1 = [y_err[i] for i in np.linspace(0,5,6).astype(int)]
		x2 = [x[i] for i in np.linspace(6,14,9).astype(int)]
		y2 = [y[i] for i in np.linspace(6,14,9).astype(int)]
		y_err_2 = [y_err[i] for i in np.linspace(6,14,9).astype(int)]

		ax.bar(x1,y1,yerr=y_err_1,align ='center',ecolor = 'k' ,edgecolor ='k',color = c_color,alpha = 1,linewidth = 1,hatch = hatch)
		ax.bar(x2,y2,yerr=y_err_2,align ='center',ecolor = 'k' ,edgecolor ='k',color = c_color,alpha = 1,linewidth = 1,hatch = hatch)
		# ax.set_xticks(x, minor = True)
		# ax.set_xticks([2,5,14])
		# ax.set_xticklabels(state_tick_list)
		ax.set_xticks(x)
		ax.set_xticklabels(state_tick_list, rotation = 'vertical')
		
		ax.set_ylim(-1,1)
		ax.set_xlim([x[0]-0.6,x[-1]+0.6])
		yticks = np.linspace(-1,1,5)
		plt.yticks(yticks)
		# rects.set_linewidth(1.5)
		# ax.hlines([0],x[0]-1,x[-1]+1,linestyles='dotted',linewidth = 2)
		ax.tick_params(axis='x', which='major', labelsize=25)
		ax.tick_params(axis='y', which='major', labelsize=25)
		ax.tick_params('both', length=5, width=3, which='major')
		ax.tick_params('x', length=3, width=1.5, which='minor')		
		# autolabel(rects)
		ax.set_ylabel('Expectation value',fontsize = 30)
		mpl.rcParams['axes.linewidth'] = 1
		ax.tick_params('both', length=3, width=1, which='major')
		y_id = np.zeros(len(y))
		if '00_init' not in figure_name:
			for ii in [6,10,14]:
				y_id[ii] = np.sign(y[ii])
		elif '00_init' in figure_name:
			for ii in [2,5,14]:
				y_id[ii] = np.sign(y[ii])
		print y_id
		ax.bar(x,y_id,align ='center',color = c_red, alpha = 0.2,ecolor = c_red, linewidth = 0)
		plt.savefig(os.path.join(folder, figure_name+'_ms_'+str(ms)+'_'+timestamps_1[0] + '.pdf'),
		format='pdf',bbox_inches='tight')
		plt.savefig(os.path.join(folder, figure_name+'_ms_'+str(ms)+'_'+timestamps_1[0] + '.png'),
		format='png',bbox_inches='tight')
		print figure_name
		print 'ms ' +str(ms)
		if '00' in figure_name:
			print 'fidelity is '+str(1/4.*(1+y[2]+y[5]+y[-1]))+'+/-'+str(1/4.*(y_err[2]**2+y_err[5]**2+y_err[-1]**2)**0.5)			
		else:
			print 'fidelity is '+str(1/4.*(1+abs(y[6])+abs(y[10])+abs(y[-1])))+'+/-'+str(1/4.*(y_err[6]**2+y_err[10]**2+y_err[-1]**2)**0.5)
		plt.show()


###########################
''' Presentation figures'''
###########################
if 0:
	for ii, figure_name in enumerate(figure_name_list):
		print figure_name

		### Carbons 5 and 1
		if figure_name == '00_init':
			timestamps_1 = ['20141230_180459','20141230_181258']
			timestamps_2 = []
			timestamps_3 = []
			det = False
		elif figure_name == 'probabilistic_entanglement':
			timestamps_1 = ['20150102_161207','20150102_161932']
			timestamps_2 = ['20150102_154713','20150102_155914']
			timestamps_3 = []	
			det = False
		elif figure_name == 'determinstic_entanglement':
			timestamps_1 = ['20141230_183900','20141230_184242']
			timestamps_2 = ['20141230_182409','20141230_182754']	
			timestamps_3 = ['20141230_183209','20141230_183420']
			det = True

		### Carbons 2 and 1
		#for supplement?

		ms_list = [None,0,1]
		hatch_list = ['','+','O']
	for jj,ms in enumerate(ms_list):
		if jj >0:
			hatch = hatch_list[jj]
		else:
			hatch = ''
		fig,ax = plt.subplots(figsize = (15,10)) 
		x, x_labels, y, y_err = append_data(ms = ms, timestamps_1 = timestamps_1,timestamps_2 = timestamps_2,timestamps_3 =timestamps_3, ssro_folder = ssro_folder, det = det)
		state_tick_list = [x_labels.tolist()[i] for i in [2,5,14]]
		x1 = [x[i] for i in np.linspace(0,5,6).astype(int)]
		y1 = [y[i] for i in np.linspace(0,5,6).astype(int)]
		y_err_1 = [y_err[i] for i in np.linspace(0,5,6).astype(int)]
		x2 = [x[i] for i in np.linspace(6,14,9).astype(int)]
		y2 = [y[i] for i in np.linspace(6,14,9).astype(int)]
		y_err_2 = [y_err[i] for i in np.linspace(6,14,9).astype(int)]

		ax.bar(x1,y1,yerr=y_err_1,align ='center',ecolor = c_color ,edgecolor =c_color,color = c_color,alpha = 0.7,linewidth = 1,hatch = hatch)
		ax.bar(x2,y2,yerr=y_err_2,align ='center',ecolor = c_color ,edgecolor =c_color,color = c_color,alpha = 1,linewidth = 1,hatch = hatch)	
		# ax.set_xticks(x, minor = True)
		# ax.set_xticks([2,5,14])
		ax.set_xticks(x)

		# ax.set_xticklabels(state_tick_list)
		ax.set_xticklabels(x_labels)
		ax.set_xticks(x, minor = True)
		ax.set_ylim(-1,1)
		ax.set_xlim([x[0]-0.6,x[-1]+0.6])
		yticks = np.linspace(-1,1,3)
		ax.set_yticks(np.arange(-1,1.1,0.2),minor = True)
		plt.yticks(yticks)

		# rects.set_linewidth(1.5)
		ax.hlines([0],x[0]-1,x[-1]+1,linestyles='dotted',linewidth = 1,color = '0.5')
		ax.tick_params(axis='x', which='major', labelsize=30)
		ax.tick_params(axis='y', which='major', labelsize=30)
		ax.tick_params('both', length=5, width=3, which='major')
		ax.tick_params('y', length=3, width=1.5, which='minor')		
		# autolabel(rects)
		ax.set_ylabel('Contrast',fontsize = 30)
		mpl.rcParams['axes.linewidth'] = 1
		ax.tick_params('both', length=4, width=1, which='major')
		plt.savefig(os.path.join(folder, figure_name+'_'+timestamps_1[0] + '.pdf'),
		format='pdf',bbox_inches='tight')
		plt.savefig(os.path.join(folder, figure_name+'_ms_'+str(ms)+'_'+timestamps_1[0] + '.png'),
		format='png',bbox_inches='tight')

	# tomo_ps.BarPlotTomoContrast(timestamps=['20141230_161531','20141230_161851'], measurement_name = ['adwindata'],folder_name ='Tomo',
	#         post_select = False, ssro_calib_timestamp ='20141230_134012')
