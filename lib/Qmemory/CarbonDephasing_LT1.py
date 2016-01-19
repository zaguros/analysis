import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common

color_list = ['Crimson','RoyalBlue','b','g','y','r','brown','s']
def fit_cos(x,y,**kw):
    a = kw.pop('a',0)
    A = kw.pop('A',0.8)
    x0 = kw.pop('x0',0)
    T = kw.pop('T',1000)
    n = kw.pop('n',1)
    f = kw.pop('a',1/200.)
    phi = kw.pop('a',0)
    fixed=kw.pop('fixed',[0,2,4])
    show_guess = kw.pop('show_guess',False)
    do_print = kw.pop('do_print',False)
    p0, fitfunc, fitfunc_str = common.fit_general_exponential_dec_cos(a,A,x0,T,n,f,phi)
    if show_guess:
        ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=fixed)
    
    return fit_result

def fit_gauss(x,y,**kw):
    a = kw.pop('a',0.1)
    A = kw.pop('A',0.7)
    x0 = kw.pop('x0',0)
    T = kw.pop('T',1000)
    n = kw.pop('n',2)
    fixed=kw.pop('fixed',[0,2,4])
    show_guess = kw.pop('show_guess',False)
    do_print = kw.pop('do_print',False)
    p0, fitfunc, fitfunc_str = common.fit_general_exponential(a,A,x0,T,n)
    if show_guess:
        ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=fixed)
    
    return fit_result
def get_tstamp_from_folder(folder):
	return folder[18:18+15]

def get_contrast_data(folders,ssro_calib_folder=None):
	'''
	General function to get positive and negative readout data and return the resulting contrast and sweep pts.
	folders = list with folder names of positive and negative readout
	'''
	if ssro_calib_folder==None:
		ssro_calib_folder=toolbox.latest_data('SSRO')

	
	for i,f in enumerate(folders):
		
		a = mbi.MBIAnalysis(f)
		a.get_sweep_pts()
		a.get_readout_results(name='adwindata')
		a.get_electron_ROC(ssro_calib_folder)

		

		swp_pts = a.sweep_pts.reshape(-1)
		# x_labels = x_labels[:-1]
		if i == 0: # this is the positive result
			pos_res = ((a.p0.reshape(-1))-0.5)*2
			pos_res_u = 2*a.u_p0.reshape(-1)
		else:
			neg_res = ((a.p0.reshape(-1))-0.5)*2
			neg_res_u = 2*a.u_p0.reshape(-1)
			
	data_dict={}
	data_dict['Signal']=pos_res/2.-neg_res/2.
	data_dict['Signal_u'] = np.sqrt(neg_res_u**2+pos_res_u**2)/2.
	data_dict['Signal_sum']=(pos_res/2.+neg_res/2.+1)
	data_dict['swp_pts'] = swp_pts
	data_dict['swp_name'] = a.sweep_name
	data_dict['y_axis_name'] = 'Contrast'

	data_dict['plt_title'] = a.timestamp+'\n'+a.measurementstring
	data_dict['msmnt_class'] = a
	#return swp_pts,Signal,Signal_u,Signal_sum
	return data_dict

def plot_and_save_contrast_data(data_dict,**kw):
	label = kw.get('label', None)
	ret = kw.get('ret', None)
	ylim = kw.get('ylim', (-1.05, 1.05))
	ax = kw.get('ax', None)
	fmt = kw.get('fmt', 'o')
	figsize=kw.get('figsize',(6,4))
	markersize = kw.get('markersize',6)
	capsize = kw.get('capsize',3)
	color=kw.get('color','RoyalBlue')
	modes=kw.get('modes',['pdf','png'])
	do_plot=kw.get('do_plot',False)
	do_save=kw.get('do_save',True)
	ret_ax=kw.get('ret_ax',False)
	if ax == None:
		fig = plt.figure()
		ax = plt.subplot()
	ax.errorbar(data_dict['swp_pts'], data_dict['Signal'], fmt='o',
	                    yerr=data_dict['Signal_u'], label=label,markersize=markersize,capsize=capsize,color=color)
	ax.set_xlabel(data_dict['swp_name'])
	ax.set_ylabel(data_dict['y_axis_name'])
	#ax.set_ylim(ylim)
	ax.set_title(data_dict['plt_title'])
	if do_save!= False:
		for mode in modes:
			fig.savefig(os.path.join(data_dict['folders'][0], 'contrast_vs_sweepparam.'+mode),
	                    format=mode)
			fig.savefig(os.path.join(data_dict['folders'][1], 'contrast_vs_sweepparam.'+mode),
	                    format=mode)
	if do_plot== False:
		plt.close(fig)
	if ret_ax:
		return ax

def plot_posneg_sum(data_dict,**kw):
	label = kw.get('labels', None)
	ret = kw.get('ret', None)
	ylim = kw.get('ylim', (-1.05, 1.05))
	ax = kw.get('ax', None)
	fmt = kw.get('fmt', 'o')
	figsize=kw.get('figsize',(6,4))
	markersize = kw.get('markersize',6)
	capsize = kw.get('capsize',3)
	color=kw.get('color','RoyalBlue')
	modes=kw.get('modes',['pdf','png'])
	do_plot=kw.get('do_plot',False)
	do_save=kw.get('do_save',True)
	if ax == None:
		fig = plt.figure()
		ax = plt.subplot()
	ax.errorbar(data_dict['swp_pts'], data_dict['Signal_sum'], fmt='o-',
	                    yerr=data_dict['Signal_u'], label=label,markersize=markersize,capsize=capsize,color=color)
	ax.set_xlabel(data_dict['swp_name'])
	ax.set_ylabel('Sum positive and negative readout')
	ax.set_ylim(ylim)
	ax.set_title(data_dict['plt_title'])
	if do_save!= False:
		for mode in modes:
			fig.savefig(os.path.join(data_dict['folders'][0], 'posneg_sum_vs_sweepparam.'+mode),
	                    format=mode)
			fig.savefig(os.path.join(data_dict['folders'][1], 'posneg_sum_vs_sweepparam.'+mode),
	                    format=mode)
	if do_plot== False:
		plt.close(fig)

def analyze_and_plot_tomo(**kw):
	pre_fix_name = kw.pop('pre_fix_name','Memory_NoOf_Repetitions_')
	post_fix_name = kw.pop('post_fix_name','')
	tomo_list = kw.pop('tomo_list',['X','Y','Z'])
	older_than = kw.pop('older_than',None)
	do_plot_and_save = kw.pop('do_plot_and_save',True)
	f_dict=make_folder_dict(pre_fix_name=pre_fix_name,post_fix_name=post_fix_name,tomo_list=tomo_list,older_than=older_than)
	d_dict=get_tomo_data(f_dict)
	if do_plot_and_save:
		for basis in d_dict.keys():
			plot_and_save_contrast_data(f_dict[basis])
			plot_posneg_sum(f_dict[basis])

	return f_dict
def get_dephasing_data(folder_dict,**kw):
	'''

	'''
	data_dict = {
	'X' : [],
	'Y' : [],
	'resX' : [],
	'resX_u' : [],
	'sum_resX' : [],
	'resY' : [],
	'resY_u': [],
	'sum_resY' : [],
	'sweep_pts': [],
	'res' : [],
	'res_u' : []
	}


	for t in ['X','Y']:
		d = get_contrast_data(folder_dict[t])
		#print x_labels
		x_labels=d['swp_pts']
		data_dict['res'+t]=d['Signal']
		data_dict['res'+t+'_u']=d['Signal_u']

		

	npY = np.array(data_dict['resY'])
	npX = np.array(data_dict['resX'])
	npY_u = np.array(data_dict['resY_u'])
	npX_u = np.array(data_dict['resX_u'])

	return x_labels,npX,npY,npX_u,npY_u

def get_tomo_data(folder_dict,**kw):
	'''
	sweeps over keys in folder dict (list with strings containing msmnt basis e.g. ['X','Y','Z'])
	gets contrast data for pos and neg RO which are supplied in folder_dict['X'] as a list.
	and adds them to folder_dict['X'] as a dictionary with all relevant RO parameters
	'''

	for key in folder_dict.keys():
		#print key

		d = get_contrast_data(folder_dict[key]['folders'],folder_dict[key]['ssro_calib_folder'])
		folder_dict[key].update(d)
		folder_dict[key].update({'y_axis_name':'<'+key+'>'})

	return folder_dict
def stitch_tomo_data(post_fix_name_list,**kw):
	tomo_list = kw.pop('tomo_list',['X','Y'])
	older_list = kw.pop('older_list',[None]*len(post_fix_name_list))
	do_plot_and_save = kw.pop('do_plot_and_save',True)
	combined_dict={}
	for t in tomo_list:
		combined_dict[t] = {'Signal':np.array([]),'Signal_sum':np.array([]),'Signal_u':np.array([]),'swp_pts':np.array([])}
	for i,older in enumerate(older_list):
		combined_dict[i]={}
		combined_dict[i].update(analyze_and_plot_tomo(older_than=older_list[i],post_fix_name=post_fix_name_list[i],tomo_list=tomo_list,do_plot_and_save=do_plot_and_save))
		for tomo in tomo_list:	
			for k in combined_dict[tomo].keys():
				#print combined_dict[tomo][k],combined_dict[i][tomo][k]
				combined_dict[tomo][k]=np.append(combined_dict[tomo][k],combined_dict[i][tomo][k])
			combined_dict[tomo]['swp_name']=combined_dict[i][tomo]['swp_name']
			combined_dict[tomo]['y_axis_name']=combined_dict[i][tomo]['y_axis_name']
			combined_dict[tomo]['plt_title']=combined_dict[i][tomo]['plt_title']
	return combined_dict    
def make_folder_dict(pre_fix_name='Memory_NoOf_Repetitions_',**kw):
	'''
	combine all tomography data of a desired dataset, based on a
	name : string that uniquely defines the measurement you want to get
	tomo_list: list of strings containing the readout basis
	carbon: string defining the carbon Number
	older_than: string 'yyyymmdd_hhmmss' to exlude all data after this time
	ssro_calib_timestamp: string with timestamp if you want to use specific ssro calib
	'''

	tomo_list = kw.pop('tomo_list',['X','Y'])
	carbon = kw.pop('carbon','1')
	older_than = kw.pop('older_than',None)
	ssro_calib_timestamp = kw.pop('ssro_calib_timestamp',None)
	post_fix_name = kw.pop('post_fix_name','')

	folder_dict={}
	
	for t in tomo_list:
		folder_dict[t]={}
		folder_dict[t]['folders']=[]
		for ro in ['positive','negative']:
			search_string = pre_fix_name+ro+'_Tomo_'+t+'_'+'C'+carbon+post_fix_name
			f=toolbox.latest_data(contains = search_string,older_than = older_than,raise_exc = False)
			folder_dict[t]['folders'].append(f)

		if ssro_calib_timestamp == None: 
		    ssro_calib_folder = toolbox.latest_data('SSRO')
		else:
			ssro_calib_folder = toolbox.latest_data(contains=ssro_calib_timestamp+'_AdwinSSRO_SSROCalibration')
		folder_dict[t]['ssro_calib_folder']=ssro_calib_folder

	return folder_dict

def plot_and_fit_dephasing_curves(post_fix_names,older_list,tomo_list,ax,label,color,fit_xy=False,fit_qsum=True):
	    comb_dict=stitch_tomo_data([post_fix_names]*len(older_list),older_list=older_list,tomo_list=tomo_list,do_plot_and_save=False)

	    xlabel='X'
	    ylabel='Y'
	    zlabel='Z'
	    blabel='Sum of pos and neg RO'
	    qlabel='quadratic sum'
	    label=label

	    color_dict={'X':'Crimson','Y':'RoyalBlue'}

	    # this code is for fitting and plotting X,Y,Z in 1 graph
	    #cd.plot_and_save_contrast_data(comb_dict['Y'],do_plot=do_plot,ax=ax,do_save=False,label=ylabel,color=color_dict['Y'])
	    #cd.plot_and_save_contrast_data(comb_dict['X'],do_plot=do_plot,ax=ax,do_save=False,label=xlabel,color=color_dict['X'])
	    #cd.plot_and_save_contrast_data(comb_dict['Z'],do_plot=do_plot,ax=ax,do_save=False,label=zlabel,color='Grey')
	    
	    if fit_xy:
	        for t in ['X','Y']:
	            x=comb_dict[t]['swp_pts']
	            y=comb_dict[t]['Signal']
	            fit_result=fit_cos(x,y)
	            x_fit=np.linspace(x[0],x[-1],500)
	            y_fit=fit_result['fitfunc'](x_fit)
	            ax.plot(x_fit,y_fit,color=color_dict[t])    

	    q_sum=      np.sqrt(
	                #np.array(comb_dict['Z']['Signal'])**2
	                +np.array(comb_dict['X']['Signal'])**2
	                +np.array(comb_dict['Y']['Signal'])**2)
	    

	    ax.errorbar(comb_dict['Y']['swp_pts'],q_sum,
	                yerr=comb_dict['Y']['Signal_u'],color=color,label=label,fmt='o')

	    ## To plot the sum of the X signal instead of the contrast (check for systematic errors like ionisation)
	    #ax.errorbar(comb_dict['X']['swp_pts'],comb_dict['X']['Signal_sum'],yerr=comb_dict['X']['Signal_u'],color=color_list[j],fmt='x',label=blabel)
	    if fit_qsum:
	        y=q_sum
	        x=comb_dict['Y']['swp_pts']
	        fit_result=fit_gauss(comb_dict['Y']['swp_pts'],q_sum,do_print=False)

	        x_fit=np.linspace(x[0],5000,500)
	        y_fit=fit_result['fitfunc'](x_fit)
	        ax.plot(x_fit,y_fit,color=color)    

	    return ax,fit_result

def Sweep_repetitions(older_than = None,
		folder_name ='Memory_NoOf_Repetitions_',
		carbon = '2',
		ssro_calib_timestamp =None, 
		plot_result= True) :

	folder_dict = {
	'X' : [],
	'Y' : [],
	'resX' : [],
	'resX_u' : [],
	'resY' : [],
	'resY_u': [],
	'sweep_pts': [],
	'res' : [],
	'res_u' : []
	}

	### search data

	for ro in ['positive','negative']:
		for t in ['X','Y']:
			search_string = folder_name+ro+'_Tomo_'+t+'_'+'C'+carbon
			folder_dict[t].append(toolbox.latest_data(contains = search_string,older_than = older_than,raise_exc = False))


	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSRO')
	else:
		ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
		ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil8'
		# print ssro_calib_folder


	x_labels,npX,npY,npX_u,npY_u = get_dephasing_data(folder_dict,ssro_calib_folder)

	folder_dict['res'] = np.sqrt(npY**2+npX**2)
	folder_dict['res_u'] = np.sqrt((npX*npX_u)**2+(npY*npY_u)**2)/np.sqrt((npX**2+npY**2))

	if folder_name == 'Memory_Sweep_repump_time_':
		plot_label = 'average repump time'
		
	elif folder_name == 'Memory_NoOf_Repetitions_':
		plot_label = 'Number of repetitions'

	elif folder_name == 'Memory_Sweep_repump_duration':
		plot_label = 'repump duration'

	elif folder_name == 'Memory_sweep_timing_':
		x_labels = np.array(x_labels)*1e6
		plot_label = 'Waiting time (us)'

	else:
		plot_label = folder_name
		

	if plot_result:
		fig = plt.figure()
		ax = plt.subplot()
		plt.errorbar(x_labels,folder_dict['res'],folder_dict['res_u'],marker='o',label='C'+carbon)


		plt.xlabel(plot_label)
		plt.ylabel('Bloch vector length')
		plt.title('Dephasing for C'+carbon)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.savefig(os.path.join(folder_dict[t][0],'CarbonDephasing.pdf'),format='pdf')
		plt.savefig(os.path.join(folder_dict[t][0],'CarbonDephasing.png'),format='png')
		plt.show()
		plt.close('all')

	else:
		return np.array(x_labels), np.array(folder_dict['res']),np.array(folder_dict['res_u']),folder_dict[t][0]

def Sweep_Rep_List(carbons = ['1','2'],older_than = None,ssro_calib_timestamp = None,**kw):


	## other key word arguments
	fit_results = kw.pop('fit_results',True)
	sequence_length = kw.pop('sequence_length',None)
	plot_ylog_scale = kw.pop('plot_ylog_scale',False)
	plot_xlog_scale = kw.pop('plot_xlog_scale',False)

	x_arr = []
	y_arr = []
	y_u_arr = []
	for c in carbons:
		x,y,y_u,folder = Sweep_repetitions(older_than = older_than, 
										carbon = c,
										ssro_calib_timestamp =ssro_calib_timestamp, 
										plot_result= False)

		x_arr.append(x)
		y_arr.append(y)
		y_u_arr.append(y_u)
	### convert to time instead of repetitions:
	if sequence_length != None:
		x_arr = [x*sequence_length for x in x_arr]
	
	fig = plt.figure()
	ax = plt.subplot()
	for x,y,y_u,carbon,jj in zip(x_arr,y_arr,y_u_arr,carbons,range(len(x_arr))):
		if fit_results:
			A0 = y[0]
			offset = 0
			decay = 50
			if sequence_length != None:
				decay = decay*sequence_length
			x0 = 0
			p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)

			fixed = [0,3]

			fit_result = fit.fit1d(x,y,None,p0 = p0, fitfunc = fitfunc, do_print = True, ret = True, fixed = fixed)
			plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax,color = color_list[jj], plot_data=False,add_txt = False, lw = 2)

		plt.errorbar(x,y,y_u,marker='o',color = color_list[jj],label='C'+carbon)

	if plot_ylog_scale:
		ax.set_yscale('log')	
	if plot_xlog_scale:
		ax.set_xscale('log')

	plt.xlabel('Repump repetitions')

	if sequence_length != None:
		plt.xlabel('elapsed time (us)')
	
	plt.ylabel('Bloch vector length')
	plt.title(get_tstamp_from_folder(folder) + ' Dephasing for C'+carbon)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.savefig(os.path.join(folder,'CarbonDephasing.pdf'),format='pdf')
	plt.savefig(os.path.join(folder,'CarbonDephasing.png'),format='png')
	plt.show()
	plt.close('all')



def Osci_period(carbon = '1',older_than = None,ssro_calib_timestamp = None,**kw):

	fit_results = kw.pop('fit_results',True)
	folder_name = kw.pop('folder_name','Memory_NoOf_Repetitions_')

	### fit parameters
	freq = kw.pop('freq',1/170.)
	offset = kw.pop('offfset',0.)
	decay = kw.pop('decay',200)
	fixed = kw.pop('fixed',[1])
	print folder_name



	folder_dict = {
	'X' : [],
	'Y' : [],
	'resX' : [],
	'resX_u' : [],
	'resY' : [],
	'resY_u': [],
	'sweep_pts': [],
	'res' : [],
	'res_u' : []
	}

	### search data

	for ro in ['positive','negative']:
		for t in ['X','Y']:
			search_string = folder_name+ro+'_Tomo_'+t+'_'+'C'+carbon
			folder_dict[t].append(toolbox.latest_data(contains = search_string,older_than = older_than,raise_exc = False))
	
	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSRO')
	else:
		ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
		ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil8'
		# print ssro_calib_folder

		### extract data
	x_labels,npX,npY,npX_u,npY_u = get_dephasing_data(folder_dict,ssro_calib_folder)


	fig = plt.figure()
	ax = plt.subplot()
	for y,y_u,jj in zip([npX,npY],[npX_u,npY_u],range(2)):
		if fit_results:
			A0 = max(y)
			phi0 = 0
			p0,fitfunc,fitfunc_str = common.fit_decaying_cos(freq,offset,A0,phi0,decay)

			# fixed = [1]

			fit_result = fit.fit1d(x_labels,y,None,p0 = p0, fitfunc = fitfunc, do_print = True, ret = True, fixed = fixed)
			plot.plot_fit1d(fit_result, np.linspace(x_labels[0],x_labels[-1],1001), ax=ax,color = color_list[jj], plot_data=False,add_txt = False, lw = 2)

		plt.errorbar(x_labels,y,y_u,marker='o',color = color_list[jj],label='C'+carbon+['X','Y'][jj])

	## define folder for data saving
	folder = folder_dict[t][0]

	plt.xlabel('Repump repetitions')
	plt.ylabel('Contrast')
	plt.title(get_tstamp_from_folder(folder) + ' Dephasing for C'+carbon)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.savefig(os.path.join(folder,'CarbonDephasing_osci.pdf'),format='pdf')
	plt.savefig(os.path.join(folder,'CarbonDephasing_osci.png'),format='png')
	plt.show()
	plt.close('all')
	print 'Results are saved in ', folder[18:18+15]

def get_PosNeg_data(name,**kw):

	ssro_calib_timestamp = kw.pop('ssro_calib_timestamp',None)
	older_than = kw.pop('older_than',None)

	data_dict = {
	'folders' : [],
	'sweep_pts': [],
	'res' : [],
	'res_u' : []
	}


	for ro in ['positive','negative']:
		search_string = ro+name
		data_dict['folders'].append(toolbox.latest_data(contains = search_string,older_than = older_than,raise_exc = False))


	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSRO')
	else:
		ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
		ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil8'
		# print ssro_calib_folder
	for i,f in enumerate(data_dict['folders']):
		a = mbi.MBIAnalysis(f)
		a.get_sweep_pts()
		a.get_readout_results(name='adwindata')
		a.get_electron_ROC(ssro_calib_folder)

		

		x_labels = a.sweep_pts.reshape(-1)
		if i == 0:
			data_dict['res'] = ((a.p0.reshape(-1))-0.5)*2
			data_dict['res_u'] = 2*a.u_p0.reshape(-1)
		else:
			y = ((a.p0.reshape(-1))-0.5)*2
			y_u = 2*a.u_p0.reshape(-1)
			data_dict['res'] = [y0/2-y[ii]/2 for ii,y0 in enumerate(data_dict['res'])]
			data_dict['res_u'] = [np.sqrt(y0**2+y_u[ii]**2)/2 for ii,y0 in enumerate(data_dict['res_u'])]


	return x_labels,data_dict['res'],data_dict['res_u'],data_dict['folders'][0]

def plot_data(x,y,**kw):

	y_u = kw.pop('y_u',None)
	if y_u != None:
		plt.errorbar(x,y,y_u)
	else: plt.plot(x,y)