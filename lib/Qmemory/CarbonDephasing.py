import numpy as np
import os,h5py
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common

color_list = ['b','g','y','r','brown','s']

def get_from_hdf5(folder,key_list):
	# gets a msmt_parameter from an hdf5 file
	Datafile=h5py.File(folder+folder[26:] + '.hdf5','r') 

	first_dict_layer = Datafile[Datafile.keys()[0]]
	return_list = []

	for key in key_list:
		return_list = return_list+[first_dict_layer.attrs[key]]
	return return_list

def get_tstamp_from_folder(folder):
	return folder[18:18+15]

def get_dephasing_data(folder_dict,ssro_calib_folder,**kw):

	tomos = kw.pop('tomos',['X','Y'])

	## contains tomo values
	data_dict = {	
	'sweep_pts': [],
	}

	# print tomos
	for t in ['X','Y','XX','XY','YX','YY']:
		data_dict.update({t:[]})
		data_dict.update({t+'_u':[]})

	for t in tomos:
		for i,f in enumerate(folder_dict[t]):
			a = mbi.MBIAnalysis(f)
			a.get_sweep_pts()
			a.get_readout_results(name='adwindata')
			a.get_electron_ROC(ssro_calib_folder)


			
			x_labels = a.sweep_pts.reshape(-1)


			if i == 0:
				data_dict[t] = ((a.p0.reshape(-1))-0.5)*2
				data_dict[t+'_u'] = 2*a.u_p0.reshape(-1)
			else:
				y = ((a.p0.reshape(-1))-0.5)*2
				y_u = 2*a.u_p0.reshape(-1)
				data_dict[t] = [y0/2-y[ii]/2 for ii,y0 in enumerate(data_dict[t])]
				data_dict[t+'_u'] = [np.sqrt(y0**2+y_u[ii]**2)/2 for ii,y0 in enumerate(data_dict[t+'_u'])]

	
	## one carbon experiment
	if len(tomos[0]) ==1:
		npY = np.array(data_dict['Y'])
		npX = np.array(data_dict['X'])
		npY_u = np.array(data_dict['Y_u'])
		npX_u = np.array(data_dict['X_u'])
		return x_labels,npX,npY,npX_u,npY_u

	elif len(tomos[0]) == 2:
		data_dict['sweep_pts'] = x_labels
		return data_dict


def Sweep_repetitions(older_than = None,
		folder_name ='Memory_NoOf_Repetitions_',
		carbon = '2',
		ssro_calib_timestamp =None, 
		plot_result= True, **kw) :

	folder_dict = {
	'sweep_pts': [],
	'res' : [],
	'res_u' : []
	}

	for t in ['X','Y','XX','XY','YX','YY']:
		folder_dict.update({t:[]})

	logicstate = kw.get('logicstate',None)
	return_fits = kw.get('return_fits',False)

	### search data
	if len(carbon) ==1:
		for ro in ['positive','negative']:
			for t in ['X','Y']:
				search_string = folder_name+ro+'_Tomo_'+t+'_'+'C'+carbon
				folder_dict[t].append(toolbox.latest_data(contains = search_string,older_than = older_than,raise_exc = False))
	### two carbons were involved
	elif len(carbon) == 2:
		for ro in ['positive','negative']:
			for t in ['XX','YY','XY','YX']:
				search_string = folder_name+ro+'_state'+logicstate+'_Tomo_'+t+'_'+'C'+carbon
				folder_dict[t].append(toolbox.latest_data(contains = search_string,older_than = older_than,raise_exc = False))

	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSRO')
	else:
		ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
		ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil8'

	if len(carbon) ==1:

		x_labels,npX,npY,npX_u,npY_u = get_dephasing_data(folder_dict,ssro_calib_folder)
		
		folder_dict['res'] = np.sqrt(npY**2+npX**2)
		folder_dict['res_u'] = np.sqrt((npX*npX_u)**2+(npY*npY_u)**2)/np.sqrt((npX**2+npY**2))

	elif len(carbon) ==2:

		tomo_dict = get_dephasing_data(folder_dict,ssro_calib_folder, tomos = ['XX','YY','XY','YX'])
		x_labels = tomo_dict['sweep_pts']
		XX,XX_u = np.array(tomo_dict['XX']),np.array(tomo_dict['XX_u'])
		YY,YY_u = np.array(tomo_dict['YY']),np.array(tomo_dict['YY_u'])
		XY,XY_u = np.array(tomo_dict['XY']),np.array(tomo_dict['XY_u'])
		YX,YX_u = np.array(tomo_dict['YX']),np.array(tomo_dict['YX_u'])

		folder_dict['res'] = np.sqrt(XX**2+YY**2+XY**2+YX**2)/np.sqrt(2) ### geometric sum of possible correlations divided by maximum --> sqrt(2)
		### calculate error bar in the same fashion

		folder_dict['res_u'] = np.sqrt((XX*XX_u)**2+(YY*YY_u)**2+(YX*YX_u)**2+(XY*XY_u)**2)/np.sqrt(XX**2+YY**2+XY**2+YX**2)
		folder_dict['res_u'] = folder_dict['res_u']/np.sqrt(2)
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

	elif not return_fits:

		return np.array(x_labels), np.array(folder_dict['res']),np.array(folder_dict['res_u']),folder_dict[t][0]

	else: 
		### fit a function tot he extracted data. return the fit results and coupling strength (extract from hdf5 file).
		A0 = folder_dict['res'][0]
		offset = 0
		decay = 50

		x0 = 0
		p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)

		fixed = [0,3]

		fit_result = fit.fit1d(x_labels,np.array(folder_dict['res']),None,p0 = p0, fitfunc = fitfunc, do_print = False, ret = True, fixed = fixed)


		### after fitting: get the coupling strength. (one could make this a subroutine at some point.)
		# for the coupling strength we need logic state ms = 0 and ms = +-1 freuqncies.
		coupling = 0
		# print fit_result['params_dict']['tau'], fit_result['error_dict']['tau']


		#extract one complete path
		for key in ['XX','X']:
			if folder_dict[key] != []:
				folder = folder_dict[key][0]

		carbon_list,logic_state = get_from_hdf5(folder,['carbon_list','2qb_logical_state'])
		### extract the expected coupling strength (works for 1 and 2 carbons)
		for ii,c in enumerate(carbon_list):
			if ii == 0:
				[coupling] = get_from_hdf5(folder,['C'+str(c)+'_freq_1_m1'])
				f_m1,f_0 = get_from_hdf5(folder,['C'+str(c)+'_freq_1_m1','C'+str(c)+'_freq_0'])
			else:
				f_m1,f_0 = get_from_hdf5(folder,['C'+str(c)+'_freq_1_m1','C'+str(c)+'_freq_0'])
				if logic_state == 'X':
					# print 'carbon freqs',f_m1,f_0
					coupling += (f_m1-f_0)
				else: 
					# print 'mX'
					# print 'carbon freqs',f_m1,f_0
					coupling +=  - f_m1-f_0

		#calculate Delta_f w.r.t. to f_0
		coupling = abs(abs(coupling)-f_0)

		return coupling, fit_result['params_dict']['tau'],fit_result['error_dict']['tau'],folder

def Sweep_Rep_List(carbons = ['1','2'],older_than = None,ssro_calib_timestamp = None,**kw):


	## other key word arguments
	fit_results = kw.pop('fit_results',True)
	sequence_length = kw.pop('sequence_length',None)
	logicstate_list = kw.pop('logicstate_list',len(carbons)*['X']) ## can be list such as ['X','mX'], used for DFS measurements.


	x_arr = []
	y_arr = []
	y_u_arr = []
	for c,logicstate in zip(carbons,logicstate_list):
		x,y,y_u,folder = Sweep_repetitions(older_than = older_than, 
										carbon = c,
										ssro_calib_timestamp =ssro_calib_timestamp, 
										plot_result= False,logicstate = logicstate,**kw)

		x_arr.append(x)
		y_arr.append(y)
		y_u_arr.append(y_u)
	### convert to time instead of repetitions:
	if sequence_length != None:
		x_arr = [x*sequence_length for x in x_arr]
	
	fig = plt.figure()
	ax = plt.subplot()
	for x,y,y_u,carbon,logicstate,jj in zip(x_arr,y_arr,y_u_arr,carbons,logicstate_list,range(len(x_arr))):
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

		label_txt = 'C'+carbon
		if len(carbon)!=1:
			label_txt = label_txt+'_'+logicstate

		plt.errorbar(x,y,y_u,marker='o',color = color_list[jj],label=label_txt)


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

def coupling_vs_repetitions(c_identifiers,**kw):

	older_than = kw.get('older_than',None)
	s = 0
	for c in c_identifiers:
		s += len(c)
	x = np.zeros(s)
	y = np.zeros(s)
	y_u = np.zeros(s)

	## acquire data
	ii = 0
	for carbon in c_identifiers:
		if len(carbon) > 1:
			for logicstate in ['X','mX']:
				x[ii],y[ii],y_u[ii],folder = Sweep_repetitions(carbon = carbon,logicstate = logicstate,return_fits=True,plot_result = False,**kw)
				ii +=1
		else:
			x[ii],y[ii],y_u[ii],folder = Sweep_repetitions(carbon = carbon,return_fits=True,plot_result = False,**kw)
			ii +=1

	return x,y,y_u,folder
	

def Osci_period(carbon = '1',older_than = None,ssro_calib_timestamp = None,**kw):

	fit_results = kw.pop('fit_results',True)
	folder_name = kw.pop('folder_name','Memory_NoOf_Repetitions_')

	### fit parameters
	freq = kw.pop('freq',1/170.)
	offset = kw.pop('offfset',0.)
	decay = kw.pop('decay',200)
	fixed = kw.pop('fixed',[1])
	show_guess = kw.pop('show_guess',False)

	auto_analysis = kw.get('auto_analysis',False)

	if auto_analysis:
		do_print  = False
		add_txt = True



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
	if auto_analysis:
		tomos = ['X']
	else:
		tomos = ['X','Y']
	for ro in ['positive','negative']:
		for t in tomos:
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
	if auto_analysis:
		res,res_u,rng = [npX],[npX_u],range(1)
	else:
		res,res_u,rng = [npX,npY],[npX_u,npY_u],range(2)
	for y,y_u,jj in zip(res,res_u,range(1)):
		if fit_results:
			A0 = max(y)
			phi0 = 0
			p0,fitfunc,fitfunc_str = common.fit_decaying_cos(freq,offset,A0,phi0,decay)

			# fixed = [1]
			print x_labels
			print y
			fit_result = fit.fit1d(x_labels,y,None,p0 = p0, fitfunc = fitfunc, do_print = do_print, ret = True, fixed = fixed)
			plot.plot_fit1d(fit_result, np.linspace(x_labels[0],x_labels[-1],1001), ax=ax,color = color_list[jj], plot_data=False,add_txt = add_txt, lw = 2)

			if show_guess:
				ax.plot(np.linspace(x_labels[0],x_labels[-1],201), fitfunc(np.linspace(x_labels[0],x_labels[-1],201)), ':', lw=2)
		
		plt.errorbar(x_labels,y,y_u,marker='o',color = color_list[jj],label='C'+carbon+['X','Y'][jj])

	## define folder for data saving
	folder = folder_dict[t][0]

	plt.xlabel('Repump repetitions')
	plt.ylabel('Contrast')
	plt.title(get_tstamp_from_folder(folder) + ' Dephasing for C'+carbon)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.savefig(os.path.join(folder,'CarbonDephasing_osci.pdf'),format='pdf')
	plt.savefig(os.path.join(folder,'CarbonDephasing_osci.png'),format='png')
	if not auto_analysis:
		plt.show()
	plt.close('all')

	print 'Results are saved in ', folder[18:18+15]
	if auto_analysis:
		return fit_result
def get_PosNeg_data(name,**kw):


	"""
	The idea is that we have two folders with the same ending: name
	searches for one folder positive_+name and for negative_+name 
	averages the data and returns the sweep points, measured contrast, uncertainty and the positive folder.
	"""
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
	label = kw.get('label',None)
	y_u = kw.pop('y_u',None)
	if y_u != None:
		plt.errorbar(x,y,y_u,label = label)
	else: plt.plot(x,y)

def fit_exp_pos_neg_data(folder_name,**kw):
	ax = kw.pop('ax',None)
	x,y1,y_u,f = get_PosNeg_data(folder_name,**kw)
	plot_data(x,y1,y_u=y_u,**kw)

	offset,A0,decay,x0 = 0,0.8,400,0
	p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)
	fixed = [0,3]
	fit_result = fit.fit1d(x,y1,None,p0 = p0, fitfunc = fitfunc, do_print = True, ret = True, fixed = fixed)
	plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001),ax=ax, plot_data=False,add_txt=False, lw = 2,**kw)
	return f

def repump_speed_singleExp(timestamp=None, measurement_name = 'adwindata', ssro_calib_timestamp =None,
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            decay_constant = 200, 
            exponent = 2, 
            plot_fit = True, do_print = False, fixed = [2], show_guess = True):
    ''' 
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''
    figsize=(6,4.7)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('ElectronRepump')

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print ssro_calib_folder


    fit_results = []
    
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder)

    ax=a.plot_results_vs_sweepparam(ret='ax',ax=None, 
                                    figsize=figsize, 
                                    ylim=(0.0,1.0)
                                    )

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    # ax.plot(x,y)
    
    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
             x0, decay_constant, exponent)

         #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    ## plot data and fit as function of total time
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    fit_results.append(fit_result)

    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')

    return fit_results

def repump_speed_doubleExp(timestamp=None, measurement_name = 'adwindata', ssro_calib_timestamp =None,
            offset = 0., 
            x0 = 0, 
            amplitude_one = 0.8,
            amplitude_two = 0.2,    
            decay_constant_one = 0.2, 
            decay_constant_two = 0.6, 
            exponent =2,
            plot_fit = True, do_print = False, fixed = [2], show_guess = True):
    ''' 
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''
    figsize=(6,4.7)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('ElectronRepump')

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print ssro_calib_folder


    fit_results = []
    
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder)

    ax=a.plot_results_vs_sweepparam(ret='ax',ax=None, 
                                    figsize=figsize, 
                                    ylim=(0.0,1.0)
                                    )

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    # ax.plot(x,y)
    
    #fitfunction: y(x) = A * exp(-x/tau)+ A2 * exp(-x/tau2) + a
    p0, fitfunc, fitfunc_str = common.fit_double_exp_decay_with_offset(offset, amplitude_one, 
             decay_constant_one, amplitude_two, decay_constant_two)

         #plot the initial guess
    if show_guess:
    	print 'i was here!!!!'
        ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    ## plot data and fit as function of total time
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    fit_results.append(fit_result)

    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')

    return fit_results


    # def fit_double_exponential(self, name='', save=True, log_plot=True, offset=0, **kw):          ret = kw.get('ret', None)
    #     ax = kw.get('ax', None)
    #     indices = kw.pop('indices',range(self.sweep_length))

    #     if ax == None:
    #         fig = self.default_fig(figsize=(6,4))
    #         ax = self.default_ax(fig)
    #     else:
    #         save = False
    #     xx=self.tail_hist_b[:-1]
        
    #     for i in indices:
            
    #         yy=self.tail_hist_h[i]+offset*i
    #         if log_plot:
    #             ax.semilogy(xx,yy,'-')
    #         else:
    #             ax.plot(xx,yy)


    #         guess_xo = 50.
    #         guess_tau = 50.
    #         guess_tau2 = 300.
    #         guess_o = 50.
    #         guess_A = 500.
    #         guess_AA = 500.
    #         AA=fit.Parameter(guess_A, 'A')
    #         A=fit.Parameter(guess_AA, 'AA')
    #         o=fit.Parameter(guess_o, 'o')
    #         xo = fit.Parameter(guess_xo, 'xo')
    #         tau = fit.Parameter(guess_tau, 'tau')
    #         tau2 = fit.Parameter(guess_tau2, 'tau2')
    #         p0 = [A, AA, o, xo, tau, tau2]
    #         #print p0

    #         def fitfunc(x):
    #             #return o() + A() *np.exp(-(x-xo())/tau())
    #             return o() + A() * np.exp(-((x-xo())/tau())) + AA()*np.exp(-((x-xo())/tau2()))
    #         fit_result = fit.fit1d(xx, yy, None, p0=p0, fitfunc=fitfunc, fixed=[], do_print=True, ret=True)
    #         x_fit=np.linspace(xx[0],xx[-1],500)
    #         y_fit=fit_result['fitfunc'](x_fit)
    #         ax.plot(x_fit,y_fit,color='k')



    #     ax.set_xlabel('Time after sync [ns]')
    #     ax.set_ylabel('Counts')

    #     if save:
    #         self.save_fig_incremental_filename(fig,'plot_tail_hist_all')
        
    #     if ret == 'ax':
    #         return ax
    #     if ret == 'fig':
    #         return fig
