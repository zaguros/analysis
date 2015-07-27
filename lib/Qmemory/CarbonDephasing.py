import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common

color_list = ['b','g','y','r','brown','s']

def get_dephasing_data(folder_dict,ssro_calib_folder,**kw):

	data_dict = {
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


	for t in ['X','Y']:
		for i,f in enumerate(folder_dict[t]):
			a = mbi.MBIAnalysis(f)
			a.get_sweep_pts()
			a.get_readout_results(name='adwindata')
			a.get_electron_ROC(ssro_calib_folder)

			

			x_labels = a.sweep_pts.reshape(-1)
			# x_labels = x_labels[:-1]
			if i == 0:
				data_dict['res'+t] = ((a.p0.reshape(-1))-0.5)*2
				data_dict['res'+t+'_u'] = 2*a.u_p0.reshape(-1)
			else:
				y = ((a.p0.reshape(-1))-0.5)*2
				y_u = 2*a.u_p0.reshape(-1)
				data_dict['res'+t] = [y0/2-y[ii]/2 for ii,y0 in enumerate(data_dict['res'+t])]
				data_dict['res'+t+'_u'] = [np.sqrt(y0**2+y_u[ii]**2)/2 for ii,y0 in enumerate(data_dict['res'+t+'_u'])]

	npY = np.array(data_dict['resY'])
	npX = np.array(data_dict['resX'])
	npY_u = np.array(data_dict['resY_u'])
	npX_u = np.array(data_dict['resX_u'])

	return x_labels,npX,npY,npX_u,npY_u


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
	if plot_result:
		fig = plt.figure()
		ax = plt.subplot()
		plt.plot(x_labels,folder_dict['res'],marker='o',label='C'+carbon)


		plt.xlabel(plot_label)
		plt.ylabel('Bloch vector length')
		plt.title('Dephasing for C'+carbon)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.savefig(os.path.join(folder_dict[t][0],'CarbonDephasing.pdf'),format='pdf')
		plt.savefig(os.path.join(folder_dict[t][0],'CarbonDephasing.png'),format='png')
		plt.show()
		plt.close('all')

	else:
		return x_labels, np.array(folder_dict['res']),np.array(folder_dict['res_u']),folder_dict[t][0]

def Sweep_Rep_List(carbons = ['1','2'],older_than = None,ssro_calib_timestamp = None,**kw):


	## other key word arguments
	fit_results = kw.pop('fit_results',True)


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
	fig = plt.figure()
	ax = plt.subplot()
	for x,y,y_u,carbon,jj in zip(x_arr,y_arr,y_u_arr,carbons,range(len(x_arr))):
		if fit_results:
			A0 = y[0]
			offset = 0
			decay = 50
			x0 = 0
			p0,fitfunc,fitfunc_str = common.fit_exp_decay_shifted_with_offset(offset,A0,decay,x0)

			fixed = [3]

			fit_result = fit.fit1d(x,y,None,p0 = p0, fitfunc = fitfunc, do_print = True, ret = True, fixed = fixed)
			plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax,color = color_list[jj], plot_data=False,add_txt = False, lw = 2)

		plt.errorbar(x,y,y_u,marker='o',color = color_list[jj],label='C'+carbon)


	plt.xlabel('Repump repetitions')
	plt.ylabel('Bloch vector length')
	plt.title('Dephasing for C'+carbon)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.savefig(os.path.join(folder,'CarbonDephasing.pdf'),format='pdf')
	plt.savefig(os.path.join(folder,'CarbonDephasing.png'),format='png')
	plt.show()
	plt.close('all')
	print 'Results are saved in ', folder


def Osci_period(carbon = '1',older_than = None,ssro_calib_timestamp = None,**kw):

	fit_results = kw.pop('fit_results',True)
	folder_name = kw.pop('folder_name','Memory_NoOf_Repetitions_')
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
			freq = 1/170.
			A0 = max(y)
			offset = 0
			decay = 170
			phi0 = 0
			p0,fitfunc,fitfunc_str = common.fit_decaying_cos(freq,offset,A0,phi0,decay)

			fixed = [1]

			fit_result = fit.fit1d(x_labels,y,None,p0 = p0, fitfunc = fitfunc, do_print = True, ret = True, fixed = fixed)
			plot.plot_fit1d(fit_result, np.linspace(x_labels[0],x_labels[-1],1001), ax=ax,color = color_list[jj], plot_data=False,add_txt = False, lw = 2)

		plt.errorbar(x_labels,y,y_u,marker='o',color = color_list[jj],label='C'+carbon+['X','Y'][jj])

	## define folder for data saving
	folder = folder_dict[t][0]

	plt.xlabel('Repump repetitions')
	plt.ylabel('Contrast')
	plt.title('Dephasing for C'+carbon)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.savefig(os.path.join(folder,'CarbonDephasing_osci.pdf'),format='pdf')
	plt.savefig(os.path.join(folder,'CarbonDephasing_osci.png'),format='png')
	plt.show()
	plt.close('all')
	print 'Results are saved in ', folder

