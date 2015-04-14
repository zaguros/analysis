"""
This file takes the Zeno data from a pickle file and tries to fit the data with a derived set of functions

2015 NK
"""

import numpy as np
import os,re
import h5py
import pickle
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import Zeno_fitting_tools as Zfits
reload(common)
reload(toolbox)
reload(Zfits)

def run():
	data = pickle.load( open( "Zeno_data.p", "rb" ) ) 

	results = [] ### we will store a string of results in this array for each measurement.

	fig=plt.figure()
	ax=plt.subplot()

	#########################################################
	#########################################################
	#########################################################

	"""
	Fit the 0 msmts curve!
	"""
	kk = 0
	evo = data[str(kk)][0]
	Fid = data[str(kk)][1]
	Fid_u = data[str(kk)][2]

	offset = 0.4
	amplitude = 0.42
	decay_constant = 20.

	p0, fitfunc, fitfunc_str = common.fit_gauss(offset, amplitude, 0., decay_constant)


	results.append('p = 0')

	#plot the initial guess
	if False:
	    ax.plot(np.linspace(evo[0],evo[-1],201), fitfunc(np.linspace(evo[0],evo[-1],201)), ':', lw=2)

	fit_result = fit.fit1d(evo,Fid, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[1]) #only fix x0 to be 0. The gaussian is centered.

	print 'fitfunction: '+fitfunc_str
	## plot data and fit as function of total time
	if True:
	    plot.plot_fit1d(fit_result, np.linspace(evo[0],evo[-1],1001), ax=ax, plot_data=False,color = 'b',add_txt = False, lw = 1)


 	#########################################################
 	#########################################################
 	#########################################################


	"""
	Fit the single qubit curve
	"""
	kk = 'single'
	evo = data[str(kk)][0]
	Fid = data[str(kk)][1]
	Fid_u = data[str(kk)][2]

	offset = 0.4
	amplitude = 0.42
	decay_constant = 21.0/np.sqrt(2)
	print 'T2* single_qubit fit: ', decay_constant

	p0, fitfunc, fitfunc_str = common.fit_gauss(offset, amplitude, 0., decay_constant)

	#plot the initial guess
	if False:
	    ax.plot(np.linspace(evo[0],evo[-1],201), fitfunc(np.linspace(evo[0],evo[-1],201)), ':', lw=2)

	fit_result_single = fit.fit1d(evo,Fid, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[1]) #only fix x0 to be 0. The gaussian is centered.

	print 'fitfunction: '+fitfunc_str
	## plot data and fit as function of total time
	if True:
	    plot.plot_fit1d(fit_result_single, np.linspace(evo[0],evo[-1],1001), ax=ax, plot_data=False,color = 'y',add_txt = False, lw = 1)


 	#########################################################
 	#########################################################
 	#########################################################


	"""
		FIT the 1 msmts curve!
		Use the parameters obtained from the 0 measurements curve to estimate the parameters for 1 measurement.
	"""


	kk = 1
	evo = data[str(kk)][0]
	Fid = data[str(kk)][1]
	Fid_u = data[str(kk)][2]


	amp0 = fit_result['params_dict']['A']
	offset0 = fit_result['params_dict']['a']
	s = fit_result['params_dict']['sigma']

	print amp0, offset0,s

	t = 21.0/np.sqrt(2)
	p = 0.05
	p0, fitfunc,fitfunc_str = Zfits.fit_1msmt_proc_fid(amp0,offset0,t, p)
	#plot the initial guess

	if False:
	    ax.plot(np.linspace(evo[1],evo[-1],201), fitfunc(np.linspace(evo[1],evo[-1],201)), ':', lw=2)


	fit_result = fit.fit1d(evo,Fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=[0,1,2]) #fix the input amplitude and offset from the 0 msmt run.

	if True:
		plot.plot_fit1d(fit_result, np.linspace(0.0,evo[-1],1001), ax=ax, plot_data=False,color = 'g',add_txt = False, lw = 1)

	p1 = str(round(fit_result['params'][0]*100,1))
	p1_u = str(round(fit_result['error'][0]*100,1))
	results.append('p = '+p1+ ' +- '+ p1_u + '(%)')
	#########################################################
	#########################################################
	#########################################################


	"""
		FIT the 2 msmts curve!
		Use the parameters obtained from the 0 measurements curve to estimate the parameters for 1 measurement.
	"""
	kk = 2
	evo = data[str(kk)][0]
	Fid = data[str(kk)][1]
	Fid_u = data[str(kk)][2]

	t = 21.0/np.sqrt(2)
	p = 0.05
	p0, fitfunc,fitfunc_str = Zfits.fit_2msmt_proc_fid(amp0,offset0,t, p)
	#plot the initial guess

	if False:
	    ax.plot(np.linspace(evo[1],evo[-1],201), fitfunc(np.linspace(evo[1],evo[-1],201)), ':', lw=2)


	fit_result = fit.fit1d(evo,Fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=[0,1,2]) #fix the input amplitude and offset from the 0 msmt run.

	if True:
		plot.plot_fit1d(fit_result, np.linspace(0.0,evo[-1],1001), ax=ax, plot_data=False,color = 'r',add_txt = False, lw = 1)
	

	p1 = str(round(fit_result['params'][0]*100,1))
	p1_u = str(round(fit_result['error'][0]*100,1))
	results.append('p = '+p1+ ' +- '+ p1_u + '(%)')
	#########################################################
	#########################################################
	#########################################################


	"""
		FIT the 3 msmts curve!
		Use the parameters obtained from the 0 measurements curve to estimate the parameters for 1 measurement.
	"""

	kk = 3
	evo = data[str(kk)][0]
	Fid = data[str(kk)][1]
	Fid_u = data[str(kk)][2]

	t = 21.0/np.sqrt(2)
	p = 0.05
	p0, fitfunc,fitfunc_str = Zfits.fit_3msmt_proc_fid(amp0,offset0,t, p)
	#plot the initial guess

	if False:
	    ax.plot(np.linspace(evo[1],evo[-1],201), fitfunc(np.linspace(evo[1],evo[-1],201)), ':', lw=2)


	fit_result = fit.fit1d(evo,Fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=[0,1,2]) #fix the input amplitude and offset from the 0 msmt run.

	if True:
		plot.plot_fit1d(fit_result, np.linspace(0.0,evo[-1],1001), ax=ax, plot_data=False,color = 'c',add_txt = False, lw = 1)

	p1 = str(round(fit_result['params'][0]*100,1))
	p1_u = str(round(fit_result['error'][0]*100,1))
	results.append('p = '+p1+ ' +- '+ p1_u + '(%)')
	#########################################################
	#########################################################
	#########################################################


	"""
		FIT the 4 msmts curve!
		Use the parameters obtained from the 0 measurements curve to estimate the parameters for 1 measurement.
	"""

	kk = 4
	evo = data[str(kk)][0]
	Fid = data[str(kk)][1]
	Fid_u = data[str(kk)][2]

	t = 21.0/np.sqrt(2)
	p = 0.05
	p0, fitfunc,fitfunc_str = Zfits.fit_4msmt_proc_fid(amp0,offset0,t, p)
	#plot the initial guess

	if False:
	    ax.plot(np.linspace(evo[1],evo[-1],201), fitfunc(np.linspace(evo[1],evo[-1],201)), ':', lw=2)


	fit_result = fit.fit1d(evo,Fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=[0,1,2]) #fix the input amplitude and offset from the 0 msmt run.

	if True:
		plot.plot_fit1d(fit_result, np.linspace(0.0,evo[-1],1001), ax=ax, plot_data=False,color = 'm',add_txt = False, lw = 1)

	p1 = str(round(fit_result['params'][0]*100,1))
	p1_u = str(round(fit_result['error'][0]*100,1))
	results.append('p = '+p1+ ' +- '+ p1_u + '(%)')
	#########################################################
	#########################################################
	#########################################################

	"""
	plot the raw data and label with error probabilities.
	"""

	for i in range(5):
		plt.errorbar(data[str(i)][0],data[str(i)][1],data[str(i)][2],fmt='o',ms=2.5,label = str(i)+' : '+results[i])

	plt.errorbar(data['single'][0],data['single'][1],data['single'][2],fmt='o',ms=2.5,label = 'single qubit')		
	plt.xlabel('Evolution time (ms)')
	plt.ylabel('Process fidelity')
	plt.title('Zeno run')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

	plt.savefig(r'D:\measuring\data\Zeno_results\Zeno_fitted.pdf',format='pdf')
	plt.savefig(r'D:\measuring\data\Zeno_results\Zeno_fitted.png',format='png')
	plt.show()

	plt.close('all')