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
import matplotlib
from matplotlib import pyplot as plt
import Zeno_fitting_tools as Zfits
reload(common)
reload(toolbox)
reload(Zfits)

def run():
	data = pickle.load( open( "Zeno_data.p", "rb" ) ) 

	results = [] ### we will store a string of results in this array for each measurement.
	fits = [] ### is used for plotting the fits in the end.
	fig=plt.figure()
	ax=plt.subplot()

	font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10}

	matplotlib.rc('font', **font)


	#########################################################
	### static for all fits --> T2* of carbon 5 		  ###
	#########################################################
	### set the average of the single qubit decay and the 0 msmts decay.
	t = (13.14+11.8)/2.
	# print 17.0/np.sqrt(2)

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
	    plot.plot_fit1d(fit_result_single, np.linspace(evo[0],evo[-1],1001), ax=ax, plot_data=False,color = 'brown',add_txt = False, lw = 1)

	### add single qubit and 0 measurement curve to the fits array.

	fits.append(fit_result_single)
	fits.append(fit_result)

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

	fits.append(fit_result)
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

	fits.append(fit_result)
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

	fits.append(fit_result)
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

	fits.append(fit_result)
	#########################################################
	#########################################################
	#########################################################



	"""
		FIT the 5 msmts curve!
		Use the parameters obtained from the 0 measurements curve to estimate the parameters for 1 measurement.
	"""

	kk = 5
	evo = data[str(kk)][0]
	Fid = data[str(kk)][1]
	Fid_u = data[str(kk)][2]


	p = 0.05
	p0, fitfunc,fitfunc_str = Zfits.fit_5msmt_proc_fid(amp0,offset0,t, p)
	#plot the initial guess

	if False:
	    ax.plot(np.linspace(evo[1],evo[-1],201), fitfunc(np.linspace(evo[1],evo[-1],201)), ':', lw=2)


	fit_result = fit.fit1d(evo,Fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=[0,1,2]) #fix the input amplitude and offset from the 0 msmt run.

	if True:
		plot.plot_fit1d(fit_result, np.linspace(0.0,evo[-1],1001), ax=ax, plot_data=False,color = 'y',add_txt = False, lw = 1)

	p1 = str(round(fit_result['params'][0]*100,1))
	p1_u = str(round(fit_result['error'][0]*100,1))
	results.append('p = '+p1+ ' +- '+ p1_u + '(%)')

	fits.append(fit_result)
	#########################################################
	#########################################################
	#########################################################


	"""
		FIT the 6 msmts curve!
		Use the parameters obtained from the 0 measurements curve to estimate the parameters for 1 measurement.
	"""

	kk = 6
	evo = data[str(kk)][0]
	Fid = data[str(kk)][1]
	Fid_u = data[str(kk)][2]


	p = 0.05
	p0, fitfunc,fitfunc_str = Zfits.fit_6msmt_proc_fid(amp0,offset0,t, p)
	#plot the initial guess

	if False:
	    ax.plot(np.linspace(evo[1],evo[-1],201), fitfunc(np.linspace(evo[1],evo[-1],201)), ':', lw=2)


	fit_result = fit.fit1d(evo,Fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=[0,1,2]) #fix the input amplitude and offset from the 0 msmt run.

	if True:
		plot.plot_fit1d(fit_result, np.linspace(0.0,evo[-1],1001), ax=ax, plot_data=False,color = 'k',add_txt = False, lw = 1)

	p1 = str(round(fit_result['params'][0]*100,1))
	p1_u = str(round(fit_result['error'][0]*100,1))
	results.append('p = '+p1+ ' +- '+ p1_u + '(%)')

	fits.append(fit_result)
	#########################################################
	#########################################################
	#########################################################


	"""
	plot the raw data and label with error probabilities.
	"""

	for i in range(7):
		ax.errorbar(data[str(i)][0],data[str(i)][1],data[str(i)][2],fmt='o',ms=2.5,label = str(i)+' : '+results[i])

	ax.errorbar(data['single'][0],data['single'][1],data['single'][2],fmt='o',color='brown', ms=2.5 ,label ='single qubit')		
	plt.xlabel('Evolution time (ms)')
	plt.ylabel('Process fidelity')
	plt.title('Zeno run')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

	plt.savefig(r'D:\measuring\data\Zeno_results\Zeno_fitted.pdf',format='pdf')
	plt.savefig(r'D:\measuring\data\Zeno_results\Zeno_fitted.png',format='png')
	plt.show()

	plt.close('all')

	#########################################################
	#########################################################
	#########################################################

	### use the fitted results and the data to plot every single run.


	### the font used for labels/ticks and legends
	font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10}

	matplotlib.rc('font', **font)

	f, axarr = plt.subplots(4,2, sharex='col', sharey='row')
	# ((axsingle, ax0, ax1, ax2), (ax3, ax4, ax5, ax6)) = axarr ### this just clears up which ax is affiliated with which measurement.
	color_array = [['brown','b','g','r'],['c','m','y','k']]

	### plot the first column.
	for i in range(3):
		axarr[i+1][0].errorbar(data[str(i)][0],data[str(i)][1],data[str(i)][2],fmt='o',ms=2.5, color = color_array[0][i+1], label = str(i)+' : '+results[i])
		plot.plot_fit1d(fits[i+1], np.linspace(0.0,evo[-1],1001), ax=axarr[i+1][0], plot_data=False,color = color_array[0][i+1], add_txt = False, lw = 1)
		axarr[i+1][0].legend()
		# axarr[i+1].tick_params(axis='both',labelsize=20)
		# axarr[i+1].set_ticks([0.4,0.7])
	### plot the second column
	for i in range(4):
		axarr[i][1].errorbar(data[str(i+3)][0],data[str(i+3)][1],data[str(i+3)][2],fmt='o',ms=2.5, color = color_array[1][i], label = str(i+3)+' : '+results[i+3])
		plot.plot_fit1d(fits[i+4], np.linspace(0.0,evo[-1],1001), ax=axarr[i][1], plot_data=False,color = color_array[1][i], add_txt = False, lw = 1)
		axarr[i][1].legend()

	### plot the single qubit result
	axarr[0][0].errorbar(data['single'][0],data['single'][1],data['single'][2],fmt='o',color='brown', ms=2.5 ,label ='single qubit')
	plot.plot_fit1d(fits[0], np.linspace(0.0,evo[-1],1001), ax=axarr[0][0], plot_data=False,color = color_array[0][0], add_txt = False, lw =1)
	axarr[0][0].legend()
	
	f.subplots_adjust(hspace=0)

	### make x-ticks disappear.
	# plt.setp([a.get_xticklabels() for a in axarr[:-1][0]], visible=False)
	# plt.setp([a.get_xticklabels() for a in axarr[:-1][1]], visible=False)

	f.set_size_inches(11.,8)


	for i in range(4):
		axarr[i][0].set_ylim([0.3,0.95])
		axarr[i][0].set_yticks(np.array([0.4,0.6,0.8]))

	# plt.yticks(np.array([0.4,0.6,0.8]))
	# plt.ylim([0.3,0.95])

	print f.get_size_inches()
	axarr[3][0].set_xlabel('Evolution time (ms)')
	axarr[3][1].set_xlabel('Evolution time (ms)')
	# plt.ylabel('Process fidelity')
	# axarr[0].set_title('Zeno run')
	plt.savefig(r'D:/measuring/data/Zeno_results/Zeno.pdf')
	plt.show()

	plt.close('all')
	
