import numpy as np
import os,re
import h5py
import pickle
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import ssro
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import analysis.scripts.Zeno.Zeno_fitting_tools as Zfits
reload(Zfits)
reload(common)
reload(toolbox)

"""
Analyze the zeno data
NK 2014
"""

def calculate_max_fidelity(fid,fid_u,fid_ort,fid_u_ort):
	"""
	calculate the maximum achievable fidelity by quadratically adding two orthogonal vectors.
	input: 	2 lists of fidelities (fid and fid_ort)
			2 lists of uncertainties for these fidelity values (fid_u and fid_u_ort)

	output:
		fid_bloch maximum obtainable fidelity (list)
		fid_u_bloch corresponding uncertainty (list)
	"""

	### add both data sets together such that the Bloch-vector length in the XY plane is obtained. (geometrically)
	### it is assumed that both measurements form a consistent data set and have the same evolution times.
	fid_bloch=[]
	fid_u_bloch=[]
	for ii in range(len(fid)):
		fid_bloch.append(np.sqrt((fid[ii]-0.5)**2+(fid_ort[ii]-0.5)**2)+0.5)

		### for a derivation of the error propagation see the Zeno onenote file. 2015-04-12 NK

		err1 = (0.25+(fid[ii]-1)*fid[ii])*fid_u[ii]**2
		err2 = (0.25+(fid_ort[ii]-1)*fid_ort[ii])*fid_u_ort[ii]**2
		denominator = 0.5 + fid[ii]*(fid[ii]-1) + fid_ort[ii]*(fid_ort[ii]-1)
		err = np.sqrt((err1+err2)/denominator)
		fid_u_bloch.append(err)

	return fid_bloch, fid_u_bloch


def get_Zeno_data(electron_RO=['positive'], 
				msmts='0',
				state='Z',ROBasis='IX',previous_evo=None, older_than=None,ssro_timestamp=None,single_qubit=False):
	"""
	this function finds timestamps according to an input which specifies the folder name
	Input: 	electron_RO 		is a list with e.g. positive or negative. 
								If the list is longer than 1 element then the 
								contrast values are returned.
			msmts & state 		both strings are used to find the correct folder
			newer_than 			Timestamp which gives the starting point for the search 

			previous_evo		Gives the last free evolution time of the previously extracted data

	Output: timestamp			specifies the latest evaluated folder
			loop_bit			Boolean which signals the end of data acquisition and if the output data should be evaluated.
			x_labels 			Tomographic bases. E.g. XI or ZZ
			y 					Read-out results
			y_err 				Read-out uncertainty
			evo_time 			Free evolution of the specific state
			folder 				The folder of the latest timestamp
			evo_was_zero		Boolean which signifies if the first evolution time was Zero. This will then stop further data acquisition.
	"""

	if not single_qubit:
		search_string=electron_RO[0]+'_logicState_'+state+'_'+msmts+'msmt__ROBasis_'+ROBasis
	
	else:# adjust the measurement string for single qubit msmts. (no entanglement)
		search_string=electron_RO[0]+'_logicState_'+state+'_'+msmts+'msmt_singleQubit_ROBasis_'+ROBasis
	loop_bit = False
	evo_was_zero = False

	timestamp=None

	if previous_evo == None:
		previous_evo = 2000 #very high value such that all evolution times are taken into account.

	#if the desired data set exists, then read the measured values.
	if toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=older_than,
									newer_than=None,
									raise_exc=False) != False and previous_evo!=0:

		timestamp,folder=toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=older_than,
									newer_than=None,
									raise_exc=False)

		evotime,y,y_err= Zeno_get_2Q_values(timestamp,ssro_calib_timestamp=ssro_timestamp)

		x_labels = folder[folder.find('ROBasis_')+8:]
	

	else: #if the searched data does not exist --> dummy results. evotime has to be larger than previous_evo
		x_labels,y,y_err,evotime= 'II',0,0,[2001]

	if evotime[0] < previous_evo or previous_evo==None:
		loop_bit = True
		
	else:
		x_labels,y,y_err,evotime= 'II',0,0,[2001]


	#if positive and negative RO are considered then adjust the search string and search for the same evo time.
	
	if len(electron_RO)>1 and (evotime[0] < previous_evo or previous_evo==None) and loop_bit:
		if not single_qubit:
			search_string=electron_RO[1]+'_logicState_'+state+'_'+msmts+'msmt__ROBasis_'+str(ROBasis)
		else:
			search_string=electron_RO[1]+'_logicState_'+state+'_'+msmts+'msmt_singleQubit_ROBasis_'+str(ROBasis) # adjust the measurement string for single qubit values.
		timestamp2,folder2=toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=older_than,
									newer_than=timestamp, #choose the right search direction.
									raise_exc=False)
		if evotime[0] < previous_evo or previous_evo==None:
			loop_bit = True

			evotime2,y2,y_err2= Zeno_get_2Q_values(timestamp2,ssro_calib_timestamp=ssro_timestamp)
		else:
			x_labels,y2,y_err2,evotime2= 'II',0,0,[2001]

		if electron_RO[0]== 'positive':	
			for i in range(len(y)):
				y[i]=(y[i]-y2[i])/2
				y_err[i]=((y_err[i]**2+y_err2[i]**2)**0.5)/2

		if electron_RO[0] == 'negative':
			for i in range(len(y)):
				y[i]=(-y[i]+y2[i])/2
				y_err[i]=((y_err[i]**2+y_err2[i]**2)**0.5)/2

	#if the evotime is set to zero then the data point needs to be shifted by the length of the applied parity measurements.
	#this duration is stored in the dictionary of the hdf group and is extracted here.
	if evotime[0]==0:
		# open the file, the folder/folder combination gives the complete filepath. 'r' for read only.
		Datafile=h5py.File(folder+folder[26:] + '.hdf5','r') 
		# you need to open the first datagroup. this can be accessed like a dictionary. 
		# that means file.keys[0] returns the first keyword of the dictionary
		first_dict_layer=Datafile[Datafile.keys()[0]]
		#the contents of msmt_params are dumped into the attributes object of this layer.
		#Attributes are again organized in dictionaries. --> Extract the desired value with a keyword
		if u'parity_duration' in first_dict_layer.attrs.keys():
			evotime[0]=first_dict_layer.attrs['parity_duration']*1e3
		Datafile.close()
		evo_was_zero=True

		


	#determine the older timestamp (for the two eRO possiblities) and return that one.
	if loop_bit:
		if len(electron_RO)==1:
			if electron_RO[0]=='negative':
				return timestamp,loop_bit,x_labels,-1*y,y_err,evotime,folder,evo_was_zero
			else:
				return timestamp,loop_bit,x_labels,y,y_err,evotime,folder,evo_was_zero

		elif toolbox.is_older(timestamp,timestamp2):
			return timestamp,loop_bit,x_labels,y,y_err,evotime,folder,evo_was_zero

		else: 
			return timestamp2,loop_bit,x_labels,y,y_err,evotime,folder2,evo_was_zero
	else:
		return older_than,loop_bit,x_labels,y,y_err,evotime,toolbox.data_from_time(older_than),evo_was_zero


def analyze_tests(older_than=None,newer_than=None,ssro_timestamp=None,N=1):
	"""
	Analyzes the test states. In a time window which is specified by older_than and newer_than in timestamp format.
	Analyzes a maximum of N states.
	"""

	N=2*N #we always perform a contrast measurement: positiv and negative --> Two folders.

	search_string='_1msmts_TESTSTATE_ZZ'

	x_arr=[]
	y_arr=[]
	y_u_arr=[]
	old_Arr=[] #collection of relevant timestamps
	fig=plt.figure()
	ax=plt.subplot		

	for i in range(N):

		#if the desired data set exists, then read the measured values.
		if toolbox.latest_data(contains=search_string,
										return_timestamp =True,
										older_than=older_than,
										newer_than=newer_than,
										raise_exc=False) != False:

			older_than,folder=toolbox.latest_data(contains=search_string,
										return_timestamp =True,
										older_than=older_than,
										newer_than=newer_than,
										raise_exc=False)

			evotime,y,y_err= Zeno_get_2Q_values(older_than,ssro_calib_timestamp=ssro_timestamp)
			x_arr.append(i)
			y_arr.extend(y)
			y_u_arr.extend(y_err)
			if i%2==0: #append only one timestamp. we always plot the contrast of positive and negative.
				old_Arr.append(older_than)

	#condense the results for positive and negative read_out
	for i in range(len(x_arr)/2+1):
		if i==0:
			pass
		else:
			y_arr[i-1]=y_arr[2*i-1]/2.-y_arr[2*i-2]/2.
			y_u_arr[i-1]=np.sqrt(y_u_arr[2*i-1]**2+y_u_arr[2*i-2]**2)/2.

	y_arr=y_arr[:len(x_arr)/2]
	y_u_arr=y_u_arr[:len(x_arr)/2]
	x_arr=x_arr[:len(x_arr)/2]
	for i in range(len(y_arr)):
		plt.errorbar(x_arr[i],y_arr[i],y_u_arr[i],marker='o',label=old_Arr[i])
	plt.xlabel('timestamp')
	plt.ylabel('ZZ contrast')
	plt.title('test states')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.xlim(-0.5,N/2-0.5)
	print folder
	plt.savefig(os.path.join(folder,'Zeno_test_states.pdf'),format='pdf')
	plt.savefig(os.path.join(folder,'Zeno_test_states.png'),format='png')
	plt.show()
	plt.close('all')



def Zeno_get_2Q_values(timestamp=None, folder=None,folder_name='Zeno',
						measurement_name = ['adwindata'], 
						ssro_calib_timestamp =None):
	"""
	Returns the relevant RO values for a given timestamp.
	"""

	if timestamp == None and folder==None:
	    timestamp, folder   = toolbox.latest_data(folder_name,return_timestamp =True)
	elif timestamp ==None and folder!=None: 
		pass
	else:
	    folder = toolbox.data_from_time(timestamp)

	if folder != None and timestamp == None:
		d,t = toolbox.get_date_time_string_from_folder(folder)
		timestamp = toolbox.timestamp_from_datetime(t)
	### SSRO calibration
	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSRO',older_than=timestamp)


	else:
	    ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
	    ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'

	a = mbi.MBIAnalysis(folder)
	a.get_sweep_pts()
	a.get_readout_results(name='adwindata')
	try: 
		a.get_electron_ROC(ssro_calib_folder)
	except IOError:
		ssro.ssrocalib(ssro_calib_folder)
		a.get_electron_ROC(ssro_calib_folder)

	x_labels = a.sweep_pts.reshape(-1)
	y= ((a.p0.reshape(-1))-0.5)*2
	y_err = 2*a.u_p0.reshape(-1)

	return x_labels,y,y_err


def Zeno_state_fidelity(older_than_tstamp=None,msmts='0',eRO_list=['positive'],
								 	state='X',
								 	plot_results=True,decoded_bit='2',
								 	ssro_timestamp=None,single_qubit=False, single_qubit_ort = False):
	"""
	Plots or returns the state fidelity for a decoded qubit as a function of time (one parity expectation value)
	"""

	loop_bit = True
	evo_time=None

	y_arr=[]
	y_err_arr=[]
	evo_time_arr=[]
	ii=0

	Tomo1Dict={'X':'ZZ','mX':'ZZ',
	    'Y':'YZ',
	    'mY':'YZ',
	    'Z':'XI',
	    'mZ':'XI'}

	Tomo2Dict={'X':'ZZ','mX':'ZZ',
	    'Y':'ZY',
	    'mY':'ZY',
	    'Z':'IX',
	    'mZ':'IX'}


	if single_qubit:
		### choose orthogonal bases to see the effect of a detuning
		if single_qubit_ort:
			Tomo1Dict={'X':'ZI','mX':'ZI',
			    'Y':'YI',
			    'mY':'YI',
			    'Z':'XI',
			    'mZ':'XI'}
			Tomo2Dict={'X':'IZ','mX':'IZ',
			    'Y':'IX',
			    'mY':'IX',
			    'Z':'IY',
			    'mZ':'IY'}
		### the regular tomography dictionaries for the expected bases.
		else:
			Tomo1Dict={'X':'ZI','mX':'ZI',
			    'Y':'YI',
			    'mY':'YI',
			    'Z':'XI',
			    'mZ':'XI'}
			Tomo2Dict={'X':'IZ','mX':'IZ',
			    'Y':'IY',
			    'mY':'IY',
			    'Z':'IX',
			    'mZ':'IX'}

	RODict={'1':Tomo1Dict,'2':Tomo2Dict}
	evo_time=[2005]
	while loop_bit:
		older_than_tstamp,loop_bit,x_labels,y,y_err,evo_time,folder,evo_was_zero=get_Zeno_data(electron_RO=eRO_list,
																					state=state,
																					older_than=older_than_tstamp,
																					previous_evo=evo_time[0],
																					msmts=msmts,
																					ssro_timestamp=ssro_timestamp,ROBasis=RODict[str(decoded_bit)][state],
																					single_qubit=single_qubit)
		#loop_bit is true as long as new data was found.
		if loop_bit:
			y_arr=np.concatenate((y_arr,y))
			y_err_arr=np.concatenate((y_err_arr,y_err))
			evo_time_arr=np.concatenate((evo_time_arr,evo_time))
			# if an evolution time of zero was reached --> stop the evaluation.
			if evo_was_zero:
				loop_bit=False



	#select the correct expectation value and the right sign for the contrast.
	sign=1

	if 'Y' in state:
		sign=-1
	elif state=='Z':
		sign=1
	elif state=='X':
		sign=1
	
	if 'm' in state:
		sign=-1*sign

	fid_arr=(sign*np.array(y_arr)+1)/2.
	fid_u_arr=np.array(y_err_arr)/2.

	if len(eRO_list)==1:
		RO_String=eRO_list[0]
	else: RO_String = 'contrast'

	fid_arr=fid_arr[np.argsort(evo_time_arr)]
	fid_u_arr=fid_u_arr[np.argsort(evo_time_arr)]
	evo_time_arr=np.sort(evo_time_arr)


	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()

		plt.errorbar(evo_time_arr,fid_arr,fid_u_arr,color='blue',marker='o')
		plt.xlabel('Evolution time (ms)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('logicState_'+state+'_stop_'+str(older_than_tstamp)+'_'+RO_String)

		plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'_'+RO_String+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'_'+RO_String+'.png'),format='png')
		plt.show()
		plt.close('all')
		print folder
	else:
		return evo_time_arr,fid_arr,fid_u_arr,older_than_tstamp,folder

def Zeno_proc_fidelity(msmts='0',
								eRO_list=['positive'],
								older_than_tstamp=None,
								plot_results=True,decoded_bit='2',
								ssro_timestamp=None,single_qubit=False):
	"""
	Plots the process fidelity for a decoded qubit as a function of time

	If plot_Results = False:
		Returns the AVERAGE state fidelity!
	"""	
	state_list=['X','mX','Y','mY','Z','mZ']


	fid_arr=[]
	fid_u_arr=[]

	#get individual state fidelities
	for state in state_list:
		evo_time,fid,fid_u,tstamp,folder=Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
										eRO_list=eRO_list,
										state=state, msmts=msmts,
										plot_results=False,decoded_bit=decoded_bit,
										ssro_timestamp=ssro_timestamp,single_qubit=single_qubit)
		if single_qubit and ('Y' in state or 'Z' in state):
			### get the orthogonal measurements in.
			evo_time_ort,fid_ort,fid_u_ort,tstamp,folder=Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
											eRO_list=eRO_list,
											state=state, msmts=msmts,
											plot_results=False,decoded_bit=decoded_bit,
											ssro_timestamp=ssro_timestamp,single_qubit=single_qubit , single_qubit_ort=True)

			### get maximum fidelity.
			fid, fid_u = calculate_max_fidelity(fid, fid_u, fid_ort, fid_u_ort)
		### append the state fidelities to the fidelity arrays.
		fid_arr.append(fid);fid_u_arr.append(fid_u)

	#calculate average state fidelity
	avg_fid=np.zeros(len(fid_arr[0]))
	avg_fid_u=np.zeros(len(fid_u_arr[0]))
	for i in range(len(fid_arr[0])):
		for ii in range(len(fid_arr)):
			avg_fid[i]= avg_fid[i] + fid_arr[ii][i]/len(fid_arr)
			avg_fid_u[i]= avg_fid_u[i] + fid_u_arr[ii][i]**2/len(state_list)**2
		avg_fid_u[i]=avg_fid_u[i]**0.5


	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot		

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		plt.errorbar(evo_time,(3*avg_fid-1)/2.,1.5*avg_fid_u,color='blue',marker='o')
		plt.xlabel('Evolution time (ms)')
		plt.ylabel('Process fidelity')
		plt.title('Process_fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)

		print folder
		plt.savefig(os.path.join(folder,'Zeno1QProcFid_decBit'+str(decoded_bit)+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QProcFid_decBit'+str(decoded_bit)+'.png'),format='png')
		plt.show()
		plt.close('all')

		fig=plt.figure()
		ax=plt.subplot()

		for i,state in enumerate(state_list):
			plt.errorbar(evo_time,fid_arr[i],fid_u_arr[i],marker='o', label=state)

		plt.xlabel('Evolution time (ms)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('State_fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)
		plt.legend()

		plt.savefig(os.path.join(folder,'Zeno1QStateFidelities_decBit'+str(decoded_bit)+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QStateFidelities_decBit'+str(decoded_bit)+'.png'),format='png')
		plt.show()
		plt.close('all')

	else:
		return evo_time,avg_fid,avg_fid_u,tstamp,folder
def fit_process_decay(msmts,ax,A0,evotime,fid):
	"""
	takes a zeno data set for a specific number of measurements and returns the a fit to the data.
	Inputs:

	msmts 			string which signifies the number of measurements
	A0 				the amplitude for the 0 measurement case.
	evotime 		np.array with the evolution times of the data points
	fid				np.array with the measured state fidelities.

	Output:

	result_string 	a string which is used for labelling the fits
	fit_result		the fitted function for plotting.
	"""

	t = 21./np.sqrt(2)
	p = 0.08
	offset0 = 0.40

	if msmts == '0':
		p0, fitfunc, fitfunc_str = common.fit_gauss(0.5, 0.43, 0., t)

	elif msmts == '1':
		p0, fitfunc,fitfunc_str = Zfits.fit_1msmt_proc_fid(A0,offset0,t, p)

		### manual option to show the intial guess in the plot.
		if False:
			ax.plot(np.linspace(0.,120.0,201), fitfunc(np.linspace(0.,120.0,201)), ':', lw=2)

	elif msmts == '2':
		p0, fitfunc,fitfunc_str = Zfits.fit_2msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '3':
		p0, fitfunc,fitfunc_str = Zfits.fit_3msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '4':
		p0, fitfunc,fitfunc_str = Zfits.fit_4msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '5':
		p0, fitfunc,fitfunc_str = Zfits.fit_5msmt_proc_fid(A0,offset0,t, p)

	elif msmts == '6':
		p0, fitfunc,fitfunc_str = Zfits.fit_6msmt_proc_fid(A0,offset0,t, p)

	### msmts = 0 is an exception
	fixed = [1]
	if msmts =='0':
		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
	else:
		fixed = [0,1] ### fixed parameter: [0,1,2] --> fix decay time, offset and amplitude, p is the only free parameter. 
										###[0,1] --> fix the amplitude and the offset for 0 measurements only.

		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)

	p1 = str(round(fit_result['params'][-1]*100,1))
	p1_u = str(round(fit_result['error'][-1]*100,1))

	result_string = p1 + ' +- ' + p1_u 

	return fit_result, result_string

def Zeno_proc_list(older_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive'],decoded_bit='2',ssro_timestamp=None,single_qubit=False,plot_results=True,fitting=False):

	"""
	Plots the process fidelities as a function of time for a list of msmts

	If plot_results == False --> Returns the arrays for evolution times, PROCESS fidelity and uncertainty.

	"""
	fid_arr,fid=[],[]
	fid_u_arr,fid_u=[],[]
	evotime_arr,evotime=[],[]
	if len(msmt_list)==0:
		print 'nothing to do here'

	else:
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder = Zeno_proc_fidelity(older_than_tstamp=older_than_tstamp,
									eRO_list=eRO_list, msmts=msmt_list[i],ssro_timestamp=ssro_timestamp,decoded_bit=decoded_bit,
									plot_results=False)
			evotime_arr.append(np.sort(evotime))
			fid_arr.append((3*fid[np.argsort(evotime)]-1)/2)
			fid_u_arr.append(1.5*fid_u[np.argsort(evotime)])
		
		### rescale evolution times to ms if it is given in seconds.

		if evotime_arr[-1][-2] < 1.0:
			for kk, timings in enumerate(evotime_arr):
				new_evo = []
				for jj in timings:
					if jj > 1.0:
						new_evo.append(jj)
					else:
						new_evo.append(jj*1e3)

				evotime_arr[kk] = np.array(new_evo)

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		if plot_results:
			fig=plt.figure()
			ax=plt.subplot()
			
			if fitting:
				color_list = ['b','g','r','c','m']

				result = ['0']*len(msmt_list) ### prepare the result strings.

				amp0 = 0.44
				offset0 = 0.40

				t = 21.0/np.sqrt(2)
				p = 0.09
				
				for ii,msmts in enumerate(msmt_list):
					fit_result, result[ii] = fit_process_decay(msmts,ax,amp0,evotime_arr[ii],fid_arr[ii])

					plot.plot_fit1d(fit_result, np.linspace(0.0,100.0,1001), ax=ax, plot_data=False,color = color_list[ii],add_txt = False, lw = 1)

					if msmts == '0':
						result[ii] = ' p = 0'

				for i in range(len(msmt_list)):
					plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],fmt='o',markersize=4,label=str(msmt_list[i]) + ' : p = ' + result[i])

			else: ### no fitting involved.
				for i in range(len(msmt_list)):
					plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],fmt='o',markersize=4,label=str(msmt_list[i])+ ' msmts')

			if single_qubit: # adds the latest single qubit measurement to the data
				evotime_single,fid_single,fid_u_single,tstamp,folder = Zeno_proc_fidelity(older_than_tstamp=older_than_tstamp,
										eRO_list=eRO_list, msmts='0',ssro_timestamp=ssro_timestamp,decoded_bit=decoded_bit,
										plot_results=False,single_qubit=True)
				plt.errorbar(np.array([t*1e3 for t in np.sort(evotime_single)]),(3*fid_single[np.argsort(evotime_single)]-1)/2,1.5*fid_u_single[np.argsort(evotime_single)],marker='o',markersize=4,label='1 qubit')

			plt.xlabel('Evolution time (ms)')
			plt.ylabel('Process fidelity')
			plt.title('Process fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)
			plt.legend()

			print 'Plots are saved in:'
			print folder
			plt.savefig(os.path.join(folder,'ZenoProc_decBit'+str(decoded_bit)+RO_String+'_combined.pdf'),format='pdf')
			plt.savefig(os.path.join(folder,'ZenoProc_decBit'+str(decoded_bit)+RO_String+'_combined.png'),format='png')
			plt.show()
			plt.close('all')
		else: 
			return evotime_arr,fid_arr,fid_u_arr,folder

def fit_State_decay(msmts,ax,A0,evotime,fid):
	"""
	takes a zeno data set for a specific state and number of measurements and returns the fitted result.
	Inputs:

	msmts 			string which signifies the number of measurements
	A0 				the amplitude for the 0 measurement case.
	evotime 		np.array with the evolution times of the data points
	fid				np.array with the measured state fidelities.

	Output:

	result_string 	a string which is used for labelling the fits
	fit_result		the fitted function for plotting.
	"""

	t = 18.2/np.sqrt(2)
	p = 0.08

	if msmts == '0':
		p0, fitfunc, fitfunc_str = common.fit_gauss(0.5, 0.43, 0., t)

	elif msmts == '1':
		p0, fitfunc,fitfunc_str = Zfits.fit_1msmt_state_fid(A0,t, p)

		### manual option to show the intial guess in the plot.
		if False:
			ax.plot(np.linspace(0.,120.0,201), fitfunc(np.linspace(0.,120.0,201)), ':', lw=2)

	elif msmts == '2':
		p0, fitfunc,fitfunc_str = Zfits.fit_2msmt_state_fid(A0,t, p)

	elif msmts == '3':
		p0, fitfunc,fitfunc_str = Zfits.fit_3msmt_state_fid(A0,t, p)

	elif msmts == '4':
		p0, fitfunc,fitfunc_str = Zfits.fit_4msmt_state_fid(A0,t, p)

	elif msmts == '5':
		p0, fitfunc,fitfunc_str = Zfits.fit_5msmt_state_fid(A0,t, p)

	elif msmts == '6':
		p0, fitfunc,fitfunc_str = Zfits.fit_6msmt_state_fid(A0,t, p)

	### msmts = 0 is an exception
	fixed = [0,1]
	if msmts =='0':
		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
	else:
		fixed = [0,1] ### fixed parameter: [0,1] --> fix time and amplitude, p is the only free parameter. 
										###[0] --> fix the amplitude for 0 measurements only.

		fit_result = fit.fit1d(evotime,fid, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)

	p1 = str(round(fit_result['params'][2-len(fixed)]*100,1))
	p1_u = str(round(fit_result['error'][2-len(fixed)]*100,1))

	result_string = p1 + ' +- ' + p1_u 

	return fit_result, result_string

def Zeno_state_list(older_than_tstamp=None,
						msmt_list=['0'],
						eRO_list=['positive'],
						state='Z',
						ssro_timestamp=None,
						decoded_bit='2',
						single_qubit=False,
						fitting = False):
	
	fid=[]
	fid_u=[]
	evotime=[]
	if len(eRO_list)==0:
		print 'nothing to do here'

	else:
		evotime_arr = []
		fid_arr = []
		fid_u_arr = []

		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder = Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
															msmts=msmt_list[i],
															eRO_list=eRO_list,decoded_bit=decoded_bit,
									plot_results=False,state=state,ssro_timestamp=ssro_timestamp)
			evotime_arr.append(evotime[np.argsort(evotime)]); fid_arr.append(fid[np.argsort(evotime)]); fid_u_arr.append(fid_u[np.argsort(evotime)])
		
		if fitting:
			color_list = ['b','g','r','c','m']

			### plot the fits to the data
			results = []
			for kk,m in enumerate(msmt_list):
				if m == '0':

					fit_result, result_string = fit_State_decay(m,ax,0.43,evotime_arr[kk],fid_arr[kk])
					plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = 1)
					results.append('p = 0')

				else:
					if 'Z' in state:
						amp0 = 0.8647
					elif 'Y' in state:
						amp0 = 0.383 * 2
					fit_result, result_string = fit_State_decay(m,ax,amp0,evotime_arr[kk],fid_arr[kk])
					plot.plot_fit1d(fit_result, np.linspace(0.0,120.0,1001), ax=ax, plot_data=False,color = color_list[kk],add_txt = False, lw = 1)
					results.append(result_string)

			### plot the data itself and incorporate the fit results in the labels.
			for i in range(len(msmt_list)):
				plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],fmt='o',label=str(msmt_list[i])+ ' : ' + results[i])
		else:
			### no fitting. plot the measured results
			for i in range(len(msmt_list)):
				plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],marker='o',label=str(msmt_list[i])+' msmts')

		#### treat the case of one qubit decaying.

		if single_qubit:
			evotime,fid,fid_u,tstamp,folder = Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
															msmts='0',
															eRO_list=eRO_list,decoded_bit=decoded_bit,
									plot_results=False,state=state,ssro_timestamp=ssro_timestamp,single_qubit=True)
			### sort arrays.
			fid = fid[np.argsort(evotime)]
			fid_u = fid_u[np.argsort(evotime)]
			### plot
			plt.errorbar(np.array([t*1e3 for t in np.sort(evotime)]),fid,fid_u,marker='o',label='1 qubit')
			if 'Y' in state or 'Z' in state:
				### acquire the data for the orthogonal state.
				evotime_ort,fid_ort,fid_u_ort,tstamp,folder = Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
																msmts='0',
																eRO_list=eRO_list,decoded_bit=decoded_bit,
										plot_results=False,state=state,ssro_timestamp=ssro_timestamp,single_qubit=True,single_qubit_ort = True)
				### sort arrays:
				fid_ort = fid_ort[np.argsort(evotime_ort)]
				fid_u_ort[np.argsort(evotime_ort)]
				plt.errorbar(np.array([t*1e3 for t in np.sort(evotime_ort)]),fid_ort,fid_u_ort,marker='o',label='1 qubit orthogonal')

				fid_bloch, fid_u_bloch = calculate_max_fidelity(fid,fid_u,fid_ort,fid_u_ort)

				### plot the maximum achievable fidelity (by quadratically adding the tomography results in X & Y)	
				plt.errorbar(np.array([t*1e3 for t in np.sort(evotime_ort)]),fid_bloch,fid_u_bloch,marker='o',label='1 qubit maximum fidelity')
		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		plt.xlabel('Evolution time (ms)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('logicState_'+state+'_stop_'+str(tstamp)+'_'+RO_String)
		plt.legend()

		print 'Plots are saved in:'
		print folder
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.png'),format='png')
		plt.show()
		plt.close('all')

def Zen_Compare_Runs(msmt='0',older_timestamp_list=[None],eRO_list=['positive'],
						decoded_bit='2',ssro_timestamp_list=[None],
						single_qubit=False,plot_results=True):
	"""
	Takes a list of 'older_than' timestamps and searches for a run with the specified number of measurements.
	Raw data is then plotted.

	If plot_results = False: 	Returns the evolution time as a 1D list and a list of lists for the process fidelities and their uncertainties. 
								Additionally returns the latest folder.

	"""

	fid_arr,fid=[],[]
	fid_u_arr,fid_u=[],[]
	evotime_arr,evotime=[],[]


	for i in range(len(older_timestamp_list)):
		evotime,fid,fid_u,folder = Zeno_proc_list(older_than_tstamp=older_timestamp_list[i],
								eRO_list=eRO_list, msmt_list=[msmt],ssro_timestamp=ssro_timestamp_list[i],decoded_bit=decoded_bit,
								plot_results=False)
		evotime_arr.append(evotime[0])
		fid_arr.append(fid[0])
		fid_u_arr.append(fid_u[0])
	
	if len(eRO_list)==1:
		RO_String=eRO_list[0]
	else: RO_String = 'contrast'

	if plot_results:
		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(older_timestamp_list)):
			plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],marker='o',markersize=4,label=str(older_timestamp_list[i]))
		
		plt.xlabel('Evolution time (ms)')
		plt.ylabel('Process fidelity')
		plt.title('Process fidelity comparison'+'_'+RO_String+'_msmts_'+msmt)
		plt.legend()

		print 'Plots are saved in:'
		print folder
		plt.savefig(os.path.join(folder,'ZenoCompareRuns'+RO_String+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'ZenoCompareRuns'+RO_String+'.png'),format='png')
		plt.show()
		plt.close('all')

	### if no plot then return the results.
	else: return evotime_arr[0], fid_arr, fid_u_arr,folder

def ShowResults():
	### this is a static function (no input parameters)
	### combines two sets (see older_timestamp_lists) of zeno measurements and plots them in a graph.
	### extendible for more than two measurement sets given that the evolution times stay the same.
	### saves a pickle file which contains the measured values as a dictionary.

	evo0,fid0,fid_u0,folder = Zen_Compare_Runs(msmt='0',older_timestamp_list=['20150404_175544','20150403_230000'],eRO_list=['positive','negative'],
	  						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)
	evo1,fid1,fid_u1,folder = Zen_Compare_Runs(msmt='1',older_timestamp_list=['20150404_175544','20150330_161704'],eRO_list=['positive','negative'],
	 						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)
	evo2,fid2,fid_u2,folder = Zen_Compare_Runs(msmt='2',older_timestamp_list=['20150404_175544','20150330_161704'],eRO_list=['positive','negative'],
	  						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)
	evo3,fid3,fid_u3,folder = Zen_Compare_Runs(msmt='3',older_timestamp_list=['20150404_175544','20150330_233000'],eRO_list=['positive','negative'],
	 						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)
	evo4,fid4,fid_u4,folder = Zen_Compare_Runs(msmt='4',older_timestamp_list=['20150404_175544','20150330_233000',],eRO_list=['positive','negative'],
	 						decoded_bit='2',ssro_timestamp_list=[None,None],plot_results=False)


	######################################
	### get the single qubit results. Beware, the function Zeno_proc_fidelity returns average state fidelities.

	evo_single, fid_single_avg, fid_u_single_avg, tstamp, folder = Zeno_proc_fidelity(older_than_tstamp='20150405_043659',
																		msmts='0',eRO_list=['positive','negative'],
																		ssro_timestamp=None,decoded_bit='2', plot_results = False, 
																		single_qubit = True)

	fid_single_proc = (3*fid_single_avg-1)/2.
	fid_u_single_proc = 1.5*fid_u_single_avg

	evo_single = evo_single*1000 ### rescale from seconds to ms.

	###########################


	fid_arr = [fid0,fid1,fid2,fid3,fid4]
	fid_u_arr = [fid_u0,fid_u1,fid_u2,fid_u3,fid_u4]
	evo_arr = [evo0,evo1,evo2,evo3,evo4]

	print fid_u_arr[0]

	### This loop fuses the measurement runs
	for i in range(len(fid_arr)):
		f_new = []
		f_u_new = []
		for j,fids in enumerate(fid_arr[i]):

			if j == 0:
				f_new = [fid/len(fid_arr[i]) for kk,fid in enumerate(fids)]
				f_u_new = [fid_u**2 for kk,fid_u in enumerate(fid_u_arr[i][j])]
			else:
				### add up fidelities
				f_new = [f_new[kk]+fid/len(fid_arr[i]) for kk,fid in enumerate(fids)]
				### start constructing the error bars. It is now a sum of squared errors.
				f_u_new = [f_u_new[kk]+fid_u**2 for kk,fid_u in enumerate(fid_u_arr[i][j])]
		### Finally, take the square root for the error entries and devide by the number of traces taken into account
		f_u_new = [np.sqrt(f_u_new[kk])/len(fid_arr[i]) for kk in range(len(f_u_new))]

		### overwrite the entries in fid_arr and fid_u_arr in order to make them ready for plotting.
		fid_arr[i] = f_new
		fid_u_arr[i] = f_u_new

	### start to plot the results
	fig=plt.figure()
	ax=plt.subplot()


	for ii in range(len(fid_arr)):

		plt.errorbar(evo_arr[ii],fid_arr[ii],fid_u_arr[ii],marker='o',markersize=4,label=str(ii))
	
	plt.xlabel('Evolution time (ms)')
	plt.ylabel('Process fidelity')
	plt.title('Process fidelity averaged over ' + str(2) + ' Zeno sets')
	plt.legend()

	print 'Plots are saved in:'
	print folder
	plt.savefig(os.path.join(folder,'ZenoRuns_averaged'+'.pdf'),format='pdf')
	plt.savefig(os.path.join(folder,'ZenoRuns_averaged'+'.png'),format='png')
	plt.show()
	plt.close('all')

	#### construct a pickle container out of the data
	### we start off with a dictionary
	pickle_dict = {'0':[evo_arr[0],fid_arr[0],fid_u_arr[0]],
					'1':[evo_arr[1],fid_arr[1],fid_u_arr[1]],
					'2':[evo_arr[2],fid_arr[2],fid_u_arr[2]],
					'3':[evo_arr[3],fid_arr[3],fid_u_arr[3]],
					'4':[evo_arr[4],fid_arr[4],fid_u_arr[4]],
					'single':[evo_single,fid_single_proc,fid_u_single_proc]}
					
	fileOut = open("Zeno_data.p","wb")
	pickle.dump(pickle_dict,fileOut)
	fileOut.close


def	Zeno_SingleQubit(older_than_tstamp=None,msmts='0',eRO_list=['positive'],
									 	state='X',
									 	plot_results=True,decoded_bit='2',
									 	ssro_timestamp=None):
		"""
		Plots or returns the state fidelity for a decoded qubit as a function of time (one parity expectation value)
		"""

		loop_bit = True
		evo_time=None

		y_arr=[]
		y_err_arr=[]
		evo_time_arr=[]
		ii=0

		evotime_dict = {
		'2' : 6.0,
		'4' : 14.0,
		'6' : 16.0,
		'8' : 16.0,
		}

		evo_time=[2005]
		while loop_bit:
			older_than_tstamp,loop_bit,x_labels,y,y_err,evo_time,folder,evo_was_zero=get_Zeno_data(electron_RO=eRO_list,
																						state=state,
																						older_than=older_than_tstamp,
																						previous_evo=evo_time[0],
																						msmts=msmts,
																						ssro_timestamp=ssro_timestamp,ROBasis='X')
			#loop_bit is true as long as new data was found.
			if loop_bit:
				y_arr=np.concatenate((y_arr,y))
				y_err_arr=np.concatenate((y_err_arr,y_err))
				evo_time_arr=np.concatenate((evo_time_arr,evo_time))
				# if an evolution time of zero was reached --> stop the evaluation.
				if evo_was_zero or evotime_dict[msmts] in evo_time:
					loop_bit=False



		#select the correct expectation value and the right sign for the contrast.
		sign=1

		if 'Y' in state:
			sign=-1
		elif state=='Z':
			sign=1
		elif state=='X':
			sign=1
		
		if 'm' in state:
			sign=-1*sign

		fid_arr=(sign*np.array(y_arr)+1)/2.
		fid_u_arr=np.array(y_err_arr)/2.

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		fid_arr=fid_arr[np.argsort(evo_time_arr)]
		fid_u_arr=fid_u_arr[np.argsort(evo_time_arr)]
		evo_time_arr=np.sort(evo_time_arr)


		if plot_results==True:
			fig=plt.figure()
			ax=plt.subplot()

			plt.errorbar(evo_time_arr,fid_arr,fid_u_arr,color='blue',marker='o')
			plt.xlabel('Evolution time (ms)')
			plt.ylabel('logical qubit state fidelity')
			plt.title('logicState_'+state+'_stop_'+str(older_than_tstamp)+'_'+RO_String)

			plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'_'+RO_String+'.pdf'),format='pdf')
			plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'_'+RO_String+'.png'),format='png')
			plt.show()
			plt.close('all')
			print folder
		else:
			return evo_time_arr,fid_arr,fid_u_arr,older_than_tstamp,folder

def Zeno_1Q_msmt_list(older_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive'],state='Z',ssro_timestamp=None,decoded_bit='2'):
	fid=[]
	fid_u=[]
	evotime=[]
	if len(eRO_list)==0:
		print 'nothing to do here'

	else:
		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder = Zeno_SingleQubit(older_than_tstamp=older_than_tstamp,
															msmts=msmt_list[i],
															eRO_list=eRO_list,decoded_bit=decoded_bit,
									plot_results=False,state=state,ssro_timestamp=ssro_timestamp)
			plt.errorbar(np.sort(evotime),fid[np.argsort(evotime)],fid_u[np.argsort(evotime)],marker='o',label=str(msmt_list[i]))
		

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		plt.xlabel('Evolution time (ms)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('logicState_'+state+'_stop_'+str(tstamp)+'_'+RO_String)
		plt.legend()

		print 'Plots are saved in:'
		print folder
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.png'),format='png')
		plt.show()
		plt.close('all')