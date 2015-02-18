import numpy as np
import os,re
import h5py
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(toolbox)

"""
Analyze the zeno data
NK 2014
"""
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
			evotime[0]=first_dict_layer.attrs['parity_duration']
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
	a.get_electron_ROC(ssro_calib_folder)

	x_labels = a.sweep_pts.reshape(-1)
	y= ((a.p0.reshape(-1))-0.5)*2
	y_err = 2*a.u_p0.reshape(-1)

	return x_labels,y,y_err


def Zeno_state_fidelity(older_than_tstamp=None,msmts='0',eRO_list=['positive'],
								 	state='X',
								 	plot_results=True,decoded_bit='2',
								 	ssro_timestamp=None,single_qubit=False):
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
		plt.xlabel('Evolution time (s)')
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
		plt.xlabel('Evolution time (s)')
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

		plt.xlabel('Evolution time (s)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('State_fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)
		plt.legend()

		plt.savefig(os.path.join(folder,'Zeno1QStateFidelities_decBit'+str(decoded_bit)+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QStateFidelities_decBit'+str(decoded_bit)+'.png'),format='png')
		plt.show()
		plt.close('all')

	else:
		return evo_time,avg_fid,avg_fid_u,tstamp,folder

def Zeno_proc_list(older_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive'],decoded_bit='2',ssro_timestamp=None,single_qubit=False,plot_results=True):
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
		
		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		if plot_results:
			fig=plt.figure()
			ax=plt.subplot()
			for i in range(len(msmt_list)):
				plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],marker='o',markersize=4,label=str(msmt_list[i])+' msmts')
			
			if single_qubit: # adds the latest single qubit measurement to the data
				evotime,fid,fid_u,tstamp,folder = Zeno_proc_fidelity(older_than_tstamp=older_than_tstamp,
										eRO_list=eRO_list, msmts='0',ssro_timestamp=ssro_timestamp,decoded_bit=decoded_bit,
										plot_results=False,single_qubit=True)
				plt.errorbar(np.sort(evotime),(3*fid[np.argsort(evotime)]-1)/2,1.5*fid_u[np.argsort(evotime)],marker='o',markersize=4,label='1 qubit')

			plt.xlabel('Evolution time (s)')
			plt.ylabel('Process fidelity')
			plt.title('Process fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)
			plt.legend()

			print 'Plots are saved in:'
			print folder
			plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.pdf'),format='pdf')
			plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.png'),format='png')
			plt.show()
			plt.close('all')
		else: return evotime_arr,fid_arr,fid_u_arr,folder

def Zeno_state_list(older_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive'],state='Z',ssro_timestamp=None,decoded_bit='2',single_qubit=False):
	fid=[]
	fid_u=[]
	evotime=[]
	if len(eRO_list)==0:
		print 'nothing to do here'

	else:
		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder = Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
															msmts=msmt_list[i],
															eRO_list=eRO_list,decoded_bit=decoded_bit,
									plot_results=False,state=state,ssro_timestamp=ssro_timestamp)
			plt.errorbar(np.sort(evotime),fid[np.argsort(evotime)],fid_u[np.argsort(evotime)],marker='o',label=str(msmt_list[i])+' msmts')
		

		if single_qubit:
			evotime,fid,fid_u,tstamp,folder = Zeno_state_fidelity(older_than_tstamp=older_than_tstamp,
															msmts='0',
															eRO_list=eRO_list,decoded_bit=decoded_bit,
									plot_results=False,state=state,ssro_timestamp=ssro_timestamp,single_qubit=True)
			plt.errorbar(np.sort(evotime),fid[np.argsort(evotime)],fid_u[np.argsort(evotime)],marker='o',label='1 qubit')
		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		plt.xlabel('Evolution time (s)')
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


	fig=plt.figure()
	ax=plt.subplot()
	for i in range(len(older_timestamp_list)):
		plt.errorbar(evotime_arr[i],fid_arr[i],fid_u_arr[i],marker='o',markersize=4,label=str(older_timestamp_list[i]))
	
	plt.xlabel('Evolution time (s)')
	plt.ylabel('Process fidelity')
	plt.title('Process fidelity comparison'+'_'+RO_String+'_msmts_'+msmt)
	plt.legend()

	print 'Plots are saved in:'
	print folder
	plt.savefig(os.path.join(folder,'ZenoCompareRuns'+RO_String+'.pdf'),format='pdf')
	plt.savefig(os.path.join(folder,'ZenoCompareRuns'+RO_String+'.png'),format='png')
	plt.show()
	plt.close('all')

 