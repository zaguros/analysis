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
				state='Z', newer_than=None, previous_evo=None):
	"""
	this function finds data according to an input which specifies the folder name
	Input: 	electron_RO 		is a list with e.g. positive and negative. 
								If the list is longer than 1 element then the 
								contrast values are returned.
			previous_evo		As soon as this reaches 0 the analysis is stopped. 
								This assumes that the evolution time is sweeped from low to high
								in actual measurements.
			msmts & state 		both strings are used to find the correct folder

	Output: timestamp			specifies the latest evaluated folder
			loop_bit			Boolean which signals the end of data acquisition and if the output data should be evaluated.
			x_labels 			Tomographic bases. E.g. XI or ZZ
			y 					Read-out results
			y_err 				Read-out uncertainty
			evo_time 			Free evolution of the specific state
			folder 				The folder of the latest timestamp
	"""

	search_string=electron_RO[0]+'_logicState_'+state+'_'+msmts+'msmt_'

	loop_bit = False

	if previous_evo == None:
		previous_evo = 2000 #very high value such that all evolution times are taken into account.

	#if the desired data set exists, then read the measured values.
	if toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=newer_than,
									newer_than=None,
									raise_exc=False) != False:

		timestamp,folder=toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=newer_than,
									newer_than=None,
									raise_exc=False)

		evotime = float(folder[folder.find('EvoTime_')+8:])
	else: #if it does not exits --> dummy results. evotime has to be larger than previous_evo
		x_labels,y,y_err,evotime= 0,0,0,2001


	if evotime < previous_evo or previous_evo==None:
		loop_bit = True
		x_labels,y,y_err= Zeno_get_2Q_values(timestamp)
	else:
		x_labels,y,y_err,evotime= 0,0,0,2001


	datafile=h5py.File(os.path.join(folder,folder[27:]+'.hdf5'),'r') # open the data read-only
	#datafile.attrs['']
	datafile.close()

	#if positive and negative RO are considered then adjust the search string and search for the same evo time.
	
	if len(electron_RO)>1 and (evotime < previous_evo or previous_evo==None):
		search_string=electron_RO[1]+'_logicState_'+state+'_'+msmts+'msmt__EvoTime_'+str(evotime)

		timestamp2,folder2=toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=newer_than,
									newer_than=None,
									raise_exc=False)

		if evotime < previous_evo or previous_evo==None:
			loop_bit = True

			x_labels,y2,y_err2= Zeno_get_2Q_values(timestamp2)
		else:
			x_labels,y2,y_err2,evotime= 0,0,0,2001

		if electron_RO[0]== 'positive':	
			for i in range(len(y)):
				y[i]=(y[i]-y2[i])/2
				y_err[i]=((y_err[i]**2+y_err2[i]**2)**0.5)/2

		if electron_RO[0] == 'negative':
			for i in range(len(y)):
				y[i]=(-y[i]+y2[i])/2
				y_err[i]=((y_err[i]**2+y_err2[i]**2)**0.5)/2

	#determine the older timestamp (for the two eRO possiblities) and return that one.
	
	if loop_bit:
		if len(electron_RO)==1:
			return timestamp,loop_bit,x_labels,y,y_err,evotime,folder

		elif toolbox.is_older(timestamp,timestamp2):
			return timestamp,loop_bit,x_labels,y,y_err,evotime,folder

		else: 
			return timestamp2,loop_bit,x_labels,y,y_err,evotime,folder2
	else:
		return newer_than,loop_bit,x_labels,y,y_err,evotime,toolbox.data_from_time(newer_than)


def Zeno_get_2Q_values(timestamp=None, folder=None,folder_name='Zeno',
						measurement_name = ['adwindata'], 
						ssro_calib_timestamp ='20150128_080328'):
	"""
	Returns the relevant 2qubit values for a given timestamp.
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

		
def Zeno_2Q_state_fidelity_decay(timestamp_start=None,
									timestamp_stop=None,
									folder_name='Zeno',
								 	state='X',
								 	plot_results=True):
	"""
	Plots or returns the state fidelity for a 2-qubit state as a function of time
	"""


	if timestamp_start==None:
		timestamp_start, folder   = toolbox.latest_data(folder_name,return_timestamp =True)

	if timestamp_stop==None:
		return Zeno_get_2Q_values(timestamp_start)

	loop_bit = True
	x_labels_arr=[]
	y_arr=[]
	x_arr=[]
	y_err_arr=[]
	evo_time_arr=[]
	ii=0

	while loop_bit:
		if toolbox.latest_data(contains='negative_logicState_'+state,
													return_timestamp =True,
													older_than=timestamp_start,
													newer_than=timestamp_stop,
													raise_exc=False) == False:
			print 'false reached'
			loop_bit =False
		else:
			timestamp_start,folder=toolbox.latest_data(contains='negative_logicState_'+state,
													return_timestamp =True,
													older_than=timestamp_start,
													newer_than=timestamp_stop,
													raise_exc=False)
			x_labels,y,y_err= Zeno_get_2Q_values(timestamp_start)
			x_labels_arr.append(x_labels); y_arr.append(y); y_err_arr.append(y_err)
			evo_time_arr.append(float(folder[folder.find('EvoTime_')+8:])) #extract free evolution time


	state_list=['X','mX','Y','mY','Z','mZ']
	sign_list=[[1,-1,1],[1,1,-1],[1,-1,-1],[1,1,1],[1,1,1],[-1,-1,1]]
	sign1,sign2,sign3=1,1,1
	#find the correct signs for the fidelity calculation of the 2-qubit fidelity
	for i,s in enumerate(state_list):
		if state==s:
			sign1=sign_list[i][0]
			sign2=sign_list[i][1]
			sign3=sign_list[i][2]

	fid_arr=np.zeros(len(y_arr))
	fid_u_arr=np.zeros(len(y_arr))
	for i,y in enumerate(y_arr):
		fid_arr[i]=(sign1*y[0]+sign2*y[1]+sign3*y[2]+1.)/4.
		fid_u_arr[i]=np.sqrt(y_err_arr[i][0]**2+y_err_arr[i][1]**2+y_err_arr[i][2]**2)/4.

	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()

		plt.errorbar(evo_time_arr,fid_arr,fid_u_arr,color='blue',marker='o')
		plt.xlabel('Free evolution time (s)')
		plt.ylabel('2 qubit state fidelity')
		plt.title('logicState_'+state+'_start_'+str(timestamp_start)+'_stop_'+str(timestamp_stop))

		plt.savefig(os.path.join(folder,'ZenoTwoQDecay.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'ZenoTwoQDecay.png'),format='png')
		plt.show()
		plt.close('all')
		print folder
	else:
		return evo_time_arr,fid_arr,fid_u_arr

def Zeno_2Q_proc_fidelity_decay(timestamp_start=None,
									timestamp_stop=None,
									folder_name='Zeno',plot_results=True):
	"""
	Plots the process fidelity for a 2-qubit state as a function of time
	"""
	state_list=['X','mX','Y','mY','Z','mZ']

	evo_time_arr=[]
	fid_arr=[]
	fid_u_arr=[]

	#get individual state fidelities
	for state in state_list:
		evo_time,fid,fid_u=Zeno_2Q_state_fidelity_decay( timestamp_start=timestamp_start,
										timestamp_stop=timestamp_stop,
										state=state,plot_results=False)
		fid_arr.append(fid);fid_u_arr.append(fid_u)

	#calculate average state fidelity
	avg_fid=np.zeros(len(fid_arr[0]))
	avg_fid_u=np.zeros(len(fid_u_arr[0]))
	for i in range(len(fid_arr[0])):
		for ii in range(len(fid_arr)):
			avg_fid[i]= avg_fid[i] + fid_arr[ii][i]/len(fid_arr)
			avg_fid_u[i]= avg_fid_u[i] + fid_u_arr[ii][i]**2/36
		avg_fid_u[i]=avg_fid_u[i]**0.5


	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()

		plt.errorbar(evo_time,avg_fid,avg_fid_u,color='blue',marker='o')
		plt.xlabel('Free evolution time (s)')
		plt.ylabel('2 qubit average fidelity')
		plt.title('Timestamps_start_'+str(timestamp_start)+'_stop_'+str(timestamp_stop))

		print toolbox.latest_data('Zeno',older_than=timestamp_start)
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=timestamp_start),'ZenoTwoQAvgDecay.pdf'),format='pdf')
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=timestamp_start),'ZenoTwoQAvgDecay.png'),format='png')

		fig=plt.figure()
		ax=plt.subplot()

		for i,state in enumerate(state_list):
			plt.errorbar(evo_time,fid_arr[i],fid_u_arr[i],marker='o', label=state)

		plt.xlabel('Free evolution time (s)')
		plt.ylabel('2 qubit state fidelity')
		plt.title('Timestamps_start_'+str(timestamp_start)+'_stop_'+str(timestamp_stop))
		plt.legend()

		print toolbox.latest_data('Zeno',older_than=timestamp_start)
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=timestamp_start),'ZenoTwoQStateFidelities.pdf'),format='pdf')
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=timestamp_start),'ZenoTwoQStateFidelities.png'),format='png')
		plt.close('all')

	else:
		return evo_time,avg_fid,avg_fid_u


def Zeno_2Q_proc_list(starts=[],
						stops=[]):
	fid=[]
	fid_u=[]
	evotime=[]
	if len(starts)==0:
		print 'nothing to do here'

	else:
		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(starts)):
			evotime,fid,fid_u = Zeno_2Q_proc_fidelity_decay(timestamp_start=starts[i],
									timestamp_stop=stops[i],
									plot_results=False)
			plt.errorbar(evotime,fid,fid_u,marker='o')
		

		plt.xlabel('Free evolution time (s)')
		plt.ylabel('2 qubit average fidelity')
		plt.title('Timestamps_start_'+str(starts[0])+'_stop_'+str(stops[-1]))

		print toolbox.latest_data('Zeno',older_than=starts[-1])
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=starts[-1]),'ZenoTwoQAvgDecays_combined.pdf'),format='pdf')
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=starts[-1]),'ZenoTwoQAvgDecays_combined.png'),format='png')
		plt.close('all')

		

def Zeno_1Q_state_fidelity_decay(newer_than_tstamp=None,msmts='0',eRO_list=['positive'],
								 	state='X',
								 	plot_results=True,
								 	decoded_bit=1):
	"""
	Plots or returns the state fidelity for a decoded qubit as a function of time (one parity expectation value)
	"""

	loop_bit = True
	evo_time=None
	x_labels_arr=[]
	y_arr=[]
	y_err_arr=[]
	evo_time_arr=[]
	ii=0

	while loop_bit:
		newer_than_tstamp,loop_bit,x_labels,y,y_err,evo_time,folder=get_Zeno_data(electron_RO=eRO_list,state=state,newer_than=newer_than_tstamp,previous_evo=evo_time,msmts=msmts)


		#loop_bit is true as long as new data was found.
		if loop_bit:
			x_labels_arr.append(x_labels); y_arr.append(y); y_err_arr.append(y_err)
			evo_time_arr.append(evo_time) #extract free evolution time

	#select the correct expectation value.
	index=0
	sign=1

	if decoded_bit == 1:
		if state=='mY' or state=='Y':
			for k,x in enumerate(x_labels_arr[-1]):
				if x=='YZ':
					index=k
					sign=-1
		elif state=='mZ' or state=='Z':
			for k,x in enumerate(x_labels_arr[-1]):
				if x=='XI':
					index=k
		elif state=='mX' or state=='X':
			for k,x in enumerate(x_labels_arr[-1]):
				if x=='ZZ':
					index=k

	if decoded_bit ==2:
		if state=='mY' or state=='Y':
			for k,x in enumerate(x_labels_arr[-1]):
				if x=='ZY':
					index=k
					sign=-1
		elif state=='mZ' or state=='Z':
			for k,x in enumerate(x_labels_arr[-1]):
				if x=='IX':
					index=k
		elif state=='mX' or state=='X':
			for k,x in enumerate(x_labels_arr[-1]):
				if x=='ZZ':
					index=k
	
	if 'm' in state:
		sign=-1*sign

	fid_arr=np.zeros(len(y_arr))
	fid_u_arr=np.zeros(len(y_arr))
	for i,y in enumerate(y_arr):
		fid_arr[i]=(sign*y[index]+1.)/2.
		fid_u_arr[i]=(y_err_arr[i][index])/2.

	if len(eRO_list)==1:
		RO_String=eRO_list[0]
	else: RO_String = 'contrast'


	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()

		plt.errorbar(evo_time_arr,fid_arr,fid_u_arr,color='blue',marker='o')
		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('logicState_'+state+'_stop_'+str(newer_than_tstamp)+'_'+RO_String)

		plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'_'+RO_String+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'_'+RO_String+'.png'),format='png')
		plt.show()
		plt.close('all')
		print folder
	else:
		return evo_time_arr,fid_arr,fid_u_arr,newer_than_tstamp,folder

def Zeno_1Q_proc_fidelity_decay(msmts='1',
								eRO_list=['positive'],
								newer_than_tstamp=None,
								plot_results=True,
								decoded_bit=1):
	"""
	Plots the process fidelity for a decoded qubit as a function of time
	"""	
	state_list=['X','mX','Y','mY','Z','mZ']


	fid_arr=[]
	fid_u_arr=[]

	#get individual state fidelities
	for state in state_list:
		evo_time,fid,fid_u,tstamp,folder=Zeno_1Q_state_fidelity_decay(newer_than_tstamp=newer_than_tstamp,
										eRO_list=eRO_list,
										state=state, msmts=msmts,
										plot_results=False,
										decoded_bit=decoded_bit)
		fid_arr.append(fid);fid_u_arr.append(fid_u)

	#calculate average state fidelity
	avg_fid=np.zeros(len(fid_arr[0]))
	avg_fid_u=np.zeros(len(fid_u_arr[0]))
	for i in range(len(fid_arr[0])):
		for ii in range(len(fid_arr)):
			avg_fid[i]= avg_fid[i] + fid_arr[ii][i]/len(fid_arr)
			avg_fid_u[i]= avg_fid_u[i] + fid_u_arr[ii][i]**2/36
		avg_fid_u[i]=avg_fid_u[i]**0.5


	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()

		plt.errorbar(evo_time,avg_fid,avg_fid_u,color='blue',marker='o')
		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit average fidelity')
		plt.title('Timestamps_start_'+str(timestamp_start)+'_stop_'+str(timestamp_stop))

		print toolbox.latest_data('Zeno',older_than=timestamp_start)
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=timestamp_start),'Zeno1QAvgDecay_decBit'+str(decoded_bit)+'.pdf'),format='pdf')
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=timestamp_start),'Zeno1QAvgDecay_decBit'+str(decoded_bit)+'.png'),format='png')

		fig=plt.figure()
		ax=plt.subplot()

		for i,state in enumerate(state_list):
			plt.errorbar(evo_time,fid_arr[i],fid_u_arr[i],marker='o', label=state)

		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('Timestamps_start_'+str(timestamp_start)+'_stop_'+str(timestamp_stop))
		plt.legend()

		print toolbox.latest_data('Zeno',older_than=timestamp_start)
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=timestamp_start),'Zeno1QStateFidelities_decBit'+str(decoded_bit)+'.pdf'),format='pdf')
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=timestamp_start),'Zeno1QStateFidelities_decBit'+str(decoded_bit)+'.png'),format='png')
		plt.close('all')

	else:
		return evo_time,avg_fid,avg_fid_u,tstamp,folder

def Zeno_1Q_proc_list(newer_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive'],decoded_bit=1):
	fid=[]
	fid_u=[]
	evotime=[]
	if len(msmt_list)==0:
		print 'nothing to do here'

	else:
		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder = Zeno_1Q_proc_fidelity_decay(newer_than_tstamp=newer_than_tstamp,
									eRO_list=eRO_list, msmts=msmt_list[i],
									plot_results=False,decoded_bit=decoded_bit)
			plt.errorbar(evotime,fid,fid_u,marker='o',label=str(msmt_list[i])+' msmts')
		
		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit average fidelity')
		plt.title('Average_fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)

		print folder
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.png'),format='png')
		plt.show()
		plt.close('all')

def Zeno_1Q_state_list(newer_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive'],decoded_bit=1,state='Z'):
	fid=[]
	fid_u=[]
	evotime=[]
	if len(eRO_list)==0:
		print 'nothing to do here'

	else:
		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder = Zeno_1Q_state_fidelity_decay(newer_than_tstamp=newer_than_tstamp,
															msmts=msmt_list[i],
															eRO_list=eRO_list,
									plot_results=False,decoded_bit=decoded_bit,state=state)
			plt.errorbar(evotime,fid,fid_u,marker='o',label=str(msmt_list[i])+' msmts')
		

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('logicState_'+state+'_stop_'+str(tstamp)+'_'+RO_String)
		plt.legend()

		print folder
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.png'),format='png')
		plt.show()
		plt.close('all')