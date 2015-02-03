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
				state='Z',ROBasis='IX',previous_evo=None, older_than=None,ssro_timestamp='20150105_155052'):
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
	"""

	search_string=electron_RO[0]+'_logicState_'+state+'_'+msmts+'msmt__RO'
	loop_bit = False

	timestamp=None

	if previous_evo == None:
		previous_evo = 2000 #very high value such that all evolution times are taken into account.

	#if the desired data set exists, then read the measured values.
	if toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=older_than,
									newer_than=None,
									raise_exc=False) != False:

		timestamp,folder=toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=older_than,
									newer_than=None,
									raise_exc=False)

		evotime,y,y_err= Zeno_get_2Q_values(timestamp)

		x_labels = folder[folder.find('ROBasis_')+8:]

	else: #if the searched data does not exist --> dummy results. evotime has to be larger than previous_evo
		x_labels,y,y_err,evotime= 'II',0,0,2001


	if evotime[0] < previous_evo or previous_evo==None:
		loop_bit = True
		
	else:
		x_labels,y,y_err,evotime= 'II',0,0,[2001]


	#if positive and negative RO are considered then adjust the search string and search for the same evo time.
	
	if len(electron_RO)>1 and (evotime < previous_evo or previous_evo==None):
		search_string=electron_RO[1]+'_logicState_'+state+'_'+msmts+'msmt__ROBasis'+str(ROBasis)

		timestamp2,folder2=toolbox.latest_data(contains=search_string,
									return_timestamp =True,
									older_than=newer_than,
									newer_than=timestamp, #choose the right search direction.
									raise_exc=False)

		if evotime[0] < previous_evo or previous_evo==None:
			loop_bit = True

			evotime2,y2,y_err2= Zeno_get_2Q_values(timestamp2,ssro_calib_timestamp =ssro_timestamp)
		else:
			x_labels,y2,y_err2,evotime2= 'II',0,0,2001

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


def Zeno_1Q_state_fidelity_decay(older_than_tstamp=None,msmts='0',eRO_list=['positive'],
								 	state='X',
								 	plot_results=True,
								 	ssro_timestamp='20150128_080328'):
	"""
	Plots or returns the state fidelity for a decoded qubit as a function of time (one parity expectation value)
	"""

	loop_bit = True
	evo_time=None
	x_labels_arr=np.array([])
	y_arr=[]
	y_err_arr=[]
	evo_time_arr=[]
	ii=0

	while loop_bit:
		older_than_tstamp,loop_bit,x_labels,y,y_err,evo_time,folder=get_Zeno_data(electron_RO=eRO_list,
																					state=state,
																					older_than=older_than_tstamp,
																					previous_evo=evo_time[-1],
																					msmts=msmts,
																					ssro_timestamp=ssro_timestamp)


		#loop_bit is true as long as new data was found.
		if loop_bit:
			x_labels_arr.extend(x_labels); y_arr.extend(y); y_err_arr.extend(y_err)
			evo_time_arr.extend(evo_time) 

	#select the correct expectation value and the right sign for the contrast.

	if state=='Y':
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


	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()

		plt.errorbar(evo_time_arr,fid_arr,fid_u_arr,color='blue',marker='o')
		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('logicState_'+state+'_stop_'+str(older_than_tstamp)+'_'+RO_String)

		plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'_'+RO_String+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'_'+RO_String+'.png'),format='png')
		plt.show()
		plt.close('all')
		print folder
	else:
		return evo_time_arr,fid_arr,fid_u_arr,older_than_tstamp,folder

def Zeno_1Q_proc_fidelity_decay(msmts='1',
								eRO_list=['positive'],
								older_than_tstamp=None,
								plot_results=True,
								ssro_timestamp='20150128_080328'):
	"""
	Plots the process fidelity for a decoded qubit as a function of time
	"""	
	state_list=['X','mX','Y','mY','Z','mZ']


	fid_arr=[]
	fid_u_arr=[]

	#get individual state fidelities
	for state in state_list:
		evo_time,fid,fid_u,tstamp,folder=Zeno_1Q_state_fidelity_decay(older_than_tstamp=older_than_tstamp,
										eRO_list=eRO_list,
										state=state, msmts=msmts,
										plot_results=False,
										ssro_timestamp=ssro_timestamp)
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
		ax=plt.subplot		

		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		plt.errorbar(evo_time,avg_fid,avg_fid_u,color='blue',marker='o')
		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit average fidelity')
		plt.title('Average_fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)

		print folder
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecay_decBit'+str(decoded_bit)+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecay_decBit'+str(decoded_bit)+'.png'),format='png')
		plt.show()
		plt.close('all')

		fig=plt.figure()
		ax=plt.subplot()

		for i,state in enumerate(state_list):
			plt.errorbar(evo_time,fid_arr[i],fid_u_arr[i],marker='o', label=state)

		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('State_fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)
		plt.legend()

		plt.savefig(os.path.join(folder,'Zeno1QStateFidelities_decBit'+str(decoded_bit)+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QStateFidelities_decBit'+str(decoded_bit)+'.png'),format='png')
		plt.show()
		plt.close('all')

	else:
		return evo_time,avg_fid,avg_fid_u,tstamp,folder

def Zeno_1Q_proc_list(older_than_tstamp=None,
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
			evotime,fid,fid_u,tstamp,folder = Zeno_1Q_proc_fidelity_decay(older_than_tstamp=older_than_tstamp,
									eRO_list=eRO_list, msmts=msmt_list[i],
									plot_results=False)
			plt.errorbar(evotime,fid,fid_u,marker='o',label=str(msmt_list[i])+' msmts')
		
		if len(eRO_list)==1:
			RO_String=eRO_list[0]
		else: RO_String = 'contrast'

		print evotime
		print fid
		print fid_u
		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit average fidelity')
		plt.title('Average_fidelity'+'_stop_'+str(tstamp)+'_'+RO_String)
		plt.legend()

		print folder
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QAvgDecays_decBit'+str(decoded_bit)+RO_String+'_combined.png'),format='png')
		plt.show()
		plt.close('all')

def Zeno_1Q_state_list(older_than_tstamp=None,
						msmt_list=['0'],eRO_list=['positive'],state='Z'):
	fid=[]
	fid_u=[]
	evotime=[]
	if len(eRO_list)==0:
		print 'nothing to do here'

	else:
		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(msmt_list)):
			evotime,fid,fid_u,tstamp,folder = Zeno_1Q_state_fidelity_decay(older_than_tstamp=older_than_tstamp,
															msmts=msmt_list[i],
															eRO_list=eRO_list,
									plot_results=False,state=state)
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