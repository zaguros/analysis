import numpy as np
import os,re
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

def Zeno_get_2Q_values(timestamp=None, folder=None,folder_name='Zeno',
						measurement_name = ['adwindata'], 
						ssro_calib_timestamp ='20150113_170834'):
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
	x = range(len(y)) 
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

def Zeno_1Q_state_fidelity_decay(timestamp_start=None,
									timestamp_stop=None,
									folder_name='Zeno',
								 	state='X',
								 	plot_results=True,
								 	decoded_bit=1):
	"""
	Plots or returns the state fidelity for a decoded qubit as a function of time (one parity expectation value)
	"""

	if timestamp_start==None:
		timestamp_start, folder   = toolbox.latest_data(folder_name,return_timestamp =True)

	if timestamp_stop==None:
		return Zeno_get_2Q_values(timestamp_start)
		#todo, convert to single q value

	loop_bit = True
	x_labels_arr=[]
	y_arr=[]
	x_arr=[]
	y_err_arr=[]
	evo_time_arr=[]
	ii=0

	while loop_bit:
		if toolbox.latest_data(contains='positive_logicState_'+state,
													return_timestamp =True,
													older_than=timestamp_start,
													newer_than=timestamp_stop,
													raise_exc=False) == False:
			loop_bit =False
		else:
			timestamp_start,folder=toolbox.latest_data(contains='positive_logicState_'+state,
													return_timestamp =True,
													older_than=timestamp_start,
													newer_than=timestamp_stop,
													raise_exc=False)
			x_labels,y,y_err= Zeno_get_2Q_values(timestamp_start)
			x_labels_arr.append(x_labels); y_arr.append(y); y_err_arr.append(y_err)
			evo_time_arr.append(float(folder[folder.find('EvoTime_')+8:])) #extract free evolution time

	#select the correct expectation value.
	index=0
	sign=1

	if decoded_bit == 1:
		if state=='mY' or state=='Y':
			for k,x in enumerate(x_labels):
				if x=='YZ':
					index=k
					sign=-1
		elif state=='mZ' or state=='Z':
			for k,x in enumerate(x_labels):
				if x=='XI':
					index=k
		elif state=='mX' or state=='X':
			for k,x in enumerate(x_labels):
				if x=='ZZ':
					index=k

	if decoded_bit ==2:
		if state=='mY' or state=='Y':
			for k,x in enumerate(x_labels):
				if x=='ZY':
					index=k
					sign=-1
		elif state=='mZ' or state=='Z':
			for k,x in enumerate(x_labels):
				if x=='IX':
					index=k
		elif state=='mX' or state=='X':
			for k,x in enumerate(x_labels):
				if x=='ZZ':
					index=k
	
	if 'm' in state:
		sign=-1*sign

	fid_arr=np.zeros(len(y_arr))
	fid_u_arr=np.zeros(len(y_arr))
	for i,y in enumerate(y_arr):
		fid_arr[i]=(sign*y[index]+1.)/2.
		fid_u_arr[i]=(y_err_arr[i][index])/2.

	if plot_results==True:
		fig=plt.figure()
		ax=plt.subplot()

		plt.errorbar(evo_time_arr,fid_arr,fid_u_arr,color='blue',marker='o')
		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('logicState_'+state+'_start_'+str(timestamp_start)+'_stop_'+str(timestamp_stop))

		plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'.pdf'),format='pdf')
		plt.savefig(os.path.join(folder,'Zeno1QDecay_decBit'+str(decoded_bit)+'.png'),format='png')
		plt.show()
		plt.close('all')
		print folder
	else:
		return evo_time_arr,fid_arr,fid_u_arr

def Zeno_1Q_proc_fidelity_decay(timestamp_start=None,timestamp_stop=None,plot_results=True,folder_name='Zeno',decoded_bit=1):
	"""
	Plots the process fidelity for a decoded qubit as a function of time
	"""	
	state_list=['X','mX','Y','mY','Z','mZ']

	evo_time_arr=[]
	fid_arr=[]
	fid_u_arr=[]

	#get individual state fidelities
	for state in state_list:
		evo_time,fid,fid_u=Zeno_1Q_state_fidelity_decay( timestamp_start=timestamp_start,
										timestamp_stop=timestamp_stop,
										state=state,plot_results=False,decoded_bit=decoded_bit)
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
		return evo_time,avg_fid,avg_fid_u

def Zeno_1Q_proc_list(starts=[],
						stops=[],decoded_bit=1):
	fid=[]
	fid_u=[]
	evotime=[]
	if len(starts)==0:
		print 'nothing to do here'

	else:
		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(starts)):
			evotime,fid,fid_u = Zeno_1Q_proc_fidelity_decay(timestamp_start=starts[i],
									timestamp_stop=stops[i],
									plot_results=False,decoded_bit=decoded_bit)
			plt.errorbar(evotime,fid,fid_u,marker='o')
		

		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit average fidelity')
		plt.title('Timestamps_start_'+str(starts[0])+'_stop_'+str(stops[-1]))

		print toolbox.latest_data('Zeno',older_than=starts[-1])
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=starts[-1]),'Zeno1QAvgDecays_decBit'+str(decoded_bit)+'_combined.pdf'),format='pdf')
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=starts[-1]),'Zeno1QAvgDecays_decBit'+str(decoded_bit)+'_combined.png'),format='png')
		plt.show()
		plt.close('all')

def Zeno_1Q_state_list(starts=[],
						stops=[],decoded_bit=1,state='Z'):
	fid=[]
	fid_u=[]
	evotime=[]
	if len(starts)==0:
		print 'nothing to do here'

	else:
		fig=plt.figure()
		ax=plt.subplot()
		for i in range(len(starts)):
			evotime,fid,fid_u = Zeno_1Q_state_fidelity_decay(timestamp_start=starts[i],
									timestamp_stop=stops[i],
									plot_results=False,decoded_bit=decoded_bit,state=state)
			plt.errorbar(evotime,fid,fid_u,marker='o')
		

		plt.xlabel('Free evolution time (s)')
		plt.ylabel('logical qubit state fidelity')
		plt.title('Timestamps_start_'+str(starts[0])+'_stop_'+str(stops[-1]))

		print toolbox.latest_data('Zeno',older_than=starts[-1])
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=starts[-1]),'Zeno1QAvgDecays_decBit'+str(decoded_bit)+'_combined.pdf'),format='pdf')
		plt.savefig(os.path.join(toolbox.latest_data('Zeno',older_than=starts[-1]),'Zeno1QAvgDecays_decBit'+str(decoded_bit)+'_combined.png'),format='png')
		plt.show()
		plt.close('all')