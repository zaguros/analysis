import os
import numpy as np
from matplotlib import pyplot as plt
#from measurement.lib.tools import toolbox


def sweep_pulse_ampl(pulse_duration=100e-9, Rabi_freq_max = 50e6, pts = 1000, detuning = 0):
	
	title= 'Simul Rabi with pulse of ' + str(pulse_duration) +' s'

	hyperfine = 2.17e6

	# Rabi frequency
	Rabi_freq=np.linspace(0.00001,Rabi_freq_max,pts)
	Proba_0=np.ones(pts)


	for i in np.linspace(-1,1,3):
		Proba_0 = Proba_0-1/3.*Rabi_freq**2/(Rabi_freq**2+ (detuning+i*hyperfine)**2)*np.sin(np.pi*np.sqrt(Rabi_freq**2+(detuning+i*hyperfine)**2)*pulse_duration)**2



	Proba_0_no_detun = 1-np.sin(np.pi*Rabi_freq*pulse_duration)**2

	Rabi_freq=Rabi_freq/1e6

	fig = plt.figure()
	ax = fig.add_subplot(111)


	ax.plot(Rabi_freq, Proba_0_no_detun,'k--',label='no $^{14}$N' )
	ax.plot(Rabi_freq, Proba_0,'r', label= '$^{14}$N' )   #, 'bo', ms = 8, label = '$Rabi$')


	ax.set_xlabel('Rabi frequency [MHz]', fontsize = 16)
	ax.set_ylabel('Proba in ms=0', fontsize = 16)
	ax.set_title(title, fontsize= 16)

	ax.tick_params(labelsize=14)

	ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)



def sweep_pulse_duration(pulse_duration_max=2000e-9, Rabi_freq = 2e6, pts = 1000, detuning = 0):
	
	title= 'Simul Rabi with Rabi freq. of ' + str(Rabi_freq/1e6) +' MHz \n detuning = ' + str(detuning/1e6) + ' MHz' 


	hyperfine = 2.17e6

	# Rabi frequency
	pulse_duration=np.linspace(1e-12,pulse_duration_max,pts)
	Proba_0=np.ones(pts)


	for i in np.linspace(-1,1,3):
		Proba_0 = Proba_0-1/3.*Rabi_freq**2/(Rabi_freq**2+ (detuning+i*hyperfine)**2)*np.sin(np.pi*np.sqrt(Rabi_freq**2+(detuning+i*hyperfine)**2)*pulse_duration)**2



	Proba_0_no_detun = 1-np.sin(np.pi*Rabi_freq*pulse_duration)**2

	Rabi_freq=Rabi_freq/1e6

	fig = plt.figure()
	ax = fig.add_subplot(111)

	pulse_duration_ns = pulse_duration*1e9

	ax.plot(pulse_duration*1e9, Proba_0_no_detun,'k--',label='no $^{14}$N' )
	ax.plot(pulse_duration*1e9, Proba_0,'r', label= '$^{14}$N' )   #, 'bo', ms = 8, label = '$Rabi$')



	ax.set_xlabel('Pulse duration [ns]', fontsize = 16)
	ax.set_ylabel('Proba in ms=0', fontsize = 16)
	ax.set_title(title, fontsize= 16)

	ax.tick_params(labelsize=14)

	ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
	#ax.set_ylim([0.,.05])
	ax.set_ylim([0.6,0.8])