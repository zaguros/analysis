from analysis.lib.nv import nvlevels
reload(nvlevels)
import numpy as np

def plot_ES_b_dependence(strain_splitting=2.):
	b_range=np.linspace(0,1000,100) #gauss
	Ex=strain_splitting/2. #Ex is approx strain_splitting divided by 2
	spectrum=np.zeros((6,))
	for Bz in b_range:
		spectrum=np.vstack((spectrum,np.sort(nvlevels.get_ES(
															E_field=[Ex,0.,0.], 
															B_field=[0.,0.,Bz],
															Ee0=-1.94,
															transitions=False,
															)[0])))
	spectrum=spectrum[1:]

	for i in range(6):
		plot(b_range,spectrum[:,i])

def plot_ES_e_dependence(Bz=20.):
	e_range=np.linspace(0,10,100) #gauss
	spectrum=np.zeros((6,))
	for Ex in e_range:
		spectrum=np.vstack((spectrum,np.sort(nvlevels.get_ES(
															E_field=[Ex,0.,0.], 
															B_field=[0.,0.,Bz],
															Ee0=-1.94,
															transitions=False,
															)[0])))
	spectrum=spectrum[1:]

	for i in range(6):
		plot(e_range,spectrum[:,i])

def plot_transitions_b_dependence(strain_splitting=2.5):
	b_range=np.linspace(0,1000,100) #gauss
	Ex=strain_splitting/2. #Ex is approx strain_splitting divided by 2
	no_transitions=8
	spectrum=np.zeros((no_transitions,))
	for Bz in b_range:
		spectrum=np.vstack((spectrum,nvlevels.get_optical_transitions(
															E_field=[Ex,0.,0.], 
															B_field=[0,0.,Bz],
															Ee0=-1.94,
															show_FB_E_transitions=False
															)))
	spectrum=spectrum[1:]

	for i in range(no_transitions):
		plot(b_range,spectrum[:,i])


if __name__=='__main__':
	plot_ES_b_dependence()
	figure()
	plot_transitions_b_dependence()
