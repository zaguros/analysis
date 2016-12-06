''' A module to calculate the 13C nuclear and electron spin dynamics
under dynamical decoupling gates. By THT '''

import numpy as np
import qutip
import analysis.lib.QEC.hyperfine_params as hf
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import warnings

### import the hyperfine parameters ###
import hyperfine_params as hf_params; reload(hf_params)
hf = hf_params.hyperfine_params
# hf = hf_params.SamSIL5_sim_params


### import the experimental values for tau and N ###
# import measurement.scripts.lt2_scripts.setup.msmt_params as msmt_params; reload(msmt_params)

### import the theoretically tuned values for tau and N ###
import gate_params as gate_params; reload(gate_params)
mp = gate_params.gp

#######################
### Basic functions ###
#######################

def pauli():
	'''Define pauli spin matrices'''
	identity = qutip.qeye(2)
	sx = qutip.sigmax()/2
	sy = qutip.sigmay()/2
	sz = qutip.sigmaz()/2
	return identity, sx, sy, sz

def basic_spin_rotations():
	''' define some simple spin rotations'''
	X = (-1j*sx*np.pi).expm();   mX = (1j*sx*np.pi).expm()
	Y = (-1j*sy*np.pi).expm();   mY = (1j*sy*np.pi).expm()
	Z = (-1j*sz*np.pi).expm();   mZ = (1j*sz*np.pi).expm()
	x = (-1j*sx*np.pi/2).expm(); mx = (1j*sx*np.pi/2).expm()
	y = (-1j*sy*np.pi/2).expm(); my = (1j*sy*np.pi/2).expm()
	z = (-1j*sz*np.pi/2).expm(); mz = (1j*sz*np.pi/2).expm()
	H = (-1j*(sx+sz)/np.sqrt(2)*np.pi).expm()

	return X,Y,Z,x,y,z,mX,mY,mZ,mx,my,mz

def basic_spin_states():
	'''define some basic spin states'''
	ket0 = qutip.basis(2,0)
	bra0 = qutip.basis(2,0).dag()
	ket1 = qutip.basis(2,1)
	bra1 = qutip.basis(2,1).dag()
	rho0 = ket0*bra0
	rho1 = ket1*bra1
	rhom = qutip.qeye(2)/2
	ketx = 1/np.sqrt(2)*(qutip.basis(2,0)+qutip.basis(2,1))
	brax = 1/np.sqrt(2)*(qutip.basis(2,0).dag()+qutip.basis(2,1).dag())
	ketmx = 1/np.sqrt(2)*(qutip.basis(2,0)-qutip.basis(2,1))
	bramx = 1/np.sqrt(2)*(qutip.basis(2,0).dag()-qutip.basis(2,1).dag())
	kety = 1/np.sqrt(2)*(qutip.basis(2,0)+1j*qutip.basis(2,1))
	bray = kety.dag()
	ketmy = 1/np.sqrt(2)*(qutip.basis(2,0)-1j*qutip.basis(2,1))
	bramy = ketmy.dag()
	rhox =ketx*brax
	rhomx = ketmx*bramx
	rhoy =kety*bray
	rhomy = ketmy*bramy
	return ket0, bra0, ket1, bra1, rho0, rho1, rhom, ketx,brax,ketmx,bramx,rhox,rhomx,kety,bray,ketmy,bramy,rhoy,rhomy

### create a set of usefull simple states and gates
Id, sx, sy, sz = pauli()                                # Electron spin operators
Id, Ix, Iy, Iz = pauli()                                # Nuclear spin operators
X,Y,Z,x,y,z,mX,mY,mZ,mx,my,mz = basic_spin_rotations()  # Basic gates
ket0, bra0, ket1, bra1, rho0, rho1, rhom, ketx,brax,ketmx,bramx,rhox,rhomx,kety,bray,ketmy,bramy,rhoy,rhomy = basic_spin_states() # Basic states

###########################
###  Helper functions   ###
###########################

def any_pure_state(alpha,beta,return_psi = False,return_rho = True):
	'''gives out your psi and if wanted your rho for a state alpha 0 + beta 1 '''
	psi = alpha*ket0+beta*ket1
	rho = psi*psi.dag()
	# print psi
	# print rho
	if return_psi == True:
		return psi
	if return_rho == True:
		return rho

def any_mixed_state(alpha,beta):
	'''gives out a mixture of rho0 and rho1'''
	rho = alpha *rho0+beta*rho1
	return rho

def evolve(system,operator):
	return operator * system * operator.dag()

def print_matrix(Qobject,div_by=100):

	print np.round(Qobject.full()*div_by)/div_by
	print type(np.round(Qobject.full()*div_by)/div_by)


class NV_system(object):
	''' Basic class to contain the parameters of the NV system, plus some functions to simulate its evolution '''

	def __init__(self,**kw):

		if kw.pop('use_msmt_params',False):
			raise Exception('Not written yet!')
		else:
			self.carbon_gate_params = kw.pop('carbon_gate_params',[])
			if kw.pop('carbon_params',False):
				self.num_carbons = len(carbon_params) 
				self.carbon_params = carbon_params
				
			else:
				warnings.warn('No carbon params passed, using dummy params!')
				self.num_carbons = 1
				self.carbon_params = [[2 * np.pi * 443141.64, 2 * np.pi * 35e3 , 2 * np.pi * 35e3]]

		self.decouple_scheme = 'simple'

		self.B_field = kw.pop('B_field',400)

		self.initial_state = self.calc_initial_state()

		self.espin_trans = '+1' # Changes sign of A_par and A_perp

	def calc_initial_state(self, NV_state = rho0):
		''' Helper function to give standard init state for system '''
		return qutip.tensor([NV_state] + [rhom] * self.num_carbons)

	def embed_e_operator(self,operator):
		''' Helper function to embed e operator with identities on carbons '''
		return qutip.tensor([operator] + [Id] * self.num_carbons)


	def NV_carbon_evolution_matrix(self,tau):
		''' Function to calculate a C13 evolution matrix, carbon params = [[omega_Larmor, A_par, A_perp]]'''

		sign = -1 if self.espin_trans == '-1' else 1

		H0s = [0.5*carbon_param[0]*Iz for carbon_param in self.carbon_params] 
		H1s = [0.5*((carbon_param[0]+ sign*carbon_param[1])*Iz + sign * carbon_param[2] * Ix) for carbon_param in self.carbon_params] 
		
		#Hamiltonians for ms=0 and ms=+/-1
		H0 = qutip.tensor(H0s)
		H1 = qutip.tensor(H1s)
	 
		#Evolution during tau for ms=0 and ms=+/-1
		expH = (-1j*(qutip.tensor(rho0,H0)+qutip.tensor(rho1,H1))*tau).expm(); 
	   
		return expH

	def decouple_unit_evolution(self,tau):
		'''Evolution during a decouple unit'''
		Xe = self.embed_e_operator(X)
		Ye = self.embed_e_operator(Y)
		
		expH = self.NV_carbon_evolution_matrix(tau)
		expHhalf = self.NV_carbon_evolution_matrix(tau/2)
		if self.decouple_scheme == 'XY4':
			expHD = Ye*expH*Xe*expH*Ye*expH*Xe*expH
		elif self.decouple_scheme == 'simple': # CHECK THIS
			expHD = expHhalf*Xe*expH*Ye*expH*Xe*expH*Xe*expHhalf
		else:
			raise Exception('Unknown scheme!')

		return expHD

	def nuclear_gate(self,number_of_pulses,tau):
		'''Gives the evolution matrix for number_of_pulses pulses'''

		return self.decouple_unit_evolution(tau) ** (number_of_pulses/2)
	   

	def c13_gate(self,carbon_nr, phase = None):
		''' THIS NEEDS WRITING ONCE DECIDE ON HOW TO DO PHASES'''
		# When phase_y is true, the phase gets the opposite sign, this is to simulate what we really do in the experiment - should be improved (JULIA)
		phase = -1* phase

	def phase_gate(self,carbon_nr, phase):
		' NEED TO DECIDE ON DECOUPLING'
		precession_freq = mp['C' + str(carbon_nr) +'_freq']*2*np.pi

		current_phase = total_time*precession_freq
		phase_dif = (phase-current_phase)%(2*np.pi)
		dec_time =  (phase_dif)/precession_freq

		tau = dec_time/4

		U0, U1 = nuclear_gate(2, tau, omega_Larmor, A_par, A_perp)

   
class NV_experiment(NV_system):

	'''Eventually this should contain code that will interpret a gate seq and make an experiment. '''


''' Here are different experiments that we commonly run on the system'''
def C13_fingerprint(NV_system,N = 32, tau_range =  np.arange(1e-6,7e-6,1e-7), show_plot = True):

		exp0 = np.zeros(np.shape(tau_range),dtype=np.complex_)
		
		proj0 = NV_system.embed_e_operator(rho0)

		for i,tau in enumerate(tau_range):
			gate_seq = NV_system.embed_e_operator(mx)*NV_system.nuclear_gate(N,tau)*NV_system.embed_e_operator(x)
			exp0[i] = (proj0 * evolve(NV_system.initial_state,gate_seq)).tr()

		### plotting
		if show_plot == True:

			plt.figure()
			plt.plot(tau_range*1e6,exp0)
			plt.title('Signal'); plt.xlabel('Tau')
			plt.ylim([0,1.1])


