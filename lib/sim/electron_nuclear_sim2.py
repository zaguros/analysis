# -*- coding: utf-8 -*-

''' A module to calculate the 13C nuclear and electron spin dynamics
under dynamical decoupling gates. By PH (inspiration from TH previous version) '''

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from importlib import reload

import numpy as np
import qutip
qutip = reload(qutip)
import analysis.lib.QEC.hyperfine_params as hf
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import warnings
import collections

import copy

# import hyperfine_params as hf_params; reload(hf_params)
# hf = hf_params.hyperfine_params

from scipy.linalg import fractional_matrix_power

#######################
### Basic functions ###
#######################

'''To do make all this a class!'''

def pauli():
	'''Define pauli spin matrices'''
	identity = qutip.qeye(2)
	sx = qutip.sigmax()/2
	sy = qutip.sigmay()/2
	sz = qutip.sigmaz()/2
	return identity, sx, sy, sz

def pauli_S1():
	'''Define pauli spin matrices for spin 1'''
	identity = qutip.qeye(3)
	sx = qutip.jmat(1,'x')
	sy = qutip.jmat(1,'y')
	sz = qutip.jmat(1,'z')

	szPseudoHalf = qutip.Qobj(-1.0*sz[1:,1:]) # Take the effective two level system bit out of the S1 sz
	return identity, sx, sy, sz, szPseudoHalf

def basic_spin_rotations():
	''' define some simple spin rotations'''
	X = (-1j*sx*np.pi).expm(); mX = (1j*sx*np.pi).expm()
	Y = (-1j*sy*np.pi).expm(); mY = (1j*sx*np.pi).expm()
	Z = (-1j*sz*np.pi).expm(); mZ = (1j*sz*np.pi).expm()
	x = (-1j*sx*np.pi/2).expm(); mx = (1j*sx*np.pi/2).expm()
	y = (-1j*sy*np.pi/2).expm(); my = (1j*sy*np.pi/2).expm()
	z = (-1j*sz*np.pi/2).expm(); mz = (1j*sz*np.pi/2).expm()
	return X,Y,Z,mX,mY,mZ,x,y,z,mx,my,mz

def spin_x_rotation(theta):
	return (-1j*sx*theta).expm()

def spin_y_rotation(theta):
	return (-1j*sy*theta).expm()

def spin_theta_rotation(phi, theta):
	return (-1j*(np.cos(phi)*sx + np.sin(phi)*sy)*theta).expm()

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

def basic_spin_states_S1():
	# Version for SU(1)
	rhom = qutip.qeye(3)/3
	ket0 = qutip.basis(3,1)
	bra0 = qutip.basis(3,1).dag()
	ket1 = qutip.basis(3,2)
	bra1 = qutip.basis(3,2).dag()
	ketm1 = qutip.basis(3,0)
	bram1 = qutip.basis(3,0).dag()
	rho0 = ket0*bra0
	rho1 = ket1*bra1
	rhom1 = ketm1*bram1

	return ket0, bra0, ket1, bra1, ketm1, bram1, rho0, rho1, rhom1, rhom

### create a set of usefull simple states and gates (should live in a class)
Id, sx, sy, sz = pauli()                                # Electron spin operators
Id_S1, sx_S1, sy_S1, sz_S1, szPseudo1_2 = pauli_S1()					# Spin 1
X,Y,Z,mX,mY,mZ,x,y,z,mx,my,mz = basic_spin_rotations()  # Basic gates
ket0, bra0, ket1, bra1, rho0, rho1, rhom, ketx,brax,ketmx,bramx,rhox,rhomx,kety,bray,ketmy,bramy,rhoy,rhomy = basic_spin_states() # Basic states

ket0_S1, bra0_S1, ket1_S1, bra1_S1, ketm1_S1, bram1_s1, rho0_S1, rho1_S1, rhom1_S1, rhom_S1 = basic_spin_states_S1()

gamma_c = 1.0705e3

def dyn_dec_signal(carbon_params,tau, N, sign = 1):
	''' Useful for quick investigation of fingerprints etc.
	'''

	if np.size(tau)!=1:
		M=np.zeros([len(carbon_params),np.size(tau)])
	else:
		M = np.zeros([len(carbon_params),np.size(N)])
	for i,carbon_param in enumerate(carbon_params):
		omega_larmor = carbon_param[0]
		HF_par = sign*carbon_param[1]
		HF_perp = sign*carbon_param[2]

		omega_tilde = np.sqrt((HF_par+omega_larmor)**2+HF_perp**2)

		alpha = omega_tilde*tau
		beta = omega_larmor*tau
		mx = HF_perp/omega_tilde
		mz = (HF_par+omega_larmor)/omega_tilde
		vec_term = mx**2 *((1-np.cos(alpha))*(1-np.cos(beta)))/(1+np.cos(alpha)*np.cos(beta)-mz*np.sin(alpha)*np.sin(beta))
		angle_term = np.sin(N*np.arccos(np.cos(alpha)*np.cos(beta)-mz*np.sin(alpha)*np.sin(beta))/2)**2

		M[i,:]= 1-(vec_term*angle_term)

	return M


###########################
### 	 Classes        ###
###########################


class NV_system(object):
	''' Basic class to contain the parameters of the NV system, plus some functions to simulate its evolution '''

	def __init__(self,**kw):

		self.B_field = kw.pop('B_field',414.1871869)
		self.gamma_c = 1.0705e3 #g-factor for C13 in Hz/G

		self.gamma_n = 0.31e3 # Nitrogen hamiltonian
		self.P_n = 5.04e6 # Nitrogen hamiltonian
		self.A_n = 2.19e6

		self.mw_ops = ['Xe','Ye','mXe','mYe','xe','ye','mxe','mye'] # list of defined rotation ops for cache and coding purposes

		self.mw_detuning = kw.pop('mw_detuning',0.0)

		self.espin_trans = kw.pop('espin_trans','+1')  # Changes sign of A_par and A_perp
		self.sign = -1 if self.espin_trans == '-1' else 1

		self.inc_nitrogen = kw.pop('inc_nitrogen',False)

		self.add_carbons(**kw)
		self.recalculate()

		self.cache_system_evn = True


	def set_mw_detuning(self,detuning):
		self.mw_detuning = detuning

	def add_carbons(self, **kw):

		if kw.pop('use_msmt_params',False):
			raise Exception('Not written yet!')
		else:
			carbon_params = kw.pop('carbon_params',False)
			if isinstance(carbon_params,list):
				self.num_carbons = len(carbon_params)
				self.carbon_params = []
				for carbon_param in carbon_params:
					self.carbon_params.append([2 * np.pi * self.B_field * self.gamma_c,2 * np.pi * carbon_param[0],2 * np.pi * carbon_param[1]])
				self.calc_c_prec_freqs()

			elif kw.pop('use_hf_library',False):

				self.carbon_params = []
				for C_key,C_val in sorted(hf.items()): # Sort to make sure alphabetical i.e. C1, C2 etc.
					if C_key == 'espin_trans':
						self.espin_trans = C_val
						self.sign = -1 if self.espin_trans == '-1' else 1
					else:
						self.carbon_params.append([2 * np.pi * self.B_field * self.gamma_c,2 * np.pi * C_val['par'],2 * np.pi * C_val['perp']])
				self.num_carbons = len(self.carbon_params)
				self.calc_c_prec_freqs()

			else:
				warnings.warn('No carbon params passed, using dummy params!')
				self.num_carbons = 1
				self.carbon_params = np.array([[2 * np.pi * 443180.29, 2 * np.pi * 35.0e3 , -2 * np.pi * 33.0e3]])
				self.calc_c_prec_freqs()

			self.recalc_Hamiltonian = True

	def calc_c_prec_freqs(self):
		''' If these arent specifed, set them from the carbon params '''
		sign = -1 if self.espin_trans == '-1' else 1

		freqs = []
		for i in range(self.num_carbons):
			omega0 = self.carbon_params[i][0]
			omega1 = np.sqrt((self.carbon_params[i][0] + sign*self.carbon_params[i][1])**2 + self.carbon_params[i][2]**2)
			omega_sup = 0.5*(omega0+omega1)
			freqs.append([omega0,omega1,omega_sup])

		self.c_prec_freqs = np.array(freqs)


	def e_op(self,operator):
		''' Helper function to embed e operator with identities on carbons '''
		return qutip.tensor([operator] + [Id] * self.num_carbons + [Id_S1] * self.inc_nitrogen)

	def c_op(self,operator,c_num):
		''' Helper function to embed c operator with identities on other carbons '''
		return qutip.tensor([Id] * c_num + [operator] + [Id] * (self.num_carbons - c_num) + [Id_S1] * self.inc_nitrogen)

	def N_op(self,N_operator):
		''' Helper function to embed N operator with identities on other carbons '''
		return qutip.tensor([Id] * (self.num_carbons + 1) + [N_operator])

	def e_C_op(self,e_operator,c_operator,c_num):
		''' Helper function to embed e operator combined with an operator on carbons '''
		return qutip.tensor([e_operator] + [Id] * (c_num-1) + [c_operator] + [Id] * (self.num_carbons - c_num) + [Id_S1] * self.inc_nitrogen)

	def e_N_op(self,e_operator,N_operator):
		''' Helper function to embed e operator combined with an operator on carbons '''
		return qutip.tensor([e_operator] + [Id] * self.num_carbons + [N_operator])

	def recalculate(self):
		self.cache_system_evn_taus = []
		self.cache_system_evn_Unitaries = []
		self.recalc_Hamiltonian = True
		self.define_useful_states()
		self._define_e_operators()

	def define_useful_states(self):
		''' standard init state for system '''
		self.NV0_carbons_mixed = qutip.tensor([rho0] + [rhom] * self.num_carbons + [rhom_S1] * self.inc_nitrogen)

	def _define_e_operators(self):
		''' Define commonly used electronic operations '''
		self._Xe = self.e_op(X)
		self._Ye = self.e_op(Y)
		self._mXe = self.e_op(mX)
		self._mYe = self.e_op(mY)
		self._xe = self.e_op(x)
		self._ye = self.e_op(y)
		self._mxe = self.e_op(mx)
		self._mye = self.e_op(my)
		self._Ide = self.e_op(Id)
		self._proj0 = self.e_op(rho0)
		self._proj1 = self.e_op(rho1)

		''' Trivial here but useful later '''
		self.Xe = lambda : self._Xe
		self.Ye = lambda : self._Ye
		self.mXe = lambda : self._mXe
		self.mYe = lambda : self._mYe
		self.xe = lambda : self._xe
		self.ye = lambda : self._ye
		self.mxe = lambda : self._mxe
		self.mye = lambda : self._mye

		self.Ide = lambda : self._Ide
		self.proj0 = lambda : self._proj0
		self.proj1 = lambda : self._proj1
		self.re = lambda theta,phi : self.e_op(spin_theta_rotation(theta, phi))


	def NV_carbon_system_Hamiltonian(self):
		''' Function to calculate the NV C13 system Hamiltonian '''

		if self.recalc_Hamiltonian == True:

			if self.num_carbons:
				self.Hsys = np.sum([self.e_C_op(rho0,carbon_param[0]*sz,i+1) \
							 + self.e_C_op(rho1,((carbon_param[0]+ self.sign*carbon_param[1])*sz + self.sign * carbon_param[2] * sx),i+1) \
					   for i,carbon_param in enumerate(self.carbon_params)]).tidyup()
			else:
				self.Hsys = 0

			if self.inc_nitrogen:

				self.Hsys += self.e_N_op(2*np.pi*self.A_n*self.sign*szPseudo1_2,sz_S1) + self.N_op(-2 * np.pi * (self.P_n*(sz_S1**2 -  1/3.0) + self.gamma_n*self.B_field*sz_S1))

			self.Hsys += self.e_op(2*np.pi*self.mw_detuning*self.sign*szPseudo1_2) # Note that funniness because NV is actually an S1 system..

			self.Hsys.tidyup()

			self.recalc_Hamiltonian = False

		return self.Hsys


	def NV_carbon_ev(self,tau):
		''' Function to calculate a C13 evolution matrix from the system Hamiltonian. Written this way so that could be overwritten'''

		''' By default, will cache evolution for a given tau, so that doesnt have to recalculate! '''
		if self.cache_system_evn:
			if tau in self.cache_system_evn_taus:
				return self.cache_system_evn_Unitaries[self.cache_system_evn_taus.index(tau)]
			else:
				unitary = (-1j*self.NV_carbon_system_Hamiltonian()*tau).expm()
				self.cache_system_evn_Unitaries.append(unitary)
				self.cache_system_evn_taus.append(tau)
				return unitary

		else:
			return (-1j*self.NV_carbon_system_Hamiltonian()*tau).expm()

# Helper function for sequences (does the actual calculation of the sequence output!)
def calc_sequence_operation(sequence):
	operation = 1.0
	for gate in sequence:
		if isinstance(gate[0],collections.deque): # Sequences can contain sequences!
			operation = calc_sequence_operation(gate[0]) ** gate[1] * operation
		else:
			operation = gate[0].gate_op() ** gate[1] * operation
	return operation

class gate(object):
	def __init__(self,gate_function,name=None,**kw):
		self.name = name
		self.gate_function = gate_function
		self.gate_properties = kw
	def gate_op(self):
		# Written this way so that could in principle mess with the properties after defined!
		# Maybe nuclear gates should have more of this functionality
		return self.gate_function(**self.gate_properties)


class basic_gate_sequence(object):

	def __init__(self,NV_system,**kw):
		self.NVsys = NV_system
		self.Ide = lambda : self.NVsys._Ide
		self.sequence = collections.deque()
		self._define_gates()

	def _reset_sequence(self):
		self.sequence = collections.deque()

	def add_gate_helper(self,gate_func,name=None,**kw):
		before = kw.pop('before', False)
		reps = kw.pop('reps',1)
		g =  gate(gate_func,name,**kw)
		self.add_gate_to_seq(g,before=before,reps=reps)

	def _define_gates(self):
		''' This is written this way so that could be overwritten for more complex behaviour'''
		self.Xe = lambda : self.add_gate_helper(self.NVsys.Xe,name ='Xe')
		self.Ye = lambda : self.add_gate_helper(self.NVsys.Ye,name ='Xe')
		self.mXe = lambda : self.add_gate_helper(self.NVsys.mXe,name ='Xe')
		self.mYe = lambda : self.add_gate_helper(self.NVsys.mYe,name ='Xe')
		self.xe = lambda : self.add_gate_helper(self.NVsys.xe,name ='Xe')
		self.ye = lambda : self.add_gate_helper(self.NVsys.ye,name ='Xe')
		self.mxe = lambda : self.add_gate_helper(self.NVsys.mxe,name ='Xe')
		self.mye = lambda : self.add_gate_helper(self.NVsys.mye,name ='Xe')

		self.proj0 = lambda **kw : self.add_gate_helper(self.NVsys.proj0,name='proj0',**kw)
		self.proj1 = lambda **kw : self.add_gate_helper(self.NVsys.proj0,name='proj1',**kw)

		self.re = lambda **kw: self.add_gate_helper(self.NVsys.re,**kw) # Note that need to pass theta and tau when calling this!

	def add_gate_to_seq(self,gate,reps=1,before=False):

		if isinstance(gate,basic_gate_sequence):
			gate = gate.sequence

		if before:
			self.sequence.appendleft([gate,reps])
		else:
			self.sequence.append([gate,reps])

		return self

	def seq_operation(self):
		return calc_sequence_operation(self.sequence)

	def apply_sequence(self,state,reps=1,norm = False):
		operation = self.seq_operation()**reps
		if not(isinstance(operation,float)):
			sysout = operation * state * operation.dag()
		else:
			sysout = state
		return sysout.unit() if norm else sysout

	def copy_seq(self):
		copied_seq = basic_gate_sequence(self.NVsys)
		copied_seq.sequence = copy.deepcopy(self.sequence)
		return copied_seq


class NV_gate_sequence(basic_gate_sequence):
# Add in more sophisticated gates (i.e. nuclear spin gates)

	def __init__(self,NV_system,**kw):
		basic_gate_sequence.__init__(self,NV_system,**kw)
		self.decouple_scheme = kw.pop('decouple_scheme','XY4')

	def copy_seq(self):
		copied_seq = NV_gate_sequence(self.NVsys)
		copied_seq.sequence = copy.deepcopy(self.sequence)
		return copied_seq

	def nuclear_ev_gate(self,in_tau,tau_factor=1.0,double_sided = False):

		if callable(in_tau):
			tau = in_tau()
		else:
			tau = in_tau

		return self.NVsys.NV_carbon_ev(self.nuclear_gate_tau(tau_factor*tau, double_sided = double_sided))

	def nuclear_gate(self,N,tau,**kw):

		'''Evolution during a decouple unit'''
		if N == 0:
			return

		# Note that these are functions, so that evaluated when the gate sequence is evaluated!
		evNV_C_tau =  gate(lambda : self.nuclear_ev_gate(tau, double_sided = True),'tau')
		evNV_C_tau_single =  gate(lambda : self.nuclear_ev_gate(tau),'tau_single')
		evNV_C_2tau =  gate(lambda : self.nuclear_ev_gate(tau, tau_factor=2.0, double_sided = True),'2_tau')

		scheme = kw.pop('scheme',self.decouple_scheme)

		seq = basic_gate_sequence(self.NVsys)

		if scheme == 'XY4':
			if N%4 != 0:
				raise Exception('Incompatible number of pulses!')

			seq.add_gate_to_seq(evNV_C_tau)
			seq.Xe()
			seq.add_gate_to_seq(evNV_C_tau_single)

			seq_repeat = basic_gate_sequence(self.NVsys)
			seq_repeat.add_gate_to_seq(evNV_C_tau_single)
			seq_repeat.Ye()
			seq_repeat.add_gate_to_seq(evNV_C_2tau)
			seq_repeat.Xe()
			seq_repeat.add_gate_to_seq(evNV_C_tau_single)
			seq.add_gate_to_seq(seq_repeat,reps=(N/2-1))

			seq.add_gate_to_seq(evNV_C_tau_single)
			seq.Ye()
			seq.add_gate_to_seq(evNV_C_tau)

		elif scheme == 'XY8':

			if N%8 != 0:
				raise Exception('Incompatible number of pulses!')

			seq.add_gate_to_seq(evNV_C_tau)
			seq.Xe()
			seq.add_gate_to_seq(evNV_C_2tau)
			seq.Ye()
			seq.add_gate_to_seq(evNV_C_2tau)
			seq.Xe()
			seq.add_gate_to_seq(evNV_C_2tau)
			seq.Ye()
			seq.add_gate_to_seq(evNV_C_tau_single)

			seq_repeat = basic_gate_sequence(self.NVsys)
			seq_repeata = basic_gate_sequence(self.NVsys)
			seq_repeata.add_gate_to_seq(evNV_C_tau_single)
			seq_repeata.Ye()
			seq_repeata.add_gate_to_seq(evNV_C_2tau)
			seq_repeata.Xe()
			seq_repeata.add_gate_to_seq(evNV_C_tau_single)

			seq_repeatb = basic_gate_sequence(self.NVsys)
			seq_repeatb.add_gate_to_seq(evNV_C_tau_single)
			seq_repeatb.Xe()
			seq_repeatb.add_gate_to_seq(evNV_C_2tau)
			seq_repeatb.Ye()
			seq_repeatb.add_gate_to_seq(evNV_C_tau_single)

			seq_repeat.add_gate_to_seq(seq_repeata,reps=2).add_gate_to_seq(seq_repeatb,reps=2)
			seq.add_gate_to_seq(seq_repeat,reps = (N/8-1))

			seq.add_gate_to_seq(evNV_C_tau)
			seq.Ye()
			seq.add_gate_to_seq(evNV_C_2tau)
			seq.Xe()
			seq.add_gate_to_seq(evNV_C_2tau)
			seq.Ye()
			seq.add_gate_to_seq(evNV_C_2tau)
			seq.Xe()
			seq.add_gate_to_seq(evNV_C_tau_single)

		elif scheme == 'simple': # CHECK THIS
			if N<2:
				seq.add_gate_to_seq(evNV_C_tau)
				seq.Xe()
				seq.add_gate_to_seq(evNV_C_tau)
			else:
				seq.add_gate_to_seq(evNV_C_tau)
				seq.Xe()
				seq.add_gate_to_seq(evNV_C_tau_single)

				seq_repeat = basic_gate_sequence(self.NVsys)
				seq_repeat.add_gate_to_seq(evNV_C_tau_single)
				seq_repeat.Xe()
				seq_repeat.add_gate_to_seq(evNV_C_tau_single)
				seq.add_gate_to_seq(seq_repeat,reps=(N-2))

				seq.add_gate_to_seq(evNV_C_tau_single)
				seq.Xe()
				seq.add_gate_to_seq(evNV_C_tau)

		else:
			raise Exception('Unknown scheme!')

		self.add_gate_to_seq(seq,**kw)

	def nuclear_gate_tau(self,tau,double_sided = False):

		'''Helper function to get tau for gate seq'''
		tau_correction_factor = self.NVsys.tau_correction_factor if hasattr(self.NVsys,'tau_correction_factor') else 0

		scale_fact = 0.5 if not(double_sided) else 1.0

		if tau_correction_factor > scale_fact*tau:
			raise Exception('mw_duration too long!')
		return (tau - scale_fact*tau_correction_factor)


	def wait_gate(self,tau,**kw):
		''' Do nothing! '''
		self.add_gate_to_seq(gate(lambda : self.nuclear_ev_gate(tau),'wait_gate'),**kw)

	def nuclear_phase_gate(self,carbon_nr, phase, state = 'sup',**kw):

		''' For small waits (no decoupling) '''
		if state == 'sup':
			precession_freq = self.NVsys.c_prec_freqs[carbon_nr-1,2]
		elif state == 0:
			precession_freq = self.NVsys.c_prec_freqs[carbon_nr-1,0]
		elif state == 1:
			precession_freq = self.NVsys.c_prec_freqs[carbon_nr-1,1]
		# NEED TO CHECK HOW PHASE IS DEFINED IN OUR EXPM
		phase = (np.pi*float(-1*phase)/180)%(2*np.pi)
		dec_time = phase/precession_freq

		self.wait_gate(dec_time,**kw)

	def mbi_sequence(self,N,tau):

		self.ye()
		self.nuclear_gate(N ,tau)
		self.mxe()



class NV_experiment(object):

	def __init__(self,NV_system,**kw):
		self.NVsys = NV_system
		self.specified_initial_state = None
		self.reset_output_state()

	def initial_state(self):
		if self.specified_initial_state is None:
			return self.NVsys.NV0_carbons_mixed
		else:
			return self.specified_initial_state

	def reset_init_state(self, state = None):
		self.specified_initial_state = copy.deepcopy(state)
		self.reset_output_state()

	def reset_output_state(self):
		self.output_state = self.initial_state()

	def gate_sequence(self): # Convenience function to quickly get a new gate sequence
		return NV_gate_sequence(self.NVsys)

	def apply_gates(self,gate_sequence,**kw):
		self.output_state = gate_sequence.apply_sequence(self.output_state, **kw)
		return self

	def measure_e(self,e_state = 0):
		if e_state == 0:
			e_state = self.NVsys.e_op(rho0)
		elif e_state == 1:
			e_state = self.NVsys.e_op(rho1)

		return np.real((e_state*self.output_state).tr())

	def measure_c(self,c_num=1,c_state = 0):
		''' Not directly accessible, but sometimes useful'''
		if c_state == 0:
			c_state = rho0
		elif c_state == 1:
			c_state = rh1

		proj = self.NVsys.c_op(c_state,c_num)

		return np.real((proj*self.output_state).tr())

	def measure_N(self,N_state = 0):
		''' Not directly accessible, but sometimes useful'''
		if N_state == 0:
			N_state = rho0_S1
		elif N_state == 1:
			N_state = rho1_S1
		elif N_state == -1:
			N_state = rhom1_S1

		proj = self.NVsys.N_op(N_state)

		return np.real((proj*self.output_state).tr())

def gaussian_envelope(t,duration):
	T_herm = 0.1667*duration

	return (1 - 0.956 * ((t- duration/2)/T_herm)**2) * np.exp(-((t - duration/2)/T_herm)**2)


class noisy_NV_system(NV_system):

	def __init__(self,**kw):

		self.mean_amp = kw.pop('mean_amp',1.0)
		self.mw_duration = kw.pop('mw_duration',10.0e-9)
		self.tau_correction_factor = self.mw_duration
		self.pulse_shape = kw.pop('pulse_shape','square')



		NV_system.__init__(self,**kw)
		self.recalculate()

	def set_mw_duration(self,mw_duration):
		self.mw_duration = mw_duration
		self.tau_correction_factor = self.mw_duration
		self.recalculate()

	def set_mw_amp(self,amp):
		self.mean_amp = amp
		self.reset_caches()

	def set_mw_detuning(self,detuning):
		self.mw_detuning = detuning
		self.reset_caches()

	def recalculate(self):
		self.recalc_Hamiltonian = True
		self.define_useful_states()
		self._define_e_operators()
		self.reset_caches()

	def amp_val(self):
		# Could do more complicated things if you want!
		return self.mean_amp

	def finite_microwave_pulse(self,duration,theta,phi,steps=20):

		duration = np.float(duration)

		if self.pulse_shape == 'square':
			Hsys = self.NV_carbon_system_Hamiltonian()
			Hint = self.e_op(phi*(np.cos(theta)*sx + np.sin(theta)*sy))

			return (-1j*(duration*Hsys+Hint)).expm()

		elif self.pulse_shape == 'Hermite':

			Hsys = self.NV_carbon_system_Hamiltonian()

			dt = duration/steps
			t = np.arange(0+dt/2,duration,dt)
			normfactor = steps/np.sum(gaussian_envelope(t,duration))
			Hint = (normfactor/duration)*phi*self.e_op((np.cos(theta)*sx + np.sin(theta)*sy))

			Utots = [(-1j*(Hsys+gaussian_envelope(ts,duration)*Hint)*dt).expm() for ts in t]
			return qutip.gate_sequence_product(Utots)

	def reset_caches(self):
		for op_string in self.mw_ops:
			setattr(self, op_string + '_cache_recalc', True)

		self.cache_system_evn_taus = []
		self.cache_system_evn_Unitaries = []

	def calc_unitary_trans(self,op_string):
		#Here is some code cleverness to only recalculate the matrices if the mw_amp changes
		if hasattr(self,op_string + '_cache'):
			if not(eval('self.' + op_string + '_cache_recalc')):
				return eval('self.' + op_string + '_cache')

		raw_gate_op = eval('self._' + op_string)

		# Add in the amplitude value
		if self.amp_val != 1.0:
			gate_op = qutip.Qobj(fractional_matrix_power(raw_gate_op.full(), self.amp_val()),dims=raw_gate_op.dims)
		else:
			gate_op = raw_gate_op
		gate_op.tidyup()

		# Add to cache
		setattr(self, op_string + '_cache', gate_op)
		setattr(self, op_string + '_cache_recalc', False)

		return gate_op

	def _define_e_operators(self):
		''' Override commonly used electronic gates '''
		self._Ide = self.e_op(Id)
		self._proj0 = self.e_op(rho0)
		self._proj1 = self.e_op(rho1)

		self._Xe = self.finite_microwave_pulse(self.mw_duration,0.0,np.pi)
		self._Ye = self.finite_microwave_pulse(self.mw_duration,0.5*np.pi,np.pi)
		self._mXe = self.finite_microwave_pulse(self.mw_duration,0.0,-np.pi)
		self._mYe = self.finite_microwave_pulse(self.mw_duration,0.5*np.pi,-np.pi)
		self._xe = self.finite_microwave_pulse(self.mw_duration,0.0,0.5*np.pi)
		self._ye = self.finite_microwave_pulse(self.mw_duration,0.5*np.pi,0.5*np.pi)
		self._mxe = self.finite_microwave_pulse(self.mw_duration,0.0,-0.5*np.pi)
		self._mye = self.finite_microwave_pulse(self.mw_duration,0.5*np.pi,-0.5*np.pi)

		for op_string in self.mw_ops:
			setattr(self,op_string, lambda op_string = op_string: self.calc_unitary_trans(op_string))  # Force eval of op_string at defn time
		self.Ide = lambda : self._Ide
		self.proj0 = lambda : self._proj0
		self.proj1 = lambda : self._proj1
		self.re = lambda theta,phi : self.finite_microwave_pulse(self.mw_duration,theta,phi*self.amp_val()) # Not cached!



#########################################
#										#
#				EXPERIMENTS				#
#										#
#########################################

''' Here are different experiments that we commonly run on the system'''


def C13_fingerprint(NV_system,N = 32, tau_range =  np.arange(1e-6,7e-6,1e-7), calc_indiv = True, quick_calc =False):
	''' Simple experiment sweeping tau for a fixed N and measuring whether e still in the same state '''
	if not(quick_calc):
		if calc_indiv:

			exp0 = np.zeros((np.shape(tau_range)[0],NV_system.num_carbons))

			carbon_params = NV_system.carbon_params
			c_prec_freqs = NV_system.c_prec_freqs


			nv_expm = NV_experiment(NV_system)
			gate_seq = nv_expm.gate_sequence()
			gate_seq.xe(), gate_seq.nuclear_gate(N ,lambda : tau), gate_seq.mxe() # can define tau later! Cool huh


			for j, carbon_param in enumerate(carbon_params):

				NV_system.carbon_params = [carbon_param]
				NV_system.num_carbons = 1
				NV_system.c_prec_freqs = [c_prec_freqs[j]]
				NV_system.recalculate()

				for i,tau in enumerate(tau_range):

					nv_expm.reset_output_state()
					nv_expm.apply_gates(gate_seq)
					exp0[i,j] = nv_expm.measure_e()
			# Reset to old values
			NV_system.carbon_params = carbon_params
			NV_system.num_carbons = len(carbon_params)
			NV_system.c_prec_freqs = c_prec_freqs
			NV_system.recalculate()

		else:

			exp0 = np.zeros((np.shape(tau_range)[0]))
			nv_expm = NV_experiment(NV_system)
			gate_seq = nv_expm.gate_sequence()
			gate_seq.xe(), gate_seq.nuclear_gate(N ,lambda : tau), gate_seq.mxe()


			for i,tau in enumerate(tau_range):

				nv_expm.reset_output_state()
				gate_seq = nv_expm.gate_sequence()
				nv_expm.apply_gates(gate_seq)
				exp0[i] = nv_expm.measure_e()


	else:
		exp0 = 0.5*(1+dyn_dec_signal(NV_system.carbon_params, tau_range, N,sign = NV_system.sign)).T

	width = 12
	height = 4
	plt.figure(figsize=(width, height))
	plt.plot(tau_range*1e6,exp0)
	plt.title('Signal for N=%s'%(N)); plt.xlabel('Tau')
	plt.ylim([-0.1,1.1])
	plt.show()
	plt.close()




def sweep_MW_amp(noisy_NV_system,N = 11, amp_range =  np.arange(0.1,2,0.05), tau = 7.5e-6,**kw):

	nv_expm = NV_experiment(noisy_NV_system)
	gate_seq = nv_expm.gate_sequence()
	gate_seq.nuclear_gate(N ,tau, scheme = 'simple')

	results = np.zeros(np.shape(amp_range))

	for i,amp in enumerate(amp_range):

		noisy_NV_system.set_mw_amp(amp)
		nv_expm.apply_gates(gate_seq)
		results[i] = nv_expm.measure_e()
		nv_expm.reset_output_state()


	plt.figure()
	plt.plot(amp_range,results)
	plt.title('Signal'); plt.xlabel('MW amplitude')
	plt.ylim(bottom = 0)
	plt.show()
	plt.close()

	noisy_NV_system.set_mw_amp(1)

	ind = np.argmin(results)
	print('Min sig. ', results[ind], ' at ', amp_range[ind])


def sweep_MW_duration(noisy_NV_system,N = 11, duration_range =  np.arange(50,200,10)*1e-9, tau = 7.5e-6,**kw):

	nv_expm = NV_experiment(noisy_NV_system)
	gate_seq = nv_expm.gate_sequence()
	gate_seq.nuclear_gate(N ,tau, scheme = 'simple')

	results = np.zeros(np.shape(duration_range))

	for i,duration in enumerate(duration_range):

		noisy_NV_system.set_mw_duration(duration)
		nv_expm.apply_gates(gate_seq)
		results[i] = nv_expm.measure_e()
		nv_expm.reset_output_state()


	plt.figure()
	plt.plot(duration_range*1e9,results)
	plt.title('Signal'); plt.xlabel('MW duration (ns)')
	plt.ylim(bottom = 0)
	plt.show()
	plt.close()


def MonteCarlo_MWFid(noisy_NV_system,N = 11, tau = 7.5e-6,N_rand = 100,mean = 1.0,sigma=0.01):
	'''Simulate doing microwave pulses with a certain standard deviation on the pulse amplitude from trial to trial '''

	nv_expm = NV_experiment(noisy_NV_system)
	gate_seq = nv_expm.gate_sequence()
	gate_seq.nuclear_gate(N ,tau, scheme = 'simple')

	rands = np.random.normal(loc = mean,scale=sigma, size=N_rand)
	infids = np.zeros(N_rand)

	for i, rand_amp in enumerate(rands):

		noisy_NV_system.set_mw_amp(rand_amp)
		nv_expm.reset_output_state()
		nv_expm.apply_gates(gate_seq)
		infids[i] = nv_expm.measure_e()


	print("Infidelity is %f \pm %f" % (np.mean(infids), np.std(infids)/np.sqrt(N_rand)))

	return infids


def dynamical_decouple(NV_system,N_range = range(0,3000,32), tau = None,**kw):

	if tau is None:
		tau = 1/(NV_system.B_field * NV_system.gamma_c)

	scheme = kw.pop("scheme", "XY8") # Note that simple means that can constructively get oscillations from Nitrogen coupling during finite pulse duration..

	nv_expm = NV_experiment(NV_system)

	results = np.zeros(np.shape(N_range))

	for i,N in enumerate(N_range):

		nv_expm.reset_output_state()
		gate_seq = nv_expm.gate_sequence()
		gate_seq.nuclear_gate(N ,tau, scheme = scheme) # Currently cant change N on the fly.
		nv_expm.apply_gates(gate_seq,norm=True)
		results[i] = nv_expm.measure_e(0)



	plt.figure()
	plt.plot(N_range,results)
	plt.title('Signal'); plt.xlabel('N')
	plt.ylim(bottom = 0,top=1.1)
	plt.show()
	plt.close()


def e_ramsey(NV_system,delay_range =  np.arange(0e-9,5e-6,50e-9)):
	''' Prepare e in X (or attempt to) and measure in X or Y '''
	results = np.zeros(np.shape(delay_range))

	nv_expm = NV_experiment(NV_system)
	ramsey_seq = nv_expm.gate_sequence()
	ramsey_seq.xe()
	ramsey_seq.wait_gate(lambda: tau)
	ramsey_seq.mxe()

	for i,tau in enumerate(delay_range):

		nv_expm.apply_gates(ramsey_seq)
		results[i] = nv_expm.measure_e()
		nv_expm.reset_output_state()

	plt.figure()
	plt.plot(delay_range*1e6,results)
	plt.title('Signal'); plt.xlabel('Tau')
	plt.show()
	plt.close()


def hahn_echo(NV_system,delay_range =  np.arange(0e-9,10e-6,10e-9)):
	''' Prepare e in X (or attempt to) and measure in X or Y '''
	results = np.zeros(np.shape(delay_range))

	nv_expm = NV_experiment(NV_system)
	ramsey_seq = nv_expm.gate_sequence()
	ramsey_seq.xe()
	ramsey_seq.wait_gate(lambda: tau)
	ramsey_seq.Ye()
	ramsey_seq.wait_gate(lambda: tau)
	ramsey_seq.mxe()

	for i,tau in enumerate(delay_range):

		nv_expm.apply_gates(ramsey_seq)
		results[i] = nv_expm.measure_e()
		nv_expm.reset_output_state()

	plt.figure()
	plt.plot(delay_range*1e6,results)
	plt.title('Signal'); plt.xlabel('Tau')
	plt.show()
	plt.close()

	ind = np.argmin(results)
	print('Min sig. ', results[ind], ' at ', delay_range[ind]*1e6)

def dark_esr(noisy_NV_system,freq_range =  np.arange(-5e6,5e6,1e5)):

	results = np.zeros(np.shape(freq_range))
	nv_expm = NV_experiment(noisy_NV_system)
	desr_seq = nv_expm.gate_sequence()
	desr_seq.re(theta = 0,phi=1.0)

	for i,freq in enumerate(freq_range):

		noisy_NV_system.set_mw_detuning(freq)
		noisy_NV_system.recalculate()

		nv_expm.apply_gates(desr_seq)
		results[i] = nv_expm.measure_e()
		nv_expm.reset_output_state()

	plt.figure()
	plt.plot(freq_range*1e-6,results)
	plt.title('Signal'); plt.xlabel('Freq')
	plt.show()
	plt.close()


def prepare_X_and_measure_XY(NV_system,N = 32, tau_range =  np.arange(1e-6,7e-6,1e-7),meas = 'eXY',**kw):
	''' Prepare carbon in X (or attempt to) and measure in X and Y '''
	X = np.zeros(np.shape(tau_range))
	Y = np.zeros(np.shape(tau_range))

	nv_expm = NV_experiment(NV_system)

	mbi_seq = nv_expm.gate_sequence()
	mbi_seq.mbi_sequence(N,lambda: tau)

	init_seq = mbi_seq.copy_seq()
	init_seq.proj0()

	mbi_seq_plus_90 = mbi_seq.copy_seq()
	mbi_seq_plus_90.nuclear_phase_gate(1,90,state=0,before=True)

	c_num = kw.pop('c_num',1)

	for i,tau in enumerate(tau_range):

		nv_expm.apply_gates(init_seq, norm = True)
		nv_expm.reset_init_state(state = nv_expm.output_state)

		if meas == 'eXY':
			nv_expm.apply_gates(mbi_seq)
			X[i] = nv_expm.measure_e()
			nv_expm.reset_output_state()

			nv_expm.apply_gates(mbi_seq_plus_90)
			Y[i] = nv_expm.measure_e()
			nv_expm.reset_output_state()

		elif meas == 'nXY':

			X[i] = nv_expm.measure_c(c_state = rhox,c_num=c_num)
			Y[i] = nv_expm.measure_c(c_state = rhoy,c_num=c_num)

		nv_expm.reset_init_state()

	Fid = (np.sqrt((X-0.5)**2 + (Y-0.5)**2)+0.5)


	plt.figure()
	plt.plot(tau_range*1e6,X,label = 'X')
	plt.plot(tau_range*1e6,Y,label = 'Y')
	plt.plot(tau_range*1e6,Fid,label = 'Fid.')
	plt.title('Signal'); plt.xlabel('Tau')
	plt.legend()
	plt.ylim([0,1.1])
	plt.show()
	plt.close()

	ind = np.argmax(Fid)
	print('Max fid. ', Fid[ind], ' at ', tau_range[ind]*1e6)


def MonteCarlo_MWAmp_CGate_fid(noisy_NV_system,N = 32, tau = 6.582e-6,N_rand = 100,mean = 0.995,sigma=0.01,meas = 'eXY'):
	'''Simulate doing a carbon gate with finite microwave durations and a certain standard deviation on the pulse amplitude from trial to trial '''

	rands = np.random.normal(loc = mean,scale=sigma, size=N_rand)
	infids = np.zeros(N_rand)

	nv_expm = NV_experiment(noisy_NV_system)

	mbi_seq = nv_expm.gate_sequence()
	mbi_seq.mbi_sequence(N,tau)

	init_seq = mbi_seq.copy_seq()
	init_seq.proj0()

	mbi_seq_plus_90 = mbi_seq.copy_seq()
	mbi_seq_plus_90.nuclear_phase_gate(1,90,state=0,before=True)

	for i, rand_amp in enumerate(rands):

		noisy_NV_system.set_mw_amp(rand_amp)

		nv_expm.apply_gates(init_seq, norm = True)
		nv_expm.reset_init_state(state = nv_expm.output_state)

		if meas == 'eXY':
			nv_expm.apply_gates(mbi_seq)
			X = nv_expm.measure_e()
			nv_expm.reset_output_state()

			nv_expm.apply_gates(mbi_seq_plus_90)
			Y = nv_expm.measure_e()
			nv_expm.reset_output_state()

		elif meas == 'nXY':

			X = nv_expm.measure_c(c_state = rhox,c_num=c_num)
			Y = nv_expm.measure_c(c_state = rhoy,c_num=c_num)

		nv_expm.reset_init_state()


		infids[i] = (np.sqrt((X-0.5)**2 + (Y-0.5)**2)+0.5)

	print("Fidelity is %f \pm %f" % (np.mean(infids), np.std(infids)))

	return infids




# 		elif scheme == 'simple': # CHECK THIS
# 			gate = self.compile([evNV_C_tau,self.Xe(),evNV_C_tau],reps = N)
# 		else:
# 			raise Exception('Unknown scheme!')

# 		return gate


# The QUTIP solvers crash on this computer :(
# This was code to do Hamiltonian evolution for the qutip stuff.

# ''' define some simple spin rotations in terms of their HAMILTONIANS'''
# 		self._hXe = self.e_op(sx*np.pi)
# 		self._hYe = self.e_op(sy*np.pi)
# 		self._hxe = self.e_op(sx*np.pi/2)
# 		self._hye = self.e_op(sy*np.pi/2)
# 		self._hmxe = self.e_op(sx*np.pi/2)
# 		self._hmye = self.e_op(sy*np.pi/2)

# def assemble_ev_mat(self,gate_comps,comp_defs):
# 		'''Helper function for gate assembly'''

# 		ev_mat = self.e_op(Id)
# 		for comp in gate_comps:
# 			ev_mat = ev_mat*comp_defs[comp]
# 		return ev_mat

# 	def assemble_Ham_list(self,gate_comps,comp_defs):
# 		'''Helper function for gate assembly'''
# 		Ham_list = []
# 		for comp in gate_comps:
# 			Ham_list.append(comp_defs[comp])
# 		return Ham_list


#
# def Hamiltonian_callback(t,args):
# 		'''To make things simple, we assume each H entry evolves the system for unit time, and so speed up the interaction strength'''

# 		return args['Ham_list']
# def nuclear_gate(self,N,tau,**kw):
# 	'''Evolution during a decouple unit'''

# 	scheme = kw.pop('scheme',self.decouple_scheme)
# 	return_option = kw.pop('return_option','ev_mat')

# 	tau_gate = self.nuclear_gate_tau(tau)

# 	if scheme == 'XY4':
# 		if N%4 != 0:
# 			raise Exception('Incompatible number of pulses!')

# 		gate_comps = ['sys_ev','Ye','sys_ev_2','Xe','sys_ev_2','Ye','sys_ev_2','Xe','sys_ev'] * (N/4)

# 	elif scheme == 'simple': # CHECK THIS
# 		dec_Hams = ['sys_ev','Xe','sys_ev'] * N

# 	else:
# 		raise Exception('Unknown scheme!')

# 	if return_option == 'Ham_list':

# 		comp_defs = {}
# 		comp_defs['sys_ev'] = tau_gate * self.NV_carbon_system_Hamiltonian()
# 		comp_defs['sys_ev_2'] = comp_defs['sys_ev']*2
# 		comp_defs['Xe'] = self.hXe
# 		comp_defs['Ye'] = self.hYe
# 		return self.assemble_Ham_list(gate_comps,comp_defs)

# 	elif return_option == 'ev_mat' or return_option == 'ev_mat_raw':

# 		comp_defs = {}
# 		comp_defs['sys_ev'] = self.NV_carbon_evolution_matrix(tau_gate)
# 		comp_defs['sys_ev_2'] = comp_defs['sys_ev']**2
# 		comp_defs['Xe'] = self.Xe
# 		comp_defs['Ye'] = self.Ye

# 		if return_option == 'ev_mat':
# 			return self.assemble_ev_mat(gate_comps,comp_defs)
# 		else:
# 			return gate_comps,comp_defs

# def nuclear_gate_tau(self,tau):
# 	tau_correction_factor = self.tau_correction_factor if hasattr(self,'tau_correction_factor') else 0
# 	if tau_correction_factor > tau:
# 		raise Exception('mw_duration too long!')
# 	return 0.5*(tau - tau_correction_factor)

# def assemble_ev_mat(self,gate_comps,comp_defs):
# 	'''Helper function for gate assembly'''

# 	ev_mat = self.e_op(Id)
# 	for comp in gate_comps:
# 		ev_mat = ev_mat*comp_defs[comp]
# 	return ev_mat

# def assemble_Ham_list(self,gate_comps,comp_defs):
# 	'''Helper function for gate assembly'''
# 	Ham_list = []
# 	for comp in gate_comps:
# 		Ham_list.append(comp_defs[comp])
# 	return Ham_list

# def C13_fingerprint_ham_ev(NV_system,N = 32, tau_range =  np.arange(1e-6,7e-6,1e-7)):

# 		exp0 = np.zeros(np.shape(tau_range))

# 		proj0 = NV_system.e_op(rho0)

# 		# tlist = np.array([float(np.shape(Ham_list)[0])])
# 		tlist = np.linspace(0,5,200)
# 		# for i,tau in enumerate(tlist):

# 		# Ham_list = [NV_system.hmxe] + NV_system.nuclear_gate(N ,tau,return_option = 'Ham_list')  + [NV_system.xe]

# 		args={'Ham_list' : qutip.sigmax().data}
# 		exp0 = qutip.mesolve(Hamiltonian_callback , proj0,tlist,[],qutip.sigmaz(),args = args).expect[0]

# 		plt.figure()
# 		plt.plot(tlist*1e6,exp0)
# 		plt.title('Signal'); plt.xlabel('Tau')
# 		plt.ylim([0,1.1])
# 		plt.show()
# 		plt.close()



# def composite_gate_calc_p_q_for_f3DD(f3DD):

# 		q1 = 4.0/(np.sqrt(5.0 + np.pi * f3DD) + -1.0)
# 		q2 = 4.0/(np.sqrt(5.0 + np.pi * f3DD) + 1.0)

#  		p = 20.0*(1/4.0 - 1.0/(2.0 * np.pi) * np.arctan(np.sqrt(q1**2 - 1.0)))
#  		q = 20.0*(1/4.0 - 1.0/(2.0 * np.pi) * np.arctan(np.sqrt(q2**2 - 1.0)))
#  		return p,q

# class composite_gate_NV_system(NV_system):

# 	def __init__(self,**kw):

# 		self.composite_scheme = kw.pop('composite_scheme','none')
# 		self.r = kw.pop('r',0.31) # For 3 pulse scheme
# 		self.f3DD = kw.pop('f3DD',0.1) # For 5 pulse scheme

# 		self.p, self.q = composite_gate_calc_p_q_for_f3DD(self.f3DD)

# 		NV_system.__init__(self,**kw)


#  	def nuclear_gate(self,N,tau,**kw):
# 		'''Evolution during a decouple unit'''

# 		if self.composite_scheme == 'none':
# 			return NV_system.nuclear_gate(self,N,tau,**kw)

# 		if self.composite_scheme == '3_pulse':

# 			if N%3 != 0:
# 				raise Exception('Incompatible number of pulses!')

# 			evNV_C_r = self.NV_carbon_ev(self.nuclear_gate_tau(6*tau*(0.5-self.r)))
# 			evNV_C_t = self.NV_carbon_ev(self.nuclear_gate_tau(6*tau*self.r,double_sided = True))

# 			gate = self.compile([evNV_C_r,self.Xe(),evNV_C_t,self.Xe(),evNV_C_t,self.Xe(),evNV_C_r],reps = N/3)
# 			return gate

# 		elif self.composite_scheme == '5_pulse_simple_pulses':

# 			if N%8 != 0:
# 				raise Exception('Incompatible number of pulses!')

# 			pulse_tau = tau/5.0
# 			evNV_C_q = self.NV_carbon_ev(self.nuclear_gate_tau(pulse_tau*self.p))
# 			evNV_C_qp = self.NV_carbon_ev(self.nuclear_gate_tau(pulse_tau*(self.q-self.p),double_sided = True))
# 			evNV_C_p = self.NV_carbon_ev(self.nuclear_gate_tau(pulse_tau*(5.0-self.q),double_sided = True))
# 			pulse30 = self.re(0,np.pi)
# 			pulse120 = self.re(0,np.pi)

# 			Xgate = self.compile([evNV_C_q,pulse30,evNV_C_qp,self.Xe(),evNV_C_p,self.Xe(),evNV_C_p,self.Xe(),evNV_C_qp,pulse30,evNV_C_q])
# 			Ygate = self.compile([evNV_C_q,pulse120,evNV_C_qp,self.Xe(),evNV_C_p,self.Xe(),evNV_C_p,self.Xe(),evNV_C_qp,pulse30,evNV_C_q])

# 			gate = self.compile([Xgate,Ygate,Xgate,Ygate,Ygate,Xgate,Ygate,Xgate],reps = N/8)
# 			return gate

# 		elif self.composite_scheme == '5_pulse':

# 			if N%8 != 0:
# 				raise Exception('Incompatible number of pulses!')

# 			pulse_tau = tau/5.0
# 			evNV_C_q = self.NV_carbon_ev(self.nuclear_gate_tau(pulse_tau*self.p))
# 			evNV_C_qp = self.NV_carbon_ev(self.nuclear_gate_tau(pulse_tau*(self.q-self.p),double_sided = True))
# 			evNV_C_p = self.NV_carbon_ev(self.nuclear_gate_tau(pulse_tau*(5.0-self.q),double_sided = True))
# 			pulse30 = self.re(np.pi/3,np.pi)
# 			pulse120 = self.re(np.pi/2 + np.pi/3,np.pi)

# 			Xgate = self.compile([evNV_C_q,pulse30,evNV_C_qp,self.Xe(),evNV_C_p,self.Ye(),evNV_C_p,self.Xe(),evNV_C_qp,pulse30,evNV_C_q])
# 			Ygate = self.compile([evNV_C_q,pulse120,evNV_C_qp,self.Ye(),evNV_C_p,self.mXe(),evNV_C_p,self.Xe(),evNV_C_qp,pulse30,evNV_C_q])

# 			gate = self.compile([Xgate,Ygate,Xgate,Ygate,Ygate,Xgate,Ygate,Xgate],reps = N/8)
# 			return gate

# class NV_experiment(NV_system):

# 	'''Eventually this should contain code that will interpret a gate seq and make an experiment. '''

# 			# current_phase = total_time*precession_freq
# 		# phase_dif = (phase-current_phase)%(2*np.pi)
# 		# dec_time =  (phase_dif)/precession_freq

# 		# tau = dec_time/4
