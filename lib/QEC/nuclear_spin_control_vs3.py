''' A module to calculate the 13C nuclear and electron spin dynamics
under dynamical decoupling gates. By THT '''

import numpy as np
import qutip
qutip = reload(qutip)
import analysis.lib.QEC.hyperfine_params as hf
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import warnings

import hyperfine_params as hf_params; reload(hf_params)
hf = hf_params.hyperfine_params

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
	X = (-1j*sx*np.pi).expm();   
	Y = (-1j*sy*np.pi).expm();   
	Z = (-1j*sz*np.pi).expm();   
	x = (-1j*sx*np.pi/2).expm(); mx = (1j*sx*np.pi/2).expm()
	y = (-1j*sy*np.pi/2).expm(); my = (1j*sy*np.pi/2).expm()
	z = (-1j*sz*np.pi/2).expm(); mz = (1j*sz*np.pi/2).expm()
	return X,Y,Z,x,y,z,mx,my,mz

def spin_theta_rotation(theta, phi):
	return (-1j*(np.cos(theta)*sx + np.sin(theta)*sy)*phi).expm()
	
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
X,Y,Z,x,y,z,mx,my,mz = basic_spin_rotations()  # Basic gates
ket0, bra0, ket1, bra1, rho0, rho1, rhom, ketx,brax,ketmx,bramx,rhox,rhomx,kety,bray,ketmy,bramy,rhoy,rhomy = basic_spin_states() # Basic states

gamma_c = 1.0705e3

def dyn_dec_signal(carbon_params,tau, N):
	''' Useful for quick investigation of fingerprints etc.
	'''
	
	if np.size(tau)!=1:
		M=np.zeros([len(carbon_params),np.size(tau)])
	else:
		M = np.zeros([len(carbon_params),np.size(N)])
	for i,carbon_param in enumerate(carbon_params):
		omega_larmor = carbon_param[0]
		HF_par = carbon_param[1]
		HF_perp = carbon_param[2]

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

		self.decouple_scheme = 'XY4'

		self.B_field = kw.pop('B_field',414.1871869) # Not currently used, could be in theory useful
		self.gamma_c = 1.0705e3 #g-factor for C13 in Hz/G

		self.espin_trans = '+1' # Changes sign of A_par and A_perp

		if kw.pop('use_msmt_params',False):
			raise Exception('Not written yet!')
		else:
			carbon_params = kw.pop('carbon_params',False)
			if isinstance(carbon_params,list):
				self.num_carbons = len(carbon_params) 
				self.carbon_params = carbon_params

			elif kw.pop('use_hf_library',False):
				self.carbon_params = []
				for C_key,C_val in sorted(hf.items()): # Sort to make sure alphabetical i.e. C1, C2 etc.
					self.carbon_params.append([2 * np.pi * self.B_field * self.gamma_c,2 * np.pi * C_val['par'],2 * np.pi * C_val['perp']])
				self.num_carbons = len(self.carbon_params)
				self.calc_c_prec_freqs()

			else:
				warnings.warn('No carbon params passed, using dummy params!')
				self.num_carbons = 1
				self.carbon_params = np.array([[2 * np.pi * 443180.29, 2 * np.pi * 35.0e3 , -2 * np.pi * 33.0e3]])
				self.calc_c_prec_freqs()

			self.recalc_Hamiltonian = True

		self.define_useful_states_and_e_operators()
		self.define_gates()

	def e_op(self,operator):
		''' Helper function to embed e operator with identities on carbons '''
		return qutip.tensor([operator] + [Id] * self.num_carbons)

	def e_C_op(self,e_operator,c_operator,c_num):
		''' Helper function to embed e operator with identities on carbons '''
		return qutip.tensor([e_operator] + [Id] * (c_num-1) + [c_operator] + [Id] * (self.num_carbons - c_num))

	def define_useful_states_and_e_operators(self):
		''' Define commonly used electronic operations '''
		self._Xe = self.e_op(X)
		self._Ye = self.e_op(Y)
		self._xe = self.e_op(x)
		self._ye = self.e_op(y)
		self._mxe = self.e_op(mx)
		self._mye = self.e_op(my)
		self._Ide = self.e_op(Id)

		''' standard init state for system '''
		self.initial_state = qutip.tensor([rho0] + [rhom] * self.num_carbons)
	
	def define_gates(self):
		''' This is written this way so that can be easily overwritten for more complex behaviour (say using Hamiltonians instead of unitaries)'''
		self.Xe = lambda : self._Xe
		self.Ye = lambda : self._Ye
		self.xe = lambda : self._xe
		self.ye = lambda : self._ye
		self.mxe = lambda : self._mxe
		self.mye = lambda : self._mye
		self.re = lambda (theta, phi) : self.e_op(spin_theta_rotation(theta, phi))

		self.proj0 = lambda : self.e_op(rho0)
		self.proj1 = lambda : self.e_op(rho1)


	def compile(self,seq,reps = 1):
		''' Again, written this way for flexibility '''
		gate = self._Ide

		for elem in seq:
			gate = elem*gate

		gate = gate**reps

		return gate

	def evolve(self,system,operator, norm = False):
		sysout = operator * system * operator.dag()

		return sysout.unit() if norm else sysout 

	def measure_e(self,sys_state,e_state = 0):
		if e_state == 0:
			e_state = self.proj0()
		elif e_state == 1:
			e_state = self.proj1()

		return np.real((e_state*sys_state).tr())

	def measure_c(self,sys_state,c_num=1,c_state = 0):
		''' Not directly accessible, but sometimes useful'''
		if c_state == 0:
			c_state = rho0
		elif c_state == 1:
			c_state = rh1

		proj = qutip.tensor([Id]*c_num + [c_state] + [Id]*(self.num_carbons - c_num))

		return np.real((proj*sys_state).tr())

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


	def NV_carbon_system_Hamiltonian(self):
		''' Function to calculate the NV C13 system Hamiltonian '''

		if self.recalc_Hamiltonian == True:
			sign = -1 if self.espin_trans == '-1' else 1

			self.Hsys = np.sum([self.e_C_op(rho0,carbon_param[0]*sz,i+1) \
						 + self.e_C_op(rho1,((carbon_param[0]+ sign*carbon_param[1])*sz + sign * carbon_param[2] * sx),i+1) \
				   for i,carbon_param in enumerate(self.carbon_params)]).tidyup()

			self.recalc_Hamiltonian = False

		return self.Hsys
		

	def NV_carbon_ev(self,tau):
		''' Function to calculate a C13 evolution matrix from the system Hamiltonian. Written this way so that could be overwritten'''

		return (-1j*self.NV_carbon_system_Hamiltonian()*tau).expm()


	def nuclear_gate(self,N,tau,**kw):
		'''Evolution during a decouple unit'''
		
		evNV_C_tau = self.NV_carbon_ev(self.nuclear_gate_tau(tau))
		evNV_C_2tau = self.NV_carbon_ev(self.nuclear_gate_tau(2*tau, double_sided = True))

		scheme = kw.pop('scheme',self.decouple_scheme)

		if scheme == 'XY4':
			if N%4 != 0:
				raise Exception('Incompatible number of pulses!')
			gate = self.compile([evNV_C_tau,self.Xe(),evNV_C_2tau,self.Ye(),evNV_C_2tau,self.Xe(),evNV_C_2tau,self.Ye(),evNV_C_tau],reps = N/4)

		elif scheme == 'simple': # CHECK THIS
			gate = self.compile([evNV_C_tau,self.Xe(),evNV_C_tau],reps = N)
		else:
			raise Exception('Unknown scheme!')

		return gate
	

	def nuclear_gate_tau(self,tau,double_sided = False):
		'''Helper function to get tau for gate seq'''
		tau_correction_factor = self.tau_correction_factor if hasattr(self,'tau_correction_factor') else 0

		scale_fact = 0.5 if not(double_sided) else 1.0

		if tau_correction_factor > scale_fact*tau:
			raise Exception('mw_duration too long!')
		return (tau - scale_fact*tau_correction_factor)

	
	def nuclear_phase_gate(self,carbon_nr, phase, state = 'sup'):
		''' For small waits (no decoupling) '''
		if state == 'sup':
			precession_freq = self.c_prec_freqs[carbon_nr-1,2]
		elif state == 0:
			precession_freq = self.c_prec_freqs[carbon_nr-1,0]
		elif state == 1:
			precession_freq = self.c_prec_freqs[carbon_nr-1,1]
		# NEED TO CHECK HOW PHASE IS DEFINED IN OUR EXPM
		phase = (np.pi*float(-1*phase)/180)%(2*np.pi)
		dec_time = phase/precession_freq
		return self.compile([self.NV_carbon_ev(dec_time)])



class noisy_NV_system(NV_system):

	def __init__(self,**kw):

		self.amp_val = kw.pop('amp_val',1.0)
		self.mw_duration = kw.pop('mw_duration',10.0e-9)
		self.tau_correction_factor = self.mw_duration

		NV_system.__init__(self,**kw)


	def set_noisy_params(self,**kw):

		amp_val = kw.pop('amp_val',None)
		if amp_val != None:
			self.amp_val = amp_val

		mw_duration = kw.pop('mw_duration',None)
		if mw_duration != None:
			self.mw_duration = mw_duration
			self.tau_correction_factor = self.mw_duration
		
		self.define_noisy_operators()
		
	
	def finite_microwave_pulse(self,duration,theta,phi):

		Hsys = duration*self.NV_carbon_system_Hamiltonian()
		Hint = self.e_op(phi*(np.cos(theta)*sx + np.sin(theta)*sy))

		return (-1j*(Hsys+Hint)).expm()

	def define_useful_states_and_e_operators(self):
		''' Override commonly used electronic gates '''
		self._Xe = self.finite_microwave_pulse(self.mw_duration,0.0,np.pi * self.amp_val)
		self._Ye = self.finite_microwave_pulse(self.mw_duration,0.5*np.pi,np.pi * self.amp_val)
		self._xe = self.finite_microwave_pulse(self.mw_duration,0.0,0.5*np.pi * self.amp_val)
		self._ye = self.finite_microwave_pulse(self.mw_duration,0.5*np.pi,0.5*np.pi * self.amp_val)
		self._mxe = self.finite_microwave_pulse(self.mw_duration,0.0,-0.5*np.pi * self.amp_val)
		self._mye = self.finite_microwave_pulse(self.mw_duration,0.5*np.pi,-0.5*np.pi * self.amp_val)
		self.re = lambda (theta, phi) : self.finite_microwave_pulse(self.mw_duration,theta,phi * self.amp_val)
		self._Ide = self.e_op(Id)

		''' standard init state for system '''
		self.initial_state = qutip.tensor([rho0] + [rhom] * self.num_carbons)




class NV_experiment(NV_system):

	'''Eventually this should contain code that will interpret a gate seq and make an experiment. '''

			# current_phase = total_time*precession_freq
		# phase_dif = (phase-current_phase)%(2*np.pi)
		# dec_time =  (phase_dif)/precession_freq

		# tau = dec_time/4




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
			
			for j, carbon_param in enumerate(carbon_params):

				NV_system.carbon_params = [carbon_param]
				NV_system.num_carbons = 1
				NV_system.c_prec_freqs = [c_prec_freqs[j]]
				NV_system.recalc_Hamiltonian = True
				NV_system.define_useful_states_and_e_operators()
				NV_system.define_gates()
				

				for i,tau in enumerate(tau_range):

					gate_seq = NV_system.compile([NV_system.xe(),NV_system.nuclear_gate(N ,tau),NV_system.mxe()]) 
					exp0[i,j] = NV_system.measure_e(NV_system.evolve(NV_system.initial_state,gate_seq))

		else:

			exp0 = np.zeros((np.shape(tau_range)[0]))

			for i,tau in enumerate(tau_range):

				gate_seq = NV_system.compile([NV_system.xe(),NV_system.nuclear_gate(N ,tau),NV_system.mxe()]) 
				exp0[i] = NV_system.measure_e(NV_system.evolve(NV_system.initial_state,gate_seq))

	else:
		exp0 = 0.5*(1+dyn_dec_signal(NV_system.carbon_params, tau_range, N)).T
	plt.figure()
	plt.plot(tau_range*1e6,exp0)
	plt.title('Signal'); plt.xlabel('Tau')
	plt.ylim([-0.1,1.1])
	plt.show()
	plt.close()


def prepare_X_and_measure_XY(NV_system,N = 32, tau_range =  np.arange(1e-6,7e-6,1e-7),meas = 'eXY',**kw):
	''' Prepare carbon in X (or attempt to) and measure in X and Y '''
	X = np.zeros(np.shape(tau_range))
	Y = np.zeros(np.shape(tau_range))
	
	for i,tau in enumerate(tau_range):
		mbi_seq = NV_system.compile([NV_system.ye(),NV_system.nuclear_gate(N ,tau),NV_system.mxe()])
		init_state = NV_system.evolve(NV_system.initial_state,NV_system.compile([mbi_seq,NV_system.proj0()]), norm = True)
			
		if meas == 'eXY':
			X[i] = NV_system.measure_e(NV_system.evolve(init_state,mbi_seq))
			Y[i] = NV_system.measure_e(NV_system.evolve(init_state,NV_system.compile([NV_system.nuclear_phase_gate(1,90,state=0),mbi_seq])))
		elif meas == 'nXY': 
			c_num = kw.pop('c_num',1)
			X[i] = NV_system.measure_c(init_state,c_state = rhox,c_num=c_num)
			Y[i] = NV_system.measure_c(init_state,c_state = rhoy,c_num=c_num)

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
	print 'Max fid. ', Fid[ind], ' at ', tau_range[ind]*1e6


def MonteCarlo_MWFid(NV_system,N = 11, tau = 6.582e-6,N_rand = 100,sigma=0.1):
	'''Simulate doing microwave pulses with a certain standard deviation on the pulse amplitude from trial to trial '''
	NV_system
	rands = np.random.normal(loc = 1.0,scale=sigma, size=N_rand)
	fids = np.zeros(N_rand)

	for i, rand_amp in enumerate(rands):

		NV_system.set_noisy_params(amp_val = rand_amp)
	
		fids[i] = 1- NV_system.measure_e(NV_system.evolve(NV_system.initial_state,NV_system.nuclear_gate(N ,tau,scheme = 'simple')))

	print "Fidelity is %f \pm %f" % (np.mean(fids), np.std(fids))

	return fids

def MonteCarlo_MWAmp_CGate_fid(NV_system,N = 32, tau = 6.582e-6,N_rand = 100,sigma=0.1,meas = 'eXY'):
	'''Simulate doing a carbon gate with finite microwave durations and a certain standard deviation on the pulse amplitude from trial to trial '''

	rands = np.random.normal(loc = 1.0,scale=sigma, size=N_rand)
	fids = np.zeros(N_rand)

	for i, rand_amp in enumerate(rands):

		NV_system.set_noisy_params(amp_val = rand_amp)

		mbi_seq = NV_system.compile([NV_system.ye(),NV_system.nuclear_gate(N ,tau),NV_system.mxe(),NV_system.proj0()])
		init_state = NV_system.evolve(NV_system.initial_state,mbi_seq, norm = True)
		if meas == 'eXY':
			X = NV_system.measure_e(NV_system.evolve(init_state,mbi_seq))
			Y = NV_system.measure_e(NV_system.evolve(init_state,NV_system.compile([NV_system.nuclear_phase_gate(1,90,state=0),mbi_seq])))
		elif meas == 'nXY': 
			X = NV_system.measure_c(init_state,c_state = rhox)
			Y = NV_system.measure_c(init_state,c_state = rhoy)

		fids[i] = (np.sqrt((X-0.5)**2 + (Y-0.5)**2)+0.5)

	print "Fidelity is %f \pm %f" % (np.mean(fids), np.std(fids))

	return fids
			
	
def prepare_X_measure_theta(NV_system,N = 32, tau = 6.582e-6,thetas = np.arange(0,360,4)):

	meas = np.zeros(np.shape(thetas))
	
	mbi_seq = NV_system.compile([NV_system.ye(),NV_system.nuclear_gate(N ,tau),NV_system.mxe(),NV_system.proj0()])
	init_state = NV_system.evolve(NV_system.initial_state,mbi_seq, norm = True)
		
	for i,theta in enumerate(thetas):
		meas[i] = NV_system.measure_e(NV_system.evolve(init_state,NV_system.compile([NV_system.nuclear_phase_gate(1,theta,state=0),mbi_seq])))
	
	plt.figure()
	plt.plot(thetas,meas)
	plt.title('Signal'); plt.xlabel('Theta')
	plt.ylim([0,1.1])
	plt.show()
	plt.close()

class composite_gate_NV_system(NV_system):

	def __init__(self,**kw):

		self.composite_scheme = kw.pop('composite_scheme','3_pulse')
		self.r = kw.pop('r',0.31) # For 3 pulse scheme
		self.p = kw.pop('p',0.209) # For 5 pulse scheme
		self.q = kw.pop('q',0.384) # For 5 pulse scheme
		
		NV_system.__init__(self,**kw)

	def nuclear_gate(self,N,tau,**kw):
		'''Evolution during a decouple unit'''

		
		if self.composite_scheme == '3_pulse':
			
			if N%3 != 0:
				raise Exception('Incompatible number of pulses!')
			
			evNV_C_r = self.NV_carbon_ev(self.nuclear_gate_tau(6*tau*(0.5-self.r)))
			evNV_C_t = self.NV_carbon_ev(self.nuclear_gate_tau(6*tau*self.r,double_sided = True))

			gate = self.compile([evNV_C_r,self.Xe(),evNV_C_t,self.Xe(),evNV_C_t,self.Xe(),evNV_C_r],reps = N/3)
			return gate

		elif self.composite_scheme == '5_pulse':

			if N%5 != 0:
				raise Exception('Incompatible number of pulses!')
			
			evNV_C_q = self.NV_carbon_ev(self.nuclear_gate_tau(10*tau*(0.5-self.q)))
			evNV_C_qp = self.NV_carbon_ev(self.nuclear_gate_tau(10*tau*(self.q-self.p),double_sided = True))
			evNV_C_p = self.NV_carbon_ev(self.nuclear_gate_tau(10*tau*self.p,double_sided = True))

			gate = self.compile([evNV_C_q,self.Xe(),evNV_C_qp,self.Xe(),evNV_C_p,self.Xe(),evNV_C_p,self.Xe(),evNV_C_qp,self.Xe(),evNV_C_q],reps=N/5)
			return gate

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