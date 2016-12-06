''' A module to calculate the 13C nuclear and electron spin dynamics
under dynamical decoupling gates. By THT '''

import numpy as np
import qutip
qutip = reload(qutip)
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

###########################
###  Helper functions   ###
###########################

def evolve(system,operator, norm = False):
	sysout = operator * system * operator.dag()

	return sysout.unit() if norm else sysout 

###########################
### 	 Classes        ###
###########################

class NV_system(object):
	''' Basic class to contain the parameters of the NV system, plus some functions to simulate its evolution '''

	def __init__(self,**kw):

		self.decouple_scheme = 'XY4'

		# self.B_field = kw.pop('B_field',400) # Not currently used, could be in theory useful

		self.espin_trans = '+1' # Changes sign of A_par and A_perp

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
				self.carbon_params = np.array([[2 * np.pi * 443180.29, 2 * np.pi * 35.0e3 , -2 * np.pi * 33.0e3]])
				self.calc_c_prec_freqs()

		self.set_initial_state()
		
		self.define_e_operators()

	def define_e_operators(self):
		''' Define commonly used electronic operations '''
		self.Xe = self.e_op(X)
		self.Ye = self.e_op(Y)
		self.xe = self.e_op(x)
		self.ye = self.e_op(y)
		self.mxe = self.e_op(mx)
		self.mye = self.e_op(my)

		''' define some simple spin rotations in terms of their HAMILTONIANS'''
		self.hXe = self.e_op(sx*np.pi)
		self.hYe = self.e_op(sy*np.pi)
		self.hxe = self.e_op(sx*np.pi/2)
		self.hye = self.e_op(sy*np.pi/2)
		self.hmxe = self.e_op(sx*np.pi/2)
		self.hmye = self.e_op(sy*np.pi/2)

	def set_initial_state(self, NV_state = rho0):
		''' Helper function to give standard init state for system '''
		self.initial_state = qutip.tensor([NV_state] + [rhom] * self.num_carbons)

	def calc_c_prec_freqs(self):
		''' If these arent specifed, set them from the carbon params '''
		sign = -1 if self.espin_trans == '-1' else 1

		freqs = []
		for i in range(self.num_carbons):
			omega0 = self.carbon_params[i,0]
			omega1 = np.sqrt((self.carbon_params[i,0] + sign*self.carbon_params[i,1])**2 + self.carbon_params[i,2]**2)
			omega_sup = 0.5*(omega0+omega1)
			freqs.append([omega0,omega1,omega_sup])

		self.c_prec_freqs = np.array(freqs)

	def e_op(self,operator):
		''' Helper function to embed e operator with identities on carbons '''
		return qutip.tensor([operator] + [Id] * self.num_carbons)

	def NV_carbon_system_Hamiltonian(self):
		''' Function to calculate the NV C13 system Hamiltonian '''

		sign = -1 if self.espin_trans == '-1' else 1

		H0s = [0.5*carbon_param[0]*sz for carbon_param in self.carbon_params] 
		H1s = [0.5*((carbon_param[0]+ sign*carbon_param[1])*sz + sign * carbon_param[2] * sx) for carbon_param in self.carbon_params] 
		
		#Hamiltonians for ms=0 and ms=+/-1
		Hsys = qutip.tensor([rho0]+H0s) + qutip.tensor([rho1]+H1s)
		
		return Hsys

	def NV_carbon_evolution_matrix(self,tau):
		''' Function to calculate a C13 evolution matrix from the system Hamiltonian'''

		return (-1j*self.NV_carbon_system_Hamiltonian()*tau).expm(); 
	

	def nuclear_gate(self,N,tau,**kw):
		'''Evolution during a decouple unit'''
		
		expH_tau = self.NV_carbon_evolution_matrix(nuclear_gate_tau(tau))
		expH_2tau = expH_tau**2

		scheme = kw.pop('scheme',self.decouple_scheme)

		if scheme == 'XY4':
			if N%4 != 0:
				raise Exception('Incompatible number of pulses!')
			expH_dec = expH_tau*self.Ye*expH_2tau*self.Xe*expH_2tau*self.Ye*expH_2tau*self.Xe*expH_tau
			expHGate = expH_dec**(N/4)
		elif scheme == 'simple': # CHECK THIS
			expH_dec = expH_tau*self.Xe*expH_tau
			expHGate = expH_dec**(N)
		else:
			raise Exception('Unknown scheme!')

		return expHGate
		
	def nuclear_gate_tau(self,tau):
		tau_correction_factor = self.tau_correction_factor if hasattr(self,'tau_correction_factor') else 0
		if tau_correction_factor > tau:
			raise Exception('mw_duration too long!')
		return 0.5*(tau - tau_correction_factor)

	def assemble_ev_mat(self,gate_comps,comp_defs):
		'''Helper function for gate assembly'''

		ev_mat = self.e_op(Id) 
		for comp in gate_comps:
			ev_mat = ev_mat*comp_defs[comp]
		return ev_mat

	def assemble_Ham_list(self,gate_comps,comp_defs):
		'''Helper function for gate assembly'''
		Ham_list = []
		for comp in gate_comps:
			Ham_list.append(comp_defs[comp])
		return Ham_list

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
		dec_time = 2*phase/precession_freq
		return self.NV_carbon_evolution_matrix(dec_time)



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
		
		self.define_e_operators()
		
	
	def finite_microwave_pulse(self,duration,theta,phi):

		Hsys = duration*self.NV_carbon_system_Hamiltonian()
		Hint = self.e_op(phi*(np.cos(theta)*sx + np.sin(theta)*sy))

		return (-1j*(Hsys+Hint)).expm()

	def define_e_operators(self):
		''' Override commonly used electronic operations '''
		self.Xe = self.finite_microwave_pulse(self.mw_duration,0.0,np.pi * self.amp_val)
		self.Ye = self.finite_microwave_pulse(self.mw_duration,0.5*np.pi,np.pi * self.amp_val)
		self.xe = self.finite_microwave_pulse(self.mw_duration,0.0,0.5*np.pi * self.amp_val)
		self.ye = self.finite_microwave_pulse(self.mw_duration,0.5*np.pi,0.5*np.pi * self.amp_val)
		self.mxe = self.finite_microwave_pulse(self.mw_duration,0.0,-0.5*np.pi * self.amp_val)
		self.mye = self.finite_microwave_pulse(self.mw_duration,0.5*np.pi,-0.5*np.pi * self.amp_val)

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


def C13_fingerprint(NV_system,N = 32, tau_range =  np.arange(1e-6,7e-6,1e-7)):

		exp0 = np.zeros(np.shape(tau_range))
		
		proj0 = NV_system.e_op(rho0)
	
		for i,tau in enumerate(tau_range):

			gate_seq = NV_system.mxe*NV_system.nuclear_gate(N ,tau)*NV_system.xe
			exp0[i] = np.real((proj0 * evolve(NV_system.initial_state,gate_seq)).tr())

		plt.figure()
		plt.plot(tau_range*1e6,exp0)
		plt.title('Signal'); plt.xlabel('Tau')
		plt.ylim([0,1.1])
		plt.show()
		plt.close()


def prepare_X_and_measure(NV_system,N = 32, tau_range =  np.arange(1e-6,7e-6,1e-7),meas = 'eXY'):

		X = np.zeros(np.shape(tau_range))
		Y = np.zeros(np.shape(tau_range))
		proj0 = NV_system.e_op(rho0)

		for i,tau in enumerate(tau_range):
			mbi_seq = NV_system.mxe*NV_system.nuclear_gate(N ,tau)*NV_system.ye
			init_state = evolve(NV_system.initial_state,proj0*mbi_seq, norm = True)
				
			if meas == 'eXY':
				X[i] = np.real((proj0 * evolve(init_state,mbi_seq)).tr())
				Y[i] = np.real((proj0 * evolve(init_state,mbi_seq*NV_system.nuclear_phase_gate(1,90,state=0))).tr())
			elif meas == 'nXY':
				init_state = init_state.ptrace(1)
				X[i] = np.real((rhox * init_state).tr())
				Y[i] = np.real((rhoy * init_state).tr())

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


def MonteCarlo_MWFid(NV_system,N = 11, tau = 15e-6,N_rand = 100,sigma=0.1):

		NV_system
		rands = np.random.normal(loc = 1.0,scale=sigma, size=N_rand)
		fids = np.zeros(N_rand)

		proj0 = NV_system.e_op(rho0)

		for i, rand_amp in enumerate(rands):

			NV_system.set_noisy_params(amp_val = rand_amp)
		
			fids[i] = 1-np.real((proj0 * evolve(NV_system.initial_state,NV_system.nuclear_gate(N ,tau,scheme = 'simple'))).tr())

		print "Fidelity is %f \pm %f" % (np.mean(fids), np.std(fids))

		return fids

def MonteCarlo_MWAmp_NGate_fid(NV_system,N = 32, tau = 6.55e-6,N_rand = 100,sigma=0.1,meas = 'eXY'):

		rands = np.random.normal(loc = 1.0,scale=sigma, size=N_rand)
		fids = np.zeros(N_rand)

		proj0 = NV_system.e_op(rho0)

		for i, rand_amp in enumerate(rands):

			NV_system.set_noisy_params(amp_val = rand_amp)

			mbi_seq = NV_system.mxe*NV_system.nuclear_gate(N ,tau)*NV_system.ye
			init_state = evolve(NV_system.initial_state,proj0*mbi_seq, norm = True)
			if meas == 'eXY':
				X = np.real((proj0 * evolve(init_state,mbi_seq)).tr())
				Y = np.real((proj0 * evolve(init_state,mbi_seq*NV_system.nuclear_phase_gate(1,90,state=0))).tr())
			elif meas == 'nXY':
				init_state = init_state.ptrace(1)
				X = np.real((rhox * init_state).tr())
				Y = np.real((rhoy * init_state).tr())

			fids[i] = (np.sqrt((X-0.5)**2 + (Y-0.5)**2)+0.5)

		print "Fidelity is %f \pm %f" % (np.mean(fids), np.std(fids))

		return fids
			
	
def prepare_X_measure_theta(NV_system,N = 32, tau = 6.51e-6,thetas = np.arange(0,360,4)):

		meas = np.zeros(np.shape(thetas))
		
		proj0 = NV_system.e_op(rho0)
		
		mbi_seq = NV_system.mxe*NV_system.nuclear_gate(N ,tau)*NV_system.ye
		init_state = evolve(NV_system.initial_state,proj0*mbi_seq, norm = True)
			
		for i,theta in enumerate(thetas):
			meas[i] = np.real((proj0 * evolve(init_state,mbi_seq*NV_system.nuclear_phase_gate(1,theta,state=0))).tr())
		
		plt.figure()
		plt.plot(thetas,meas)
		plt.title('Signal'); plt.xlabel('Theta')
		plt.ylim([0,1.1])
		plt.show()
		plt.close()



# The QUTIP solvers crash on this computer :( 	
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
