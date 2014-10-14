
import numpy as np
from analysis.lib.fitting import fit, ramsey, common
from analysis.lib.tools import plot
import random
from matplotlib import rc, cm
import os, sys
import h5py
import logging, time

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot
from analysis.lib.tools import compare_functions as compare
from matplotlib import rc, cm

reload(sequence)
reload(compare)

	
class RamseySequence():

	def __init__ (self, N_msmnts, reps, tau0):
		self.N = N_msmnts
		self.points = 2**(self.N+1)+3
		self.discr_steps = 2*self.points+1
		self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
		self.msmnt_results = None
		self.msmnt_phases = None
		self.msmnt_times = None
		self.table_elements = None #Only used by table-based protocols
		self.set_detuning = None
		self.T2 = 96e-6
		self.fid0 = 0.9
		self.fid1 = 0.015
		self.theta = 0*np.pi/180.
		self.t0 = tau0
		self.B_max = 1./(2*tau0)
		self.n_points = 2**(self.N+3)
		self.curr_rep = 0
		self.curr_msmnt = 1
		self.reps = reps
		self.N_total = self.N
		self.M=1
		self.p_k[self.points] = 1/(2.*np.pi)
		self.verbose = True
		self.use_ROfid_in_update = False
		self.renorm_ssro = False
		#parameters for majority vote
		self.maj_reps = None
		self.maj_thr = None
		#parameters for variable-M protocols
		self.G = 0
		self.K = 0
		self.F = 0
		
	def set_ideal (self):
		self.T2 = 1000.
		self.fid0 = 1.
		self.fid1 = 0.
	
	def field_values_MHz (self):
		return 1e-6*np.linspace (-self.B_max, self.B_max, self.n_points)
		
	def reset (self):
		self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
		self.p_k[self.points] = 1/(2.*np.pi)
			
	def reset_rep_counter(self):
		self.curr_rep = 0
		self.curr_msmnt = 1
		
	def plot_p_k (self):	
		x = np.arange(self.discr_steps)-self.points
		plt.plot (x, np.real (self.p_k))
		plt.plot (x, np.imag (self.p_k))
		plt.xlabel ('k')
		plt.show()

	def bayesian_update(self, m_n, phase_n, t_n,repetition = None):
			
		if (repetition == None):
			repetition = self.curr_rep

		p_old = np.copy(self.p_k)
		
		if self.use_ROfid_in_update:
			p0 = (0.5*self.fid0 + m_n-m_n*self.fid0)*p_old 
			p1 = 0.25*self.fid0*(np.exp(1j*(m_n*np.pi+phase_n))*np.roll(p_old, shift = -t_n)) 
			p2 = 0.25*self.fid0*(np.exp(-1j*(m_n*np.pi+phase_n))*np.roll(p_old, shift = +t_n)) 		
		else:
			p0 = 0.5*p_old 
			p1 = 0.25*(np.exp(1j*(m_n*np.pi+phase_n))*np.roll(p_old, shift = -t_n)) 
			p2 = 0.25*(np.exp(-1j*(m_n*np.pi+phase_n))*np.roll(p_old, shift = +t_n)) 
		p = p0+p1+p2
		p = (p/np.sum(np.abs(p)**2)**0.5)
		self.p_k = np.copy (p)
			
	def inc_rep (self):
		self.curr_rep = self.curr_rep + 1
		if (self.curr_rep>self.reps):
			print 'Maximum repetition reached... Setting counter back to zero!'
			self.curr_rep = 0
			
	def optimal_phase (self):
		ttt = -2**(self.N-self.curr_msmnt+1)
		optimal_phase = 0.5*np.angle (self.p_k[ttt+self.points])
		return optimal_phase
		
	def holevo_variance_pk (self):
		return (2*np.pi*self.p_k[-1+self.points])**(-2)-1.

	def print_msmnt(self):
		print 'Measurement results: ', self.msmnt_results
		print 'Phases: ', self.msmnt_phases
		print 'Ramsey time: ', self.msmnt_times*1e9, 'ns'
		
	def convert_to_dict (self):
		
		self.msmnt_dict = {}
		self.msmnt_multiplicity = {}
		self.phases_dict = {}
		self.table_el_dict = {}
		
		ind = 0		
		for j in np.arange(self.reps):
			curr_msmnt = self.msmnt_results [j, :]
			curr_phases = self.msmnt_phases [j, :]
			
			found = 0
			for k in self.msmnt_dict:
				if (np.sum(np.abs(self.msmnt_dict[k]-curr_msmnt))<1e-5):
					self.msmnt_multiplicity[k] = self.msmnt_multiplicity[k]+1
					found = 1
					
			if (found == 0):
				self.msmnt_dict[ind] = curr_msmnt
				self.phases_dict[ind] = curr_phases
				self.msmnt_multiplicity[ind] = 1
				if not(self.table_elements == None):
					self.table_el_dict [ind] = self.table_elements[ind]
				ind = ind+1

	def analysis (self, corrected=True, N_max = None, repetition = None):
		#assumes in each adptv step, we store the number of ones, out of M msmnts (M>1, potentially)
		if (N_max==None):
			N_max = self.N
								
		if (repetition == None):
			repetition = self.curr_rep
			
		prob = np.ones(self.n_points)
		beta = np.linspace (-self.B_max, self.B_max, self.n_points)

		print 'Analysis msmnt results:'
		print self.msmnt_results[repetition, :]

		for n in np.arange(N_max) +(self.N-N_max):
			q = 2*np.pi*beta*self.msmnt_times[n]*self.t0+self.msmnt_phases[repetition, n]
			dec = np.exp(-(self.msmnt_times[n]*self.t0/self.T2)**2)
			nr_ones = self.msmnt_results[repetition, n]
			nr_zeros = self.M-nr_ones
			
			for j in np.arange(nr_ones):
				prob = prob*(1-dec*np.cos(q))
			for j in np.arange(nr_zeros):
				prob = prob*(1+dec*np.cos(q))

			#for j in np.arange(nr_ones):
			#	prob = prob*self.fid0*(1-dec*np.cos(q))
			#for j in np.arange(nr_zeros):
			#	prob = prob*((1-0.5*self.fid0)+0.5*self.fid0*dec*np.cos(q))

			
		prob = prob/np.sum(np.abs(prob))
		return beta, prob


	def analysis_dict (self, phase=[], msmnt_results=[], times=[]):
		#assumes in each adptv step, we store the number of ones, out of M msmnts (M>1, potentially)
			
		prob = np.ones(self.n_points)
		beta = np.linspace (-self.B_max, self.B_max, self.n_points)
		N_max = len(msmnt_results)
		for n in np.arange(N_max) +(self.N-N_max):
			q = 2*np.pi*beta*times[n]*self.t0+phase[n]
			dec = np.exp(-(times[n]*self.t0/self.T2)**2)
			nr_ones = msmnt_results[n]
			nr_zeros = self.M-nr_ones

			for j in np.arange(nr_ones):
				prob = prob*(1-dec*np.cos(q))
			for j in np.arange(nr_zeros):
				prob = prob*(1+dec*np.cos(q))
		prob = prob/np.sum(np.abs(prob))
		return beta, prob


	def analyse_ramsey (self):

		if not(np.all(self.msmnt_phases)):
			print 'Attentions! Phases are not all equal'

		total_reps = 0
		for k in self.msmnt_dict:
			curr_phase = np.rint(self.phases_dict[k])
			curr_msmnt = np.squeeze(np.array(np.rint(self.msmnt_dict[k])))
			mult = np.rint(self.msmnt_multiplicity[k])

			if (total_reps==0):
				msmnt_counts = mult*curr_msmnt
			else:
				msmnt_counts = msmnt_counts + mult*curr_msmnt
			total_reps = total_reps + mult

		ramsey = msmnt_counts/float(total_reps*self.M)
		plt.plot (self.msmnt_times*self.t0*1e9, ramsey, 'ob')
		plt.ylim([0,1])
		plt.title ('Ramsey')
		plt.xlabel ('[ns]')
		plt.show()




	def plot_phase_distribution (self, repetition = 0):
		beta, prob = self.analysis (corrected=False, repetition = repetition)
		plt.plot (beta*1e-6, prob)
		plt.xlim([-20,20])
		plt.show()
		
	def mean_square_error (self, set_value = None, do_plot=False, save_plot=False, xlim = None):
		
		msqe = 0
		total_reps = 0

		if (set_value==None):
			set_value = self.set_detuning
		
		for k in self.msmnt_dict:
			curr_phase = np.rint(self.phases_dict[k])
			curr_msmnt = np.rint(self.msmnt_dict[k])
			mult = np.rint(self.msmnt_multiplicity[k])
			beta, prob = self.analysis_dict (phase = curr_phase, msmnt_results = curr_msmnt, times = np.rint(self.msmnt_times))

			fase = np.exp(1j*2*np.pi*beta*self.t0)
			phi_m = np.sum(fase*prob)

			if (total_reps==0):
				avg_prob = mult*prob
			else:
				avg_prob = avg_prob + mult*prob
			
			phi_set = np.exp(1j*2*np.pi*set_value*self.t0)
			msqe = msqe + mult*(phi_m/phi_set)
			total_reps = total_reps + mult
		avg_prob = avg_prob/np.sum(avg_prob)
		msqe = msqe/float(total_reps)
		msqe = np.abs(msqe)**(-2)-1
		msqe_fB = msqe/((2*np.pi*self.t0)**2)
		sigma_fB = 1e-6*msqe_fB**0.5
		fase = np.exp(1j*2*np.pi*beta*self.t0)
		phi_m = np.sum(fase*avg_prob)
		mean_fB = 1e-6*np.angle(phi_m)/(2*np.pi*self.t0)


		if do_plot:
			f1 = plt.figure()
			plt.plot (beta/1e6, avg_prob)
			plt.xlabel ('magnetic field detuning [MHz]')
			plt.ylabel ('prob distrib')
			plt.title('(B = '+str(mean_fB)+' +- '+str(sigma_fB) + ') MHz')
			if not(xlim==None):
				plt.xlim(xlim)

			if save_plot:
				f_name = 'prob_distr.png'
				savepath = os.path.join(self.folder, f_name)
				f1.savefig(savepath)
			plt.show()

		return beta, avg_prob, msqe, mean_fB, sigma_fB

	def compare_to_simulations(self):

		beta_exp, p_exp, err_exp, mB, sB = self.mean_square_error(set_value=self.set_detuning, do_plot=False)

		s = RamseySequence_Simulation (N_msmnts = self.N, reps=self.reps, tau0=self.t0)
		s.setup_simulation (magnetic_field_hz = self.set_detuning, M=self.M, lab_pc = True)
		s.T2 = self.T2
		s.fid0 = self.fid0
		s.fid1 = self.fid1
		s.maj_reps = self.maj_reps
		s.maj_thr = self.maj_thr

		s.table_based_simulation()
		s.convert_to_dict()
		s.print_table_positions()		
		beta_sim, p_sim, err_sim, a, b = s.mean_square_error(set_value=self.set_detuning, do_plot=False)

		plt.plot (beta_exp*1e-6, p_exp, 'b', label = 'exp')
		plt.plot (beta_sim*1e-6, p_sim, 'r', label = 'sim')\
		#plt.title('(B_exp = '+'{0:.2f}'.format(mB)+' +- '+str(sB) + ') MHz')
		plt.xlabel ('magnetic field detuning [MHz]')
		plt.ylabel ('probability distribution')
		plt.legend()
		plt.show()


	def print_results (self):

		print '---    Measurement results  - SUMMARY    ---'
		print '--------------------------------------------'
		
		for k in self.msmnt_dict:
			curr_phase = self.phases_dict[k]
			curr_msmnt = self.msmnt_dict[k]
			mult = self.msmnt_multiplicity[k]

			print '(*)', curr_msmnt, ' - ', mult, ' times --- phases: ', np.round(curr_phase*180/(np.pi)), ' deg'
		print '--------------------------------------------'


	def print_table_positions (self):

		print '---    Measurement results  - SUMMARY    ---'
		print '--------------------------------------------'
		np.set_printoptions(suppress=True)
		for k in self.msmnt_dict:
			curr_phase = self.phases_dict[k]
			curr_msmnt = self.msmnt_dict[k]
			mult = self.msmnt_multiplicity[k]
			print curr_msmnt, ' - X', mult, '- phases: ', np.round(curr_phase*180/(np.pi)), ' - pos: ', np.rint(self.table_el_dict[k])
		print '--------------------------------------------'


	def hist_phases (self):

		for i in np.arange (self.N_total):
			a = self.msmnt_phases [:,i]*180/np.pi
			plt.hist (a, 30)
			plt.title ('msmnt step nr. ' +str(i+1))
			plt.xlabel ('adaptive phase [deg]')
			plt.show()

	def print_phases(self):
		for i in np.arange (self.N_total):
			print '------------------------------------------'
			print ' Msmnt nr. '+str(i+1)
			print np.unique(self.msmnt_phases [:,i]*180/np.pi)


class RamseySequence_Simulation (RamseySequence):

	def setup_simulation(self, magnetic_field_hz = 0., M = 1, lab_pc=True):
		self.beta_sim = magnetic_field_hz
		self.M = M
		if lab_pc:
			self.root_folder = 'D:/measuring/'
		else:
			self.root_folder = '/home/cristian/Work/Research/teamdiamond/'

	def save_folder (self, folder = '/home/cristian/Work/Research/adaptive magnetometry/'):
		self.save_folder = folder
	
	def ramsey (self, t=0., theta=0.):

		n0 = 0.5*(1-np.cos(2*np.pi*self.beta_sim*t+theta))
		n1 = 0.5*(1+np.cos(2*np.pi*self.beta_sim*t+theta))
	
		pp0 = self.fid0*n0+self.fid1
		pp1 = n1+(1-self.fid0*n0)
		p0 = pp0/(pp0+pp1)
		p1 = pp1/(pp0+pp1)
		np.random.seed()
		result = 1-np.random.choice (2, 1, p=[p0, p1])
		return result[0]		
		
	def majority_vote_msmnt (self, theta_n, t_n):
	
		m_total = 0
		for m in np.arange (self.maj_reps):
			res = self.ramsey (theta=theta_n, t = t_n*self.t0)
			m_total = m_total + res
		m_total = round(m_total/self.fid0)
		if (m_total>self.maj_thr):
			m_res = 1
		else:
			m_res = 0	
		return m_res

	def M_ssro_msmnt (self, theta_n, t_n):
		m_total = 0
		for m in np.arange (self.M)+1:
			res = self.ramsey (theta=theta_n, t = t_n*self.t0)
			m_total = m_total + res
			
		n_ones = m_total
		if self.renorm_ssro:
			n_ones = round(m_total/self.fid0)
			if (n_ones>self.M):
				n_ones = self.M
		return n_ones, self.M-n_ones

	def sim_cappellaro (self):
		
		if self.verbose:				
			print '----------------------------------'
			print 'Simulating Cappellaro protocol'
			print '----------------------------------'
			print '- N = '+str(self.N)+ ', M = '+str(self.M)

		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))
		
		nn = np.arange(self.N)+1
		tau = 2**(self.N-nn)
		self.msmnt_times = tau
		self.reset_rep_counter()
		
		for r in np.arange(self.reps):

			msmnt_results = np.zeros (self.N)
			t = np.zeros (self.N)
			phase = np.zeros(self.N)
			self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
			self.p_k [self.points] = 1/(2.*np.pi)
	
			for n in np.arange(self.N)+1:

				t[n-1] = int(2**(self.N-n))
				ttt = -2**(self.N-n+1)
				phase[n-1] = 0.5*np.angle (self.p_k[ttt+self.points])

				n_ones, n_zeros = self.M_ssro_msmnt (theta_n = phase[n-1], t_n = t[n-1])
				
				for jj in np.arange(n_ones):
					self.bayesian_update (m_n = 1, phase_n = phase[n-1], t_n = 2**(self.N-n))

				for jj in np.arange(n_zeros):
					self.bayesian_update (m_n = 0, phase_n = phase[n-1], t_n = 2**(self.N-n))
			
				msmnt_results[n-1] = n_ones
				
			self.msmnt_results [r, :] = np.copy(msmnt_results)
			self.msmnt_phases [r, :] = np.copy(phase)
			self.inc_rep()

	def sim_cappellaro_majority (self):
		
		if self.verbose:				
			print '-----------------------------------------'
			print 'Simulating Cappellaro protocol (maj vote)'
			print '-----------------------------------------'
			print '- N = '+str(self.N)+ ', M = '+str(self.M)

		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))
		
		nn = np.arange(self.N)+1
		tau = 2**(self.N-nn)
		self.msmnt_times = tau
		self.reset_rep_counter()
		
		for r in np.arange(self.reps):

			msmnt_results = np.zeros (self.N)
			t = np.zeros (self.N)
			phase = np.zeros(self.N)
			self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
			self.p_k [self.points] = 1/(2.*np.pi)
	
			for n in np.arange(self.N)+1:

				t[n-1] = int(2**(self.N-n))
				ttt = -2**(self.N-n+1)
				phase[n-1] = 0.5*np.angle (self.p_k[ttt+self.points])

				m_total = 0
				for m in np.arange(self.M):
					m_res = self.majority_vote_msmnt (theta_n = phase[n-1], t_n = t[n-1])				
					self.bayesian_update (m_n = m_res, phase_n = phase[n-1], t_n = 2**(self.N-n))
					m_total = m_total + m_res
				msmnt_results[n-1] =m_total
				
			self.msmnt_results [r, :] = np.copy(msmnt_results)
			self.msmnt_phases [r, :] = np.copy(phase)
			self.inc_rep()

	def sim_cappellaro_variable_M (self):
		
		if self.verbose:				
			print '-------------------------------------------'
			print 'Simulating Cappellaro protocol (variable M)'
			print '-------------------------------------------'
			print '- N = '+str(self.N)+ ', M = '+str(self.M)

		if (self.G+self.F+self.K==0):
			print 'Simulation parameters G, K, F not set!!'
		else:			
			self.msmnt_phases = np.zeros((self.reps,self.N))
			self.msmnt_times = np.zeros(self.N)
			self.msmnt_results = np.zeros((self.reps,self.N))
		
			nn = np.arange(self.N)+1
			tau = 2**(self.N-nn)
			self.msmnt_times = tau
			self.reset_rep_counter()
			
			for r in np.arange(self.reps):

				msmnt_results = np.zeros (self.N)
				t = np.zeros (self.N)
				phase = np.zeros(self.N)
				self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
				self.p_k [self.points] = 1/(2.*np.pi)
		
				for n in np.arange(self.N)+1:

					t[n-1] = int(2**(self.N-n))
					ttt = -2**(self.N-n+1)
					phase[n-1] = 0.5*np.angle (self.p_k[ttt+self.points])

					m_total = 0
					for m in np.arange(self.M):
						m_res = self.majority_vote_msmnt (theta_n = phase[n-1], t_n = t[n-1])				
						self.bayesian_update (m_n = m_res, phase_n = phase[n-1], t_n = 2**(self.N-n))
						m_total = m_total + m_res
					msmnt_results[n-1] =m_total
					
				self.msmnt_results [r, :] = np.copy(msmnt_results)
				self.msmnt_phases [r, :] = np.copy(phase)
				self.inc_rep()


	def sim_nonadaptive_variable_M (self):
		
		if self.verbose:				
			print '---------------------------------------------'
			print 'Simulating non-adaptive protocol (variable M)'
			print '---------------------------------------------'

		if (self.G+self.F+self.K==0):
			print 'Simulation parameters G, K, F not set!!'
		else:			
			self.total_nr_msmnts = self.G*(2**(self.K+1)-1) + self.F*(2**(self.K+1)-2-self.K)
			self.msmnt_phases = np.zeros((self.reps,self.total_nr_msmnts))
			self.msmnt_times = np.zeros(self.total_nr_msmnts)
			self.msmnt_results = np.zeros((self.reps,self.total_nr_msmnts))
		
			k_array = self.K-np.arange(self.K+1)
			tau = 2**(k_array)
			self.msmnt_times = tau
			self.reset_rep_counter()
			
			for r in np.arange(self.reps):

				msmnt_results = np.zeros (self.total_nr_msmnts)
				t = np.zeros (self.total_nr_msmsnts)
				phase = np.zeros(self.total_nr_msmnts)
				n = 0
				for k in k_array:
					t[n] = 2**(k_array)
					if (n>0):
						phase[n] = phase[n-1] + np.pi/2
					else:
						phase[n] = 0

					MK = self.G+self.F*(self.K-k)
					for m in np.arange(MK):
						m_res = self.ramsey (theta=phase[n], t = t[n]*self.t0)
						self.bayesian_update (m_n = m_res, phase_n = phase[n], t_n = t[n])
					msmnt_results[n] = m_res
					
				self.msmnt_results [r, :] = np.copy(msmnt_results)
				self.msmnt_phases [r, :] = np.copy(phase)
				self.inc_rep()



			
	def load_table (self, N, M):
		ttt = int(np.round(self.t0*1e9))
		self.table_folder = self.root_folder+'/measurement/scripts/Magnetometry/adaptive_tables_lt1/tau0='+str(ttt)+'ns/'		
		name = 'adptv_table_cappellaro_N='+str(self.N)+'_M='+str(self.M)+'.npz'
		a = np.load(self.table_folder+name)
		self.table = a['table']


	def table_based_simulation (self):
	
		if (self.maj_reps == None):
			print 'Majority reps value not set!'
			a = raw_input ('Set it now...')
			self.maj_reps = a
			self.maj_thr = 0
		
		if self.verbose:				
			print '--------------------------------------------'
			print 'Simulating Cappellaro protocol (table-based)'
			print '--------------------------------------------'
			print '- N = '+str(self.N)+ ', M = '+str(self.M)

		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))
		self.table_elements = np.zeros((self.reps,self.N))
		
		nn = np.arange(self.N)+1
		tau = 2**(self.N-nn)
		self.msmnt_times = tau
		self.reset_rep_counter()
		self.load_table (N=self.N, M=self.M)
		
		for r in np.arange(self.reps):

			msmnt_results = np.zeros (self.N)
			t = np.zeros (self.N)
			phase = np.zeros(self.N)
			array_el = np.zeros (self.N)
	
			for n in np.arange(self.N)+1:
						
				if (n==1):
					pp0 = 1
				else:
					pp0 = 1+((self.M+1)**(n-1))/self.M
				pp1 = 0
				for j in np.arange(n-1):
					pp1 = pp1 + msmnt_results[j]*((self.M+1)**j)
				pos = pp0 + pp1
				phase[n-1] = self.table[pos-1]*np.pi/180.
				array_el[n-1] = int(pos)
				phase[0]=0

				m_total = 0
				for m in np.arange(self.M):
					m_res = self.majority_vote_msmnt (theta_n = phase[n-1], t_n = 2**(self.N-n))				
					self.bayesian_update (m_n = m_res, phase_n = phase[n-1], t_n = 2**(self.N-n))
					m_total = m_total + m_res
				msmnt_results[n-1] = m_total
				
			self.msmnt_results [r, :] = np.copy(msmnt_results)
			self.msmnt_phases [r, :] = np.copy(phase)
			self.table_elements[r,:] = np.copy(array_el)
			self.inc_rep()


					
	def save (self, name):
		
		np.savez (self.folder+name+'.npz', msmnt_results = self.msmnt_results, msmnt_times = self.msmnt_times, msmnt_phases=self.msmnt_phases,
					T2=self.T2, fid0=self.fid0, fid1=self.fid1, B_sim = self.beta_sim, tau0 =self.t0, N=self.N, reps = self.reps)
		
	def load (self, name):
		a = np.load (self.folder+name+'.npz')
		self.msmnt_results = a['msmnt_results']
		self.msmnt_times = a['msmnt_times']
		self.msmnt_phases = a['msmnt+phases']
		self.T2 = a['T2']
		self.fid0 = a['fid0']
		self.fid1 = a['fid1']
		self.beta_sim = a['B_sim']
		self.t0 = a['tau0']
		self.N = a['N']
		self.reps = a[ 'reps']
		self.points = 2**(self.N)+3
		self.discr_steps = 2*self.points+1
		self.p_k = np.zeros ((reps,self.discr_steps))+1j*np.zeros ((reps, self.discr_steps))
		self.theta = 0*np.pi/180.
		self.B_max = 1./(2*self.t0)
		self.curr_rep = 0
		self.n_points = 50000
		self.curr_msmnt = 1
		
		for i in np.arange(reps):
			self.p_k[i, self.points] = 1/(2.*np.pi)
		
		print 'Data loaded!'		


class RamseySequence_Exp (RamseySequence):

	def __init__ (self, folder = '', sub_string = ''):
		self.folder = folder
		self.sub_string = sub_string
		self.t0 = 1e-9
		self.B_max = 1./(2*self.t0)
		self.table_elements=None
		
	def set_exp_pars(self, T2, fid0, fid1):
		self.T2=T2
		self.fid0 = fid0
		self.fid1 = fid1

	def load_exp_data (self):

		a = sequence.SequenceAnalysis(self.folder)
		a.get_sweep_pts()
		a.get_magnetometry_data(name='adwindata', ssro = False)
		self.msmnt_results = a.clicks
		self.reps, self.N = np.shape (a.clicks)
		self.n_points =  2**(self.N+3)
		self.t0 = a.t0
		self.B_max = 1./(2*self.t0)
		self.msmnt_times = np.zeros(len(a.ramsey_time))
		self.set_detuning = a.set_detuning
		self.maj_reps = a.maj_reps
		self.maj_thr = a.maj_thr
		for j in np.arange(len(a.ramsey_time)):
			self.msmnt_times[j] = a.ramsey_time[j]/self.t0
		self.msmnt_phases = 2*np.pi*a.set_phase/255.
		self.M=a.M
	
		phases_detuning = 2*np.pi*a.phases_detuning/360.
		b = np.ones(self.reps)
		self.msmnt_phases = np.mod(self.msmnt_phases - np.outer (b, phases_detuning), 2*np.pi)



class AdaptiveTable ():

	def __init__(self, N, M):
		self.N = N
		self.M = M 
		self.T2 = 96e-6
		self.t0 = 1e-9
		
		self.table = None
		self.test_pp1 = None
		self.test_pp0 = None
		self.max_element = None
		
		self.verbose = False
		self.save_folder = None

	def set_tau0(self, tau0):
		self.t0 = tau0

	def M_conv (self, element = 1):
	
		msmnt_results = np.zeros(self.N)
		a = element
		for k in np.arange(self.N-1)+1:
			msmnt_results [self.N-k] = np.mod(a,self.M+1)
			a = a/(self.M+1)
		return msmnt_results[1:]

	
	def generate (self):
	
		total_length = int(((self.M+1)**(self.N)-1)/self.M)+1

		self.table = np.zeros(total_length)
		self.test_pp1 = np.zeros(total_length)
		self.test_pp0 = np.zeros(total_length)
		self.max_element = (self.M+1)**(self.N-1)
		print 'Total length array: ', total_length
				
		for j in np.arange (self.max_element):
		
			self.curr_j = j
			msmnt_results = self.M_conv (element=j)
			opt_phase, pos = self.msmnt_to_position (msmnt_results = msmnt_results)

			for i in np.arange(len(opt_phase)):
				self.table[pos[i]] = opt_phase[i]
				
			if self.verbose:
				print '   ----- Analyzing round ' +str(j)
				print '      msmnt results: ', msmnt_results
				print '      optimal phases: ', opt_phase
				print '      stored in positions: pp0 = ', pos
				
				
		#adwin basic counts arrays elements starting from 1:			
		self.table = self.table[1:]



	def msmnt_to_position (self, msmnt_results = []):
	
		seq = RamseySequence_Simulation (N_msmnts = self.N, reps=1, tau0=self.t0)
		seq.M = self.M

		phase = np.zeros(self.N)
		seq.p_k = np.zeros (seq.discr_steps)+1j*np.zeros (seq.discr_steps)
		seq.p_k [seq.points] = 1/(2.*np.pi)

		#n is the msmnt currently being performed
		opt_phase = 0.
		optimal_phases = []
		summary_pp0 = []
		summary_pp1 = []
		
		for n in np.arange(self.N)+1:
			ttt = -2**(self.N-n+1)
			opt_phase = 0.5*np.angle (seq.p_k[ttt+seq.points])
			optimal_phases.append(np.round(opt_phase*180/np.pi))

			pp0 = 1+((((self.M+1)**(n-1))-1)/self.M)
			pp1 = 0
			for i in np.arange(n-1):
				pp1 = pp1+msmnt_results[i]*(self.M+1)**(i)			
			position = pp0+pp1
			summary_pp0.append(pp0)
			summary_pp1.append(pp1)	
			pos_array = n-1	

			if n<self.N:		
				for m in np.arange(msmnt_results[pos_array]):
					seq.bayesian_update (m_n = 1, phase_n = opt_phase, t_n = 2**(self.N-n))
				for m in np.arange(seq.M - msmnt_results[pos_array]):
					seq.bayesian_update (m_n = 0, phase_n = opt_phase, t_n = 2**(self.N-n))
				
		if self.verbose:
			print '      msmnt results: ', msmnt_results
			print '      optimal phases: ', optimal_phases
			print '      stored in positions: pp0 = ', summary_pp0
			print '                           pp1 = ', summary_pp1
			print '                           pos = ', np.array(summary_pp1)+np.array(summary_pp0)
			
		return np.array(optimal_phases), np.array(summary_pp1)+np.array(summary_pp0)
				

	def save_table(self, label=''):
		if (self.save_folder==None):
			print 'Please set save folder!!'
		else:
			name = 'adptv_table_cappellaro_N='+str(self.N)+'_M='+str(self.M)+label+'.npz'
			np.savez (self.save_folder+name, table = self.table, N=self.N, M=self.M, T2=self.T2, t0=self.t0)
		

class AdaptiveMagnetometry ():

	def __init__(self, N, tau0):
		self.N = N 				#here it corresponds to the maximum N
		self.t0 = tau0
		self.reps = None
		self.B_max = 1./(2*self.t0)
		self.n_points = 2**(self.N+3)
		self.nr_B_points = 2**(self.N+2)/2
		self.results_dict = {} 
		self.prob_density_dict = {}
		self.analyzed_N = []

		self.nr_periods = None
		self.nr_points_per_period = None
		
		self.scaling_variance = np.zeros(self.N+1)
		self.total_time = np.zeros(self.N+1)

		self.simulated_data = None

		
	def set_exp_params(self, T2 = 96e-6, fid0 = 0.88, fid1 = 0.015):
		self.T2 = T2
		self.fid0 = fid0
		self.fid1 = fid1
		
	def set_protocol (self, M=1, renorm_ssro = False, maj_reps = None, maj_thr = None):
		self.M=M
		self.theta = 0*np.pi/180.
		self.use_ROfid_in_update = False
		self.renorm_ssro = renorm_ssro
		self.maj_reps = maj_reps
		self.maj_thr = maj_thr

	def set_sweep_params (self, nr_periods = 3, nr_points_per_period=11):
		self.nr_periods = nr_periods
		self.nr_points_per_period = nr_points_per_period

		
	def sweep_field_fixedN_simulation (self, N = None):
	
		if (N==None):
			N = self.N

		self.simulated_data = True			
		self.B_values = np.linspace(0, self.B_max/2., self.nr_B_points)
		ind = 0
		print 'Processing: N = '+str(N), ' (M='+str(self.M)+')'
		for b in self.B_values:		
			s = RamseySequence_Simulation (N_msmnts = N, reps=self.reps, tau0=self.t0)
			s.setup_simulation (magnetic_field_hz = b, M=self.M)
			s.verbose = False
			s.T2 = self.T2
			s.fid0 = self.fid0
			s.fid1 = self.fid1
			s.renorm_ssro = self.renorm_ssro
			s.maj_reps = self.maj_reps
			s.maj_thr = self.maj_thr
			s.table_based_simulation ()
			s.convert_to_dict()
			beta, p, err = s.mean_square_error(set_value=b)
			self.results[N-1, ind] = err**0.5
			ind =ind+1

	def load_sweep_field_data (self, N):

		self.simulated_data = False
		self.analyzed_N.append(N)
		msqe = np.zeros(self.nr_points_per_period*self.nr_periods)
		B_field = np.zeros(self.nr_points_per_period*self.nr_periods)
		ind = 0
		check_params_labels = ['tau0', 'M', 'maj_reps', 'maj_thr']

		for per in np.arange(self.nr_periods):
			for pt in np.arange (self.nr_points_per_period):
				label = 'N='+str(N)+'_M='+str(self.M)+'_majReps='+str(self.maj_reps)+'_majThr='+str(self.maj_thr)+'_p'+str(per)+'b'+str(pt)
				#print label
				f = toolbox.latest_data(contains=label)
				s = RamseySequence_Exp (folder = f)
				s.set_exp_pars (T2=96e-6, fid0=0.876, fid1=1-.964)
				s.load_exp_data()
				check_params = [(s.t0 == self.t0), (s.M == self.M), (s.maj_reps==self.maj_reps), (s.maj_thr==self.maj_thr)]
				if np.all(check_params):
					s.convert_to_dict()
					beta, prob, err, mB, sB = s.mean_square_error(do_plot=False, save_plot=True)
					self.prob_density_dict[label] = prob
					#print s.set_detuning, mB
					msqe [ind] = err
					B_field [ind] = mB
				else:
					msg = []
					for i in np.arange (4):
						if (check_params[i]==False): msg.append(check_params_labels[i])
					print 'Non matching parameters: ', msg, ' --- ', label
				ind +=1
		self.results_dict[str(N)] = {'B_field':B_field, 'msqe':msqe, 'M':self.M, 'maj_reps':self.maj_reps, 'maj_thr':self.maj_thr}

	def plot_msqe_dictionary(self):

		C = compare.compare_functions()
		C.xlabel = 'magnetic field detuning [MHz]'
		C.ylabel = 'mean square error [MHz^2]'
		for n in self.analyzed_N:	
			C.add (x =self.results_dict[str(n)]['B_field']*1e-6, y=self.results_dict[str(n)]['msqe'], label=str(n))
		C.plot()

	def sweep_field (self, do_simulate = True):
	
		for n in np.arange(self.N)+1:
			self.sweep_field_fixedN (N=n, do_simulate = do_simulate)
			
	def plot_msqe (self, N=[]):
		r = compare_functions()	
		r.xlabel = 'field [MHz]'
		r.ylabel = 'mean square error'
		r.title = 'N = '+str(self.N)+', M = '+str(self.M)+', sweep magnetic field'
		for n in []:
			r.add (x=self.B_values[n-1,:]/1e6, y=self.results[n-1, :], legend = str(n))
		r.plot()


	def plot_scaling (self):

		for n in self.analyzed_N:
			msqe = self.results_dict[str(n)]['msqe']
			self.scaling_variance [n] = np.mean(msqe)
			self.total_time [n] = self.M*self.maj_reps*(2**(n+1)-1)
		
		self.plot_scaling()
			
		plt.figure()
		plt.loglog (self.total_time*self.t0/1000., self.scaling_variance*self.total_time, 'b')
		plt.loglog (self.total_time*self.t0/1000., self.scaling_variance*self.total_time, 'ob')
		plt.xlabel ('total ramsey time')
		plt.ylabel ('$V_H*T$')
		plt.show()
		
	def save(self, folder = None):
		
		if (folder==None):
			self.folder = '/home/cristian/Work/Research/adaptive magnetometry/data_analysis/'
		else:
			self.folder = folder

		fName = time.strftime ('%Y%m%d_%H%M%S')+'_adaptive_magnetometry_M='+str(self.M)+'_maj=('+str(self.maj_reps)+', '+str(maj_thr)+').hdf5'
		print fName
		#np.savez(self.folder+fName, B_field=self.B_values, results =self.results, total_time = self.total_time, scaling_variance = self.scaling_variance) 
	
