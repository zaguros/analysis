
import numpy as np
from analysis.lib.fitting import fit, ramsey, common
from analysis.lib.tools import plot
import random
from matplotlib import rc, cm
import os, sys
import h5py
import logging, time, timeit

import matplotlib.pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import  sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr, common
from analysis.lib.tools import plot
from analysis.lib.tools import compare_functions as compare
from analysis.lib.m2 import m2
from matplotlib import rc, cm
from analysis.lib.math import statistics as stat
from compiler.ast import flatten

reload(sequence)
reload(compare)
reload(toolbox)
reload (stat)
reload (fit)

	
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
		self.theta = 0*np.pi/2#0*np.pi/180.
		self.t0 = tau0
		self.B_max = 1./(2*tau0)
		self.n_points = 2**(self.N+3)
		self.curr_rep = 0
		self.curr_msmnt = 1
		self.reps = reps
		self.N_total = self.N
		self.p_k[self.points] = 1/(2.*np.pi)
		self.verbose = True
		self.use_ROfid_in_update = False
		self.renorm_ssro = False
		self.phase_update = False

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
		plt.figure(figsize=(8,4))
		plt.plot (x-0.05, 2*np.pi*np.real(self.p_k), 'ok', markersize=6, markeredgecolor='k')
		plt.plot (x+0.05, 2*np.pi*np.imag(self.p_k), 'or', markersize=6, markeredgecolor='r')
		for i in np.arange (len(self.p_k)):
			plt.vlines (x=x[i]-0.05, ymin = min(0, 2*np.pi*np.real (self.p_k[i])), ymax =  max(0, 2*np.pi*np.real (self.p_k[i])), linewidth =1, color = 'k')
			plt.vlines (x=x[i]+0.05, ymin = min(0, 2*np.pi*np.imag (self.p_k[i])), ymax = max(0, 2*np.pi*np.imag (self.p_k[i])), linewidth =1, color = 'r')
		plt.xlabel ('k', fontsize = 15)
		rc('xtick', labelsize=15) 
		rc('ytick', labelsize=15)
		plt.xlim([-20,20])
		plt.ylim([-1,1])
		plt.show()

	def bayesian_update(self, m_n, phase_n, t_n,repetition = None):

		if (repetition == None):
			repetition = self.curr_rep

		p_old = np.copy(self.p_k)
		

		p0 = p_old*((1-m_n)-((-1)**m_n)*(self.fid0+self.fid1)/2.) 
		p1 = ((-1)**m_n)*(self.fid0-self.fid1)*0.25*(np.exp(1j*(phase_n))*np.roll(p_old, shift = -t_n)) 
		p2 = ((-1)**m_n)*(self.fid0-self.fid1)*0.25*(np.exp(-1j*(phase_n))*np.roll(p_old, shift = +t_n)) 
		p = p0+p1+p2
		p = (p/np.sum(np.abs(p)**2)**0.5)
		p = p/(2*np.pi*np.real(p[self.points]))
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
		self.index = {}
		
		ind = 0		
		for j in np.arange(self.reps):

			if (self.N>1):
				curr_msmnt = self.msmnt_results [j, :]
				curr_phases = self.msmnt_phases [j, :]
			else:
				curr_msmnt = np.array(self.msmnt_results [j])
				curr_phases = np.array(self.msmnt_phases [j])
			
			found = 0
			for k in self.msmnt_dict:
				if (np.sum(np.abs(self.msmnt_dict[k]-curr_msmnt))<1e-5):
					self.msmnt_multiplicity[k] = self.msmnt_multiplicity[k]+1
					self.index[k].append(j)
					found = 1
					
			if (found == 0):
				self.msmnt_dict[ind] = curr_msmnt
				self.phases_dict[ind] = curr_phases
				self.msmnt_multiplicity[ind] = 1
				self.index[ind]=[j]
				if not(self.table_elements == None):
					self.table_el_dict [ind] = self.table_elements[ind]
				ind = ind+1

	def reshape_SQL (self, G):		
		total_data_size = len(self.stored_msmnt_results)
		eff_reps = round(total_data_size/G)
		eff_data_size = eff_reps*G

		self.msmnt_results = np.copy (np.reshape(self.stored_msmnt_results[:eff_data_size], (eff_reps, G)))
		self.msmnt_phases = np.copy (np.reshape(self.stored_msmnt_phases[:eff_data_size], (eff_reps, G)))
		self.msmnt_times = np.ones(G)*self.stored_msmnt_times
		self.G = G
		self.reps = eff_reps
		print  'Reshaping: G=', self.G, '   reps = ', self.reps 


	def plot_updates (self, repetition = None):
								
		if (repetition == None):
			repetition = 0
			
		prob = np.ones(self.n_points)
		beta = np.linspace (-self.B_max, self.B_max, self.n_points)
		prob = prob/np.sum(prob)

		plt.plot (beta*1e-6, prob, 'b', linewidth = 3)
		plt.ylim([0,max(prob)+0.01])
		plt.xlim([-25, 25])
		plt.show()

		for n in np.arange(self.N):
			q = 2*np.pi*beta*self.msmnt_times[n]*self.t0+self.msmnt_phases[repetition, n]
			dec = np.exp(-(self.msmnt_times[n]*self.t0/self.T2)**2)
			nr_ones = self.msmnt_results[repetition, n]
			MK=self.G+self.F*n
			nr_zeros = MK-nr_ones
			rams_1 = ((self.fid0+self.fid1)-(self.fid0-self.fid1)*dec*np.cos(q))/2.
			rams_0 = (1-((self.fid0+self.fid1)-(self.fid0-self.fid1)*dec*np.cos(q))/2.)

			for j in np.arange(nr_ones):
				prob = prob*rams_1
			for j in np.arange(nr_zeros):
				prob=prob*rams_0

			print 'result: ', nr_ones, ' --- phase: ', self.msmnt_phases[repetition, n]*180/np.pi
			prob = prob/np.sum(np.abs(prob))
			plt.plot (beta*1e-6, prob, 'b', linewidth = 4)
			plt.plot (beta*1e-6, rams_0*max(prob), '--r')
			plt.plot (beta*1e-6, rams_1*max(prob), '--k')

			plt.ylim([0,max(prob)+0.01])
			plt.xlim([-25, 25])
			plt.show()

		return beta, prob




	def analysis_dict (self, phase=[], msmnt_results=[], times=[]):
		#assumes in each adptv step, we store the number of ones, out of M msmnts (M>1, potentially)
			
		prob = np.ones(self.n_points)
		beta = np.linspace (-self.B_max, self.B_max, self.n_points)
		try:
			N_max = len(msmnt_results)
		except:
			N_max=1
			msmnt_results = np.array([msmnt_results])
			times = np.array([times])
			phase = np.array([phase])
		for n in np.arange(N_max) +(self.N-N_max):
			q = 2*np.pi*beta*times[n]*self.t0+phase[n]
			dec = np.exp(-(times[n]*self.t0/self.T2)**2)
			nr_ones = msmnt_results[n]
			MK=self.G+self.F*n
			nr_zeros = MK-nr_ones
			for j in np.arange(nr_ones):
					prob = prob*((self.fid0+self.fid1)-(self.fid0-self.fid1)*dec*np.cos(q))/2.
			for j in np.arange(nr_zeros):
					prob = prob*(1-((self.fid0+self.fid1)-(self.fid0-self.fid1)*dec*np.cos(q))/2.)
		prob = prob/np.sum(np.abs(prob))

		return beta, prob

	def analysis_dict_phase_update (self, phase=[], msmnt_results=[], times=[]):
			
		prob = np.ones(self.n_points)
		beta = np.linspace (-self.B_max, self.B_max, self.n_points)
		try:
			N_max = len(msmnt_results)
		except:
			N_max=1
			msmnt_results = np.array([msmnt_results])
			times = np.array([times])
			phase = np.array([phase])

		for n in np.arange(N_max):
			q = 2*np.pi*beta*times[n]*self.t0+phase[n]
			dec = np.exp(-(times[n]*self.t0/self.T2)**2)
			if (msmnt_results[n]==1):
				prob = prob*((self.fid0+self.fid1)-(self.fid0-self.fid1)*dec*np.cos(q))/2.
			else:
				prob = prob*(1-((self.fid0+self.fid1)-(self.fid0-self.fid1)*dec*np.cos(q))/2.)

		prob = prob/np.sum(np.abs(prob))

		return beta, prob



	def B_vs_index (self):
		B_dict = {}
		for k in self.msmnt_dict:
			curr_phase = self.phases_dict[k]
			curr_msmnt = np.rint(self.msmnt_dict[k])
			beta, prob = self.analysis_dict (phase = curr_phase, msmnt_results = curr_msmnt, times = np.rint(self.msmnt_times))
			fase = np.exp(1j*2*np.pi*beta*self.t0)
			phi_m = np.sum(fase*prob)
			mean_fB = 1e-6*np.angle(phi_m)/(2*np.pi*self.t0)
			B_dict[k] = mean_fB
		return B_dict, self.index
						

	def analyse_ramsey (self):

		if not(np.all(self.msmnt_phases)):
			print 'Attentions! Phases are not all equal'

		total_reps = 0
		for k in self.msmnt_dict:
			curr_phase = self.phases_dict[k]
			curr_msmnt = np.squeeze(np.array(np.rint(self.msmnt_dict[k])))
			mult = np.rint(self.msmnt_multiplicity[k])

			if (total_reps==0):
				msmnt_counts = mult*curr_msmnt
			else:
				msmnt_counts = msmnt_counts + mult*curr_msmnt
			total_reps = total_reps + mult
		# Note: does not work for F>0
		ramsey = msmnt_counts/float(total_reps*self.G)
		plt.plot (self.msmnt_times*self.t0*1e9, ramsey, 'ob')
		plt.ylim([0,1])
		plt.title ('Ramsey')
		plt.xlabel ('[ns]')
		plt.show()

	def mean_square_error (self, set_value = None, do_plot=False, save_plot=False, show_plot = False, xlim = None,y_log=False, return_all_estimated_phases=False):
		
		msqe = 0
		total_reps = 0
		ave_exp=0
		if (set_value==None):
			set_value = self.set_detuning
		
		list_phase_values = []

		for k in self.msmnt_dict:
			
			curr_phase = self.phases_dict[k]
			curr_msmnt = np.rint(self.msmnt_dict[k])
			#print len(curr_phase)
			mult = np.rint(self.msmnt_multiplicity[k])

			if (self.phase_update or self.swarm_opt):
				beta, prob = self.analysis_dict_phase_update (phase = curr_phase, msmnt_results = curr_msmnt, times = np.rint(self.msmnt_times))
			else:
				beta, prob = self.analysis_dict (phase = curr_phase, msmnt_results = curr_msmnt, times = np.rint(self.msmnt_times))
			
			fase = np.exp(1j*2*np.pi*beta*self.t0)
			phi_m = np.sum(fase*prob)
			phi_set = np.exp(1j*2*np.pi*set_value*self.t0)

			if return_all_estimated_phases:
				list_phase_values.append(mult*[phi_m/phi_set])

			if (total_reps==0):
				avg_prob = mult*prob
			else:
				avg_prob = avg_prob + mult*prob

			msqe = msqe + mult*(np.sum(prob*beta)-set_value)**2
			ave_exp = ave_exp + mult*(phi_m/phi_set)
			total_reps = total_reps + mult
		avg_prob = avg_prob/np.sum(avg_prob)
		ave_exp = ave_exp/float(total_reps)
		H = np.abs(ave_exp)**(-2)-1
		sigma_fB = 1e-6*(msqe/float(total_reps))**0.5
		
		fase = np.exp(1j*2*np.pi*beta*self.t0)
		phi_m = np.sum(fase*avg_prob)
		mean_fB = 1e-6*np.angle(phi_m)/(2*np.pi*self.t0)

		if do_plot:

			if show_plot:
				plt.ion()
			else:
				plt.ioff()
			f1 = plt.figure()
			plt.plot (beta/1e6, avg_prob)
			if y_log:
				plt.yscale('log')
				plt.ylim((1e-10,1e-1))
			plt.xlabel ('magnetic field detuning [MHz]')
			plt.ylabel ('prob distrib')
			plt.title('(B_exp = '+str('{0:.4f}'.format(mean_fB))+' +- '+str('{0:.4f}'.format(sigma_fB)) + ') MHz')
			if not(xlim==None):
				plt.xlim(xlim)

			if save_plot:
				f_name = 'probability_distribution.png'
				savepath = os.path.join(self.folder, f_name)
				f1.savefig(savepath)
			if show_plot:
				plt.show()
			plt.ion()

		if return_all_estimated_phases:
			list_phase_values = flatten(list_phase_values)
			return beta, avg_prob, ave_exp,H, mean_fB, sigma_fB, list_phase_values
		else:
			return beta, avg_prob, ave_exp,H, mean_fB, sigma_fB

	def compare_to_simulations(self, do_save = False, show_plot = False, verbose=False,plot_log=False):

		if show_plot:
			plt.ion()
		else:
			plt.ioff()

		f1 = plt.figure()
		
		beta_exp, p_exp, ave_exp,err_exp, mB, sB = self.mean_square_error(set_value=self.set_detuning, do_plot=False, show_plot=False, save_plot=False)
		plt.plot (beta_exp*1e-6, p_exp, color='RoyalBlue', label = 'exp')

		
		s = RamseySequence_Simulation (N_msmnts = self.N, reps=self.reps, tau0=self.t0)
		s.setup_simulation (magnetic_field_hz = self.set_detuning, G=self.G,F=self.F,K=self.N-1)
		s.verbose=verbose
		s.T2 = self.T2
		print 'T2 = ', self.T2
		s.fid0 = self.fid0
		s.fid1 = self.fid1
		if self.swarm_opt:
			s.sim_swarm_optim()
		else:
			s.sim_cappellaro_variable_M()
		s.convert_to_dict()
		if s.verbose:
			s.print_results()		
		beta_sim, p_sim, ave_exp_sim,err_sim, a, b = s.mean_square_error(set_value=self.set_detuning, do_plot=False, show_plot=False, save_plot=False)

		plt.plot (beta_sim*1e-6, p_sim, color='Crimson',label = 'sim')
		B_sim_string='\n (B_sim = '+str('{0:.4f}'.format(a))+' +- '+str('{0:.4f}'.format(b)) + ') MHz' + ';  H = ' + str('{0:.4f}'.format(err_sim))

		if plot_log:
			plt.yscale('log')
			plt.ylim((1e-10,0.5))
		plt.title('(B_exp = '+str('{0:.4f}'.format(mB))+' +- '+str('{0:.4f}'.format(sB)) + ') MHz' + ';  H = ' + str('{0:.4f}'.format(err_exp)) \
             + B_sim_string)		
		plt.xlabel ('magnetic field detuning [MHz]')
		plt.ylabel ('probability distribution')
		plt.legend()
		if do_save:
			f_name = 'probability_distribution.png'
			savepath = os.path.join(self.folder, f_name)
			print 'saving figure in ' , savepath
			f1.savefig(savepath)
		if show_plot:
			plt.show()
		plt.ion()
		return beta_exp, p_exp, ave_exp,err_exp, mB, sB



	def print_results (self):

		print '---    Measurement results  - SUMMARY    ---'
		print '--------------------------------------------'
		
		for k in self.msmnt_dict:
			curr_phase = self.phases_dict[k]
			curr_msmnt = self.msmnt_dict[k]
			mult = self.msmnt_multiplicity[k]

			print '(*)', curr_msmnt, ' - ', mult, ' times --- phases: ', np.round(curr_phase*180/(np.pi)), ' deg'
			print '     at indexes: ', self.index[k]			
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

	def setup_simulation(self, magnetic_field_hz = 0., G=0,F=0,K=0):
		self.beta_sim = magnetic_field_hz
		self.G = G
		self.F = F
		self.K = K
		if (os.name=='posix'):
			self.root_folder = '/home/cristian/Work/Research/teamdiamond/'
		else:
			self.root_folder = 'D:/measuring/'
	

	def save_folder (self, folder = '/home/cristian/Work/Research/adaptive magnetometry/'):
		self.save_folder = folder
	
	def ramsey (self, t=0., theta=0.):

		A = 0.5*(self.fid0 + self.fid1)
		B = 0.5*(self.fid1 - self.fid0)		

		p0 = (1-A)-B*np.exp(-(t/self.T2)**2)*np.cos(2*np.pi*self.beta_sim*t+theta)
		p1 = A+B*np.exp(-(t/self.T2)**2)*np.cos(2*np.pi*self.beta_sim*t+theta)

		np.random.seed()
		result = np.random.choice (2, 1, p=[p0, p1])
		return result[0]	

	def plot_ramsey(self, nr_datapoints=20, max_tau = None, theta=0):

		if (max_tau == None):
			max_tau = 2*self.T2
		t = np.linspace (0, max_tau, nr_datapoints)
		r = np.zeros(nr_datapoints)
		for i in np.arange(nr_datapoints):
			for j in np.arange(self.reps):
				r[i] = r[i] + self.ramsey(t=t[i], theta = theta)

		r = r/(float(self.reps))
		plt.plot (t*1e6, r, 'ob')
		plt.xlabel ('free evol time [us]')
		plt.ylabel ('ramsey')
		plt.ylim([0,1])
		plt.show()


		result = np.random.choice (2, 1, p=[p0, p1])
		
		return result[0]		


	def sim_cappellaro_variable_M (self):
		
		print 'Cappellaro!!!!!'		
		if self.verbose:				
			print '-------------------------------------------'
			print 'Simulating Cappellaro protocol (variable M)'
			print '-------------------------------------------'
			print '- N = '+str(self.N)+ ', G = '+str(self.G)+', F = '+str(self.F)

		if (self.G+self.F+self.K==0):
			print 'Simulation parameters G, K, F not set!!'
		else:			

			self.total_nr_msmnts = self.G*(2**(self.K+1)-1) + self.F*(2**(self.K+1)-2-self.K)
			self.msmnt_phases = np.zeros((self.reps,self.K+1))
			self.msmnt_times = np.zeros(self.K+1)
			self.msmnt_results = np.zeros((self.reps,self.K+1))
		
			k_array = self.K-np.arange(self.K+1)
			tau = 2**(k_array)
			self.msmnt_times = tau
			self.reset_rep_counter()

			for r in np.arange(self.reps):

				msmnt_results = np.zeros (self.K+1)
				t = np.zeros (self.K+1)
				phase = np.zeros(self.K+1)
				self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
				self.p_k [self.points] = 1/(2.*np.pi)
		
				for i,k in enumerate(k_array):

					t[i] = int(2**(k))
					ttt = -2**(k+1)
					
					phase[i] = 0.5*np.angle (self.p_k[ttt+self.points])
					m_total = 0
					MK = self.G+self.F*(self.K-k)
					#print 'MK=',MK,'time= ', t[i]
					for m in np.arange(MK):
						m_res = self.ramsey (theta=phase[i], t = t[i]*self.t0)#self.majority_vote_msmnt (theta_n = phase[i], t_n = t[i])					
						self.bayesian_update (m_n = m_res, phase_n = phase[i], t_n = 2**(k))
						m_total = m_total + m_res
					msmnt_results[i] =m_total
					
				self.msmnt_results [r, :] = np.copy(msmnt_results)
				self.msmnt_phases [r, :] = np.copy(phase)
				self.inc_rep()

	def sim_berry_protocol (self, do_adaptive):
		if self.verbose:				
			if do_adaptive:
				adptv = '(adaptive)'
			else:
				adptv = '(non adaptive)'
			print '-------------------------------------------'
			print 'Simulating Berry protocol '+adptv
			print '-------------------------------------------'
			print '- N = '+str(self.N)+ ', G = '+str(self.G)+', F = '+str(self.F)

		if (self.G+self.F+self.K==0):
			print 'Simulation parameters G, K, F not set!!'
		else:			
			self.phase_update = True
			self.total_nr_msmnts = self.G*(2**(self.K+1)-1) + self.F*(2**(self.K+1)-2-self.K)
			nr_results = int((self.K+1)*(self.G + self.F*self.K/2.))
			self.msmnt_phases = np.zeros((self.reps, nr_results))
			self.phase_upd_values = np.zeros((self.reps, nr_results))
			self.msmnt_times = np.zeros(nr_results)
			self.msmnt_results = np.zeros((self.reps,nr_results))
		
			k_array = self.K-np.arange(self.K+1)
			tau = 2**(k_array)
			self.reset_rep_counter()

			for r in np.arange(self.reps):
				msmnt_results = np.zeros (nr_results)
				t = np.zeros (self.K+1)
				phase = np.zeros(self.K+1)
				self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
				self.p_k [self.points] = 1/(2.*np.pi)
				res_idx = 0

				for i,k in enumerate(k_array):

					t[i] = int(2**(k))
					ttt = -2**(k+1)
					
					if do_adaptive:
						phase[i] = 0.5*np.angle (self.p_k[ttt+self.points])
					m_total = 0
					MK = self.G+self.F*(self.K-k)
					for m in np.arange(MK):
						phase_upd = m*np.pi/float(self.G + self.F*(self.K - k))
						m_res = self.ramsey (theta=phase[i]+phase_upd, t = t[i]*self.t0)					
						if do_adaptive:
							self.bayesian_update (m_n = m_res, phase_n = phase[i]+phase_upd, t_n = 2**(k))
						m_total = m_total + m_res
						self.msmnt_results[r, res_idx] = m_res
						self.msmnt_phases[r, res_idx] = (phase[i]+phase_upd)
						self.msmnt_times[res_idx] = t[i]
						res_idx = res_idx + 1
				#print 'Msmsnt results: ', self.msmnt_results
				#print 'Cappellaro phases: ', np.round(self.msmnt_phases*180/np.pi)
				#print 'Msmnt times: ', self.msmnt_times
				self.inc_rep()


	def sim_recalculate_optimal_phase (self, always_recalculate_phase, N1_sweep=False):
		#PARAMETERS:
		#always_recalculate_phase: calculate Cappellaro phase before each msmnt, not only when updating tau
		#N1_sweep: ONLY for N=1, update phase as m*pi/G

		if (self.G+self.F+self.K==0):
			print 'Simulation parameters G, K, F not set!!'
		else:			
			self.phase_update = True
			self.total_nr_msmnts = self.G*(2**(self.K+1)-1) + self.F*(2**(self.K+1)-2-self.K)
			nr_results = int((self.K+1)*(self.G + self.F*self.K/2.))
			self.msmnt_phases = np.zeros((self.reps, nr_results))
			self.phase_upd_values = np.zeros((self.reps, nr_results))
			self.msmnt_times = np.zeros(nr_results)
			self.msmnt_results = np.zeros((self.reps,nr_results))
		
			k_array = self.K-np.arange(self.K+1)
			tau = 2**(k_array)
			self.reset_rep_counter()

			for r in np.arange(self.reps):
				msmnt_results = np.zeros (nr_results)
				t = np.zeros (self.K+1)
				self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
				self.p_k [self.points] = 1/(2.*np.pi)
				res_idx = 0

				for i,k in enumerate(k_array):

					t[i] = int(2**(k))
					ttt = -2**(k+1)					
					m_total = 0
					MK = self.G+self.F*(self.K-k)

					phase = 0.5*np.angle (self.p_k[ttt+self.points])

					for m in np.arange(MK):
						if always_recalculate_phase:
							phase = 0.5*np.angle (self.p_k[ttt+self.points])
							if N1_sweep:
								if i==0:
									phase = m*np.pi/float(self.G)
						m_res = self.ramsey (theta=phase, t = t[i]*self.t0)					
						self.bayesian_update (m_n = m_res, phase_n = phase, t_n = 2**(k))
						self.msmnt_results[r, res_idx] = m_res
						self.msmnt_phases[r, res_idx] = phase
						self.msmnt_times[res_idx] = t[i]
						res_idx = res_idx + 1
				#print 'Msmsnt results: ', self.msmnt_results
				#print 'Cappellaro phases: ', np.round(self.msmnt_phases*180/np.pi)
				#print 'Msmnt times: ', self.msmnt_times
				self.inc_rep()



	def sim_swarm_optim (self):
		
		if (self.G+self.F+self.K==0):
			print 'Simulation parameters G, K, F not set!!'
		else:			
			self.phase_update = True #Need this to be set to True to perform analysis considering data sotred as each msmnt_result in one separate array location!
			self.always_recalculate_phase = True
			self.swarm_opt = True
			self.total_nr_msmnts = self.G*(2**(self.K+1)-1) + self.F*(2**(self.K+1)-2-self.K)
			nr_results = int((self.K+1)*(self.G + self.F*self.K/2.))
			self.msmnt_phases = np.zeros((self.reps, nr_results))
			self.phase_upd_values = np.zeros((self.reps, nr_results))
			self.msmnt_times = np.zeros(nr_results)
			self.msmnt_results = np.zeros((self.reps,nr_results))
		
			k_array = self.K-np.arange(self.K+1)
			tau = 2**(k_array)
			self.reset_rep_counter()

			file_old = 'D:/measuring/analysis/scripts/magnetometry/swarm_optimization/phases_G'+str(self.G)+'_F'+str(self.F)+'/swarm_opt_G='+str(self.G)+'_F='+str(self.F)+'_K='+str(self.K)+'.npz'
			round_fid = int(round(self.fid0*100))
			file_new = 'D:/measuring/analysis/scripts/magnetometry/swarm_optimization/incr_fid'+str(round_fid)+'_G'+str(self.G)+'/incr_fid'+str(round_fid)+'_G'+str(self.G)+'F'+str(self.F)+'_K='+str(self.K)+'.npz'
			if os.path.exists (file_new):
				swarm_incr_file = file_new
			elif os.path.exists (file_old):
				swarm_incr_file = file_old
			else:
				print 'ATTENTION!!! No file found!'
				print file_new
				print file_old
			swarm_opt_pars = np.load (swarm_incr_file)
			self.u0 = swarm_opt_pars['u0']
			self.u1 = swarm_opt_pars['u1']			

			for r in np.arange(self.reps):
				msmnt_results = np.zeros (nr_results)
				t = np.zeros (self.K+1)
				self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
				self.p_k [self.points] = 1/(2.*np.pi)
				res_idx = 0
				m_res = 0

				for i,k in enumerate(k_array):

					t[i] = int(2**(k))
					ttt = -2**(k+1)					
					m_total = 0
					MK = self.G+self.F*(self.K-k)

					for m in np.arange(MK):
						if m_res == 0:
							phase_inc_swarm = self.u0 [res_idx]
						else:
							phase_inc_swarm = self.u1 [res_idx]

						phase_cappellaro = 0.5*np.angle (self.p_k[ttt+self.points])
						phase = phase_cappellaro + phase_inc_swarm
						m_res = self.ramsey (theta=phase, t = t[i]*self.t0)					
						self.bayesian_update (m_n = m_res, phase_n = phase, t_n = 2**(k))
						self.msmnt_results[r, res_idx] = m_res
						self.msmnt_phases[r, res_idx] = phase
						self.msmnt_times[res_idx] = t[i]
						res_idx = res_idx + 1
				self.inc_rep()

	def sim_SQL(self, phase_deg):
		#how do we measure the SQL data, just hacking the swarm opt code??
		self.phase_update = True #Need this to be set to True to perform analysis considering data sorted as each msmnt_result in one separate array location!
		self.always_recalculate_phase = True
		self.swarm_opt = True
		self.total_nr_msmnts = self.G*(2**(self.K+1)-1) + self.F*(2**(self.K+1)-2-self.K)
		nr_results = int((self.K+1)*(self.G + self.F*self.K/2.))
		self.msmnt_phases = np.zeros((self.reps, nr_results))
		self.phase_upd_values = np.zeros((self.reps, nr_results))
		self.msmnt_times = np.zeros(nr_results)
		self.msmnt_results = np.zeros((self.reps,nr_results))
	
		k_array = self.K-np.arange(self.K+1)
		tau = 2**(k_array)
		self.reset_rep_counter()

		for r in np.arange(self.reps):
			msmnt_results = np.zeros (nr_results)
			t = np.zeros (self.K+1)
			self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
			self.p_k [self.points] = 1/(2.*np.pi)
			res_idx = 0
			m_res = 0

			for i,k in enumerate(k_array):

				t[i] = int(2**(k))
				ttt = -2**(k+1)					
				m_total = 0
				MK = self.G+self.F*(self.K-k)

				for m in np.arange(MK):

					phase = phase_deg*np.pi/180. #HACK for SQL!
					m_res = self.ramsey (theta=phase, t = t[i]*self.t0)					
					self.bayesian_update (m_n = m_res, phase_n = phase, t_n = 2**(k))
					self.msmnt_results[r, res_idx] = m_res
					self.msmnt_phases[r, res_idx] = phase
					self.msmnt_times[res_idx] = t[i]
					res_idx = res_idx + 1
			self.inc_rep()

		rows, cols = np.shape(self.msmnt_results)
		self.stored_msmnt_results = np.copy (np.reshape(self.msmnt_results, rows*cols))
		self.stored_msmnt_phases = np.copy (np.reshape(self.msmnt_phases, rows*cols))
		self.stored_msmnt_times = np.copy (self.msmnt_times[0])



	def sim_nonadaptive_variable_M (self):
		
		if self.verbose:				
			print '---------------------------------------------'
			print 'Simulating non-adaptive protocol (variable M)'
			print '---------------------------------------------'
			print '- N = '+str(self.N)+ ', G = '+str(self.G)+', F = '+str(self.F)
		
		if (self.G+self.F+self.K==0):
			print 'Simulation parameters G, K, F not set!!'
		else:			

			self.total_nr_msmnts = self.G*(2**(self.K+1)-1) + self.F*(2**(self.K+1)-2-self.K)
			self.msmnt_phases = np.zeros((self.reps,self.K+1))
			self.msmnt_times = np.zeros(self.K+1)
			self.msmnt_results = np.zeros((self.reps,self.K+1))
		
			k_array = self.K-np.arange(self.K+1)
			tau = 2**(k_array)
			self.msmnt_times = tau
			self.reset_rep_counter()

			for r in np.arange(self.reps):

				msmnt_results = np.zeros (self.K+1)
				t = np.zeros (self.K+1)
				phase = np.zeros(self.K+1)
				self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
				self.p_k [self.points] = 1/(2.*np.pi)
		
				for i,k in enumerate(k_array):

					t[i] = int(2**(k))
					ttt = -2**(k+1)
					
					phase[i] = 0.5*np.angle (self.p_k[ttt+self.points])
					m_total = 0
					MK = self.G+self.F*(self.K-k)
					#print 'MK=',MK,'time= ', t[i]
					for m in np.arange(MK):
						m_res = self.ramsey (theta=phase[i], t = t[i]*self.t0)#self.majority_vote_msmnt (theta_n = phase[i], t_n = t[i])					
						m_total = m_total + m_res
					msmnt_results[i] =m_total
					
				self.msmnt_results [r, :] = np.copy(msmnt_results)
				self.msmnt_phases [r, :] = np.copy(phase)
				self.inc_rep()


			
	def load_table (self, N, M):
		# could probably remove this function
		ttt = int(np.round(self.t0*1e9))
		self.table_folder = self.root_folder+'/measurement/scripts/Magnetometry/adaptive_tables_lt1/tau0='+str(ttt)+'ns/'		
		name = 'adptv_table_cappellaro_N='+str(self.N)+'_M='+str(self.M)+'.npz'
		a = np.load(self.table_folder+name)
		self.table = a['table']


	def table_based_simulation (self):
		# could probably remove this function
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
					T2=self.T2, fid0=self.fid0, fid1=self.fid1, B_sim = self.beta_sim, tau0 =self.t0, N=self.N,K=self.K,G=self.G,F=self.F, reps = self.reps)
		
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
		self.K=a['K']
		self.F=a['F']
		self.G=a['G']
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

		a = sequence.MagnetometrySequenceAnalysis(self.folder)
		#a.get_sweep_pts()
		a.get_magnetometry_data(name='adwindata', ssro = False)
		print 'a keys ', a.N
		self.msmnt_results = a.clicks
		#print 'msmnt_results (load_exp_data): ', self.msmnt_results

		#if ((np.shape(np.shape(a.clicks)))[0]==1):
		#	self.reps = len(a.clicks)
		#	self.N = 1
		#else:
		#	self.reps, self.N = np.shape (a.clicks)
		self.N=a.N
		self.M = a.M
		self.F=a.F
		self.G=a.G
		self.K=a.K
		self.repetitions=a.repetitions
		self.reps=a.reps
		self.n_points =  2**(self.N+3)
		self.t0 = a.t0
		self.B_max = 1./(2*self.t0)
		self.msmnt_times_tmp = np.zeros(len(a.ramsey_time))
		self.set_detuning = a.set_detuning
		self.CR_after = a.CR_after
		try:
			self.do_add_phase=a.do_add_phase
			if self.do_add_phase==0:
				self.phase_update=False
			else:
				self.phase_update=True	
		except:	
			self.phase_update=False

		self.swarm_opt = a.swarm_opt
		self.p_tn = a.p_tn
		self.p_2tn = a.p_2tn

		for j in np.arange(len(a.ramsey_time)):
			self.msmnt_times_tmp[j] = a.ramsey_time[j]/self.t0
		self.theta_rad = 2*np.pi*a.theta/360.
		print self.phase_update
		if self.phase_update:
			self.msmnt_phases=self.msmnt_results*0.
			self.msmnt_times=self.msmnt_results[1]*0.
			for rep in np.arange(self.reps):
				i=0
				for n in np.arange(self.N):
					MK=self.G+self.F*n
					for m in np.arange(MK):
						added_phase=m*np.pi/float(self.G+self.F*n)
						self.msmnt_phases[rep,i]=self.theta_rad[rep,n]+added_phase
						if rep==0:
							#print 'added phase = ',added_phase, '; theta rad = ', self.theta_rad[rep,n]
							#print self.msmnt_phases[0,:]
							self.msmnt_times[i]=self.msmnt_times_tmp[n]
						i+=1
		elif self.swarm_opt:
			print 'Loading exp data: formatting time array for swarm opt!'
			self.msmnt_phases=self.theta_rad
			i = 0 
			if self.swarm_opt:
				self.msmnt_times=self.msmnt_results[1]*0.
				for n in np.arange(self.N):
					MK=self.G+self.F*n
					for m in np.arange(MK):
						self.msmnt_times[i]=self.msmnt_times_tmp[n]
						i = i + 1
		else:
			self.msmnt_phases=self.theta_rad
			self.msmnt_times=self.msmnt_times_tmp

		#print a.theta

		self.discarded_elements = []
		phases_detuning = 2*np.pi*a.phases_detuning/360.
		b = np.ones(self.reps)
		#self.msmnt_phases = np.mod(self.msmnt_phases - np.outer (b, phases_detuning), 2*np.pi)

		self.msmnt_type = a.msmnt_type
		self.timer = a.timer

		if (a.msmnt_type=='realtime'):
			self.test_adwin = a.test_adwin
			self.G = a.G
			self.F = a.F
			self.renorm_ssro = a.renorm_ssro
			self.exp_fid0 = a.exp_fid0
			self.exp_fid1 = a.exp_fid1
			self.opt_phase = a.theta_opt #to be removed, only for testing yesterday's data!!! (should load a.theta_opt)

		self.T2 = a.T2_mult_t0*self.t0
		self.save_pk_n = a.save_pk_n
		self.save_pk_m = a.save_pk_m
		self.real_pk_adwin = a.real_pk_adwin
		self.imag_pk_adwin = a.imag_pk_adwin



	def CR_after_postselection(self,threshold=15):

		if (self.N>1):
			res = np.copy(self.msmnt_results)
			phases = np.copy(self.msmnt_phases)
			#print np.shape(phases)
			#print np.shape(self.msmnt_results)
			self.discarded_elements = []
			new_results = self.msmnt_results*0.
			new_phases = self.msmnt_results*0.
			rep = 0
			self.CR_after_threshold = threshold
			for j in np.arange(self.reps):
				if any(t<self.CR_after_threshold for t in self.CR_after[j,:]):
					self.discarded_elements.append(j)
				else:
				#print 'for i = ', j , 'CR array',self.CR_after[j,:]
					new_results[rep,:] = np.copy(res[j,:])
					new_phases[rep,:] = np.copy(phases[j,:])
					rep = rep + 1

				#if (len(self.CR_after[j,:])==np.count_nonzero(self.CR_after[j,:])):
					#print 'for i = ', j , 'CR array',self.CR_after[j,:]
				#	new_results[rep,:] = np.copy(res[j,:])
				#	new_phases[rep,:] = np.copy(phases[j,:])
				#	rep = rep + 1
				#else:
				#	self.discarded_elements.append(j)
			self.reps = rep
			self.msmnt_results =new_results[:self.reps,:]
			self.msmnt_phases =  new_phases[:self.reps,:]
			if len(self.discarded_elements)>0:
				print 'CR post-selection, discard elements: ', self.discarded_elements


	def check_realtime_phases(self):
		adwin_phase = self.opt_phase[:self.N]
		adwin_results = self.msmnt_results[0,:]
		sys.stdout.write ('Msmnt type: '+self.msmnt_type+ '   --- N = '+str(self.N)+ ', M = '+str(self.M))
		if self.renorm_ssro:
			print '	- ssro renorm: = ',self.exp_fid0, ', fid1 = ', self.exp_fid1,')'
		else:
			print '	- no ssro renorm'

		print 'Adwin: results = ', adwin_results, ' - phases: ', adwin_phase

		t = AdaptiveTable (N=self.N,M=self.M)
		t.set_tau0(tau0=self.t0)
		t.verbose = False
		sim_ph, pos = t.msmnt_to_position (msmnt_results = adwin_results)
		print 'Python simulated phases: ', sim_ph

		diff = np.sum(np.abs(adwin_phase-sim_ph))/self.N
		print 'MEAN ERROR: ', diff, ' deg'


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
		self.gamma_e = 28e9 #~28 GHz/T

		self.nr_periods = None
		self.nr_points_per_period = None
		
		self.scaling_variance=[]
		self.total_time=[]
		self.sensitivity = []

		self.scaling_factor = None
		self.error_scaling_factor = None
		self.simulated_data = None
		self.error_bars = False

		self.protocols = {}
		capp_ptcl = {'do_adaptive':True, 'phase_update':False, 'always_recalculate_phase':False, 'swarm_optimization': False, 'code':'capp'}
		mod_capp_ptcl = {'do_adaptive':True, 'phase_update':False, 'always_recalculate_phase':True, 'swarm_optimization': False, 'code':'modCapp'}
		nonadptv_ptcl = {'do_adaptive':False, 'phase_update':True, 'always_recalculate_phase':False, 'swarm_optimization': False, 'code':'nnAdptv'}
		phaseUpdate_ptcl = {'do_adaptive':True, 'phase_update':True, 'always_recalculate_phase':False, 'swarm_optimization': False, 'code':'adptvPhUpdate'}
		swarmOptim_ptcl = {'do_adaptive':True, 'phase_update':True, 'always_recalculate_phase':True, 'swarm_optimization': True, 'code':'swarmOpt'}
		self.protocols['cappellaro'] = capp_ptcl
		self.protocols['modified_cappellaro'] = mod_capp_ptcl
		self.protocols['non_adaptive'] = nonadptv_ptcl
		self.protocols['cappellaro_phase_update'] = phaseUpdate_ptcl
		self.protocols['swarm_optimization'] = swarmOptim_ptcl
		self.do_sql = False

		if (os.name =='posix'):
			self.folder = '/home/cristian/Work/Research/adaptive magnetometry/data_analysis/'
		else:
			folder = r'M:/tnw/ns/qt/Diamond/Projects/Magnetometry with adaptive measurements/Data/analyzed data'
			if (os.path.exists(folder)):
				self.folder = folder
			else:
				self.folder = r'V:/tnw/ns/qt/Diamond/Projects/Magnetometry with adaptive measurements/Data/analyzed data'

	def set_exp_params(self, T2 = 96e-6, fid0 = 0.88, fid1 = 0.015):
		self.T2 = T2
		self.fid0 = fid0
		self.fid1 = fid1
		
	def set_protocol (self, G=0,K=0,F=0):
		self.G=G
		self.F=F
		self.K=K
		self.theta = 0*np.pi/180.

	def set_sweep_params (self, nr_periods = 3, nr_points_per_period=11, reps=200):
		self.nr_periods = nr_periods
		self.nr_points_per_period = nr_points_per_period
		self.reps = reps

	def sample_B_space(self, N):

		delta_f = 1./(2*self.t0*(2**N))
		nr_available_periods = 2**(N)

		if (self.nr_periods == None):
			self.nr_periods = int(N^(3/2))

		if (self.nr_periods>nr_available_periods):
			nr_periods = nr_available_periods
		else:
			nr_periods = self.nr_periods

		periods = (np.unique(np.random.randint(0, nr_available_periods, size=nr_periods*3)))-nr_available_periods/2

		if (len(periods)>nr_periods):
			idxs = np.random.randint (0, len(periods), size = nr_periods)
			periods = periods [idxs]
		periods = np.sort(periods)
		B_values = np.array([])
		label_array = []
		for per in periods:
		    B = np.linspace(per*delta_f, (per+1)*delta_f, self.nr_points_per_period)
		    B_values = np.hstack((B_values, B))
		    for l in np.arange(self.nr_points_per_period):
		    	label_array.append('N='+str(N)+'_p'+str(per)+'_'+str(l))
		return label_array, B_values

	def specific_B_fields_for_sweep_sim(self,N,pts):
		B_dict={'1':np.linspace(-25.e6,25e6,pts),
				'2':np.linspace(-12.5e6,25e6,pts),
				'3':np.linspace((12.5-6.25)*1e6,(18.75+6.25)*1e6,pts),
				'4':np.linspace((3.125-3.125)*1e6,(6.25+3.125)*1e6,pts),
				'5':np.linspace((-14.0625-1.5625)*1e6,(-12.5+1.5625)*1e6,pts),
				'6':np.linspace((0.78125-0.78125)*1e6,(1.5625+0.78125)*1e6,pts),
				'7':np.linspace((-21.09375-0.390625)*1e6,(-20.703125+0.390625)*1e6,pts),
				'8':np.linspace((2.9296875-0.1953125)*1e6,(3.125+0.1953125)*1e6,pts),
				'9':np.linspace((23.92578125-0.09765625)*1e6,(24.0234375+0.09765625)*1e6,pts),
				'10':np.linspace((-1.953125-0.048828125)*1e6,(-1.904296875+0.048828125)*1e6,pts),
				'11':np.linspace((0.29296875-0.0244140625)*1e6,(0.3173828125+0.0244140625)*1e6,pts),
				'12':np.linspace((-21.7407226563-0.012207031300000892)*1e6,(-21.728515625+0.012207031300000892)*1e6,pts),
				'13':np.linspace((1.30004882813-0.00610351562)*1e6,(1.30615234375+0.00610351562)*1e6,pts),
				'14':np.linspace((-7.58972167969-0.003051757810000666)*1e6,(-7.58666992188+0.003051757810000666)*1e6,pts),
				
		}

		return B_dict[str(N)]
	
	def sweep_field_simulation (self, N, protocol, table_based=False, print_results=False, specific_B=False):
		self.simulated_data = True		
		self.analyzed_N.append(N)	
		
		B_values = np.array([])
		label_array = []

		per=0
		delta_f = 1./(self.t0*(2**N))

		B = np.linspace(per*delta_f, (per+1)*delta_f, self.nr_points_per_period)
		if specific_B:
			B=self.specific_B_fields_for_sweep_sim(N,self.nr_points_per_period)
		self.B_values = np.hstack((B_values, B))
		for l in np.arange(self.nr_points_per_period):
			label_array.append('N='+str(N)+'G='+str(self.G)+'F='+str(self.F)+'_p'+str(0)+'_'+str(l))

		label_array, self.B_values = self.sample_B_space (N=N)		
		self.B_values=np.unique(self.B_values)
		nr_points = len(self.B_values)
		msqe = np.zeros(nr_points)
		ave_exps = [0]*nr_points
		B_field = np.zeros(nr_points)

		ind = 0
		print 'PROTOCOL: ', self.protocols[protocol]
		print "Simulating N="+str(N)+', '+str(len(self.B_values))+" instances of magnetic field"

		self.phase_update = self.protocols[protocol]['phase_update']
		self.do_adaptive = self.protocols[protocol]['do_adaptive']
		self.always_recalculate_phase = self.protocols[protocol]['always_recalculate_phase']
		self.swarm_optim = self.protocols[protocol]['swarm_optimization']
		self.code = self.protocols[protocol]['code']
		self.protocol = protocol

		list_estim_phases = []
		self.N1_sweep = False
		for b in np.arange(nr_points):
			sys.stdout.write(str(ind)+', ')	
			s = RamseySequence_Simulation (N_msmnts = N, reps=self.reps, tau0=self.t0)
			s.setup_simulation (magnetic_field_hz = self.B_values[b], F=self.F,G=self.G,K=self.K) 
			s.verbose = False
			s.T2 = self.T2
			s.fid0 = self.fid0
			s.fid1 = self.fid1
	
			if protocol == 'cappellaro':
				s.sim_cappellaro_variable_M()
			elif protocol == 'modified_cappellaro':
				s.sim_recalculate_optimal_phase(always_recalculate_phase=True)
			elif protocol == 'non_adaptive':
				s.sim_berry_protocol(do_adaptive=False)
			elif protocol == 'cappellaro_phase_update':
				s.sim_berry_protocol(do_adaptive=True)
			elif protocol == 'swarm_optimization':
				s.sim_swarm_optim()
			else:
				print 'Unknown protocol: ', protocol

			s.convert_to_dict()

			if print_results:
				s.print_results()

			if self.error_bars:
				beta, p, ave_exp,H, mB, sB, list_phase_values = s.mean_square_error(set_value=self.B_values[b], do_plot=False, return_all_estimated_phases = True)
				list_estim_phases.append(list_phase_values)
			else:				
				beta, p, ave_exp,H, mB, sB = s.mean_square_error(set_value=self.B_values[b], do_plot=False)
			self.prob_density_dict[label_array[ind]] = p
			ave_exps[ind]=ave_exp
			msqe [ind] = H
			B_field [ind] = self.B_values[b]
			ind =ind+1
		list_estim_phases = flatten(list_estim_phases)
		self.results_dict[str(N)] = {'B_field':B_field, 'ave_exp':ave_exps,'msqe':msqe, 'G':self.G,'K':self.K,'F':self.F, 'estimated_phase_values':list_estim_phases}

	def load_sweep_field_data (self, N, compare_to_simulations=False,older_than=None,newer_than=None,CR_after_threshold=2):

		self.simulated_data = False
		self.analyzed_N.append(N)
		msqe = np.zeros(self.nr_points_per_period*self.nr_periods)
		ave_exps = [0]*self.nr_points_per_period*self.nr_periods
		B_field = np.zeros(self.nr_points_per_period*self.nr_periods)
		ind = 0
		check_params_labels = ['tau0','F','G']
		self.nr_discarded_elements=[]
		list_estim_phases=[]
		for per in np.arange(self.nr_periods):
			for pt in np.arange (self.nr_points_per_period):
				if self.phase_update:
					label = '_N = '+str(N)+'_'+'M=('+str(self.G)+', '+str(self.F)+')'+'_addphase_rtAdwin_'+'_p'+str(per)+'b'+str(pt)
				elif self.swarm_opt:
					label = '_N = '+str(N)+'_'+'M=('+str(self.G)+', '+str(self.F)+')'+'_swarm'+'_p'+str(per)+'b'+str(pt)	
				else:
					label = '_N = '+str(N)+'_'+'M=('+str(self.G)+', '+str(self.F)+')'+'_rtAdwin_'+'_p'+str(per)+'b'+str(pt)
				print 'Processing...', label
				if older_than:
					f = toolbox.latest_data(contains=label,older_than=older_than,newer_than=newer_than)
				else:
					f=toolbox.latest_data(contains=label)	
				s = RamseySequence_Exp (folder = f)
				s.set_exp_pars (T2=96e-6, fid0=0.876, fid1=1-.964)
				s.load_exp_data()
				#print np.shape(s.msmnt_results)
				#print s.msmnt_results
				self.repetitions = s.repetitions
				s.CR_after_postselection(CR_after_threshold)
				self.CR_after_threshold=s.CR_after_threshold
				self.nr_discarded_elements.append(len(s.discarded_elements))
				check_params = [(s.t0 == self.t0), (s.F == self.F),(s.G == self.G)]
				if np.all(check_params):
					s.convert_to_dict()
					if compare_to_simulations:
						beta, prob, err, mB, sB = s.compare_to_simulations (show_plot=False, do_save=True, verbose=False)
					else:
						#print 'ERROR BARS: ', self.error_bars
						if self.error_bars:
							beta, p, ave_exp,H, mB, sB, list_phase_values = s.mean_square_error(set_value=s.set_detuning, do_plot=False, return_all_estimated_phases = True)
							list_estim_phases.append(list_phase_values)
						else:				
							beta, p, ave_exp,H, mB, sB = s.mean_square_error(set_value=s.set_detuning, do_plot=False)

						beta, prob, ave_exp,err, mB, sB = s.mean_square_error(show_plot=False, save_plot=True, do_plot=False)
					self.prob_density_dict[label] = prob
					#print s.set_detuning, mB
					ave_exps[ind]=ave_exp
					msqe [ind] = err
					B_field [ind] = s.set_detuning
				else:
					msg = []
					for i in np.arange (4):
						if (check_params[i]==False): msg.append(check_params_labels[i])
					print 'Non matching parameters: ', msg, ' --- ', label
				ind +=1
		list_estim_phases = flatten(list_estim_phases)
		self.results_dict[str(N)] ={'B_field':B_field, 'ave_exp':ave_exps,'msqe':msqe, 'G':self.G,'K':self.K,'F':self.F, 'estimated_phase_values':list_estim_phases,'nr_discarded_elements':self.nr_discarded_elements}

	def plot_msqe_dictionary(self,y_log=False, save_plot=False):

		if self.simulated_data:
			fName = time.strftime ('%Y%m%d_%H%M%S')+'_simulated_adaptive_magnetometry_'+'N='+str(self.N)+'G='+str(self.G)+'F='+str(self.F)+'_fid0='+str(self.fid0)+'.png'
		else:
			fName = time.strftime ('%Y%m%d_%H%M%S')+'_adaptive_magnetometry_'+'N='+str(self.N)+'G='+str(self.G)+'F='+str(self.F)+'_fid0='+str(self.fid0)+'.png'

		C = compare.compare_functions()
		C.xlabel = 'magnetic field detuning [MHz]'
		#C.ylabel = 'mean square error [MHz^2]'
		C.ylabel = '$V_{H}$'
		for n in self.analyzed_N:	
			C.add (x =self.results_dict[str(n)]['B_field']*1e-6, y=self.results_dict[str(n)]['msqe'], label=str(n))
		fig = C.plot(y_log=y_log)
		if save_plot:
			savepath = os.path.join(self.folder, fName)
			fig.savefig(savepath)


	def sweep_field (self, do_simulate = True):
	
		for n in np.arange(self.N)+1:
			self.sweep_field_fixedN (N=n, do_simulate = do_simulate)
			

	def calculate_scaling (self):
		print 'Calculating scaling ... '
		self.scaling_variance=[]
		self.total_time=[]
		self.std_H = []
		for i,n in enumerate(self.analyzed_N):
			
			msqe_phi = self.results_dict[str(n)]['ave_exp']
			self.scaling_variance.append(np.abs(np.mean(msqe_phi))**(-2)-1)

			if self.error_bars:
				bs = stat.BootStrap (n_boots = 1000)
				estim_phases = np.array(self.results_dict[str(n)]['estimated_phase_values'])
				bs.set_y (y=estim_phases)
				bs.run_bootstrap_holevo()
				self.std_H.append(bs.errH_bootstrap)

			if (self.G+self.F+self.K==0):
				print 'Error G F and K are not set!'
				break
			else:	
				if self.do_sql:
					self.total_time.append(self.t0*n)
					print 'current G = ', n, '   ---- time: ', self.t0*n*1e6, 'us'
				else:
					self.total_time.append(self.t0*(self.G*(2**(n)-1)+self.F*(2**(n)-1-n)))

			#print np.array(self.total_time)/self.t0
		self.total_time = np.array(self.total_time)
		self.scaling_variance=np.array(self.scaling_variance)		
		self.sensitivity = (self.scaling_variance*self.total_time)#/((2*np.pi*self.gamma_e*self.t0)**2)
		if self.error_bars:
			self.std_H  = np.array(self.std_H)
			self.err_sensitivity = self.std_H*self.total_time
		else:
			self.std_H  = self.sensitivity*0
			self.err_sensitivity = self.sensitivity*0

	def plot_sensitivity_scaling (self, do_fit = True, save_plot=False):
		if (self.scaling_variance == []):
			self.calculate_scaling()

		if self.simulated_data:
			if self.always_recalculate_phase:
				add = '_alwRecPhase'
			else:
				add = ''
			fName = time.strftime ('%Y%m%d_%H%M%S')+'_simulated_adaptive_magnetometry_'+'N='+str(self.N)+'G='+str(self.G)+'F='+str(self.F)+'_fid0='+str(self.fid0)+add+'.png'
		else:
			fName = time.strftime ('%Y%m%d_%H%M%S')+'_adaptive_magnetometry_'+'N='+str(self.N)+'G='+str(self.G)+'F='+str(self.F)+'_fid0='+str(self.fid0)+'.png'


		print 'PLOTTING SCALING!!!!!!'
		plt.figure()


		#NOTE!!!

		plt.loglog (self.total_time, self.sensitivity, 'ob')
		if self.error_bars:
			#plt.loglog (self.total_time, self.sensitivity+self.err_sensitivity, ':r')
			#plt.loglog (self.total_time, self.sensitivity-self.err_sensitivity, ':r')		
			plt.fill_between (self.total_time, self.sensitivity-self.err_sensitivity, self.sensitivity+self.err_sensitivity, color='b', alpha=0.1)
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel ('total ramsey time [$\mu$s]', fontsize=15)
		plt.ylabel ('sensitivity [$\mu$T$^2$*Hz$^{-1}$]', fontsize=15)
		plt.ylabel ('$V_{H}$ T ')
		plt.show()


		
		#x0 = np.log10(self.total_time*1e6)
		#y0 = np.log10(self.sensitivity*1e12)
		#NOTE!!!!!!!!
		x_full = np.log10(self.total_time/self.t0)
		y_full = np.log10(self.sensitivity/self.t0)
		err_y = self.err_sensitivity/self.sensitivity
		#x0 = self.total_time*1e6
		#y0 = self.sensitivity*1e12
		a='y'
		#a = raw_input('Do you want to use a sub-set of points for the fit? [y/n]')
		if (a=='y'):
			n0 = 1#raw_input ('First [1-'+str(len(self.analyzed_N))+'] = ?')
			n1 = len(self.analyzed_N)#raw_input ('Last [1-'+str(len(self.analyzed_N))+'] = ?')
			print len(self.analyzed_N)
			n0 = int(n0)
			n1 = int(n1)

			#Note!!!
			#x0 = np.log10(self.total_time[n0-1:n1-1]*1e6)
			#y0 = np.log10(self.sensitivity[n0-1:n1-1]*1e12)
			
			x0 = np.log10(self.total_time[n0:n1]/self.t0)
			y0 = np.log10(self.sensitivity[n0:n1]/self.t0)

			#x = self.total_time*1e6
			#NOTE!!!
			x = x0


			y = y0
		

		try:
			guess_b = -1#np.abs((y0[1]-y0[0]))/(x0[1]-x0[0])
			guess_a = 0*y0[-1]+1
			print 'aguess = ',guess_a, ' guess_b = ', guess_b
			
			a = fit.Parameter(guess_a, 'a')
			b = fit.Parameter(guess_b, 'b')
			
			p0 = [a, b]
			fitfunc_str = ''

			def fitfunc(x):
				return a()+b()*x

			if self.error_bars:
				error_array = err_y[n0:n1]
			else:
				error_array = []

			fit_result = fit.fit1d(x0,y0, None, p0=p0, fitfunc=fitfunc, fixed=[],
                	do_print=False, ret=True, err_y = error_array)
			a_fit = fit_result['params_dict']['a']
			b_fit = fit_result['params_dict']['b']
			b_err = fit_result['error_dict']['b']
			print 'a= ',a_fit
			print 'b=',b_fit
			#NOTE!!!
			#x0 = 1e6*self.total_time[0]/2.
			#x_end = 1e6*self.total_time[-1]*2.
			

			x_fit = np.linspace (x0[0], x0[-1], 100)
			y_fit = a_fit+b_fit*x_fit
			y_guess=guess_a+guess_b*x_fit
			do_fit = True
		except:
			print 'Fit failed!'
			do_fit = False


		fig = plt.figure(figsize=(8,6))
		p = fig.add_subplot(111)
		p.tick_params(axis='both', which='major', labelsize=15)
		
		if do_fit:
			#y_SQL=a_fit*(x_fit**(0))
			#y_heis=a_fit*(x_fit**(-1))
			#p.plot(x_fit, y_heis, 'Grey')
			#p.plot (x_fit, y_SQL, 'Grey')
			self.scaling_factor = b_fit
			self.error_scaling_factor = b_err

			p.plot (x_fit, y_fit, 'r')
			#p.plot (x_fit, y_guess, 'b')
			plt.title('scaling:  '+str('{0:.2f}'.format(b_fit))+' +- '+str('{0:.2f}'.format(b_err)) + '$\mu$T*Hz$^{1/2}$', fontsize=15)
		else:
			self.scaling_factor = 1
			self.error_scaling_factor = 1
		#NOTE!!!!!!!!!!
		#p.loglog (self.total_time*1e6, self.sensitivity*1e12, 'o', markersize=10, markeredgecolor = 'k', markerfacecolor='b')
		if self.error_bars:
			p.errorbar (x_full, y_full, yerr= err_y, fmt='o', markersize=10, markeredgecolor = 'k', markerfacecolor='b')
		else:
			p.plot (x_full, y_full, 'o', markersize=10, markeredgecolor = 'k', markerfacecolor='b')			
		plt.xlabel ('Log(total ramsey time$ * tau_{0}^{-1}$)', fontsize=15)
		plt.ylabel ('sensitivity [$\mu$T$^2$*Hz$^{-1}$]', fontsize=15)
		plt.ylim ([min(y_full)-0.2, max(y_full)+0.2])
		plt.ylabel ('Log($V_{H}$ T)')

		fig = plt.figure(figsize=(8,6))
		p = fig.add_subplot(111)
		p.tick_params(axis='both', which='major', labelsize=15)
		#p.loglog(x_fit, y_heis, 'Grey')
		#p.loglog(x_fit, y_SQL, 'Grey')
		if do_fit:
			p.plot (x_fit, y_fit, 'r')
			#p.loglog (x_fit, y_guess, 'b')
			pars = '{F='+str(self.F)+ ', G='+str(self.G)+'}'
			plt.title(pars+'scaling:  '+str('{0:.2f}'.format(b_fit))+' +- '+str('{0:.2f}'.format(b_err)) + '$\mu$T*Hz$^{1/2}$', fontsize=15)
		
		#NOTE!!!!!!!!!!
		#p.loglog (self.total_time*1e6, self.sensitivity*1e12, 'o', markersize=10, markeredgecolor = 'k', markerfacecolor='b')
		if self.error_bars:
			p.errorbar (x_full, y_full, yerr= err_y, fmt='o', markersize=10, markeredgecolor = 'k', markerfacecolor='b')
		else:
			p.plot (x_full, y_full, 'o', markersize=10, markeredgecolor = 'k', markerfacecolor='b')			
		

		plt.xlabel ('total ramsey time$ * tau_{0}^{-1}$', fontsize=15)
		plt.ylabel ('sensitivity [$\mu$T$^2$*Hz$^{-1}$]', fontsize=15)
		plt.ylabel ('$V_{H}$ T')


		if save_plot:
			savepath = os.path.join(self.folder, fName)
			fig.savefig(savepath)

		plt.show()
		

	def save(self, name = '', folder = None):
		
		if (folder==None):
			folder = self.folder 

		if self.simulated_data:

			try:
				fName = time.strftime ('%Y%m%d_%H%M%S')+ '_'+self.code+'_N='+str(self.N)+'G='+str(self.G)+'F='+str(self.F)+'_fid0='+str(self.fid0)
			except:

				if self.do_adaptive:
					fName = time.strftime ('%Y%m%d_%H%M%S')+'_simulated_adaptive_magnetometry_'+'N='+str(self.N)+'G='+str(self.G)+'F='+str(self.F)+'_fid0='+str(self.fid0)
					if self.always_recalculate_phase:
						fName = fName + '_alwaysRecPhase'
				else:
					fName = time.strftime ('%Y%m%d_%H%M%S')+'_simulated_NON_adaptive_magnetometry_'+'N='+str(self.N)+'G='+str(self.G)+'F='+str(self.F)+'_fid0='+str(self.fid0)

		else:
			fName = time.strftime ('%Y%m%d_%H%M%S')+'_adaptive_magnetometry_'+'N='+str(self.N)+'G='+str(self.G)+'F='+str(self.F)+'_fid0='+str(self.fid0)

		try:
			if self.N1_sweep:
				fName = fName +"_N1sweep"
		except:
			pass

		fName = fName+name
		if not os.path.exists(os.path.join(folder, fName+'.hdf5')):
			mode = 'w'
		else:
			mode = 'r+'
			print 'Output file already exists!'

		f = h5py.File(os.path.join(folder, fName+'.hdf5'), mode)
		f.attrs ['N'] = self.N
		f.attrs ['K'] = self.K
		f.attrs ['F'] = self.F
		f.attrs ['G'] = self.G
		f.attrs ['tau0']= self.t0
		f.attrs['analyzed_N'] = self.analyzed_N
		f.attrs ['fid0']=self.fid0
		f.attrs ['fid1']=self.fid1
		f.attrs ['T2']=self.T2
		try:
			f.attrs ['repetitions']=self.repetitions
		except:
			f.attrs ['repetitions']=-1
			print 'param Repetitions missing...'
		try:	
			f.attrs['CR_after_threshold'] = self.CR_after_threshold
		except:
			f.attrs['CR_after_threshold'] = -1
		f.attrs['reps'] = self.reps

		pr_grp = f.create_group('probability_densities')
		msqe_grp = f.create_group('mean_square_error')
		scaling = f.create_group('scaling')

		for i in self.prob_density_dict:
			pr_grp.create_dataset(i, data = self.prob_density_dict[i])

		for n in self.analyzed_N:
			n_grp = msqe_grp.create_group('N = '+str(n))
			n_grp.create_dataset ('B_field', data = self.results_dict[str(n)] ['B_field'])
			n_grp.create_dataset ('msqe', data = self.results_dict[str(n)] ['msqe'])
			n_grp.create_dataset ('ave_exp', data = self.results_dict[str(n)] ['ave_exp'])
			n_grp.create_dataset ('estimated_phase_values', data = self.results_dict[str(n)] ['estimated_phase_values'])
			try:
				n_grp.create_dataset ('nr_discarded_elements', data = self.results_dict[str(n)] ['nr_discarded_elements'])
			except:
				pass
			n_grp.attrs['N'] = n

		scaling.create_dataset ('total_time', data = self.total_time)
		scaling.create_dataset ('scaling_variance_phi', data = self.scaling_variance)
		scaling.create_dataset ('scaling_sensitivity', data = self.sensitivity)		
		scaling.attrs['scaling factor'] = self.scaling_factor
		scaling.attrs['error scaling factor'] = self.error_scaling_factor

		f.close()
		print 'Data saved!'

	def load_analysis (self, timestamp):

		name = toolbox.file_in_folder(folder=self.folder, timestamp = timestamp)
		print name
		f = h5py.File(os.path.join(self.folder, name),'r')

		self.F = f.attrs['F']
		self.G = f.attrs['G']
		self.K=f.attrs['K']
		self.t0 = f.attrs ['tau0']
		self.analyzed_N = f.attrs['analyzed_N']
		self.fid0 = f.attrs ['fid0']
		self.fid1 = f.attrs ['fid1']
		self.T2 = f.attrs ['T2']
		try:
			self.repetitions=f.attrs['repetitions']
		except:
			print 'Nr of repetitions: parameter not found!'
		pr_grp = f['/probability_densities']
		msqe_grp = f['/mean_square_error']

		for k in pr_grp.keys():
			self.prob_density_dict [k] = pr_grp[k].value

		for k in msqe_grp.keys():
			if type(msqe_grp[k])==h5py.Group:
				curr_subgrp = msqe_grp[k]
				curr_n = curr_subgrp.attrs['N']
				B_field = curr_subgrp['B_field'].value
				msqe = curr_subgrp['msqe'].value
				ave_exp = curr_subgrp['ave_exp'].value
				if self.error_bars:
					estimated_phase_values=curr_subgrp['estimated_phase_values'].value
				else:	
					estimated_phase_values=0*curr_subgrp['ave_exp'].value
				try:
					nr_discarded_elements=curr_subgrp['nr_discarded_elements'].value
				except:
					nr_discarded_elements=None	
				self.results_dict[str(curr_n)] =  {'B_field':B_field,'ave_exp':ave_exp, 'msqe':msqe, 'F':self.F,'G':self.G,'estimated_phase_values':estimated_phase_values,'nr_discarded_elements':nr_discarded_elements}
		#for k in scaling_grp.keys():
		#	self.scaling_grp[k]=scaling_grp[k].value
		f.close()


class magnetometrySQL(AdaptiveMagnetometry):

	def __init__(self, tau0, max_G):
		self.N = 1			#here it corresponds to the maximum N
		self.t0 = tau0
		self.reps = None
		self.B_max = 1./(4*self.t0) #this is different than the usual!
		self.n_points = round(100*max_G**0.5)
		#self.nr_B_points = 2**(self.N+2)/2
		self.results_dict = {} 
		self.prob_density_dict = {}
		self.analyzed_N = []
		self.gamma_e = 28e9 #~28 GHz/T

		self.nr_periods = None
		self.nr_points_per_period = None
		
		self.scaling_variance=[]
		self.total_time=[]
		self.sensitivity = []
		self.analyzed_N = []

		self.scaling_factor = None
		self.error_scaling_factor = None
		self.simulated_data = None
		self.error_bars = False

		self.always_recalculate_phase = False
		self.do_adaptive = False

		self.F = 0
		self.K = 0

		self.protocols = 'sql'
		self.do_sql = True

		if (os.name =='posix'):
			self.folder = '/home/cristian/Work/Research/adaptive magnetometry/data_analysis/'
		else:
			self.folder = r'M:/tnw/ns/qt/Diamond/Projects/Magnetometry with adaptive measurements/Data/analyzed data'


	def sweep_field_simulation_SQL (self, print_results=False, specific_B=False):
		self.simulated_data = True	
		B_values = np.array([])
		label_array = []
		N=1

		per=0
		delta_f = 1./(self.t0*(2**N))

		B = np.linspace(per*delta_f, (per+1)*delta_f, self.nr_points_per_period)
		if specific_B:
			B=self.specific_B_fields_for_sweep_sim(N,self.nr_points_per_period)
		self.B_values = np.hstack((B_values, B))
		for l in np.arange(self.nr_points_per_period):
			label_array.append('N='+str(N)+'G='+str(self.G)+'F='+str(self.F)+'_p'+str(0)+'_'+str(l))

		label_array, self.B_values = self.sample_B_space (N=N)		
		self.B_values=np.unique(self.B_values)*0.25
		nr_points = len(self.B_values)

		list_estim_phases = []
		G = 100*(np.arange(10)+1)
		msqe = np.zeros((len(G), nr_points))
		ave_exps = 0*msqe+1j*0*msqe
		B_field = np.zeros((len(G), nr_points))

		for i in np.arange(len(G)):
			list_estim_phases.append([])
			self.analyzed_N.append(G[i])	



		ind = 0
		for b in np.arange(nr_points):

			sys.stdout.write(str(ind)+', ')	
			s = RamseySequence_Simulation (N_msmnts = 1, reps=self.reps, tau0=self.t0)
			s.setup_simulation (magnetic_field_hz = self.B_values[b], F=0,G=self.G,K=0) 
			s.verbose = False
			s.T2 = self.T2
			s.fid0 = self.fid0
			s.fid1 = self.fid1
			s.B_max = 1/(4*s.t0)
			s.sim_SQL(phase_deg = 90)
			s.n_points = 10000
			g_idx = 0

			for g in G:
				print 'G  = ', g
				s.reshape_SQL(G=g)
				s.convert_to_dict()
				beta, p, ave_exp,H, mB, sB, list_phase_values = s.mean_square_error(set_value=self.B_values[b], do_plot=True, return_all_estimated_phases = True)
				list_estim_phases[g_idx].append(list_phase_values)
				ave_exps[g_idx, ind]=ave_exp
				msqe [g_idx, ind] = H
				B_field [g_idx, ind] = self.B_values[b]
				g_idx = g_idx +1
			ind =ind+1
		for i in np.arange(len(G)):
			list_estim_phases[i] = flatten(list_estim_phases[i])
			self.results_dict[str(G[i])] = {'B_field':B_field[i, :], 'ave_exp':ave_exps[i,:],'msqe':msqe[i,:], 'G':i,'K':self.K,'F':self.F, 'estimated_phase_values':list_estim_phases[i]}
		print 'Simulations done!!'

