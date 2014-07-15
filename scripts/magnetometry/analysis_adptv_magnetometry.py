

import numpy as np
from analysis.lib.fitting import fit, ramsey, common
from analysis.lib.tools import plot
import random
from matplotlib import rc, cm
import os, sys
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot
from matplotlib import rc, cm

reload(sequence)

class PhaseDistributionArray ():

	def __init__ (self, n_points, N_max, B_max, t0):
		self.n_points = n_points
		self.B_max = B_max
		self.N_max = N_max
		self.data = np.zeros ((self.N_max, self.n_points+1))
		self.beta = np.linspace (-self.B_max, self.B_max, self.n_points)
		self.curr = 0
		self.t0 = t0
		self.n = None
		self.mean = None
		self.var = None
		self.holevo_var = None
		
	def add (self, n, prob_distr):
		self.data[self.curr, 0] = n
		self.data[self.curr, 1:] = prob_distr
		self.curr=self.curr+1
		
	def statistics(self):
		m = []
		n = []
		v = []
		vH = []
		for i in np.arange(self.N_max):
			if (self.data[i,0] > 0):
				p = self.data[i, 1:]/np.sum(self.data[i,1:])
				m0 = np.sum(p*self.beta)
				m.append(m0)
				n.append(self.data[i, 0])
				v0 = np.sum(self.beta*self.beta*p)-m0**2
				v.append(v0)
				vH.append ((np.abs(np.sum(p*np.exp(1j*(2*np.pi*self.beta*self.t0)))))**(-2)-1)
		self.n = np.squeeze(np.array(n))
		self.mean = np.squeeze(np.array(m))
		self.var = np.squeeze(np.array(v))
		self.holevo_var = np.squeeze(np.array(vH))
		return n, m, v, vH
		
	def plot_distribution (self):
		ccc = cm.gist_heat(np.linspace(0., 0.8, self.N_max))
		for i in np.arange(self.N_max):
			if (self.data[i,0] > 0):
				plt.plot (self.beta*1e6, self.data[i, 1:], label =  str(self.data[i,0]), color = ccc[i])
		plt.legend()
		plt.show()
		
	def plot_scaling (self,M=None):
		self.statistics()
		if M==None:
			M=self.M
		self.T = M*self.t0*(2**(self.n+1)-1)
		plt.loglog (self.T, self.var*self.T, 'ob')
		plt.xlabel ('total estimation time')
		plt.ylabel ('variance*T ((B/g)**2*Hz)')
		plt.show()
		
		plt.loglog (self.T, self.holevo_var*self.T, 'ob')
		plt.xlabel ('total estimation time')
		plt.ylabel ('holevo_variance*T ((B/g [rad])**2*Hz)')
		plt.show()
		
	def save(self, name):
		np.savez (name, n_points= self.n_points, B_max = self.B_max, N_max = self.N_max, data = self.data,
				beta = self.beta, tau0 = self.t0, n = self.n, mean = self.mean, var = self.var, holevo_var = self.holevo_var)

	
		
				
	
class RamseySequence():

	def __init__ (self, N_msmnts, reps, tau0):
		self.N = N_msmnts
		self.points = 2**(self.N)+3
		self.discr_steps = 2*self.points+1
		self.p_k = np.zeros ((reps,self.discr_steps))+1j*np.zeros ((reps, self.discr_steps))
		self.msmnt_results = None
		self.msmnt_phases = None
		self.msmnt_times = None
		self.T2 = 4452e-9
		self.fid0 = 0.85
		self.fid1 = 0.015
		self.theta = 0*np.pi/180.
		self.t0 = tau0
		self.B_max = 1./(2*tau0)
		self.curr_rep = 0
		self.n_points = 50000
		self.curr_msmnt = 1
		self.reps = reps
		self.N_total = self.N
		self.majority_vote = False
		self.M=1
		self.M_eff = self.M
		
		for i in np.arange(reps):
			self.p_k[i, self.points] = 1/(2.*np.pi)
		
	def set_ideal (self):
		self.T2 = 1000.
		self.fid0 = 1.
		self.fid1 = 0.
	
	def field_values_MHz (self):
		return 1e-6*np.linspace (-self.B_max, self.B_max, self.n_points)
		
	def reset (self):
		self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
		self.p_k[self.points] = 1/(2.*np.pi)
		self.msmnt_results = None
		self.msmnt_phases = None
		self.msmnt_times = None
		self.curr_msmnt = 1
		
	def reset_rep_counter(self):
		self.curr_rep = 0
		self.curr_msmnt = 1


	def bayesian_update(self, m_n, phase_n, t_n,repetition = None):
			
		if (repetition == None):
			repetition = self.curr_rep
		p0 = np.copy (self.p_k[repetition, :])
		t_total = 2**(self.N)+1

		for k in np.arange(-t_total-1, t_total+1):
			if (k-t_n+self.points<0):
				p0_a = 0
			else:
				p0_a = p0 [k-t_n+self.points]
			
			if (k+t_n+self.points>self.discr_steps-1):
				p0_b = 0
			else:
				p0_b = p0 [k+t_n+self.points]
			q = m_n*np.pi+phase_n
			self.p_k [repetition, k+self.points] = 0.5*p0[k+self.points] + 0.25*(np.exp(1j*q)*p0_a + np.exp(-1j*q)*p0_b)
			norm = float(np.sum(np.abs(self.p_k[repetition, :])**2)**0.5)
			self.p_k[repetition, :] = self.p_k[repetition, :]/norm
		#plt.plot (np.abs(self.p_k[repetition,:]))
		#plt.show()
		#self.curr_msmnt = self.curr_msmnt+1		

	def inc_rep (self):
		self.curr_rep = self.curr_rep + 1
		self.curr_msmnt = 1
		if (self.curr_rep>self.reps):
			print 'Maximum repetition reached... Setting counter back to zero!'
			self.curr_rep = 0
			
	def phase_cappellaro_protocol (self):
		ttt = -2**(self.N-self.curr_msmnt+1)
		optimal_phase = -0.5*np.angle (self.p_k[self.curr_rep, ttt+self.points])
		return optimal_phase
		
	def holevo_variance (self):
		return (2*np.pi*self.p_k[-1+self.points])**(-2)-1.

	def print_msmnt(self):
		print 'Measurement results: ', self.msmnt_results
		print 'Phases: ', self.msmnt_phases
		print 'Ramsey time: ', self.msmnt_times*1e9, 'ns'
		

	def analysis (self, corrected=True, N_max = None, repetition = None):
		
		if (N_max==None):
			N_max = self.N
								
		if (repetition == None):
			repetition = self.curr_rep
			
		prob = np.ones(self.n_points)
		beta = np.linspace (-self.B_max, self.B_max, self.n_points)

		for n in np.arange(N_max*self.M_eff) +(self.N_total-N_max*self.M_eff):
			q = 2*np.pi*beta*self.msmnt_times[n]+self.msmnt_phases[repetition, n]
			dec = np.exp(-(self.msmnt_times[n]/self.T2))
			m = 1-self.msmnt_results[repetition, n]
			if corrected:
				prob = prob*(((1-m)*self.fid0+m)*(1-dec*np.cos(q)*(-1)**(m))+((m)*(1-self.fid0))*(1-((-1)**(1-m))*dec*np.cos(q)))
			else:
				prob = prob*(1-((-1)**m)*dec*np.cos(q))

		prob = prob/np.sum(np.abs(prob))
		return beta, prob

	def analysis_M (self, corrected=True, N_max = None, repetition = None):
		
		if (N_max==None):
			N_max = self.N
								
		if (repetition == None):
			repetition = self.curr_rep
			
		prob = np.ones(self.n_points)
		beta = np.linspace (-self.B_max, self.B_max, self.n_points)

		for n in np.arange(N_max) +(self.N-N_max):
			q = 2*np.pi*beta*self.msmnt_times[n]+self.msmnt_phases[repetition, n]
			dec = np.exp(-(self.msmnt_times[n]/self.T2))
			nr_ones = self.msmnt_results[repetition, n]
			nr_zeros = self.M-nr_ones

			#print 'Bayesian update: '+str(nr_zeros)+' zeros and '+str(nr_ones)+' ones.'  

			for j in np.arange(nr_ones):
				prob = prob*(1-dec*np.cos(q))
			for j in np.arange(nr_zeros):
				prob = prob*(1+dec*np.cos(q))

		prob = prob/np.sum(np.abs(prob))
		return beta, prob


		
	def plot_phase_distribution (self, repetition = 0):
		beta, prob = self.analysis_M (corrected=False, repetition = repetition)
		plt.plot (beta*1e-6, prob)
		
	def plot_avg_phase_distribution (self, max_rep = None):
		if (max_rep==None):
			max_rep = self.reps
		for j in np.arange(max_rep):
			if ( np.mod(10000*j/max_rep, 1000) == 0):
				print 100*j/max_rep, '%'
			beta, prob = self.analysis_M (corrected=False, repetition = j)
			if (j==0):
				prob_avg = prob
			else:
				prob_avg = prob_avg + prob
				
		plt.plot (beta*1e-6, prob_avg)
		plt.xlabel ('magnetic field [MHz]')
		plt.show()	
		return prob_avg
		

		
	def phase_distribution_scaling (self, max_rep = None, do_plot=False):
	
		print 'Analyzing scaling of estimation...'
		
		if (max_rep == None):
			max_rep = self.reps
		
		phase_distr = PhaseDistributionArray (N_max = self.N, n_points = self.n_points, B_max = self.B_max, t0=self.t0)

		for n in np.arange(self.N-1)+2:
			print '############', n
			for j in np.arange(max_rep):
				beta, prob = self.analysis_M (corrected = False, N_max = n, repetition = j)
				if (j==0):
					p_avg = prob
				else:
					p_avg = p_avg+prob
			p_avg = p_avg/np.sum(np.abs(p_avg))
			phase_distr.add (n = n, prob_distr = p_avg)
			mean = np.sum(p_avg*beta)
			print 'mean: ', mean*1e-6, 'MHz'
			variance = np.sum (p_avg*beta*beta)-mean**2
			print 'std: ', (variance**0.5)*1e-6, 'MHz'
			phi = 2*np.pi*beta*self.t0
			vH = (np.abs(np.sum(p_avg*np.exp(1j*phi))))**(-2)-1
			print np.abs(np.sum(p_avg*np.exp(1j*phi)))
			print 'holevo std: ', (vH/(2*np.pi*self.t0**2))**0.5
			if do_plot:
				plt.plot (beta*1e-6, p_avg)
				plt.xlabel ('magnetic field [MHz]')
				plt.show()
		
		return phase_distr
		
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

	def setup_simulation(self, magnetic_field_hz = 0., M = 1):
		self.beta_sim = magnetic_field_hz
		self.M = M

	def save_folder (self, folder = '/home/cristian/Work/Research/adaptive magnetometry/'):
		self.save_folder = folder
	
	def ramsey (self, t=0., theta=0.):
		n0 = 0.5*(1-np.exp(-(t/self.T2))*np.cos(2*np.pi*self.beta_sim*t+theta))
		n1 = 0.5*(1+np.exp(-(t/self.T2))*np.cos(2*np.pi*self.beta_sim*t+theta))
		pp0 = self.fid0*n0+self.fid1
		pp1 = n1+(1-self.fid0*n0)
		p0 = pp0/(pp0+pp1)
		p1 = pp1/(pp0+pp1)
		np.random.seed()
		result = 1-np.random.choice (2, 1, p=[p0, p1])
		return result[0]		
		

	def sim_non_adaptive_trivial (self, phase = [], tau = []):
		
		if (phase == []):
			phase = np.zeros(self.N)
		
		if (tau == []):
			nn = np.arange(self.N)+1
			tau = 2**(self.N-nn)
	
		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros((self.reps,self.N))
		self.msmnt_results = np.zeros((self.reps,self.N))
		
		for r in np.arange(self.reps):
			for i in np.arange(self.N):
				self.msmnt_phases[r, i] = phase[i]
				self.msmnt_times[i] = tau[i]*self.t0
				self.msmnt_results[r, i] = self.ramsey(t = tau[i]*self.t0, theta=phase[i])

	def sim_non_adaptive (self, K = 0, F=0, MK=0):
		
		total = 200#MK*(2**(K)-1)+F*(2**(K)-2-K)
		print 'total nr of msmnts: ', total
		
		self.msmnt_phases = np.zeros((self.reps, total))
		self.msmnt_times = np.zeros(total)
		self.msmnt_results = np.zeros((self.reps, total))
		
		for r in np.arange(self.reps):
			i = 0
			print 'Rep: ',r
			for l in np.arange(K+1):
				k = K-l
				self.msmnt_times[i] = (2**i)*self.t0
				phi = 0.
				M = MK+F*(K-k)	
				print 'M: ', M
				for m in np.arange (M)+1:
					self.msmnt_phases[r, i] = phi
					self.msmnt_results[r, i] = self.ramsey(t = (2**i)*self.t0, theta=phi)
					phi = np.mod(phi+np.pi/2, 2*np.pi)
					i = i+1
					
		
	def sim_cappellaro (self, debug = False, majority_vote = False):
		self.majority_vote = majority_vote
		
		print '----------------------------------'
		print 'Simulating Cappellaro protocol'
		print '----------------------------------'
		if majority_vote:
			self.N_total = self.N
			m_thr = 1+int(self.fid0*self.M/2.)
			self.M_eff = 1
		else:
			self.N_total = self.N*self.M
			self.M_eff = self.M
		self.msmnt_phases = np.zeros((self.reps,self.N_total))
		self.msmnt_times = np.zeros(self.N_total)
		self.msmnt_results = np.zeros((self.reps,self.N_total))
		
		nn = np.arange(self.N)+1
		tau = 2**(self.N-nn)
		
		self.reset_rep_counter()
		for r in np.arange(self.reps):
			print '- rep: ', self.curr_rep
			fff=0.
			ind = 0
			for i in np.arange(self.N):
				self.curr_msmnt = i
				
				total_counts = 0
				
				for m_step in np.arange(self.M):
					m = self.ramsey(t = tau[i]*self.t0, theta=fff)

				
					if debug:
						print '------ msmnt_step: (', self.curr_msmnt, ', ', m_step, ')'
						print '-- tau: ', tau[i], ' -- phase: ', fff*180/np.pi, ' --> ', m
					if (majority_vote == False):
						self.bayesian_update (m_n = m, phase_n = fff, t_n = tau[i])
						self.msmnt_results[r, ind] = m
						self.msmnt_phases[r, ind] = fff
						self.msmnt_times[ind] = tau[i]*self.t0
						ind = ind + 1
					else:
						total_counts = total_counts + m
				
				if majority_vote:
					if (total_counts>m_thr-1):
						m = 1
					else:
						m = 0
						
					if debug:
						print '------ Majority-vote nr. : ', self.curr_msmnt
						print '-- tau: ', tau[i], ' -- phase: ', fff*180/np.pi, ' - ', total_counts, ' --> ', m
					
					self.bayesian_update (m_n = m, phase_n = fff, t_n = tau[i])
					self.msmnt_results[r, i] = m
					self.msmnt_phases[r, i] = fff
					self.msmnt_times[i] = tau[i]*self.t0
					ind = ind + 1
					
				fff = self.phase_cappellaro_protocol()
				
				if debug:
					print '-- new phase: ', fff*180/np.pi
					

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
		self.n_points = 50000
		self.B_max = 1./(2*self.t0)
		
	def set_exp_pars(self, T2, fid0, fid1):
		self.T2=T2
		self.fid0 = fid0
		self.fid1 = fid1

	def load_data (self):

		a = sequence.SequenceAnalysis(self.folder)
		a.get_sweep_pts()
		a.get_magnetometry_data(name='adwindata', ssro = False)
		self.msmnt_results = a.clicks
		self.reps, self.N = np.shape (a.clicks)
		self.msmnt_times = np.zeros(len(a.ramsey_time))
		for j in np.arange(len(a.ramsey_time)):
			self.msmnt_times[j] = a.ramsey_time[j]
		self.msmnt_phases = 2*np.pi*a.set_phase/255.
		self.M=a.M
		print ' ---- phases_detuning' 
		print a.phases_detuning	
		print ' ---- adptv_phases:'
		print np.unique (self.msmnt_phases*180/np.pi)


		phases_detuning = 2*np.pi*a.phases_detuning/360.
		b = np.ones(self.reps)
		self.msmnt_phases = np.mod(self.msmnt_phases - np.outer (b, phases_detuning), 2*np.pi)

		print 'Some examples: '
		print self.msmnt_phases[0,:]*180/np.pi
		print self.msmnt_phases[1,:]*180/np.pi
		print self.msmnt_phases[10,:]*180/np.pi

		#print "Phases detuning", self.phases_detuning
					

	def old_real_space_analysis (self, log_plot=False, corrected = False, n_points=10000):
		self.reps = 30
		if (self.N_values == None):
			self.N_values = np.arange(self.N-1)+2
		
		ccc = cm.gist_heat(np.linspace(0., 0.8, len(self.N_values)))


		ind = 0
		for N in self.N_values:
			print 'Evaluating N = ', N, ' msmnts'
			p0 = ProbabilityDistribution (B_max=B_max, n_points=n_points)
			for j in np.arange(self.reps):
				r = AdaptiveRamsey(N_msmnts=N, B_max=B_max)
				r.msmnt_results = self.data [j,:N] 
				r.msmnt_phases = np.zeros(N)
				r.msmnt_times = self.ramsey_times [:N]

				p = r.analysis()
				p0.add(prob=p.data())
			p0.normalize()
			self.mean.append(p0.mean())
			self.variance.append((p0.variance()/B_max))
			plt.plot (p0.beta*1e-6, p0.data()/max(p0.data()), linewidth = 1, label = str(N), color = ccc[ind])
			ind = ind+1

		plt.xlabel ('magnetic field [MHz]')
		#plt.xlim([-3, 5])
		plt.legend()
		plt.show()
		
		total_time = 1e-9*(2**(self.N_values+1)-1)
		plt.plot (self.N_values, self.mean, 'ob')
		plt.show()

		plt.loglog (total_time*1e6, self.variance*total_time, 'ob')
		plt.xlabel ('total msmnt time [us]')
		plt.show()



def generate_table_old (N=10, debug=False):

	table = np.zeros(2**(N)-1)
	map_table = np.zeros(2**(N+1)-1)
	for n in np.arange(N)+1:
		start = 2**(n-1)
		tn = 2**(N-n)
		print '######## mmsnt step nr.', n, 'starting in array @ ', start
		print 'Ramsey time: ', tn
		
		#print (n-1)
		p_estim = AdaptiveRamsey (N_msmnts=n-1, B_max = B_max)
		for i in np.arange(2**(n-1)):
			msmnt_results = np.zeros(n-1)
			a = i
			for k in np.arange(n-1)+1:
				msmnt_results [n-1-k] = np.mod(a,2)
				a = a/2
			#print start+i, msmnt_results
			fff = 0.
			for k in np.arange(n-1):

				p_estim.bayesian_update (msmnt_result = msmnt_results[k], phase=fff, tn=tn) 
				fff= p_estim.cappellaro_protocol()
				print p_estim.curr_msmnt, '-', fff*180/np.pi
			opt_phase = fff
			table[start+i]=(np.mod(opt_phase*180/np.pi, 360)) 
			map_table[start+i] = i
			print n, i, msmnt_results, table[start+i]
			p_estim.reset()		

			
	#return map_table, table

def generate_table (N=10):
	table = np.zeros(2**N+2)
	points = 2**(N+3)
	discr_steps = 2*points+1
	
	x = np.linspace (0, 2*np.pi, discr_steps)
	m = np.zeros (N+1)
	t = np.zeros (N+1)
	p = np.zeros ((N+1, discr_steps))+1j*np.zeros ((N+1, discr_steps))
	p [0, points] = 1/(2.*np.pi)
	t_total = 2**(N+2)+1

	for j in np.arange (2**(N-1)):

		th = np.zeros(N+1)
		p = np.zeros ((N+1, discr_steps))+1j*np.zeros ((N+1, discr_steps))
		p [0, points] = 1/(2.*np.pi)
		msmnt_results = np.zeros(N)
		a = j
		for k in np.arange(N-1)+1:
			msmnt_results [N-k-1] = np.mod(a,2)
			a = a/2
		print '#############'
		print 'msmnt results: ', msmnt_results
		
		for n in np.arange(N)+1:
	
			t[n] = 2**(N-n)
			ttt = -2**(N-n+1)
			th[n] = -0.5*np.angle (p[n-1, ttt+points])
			m[n] = msmnt_results[n-1]
			pp0 = 2**(n-1)
			pp1 = 0
			for i in np.arange (n-1)+1:
				pp1 = pp1+msmnt_results[i-1]*2**(i-1)
			pos0 = pp0+pp1
			#update rule:
			for k in np.arange(-t_total-1, t_total+1):
				p [n, k+points] = 0.5*p[n-1, k+points] + 0.25*np.exp(-(t[n]/4400.))*(np.exp(1j*(m[n]*np.pi+th[n]))*p [n-1, k-t[n]+points] + 
						np.exp(-1j*(m[n]*np.pi+th[n]))*p [n-1, k+t[n]+points])

			p [n, :] = p[n, :]/np.sum(np.abs(p[n,:]))
			fff = -0.5*np.angle (p[n, ttt+points])	
			#print int(pos0), t[n], msmnt_results[:n-1], fff*180/np.pi
			table[int(pos0)] = fff*180/np.pi
		#print j, msmnt_results, np.mod(th*180/np.pi, 360)
			
	return table[1:]

	
def generate_M_ary_table (N=7, M = 2):
	
	emme = M+1
	total_length = int(((M+1)**(N)-1)/M)+1
	print 'Total length array: ', total_length
	table = np.zeros(total_length)
	points = 2**(N+3)
	discr_steps = 2*points+1
	
	m = np.zeros (N*M+1)
	t = np.zeros (N*M+1)
	#p = np.zeros ((N*M+1, discr_steps))+1j*np.zeros ((N*M+1, discr_steps))
	#p [0, points] = 1/(2.*np.pi)
	t_total = 2**(N+2)+1

	for j in np.arange ((emme)**(N-1)):

		th = np.zeros(N+1)
		p = np.zeros (discr_steps)+1j*np.zeros (discr_steps)
		p [points] = 1/(2.*np.pi)
		p0 = np.copy (p)
		msmnt_results = np.zeros(N)
		a = j
		for k in np.arange(N-1)+1:
			msmnt_results [N-k-1] = np.mod(a,emme)
			a = a/emme
		print '#############'
		print 'msmnt results: ', msmnt_results
		

		for n in np.arange(N)+1:
	
			t[n] = 2**(N-n)
			ttt = -2**(N-n+1)
			th[n] = -0.5*np.angle (p[ttt+points])
			
			n_ones = msmnt_results[n-1]
			n_zeros = M-n_ones

			pp0 = ((M+1)**(n-1)-1)/M+1
			pp1 = 0
			for i in np.arange (n-1)+1:
				pp1 = pp1+msmnt_results[i]*(emme)**(i-1)
			pos0 = pp0+pp1

			#update rule (zeros):
			for ss in np.arange(n_zeros):
				for k in np.arange(-t_total-1, t_total+1):
					p [k+points] = 0.5*p0[k+points] + 0.25*np.exp(-(t[n]/4400.))*(np.exp(1j*(th[n]))*p0 [k-t[n]+points] + 
							np.exp(-1j*(th[n]))*p0 [k+t[n]+points])
				p = p/np.sum(np.abs(p))
				p0 = np.copy (p)

			#update rule (ones):
			for ss in np.arange(n_ones):
				for k in np.arange(-t_total-1, t_total+1):
					p [k+points] = 0.5*p0[k+points] + 0.25*np.exp(-(t[n]/4400.))*(np.exp(1j*(np.pi+th[n]))*p0 [k-t[n]+points] + 
							np.exp(-1j*(np.pi+th[n]))*p0 [k+t[n]+points])
				p = p/np.sum(np.abs(p))
				p0 = np.copy (p)
			
			fff = -0.5*np.angle (p[ttt+points])	
			#print int(pos0), t[n], msmnt_results[:n-1], fff*180/np.pi
			table[int(pos0)] = fff*180/np.pi
		#print j, msmnt_results, np.mod(th*180/np.pi, 360)

	return table[1:]


f = toolbox.latest_data(contains='220832')
s = RamseySequence_Exp (folder = f)
s.set_exp_pars (T2=4000e-9, fid0=0.85, fid1=0.015)
#s.M=1
#s.set_ideal()
s.load_data()
p = s.plot_avg_phase_distribution()
#s.print_phases()
p = s.phase_distribution_scaling (do_plot=True)
p.plot_scaling()

