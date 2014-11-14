

import numpy as np
from analysis.lib.fitting import fit, ramsey, common
from analysis.lib.tools import plot
import random
from matplotlib import rc, cm
import os, sys
import h5py
import logging, time, timeit

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import  sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr, common
from analysis.lib.tools import plot
from analysis.lib.tools import compare_functions as compare
from analysis.lib.m2 import m2
from matplotlib import rc, cm
from analysis.lib.magnetometry import adaptive_magnetometry as adptv_mgnt


#reload(sequence)
reload(compare)
reload(toolbox)
reload(adptv_mgnt)


class RamseySequence_Adwin (adptv_mgnt.RamseySequence_Simulation):


	def f_pr (self, value):

		if (np.real(value)==0):
			return '{0:.5f}'.format(value)
		else:
			if (float(np.imag(value))>0):
				s = '+'
			else:
				s = '-'
			return '{0:.5f}'.format(float(np.real(value)))+s+'{0:.5f}'.format(abs(float(np.imag(value))))+'j'



	def sim_cappellaro_majority_with_plots (self, N_sim = []):
		
		if self.verbose:				
			print '-----------------------------------------'
			print 'Simulating Cappellaro protocol (maj vote)'
			print '-----------------------------------------'
			print '- N = '+str(self.N)+ ', M = '+str(self.M)

		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))

		actual_points = self.points  
		
		self.points = 2**(self.N+2)+1
		self.discr_steps = 2*self.points+1
		
		nn = np.arange(self.N)+1
		tau = 2**(self.N-nn)
		self.msmnt_times = tau
		self.reset_rep_counter()

		msmnt_results = np.zeros (self.N)
		t = np.zeros (self.N)
		phase = np.zeros(self.N)
		self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
		self.p_k [self.points] = 1/(2.*np.pi)

		x_axis = np.arange(self.discr_steps)-self.points

		for n in N_sim:

			t[n-1] = int(2**(self.N-n))
			ttt = -2**(self.N-n+1)
			phase[n-1] = 0.5*np.angle (self.p_k[ttt+self.points])

			m_total = 0
			print '### N  = ', n
			for m in np.arange(self.M):
				m_res = self.majority_vote_msmnt (theta_n = phase[n-1], t_n = t[n-1])				
				self.bayesian_update (m_n = m_res, phase_n = phase[n-1], t_n = 2**(self.N-n))
				plt.semilogy (x_axis, np.real(self.p_k), 'o', label = 'M = '+str(m+1))
				m_total = m_total + m_res
				print '--- m = ', m, '. Non-zero elements:'
				print ' Re[p_k], k = ', np.nonzero(np.real(self.p_k))-(self.points)*np.ones(len(np.nonzero(np.real(self.p_k))))

			msmnt_results[n-1] =m_total
			plt.title ('N = '+str(self.N)+', n = '+str(n), fontsize=15)
			plt.legend()
			plt.show()
				
		self.points = actual_points
		self.discr_steps = 2*self.points+1
		self.msmnt_results [0, :] = np.copy(msmnt_results)
		self.msmnt_phases [0, :] = np.copy(phase)
		self.inc_rep()


	def sim_cappellaro_track_coefficients (self, msmnt_result = []):
		

		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))

		actual_points = self.points  
		
		self.points = 2**(self.N+2)+1
		self.discr_steps = 2*self.points+1
		
		nn = np.arange(self.N)+1
		tau = 2**(self.N-nn)
		self.msmnt_times = tau
		self.reset_rep_counter()

		msmnt_results = np.zeros (self.N)
		t = np.zeros (self.N)
		phase = np.zeros(self.N)
		self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
		self.p_k [self.points] = 1/(2.*np.pi)

		pi_tn = np.zeros(self.N)+1j*np.zeros(self.N)
		pi_2tn = np.zeros(self.N)+1j*np.zeros(self.N)

		for n in np.arange(self.N)+1:

			t[n-1] = int(2**(self.N-n))
			ttt = -2**(self.N-n+1)
			phase[n-1] = 0.5*np.angle (self.p_k[ttt+self.points])
			pi_tn[n-1] = self.p_k[t[n-1]+self.points]
			pi_2tn[n-1] = self.p_k[-ttt+self.points]

			m_total = 0
			nr_ones = msmnt_result[n-1]
			nr_zeros = self.M-nr_ones
			for m in np.arange(nr_ones):
				self.bayesian_update (m_n = 1, phase_n = phase[n-1], t_n = 2**(self.N-n))
			for m in np.arange(nr_zeros):
				self.bayesian_update (m_n = 0, phase_n = phase[n-1], t_n = 2**(self.N-n))
			
		self.points = actual_points
		self.discr_steps = 2*self.points+1
		self.msmnt_results [0, :] = np.copy(msmnt_results)
		self.msmnt_phases [0, :] = np.copy(phase)
		return phase, pi_tn, pi_2tn

	def basic_adwin_algorithm (self, debug = False, exec_speed=False):

		discr_steps = 2**(self.N+2)+1
		m = np.zeros (self.N+1)
		t = np.zeros (self.N+1)
		th = np.zeros(self.N+1)
		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))
		exec_time = None
		
		if exec_speed:
			start = time.clock()

		for rep in np.arange(self.reps):
			p_real = np.zeros (discr_steps)
			p_imag = np.zeros (discr_steps)
			p_real [2**(self.N+1)] = 1./(2*np.pi)

			curr_phase = 0
			for n in np.arange(self.N)+1:
				t[n] = 2**(self.N-n)
				k_opt = -2**(self.N-n+1)+2**(self.N+1)
				th[n] = -0.5*np.angle(1j*p_imag[k_opt]+p_real[k_opt])
				m[n] = self.ramsey (theta=th[n], t = t[n]*self.t0)

				if debug:
					print '### n =', n
					print 't = ', t[n], '  theta = ', th[n]*180./np.pi, '  --- res: ', m[n] 
				#update rule:
				cn = m[n]*np.pi+th[n]
				p0_real = np.copy (p_real)
				p0_imag = np.copy (p_imag)

				for k in np.arange(2**self.N, 3*2**self.N):
					p_real [k] = 0.5*p0_real[k] + 0.25*(np.cos(cn)*(p0_real [k-t[n]] + p0_real [k+t[n]]) - np.sin(cn)*(p0_imag [k-t[n]] - p0_imag [k+t[n]])) 
					p_imag [k] = 0.5*p0_imag[k] + 0.25*(np.cos(cn)*(p0_imag [k-t[n]] + p0_imag [k+t[n]]) + np.sin(cn)*(p0_real [k-t[n]] - p0_real [k+t[n]])) 

			self.msmnt_results [rep, :] = m[1:]
			self.msmnt_phases [rep, :] = th[1:]
			self.msmnt_times = t [1:]
		if exec_speed:
			stop = time.clock()
			exec_time = stop-start
		return exec_time


	def adwin_algorithm (self, outcomes=[], do_plot=False, do_print=True):

		discr_steps = 2**(self.N+2)+1
		m = np.zeros (self.N+1)
		t = np.zeros (self.N+1)
		th = np.zeros(self.N+1)
		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))
		p_tn = np.zeros(self.N+1)+1j*np.zeros(self.N+1)
		p_2tn = np.zeros(self.N+1)+1j*np.zeros(self.N+1)


		for rep in np.arange(self.reps):
			p_real = np.zeros (discr_steps+1)
			p_imag = np.zeros (discr_steps+1)
			p_real [1] = self.N*self.M**2

			curr_phase = 0

			if do_plot:
				f, axarr = plt.subplots(2, sharex=True, figsize=(10,10))
				x = np.arange(2**self.N+1)
			for n in np.arange(self.N)+1:
				t[n] = 2**(self.N-n)
				k_opt = 2*t[n]
				th[n] = -0.5*np.angle(-1j*p_imag[1+k_opt]+p_real[1+k_opt])
				nr_ones = outcomes[n-1]
				nr_zeros = self.M-nr_ones

				if do_plot:
					axarr[0].plot (p_real[1:], ':k')
					axarr[1].plot (p_imag[1:], ':k')
					axarr[0].plot (p_real[1:], 'o', label=str(n))
					axarr[1].plot (p_imag[1:], 'o', label=str(n))

				p_tn [n] = p_real[1+t[n]]+1j*p_imag[1+t[n]]
				p_2tn [n] = p_real[1+2*t[n]]+1j*p_imag[1+2*t[n]]

				for mmm in np.arange(nr_ones):
					p0_real = np.copy (p_real)
					p0_imag = np.copy (p_imag)
					cn = np.pi - th[n]
					for k in np.arange(0, 2**(self.N+1)+1):
						if (k<t[n]):
							p_im_min = -p0_imag[1+np.abs(k-t[n])]
						else:
							p_im_min = p0_imag[1+np.abs(k-t[n])]
						p_real [1+k] = 0.5*p0_real[1+k] + 0.25*(np.cos(cn)*(p0_real [1+np.abs(k-t[n])] + p0_real [1+k+t[n]]) - np.sin(cn)*(p_im_min - p0_imag [1+k+t[n]])) 
						p_imag [1+k] = 0.5*p0_imag[1+k] + 0.25*(np.cos(cn)*(p_im_min + p0_imag [1+k+t[n]]) + np.sin(cn)*(p0_real [1+np.abs(k-t[n])] - p0_real [1+k+t[n]])) 

					norm = p_real[1]*2*np.pi
					p_real = p_real/norm
					p_imag = p_imag/norm

				for mmm in np.arange(nr_zeros):
					p0_real = np.copy (p_real)
					p0_imag = np.copy (p_imag)
					cn = 0*np.pi - th[n]
					for k in np.arange(0, 2**(self.N+1)+1):
						if (k<t[n]):
							p_im_min = -p0_imag[1+np.abs(k-t[n])]
						else:
							p_im_min = p0_imag[1+np.abs(k-t[n])]
						p_real [1+k] = 0.5*p0_real[1+k] + 0.25*(np.cos(cn)*(p0_real [1+np.abs(k-t[n])] + p0_real [1+k+t[n]]) - np.sin(cn)*(p_im_min - p0_imag [1+k+t[n]])) 
						p_imag [1+k] = 0.5*p0_imag[1+k] + 0.25*(np.cos(cn)*(p_im_min + p0_imag [1+k+t[n]]) + np.sin(cn)*(p0_real [1+np.abs(k-t[n])] - p0_real [1+k+t[n]])) 

					norm = p_real[1]*2*np.pi
					p_real = p_real/norm
					p_imag = p_imag/norm

			if do_plot:
				axarr[0].set_title('real part')
				axarr[1].set_title('imaginary part')
				axarr[0].set_xlim([0, 2**self.N+1])
				axarr[1].set_xlim([0, 2**self.N+1])
				plt.legend()
				plt.show()
			self.msmnt_results [rep, :] = np.array(outcomes)
			self.msmnt_phases [rep, :] = th[1:]
			self.msmnt_times = t [1:]

			if do_print:
				print 'Measurement outcomes: ', outcomes
				print 'Optimal phases: ', np.round(th[1:]*180/np.pi)

				for n in np.arange(self.N)+1:
					print '*step-'+str(n)+' p[tn] = ', self.f_pr(p_tn[n]), '   p[2tn] = ', self.f_pr(p_2tn[n])


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
		norm = np.real(p[self.points])*2*np.pi
		p = p/norm
		self.p_k = np.copy (p)




	def compare_algorithms (self, outcomes = [], do_plot=True, do_print = False, do_plot_adwin = True):
		
		#compare results of "adwin"-aimed protocol to full-scale python simulations

		self.points = 2**(self.N+6)+1
		self.discr_steps = 2*self.points+1

		
		self.p_k = np.zeros (self.discr_steps)+1j*np.zeros (self.discr_steps)
		self.p_k [self.points] = 1/(2.*np.pi)

		pi_tn = np.zeros(self.N+1)+1j*np.zeros(self.N+1)
		pi_2tn = np.zeros(self.N+1)+1j*np.zeros(self.N+1)
		t = np.zeros(self.N+1)
		phase_python = np.zeros(self.N+1)
		phase_adwin = np.zeros(self.N+1)
		range_k_sweep = 2**(self.N+4)+1
		p_real = np.zeros (range_k_sweep+1)
		p_imag = np.zeros (range_k_sweep+1)
		p_real [1] = 1/(2*np.pi)
		p_2tn_adwin = np.zeros(self.N+1)+1j*np.zeros(self.N+1)
		p_2tn_python = np.zeros(self.N+1)+1j*np.zeros(self.N+1)

		diff = None
		for n in np.arange(self.N)+1:

			t[n] = int(2**(self.N-n))
			k_opt = -2*t[n]
			phase_python[n] = 0.5*np.angle (self.p_k[k_opt+self.points])
			phase_adwin[n] = 0.5*np.angle(-1j*p_imag[int(1+2*t[n])]+p_real[int(1+2*t[n])])

			'''
			if np.abs(self.p_k[k_opt+self.points])<1e-7:
				phase_python[n] = 0
			
			if np.abs(-1j*p_imag[int(1+2*t[n])]+p_real[int(1+2*t[n])])<1e-7:
				phase_adwin[n] = 0
			'''
			if do_print:
				print ' #### n = ', n
				print 'Optimal phase els: pyth = ', k_opt+self.points, '(',k_opt,')) - adwin = ', int(1+2*t[n])
				print 'Values : ', self.p_k[k_opt+self.points], ' --- ', -1j*p_imag[int(1+2*t[n])]+p_real[int(1+2*t[n])]
				print 'Phases : ', phase_python[n], ' --- ', phase_adwin[n]

			p_2tn_adwin[n] = p_real[1+2*t[n]]+1j*p_imag[1+2*t[n]]
			p_2tn_python[n] = self.p_k[2*t[n]+self.points]


			#print '--- p[2*tn]:  adwin = ', self.f_pr(p_2tn_adwin[n]), '    python = ', self.f_pr(p_2tn_python[n])


			nr_ones = outcomes[n-1]
			nr_zeros = self.M-nr_ones

			max_k = int(t[n])*(self.M**2+10)
			if (max_k>range_k_sweep):
				nr_k  = int(range_k_sweep)/int(t[n])
				max_k = int(t[n]+1)*nr_k
			k_space = np.arange(0, max_k, int(t[n]))
			print 'max_k: ', max_k, '     --- range_k_sweep: ', range_k_sweep


			#print 'k_space: ', k_space

			for m in np.arange(nr_ones):
				self.bayesian_update (m_n = 1, phase_n = phase_python[n], t_n = int(t[n]))

				p0_real = np.copy (p_real)
				p0_imag = np.copy (p_imag)
				cn = np.pi + phase_adwin[n]

	

				for k in k_space:
					if (k+t[n]>range_k_sweep):
						p0r_pl = 0
						p0i_pl = 0
					else:
						p0r_pl = p0_real [1+k+t[n]]
						p0i_pl = p0_imag [1+k+t[n]]
					if (k<t[n]):
						p_im_min = -p0_imag[1+np.abs(k-t[n])]
					else:
						p_im_min = p0_imag[1+np.abs(k-t[n])]
					p_real [1+k] = 0.5*p0_real[1+k] + 0.25*(np.cos(cn)*(p0_real [1+np.abs(k-t[n])] + p0r_pl) + np.sin(cn)*(p_im_min - p0i_pl)) 
					p_imag [1+k] = 0.5*p0_imag[1+k] + 0.25*(np.cos(cn)*(p_im_min + p0i_pl) - np.sin(cn)*(p0_real [1+np.abs(k-t[n])] - p0r_pl)) 

				#norm = p_real[1]*2*np.pi
				#p_real = p_real/norm
				#p_imag = p_imag/norm

			for m in np.arange(nr_zeros):
				self.bayesian_update (m_n = 0, phase_n = phase_python[n], t_n = int(t[n]))

				p0_real = np.copy (p_real)
				p0_imag = np.copy (p_imag)
				cn = 0*np.pi + phase_adwin[n]
				for k in k_space:
					if (k+t[n]>range_k_sweep):
						p0r_pl = 0
						p0i_pl = 0
					else:
						p0r_pl = p0_real [1+k+t[n]]
						p0i_pl = p0_imag [1+k+t[n]]
					if (k<t[n]):
						p_im_min = -p0_imag[1+np.abs(k-t[n])]
					else:
						p_im_min = p0_imag[1+np.abs(k-t[n])]
					p_real [1+k] = 0.5*p0_real[1+k] + 0.25*(np.cos(cn)*(p0_real [1+np.abs(k-t[n])] + p0r_pl) + np.sin(cn)*(p_im_min - p0i_pl)) 
					p_imag [1+k] = 0.5*p0_imag[1+k] + 0.25*(np.cos(cn)*(p_im_min + p0i_pl) - np.sin(cn)*(p0_real [1+np.abs(k-t[n])] - p0r_pl)) 

				#norm = p_real[1]*2*np.pi
				#p_real = p_real/norm
				#p_imag = p_imag/norm

			if do_plot:
				len_array = len(p_real)
				th_adw = int(round(phase_adwin[n]*180/np.pi))
				th_pyt = int(round(phase_python[n]*180/np.pi))

				print '######  n = ', n
				print '   th_adw = ', th_adw, '  --- th_pyth = ', th_pyt
				print '	  msmnt outcome: ', nr_ones
				diff = (((p_real[1:]-np.real(self.p_k[self.points:self.points+len_array-1])))**2+((p_imag[1:]-np.imag(self.p_k[self.points:self.points+len_array-1])))**2)**0.5
				avg_diff = np.sum(diff)/len(diff)
				f, axarr = plt.subplots(3, sharex=True, figsize=(10,10))
				x = np.arange(2**self.N+1)
				axarr[0].plot (p_real[1:], ':k')
				axarr[0].plot (np.real(self.p_k[self.points:]), ':b')
				axarr[1].plot (p_imag[1:], ':k')
				axarr[1].plot (np.imag(self.p_k[self.points:]), ':b')
				axarr[0].plot (p_real[1:], 'o', label='adw_'+str(n))
				axarr[1].plot (p_imag[1:], 'o', label='adw_'+str(n))
				axarr[0].plot (np.real(self.p_k[self.points:]), 'o', label='pyth_'+str(n))
				axarr[1].plot (np.imag(self.p_k[self.points:]), 'o', label='pyth_'+str(n))
				axarr[2].semilogy (diff, ':k')
				axarr[2].semilogy (diff, 'ok')
				axarr[0].set_title('real part, n = '+str(n))
				axarr[1].set_title('imaginary part, n = '+str(n))
				axarr[2].set_title('difference, n = '+str(n)+', --- avg_diff = '+str(avg_diff)			)
				axarr[0].set_xlim([0, 2**self.N+1])
				axarr[1].set_xlim([0, 2**self.N+1])
				axarr[2].set_xlim([0, 2**self.N+1])
				axarr[0].set_ylim([-0.2, 0.2])
				axarr[1].set_ylim([-0.2, 0.2])
				axarr[0].legend()
				axarr[1].legend()
				axarr[2].legend()
				plt.show()

			if do_plot_adwin:

				plt.figure()
				plt.plot (p_real[1:], ':b')
				plt.plot (p_real[1:], 'ob', label='real_part')
				plt.plot (p_imag[1:], ':r')
				plt.plot (p_imag[1:], 'or', label='imag_part')
				plt.xlabel ('k', fontsize=12)
				plt.xlim([0, 16])
				plt.legend()
				plt.show()
		
		return phase_adwin, phase_python, diff, p_2tn_adwin, p_2tn_python

	def compare_adwin_python (self, do_plot=False, newer_than = None, use_fid_bayesian_update = False):
		
		#compares {p_k} distributions obtaind by the adwin with 
		#the ones calculatd in python

		range_k_sweep = 2**self.N*(self.M+3)
		p_real = np.zeros (range_k_sweep*2+1)
		p_imag = np.zeros (range_k_sweep*2+1)
		p_real [1] = 1/(2*np.pi)

		t = np.zeros(self.N+1)
		outcomes = np.zeros(self.N+1)
		phase_adwin = np.zeros(self.N+1)
		phase_simulation = np.zeros(self.N+1)
		outcomes = np.zeros(self.N+1)
		ext_outcomes = self.outcomes

		label = ''
		for i in np.arange(self.N):
			label=label+str(ext_outcomes[i])
		print 'Label: ', label

		for n in np.arange(self.N)+1:
			t[n] = int(2**(self.N-n))
			k_opt = -2*t[n]
			phase_simulation[n] = 0.5*np.angle(-1j*p_imag[1+2*t[n]]+p_real[1+2*t[n]])
			print 'n = ', n

			max_k = int(t[n])*(self.M+3)
			if (max_k>range_k_sweep):
				nr_k  = int(range_k_sweep)/int(t[n])
				max_k = int(t[n])*nr_k
			k_space = np.arange(0, max_k, int(t[n]))
			print 'max_k: ', max_k, '     --- range_k_sweep: ', range_k_sweep

			for m in np.arange(self.M)+1:				

				exp_exists = True
				try:
					f_name = 'N='+str(self.N)+'_M=('+str(self.M)+', 0)_rtAdwin_'+label+'_test_pk_(n='+str(n)+'_m='+str(m)+')'
					f = toolbox.latest_data (contains = f_name, return_all=False, newer_than = newer_than)
					exp = adptv_mgnt.RamseySequence_Exp (folder = f)
					exp.load_exp_data()
					phase_adwin = exp.opt_phase[:self.N]
					print 'Phase adwin (exp): ', phase_adwin
					self.fid0 = exp.exp_fid0
					self.fid1 = exp.exp_fid1
					#self.T2 = exp.T2
					exp_exists = True
				except:
					exp_exists = False
					print f_name+' does not exist!'


				if (ext_outcomes == []):
					m_res = self.ramsey (t = t[n]*self.t0, theta = phase_simulation[n])
				else:
					ext_outcomes[n-1] = ext_outcomes[n-1] - 1
					if (ext_outcomes[n-1]<0):
						m_res = 0
					else:
						m_res = 1

				outcomes[n] = outcomes[n] + m_res

				p0_real = np.copy (p_real)
				p0_imag = np.copy (p_imag)
				cn = m_res*np.pi + phase_simulation[n]


				for k in k_space:
					if (k+t[n]>max(k_space)):
						p0r_pl = 0
						p0i_pl = 0
					else:
						p0r_pl = p0_real [1+k+t[n]]
						p0i_pl = p0_imag [1+k+t[n]]
					if (k<t[n]):
						p_im_min = -p0_imag[1+np.abs(k-t[n])]
					else:
						p_im_min = p0_imag[1+np.abs(k-t[n])]

					if (int(m_res) == 0):
						A = 1-0.5*(self.fid0+self.fid1)
					elif (int(m_res) == 1):
						A = 0.5*(self.fid0+self.fid1)
					B = 0.5*(self.fid1-self.fid0)*np.exp(-(t[n]*self.t0/self.T2)**2)
					p_real [1+k] = A*p0_real[1+k] - 0.5*B*np.cos(cn)*(p0_real [1+np.abs(k-t[n])] + p0r_pl) - 0.5*B* np.sin(cn)*(p_im_min - p0i_pl)
					p_imag [1+k] = A*p0_imag[1+k] - 0.5*B*np.cos(cn)*(p_im_min + p0i_pl) +0.5*B*np.sin(cn)*(p0_real [1+np.abs(k-t[n])] - p0r_pl)

				norm = p_real[1]*2*np.pi
				p_real = p_real/norm
				p_imag = p_imag/norm


				if exp_exists:
					len_array = len(p_real)
					th_sim = int(round(phase_simulation[n]*180/np.pi))
					th_adw = phase_adwin[n-1]

					print '######  n = ', n, ' ---- m = ', m
					print '   phase_th = ', th_sim, '  --- phase_exp = ', th_adw
					print 'cos(cn) = ', np.cos(cn), '    sin(cn) = ', np.sin(cn)

					if do_plot:
						diff = (((p_real[1:2**self.N+2]-exp.real_pk_adwin[:2**self.N+1]))**2+((p_imag[1:2**self.N+2]-exp.imag_pk_adwin[:2**self.N+1]))**2)**0.5
						rel_diff = 100*diff/((p_real[1:2**self.N+2]**2+p_imag[1:2**self.N+2]**2)**0.5)
						avg_diff = np.sum(diff)/len(diff)
						f, axarr = plt.subplots(3, sharex=True, figsize=(10,10))
						x = np.arange(2**self.N+1)

						axarr[0].plot (p_real[1:], ':k')
						axarr[0].plot (exp.real_pk_adwin, ':b')
						axarr[1].plot (p_imag[1:], ':k')
						axarr[1].plot (exp.imag_pk_adwin, ':b')
						axarr[0].plot (p_real[1:], 'o', label='sim_'+str(n)+','+str(m))
						axarr[1].plot (p_imag[1:], 'o', label='sim_'+str(n)+','+str(m))
						axarr[0].plot (exp.real_pk_adwin, 'o', label='exp_'+str(n)+','+str(m))
						axarr[1].plot (exp.imag_pk_adwin, 'o', label='exp_'+str(n)+','+str(m))
						axarr[2].plot (rel_diff, ':k')
						axarr[2].plot (rel_diff, 'ok')
						axarr[0].set_title('real part, n = '+str(n)+', m = '+str(m))
						axarr[1].set_title('imaginary part, n = '+str(n)+', m = '+ str(m))
						axarr[2].set_title('difference, n = '+str(n)+', m = '+ str(m))
						axarr[0].set_xlim([0, self.M*2**(self.N-1)+1])
						axarr[1].set_xlim([0, self.M*2**(self.N-1)+1])
						axarr[2].set_xlim([0, self.M*2**(self.N-1)+1])
						axarr[0].set_ylim([-0.2, 0.2])
						axarr[1].set_ylim([-0.2, 0.2])
						axarr[2].set_ylim([0, 100])					
						axarr[2].set_ylabel('rel diff %', fontsize=15)
						axarr[0].legend()
						axarr[1].legend()
						axarr[2].legend()
						time_label = time.strftime ('%Y%m%d_%H%M%S')
						fname = 'M:/tnw/ns/qt/Diamond/Projects/Magnetometry with adaptive measurements/Data/adwin_tests/'+time_label+'_N='+str(self.N)+'_M='+str(self.M)+'_pk_'+str(n)+'_'+str(m)+'.png'
						print fname
						f.savefig (fname)
						plt.show()
		print 'Phase simulation: ', phase_simulation[1:]*180/np.pi


	def adwin_optimal_looping (self, use_fid_bayesian_update = True, ext_outcomes = []):
		
		#optimal-looping code (designed for adwin)

		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))

		self.reset_rep_counter()
		range_k_sweep = 2**self.N*(self.M/2+1)#2**(self.N+4)+1
		#print 'range_k_sweep: ', range_k_sweep,  ' -- ', 2**(self.N+4)+1
		print 'range_k_sweep: ', range_k_sweep

		if (ext_outcomes != []):
			self.reps = 1

		for r in np.arange(self.reps):

			if np.mod(r, 100)==0:
				print str(r)+'/'+str(self.reps)
			p_real = np.zeros (range_k_sweep)
			p_imag = np.zeros (range_k_sweep)
			p_real [1] = 1/(2*np.pi)

			t = np.zeros(self.N+1)
			outcomes = np.zeros(self.N+1)
			phase_adwin = np.zeros(self.N+1)

			for n in np.arange(self.N)+1:

				t[n] = int(2**(self.N-n))
				k_opt = -2*t[n]
				phase_adwin[n] = 0.5*np.angle(-1j*p_imag[1+2*t[n]]+p_real[1+2*t[n]])

				#nr_ones, nr_zeros = self.M_ssro_msmnt (theta_n = phase_adwin[n], t_n = int(t[n]))	

				max_k = int(t[n])*(self.M+3)
				if (max_k>range_k_sweep):
					nr_k  = int(range_k_sweep)/int(t[n])
					max_k = int(t[n])*nr_k
				k_space = np.arange(0, max_k, int(t[n]))

				for m in np.arange(self.M):				

					if (ext_outcomes == []):
						m_res = self.ramsey (t = t[n]*self.t0, theta = phase_adwin[n])
					else:
						ext_outcomes[n-1] = ext_outcomes[n-1] - 1
						if (ext_outcomes[n-1]<0):
							m_res = 0
						else:
							m_res = 1
					outcomes[n] = outcomes[n] + m_res

					p0_real = np.copy (p_real)
					p0_imag = np.copy (p_imag)
					cn = m_res*np.pi + phase_adwin[n]

					for k in k_space:
						if (k+t[n]>max(k_space)):
							p0r_pl = 0
							p0i_pl = 0
						else:
							p0r_pl = p0_real [1+k+t[n]]
							p0i_pl = p0_imag [1+k+t[n]]
						if (k<t[n]):
							p_im_min = -p0_imag[1+np.abs(k-t[n])]
						else:
							p_im_min = p0_imag[1+np.abs(k-t[n])]

						if use_fid_bayesian_update:
							if (int(m_res) == 0):
								A = 1-0.5*(self.fid0+self.fid1)
							elif (int(m_res) == 1):
								A = 0.5*(self.fid0+self.fid1)
							B = 0.5*(self.fid1-self.fid0)*np.exp(-(t[n]*self.t0/self.T2)**2)
							p_real [1+k] = A*p0_real[1+k] - 0.5*B*(np.cos(cn)*(p0_real [1+np.abs(k-t[n])] + p0r_pl) + np.sin(cn)*(p_im_min - p0i_pl)) 
							p_imag [1+k] = A*p0_imag[1+k] - 0.5*B*(np.cos(cn)*(p_im_min + p0i_pl) - np.sin(cn)*(p0_real [1+np.abs(k-t[n])] - p0r_pl)) 
						else:
							p_real [1+k] = 0.5*p0_real[1+k] + 0.25*(np.cos(cn)*(p0_real [1+np.abs(k-t[n])] + p0r_pl) + np.sin(cn)*(p_im_min - p0i_pl)) 
							p_imag [1+k] = 0.5*p0_imag[1+k] + 0.25*(np.cos(cn)*(p_im_min + p0i_pl) - np.sin(cn)*(p0_real [1+np.abs(k-t[n])] - p0r_pl)) 

					norm = p_real[1]*2*np.pi
					p_real = p_real/norm
					p_imag = p_imag/norm

			self.msmnt_results [r, :] = np.copy(outcomes[1:])
			self.msmnt_phases [r, :] = np.copy(phase_adwin[1:])
			self.msmnt_times = np.copy(t[1:])
			self.inc_rep()

		if (ext_outcomes != []):
			print ' ------ Adwin simulation: ------'
			print 'Outcomes: ', outcomes [1:]
			print 'Phases: ', np.round(phase_adwin [1:]*180/np.pi)



	def adwin_optimal_looping_storage (self, ext_outcomes = [], do_plot = False):
			
		#adwin-designed code with optimized looping and optimized 
		#storage of relevant coefficients {p[k]}

		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))

		self.reset_rep_counter()
		range_k_sweep = self.G+self.F*self.N+10
		print 'range_k_sweep: ', range_k_sweep

		if (ext_outcomes != []):
			self.reps = 1

		for r in np.arange(self.reps):

			if np.mod(r, 100)==0:
				print str(r)+'/'+str(self.reps)
			p_real = np.zeros (range_k_sweep)
			p_imag = np.zeros (range_k_sweep)
			p_real [1] = 1/(2*np.pi)
			p0_real = np.copy (p_real)
			p0_imag = np.copy (p_imag)

			t = np.zeros(self.N+1)
			outcomes = np.zeros(self.N+1)
			phase_adwin = np.zeros(self.N+1)

			for n in np.arange(self.N)+1:

				t[n] = int(2**(self.N-n))
				phase_adwin[n] = 0.5*np.angle(-1j*p_imag[3]+p_real[3])
				
				M = self.G+self.F*(self.N-n)
				max_k = M+3
				k_space = np.arange(max_k)
				B = 0.5*(self.fid1-self.fid0)*np.exp(-(t[n]*self.t0/self.T2)**2)

				if do_plot:
					print '### n = ', n
					f, axarr = plt.subplots(2, sharex=True, figsize=(10,10))

				for m in np.arange(M):				

					if do_plot:
						axarr[0].plot (np.arange (range_k_sweep-1)*t[n], p_real[1:])
						axarr[1].plot (np.arange (range_k_sweep-1)*t[n], p_imag[1:])
						#axarr[0].plot (np.arange (range_k_sweep-1)*t[n], p_real[1:], 'o')
						#axarr[1].plot (np.arange (range_k_sweep-1)*t[n], p_imag[1:], 'o')

					if (ext_outcomes == []):
						m_res = self.ramsey (t = t[n]*self.t0, theta = phase_adwin[n])
					else:
						ext_outcomes[n-1] = ext_outcomes[n-1] - 1
						if (ext_outcomes[n-1]<0):
							m_res = 0
						else:
							m_res = 1
					outcomes[n] = outcomes[n] + m_res
					cn = m_res*np.pi + phase_adwin[n]

					if (int(m_res) == 0):
						A = 1-0.5*(self.fid0+self.fid1)
					elif (int(m_res) == 1):
						A = 0.5*(self.fid0+self.fid1)

					#update rule:
					p_real [1] = A*p0_real[1] - B*(np.cos(cn)*p0_real [2] - np.sin(cn)*p0_imag[2]) 
					p_imag [1] = A*p0_imag[1]
					for k in np.arange(max_k)+1:
						p_real [1+k] = A*p0_real[1+k] - 0.5*B*np.cos(cn)*(p0_real [k] + p0_real [k+2]) - 0.5*B*np.sin(cn)*(p0_imag [k] - p0_imag [k+2])
						p_imag [1+k] = A*p0_imag[1+k] - 0.5*B*np.cos(cn)*(p0_imag [k] + p0_imag [k+2]) + 0.5*B*np.sin(cn)*(p0_real [k] - p0_real [k+2])

					norm = p_real[1]*2*np.pi
					p_real = p_real/norm
					p_imag = p_imag/norm

					p0_real = np.copy (p_real)
					p0_imag = np.copy (p_imag)

				#Strech k_space operation!
				for k in np.arange(range_k_sweep/2):
					p_real[1+2*k] = p0_real[1+k]
					p_imag[1+2*k] = p0_imag[1+k]
				p0_real = np.copy (p_real)
				p0_imag = np.copy (p_imag)

				if do_plot:
					plt.show()

			self.msmnt_results [r, :] = np.copy(outcomes[1:])
			self.msmnt_phases [r, :] = np.copy(phase_adwin[1:])
			self.msmnt_times = np.copy(t[1:])
			self.inc_rep()

		if (ext_outcomes != []):
			print ' ------ Adwin simulation: ------'
			print 'Outcomes: ', outcomes [1:]
			print 'Phases: ', np.round(phase_adwin [1:]*180/np.pi)



<<<<<<< HEAD

		
=======
	def compare_adwin_python_optimal_looping_storage (self, do_plot=False, newer_than = None, ext_outcomes = []):
		
		#compares {p_k} distributions obtaind by the adwin with 
		#the ones calculatd in python ("optimal" looping and storage)

		range_k_sweep = self.G+self.F*self.N+10
		print '## range_k_sweep: ', range_k_sweep
		p_real = np.zeros (range_k_sweep)
		p_imag = np.zeros (range_k_sweep)
		p_real [1] = 1/(2*np.pi)
		p0_real = np.copy (p_real)
		p0_imag = np.copy (p_imag)

		t = np.zeros(self.N+1)
		outcomes = np.zeros(self.N+1)
		phase_adwin = np.zeros(self.N+1)
		phase_simulation = np.zeros(self.N+1)
		outcomes = np.zeros(self.N+1)
		#ext_outcomes = self.outcomes

		label = ''
		for i in np.arange(self.N):
			label=label+str(ext_outcomes[i])
		print 'Label: ', label


		t = np.zeros(self.N+1)
		outcomes = np.zeros(self.N+1)
		phase_adwin = np.zeros(self.N+1)

		for n in np.arange(self.N)+1:

			t[n] = int(2**(self.N-n))
			phase_adwin[n] = 0.5*np.angle(-1j*p_imag[3]+p_real[3])
			
			M = self.G+self.F*(self.N-n)
			max_k = M+3
			k_space = np.arange(max_k)
			B = 0.5*(self.fid1-self.fid0)*np.exp(-(t[n]*self.t0/self.T2)**2)

			if do_plot:
				print '### n = ', n
				f, axarr = plt.subplots(2, sharex=True, figsize=(10,10))

			for m in np.arange(M)+1:				
				exp_exists = True
				f_name = 'N='+str(self.N)+'_M=(5, 0)_rtAdwin_'+label+'_test_pk_(n='+str(n)+'_m='+str(m)+')'

				try:
					f = toolbox.latest_data (contains = f_name, return_all=False, newer_than = newer_than)
					exp = adptv_mgnt.RamseySequence_Exp (folder = f)
					exp.load_exp_data()
					phase_adwin = exp.opt_phase[:self.N]
					print 'Phase adwin (exp): ', phase_adwin
					self.fid0 = exp.exp_fid0
					self.fid1 = exp.exp_fid1
					#self.T2 = exp.T2
					exp_exists = True
				except:
					exp_exists = False
					print f_name+' does not exist!'

				if (ext_outcomes == []):
					m_res = self.ramsey (t = t[n]*self.t0, theta = phase_adwin[n])
				else:
					ext_outcomes[n-1] = ext_outcomes[n-1] - 1
					if (ext_outcomes[n-1]<0):
						m_res = 0
					else:
						m_res = 1
				outcomes[n] = outcomes[n] + m_res
				cn = m_res*np.pi + phase_adwin[n]

				if (int(m_res) == 0):
					A = 1-0.5*(self.fid0+self.fid1)
				elif (int(m_res) == 1):
					A = 0.5*(self.fid0+self.fid1)

				#update rule:
				p_real [1] = A*p0_real[1] - B*(np.cos(cn)*p0_real [2] - np.sin(cn)*p0_imag[2]) 
				p_imag [1] = A*p0_imag[1]
				for k in np.arange(max_k)+1:
					p_real [1+k] = A*p0_real[1+k] - 0.5*B*np.cos(cn)*(p0_real [k] + p0_real [k+2]) - 0.5*B*np.sin(cn)*(p0_imag [k] - p0_imag [k+2])
					p_imag [1+k] = A*p0_imag[1+k] - 0.5*B*np.cos(cn)*(p0_imag [k] + p0_imag [k+2]) + 0.5*B*np.sin(cn)*(p0_real [k] - p0_real [k+2])

				norm = p_real[1]*2*np.pi
				p_real = p_real/norm
				p_imag = p_imag/norm

				p0_real = np.copy (p_real)
				p0_imag = np.copy (p_imag)

				if exp_exists:
					len_array = len(p_real)
					th_sim = int(round(phase_simulation[n]*180/np.pi))
					th_adw = phase_adwin[n-1]

					print '######  n = ', n, ' ---- m = ', m
					print '   phase_th = ', th_sim, '  --- phase_exp = ', th_adw
					print 'cos(cn) = ', np.cos(cn), '    sin(cn) = ', np.sin(cn)

					if do_plot:
						#diff = (((p_real[1:2**self.N+2]-exp.real_pk_adwin[:2**self.N]))**2+((p_imag[1:2**self.N+2]-exp.imag_pk_adwin[:2**self.N]))**2)**0.5
						#rel_diff = 100*diff/((p_real[1:2**self.N+2]**2+p_imag[1:2**self.N+2]**2)**0.5)
						#avg_diff = np.sum(diff)/len(diff)
						f, axarr = plt.subplots(3, sharex=True, figsize=(10,10))
						x = np.arange(2**self.N+1)

						print exp.real_pk_adwin
						axarr[0].plot (np.arange (range_k_sweep-1)*t[n],p_real[1:], ':k')
						axarr[0].plot (np.arange (range_k_sweep-1)*t[n], exp.real_pk_adwin[:-1], ':b')
						axarr[1].plot (np.arange (range_k_sweep-1)*t[n], p_imag[1:], ':k')
						axarr[1].plot (np.arange (range_k_sweep-1)*t[n], exp.imag_pk_adwin[:-1], ':b')
						axarr[0].plot (np.arange (range_k_sweep-1)*t[n], p_real[1:], 'o', label='sim_'+str(n)+','+str(m))
						axarr[1].plot (np.arange (range_k_sweep-1)*t[n], p_imag[1:], 'o', label='sim_'+str(n)+','+str(m))
						axarr[0].plot (np.arange (range_k_sweep-1)*t[n], exp.real_pk_adwin[:-1], 'o', label='exp_'+str(n)+','+str(m))
						axarr[1].plot (np.arange (range_k_sweep-1)*t[n], exp.imag_pk_adwin[:-1], 'o', label='exp_'+str(n)+','+str(m))
						#axarr[2].plot (np.arange (range_k_sweep-1)*t[n], rel_diff, ':k')
						#axarr[2].plot (np.arange (range_k_sweep-1)*t[n], rel_diff, 'ok')
						axarr[0].set_title('real part, n = '+str(n)+', m = '+str(m))
						axarr[1].set_title('imaginary part, n = '+str(n)+', m = '+ str(m))
						axarr[2].set_title('difference, n = '+str(n)+', m = '+ str(m))
						axarr[0].set_xlim([0, (range_k_sweep-2)*t[n]])
						axarr[1].set_xlim([0, (range_k_sweep-2)*t[n]])
						axarr[2].set_xlim([0, range_k_sweep*t[n]])
						axarr[0].set_ylim([-0.2, 0.2])
						axarr[1].set_ylim([-0.2, 0.2])
						axarr[2].set_ylim([0, 100])					
						axarr[2].set_ylabel('rel diff %', fontsize=15)
						axarr[0].legend()
						axarr[1].legend()
						axarr[2].legend()
						#time_label = time.strftime ('%Y%m%d_%H%M%S')
						#fname = 'M:/tnw/ns/qt/Diamond/Projects/Magnetometry with adaptive measurements/Data/adwin_tests/'+time_label+'_N='+str(self.N)+'_M='+str(self.M)+'_pk_'+str(n)+'_'+str(m)+'.png'
						#print fname
						#f.savefig (fname)
						plt.show()

			#Strech k_space operation!
			for k in np.arange(range_k_sweep/2):
				p_real[1+2*k] = p0_real[1+k]
				p_imag[1+2*k] = p0_imag[1+k]
			p0_real = np.copy (p_real)
			p0_imag = np.copy (p_imag)


	
>>>>>>> e1d74dbc3c6181d681a6a5c611773b57b7538bde
	def adwin_only_positive (self, debug = False, exec_speed = False):

		#use only {p[k], k>=0}, since p[-k]=p*[k] for real probability distribution

		discr_steps = 2**(self.N+1)+1
		m = np.zeros (self.N+1)
		t = np.zeros (self.N+1)
		th = np.zeros(self.N+1)
		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))
		exec_time = None
		
		if exec_speed:
			start = time.clock()
		for rep in np.arange(self.reps):
			p_real = np.zeros (discr_steps)
			p_imag = np.zeros (discr_steps)
			p_real [0] = 1./(2*np.pi)

			curr_phase = 0
			for n in np.arange(self.N)+1:
				t[n] = 2**(self.N-n)
				k_opt = 2**(self.N-n+1)
				#print k_opt, p_real[k_opt]-1j*p_imag[k_opt]
				th[n] = -0.5*np.angle(-1j*p_imag[k_opt]+p_real[k_opt])
				m[n] = self.ramsey (theta=th[n], t = t[n]*self.t0)

				if debug:
					print '### n =', n
					print 't = ', t[n], '  theta = ', th[n]*180./np.pi, '  --- res: ', m[n] 
				#update rule:
				cn = m[n]*np.pi+th[n]
				p0_real = np.copy (p_real)
				p0_imag = np.copy (p_imag)

				for k in np.arange(t[n]+1, 2**(self.N)+1):
					p_real [k] = 0.5*p0_real[k] + 0.25*(np.cos(cn)*(p0_real [k-t[n]] + p0_real [k+t[n]]) - np.sin(cn)*(p0_imag [k-t[n]] - p0_imag [k+t[n]])) 
					p_imag [k] = 0.5*p0_imag[k] + 0.25*(np.cos(cn)*(p0_imag [k-t[n]] + p0_imag [k+t[n]]) + np.sin(cn)*(p0_real [k-t[n]] - p0_real [k+t[n]])) 
		
				for k in np.arange(0, t[n]+1):
					p_real [k] = 0.5*p0_real[k] + 0.25*(np.cos(cn)*(p0_real [k-t[n]] + p0_real [k+t[n]]) - np.sin(cn)*(-p0_imag [t[n]-k] - p0_imag [k+t[n]])) 
					p_imag [k] = 0.5*p0_imag[k] + 0.25*(np.cos(cn)*(-p0_imag [t[n]-k] + p0_imag [k+t[n]]) + np.sin(cn)*(p0_real [k-t[n]] - p0_real [k+t[n]])) 

			self.msmnt_results [rep, :] = m[1:]
			self.msmnt_phases [rep, :] = th[1:]
			self.msmnt_times = t [1:]
		if exec_speed:
			stop = time.clock()
			exec_time = stop-start
		return exec_time

	def adwin_positive_nonzero (self, debug = False, exec_speed=False):

		#use only {p[k], k>=0}, since p[-k]=p*[k] for real probability distribution

		discr_steps = 2**(self.N+1)+1
		m = np.zeros (self.N+1)
		t = np.zeros (self.N+1)
		th = np.zeros(self.N+1)
		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))
		exec_time = None
		
		#print 'Create \'relevant-elements\' sequence: '
		elems1 = [0]
		elems2 = [2**(self.N-1)]
		for i in np.arange(self.N-1)+2:
			tn = 2**(self.N-i)
			arr = np.ndarray.tolist(np.array(elems1)+tn)
			elems1 = elems1+arr
			arr = np.ndarray.tolist(np.array(elems2)+tn)
			elems2 = elems2+arr
		elems1 = np.squeeze(np.array(elems1))
		elems2 = np.squeeze(np.array(elems2))
		#print elems1
		#print elems2

		if exec_speed:
			start = time.clock()

		for rep in np.arange(self.reps):
			p_real = np.zeros (discr_steps)
			p_imag = np.zeros (discr_steps)
			p_real [0] = 1./(2*np.pi)

			curr_phase = 0
			for n in np.arange(self.N)+1:
				t[n] = 2**(self.N-n)
				k_opt = 2**(self.N-n+1)
				th[n] = -0.5*np.angle(-1j*p_imag[k_opt]+p_real[k_opt])
				m[n] = self.ramsey (theta=th[n], t = t[n]*self.t0)

				if debug:
					print '### n =', n
					print 't = ', t[n], '  theta = ', th[n]*180./np.pi, '  --- res: ', m[n] 
				#update rule:
				cn = m[n]*np.pi+th[n]
				p0_real = np.copy (p_real)
				p0_imag = np.copy (p_imag)

				tt2 = elems2[:2**(n-1)]
				for k in tt2:
					p_real [k] = 0.5*p0_real[k] + 0.25*(np.cos(cn)*(p0_real [k-t[n]] + p0_real [k+t[n]]) - np.sin(cn)*(p0_imag [k-t[n]] - p0_imag [k+t[n]])) 
					p_imag [k] = 0.5*p0_imag[k] + 0.25*(np.cos(cn)*(p0_imag [k-t[n]] + p0_imag [k+t[n]]) + np.sin(cn)*(p0_real [k-t[n]] - p0_real [k+t[n]])) 

				tt1 = elems1[:2**(n-1)]
				for k in tt1:
					p_real [k] = 0.5*p0_real[k] + 0.25*(np.cos(cn)*(p0_real [k-t[n]] + p0_real [k+t[n]]) - np.sin(cn)*(p0_imag [k-t[n]] - p0_imag [k+t[n]])) 
					p_imag [k] = 0.5*p0_imag[k] + 0.25*(np.cos(cn)*(p0_imag [k-t[n]] + p0_imag [k+t[n]]) + np.sin(cn)*(p0_real [k-t[n]] - p0_real [k+t[n]])) 

			self.msmnt_results [rep, :] = m[1:]
			self.msmnt_phases [rep, :] = th[1:]
			self.msmnt_times = t [1:]
		if exec_speed:
			stop = time.clock()
			exec_time = stop-start
		return exec_time

	def adwin_ultrafast (self, debug=False, exec_speed = False):

		discr_steps = 2**(self.N)+1
		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))
		exec_time = None

		if exec_speed:
			start = time.clock()

		for rep in np.arange(self.reps):
			m = np.zeros (self.N+1)
			t = np.zeros (self.N+1)
			th = np.zeros(self.N+1)

			p_real = np.zeros (discr_steps)
			p_imag = np.zeros (discr_steps)
			p_real [0] = 1./(2*np.pi)
			t[0] = 2**self.N

			for n in np.arange(self.N)+1:
				t[n] = int(2**(self.N-n))
				k_opt = int(2**(self.N-n+1))
				th[n] = -0.5*np.angle(-1j*p_imag[k_opt]+p_real[k_opt])

				for mmm in np.arange(self.M):
					meas_res = self.ramsey (theta=th[n], t = 2**(self.N-n)*self.t0)
					m[n] = m[n] + meas_res
					cn = meas_res*np.pi+th[n]
					p0_real = np.copy (p_real)
					p0_imag = np.copy (p_imag)
					k = t[n]
					p_real [k] = 0.5*p0_real[k] + 0.25*(np.cos(cn)*(p0_real [0] + p0_real [2**(self.N-n+1)]) - np.sin(cn)*(p0_imag [0] - p0_imag [2**(self.N-n+1)])) 
					p_imag [k] = 0.5*p0_imag[k] + 0.25*(np.cos(cn)*(p0_imag [0] + p0_imag [2**(self.N-n+1)]) + np.sin(cn)*(p0_real [0] - p0_real [2**(self.N-n+1)])) 
			self.msmnt_results [rep, :] = m[1:]
			self.msmnt_phases [rep, :] = th[1:]
			self.msmnt_times = t [1:]
		if exec_speed:
			stop = time.clock()
			exec_time = stop-start
		return exec_time


	def adwin_ultrafast_print_steps (self, msmnt_results = []):

		if (len(msmnt_results) != self.N):
			print 'Incorrect input!'
		else:
			discr_steps = 2**(self.N)+1
			m = np.zeros (self.N+1)
			m[1:] = msmnt_results
			t = np.zeros (self.N+1)
			th = np.zeros(self.N+1)

			p_real = np.zeros (discr_steps)
			p_imag = np.zeros (discr_steps)
			p_real [0] = 1./(2*np.pi)
			t[0] = 2**self.N

			if self.verbose:
				print '###############################################'
				print 'Measurement result sequence: ', msmnt_results
			for n in np.arange(self.N)+1:
				t[n] = int(2**(self.N-n))
				k_opt = int(2**(self.N-n+1))
				th[n] = -0.5*np.angle(-1j*p_imag[k_opt]+p_real[k_opt])
				k = t[n]

				if self.verbose:
					print ' * n = ', n, ' ----> tn = ', t[n],  ' - k_opt + 1 = ', k_opt+1, ' - th_opt = ', int(th[n]*180/(np.pi))
					print '  		p[k_opt+1] = ' +str('{0:.4f}'.format(-1j*p_imag[k_opt]+p_real[k_opt]))
					print '			angle = ', int(np.angle(-1j*p_imag[k_opt]+p_real[k_opt])*180/np.pi)
					print '			Important coefficients: '
					print '				p[0+1] = '+ str('{0:.4f}'.format(p_real[0]+1j*p_imag[0]))
					print '				p[tn+1] = '+str('{0:.4f}'.format(p_real[t[n]]+1j*p_imag[t[n]]))
					print '				p[2*tn+1] = '+str('{0:.4f}'.format(p_real[2*t[n]]+1j*p_imag[2*t[n]]))

				nr_ones = m[n]
				nr_zeros = self.M - nr_ones
				for i in np.arange(nr_ones):
					cn = np.pi+th[n]
					p0_real = np.copy (p_real)
					p0_imag = np.copy (p_imag)
					p_real [k] = 0.5*p0_real[k] + 0.25*(np.cos(cn)*(p0_real [0] + p0_real [2*k]) - np.sin(cn)*(p0_imag [0] - p0_imag [2*k])) 
					p_imag [k] = 0.5*p0_imag[k] + 0.25*(np.cos(cn)*(p0_imag [0] + p0_imag [2*k]) + np.sin(cn)*(p0_real [0] - p0_real [2*k])) 
				for i in np.arange(nr_zeros):
					cn = th[n]
					p0_real = np.copy (p_real)
					p0_imag = np.copy (p_imag)
					p_real [k] = 0.5*p0_real[k] + 0.25*(np.cos(cn)*(p0_real [0] + p0_real [2*k]) - np.sin(cn)*(p0_imag [0] - p0_imag [2*k])) 
					p_imag [k] = 0.5*p0_imag[k] + 0.25*(np.cos(cn)*(p0_imag [0] + p0_imag [2*k]) + np.sin(cn)*(p0_real [0] - p0_real [2*k])) 


		return np.round(th[1:]*180/(np.pi))


	def adwin_update (self, p_real, p_imag, meas_res, phase, k, tn):
		cn = meas_res*np.pi+phase
		p0_real = np.copy (p_real)
		p0_imag = np.copy (p_imag)


		
		try:
			if (k>=tn):	
				p_real [k] = 0.5*p0_real[k] + 0.25*(np.cos(cn)*(p0_real [k-tn] + p0_real [k+tn]) - np.sin(cn)*(p0_imag [k-tn] - p0_imag [k+tn])) 
				p_imag [k] = 0.5*p0_imag[k] + 0.25*(np.cos(cn)*(p0_imag [k-tn] + p0_imag [k+tn]) + np.sin(cn)*(p0_real [k-tn] - p0_real [k+tn])) 
			else:
				p_real [k] = 0.5*p0_real[k] + 0.25*(np.cos(cn)*(p0_real [tn-k] + p0_real [k+tn]) - np.sin(cn)*(p0_imag [tn-k] - p0_imag [k+tn])) 
				p_imag [k] = 0.5*p0_imag[k] + 0.25*(np.cos(cn)*(-p0_imag [tn-k] + p0_imag [k+tn]) + np.sin(cn)*(p0_real [tn-k] - p0_real [k+tn])) 
		except:
			pass
		return p_real, p_imag



	def adwin_ultrafast_M (self, debug=False, nr_coeff = 1, exec_speed=False):

		self.max_k = 2**self.N
		discr_steps = 2*self.max_k+1
		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))
		exec_time = None

		if exec_speed:
			start = time.clock()

		for rep in np.arange(self.reps):
			m = np.zeros (self.N+1)
			t = np.zeros (self.N+1)
			th = np.zeros(self.N+1)

			p_real = np.zeros (discr_steps)
			p_imag = np.zeros (discr_steps)
			p_real [0] = 1./(2*np.pi)
			t[0] = 2**self.N

			for n in np.arange(self.N)+1:
				t[n] = int(2**(self.N-n))
				k_opt = int(2**(self.N-n+1))
				th[n] = -0.5*np.angle(-1j*p_imag[k_opt]+p_real[k_opt])

				for mmm in np.arange(self.M):
					meas_res = self.ramsey (theta=th[n], t = 2**(self.N-n)*self.t0)
					m[n] = m[n] + meas_res
					for j in np.arange(nr_coeff)-nr_coeff/2:
						p_real, p_imag = self.adwin_update (p_real = p_real, p_imag = p_imag, meas_res = meas_res, phase = th[n], tn = t[n], k=k_opt+j*t[n])
					if debug:
						plt.plot (p_real, label = str(mmm+1))
				if debug:
					plt.title ('N = '+str(n), fontsize=15)
					plt.legend()
					plt.show()

			self.msmnt_results [rep, :] = m[1:]
			self.msmnt_phases [rep, :] = th[1:]
			self.msmnt_times = t [1:]
		if exec_speed:
			stop = time.clock()
			exec_time = stop-start
		return exec_time

	def adwin_ultrafast_M_set_input (self, msmnt_result = [],  nr_coeff = 1):

		self.max_k = 2**self.N
		discr_steps = 2*self.max_k+1
		self.msmnt_phases = np.zeros((self.reps,self.N))
		self.msmnt_times = np.zeros(self.N)
		self.msmnt_results = np.zeros((self.reps,self.N))

		m = np.zeros (self.N+1)
		t = np.zeros (self.N+1)
		th = np.zeros(self.N+1)

		p_real = np.zeros (discr_steps)
		p_imag = np.zeros (discr_steps)
		p_real [0] = 1./(2*np.pi)
		t[0] = 2**self.N

		for n in np.arange(self.N)+1:
			t[n] = int(2**(self.N-n))
			k_opt = int(2**(self.N-n+1))
			th[n] = -0.5*np.angle(-1j*p_imag[k_opt]+p_real[k_opt])

			nr_ones = msmnt_result [n-1]
			nr_zeros = self.M - nr_ones

			for mmm in np.arange(nr_ones):
				for j in np.arange(nr_coeff)-nr_coeff/2:
					if (k_opt+j*t[n]>=0):
						p_real, p_imag = self.adwin_update (p_real = p_real, p_imag = p_imag, meas_res = 1, phase = th[n], tn = t[n], k=k_opt+j*t[n])
			for mmm in np.arange(nr_zeros):
				for j in np.arange(nr_coeff)-nr_coeff/2:
					if (k_opt+j*t[n]>=0):
						p_real, p_imag = self.adwin_update (p_real = p_real, p_imag = p_imag, meas_res = 0, phase = th[n], tn = t[n], k=k_opt+j*t[n])

		return th[1:]



	def check_adwin_phases(self, nr_coeff):

		sim_result = np.random.randint(self.M+1, size=self.N)
		print "Simulated outcome: ", sim_result

		th_adwin = np.round(self.adwin_ultrafast_M_set_input (nr_coeff=nr_coeff, msmnt_result = sim_result)*180/np.pi)
		print 'Adwin simualted phases: ', th_adwin

		t = AdaptiveTable (N=self.N,M=self.M)
		t.set_tau0(tau0=self.t0)
		t.verbose = False
		sim_ph, pos = t.msmnt_to_position (msmnt_results = sim_result)
		print 'Python simulated phases: ', sim_ph

		diff = np.sum(np.abs(th_adwin-sim_ph))/self.N
		print 'MEAN ERROR: ', diff, ' deg'


