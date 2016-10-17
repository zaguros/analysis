

import numpy as np
from matplotlib import rc, cm

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot
from matplotlib import rc, cm
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
from analysis.lib.magnetometry import adwin_debug_magnetometry as adwin_mgnt

reload (sequence)
reload(magnetometry)
reload(adwin_mgnt)


class Didactic(magnetometry.RamseySequence_Simulation):


	def didactic_cappellaro (self):

		if (self.G+self.F+self.K==0):
			print 'Simulation parameters G, K, F not set!!'
		else:			
			self.reps = 1
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
				self.plot_p_k()
		
				for i,k in enumerate(k_array):

					t[i] = int(2**(k))
					ttt = -2**(k+1)
					
					phase[i] = 0.5*np.angle (self.p_k[ttt+self.points])
					m_total = 0
					MK = self.G+self.F*(self.K-k)

					for m in np.arange(MK):
						m_res = self.ramsey (theta=phase[i], t = t[i]*self.t0)#self.majority_vote_msmnt (theta_n = phase[i], t_n = t[i])					
						self.bayesian_update (m_n = m_res, phase_n = phase[i], t_n = 2**(k))
						m_total = m_total + m_res
					msmnt_results[i] =m_total
					self.plot_p_k()

					
				self.msmnt_results [r, :] = np.copy(msmnt_results)
				self.msmnt_phases [r, :] = np.copy(phase)
				self.inc_rep()


s = Didactic(N_msmnts = 9, tau0=20e-9, reps=101)
s.setup_simulation (magnetic_field_hz = -20e6, G=10,F=5,K=8)
s.fid0=1
s.fid1=0
s.didactic_cappellaro()
s.plot_updates()

