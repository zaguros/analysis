
import numpy as np
import math
import pylab as plt
import matplotlib
from matplotlib import rc, cm


class CavitySims ():

	def __init__(self):

		#NV parameters
		#experimental values taken from PRL 110, 243602
		self.freq = np.array([469.6e12, 461.9e12,452.6e12, 439.7e12, 429.0e12, 416.2e12, 403.2e12, 392.4e12])
		self.A = np.array([1520., 5260., 18600., 16400., 14000., 9180., 6570., 3270.])
		self.d = self.freq - self.freq[0] 
		#self.freqC = freq[0]

		self.A_tot = np.sum(self.A)
		self.epsilon = self.A/self.A_tot #relative strength transitions
		self.gm_tot = 35e6 #deduced by g2 msmnts
		self.g_relative = self.gm_tot*self.epsilon

		#CAVITY parameter
		self.cavity_length = 1.1e-6 # only air length
		self.Q = None
		self.T = None
		self.diamond_thickness = 3e-6
		self.radius_curvature = 15e-6
		self.wavelength_ZPL = 3e8/float(self.freq[0])

		self._linewidth = False

		#simulation parameters
		#self.N_sim = 50000

	def set_mirror_curvature (self, value):
		self.radius_curvature = value

	def set_diamond_thickness (self, value):
		self.diamond_thickness = value

	def set_cavity_length (self, value):
		ret = 0
		self.cavity_length = value
		d = self.diamond_thickness*2.5 + self.cavity_length
		if (d>self.radius_curvature-2e-6):
			print "Cavity is unstable!"
			ret = 1
		#self.waist = ((self.wavelength_ZPL/np.pi)**0.5)*(d*(self.radius_curvature-d))**(1/4.)
		#self.cavity_volume = 0.5*np.pi*self.waist**2*d
		self.waist = ((self.wavelength_ZPL/np.pi)**0.5)*(d*(self.radius_curvature))**(1/4.) # ref New.J.Phys
		self.cavity_volume = np.pi*((0.5*self.waist)**2)*d
		return ret

	def set_temperature(self, value):
		self.T = value
		dephZPL = (6.28*16.2+9.2e-7*(1/0.0125)*self.T**5)*1e6/(2*math.pi) #in Hz. Only for T < 100 K!, do you need 2pi again?! 
		self.linewidths = np.array([dephZPL, 15.9e12, 15.5e12, 15.0e12, 16.5e12, 12.7e12, 13.7e12, 16.1e12])
		self.gm_star = 2*np.pi*self.linewidths[0] 

	def set_ZPL_linewidth (self, value):
		self.linewidths = np.array([value, 15.9e12, 15.5e12, 15.0e12, 16.5e12, 12.7e12, 13.7e12, 16.1e12])
		self.gm_star = 2*np.pi*self.linewidths[0]		
		self._linewidth = True

	def set_Q (self, value,verbose=False):
		self.Q = value
		self.k = np.pi*3e8/(self.wavelength_ZPL*self.Q)
		if verbose:
			print 'k',self.k/1.e9,'GHz'

	def calculate_params (self, do_plot = False, verbose=False):
		
		self.FSR = 3e8/(2*(self.diamond_thickness*2.5+self.cavity_length)) #FSR
		self.resonant_frq = 3e8/self.wavelength_ZPL

		#longitudinal cavity modes
		self.lambda_modes_nm=[]
		self.freq_modes = []
		l=0
		i=0
		while ((l<800) and (i< 1000)):
			nu = self.resonant_frq-i*self.FSR
			l = 1e9*3e8/nu 
			self.lambda_modes_nm.append(l)
			self.freq_modes.append(i*self.FSR)
			i = i+1

		if do_plot:
			plt.figure()
			plt.plot (self.lambda_modes_nm, 'o')
			plt.ylabel ('wavelength [nm]')
			plt.show()

		#mode volume
		d = self.diamond_thickness*2.5 + self.cavity_length
		#self.waist = ((self.wavelength_ZPL/np.pi)**0.5)*(d*(self.radius_curvature-d))**(1/4.)
		self.waist = ((self.wavelength_ZPL/np.pi)**0.5)*(d*(self.radius_curvature))**(1/4.) # ref New.J.Phys

		#self.cavity_volume = 0.5*np.pi*self.waist**2*d
		self.cavity_volume = np.pi*((0.5*self.waist)**2)*d

		if verbose:
			print 'lambda_ZPL = ', self.wavelength_ZPL*1e9, '[nm]'
			print 'FSR = ', 1e9*((self.wavelength_ZPL**2)/3e8)*self.FSR, '[nm]'
			print "mode volume: ", self.cavity_volume*1e18,  "micron^3"

		#T^5 dephasing ZPL, linewidths = w
		if not(self._linewidth):
			dephZPL = (16.2+9.2e-7*(1/0.0125)*self.T**5)*1e6 #in Hz. This is w - thus differs a factor 2pi from Jahn-Teller paper
			self.linewidths = np.array([dephZPL, 15.9e12, 15.5e12, 15.0e12, 16.5e12, 12.7e12, 13.7e12, 16.1e12])
		self.gm_star = 2*np.pi*self.linewidths[0]
		self.g = np.sqrt((3*3e8*(637e-9**2)*self.gm_tot/2)/(4*np.pi*self.cavity_volume))
		self.gm_relax = 2*np.pi*self.linewidths-2*np.pi*self.linewidths[0]
		self.g_all = self.g*np.sqrt(self.epsilon)
		self.gamma=(self.k+self.gm_tot+self.gm_relax+self.gm_star)

	def simulate (self, do_plot = True, do_save = False, verbose=False):
		
		delta = np.zeros([len(self.linewidths), len(self.lambda_modes_nm)])

		i=-1
		for f_line in self.d:
			j = -1
			i = i+1
			for f_mode in self.freq_modes:
				j = j+1
				delta [i, j] = f_line - f_mode

		G_all = self.g_all
		Gamma = self.gamma
		for j in np.arange(len(self.freq_modes)-1):
			Gamma = np.vstack([Gamma, self.gamma])
			G_all = np.vstack([G_all, self.g_all])
		self.Gamma = Gamma.transpose()
		self.G_all = G_all.transpose()
		if verbose:
			print 'g',self.g_all
		self.R = (4*self.G_all**2/self.Gamma)*(1/(1+(2*delta/(self.Gamma))**2))
		self.F = self.R/float(self.gm_tot)
		self.F0 = self.F[0,0]
		if verbose:
			print 'Rij',self.R[:,0]
		self.P = self.R/ float(np.sum(np.sum(self.R))+self.gm_tot)
		self.ZPL_in_0 = self.P[0,0]
		self.P_tot = np.sum(self.P, axis=0)

		if do_plot:	
			plt.figure()
			plt.plot (self.lambda_modes_nm, self.P_tot, 'ob', markersize=2)
			plt.ylabel ('probability of emission in the mode')
			plt.xlabel ('wavelength [nm]')
			plt.show()


	def emission_in_ZPL (self, sweep_param = 'cavity_length', min_val=0, max_val=15e-6, nr_points=50, xlogscale=False):
		sweep_vals = np.linspace (min_val, max_val, nr_points)
		if xlogscale:
			sweep_vals = np.logspace(np.log10(min_val), np.log10(max_val), nr_points)

		emission_prob = np.zeros (nr_points)
		emission_prob_others = np.zeros(nr_points)
		emission_prob_ZPL_in_0 = np.zeros(nr_points)
		purcell = np.zeros(nr_points)
		ind = 0
		mode_vol = np.zeros(nr_points)

		print 'Q',self.Q
		for p in sweep_vals:
			err = 0
			if (sweep_param == 'cavity_length'):
				err = self.set_cavity_length (p) #err is 1 if cavity is unstable
			elif (sweep_param == 'temperature'):
				self.set_temperature (p)
			elif (sweep_param == 'gm_star'):
				self.set_ZPL_linewidth (p/(2*math.pi))
			elif (sweep_param == 'Q'):
				self.set_Q (p)
			
			self.calculate_params(do_plot = False)
			self.simulate(do_plot=False,verbose=False)
			if (err==0):
				emission_prob[ind] = self.P_tot[0]
				emission_prob_others[ind] = np.sum(self.P_tot[1:])
				emission_prob_ZPL_in_0[ind] = self.ZPL_in_0
				mode_vol[ind] = self.cavity_volume
				purcell[ind] = self.F0
				ind = ind + 1

		return sweep_vals[:ind], emission_prob[:ind], mode_vol[:ind], purcell[:ind],emission_prob_others[:ind],emission_prob_ZPL_in_0[:ind]


	def plot_emission_spectrum (self):

		ll = np.linspace (580, 1000, 10000)*1e-9
		nu = 3e8/ll

		y = np.zeros(len(nu))

		#devo fare un for-loop sui modi e plottarli come lorentziane
		for i in np.arange(len(self.freq)):

			y = y + 2*self.A[i]*(np.pi*self.linewidths[i])/(4*(nu-self.freq[i])**2+self.linewidths[i]**2)

		y = y/np.sum(y)
		plt.figure(figsize=(10,5))
		plt.plot (ll*1e9, y, linewidth=2)
		plt.xlim ([600, 800])
		plt.show()


matplotlib.rc ('xtick', labelsize=15)
matplotlib.rc ('ytick', labelsize=15)



c = CavitySims()
c._linewidth=False

F=10000
L=2e-6
lambdaa=0.637e-6#0.6393e-6
#L=9./2*lambdaa
Q = F*L/lambdaa

c.set_diamond_thickness (3e-6)
c.set_mirror_curvature (15e-6)
c.set_cavity_length (L)

c.set_Q (Q)
c.set_temperature (5)

#c.set_diamond_thickness (2.e-6)
# c.set_cavity_length (1.1e-6)
# c.set_mirror_curvature (7.6e-6)
c.set_ZPL_linewidth (2.44e12)

#c.calculate_params()
#c.simulate()

sweep_param1='temperature'
sweep_param2='Q'

x1, y10, V, F10,y10b,y10c = c.emission_in_ZPL(sweep_param = sweep_param1, min_val = 0, max_val=100, nr_points = 500)

c.set_temperature (5)
c._linewidth = False
# c.set_cavity_length (5.5e-6)
# x2, y20, V, F20,y20b,y20c = c.emission_in_ZPL(sweep_param = 'gm_star', min_val = 1e9, max_val= 100000e9, nr_points = 500, xlogscale=True)
# x2=x2/1.e9
# x2, y20, V, F20,y20b,y20c = c.emission_in_ZPL(sweep_param = sweep_param2, min_val = 0, max_val=100, nr_points = 500)
# x2, y20, V, F20,y20b,y20c = c.emission_in_ZPL(sweep_param = sweep_param2, min_val = 0.e-6, max_val=5.5e-6, nr_points = 500)
# x2=x2/1.e-6#convert to um
x2, y20, V, F20,y20b,y20c = c.emission_in_ZPL(sweep_param = sweep_param2, min_val = 0, max_val=200000, nr_points = 500)
#x2=x2/1.e-6#convert to um

#x, y_10000, V = c.emission_in_ZPL(max_val=15e-6)
fig = plt.figure (figsize=(16,6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

new_tick_locations = np.array([ 1000, 2500, 5000, 10000])

def tick_function(Q):
	linewidth = 1e-9*3e8/(637e-9*Q)
	return ["%.1f" % z for z in linewidth]

# ax2.set_xticks(new_tick_locations)
# ax2.set_xticklabels(tick_function(new_tick_locations))
ax1.plot (x1, y10, color = 'DarkViolet', linewidth = 4, label = 'zero-mode') #0.5 factor from losses
ax1.plot (x1, y10b, color = 'Green', linewidth = 4, label = 'all other modes') 
ax2.plot (x2, y20, color = 'DarkViolet', linewidth = 4, label = 'zero-mode') #0.5 factor from losses
ax2.plot (x2, y20b, color = 'Green', linewidth = 4, label = 'all other modes') 
# ax1.plot (x1, y10c, color = 'Green', linewidth = 4, label = 'd = 1 $\mu$m') #0.5 factor from losses
# ax2.plot (x2, y20c, color = 'crimson', linewidth = 4, label = 'd = 3 $\mu$m')
#ax1.plot (x1, 1+F10, '--', color = 'RoyalBlue', linewidth = 4, label = 'd = 1 $\mu$m') #0.5 factor from losses
#ax2.plot (x2, 1+F20, '--', color = 'crimson', linewidth = 4, label = 'd = 3 $\mu$m')
ax1.set_xlabel('Temperature (K)', fontsize = 20)
ax1.set_ylabel ('Probability of NV emission\ninto the cavity mode', fontsize = 20)
# ax2.set_xlim([0,5.5])
ax2.set_xlabel('Quality factor', fontsize = 20)
#ax2.set_ylabel ('Probability of NV emission\n(from ZPL & PSB)\ninto the cavity mode', fontsize = 20)

# ax2.set_xlabel('gamma * (GHz)', fontsize = 15)
# ax2.set_ylabel ('emission fraction', fontsize = 16)
# ax2.set_xlim([x2[0],x2[-1]])
ax1.set_ylim([0,1])
ax2.set_ylim([0,1])
# ax2.set_xscale('log')
#ax2.set_xlabel ('cavity linewidth [GHz]', fontsize = 15)
ax1.legend(loc = 'upper right', title = 'd = 3 $\mu$m; L = 2 $\mu$m; Q = 31000 ',fontsize=14)
ax2.legend(loc = 'upper left', title = 'd = 3 $\mu$m; L = 2 $\mu$m; T = 5 K', fontsize=14)
fig.savefig ('K://ns/qt/Diamond/Projects/Cavities/simulations/simulResults/emission_frac_vs_'+sweep_param1+sweep_param2+'.pdf')
fig.savefig ('K://ns/qt/Diamond/Projects/Cavities/simulations/simulResults/emission_frac_vs'+sweep_param1+sweep_param2+'.png')
plt.show()

#plt.plot(x, 1+F)
#plt.title('Freespace/cavity emitter lifetime reduction', fontsize = 17)
#plt.xlabel('Quality Factor Q', fontsize = 15)
#plt.ylabel('Purcell factor', fontsize = 15 )
#plt.show()