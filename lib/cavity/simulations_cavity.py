"""Simulations of a cavity, based on Albrecht, ..., Becher et al.  Author: Cristian Bonato
Improvements to be made:
- parameter initialization more clear.
- save parameters with a figure.
- get the right ZPL - PSB branching ratio, lifetime data. (not for nanodiamonds)
"""


import numpy as np
import math
import pylab as plt
import matplotlib
from analysis.lib.tools import toolbox as tb
from matplotlib import rc, cm

n_diamond=2.419
c = 3.e8

class CavitySims ():

    def __init__(self):

        #NV parameters
        #experimental values taken from PRL 110, 243602
        #self.freq = np.array([469.6e12, 461.9e12,452.6e12, 439.7e12, 429.0e12, 416.2e12, 403.2e12, 392.4e12])
        #self.A = np.array([1520., 5260., 18600., 16400., 14000., 9180., 6570., 3270.]) ##branching ratios from Albrecht et al.
        #self.A = np.array([0.03, 0.1, 0.1, 0.2, 0.17, 0.2, 0.1, 0.1]) ##fake branching ratio, but with ZPL/(all)=0.03
        self.freq = np.array([471.1e12,456.3e12,440.6e12,432.4e12,423.4e12])
        self.A = np.array([119.,1648.,646.,439.,488.])
        self.linewidths = np.array([0.4e12,19.8e12,10.8e12,10.6e12,9.1e12])
        # self.freq = np.array([471.1e12,456.9e12,438.4e12,423.7e12])
        # self.A = np.array([120.,1463.,1560.,496.])
        # self.linewidths = np.array([0.4e12,19.2e12,19.8e12,10.0e12])

        self.d = self.freq - self.freq[0] 
        #self.freqC = freq[0]

        self.A_tot = np.sum(self.A)
        self.epsilon = self.A/self.A_tot #relative strength transitions
        print 'brancing ratio ZPL to rest:',self.epsilon[0]
        #self.gm_tot = 35.e6 #deduced by g2 msmnts on nanodiamonds (Albrecht)
        self.gm_tot = 1/(12.e-9) ####the above value is from the Albrecht  (becher) paper, but corresponds to t = 28 ns, in nanodiamond.. this is typically longer.. 
        print 'lifetime excited state:',str(round(1/self.gm_tot*1.e9,1)),'ns'
        self.g_relative = self.gm_tot*self.epsilon
        #CAVITY parameter
        self.cavity_length = 1.1e-6 # only air length
        self.Q = None
        self.T = None
        self.diamond_thickness = 3.e-6
        self.radius_curvature = 15.e-6
        self.wavelength_ZPL = c/float(self.freq[0])

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
        d = self.diamond_thickness*n_diamond + self.cavity_length
        if (d>self.radius_curvature-2e-6):
            print "Cavity is unstable!"
            ret = 1
        self.waist = ((self.wavelength_ZPL/np.pi)**0.5)*(d*(self.radius_curvature-d))**(1/4.)
        self.cavity_volume = np.pi*((0.5*self.waist)**2)*d
        return ret

    def set_temperature(self, value):
        self.T = value
        dephZPL = (16.2+9.2e-7*(1/0.0125)*self.T**5)*1e6 #in Hz. Only for T < 100 K!, removed 2*pi, compared to Cristian's code. 
        #self.linewidths = np.array([dephZPL, 15.9e12, 15.5e12, 15.0e12, 16.5e12, 12.7e12, 13.7e12, 16.1e12])
        self.linewidths[0] = dephZPL
        self.gm_star = 2*np.pi*self.linewidths[0] 

    def set_ZPL_linewidth (self, value):
        #self.linewidths = np.array([value, 15.9e12, 15.5e12, 15.0e12, 16.5e12, 12.7e12, 13.7e12, 16.1e12])
        self.linewidths[0] = dephZPL       
        self.gm_star = 2*np.pi*self.linewidths[0]
        self._linewidth = True

    def set_Q (self, value,verbose=False):
        self.Q = value
        self.k = math.pi*c/(self.wavelength_ZPL*self.Q)
        if verbose:
            print 'k',self.k/1.e9,'GHz'
        #print self.Q, self.k

    def calculate_params (self, do_plot = False, verbose=False):
        
        self.FSR = c/(2*(self.diamond_thickness*2.5+self.cavity_length)) #FSR
        self.resonant_frq = c/self.wavelength_ZPL

        #longitudinal cavity modes
        self.lambda_modes_nm=[]
        self.freq_modes = []
        l=0
        i=0
        while ((l<800) and (i< 1000)):
            nu = self.resonant_frq-i*self.FSR
            l = 1e9*c/nu 
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
        self.waist = ((self.wavelength_ZPL/np.pi)**0.5)*(d*(self.radius_curvature-d))**(1/4.)
        self.cavity_volume = np.pi*((0.5*self.waist)**2)*d

        if verbose:
            print 'lambda_ZPL = ', self.wavelength_ZPL*1e9, '[nm]'
            print 'FSR = ', 1e9*((self.wavelength_ZPL**2)/3e8)*self.FSR, '[nm]'
            print "mode volume: ", self.cavity_volume*1e18,  "micron^3"

        #T^5 dephasing ZPL, linewidths = w
        if not(self._linewidth):
            dephZPL = (16.2+9.2e-7*(1/0.0125)*self.T**5)*1e6 #in Hz. This is w - thus differs a factor 2pi from Jahn-Teller paper
            self.linewidths[0] = dephZPL 
            #self.linewidths = np.array([dephZPL, 15.9e12, 15.5e12, 15.0e12, 16.5e12, 12.7e12, 13.7e12, 16.1e12])
        self.gm_star = 2*np.pi*self.linewidths[0]
        self.g = np.sqrt((3.*c*(self.wavelength_ZPL**2)*self.gm_tot/2.)/(4.*math.pi*self.cavity_volume))
        self.gm_relax = 2*math.pi*self.linewidths-2*np.pi*self.linewidths[0]
        self.g_all = self.g*np.sqrt(self.epsilon)
        self.gamma=(self.k+self.gm_tot+self.gm_relax+self.gm_star)
        
        self.naive_FpA = 3. * self.Q * (self.wavelength_ZPL)**3. / n_diamond**3. /(4.*math.pi**2 * self.cavity_volume) 
        self.naive_FpD = 3. * self.Q * (self.wavelength_ZPL)**3. / (n_diamond**3 - n_diamond) /(4.*math.pi**2 * self.cavity_volume) 

        #print 'Q',self.Q
        #print self.g
        #print self.k
        #print 'g_relax',self.linewidths[0]
        #print 'lambdaZPL',self.wavelength_ZPL
        
        self.naive_Fp2= 4.*((self.g)**2)/(self.k*self.gm_tot)
        self.naive_Fp3 = 3. * self.Q * (self.wavelength_ZPL)**3. / (4.*(math.pi**2) * self.cavity_volume) 

        
    def simulate (self, do_plot = True, do_save = False, verbose=False):
        
        delta = np.zeros([len(self.linewidths), len(self.lambda_modes_nm)])
        
        i=-1
        for f_line in self.d:
            j = -1
            i = i+1
            for f_mode in self.freq_modes:
                j = j+1
                delta [i, j] = f_line - f_mode ###get delta_ij = det between NV-line i and cav-mode j

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
        self.P_tot = np.sum(self.P, axis=0) # sum over NV lines
        self.R_tot = np.sum(self.R, axis=1) ## the emission rate from the ZPL, PSBs in all modes (sum over cavity modes)
        self.R_plus_gm = np.add(self.R_tot, self.g_relative)
        self.P_PSB = np.sum(self.P[1:,:],axis=0) #sum over the PSB lines
        self.P_ZPL = self.P[0,:]
        self.P_PSB_zero = self.P_PSB[0]
        self.P_PSB_nonzero = np.sum(self.P_PSB[1:])
        self.P_ZPL_zero  = self.P_ZPL[0]
        self.P_ZPL_nonzero = np.sum(self.P_ZPL[1:])       
        self.lifetime = 1./np.sum(self.R_plus_gm)

        if do_plot:
            plt.figure()
            plt.plot (self.lambda_modes_nm, self.P_tot, 'ob', markersize=2)
            plt.ylabel ('probability of emission in the mode')
            plt.xlabel ('wavelength [nm]')
            plt.show()
            
    
    def emission_in_ZPL (self, sweep_param = 'Cavity length (um)', min_val=0, max_val=15e-6, nr_points=50, xlogscale=False):
        sweep_vals = np.linspace (min_val, max_val, nr_points)
        if xlogscale:
            sweep_vals = np.logspace(np.log10(min_val), np.log10(max_val), nr_points)

        emission_prob = np.zeros (nr_points)
        emission_prob_others = np.zeros(nr_points)
        emission_prob_ZPL_in_0 = np.zeros(nr_points)
        emission_prob_PSB_in_0 = np.zeros(nr_points)
        emission_prob_ZPL_nonzero = np.zeros(nr_points)
        emission_prob_PSB_nonzero = np.zeros(nr_points)
        purcell = np.zeros(nr_points)
        purcellA = np.zeros(nr_points)
        purcellD = np.zeros(nr_points)
        ind = 0
        mode_vol = np.zeros(nr_points)
        rate_ZPL = np.zeros(nr_points)
        rate_PSB = np.zeros(nr_points)
        gamma_ZPL = np.zeros(nr_points)
        gamma_PSB = np.zeros(nr_points)
        lifetime = np.zeros(nr_points)

        print 'Q',int(self.Q)
        print 'L_a',self.cavity_length
        print 'L_d',self.diamond_thickness
        for p in sweep_vals:
            err = 0
            if (sweep_param == 'Cavity length (um)'):
                err = self.set_cavity_length (p) #err is 1 if cavity is unstable
            elif (sweep_param == 'Temperature (K)'):
                self.set_temperature (p)
            elif (sweep_param == 'dephasing rate (GHz)'):
                self.set_ZPL_linewidth (p*1e9/(2*math.pi))
            elif (sweep_param == 'Quality factor'):
                self.set_Q (p)
        
            self.calculate_params(do_plot = False)
            self.simulate(do_plot=False,verbose=False)
            if (err==0):
                emission_prob[ind] = self.P_tot[0]
                emission_prob_others[ind] = np.sum(self.P_tot[1:])
                emission_prob_ZPL_in_0[ind] = self.P_ZPL_zero
                emission_prob_PSB_in_0[ind] = self.P_PSB_zero
                emission_prob_ZPL_nonzero[ind] = self.P_ZPL_nonzero
                emission_prob_PSB_nonzero[ind] = self.P_PSB_nonzero
                mode_vol[ind] = self.cavity_volume
                purcell[ind] = self.F0
                purcellA[ind] = self.naive_FpA
                purcellD[ind] = self.naive_FpD
                rate_ZPL[ind] = self.R_plus_gm[0]
                rate_PSB[ind] = np.sum(self.R_plus_gm[1:])
                gamma_ZPL[ind] = self.g_relative[0]/float(np.sum(np.sum(self.R))+self.gm_tot)
                gamma_PSB[ind] = np.sum(self.g_relative[1:])/float(np.sum(np.sum(self.R))+self.gm_tot)
                lifetime[ind] = self.lifetime
                ind = ind + 1

        return sweep_vals[:ind], emission_prob[:ind], mode_vol[:ind], purcell[:ind], purcellA[:ind],purcellD[:ind],emission_prob_others[:ind],emission_prob_ZPL_in_0[:ind],emission_prob_PSB_in_0[:ind],emission_prob_ZPL_nonzero[:ind],emission_prob_PSB_nonzero[:ind],rate_ZPL[:ind],rate_PSB[:ind], gamma_ZPL[:ind],gamma_PSB[:ind],lifetime[:ind]


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
    
    def plot_vs_sweepparam(self, sweep_param, min_val, max_val, nr_points, xlogscale=False, plotmode='all_to_zero'):
        x, y1, V, F_p, F_pA, F_pD, y2, y3, y4,y5,y6, y7, y8,y9,y10,y11 = self.emission_in_ZPL(sweep_param, min_val, max_val, nr_points, xlogscale)
            #y1 : emission all - into zero mode
            #y2 : emission all - into others modes
            #y3 : emission ZPL - into zero mode
        fig = plt.figure (figsize=(8,6))
        ax = fig.add_subplot(111)

        print plotmode
        alternative_ylabel=False
        if plotmode == 'all_to_zero':
            ax.plot(x, y1, color = 'DarkViolet', linewidth = 4, label = 'resonant mode') #0.5 factor from losses
        elif plotmode == 'all_to_rest':
            ax.plot(x, y2, color = 'Green', linewidth = 4, label = 'all other modes') 
        elif plotmode == 'ZPL_to_zero':
            ax.plot(x, y3, color = 'RoyalBlue', linewidth = 4, label = 'ZPL to resonant mode') 
        elif plotmode == 'ZPL_PSB_inout':
            ax.plot(x,y3, color = 'DarkViolet', linewidth = 4, label = 'ZPL to resonant mode')
            ax.plot(x,y4, color = 'Green', linewidth = 4, label = 'PSB to resonant mode')
            ax.plot(x,y5,'--', color = 'Darkviolet', linewidth = 4, label = 'ZPL to other modes')
            ax.plot(x,y6, '--',color = 'Green', linewidth = 4, label = 'PSB to other modes')
            ax.plot(x,y9, '.',color = 'DarkViolet', linewidth = 4, label = 'ZPL outside cavity')
            ax.plot(x,y10, '.',color = 'Green', linewidth = 4, label = 'PSB outside cavity')
            #ax.text(x[0]+1,y3[0]+0.02,'P00='+str(round(y3[0],2)),color='DarkViolet')
            #ax.text(x[-1]-12,y3[-1]+0.08,'P00='+str(round(y3[-1],2)),color='DarkViolet')
        elif plotmode == 'rates':
            ax.plot(x, y7/1.e6, color = 'DarkViolet', linewidth = 4, label = 'Rate ZPL to all') 
            ax.plot(x, y8/1.e6, color = 'Green', linewidth = 4, label = 'Rate PSB to all') 
            alternative_ylabel=True
            ylabel = 'Rate (MHz)'
        elif plotmode == 'branching_ratio':
            branchingratio = np.divide(y7,(np.add(y7,y8)))
            ax.plot(x, branchingratio, color = 'DarkViolet', linewidth = 4, label = 'branching ratio') 
            ax.text(x[0]+4000,branchingratio[0],str(round(branchingratio[0],2)),fontsize=15,color='DarkViolet')
            ax.text(x[-1]-30000,branchingratio[-1]-0.1,str(round(branchingratio[-1],2)),fontsize=15,color='DarkViolet')
            alternative_ylabel=True
            ylabel = 'probability of emission into ZPL'
            ax.set_ylim([0,1])
        elif plotmode == 'Fpurcell':
            ax.plot(x, F_p, color = 'DarkViolet', linewidth = 4, label = '$R_{0,0}/\gamma$') 
            ax.plot(x, F_pA, color = 'Green', linewidth = 4, label = '$F_p^{(A)}$') 
            ax.plot(x, F_pD, color = 'Orange', linewidth = 4, label = '$F_p^{(D)}$') 
            alternative_ylabel=True
            ylabel = 'Purcell factor'
        elif plotmode == 'lifetime':
            ax.plot(x, y11*1.e9, color = 'Orange', linewidth = 4, label = 'lifetime') 
            alternative_ylabel=True
            ylabel = 'lifetime(ns)'
            
        if xlogscale:
            ax.set_xscale('log')
        
        ax.set_xlabel(sweep_param, fontsize = 20)
        if not alternative_ylabel:
            ax.set_ylabel ('Probability of emission into channel', fontsize = 20)
            ax.set_ylim([0,1])
        else:
            ax.set_ylabel (ylabel, fontsize = 20)
        # ax.text(x[-1]-8000,y3[-1]-0.1,str(round(y3[-1],2)),fontsize=15,color='DarkViolet')
        ax.tick_params(axis='both', which='major', labelsize=16)
        # ax.legend(loc = 'right')
        ax.set_title('d = '+str(round(self.diamond_thickness*1e6,1))+' $\mu$m; L = ' +str(round(self.cavity_length*1e6,1)) +' $\mu$m; ROC = '+str(int(self.radius_curvature*1.e6))+' $\mu$m; Q ='+str(int(self.Q)) +'\n',fontsize=16)
        sweepparam = sweep_param.replace(" ","")
        
        ax.yaxis.set_ticks(np.arange(0, 1.2, 0.5))
        ax.tick_params(axis='both', which='major',   top='off',right = 'off',labelsize=14)
        ax.set_ylim(-0.02,1.02)
        plt.tight_layout()


        try:
            fig.savefig('H://My Documents/Cavities/Simulations/'+tb.get_timestamp_from_now()+plotmode+'_vs_'+sweepparam+'.pdf')
        except: 
            print 'Figure not saved'
        plt.show()

plt.close('all')

s = CavitySims()
s._linewidth=False

F=10000
L=2e-6
lambdaa=0.637e-6#0.6393e-6
#L=9./2*lambdaa

dnu = 10.e9
nu = c/637.e-9
Q = nu/dnu
d_diamond=3.e-6
ROC = 15.e-6

s.set_diamond_thickness (d_diamond)
s.set_mirror_curvature (ROC)
s.set_cavity_length (L)
s.set_Q (Q)

#s.calculate_params()
#s.simulate()

sweep_param1='Temperature (K)'
sweep_param2='Quality factor'
sweep_param3='dephasing rate (GHz)'

#s.plot_vs_sweepparam(sweep_param1, 0, 100, 100, plotmode='ZPL_PSB_inout')
#s.plot_vs_sweepparam(sweep_param1, 0, 100, 100, plotmode='branching_ratio')
#s.plot_vs_sweepparam(sweep_param1, 0, 100, 20, plotmode='rates')
#s.plot_vs_sweepparam(sweep_param1, 0, 100, 20, plotmode='Fpurcell')
#s.plot_vs_sweepparam(sweep_param1, 0, 100, 21, plotmode='ZPL_to_zero')
#s.plot_vs_sweepparam(sweep_param1, 0, 100, 51, plotmode='lifetime')

s.set_temperature (5) #after sweeping the temperature, set T to 5.
s.plot_vs_sweepparam(sweep_param2, 0, 400000, 100, plotmode='branching_ratio')
# s.plot_vs_sweepparam(sweep_param2, 44999, 45000, 3, plotmode='lifetime')
# print "t",s.lifetime
# print "p",s.P_ZPL_zero
# s.plot_vs_sweepparam(sweep_param2, 0, 100000, 100, plotmode='ZPL_PSB_inout')
#s.plot_vs_sweepparam(sweep_param2, 0, 200000, 500, plotmode='Fpurcell')
#s.plot_vs_sweepparam(sweep_param3, 0.00001, 10, 500, plotmode='Fpurcell',xlogscale=True)

