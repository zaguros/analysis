##############################################
# SvD, May 2017
##############################################
# Simulation of the resonant modes of a cavity with a diamond membrane

import numpy as np
import scipy as sc
import scipy.constants
import scipy.linalg
import math 
import os

from analysis.lib.tools import toolbox as tb
from matplotlib import pyplot as plt
import analysis.scripts.cavity.simulations_gaussian_beams as sim_gb; reload(sim_gb)
import analysis.scripts.cavity.simulations_vibrations as sim_vib; reload(sim_vib)
import analysis.lib.cavity.simulations_cavity as sim_cav; reload(sim_cav)



data_dir = 'K://ns/qt/Diamond/Projects/Cavities/simulations/transfer_matrix_model/'
if os.path.isdir(data_dir):   
    data_folder = os.path.join(data_dir,tb.get_timestamp_from_now()[0:8])
    if os.path.isdir(data_folder):
        print 'saving data in ', data_folder
    else:
        print 'creating directory', data_folder
        os.mkdir(data_folder)
else:
    print 'data_dir:  '+data_dir+ ' does not exist! Saving disabled.'


n_air = 1.0
n_diamond = 2.41


class Cavity():

    def __init__(self, t_d, t_a_g, R, **kw):
        self.t_d = t_d
        self.t_a_g = t_a_g
        self.R = R
        self.lambda_cav = kw.pop('lambda_cav', 637.e-9)
        self.n_1 = kw.pop('n_1', 2.15) #(Ta2O5)
        self.n_2 = kw.pop('n_2', 1.46) #SiO2
        self.M = kw.pop('M', 10) #number of mirror layers
        self.cav_type = kw.pop('cav_type','hybrid') #cavity type. can be hybrid or air
        self.lambda_i = kw.pop('lambda_i', 637.e-9)
        self.res_search_range = kw.pop('res_search_range',100e-9)
        self.res_wl_search_range = kw.pop('res_wl_search_range', 0.5e-9)
        self.dnu = kw.pop('dnu',1.e9)
        self.beta0 = kw.pop('beta0',0.03) #free space branching ratio into ZPL
        self.tau0 = kw.pop('tau0',12) #excited state lifetime in ns

        if self.cav_type == 'hybrid':
            self.n_z0 = n_diamond
        elif self.cav_type == 'air':
            self.n_z0 = n_air
            self.t_d= 0
        elif self.cav_type == 'test1':
            print 'cav type is tst1'
            self.n_z0 = n_air      
            self.M1 = kw.pop('M1', 10) #number of mirror layers, first mirror 
            self.M2 = kw.pop('M2', 10) #number of mirror layers, second mirror
            self.t_s = kw.pop('t_s',5e-6)
            self.t_d = 0
        else:
            self.n_z0 = 0
            'specify valid cavity type (hybrid or air)'

        self.ns_in_cavity()

    def boundary_matrix(self,n1,n2):
        rho = (n1-n2)/(n1+n2)
        tau = 1+rho
        matrix = (1./tau)*np.array(((1,rho),(rho,1)))
        return matrix

    def propagation_matrix(self,n,t):
        k_c = 2*math.pi*n/self.lambda_i
        matrix = np.array(((np.exp(1j*k_c*t),0),(0,np.exp(-1j*k_c*t))))
        return matrix

    def electric_field_from_transfer_matrix(self,matrix,E):
        out = np.dot(matrix,E)
        return out  

    def ns_in_cavity(self):
        if self.cav_type == 'hybrid':
            ns_mirror = np.append(np.array((self.n_1,self.n_2)*self.M),np.array((self.n_1)))
            ns_cav = np.array((n_diamond, n_air))
            self.ns = np.concatenate((ns_mirror,ns_cav,ns_mirror))
        elif self.cav_type == 'air':
            ns_mirror = np.append(np.array((self.n_1,self.n_2)*self.M),np.array((self.n_1)))
            ns_cav = np.array([n_air])
            self.ns = np.concatenate((ns_mirror,ns_cav,ns_mirror)) 
        elif self.cav_type == 'test1':
            ts_substrate = np.array([n_diamond])
            ns_mirror1 = np.append(np.array((self.n_1,self.n_2)*self.M1),np.array((self.n_1)))
            ns_mirror2 = np.append(np.array((self.n_1,self.n_2)*self.M2),np.array((self.n_1)))
            ns_cav = np.array([n_air])
            self.ns = np.concatenate((ts_substrate,ns_mirror1,ns_cav,ns_mirror2))  
        else:
            'specify valid cavity type (hyrbid or air)!'
            return
        return self.ns

    def ts_in_cavity(self):
        if self.cav_type == 'hybrid':
            ts_mirror = np.append(np.array((self.lambda_cav/(4*self.n_1),self.lambda_cav/(4*self.n_2))*self.M),np.array(self.lambda_cav/(4*self.n_1)))
            ts_cav = np.array((self.t_d,self.t_a))
            self.ts = np.concatenate((ts_mirror,ts_cav,ts_mirror))
        elif self.cav_type == 'air':
            ts_mirror = np.append(np.array((self.lambda_cav/(4*self.n_1),self.lambda_cav/(4*self.n_2))*self.M),np.array(self.lambda_cav/(4*self.n_1)))
            ts_cav = np.array([self.t_a])
            self.ts = np.concatenate((ts_mirror,ts_cav,ts_mirror))
        elif self.cav_type == 'test1':
            ts_substrate = np.array([self.t_s])
            ts_mirror1 = np.append(np.array((self.lambda_cav/(4*self.n_1),self.lambda_cav/(4*self.n_2))*self.M1),np.array(self.lambda_cav/(4*self.n_1)))
            ts_mirror2 = np.append(np.array((self.lambda_cav/(4*self.n_1),self.lambda_cav/(4*self.n_2))*self.M2),np.array(self.lambda_cav/(4*self.n_1)))
            ts_cav = np.array([self.t_a])
            self.ts = np.concatenate((ts_substrate,ts_mirror1,ts_cav,ts_mirror2))        

        else:
            'specify valid cavity type (hyrbid or air)!'
            return
        return self.ts


    def cavity_reflectivity(self,na=n_air,nb=n_air):
        """
        calculates the cavity reflectivity, assuming there is no outgoing field.
        input:
        ts          lengths of all elements. right-most element first
        ns          refractive indices of all elements. right-most element first 
        na          refractive index of left-most element
        nb          refractive index of right-most element
        lambdas     wavelengths at which to evaluate
        """
        m1 = self.boundary_matrix(self.ns[0],na)
        mf = self.boundary_matrix(nb,self.ns[-1])

        Er = np.transpose(np.array((1,0)))

        for i,(t,n) in enumerate(zip(self.ts,self.ns)):
            if i==0:
                c = np.identity(2)
                m = m1
                E=Er
            else:
                m = self.boundary_matrix(self.ns[i],self.ns[i-1])
            c = np.dot(m,c)
     
            pi = self.propagation_matrix(n,t)
            c = np.dot(pi,c)

        c = np.dot(mf,c)
        r = c[1,0]/c[0,0] #reflectivity

        return r

    def find_res_wavelength(self,nr_pts=161,plot_r=False):
        """
        find the resonance condition by varying lambda around lambda_i, given cavity lengths.
        uses two consecutive optimisation steps
        """
        self.t_a = self.t_a_g 
        ts = self.ts_in_cavity()

        lambda_is = np.linspace(self.lambda_i-(self.res_wl_search_range),self.lambda_i+self.res_wl_search_range,nr_pts)
        r_is = np.zeros(len(lambda_is))
        for i,lambda_ii in enumerate(lambda_is):
            self.lambda_i = lambda_ii
            r = self.cavity_reflectivity()
            r_is[i]=np.abs(r)**2

        lambda_ii=lambda_is[np.argmin(r_is)]
        r_i = r_is[np.argmin(r_is)]
        
        if plot_r:
            fig,ax = plt.subplots()
            ax.plot(lambda_is*1.e6,r_is)
            plt.show()
            plt.close()
        
        lambda_is = np.linspace(lambda_ii-self.res_wl_search_range/100.,lambda_ii+self.res_wl_search_range/100.,2*nr_pts)
        r_iis = np.zeros(len(lambda_is))
        for i,lambda_iii in enumerate(lambda_is):
            self.lambda_i = lambda_iii
            r = self.cavity_reflectivity()
            r_iis[i]=np.abs(r)**2

        self.lambda_i = lambda_is[np.argmin(r_iis)]
        self.r = r_iis[np.argmin(r_iis)]
        if plot_r:
            print 'air gap',self.t_a*1.e6, ' um; r = ', self.r
            fig,ax = plt.subplots()
            ax.plot(lambda_is*1.e9,r_iis)
            plt.show()
            plt.close()
        
        self.optical_length = self.t_a + self.t_d*n_diamond
        return self.lambda_i,self.r,r_iis#,t_ai,r_i,r_is#,


    def find_res_condition(self,nr_pts=161,plot_r=False):
        """
        find the resonance condition by varying t_a around t_a_g
        uses two consecutive optimisation steps
        """
        
        t_as = np.linspace(self.t_a_g-(self.res_search_range),self.t_a_g+self.res_search_range,nr_pts)
        r_is = np.zeros(len(t_as))
        for i,t_ai in enumerate(t_as):
            self.t_a = t_ai
            ts = self.ts_in_cavity()
            r = self.cavity_reflectivity()
            r_is[i]=np.abs(r)**2

        t_ai=t_as[np.argmin(r_is)]
        r_i = r_is[np.argmin(r_is)]
        
        if plot_r:
            fig,ax = plt.subplots()
            ax.plot(t_as*1.e6,r_is)
            plt.show()
            plt.close()
        
        t_as = np.linspace(t_ai-self.res_search_range/100.,t_ai+self.res_search_range/100.,2*nr_pts)
        r_iis = np.zeros(len(t_as))
        for i,t_aii in enumerate(t_as):
            self.t_a = t_aii
            ts = self.ts_in_cavity()
            r = self.cavity_reflectivity()
            r_iis[i]=np.abs(r)**2

        self.t_a = t_as[np.argmin(r_iis)]
        self.r = r_iis[np.argmin(r_iis)]
        if plot_r:
            print 'air gap',self.t_a*1.e6, ' um; r = ', self.r
            fig,ax = plt.subplots()
            ax.plot(t_as*1.e6,r_iis)
            plt.show()
            plt.close()
        
        self.optical_length = self.t_a + self.t_d*n_diamond
        return self.t_a,self.r,r_iis#,t_ai,r_i,r_is#,

    def calculate_w0(self):
        """
        calculate the beam waist of the beam in diamond. Use Gaussian beam optics.
        parameters:
        La  optical cavity length
        d   diamond thickness
        R  radius of curvature
        lambda_i  wavelength in cavity
        n_diamond
        ---> rather: assume gaussian beams, calculated currently in my mathematica script. Use tables created with this.
        """
        #print 'looking up w0 from : La,d,R = ', self.t_a,self.t_d,self.R
        self.w0 = sim_gb.get_w0_from_table(La=self.t_a,d=self.t_d,R=self.R)
        #self.w0 =  1.04e-6#*math.from my calculations of beam waist in cavity Riedel
        #self.w0 = 0.704e-6*2#from FWHM quted in Riedels paper
        # self.w0 = 1.61e-6#from calculations for d=12, L = 4, R=39
        return self.w0

    def electric_field_distribution(self,na=n_air,nb=n_air,lambda_i = 637.e-9,nr_points=50001):
        """
        calculates the electric field disrtibution, assuming there is no leftmoving field in the transmission direction of the cavity.
        input:
        ts          lengths of all elements. right-most element first
        ns          refractive indices of all elements. right-most element first 
        na          refractive index of left-most element
        nb          refractive index of right-most element
        lambdas     wavelengths at which to evaluate
        """
        L=np.sum(self.ts)
        self.zs = np.linspace(0,L,nr_points)
        m1 = self.boundary_matrix(self.ns[0],na)
        mf = self.boundary_matrix(nb,self.ns[-1])
        self.ns_vs_z = np.zeros(nr_points)
        Efw_vs_z = np.zeros((nr_points),dtype='complex')
        Ebw_vs_z= np.zeros((nr_points),dtype='complex')
        
        t_tot=0

        Er = np.transpose(np.array((1,0),dtype='complex'))

        for i,(t,n) in enumerate(zip(self.ts,self.ns)):
            ki =  2*math.pi*n/self.lambda_i
            if i==0:
                c = np.identity(2)
                m = m1
                E=Er
            else:
                m = self.boundary_matrix(self.ns[i],self.ns[i-1])
            c = np.dot(m,c)
            E = self.electric_field_from_transfer_matrix(m,E)

            i_ts = np.where(((self.zs <= L-t_tot) & (self.zs > L-t_tot-t)))

            self.ns_vs_z[i_ts] = n
            Efw_vs_z[i_ts] = E[0]*np.exp(1j*(L-t_tot-self.zs[i_ts])*2*math.pi*n/self.lambda_i)#E2_fw*exp(1j*(L-t_tot-zs[i_t])*2.*math.pi*n/lambda_i)
            Ebw_vs_z[i_ts] = E[1]*np.exp(-1j*(L-t_tot-self.zs[i_ts])*2*math.pi*n/self.lambda_i)#E2_bw*exp(-1j*(L-t_tot-zs[i_t])*2.*math.pi*n/lambda_i)

            p_i = self.propagation_matrix(n,t)
            c = np.dot(p_i,c)

            E = self.electric_field_from_transfer_matrix(p_i,E)


            t_tot=t_tot+t

        c = np.dot(mf,c)
        
        self.Etot_vs_z = Efw_vs_z+Ebw_vs_z #E
        H_vs_z = (Efw_vs_z-Ebw_vs_z)*n #H
        Z_vs_z = self.Etot_vs_z/ H_vs_z# 

        self.r = c[1,0]/c[0,0] #reflectivity
        self.Etot_vs_z=np.abs(self.Etot_vs_z/c[0,0])
        H_vs_z=np.abs(H_vs_z/c[0,0])

    def calculate_E_max(self):
        """
        calculate the maximum electric field density and its location, in the cavity region with refractive index n_z0
        """
        dz = np.abs(self.zs[0]-self.zs[1])
        arg_zs_in_cavity = np.where((self.ns_vs_z<self.n_z0+0.001)&(self.ns_vs_z>self.n_z0-0.001))
        zs_in_cavity = self.zs[arg_zs_in_cavity]
        arg_shallow_zs_in_cavity = arg_zs_in_cavity[0][-int(self.lambda_i/2/self.n_z0/dz):]
        z0_i = arg_shallow_zs_in_cavity[0] + np.argmax(self.Etot_vs_z[arg_shallow_zs_in_cavity])
        self.z0 = self.zs[z0_i]
        self.E_z0 = self.Etot_vs_z[z0_i]
        self.dz0 = np.abs(zs_in_cavity[-1]-self.z0)
        
    def calculate_energy_dist_length(self):
        """
        calculate energy distribution length int((n(z)^2)*(E(z))^2 dz)/(n(z0))^2*E(z0)^2
        """
        dz = np.abs(self.zs[0]-self.zs[1])
        self.effective_length = np.real(np.sum((self.ns_vs_z**2)*(self.Etot_vs_z**2))*dz/((self.n_z0**2)*(self.E_z0**2)))


    def calculate_mode_volume(self):
        self.mode_volume = math.pi*self.w0**2/2.*self.effective_length #note that since we are using the effective length, we here have to use pi*wo^2/2
        #print 'mode vol',self.mode_volume
        
    def calculate_max_Evac(self):
        """
        Uses the mode volume to determine the maximum vacuum electric field (hbar omega/(2epsilon*V)) (epsilon = n^2*epsilon_0)
        """
        self.Evac_max = np.sqrt(sc.constants.hbar*2*math.pi*(sc.constants.c/self.lambda_i)/(2*sc.constants.epsilon_0*self.n_z0**2*self.mode_volume))

    def calculate_energy_in_m(self,n_m):
        """
        calculate the electric field energy (J/m^2) in cavity material with refractive index n_m
        """
        arg_zs_in_m = np.where((self.ns_vs_z<n_m+0.01)&(self.ns_vs_z>n_m-0.01))
        ns_vs_z_in_m = self.ns_vs_z[arg_zs_in_m]
        Etot_vs_z_in_m = self.Etot_vs_z[arg_zs_in_m]
        dz = np.abs(self.zs[0]-self.zs[1])
        Etot = 1./2*sc.constants.epsilon_0*np.sum((ns_vs_z_in_m**2)*(Etot_vs_z_in_m**2))*dz 
        return Etot #in J/m^2 

    def calculate_Purcell(self):
        """
        Calculate the purcell factor.
        """
        self.Fp = 3*scipy.constants.c*self.lambda_i**2/(4*math.pi**2*self.n_z0**3*self.dnu)/self.mode_volume
        return self.Fp

    def calculate_into_ZPL(self):
        self.intoZPL = self.beta0*(self.Fp+1)/(self.beta0*(self.Fp)+1) #into ZPL AND into cavity (into ZPL alone: beta*(F+1)/(beta*F+1))
        return self.intoZPL

    def calculate_lifetime(self):
        """
        Calculate excited state lifetime after Purcell enhancement of ZPL. 
        Use: gamma_new = (1+Fp)gamma_ZPL+gamma_PSB  = (1+Fp)beta0*gamma_tot+(1-beta0)*gamma_tot=(Fp*beta0+1)gamma_tot
        since:
        gamma_tot = gamma_ZPL+gamma_PSB = 1/(tau_tot)
        gamma_ZPL = beta0*gamma_tot
        """
        self.tau = 1./(self.Fp*self.beta0+1)*self.tau0
        return self.tau

    def analyse_cavity(self,nr_points=30001,plot_Evac=False,save_plot=True):

        self.ts_in_cavity()
        self.electric_field_distribution()
        self.calculate_E_max()
        self.calculate_energy_dist_length()

        self.calculate_mode_volume()#t_a+t_d*2.4)
        self.calculate_max_Evac()

        self.Etot_vs_z = self.Etot_vs_z*self.Evac_max/self.E_z0

        self.Etot_diamond = self.calculate_energy_in_m(n_diamond)
        self.Etot_air = self.calculate_energy_in_m(n_air)

        self.calculate_Purcell()
        self.calculate_into_ZPL()

        if plot_Evac:
            print 'energy distribution length = ',self.effective_length*1.e6, 'um'
            print 'beam waist', self.w0*1.e6, 'um'
            print 'mode volume', self.mode_volume*(1.e6)**3, 'um^3 = ', self.mode_volume/(self.lambda_i**3), 'lambda^3'
            print 'max Evac',self.Evac_max/1000.,'kV/m' 
            today = tb.get_timestamp_from_now()[:8]
            title_string = today+'_td_%.2f_ta_%.2f_lambdai_%.1f_R_%.1f'%(self.t_d*1.e6,self.t_a*1.e6,self.lambda_i*1.e9,self.R*1.e6)

            fig,ax = plt.subplots(figsize=(12,6))
            ax.plot(self.zs*1.e6,self.ns_vs_z)
            ax.set_ylabel('n')
            ax.set_xlabel('z (um)')
            ax.set_title(title_string)
            # ax.set_ylim((1.,2.6))
            ax2=ax.twinx()
            # ax2.set_xlim((6.2,6.5))
            # ax2.set_ylim((0,1))
            ax2.plot(self.zs*1.e6,self.Etot_vs_z/1000.,'g')
            ax2.plot([self.z0*1.e6],[self.Evac_max/1000.],'ok',label='dz = %.1f nm, Emax = %.1f kV/m'%(self.dz0*1.e9,self.Evac_max/1000))
            ax2.set_ylabel('E (kV/m)')
            ax2.legend()

            if save_plot:
                plt.savefig(data_folder+'/'+title_string+'.png')
            plt.show()
            plt.close()

        return self.zs,self.ns_vs_z,self.Etot_vs_z,self.Evac_max,self.z0,self.dz0,self.Etot_diamond,self.Etot_air


def calculate_Fp_vs_finesse(t_d=4.e-6,t_a_g=1.2e-6,R=18.e-6,finesses=np.linspace(100,20000,199)):
    s = Cavity(t_d, t_a_g, R)
    s.find_res_condition(plot_r=True)
    s.calculate_w0()
    s.analyse_cavity(plot_Evac=True)

    dnus = np.zeros(len(finesses))
    Fps = np.zeros(len(finesses))
    intoZPLs = np.zeros(len(finesses))
    taus = np.zeros(len(finesses))

    for i,f in enumerate(finesses):
        dnu = scipy.constants.c/(2*s.optical_length)/f #use the optical length here, since this will give the visible FSR.
        s.dnu = dnu
        if f == 5000:
            print 'du',dnu*1.e-9
        dnus[i] = s.dnu
        Fps[i] = s.calculate_Purcell()
        intoZPLs[i] = s.calculate_into_ZPL()
        taus[i] = s.calculate_lifetime()

    fig,ax = plt.subplots()
    ax.plot(finesses,taus,'magenta',label='tau (ns)')
    ax.set_ylabel('tau (ns)')
    ax.yaxis.label.set_color('magenta')
    ax.tick_params(axis='y', colors='magenta')
    ax2 = ax.twinx()
    ax2.plot(finesses,intoZPLs,'orange',label='%% into ZPL')
    ax2.yaxis.label.set_color('orange')
    ax2.tick_params(axis='y', colors='orange')
    ax2.spines['right'].set_color('orange')
    ax2.spines['left'].set_color('magenta')
    ax.set_xlabel('Finesse')
    ax2.set_ylabel('% into ZPL into cavity')
    ax.set_xscale('log')
    ax2.set_xscale('log')
    ax2.set_xlim((min(finesses),max(finesses)))

    plt.show()

    fig,ax = plt.subplots()
    ax.plot(dnus*1.e-9,Fps,'magenta',label='Fp')
    ax.set_ylabel('Fp')
    ax.set_xlabel('linewidth (GHz)')
    ax.yaxis.label.set_color('magenta')
    ax.tick_params(axis='y', colors='magenta')
    ax2 = ax.twinx()
    ax2.plot(dnus*1.e-9,intoZPLs,'orange',label='%% into ZPL')
    ax2.yaxis.label.set_color('orange')
    ax2.tick_params(axis='y', colors='orange')
    ax2.spines['right'].set_color('orange')
    ax2.spines['left'].set_color('magenta')
    ax2.set_ylabel('% into ZPL')
    ax2.set_xlim((min(dnus*1.e-9),6))

    plt.show()
    return dnus,Fps,intoZPLs

#this doesn't work, since the cavity resonance is more subtley dependent on the cavity length:
#it depends on the phase of the Efield at the diamond-air transition. As of course it does: 
#this gives us the wobbly behaviour of the cavity modes.
# def spectral_overlap(p_cav_length,cav_lengths,lambda_ZPL,dnu):
#     cav = sim_cav.CavitySims()
#     lambda_cavs = []
#     Qs=[]
#     for l in cav_lengths:
#         cav.optical_length = l
#         cav.wavelength_ZPL = lambda_ZPL
#         cav.calc_longitudinal_cav_modes(reference = 'cavity_length')
#         i_cav,lambda_cav,freq_cav = cav.find_nearest_cav_mode()
#         lambda_cavs.append (lambda_cav)
#         Qs.append(freq_cav/dnu)

#     print lambda_cavs
#     spectral_overlap = 1./(1.+4*(Qs**2)*((lambda_ZPL/lambda_cavs)-1.)**2)

def spectral_overlap(lambda_ZPL,lambda_cav,dnu):
    Q = sc.constants.c/lambda_cav/dnu#nu/dnu
    spectral_overlap = 1./(1.+4*(Q**2)*((lambda_ZPL/lambda_cav)-1.)**2)
    return spectral_overlap

def calculate_Fp_w_vibrations(t_d=4.e-6,t_a_g=1.2e-6,R=18.e-6,vib_dL=0.4e-9,dnu=3.5e9,lambda_ZPL=637.e-9,dLmax=1.e-9,nr_pts=51,show_plots=False,fast_approx=False): #dL in nm
    """
    Calculates the average emission into the ZPL, when taking vibrations into account.
    Assumes a gaussian distribution of cavity length with FWHM of vib_dL as a result of vibrations

    """
    s = Cavity(t_d, t_a_g, R,dnu=dnu,lambda_i=lambda_ZPL)
    s.find_res_condition(plot_r=False)
    t_a0 = s.t_a


    s.calculate_w0()
    s.analyse_cavity(plot_Evac=False)    
    if fast_approx:
        s.calculate_Purcell()
        intoZPL0 = s.calculate_into_ZPL()

    dLs = np.linspace(-dLmax,dLmax,nr_pts)
    p_cav_length = sim_vib.gaussian(t_a0+dLs,t_a0,vib_dL)

    into_ZPLs = np.zeros(len(dLs))
    lambda_ress = np.zeros(len(dLs))
    Purcells=np.zeros(len(dLs))

    for i,dL in enumerate(dLs):
        cav_lengths = np.linspace(s.optical_length-4*dL,s.optical_length+4*dL,nr_pts)
        
        si = Cavity(t_d, t_a0+dL, R,dnu=dnu)

        si.find_res_wavelength(plot_r=False)
        lambda_ress[i] = si.lambda_i

        if fast_approx:
            intoZPLonR = intoZPL0
        else:
            si.calculate_w0()
            si.analyse_cavity()
            Purcells[i] = si.calculate_Purcell()
            intoZPLonR = si.calculate_into_ZPL()

        into_ZPL = intoZPLonR  *spectral_overlap(lambda_ZPL,si.lambda_i,si.dnu)
        into_ZPLs[i] = into_ZPL

    avg_p_ZPL = sim_vib.avg_p_ZPL_to_zero(into_ZPLs,p_cav_length)

    if show_plots:
        fig,ax = plt.subplots()
        ax.plot(dLs, p_cav_length/np.sum(p_cav_length))
        ax.set_xlabel('dL w.r.t. resonance (nm)')
        ax.set_ylabel('probability distribution')
        plt.show()
        plt.close()

        fig,ax = plt.subplots()
        ax.plot(dLs*1.e9,into_ZPLs)
        ax.set_xlabel('dL w.r.t. resonance (nm)')
        ax.set_ylabel('into ZPL')
        plt.show()
        plt.close()

        fig,ax = plt.subplots()
        ax.plot(dLs*1.e9,lambda_ress*1.e9)
        ax.set_xlabel('dL w.r.t. resonance (nm)')
        ax.set_ylabel('resonant wavelength')
        plt.show()
        plt.close()

        fig,ax = plt.subplots()
        ax.plot(dLs*1.e9,Purcells)
        ax.set_xlabel('dL w.r.t. resonance (nm)')
        ax.set_ylabel('Purcell factor')
        plt.show()
        plt.close()

    return into_ZPLs,Purcells,lambda_ress,avg_p_ZPL

def calc_avg_ZPL_vs_vibrations(t_d=4.e-6,t_a_g=1.2e-6,R=18.e-6,vib_dLs=np.linspace(0.05,1.,20)*1.e-9,dnu=3.5e9,lambda_ZPL=637.e-9,dLmax=4.e-9,nr_pts=3):
    """
    calculate average emission into ZPL into cavity vs FWHM of vibrations
    """
    avg_ZPLs = np.zeros(len(vib_dLs))
    for i,vib_dL in enumerate(vib_dLs):
        print i, 'out of ', len(vib_dLs), 'done'
        print vib_dL
        a,b,c,avg_ZPL = calculate_Fp_w_vibrations(t_d=t_d,t_a_g=t_a_g,R=R,vib_dL=vib_dL,dnu=dnu,lambda_ZPL=lambda_ZPL,dLmax=dLmax,nr_pts=nr_pts,fast_approx=True)
        print avg_ZPL
        avg_ZPLs[i]=avg_ZPL

    fig,ax = plt.subplots()
    ax.plot(vib_dLs,avg_ZPLs)
    ax.set_xlabel('FWHM vibrations (nm)')
    ax.set_ylabel('avg emission into ZPL')
    plt.show()
    plt.close()

    return avg_ZPLs

if __name__ == '__main__':
    R=18.e-6
    t_d = 4.e-6
    t_a_g = 1.2e-6
    s = Cavity(t_d, t_a_g, R)
    s.find_res_condition(plot_r=True)
    s.calculate_w0()
    s.analyse_cavity(plot_Evac=True)

 





def dielmirror_recursive_refl_coeff(Gamma_ip1, k_i, l_i,rho_i):
    Gamma_i = (rho_i+Gamma_ip1*np.exp(-2j*k_i*l_i) )/(1+rho_i*Gamma_ip1*np.exp(-2j*k_i*l_i) )
    return Gamma_i

def dielmirror_reflectivity(n_1,n_2,M,lambda_design=637.e-9,na=n_air,nb=n_air,lambdas = np.array([637.e-9])):
    """
    calculates dielectric mirror reflectivity optimised for wavelength 'lambda', using recursive application of reflection coefficient.
    Inputs: 
    n_1         refractive index start & end layer
    n_2         refractive index alternating layer
    M           number of alternating layers - total # layers is 2M+1
    lambda_design   wavelength for which mirror is designed (layers are lambda/4/n thick)
    na          refractive index outside layer (default: n_air)    
    nb          refractive index outside layer (default: n_air)
    lambdas     array of wavelenths at which to evaluate the reflectivity
    """

    rho_Mp1 = (n_1 - nb)/(n_1+nb)
    rho_1 = (na - n_1) /(na+n_1)
    rho_odd = (n_2 - n_1)/(n_2+n_1)
    rho_even = (n_1 - n_2)/(n_2+n_1)

    L_1 = lambda_design / (4*n_1)
    L_2 = lambda_design / (4*n_2)

    ks_1 = 2*math.pi*n_1/lambdas
    ks_2 = 2*math.pi*n_2/lambdas

    # initiliaze the recursion using the right-most interface (E'-_Mp1 = 0)
    Gamma_Mp1 = np.ones(len(lambdas))*rho_Mp1

    Gamma = Gamma_Mp1  

    for n in np.arange(M):
        Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_1, L_1,rho_odd)
        Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_2, L_2,rho_even)

    # and the final layer, back to material na
    Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_1,L_1, rho_1)
    Refl = np.conjugate(Gamma)*Gamma
    return Gamma,Refl


def dielmirror_narrowsb_reflectivity(n_1,n_2,M,lambda_design=637.e-9,na=n_air,nb=n_air,lambdas = np.array([637.e-9])):
    """
    calculates dielectric mirror reflectivity optimised for wavelength 'lambda', using recursive application of reflection coefficient.
    Inputs: 
    n_1         refractive index start & end layer
    n_2         refractive index alternating layer
    M           number of alternating layers - total # layers is 2M+1+4 
    lambda_design   wavelength for which mirror is designed (layers are lambda/4/n thick)
    na          refractive index outside layer (default: n_air)    
    nb          refractive index outside layer (default: n_air)
    lambdas     array of wavelenths at which to evaluate the reflectivity
    """

    rho_Mp1 = (n_1 - nb)/(n_1+nb)
    rho_1 = (na - n_1) /(na+n_1)
    rho_odd = (n_2 - n_1)/(n_2+n_1)
    rho_even = (n_1 - n_2)/(n_2+n_1)

    L_1 = lambda_design / (4*n_1)
    L_2 = lambda_design / (4*n_2)

    ks_1 = 2*math.pi*n_1/lambdas
    ks_2 = 2*math.pi*n_2/lambdas

    # initiliaze the recursion using the right-most interface (E'-_Mp1 = 0)
    Gamma_Mp1 = np.ones(len(lambdas))*rho_Mp1

    Gamma = Gamma_Mp1  

    Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_1, 1.27335*L_1,rho_odd)
    Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_1, 2.94525*L_2,rho_even)
    Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_1, 0.94808*L_1,rho_odd)
    Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_1, 2.99742*L_1,rho_even)


    for n in np.arange(M):
        Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_1, L_1,rho_odd)
        Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_2, 3*L_2,rho_even)

    # and the final layer, back to material na
    Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_1,L_1, rho_1)
    Refl = np.conjugate(Gamma)*Gamma
    return Gamma,Refl


def calculate_transfer_matrices_mirrors(wavelength, lambda_design=637e-9,n_air=n_air,n_H=2.15,n_L=1.46,na=n_air,nb=n_air,number_of_layers=10):
    k_H = 2.*math.pi*n_H/wavelength
    k_L = 2.*math.pi*n_L/wavelength
    t_H = lambda_design/(4*n_H)
    t_L = 3*lambda_design/(4*n_L)
    v_H =np.exp(1j*k_H*t_H)


    rho = (n_H-n_L)/(n_H+n_L)
    rho_in = (na-n_H)/(na+n_H)
    rho_out = (n_H-nb)/(n_H+nb)

    B1 = 1./(1+rho)*np.array(((1,rho),(rho,1)))
    P1 = np.array(((np.exp(1j*k_H*t_H),0),(0,np.exp(-1j*k_H*t_H))))
    B2 = 1./(1-rho)*np.array(((1,-rho),(-rho,1)))
    P2 = np.array(((np.exp(1j*k_L*t_L),0),(0,np.exp(-1j*k_L*t_L))))

    F = np.dot(B2,P2)
    F = np.dot(P1,F)
    F = np.dot(B1,F)

    B_in = 1./(1+rho_in)*np.array(((1,rho_in),(rho_in,1)))
    P_in = np.array(((np.exp(1j*k_H*t_H),0),(0,np.exp(-1j*k_H*t_H))))
    F_in = np.dot(B_in,P_in)
    F_out = (1./(1.+rho_out))*np.array(((1.,rho_out),(rho_out,1.)))


    N = number_of_layers
    # add all layers of the mirror. 
    H = F_in
    for i in range(N):
        H = np.dot(H,F)

    M_gma = np.dot(H,F_out)
    M_amg = scipy.linalg.inv(M_gma)

    return(M_gma,M_amg)

def mirror_reflectivity_from_matrix(lambdas,number_of_layers=10):
    """
    Alternative method to get mirror reflectivity from total transfer matrix. 
    """

    r_coeff = np.zeros(len(lambdas),dtype='complex')
    R = np.zeros(len(lambdas),dtype='complex')

    for i,wavelength in enumerate(lambdas):
        M_gma,M_amg = calculate_transfer_matrices_mirrors(wavelength=wavelength,number_of_layers=number_of_layers)

        r_coeff[i] = M_gma[1,0]/M_gma[0,0]

        R[i] = np.abs(r_coeff[i])**2
        R[i] = np.conjugate(r_coeff[i])*r_coeff[i]
    return R,r_coeff





