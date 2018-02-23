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
import time
import json
import h5py


from analysis.lib.tools import toolbox as tb
from matplotlib import pyplot as plt
import analysis.scripts.cavity.simulations_gaussian_beams as sim_gb; reload(sim_gb)
import analysis.scripts.cavity.simulations_vibrations as sim_vib; reload(sim_vib)
import analysis.lib.cavity.simulations_cavity as sim_cav; reload(sim_cav)

from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common

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




def default(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError('Not serializable')
    
def save_to_json_file(folder,filename,analysis_dict):
    f = open(os.path.join(folder,filename+'.json'), 'w')
    json.dump(analysis_dict,f,default=default)
    f.close()

n_air = 1.0
n_diamond = 2.41
n_Al2O3 = 1.77


class Cavity():

    def __init__(self, t_d, t_a_g, R, **kw):
        self.t_d = t_d
        self.t_a_g = t_a_g
        self.R = R
        self.lambda_cav = kw.pop('lambda_cav', 637.e-9)
        self.n_1 = kw.pop('n_1', 2.1349) #(Ta2O5)
        self.n_2 = kw.pop('n_2', 1.476) #SiO2
        self.M = kw.pop('M', 10) #number of mirror layers
        self.cav_type = kw.pop('cav_type','hybrid') #cavity type. can be hybrid or air
        self.lambda_i = kw.pop('lambda_i', 637.e-9)
        self.res_search_range = kw.pop('res_search_range',100e-9)
        self.res_wl_search_range = kw.pop('res_wl_search_range', 0.5e-9)
        self.dnu = kw.pop('dnu',3.e9)
        self.calculate_dnu =kw.pop('calculate_dnu',False)
        self.beta0 = kw.pop('beta0',0.03) #free space branching ratio into ZPL
        self.tau0 = kw.pop('tau0',12) #excited state lifetime in ns
        self.AR_coating = kw.pop('AR_coating',False)
        self.AR_type = kw.pop('AR_type','ideal')
        self.bond_gap = kw.pop('bond_gap', False)
        self.t_bondgap = kw.pop('t_bondgap',0.e-6)
        self.N_n2 = kw.pop('N_n2',1) #number of lambda/4 thicknesses of n2 (3 for narrow stopband, 1 for broad stopband)
        self.realistic_mirrors = kw.pop('realistic_mirrors',False)
        self.include_losses = kw.pop('include_losses',False)
        # print 'cavity with design wavelength: ',self.lambda_cav*1.e9,' nm'

        if self.dnu != 0.e9 and self.calculate_dnu:
            'warning: linewidth will be calculated, not externally set'

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

        if self.realistic_mirrors:
            self.N_n2_2 = kw.pop('N_n2_2',3)
            self.N_n2_1 = kw.pop('N_n2_1',1)
            self.M_2 = kw.pop('M_2',7)
            self.M_1 = kw.pop('M_1',10)

        if self.include_losses:
            print 'including losses'
            self.n_1_i = kw.pop('n_1_i',0.001)
            self.n_2_i = kw.pop('n_2_i',0.001)
            # self.n_1 = self.n_1 + 1j*self.n_1_i
            # self.n_2 = self.n_2 + 1j*self.n_1_i

        if self.AR_coating:
            if self.AR_type == 'ideal':
                # print 'using an ideal AR coating'
                self.n_AR = math.sqrt(n_diamond*n_air)
            elif self.AR_type == 'Al2O3':
                # print 'using an Al2O3 AR coating'
                self.n_AR = n_Al2O3
            else:
                self.n_AR = 1.
                'Specify valid AR_type (ideal or Al2O3)'

        if self.bond_gap:
            if self.AR_coating:
                'WARNING: setting bond gap and AR coating cannot be combined. only evaluating AR coating'
            else:
                'Setting bond gap to ',self.t_bondgap*1.e6,'um.'

        self.ns_in_cavity()

    def boundary_matrix(self,n1,n2,sigma=0):
        rho12 = (n1-n2)/(n1+n2)
        rho21 = (n2-n1)/(n1+n2)
        tau12 = 2*n1/(n1+n2)
        tau21 = 2*n2/(n1+n2)
        # rho = rho12
        # tau = 1+rho
        if not self.include_losses:
            rho = rho12 
            tau = tau12 #=1+rho
            matrix = (1./tau)*np.array(((1,rho),(rho,1)))
        else:
            if n1 == n_air and n2 == n_diamond:
                # print 'air-diamond interface'
                sigma=1.0e-9
                # print sigma
            if abs(n1) == n_diamond and abs(n2) == self.n_1:
                # print 'diamond-mirror interface'
                sigma=0.4e-9
                # print sigma
            if abs(n1) == self.n_1 and abs(n2) == n_air:
                # print 'air-mirror interface'
                sigma=0.2e-9
                # print sigma
            else:
                sigma=sigma
            rho12 = rho12*np.exp(-2.*(2*math.pi*sigma*n1/self.lambda_i)**2)
            rho21 = rho21*np.exp(-2.*(2*math.pi*sigma*n2/self.lambda_i)**2)
            tau12 = tau12*np.exp(-1./2*(2*math.pi*sigma*(n2-n1)/self.lambda_i)**2)
            tau21 = tau21*np.exp(-1./2*(2*math.pi*sigma*(n1-n2)/self.lambda_i)**2)
            # print rho,tau
            matrix = (1./tau12)*np.array(((1,-rho21),(rho12, (tau12*tau21)-(rho12*rho21) )) )
            # print matrix

        return matrix

    def propagation_matrix(self,n,t):
        k_c = 2*math.pi*n/self.lambda_i
        matrix = np.array(((np.exp(1j*k_c*t),0),(0,np.exp(-1j*k_c*t))))
        return matrix

    def electric_field_from_transfer_matrix(self,matrix,E):
        out = np.dot(matrix,E)
        return out  

    def ns_in_cavity(self):
        if self.realistic_mirrors:
            ns_mirror1 = np.append(np.array((self.n_1,self.n_2)*self.M_1),np.array((self.n_1)))
            ns_mirror2 = np.append(np.array((self.n_1,self.n_2)*self.M_2),np.array((self.n_1)))   
        else:
            ns_mirror = np.append(np.array((self.n_1,self.n_2)*self.M),np.array((self.n_1)))
            ns_mirror1=ns_mirror
            ns_mirror2=ns_mirror

        if self.cav_type == 'hybrid':
            if self.AR_coating:
                ns_cav = np.array((n_diamond, self.n_AR, n_air))
            elif self.bond_gap:
                ns_cav = np.array((n_air, n_diamond, n_air))
            else:
                ns_cav = np.array((n_diamond, n_air))

        elif self.cav_type == 'air':
            ns_cav = np.array([n_air])
        else:
            'specify valid cavity type (hyrbid or air)!'
            return

        self.ns = np.concatenate((ns_mirror2,ns_cav,ns_mirror1))  

        return self.ns

    def ts_in_cavity(self):
        if self.realistic_mirrors:
            ts_mirror_1 = np.append(np.array((self.lambda_cav/(4*self.n_1),self.N_n2_1*self.lambda_cav/(4*self.n_2))*self.M_1),np.array(self.lambda_cav/(4*self.n_1)))
            ts_mirror_2 = np.append(np.array((self.lambda_cav/(4*self.n_1),self.N_n2_2*self.lambda_cav/(4*self.n_2))*self.M_2),np.array(self.lambda_cav/(4*self.n_1)))
        else:
            ts_mirror = np.append(np.array((self.lambda_cav/(4*self.n_1),self.N_n2*self.lambda_cav/(4*self.n_2))*self.M),np.array(self.lambda_cav/(4*self.n_1)))
            ts_mirror_1=ts_mirror
            ts_mirror_2=ts_mirror

        if self.cav_type == 'hybrid':
            if self.AR_coating:
                ts_cav = np.array((self.t_d,self.lambda_cav/(4*self.n_AR),self.t_a))
            elif self.bond_gap:
                ts_cav = np.array((self.t_bondgap,self.t_d,self.t_a))
            else:
                ts_cav = np.array((self.t_d,self.t_a))

        elif self.cav_type == 'air':
            ts_cav = np.array([self.t_a])
        else:
            'specify valid cavity type (hyrbid or air)!'
            return

        self.ts = np.concatenate((ts_mirror_2,ts_cav,ts_mirror_1))
        self.tmirror1 = sum(ts_mirror_1)
        self.tmirror2 = sum(ts_mirror_2)

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

    def find_linewidth(self,nr_pts=21,sweep_range =4e9,plot_data=False):
        """
        find the linewidth by sweeping around lambda_i
        """
        ts = self.ts_in_cavity()
        self.freq_i = scipy.constants.c/self.lambda_i
        freq_is = np.linspace(-sweep_range,+sweep_range,nr_pts)
        lambda_was = self.lambda_i
        r_iis= np.zeros(len(freq_is))
        for i,freq_ii in enumerate(freq_is):
            self.lambda_i = scipy.constants.c/(freq_ii+self.freq_i)
            r = self.cavity_reflectivity()
            r_iis[i]=np.abs(r)**2
        

        self.lambda_i=lambda_was
        self.r=r_iis[np.argmin(r_iis)] #should be the same as initially, though seems not (always) to be
        halfmax = ((1-self.r)/2+self.r)
        g_offset=1
        g_x0 = self.freq_i
        g_gamma = abs(2*(freq_is[np.argmin(abs(r_iis - halfmax))]))
        g_A = -g_gamma

        fixed=[0]
        p0, fitfunc, fitfunc_str = common.fit_lorentz(g_offset, g_A, g_x0, g_gamma)
        fit_result = fit.fit1d((self.freq_i+freq_is),r_iis, None, p0=p0, 
            fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
        try:
            self.linewidth=fit_result['params_dict']['gamma']
            if plot_data:
                plot.plot_fit1d(fit_result,fit_xvals = self.freq_i+freq_is)
        except:
            'failed to fit linewidth. setting to guess'
            self.linewidth=g_gamma

        # if plot_data:
        #     fig,ax = plt.subplots()
        #     ax.set_xlabel('frequency THz')
        #     ax.set_ylabel('reflectivity')
        #     ax.plot((self.freq_i+freq_is)*1.e-12,r_iis,'o')
        #     ax.plot([self.freq_i*1.e-12-self.linewidth*1.e-12/2.,self.freq_i*1.e-12+self.linewidth*1.e-12/2.],[halfmax,halfmax],'c',linewidth=4)
        #     plt.show()
        #     plt.close()

        return self.linewidth

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
            ax.plot(lambda_is*1.e9,r_is, 'o')
            plt.show()
            plt.close()
        
        lambda_is = np.linspace(lambda_ii-self.res_wl_search_range/50.,lambda_ii+self.res_wl_search_range/50.,nr_pts/4)
        r_iis = np.zeros(len(lambda_is))
        for i,lambda_iii in enumerate(lambda_is):
            self.lambda_i = lambda_iii
            r = self.cavity_reflectivity()
            r_iis[i]=np.abs(r)**2

        self.lambda_i = lambda_is[np.argmin(r_iis)]
        self.r = r_iis[np.argmin(r_iis)]
        # print ((1-self.r)/2+self.r)
        # print np.argmin(abs(r_iis - ((1-self.r)/2+self.r)))
        # print lambda_is[np.argmin(abs(r_iis - ((1-self.r)/2+self.r)))]
        # print scipy.constants.c/lambda_is[np.argmin(abs(r_iis - ((1-self.r)/2+self.r)))]*1.e-12


        # self.linewidth = 2*(scipy.constants.c/self.lambda_i - scipy.constants.c/lambda_is[np.argmin(abs(r_iis - ((1-self.r)/2+self.r)))])

        if plot_r:
            print 'air gap',self.t_a*1.e6, ' um; r = ', self.r
            fig,ax = plt.subplots()
            ax.plot(lambda_is*1.e9,r_iis,'o')
            ax.set_xlabel('frequency (THz)')
            plt.show()
            plt.close()
        
        self.optical_length = self.t_a + self.t_d*n_diamond
        return self.lambda_i,self.r,r_iis#,t_ai,r_i,r_is#,

    def find_res_condition(self,nr_pts=161,plot_r=False):
        """
        find the resonance condition by varying t_a around t_a_g
        uses two consecutive optimisation steps
        could replace the second search for a lorentzian fit.
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
            ax.plot(t_as*1.e6,r_is,'o')
            plt.show()
            plt.close()
        
        t_as = np.linspace(t_ai-self.res_search_range/100.,t_ai+self.res_search_range/100.,1.5*nr_pts)
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
            ax.plot(t_as*1.e6,r_iis,'o')
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
        # self.w0 =  1.04e-6#*math.from my calculations of beam waist in cavity Riedel -> corresponds with both predictions from mathematica, and outcomes!! 
        # self.w0 =0.83e-6#0.705e-6# from 'FWHM beam waist' quoted in Riedels paper
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

    def calculate_E_max_at_NV(self):
        """
        calculate the maximum electric field density and its location, in the cavity region with refractive index n_z0
        """
        z0,E_z0,dz0 = self.calculate_E_max_in_n()
        self.z0 = z0
        self.E_z0 = E_z0
        self.dz0 = dz0

    def calculate_E_max_in_n(self, **kw):
        n_0 = kw.pop('n_0',self.n_z0)
        dz = np.abs(self.zs[0]-self.zs[1])
        arg_zs_in_cavity = np.where((self.ns_vs_z<n_0+0.001)&(self.ns_vs_z>n_0-0.001))
        zs_in_cavity = self.zs[arg_zs_in_cavity]
        arg_shallow_zs_in_cavity = arg_zs_in_cavity[0][-int(self.lambda_i/2/n_0/dz):]

        z0_i = arg_shallow_zs_in_cavity[0] + np.argmax(self.Etot_vs_z[arg_shallow_zs_in_cavity])
        z0 = self.zs[z0_i]
        E_z0 = self.Etot_vs_z[z0_i]
        dz0 = np.abs(zs_in_cavity[-1]-z0)
        return z0,E_z0,dz0

    def calculate_energy_dist_length(self):
        """
        calculate energy distribution length int((n(z)^2)*(E(z))^2 dz)/(n(z0))^2*E(z0)^2
        """
        dz = np.abs(self.zs[0]-self.zs[1])
        self.effective_length = np.real(np.sum((self.ns_vs_z**2)*(self.Etot_vs_z**2))*dz/((self.n_z0**2)*(self.E_z0**2)))
        return self.effective_length

    def calculate_energy_dist_length_for_args(self,args):
        """
        calculate energy distribution length int((n(z)^2)*(E(z))^2 dz)/(n(z0))^2*E(z0)^2
        """
        ns_vs_z_in_m = self.ns_vs_z[args]
        Etot_vs_z_in_m = self.Etot_vs_z[args]
        # print ns_vs_z_in_m
        dz = np.abs(self.zs[0]-self.zs[1])
        effective_length_in_m = np.real(np.sum((ns_vs_z_in_m**2)*(Etot_vs_z_in_m**2))*dz/((self.n_z0**2)*(self.Evac_max**2)))
        return effective_length_in_m

    def calculate_energy_dist_length_in_m(self,n_m):
        """
        calculate the electric field energy (J/m^2) in cavity material with refractive index n_m
        """
        arg_zs_in_m = np.where((self.ns_vs_z<n_m+0.01)&(self.ns_vs_z>n_m-0.01))
        effective_length_in_m = self.calculate_energy_dist_length_for_args(arg_zs_in_m)
        return effective_length_in_m #in J/m^2 

    def calculate_mode_volume(self):
        self.mode_volume = math.pi*self.w0**2/2.*self.effective_length #note that since we are using the effective length, we here have to use pi*wo^2/2. this is int(|E(x,y,z0)|^2dxdy) = int(e^(-2r^2/w0^2)rdr = pi*w0^2/2)
        #print 'mode vol',self.mode_volume
        return self.mode_volume
        
    def calculate_max_Evac(self):
        """
        Uses the mode volume to determine the maximum vacuum electric field (hbar omega/(2epsilon*V)) (epsilon = n^2*epsilon_0)
        """
        self.Evac_max =np.sqrt(sc.constants.hbar*2*math.pi*(sc.constants.c/self.lambda_i)/(2*sc.constants.epsilon_0*self.n_z0**2*self.mode_volume))

    def calculate_energy_in_cavity(self):
        dz = np.abs(self.zs[0]-self.zs[1])
        Etot_L = 1./2*sc.constants.epsilon_0*np.sum((self.ns_vs_z**2)*(self.Etot_vs_z**2))*dz 
        Etot_V = Etot_L*math.pi*self.w0**2/2
        return Etot_L,Etot_V

    def calculate_energy_in_args(self,args):
        ns_vs_z_in_m = self.ns_vs_z[args]
        # arg=args[0]
        # print arg
        # print args
        # if len(arg)>1:
        #     plt.figure()
        #     plt.plot(self.zs,self.ns_vs_z)

        #     plt.plot(self.zs, np.concatenate((np.zeros(arg[0]),np.ones(len(arg)),np.zeros(len(self.ns_vs_z)-arg[-1]-1))))
        #     plt.show()
        #     plt.close()

        Etot_vs_z_in_m = self.Etot_vs_z[args]
        dz = np.abs(self.zs[0]-self.zs[1])
        # print ns_vs_z_in_m
        # print 'length', ns_vs_z_in_m[1] ,np.sum(np.ones(len(ns_vs_z_in_m)))*dz *1.e6

        Etot = 1./2*sc.constants.epsilon_0*np.sum((ns_vs_z_in_m**2)*(Etot_vs_z_in_m**2))*dz 
        return Etot #in J/m^2 

    def calculate_energy_in_m(self,n_m):
        """
        calculate the electric field energy (J/m^2) in cavity material with refractive index n_m
        """
        arg_zs_in_m = np.where((self.ns_vs_z<n_m+0.01)&(self.ns_vs_z>n_m-0.01))
        # ns_vs_z_in_m = self.ns_vs_z[arg_zs_in_m]
        # Etot_vs_z_in_m = self.Etot_vs_z[arg_zs_in_m]
        # dz = np.abs(self.zs[0]-self.zs[1])
        # Etot = 1./2*sc.constants.epsilon_0*np.sum((ns_vs_z_in_m**2)*(Etot_vs_z_in_m**2))*dz 
        Etot = self.calculate_energy_in_args(arg_zs_in_m)
        return Etot #in J/m^2 

    def calculate_pen_depth(self):
        # Eair = self.Etot_air
        # Ediamond = self.Etot_diamond
        # Etot = self.Etot
        # if self.AR_coating:
        #     EAR = self.calculate_energy_in_m(n_m=self.n_AR)F
        #     Emirrors = Etot - Ediamond - Eair - EAR
        # else:
        #     Emirrors = Etot - Ediamond - Eair 


        argsmirror1 = np.where(abs(self.zs < self.tmirror1))
        Emirror1 = self.calculate_energy_in_args(argsmirror1)
        argsmirror2 = np.where(abs(self.zs>(self.zs[-1]-self.tmirror2)))
        Emirror2 = self.calculate_energy_in_args(argsmirror2)

        Emax_air = self.calculate_E_max_in_n(n_0=1)[1]
        # t_pen_mirror = Emirror/(Eair+Ediamond)*(self.t_a+n_diamond*self.t_d)
        if self.cav_type == 'hybrid':
            Emax_diamond = self.calculate_E_max_in_n(n_0=n_diamond)[1]

            # self.t_pen_mirror1 = Emirror1/(Ediamond)*(n_diamond*self.t_d)
            # self.t_pen_mirror2 = Emirror2/(Ediamond)*(n_diamond*self.t_d)
            #     # print '1',self.t_pen_mirror1*1.e6

            #compare both to the max field in diamond such that get close to effective energy distribution length
            self.t_pen_mirror1 = Emirror1/(1/4.*scipy.constants.epsilon_0*n_diamond*Emax_diamond**2)
            self.t_pen_mirror2 = Emirror2/(1/4.*scipy.constants.epsilon_0*n_diamond*Emax_diamond**2)

            # print '2',self.t_pen_mirror1*1.e6
            # if self.AR_coating:
            #     # self.t_pen_mirror1 = Emirror1/(Eair+EAR+Ediamond)*(self.t_a+n_diamond*self.t_d+self.lambda_cav/4)
            #     # self.t_pen_mirror2 = Emirror2/(Eair+EAR+Ediamond)*(self.t_a+n_diamond*self.t_d+self.lambda_cav/4)

            # else:
            #     self.t_pen_mirror1 = Emirror1/(Eair+Ediamond)*(self.t_a+n_diamond*self.t_d)
            #     self.t_pen_mirror2 = Emirror2/(Eair+Ediamond)*(self.t_a+n_diamond*self.t_d)    

        elif self.cav_type == 'air':
            self.t_pen_mirror1 =  Emirror1/(1/4.*scipy.constants.epsilon_0*Emax_air**2)
            self.t_pen_mirror2 = Emirror2/(1/4.*scipy.constants.epsilon_0*Emax_air**2)
        else: 
            print 'specify valid cavity type!'

        # print Emirror1, Emirror2
        # print 'Eair', Eair/self.t_a, 'Ediamond', Ediamond/(n_diamond*self.t_d) ,
        # if self.AR_coating:
        #     print 'EAR', EAR/(self.n_AR*self.lambda_cav/(4*self.n_AR)) 
        # print 'Eair/diamond', Eair/self.t_a/(Ediamond/(n_diamond*self.t_d))
        # print Eair/Etot,Ediamond/Etot,Emirrors/Etot

        # print n_diamond*self.t_d/(self.t_a+n_diamond*self.t_d),self.t_a/(self.t_a+n_diamond*self.t_d)
        # print 't_pen_mirror',t_pen_mirror
        # print 't_pen_mirror1',t_pen_mirror1
        # print 't_pen_mirror2',t_pen_mirror2
        return self.t_pen_mirror1,self.t_pen_mirror2

    def calculate_dnu_from_LM1_LM2(self,LM1=100e-6,LM2=2200e-6):
        """
        calculate the linewidth based on losses at the mirrors only (no diamond losses), using the unequal distribution of E-field in the cavity
        """
        self.calculate_pen_depth()
        if self.cav_type == 'hybrid':
            Emax_air = self.calculate_E_max_in_n(n_0=1)[1]
            Emax_diamond = self.calculate_E_max_in_n(n_0=n_diamond)[1]
            f = Emax_air**2/(n_diamond*Emax_diamond**2) #ratio of time-averaged poynting vectors!
            dnu = sc.constants.c/(4*math.pi)*(LM1*f+LM2)/(2*n_diamond*self.effective_length)
        else:
            dnu = sc.constants.c/(4*math.pi)*(LM1+LM2)/(2*self.effective_length)
            # if self.AR_coating:
            #     dnus =  sc.constants.c/(4*math.pi)*(LM2 + LM1)/(self.effective_length)
            #     dnu = sc.constants.c/(4*math.pi)*(LM2 + LM1)/(n_diamond*self.t_d + self.t_a+self.t_pen_mirror1+self.t_pen_mirror2+self.lambda_cav/4)
            # else:

                
            #     dside = n_diamond*self.t_d+self.t_pen_mirror2 #~optical cycles on the diamond side
            #     aside = self.t_a + self.t_pen_mirror1 #~optical cycles on the air side
            #     print 'use length',(f/(1-f)*aside+dside)
            #     dnu = sc.constants.c/(4*math.pi)*(LM1*f+LM2)/(f/(1-f)*aside+dside)
                
        return dnu

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

    def analyse_cavity(self,nr_points=50001,plot_Evac=False,save_plot=True):
        self.ts_in_cavity()
        self.electric_field_distribution(nr_points=nr_points)

        self.calculate_E_max_at_NV()
        self.calculate_energy_dist_length()
        self.calculate_mode_volume()#t_a+t_d*2.4)
        self.calculate_max_Evac()
        self.Etot_vs_z = self.Etot_vs_z*self.Evac_max/self.E_z0


        self.Etot,self.Etot_3d = self.calculate_energy_in_cavity()
        self.Etot_diamond = self.calculate_energy_in_m(n_diamond)
        self.Etot_air = self.calculate_energy_in_m(n_air)

        if self.calculate_dnu:
            self.find_linewidth()
            self.dnu = self.linewidth

        self.calculate_Purcell()
        self.calculate_into_ZPL()

        if plot_Evac:
            # print 'energy distribution length = ',self.effective_length*1.e6, 'um'
            # print 'beam waist', self.w0*1.e6, 'um'
            # print 'mode volume', self.mode_volume*(1.e6)**3, 'um^3 = ', self.mode_volume/(self.lambda_i**3), 'lambda^3'
            # print 'max Evac',self.Evac_max/1000.,'kV/m' 
            today = tb.get_timestamp_from_now()[:8]
            title_string = today+'_td_%.2f_ta_%.2f_lambdai_%.1f_R_%.1f'%(self.t_d*1.e6,self.t_a*1.e6,self.lambda_i*1.e9,self.R*1.e6)

            fig,ax = plt.subplots(figsize=(12,6))
            ax.plot(self.zs*1.e6,self.ns_vs_z)
            ax.set_ylabel('n')
            ax.set_xlabel('z (um)')
            ax.set_title(title_string)
            # ax.set_ylim((1.,2.6))
            ax2=ax.twinx()
            # ax2.set_xlim(4,6)#((1.8,2.3))
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
    ax2.set_xlim((min(dnus*1.e-9),10))

    plt.show()


    fig,ax = plt.subplots()
    ax.plot(dnus*1.e-9,taus,'magenta',label='Fp')
    ax.set_ylabel('tau (ns)',fontsize=14)
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
    ax2.set_xlim((min(dnus*1.e-9),10))

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

def calculate_avgintoZPL_vs_vibrations(t_d=4.e-6,t_a_g=1.2e-6,Rs=[18.e-6],vib_dLs=[0.4e-9], dLmax=1.0e-9, Ltots = [1000e-6],LM1s = [242e-6],lambda_ZPL=637.e-9,
        nr_pts = 51, show_plots=False, plot_r=False,AR_coating=False,realistic_mirrors=False, method='numeric',save_data=False,tag=''):
    """
    calculate the average emission into the ZPL, and outcoupling throught the top mirror, under the influence of vibrations and varying division of losses over the mirrors.
    For the numeric case, always assume that the purcell factor on-resonance does not change for the small fluctuations in air gap due oto vibrations. 
    This is a valid approximation, as can be verified using the function 'calculate_Fp_w_vibrations'.
    """

    s = Cavity(t_d, t_a_g, Rs[0],lambda_i=lambda_ZPL,realistic_mirrors=realistic_mirrors,AR_coating=AR_coating,calculate_dnu=True)
    s.find_res_condition(plot_r=plot_r)
    t_a0 = s.t_a
    s.calculate_w0()
    s.analyse_cavity(plot_Evac=False,nr_points=20001)    
    linewidth_for_seach_range = s.linewidth #we use this value ONLY to determine the res_wl_search_range of the cavity objects below. It is NOT the linewidth of our cavity, that is specified by the different losses.
    res_wl_search_range = 10*linewidth_for_seach_range*scipy.constants.c/((scipy.constants.c/lambda_ZPL)**2)#transform change in frequency to change in wavelength. multiply by 20 to do wide range search initially.
    dLmax = max(2*res_wl_search_range,1.8*max(vib_dLs))

    #In the high-finesse limit, the resonance condition does not end up at 637nm. The cavity object below is used to calculate the resonant wavelength.
    sr = Cavity(t_d, t_a0, Rs[0],res_wl_search_range=res_wl_search_range,AR_coating=AR_coating,realistic_mirrors=realistic_mirrors)
    sr.find_res_wavelength(plot_r=False,nr_pts=101)
    lambda_ref = sr.lambda_i #use this for simple_analytics and numeric methods as reference lambda in the spectral overlap.

    dLs = np.linspace(-dLmax,dLmax,nr_pts)

    lambda_ress = np.zeros(len(dLs))

    if method == 'simple_analytic': #use dnu=dL/L*nu
        dnus = -dLs*(scipy.constants.c/lambda_ZPL)/(2*n_diamond*s.effective_length)
        lambda_ress = scipy.constants.c/((scipy.constants.c/lambda_ZPL)+dnus)

    elif method == 'analytic':
        Nmin =int(2*(t_a0+t_d*n_diamond)*450e12/scipy.constants.c)
        N= np.arange(Nmin,Nmin+10)
        r_diamond=(1-n_diamond)/(n_diamond+1)

        nu0_analytics = scipy.constants.c/(2*math.pi*(t_a0+n_diamond*t_d))*(math.pi*(N) - (-1)**(N)*np.arcsin( -r_diamond*np.sin((N)*math.pi*(t_a0 - n_diamond*t_d)/(t_a0 +n_diamond*t_d) )  )  )
        nu0_res = nu0_analytics[np.argmin(abs(nu0_analytics-scipy.constants.c/lambda_ZPL))]
        N0 = N[np.argmin(abs(nu0_analytics-scipy.constants.c/lambda_ZPL))]
        lambda_ref=scipy.constants.c/nu0_res

        t_as = t_a0+dLs
        nu_analytics = scipy.constants.c/(2*math.pi*(t_as+n_diamond*t_d))*(math.pi*(N0) - (-1)**(N0)*np.arcsin( -r_diamond*np.sin((N0)*math.pi*(t_as - n_diamond*t_d)/(t_as +n_diamond*t_d) )  )  )
        lambda_ress = scipy.constants.c/nu_analytics

    elif method == 'numeric':

        for i,dL in enumerate(dLs):
            si = Cavity(t_d, t_a0+dL, Rs[0],res_wl_search_range=res_wl_search_range,AR_coating=AR_coating,realistic_mirrors=realistic_mirrors)
            si.find_res_wavelength(plot_r=False,nr_pts=101)
            lambda_ress[i] = si.lambda_i
        

    else:
        'WARNING: define a valid method!'
   
    t0 = time.time()


    #Incorporating the option to sweep: Ltots, LM1s, and vib_dLs
    dnu_from_losses = np.zeros((len(Ltots),len(LM1s)))
    intoZPL0s = np.zeros((len(Rs),len(Ltots),len(LM1s)))
    intoZPLs = np.zeros((len(Rs),len(Ltots),len(LM1s),len(dLs)))
    avg_p_ZPLs = np.zeros((len(Rs),len(Ltots),len(LM1s),len(vib_dLs) ))
    out_through_M2s = np.zeros((len(Rs),len(Ltots),len(LM1s),len(vib_dLs) ))
    w0s = np.zeros(len(Rs))
    mode_volumes=np.zeros(len(Rs))

    for i, R in enumerate(Rs):
        s.R = R
        s.calculate_w0()
        s.calculate_mode_volume()
        w0s[i] = s.w0
        mode_volumes[i]=s.mode_volume
        for j, Ltot in enumerate(Ltots):
            for k,LM1 in enumerate(LM1s):
                dnu_from_losses[j,k] = s.calculate_dnu_from_LM1_LM2(LM1 = LM1,LM2 = Ltot-LM1)#make sure that max(LM1s),min(Ltots)
                s.dnu = dnu_from_losses[j,k] #need this to get the right into ZPL calculation
                s.calculate_Purcell()
                intoZPL0s[i,j,k] = s.calculate_into_ZPL()
                intoZPLs[i,j,k,:] = intoZPL0s[i,j,k] * spectral_overlap(lambda_ref,lambda_ress,dnu_from_losses[j,k] )
                for l,vib_dL in enumerate(vib_dLs):    #here we are really taking the size of the vibrations into account
                    p_cav_length = sim_vib.gaussian(dLs,0,vib_dL)
                    avg_p_ZPLs[i,j,k,l] = sim_vib.avg_p_ZPL_to_zero(intoZPLs[i,j,k,:],p_cav_length)
                    # out_through_M2s[i,j,k,l] = avg_p_ZPLs[i,j,k,l]*(Ltot-LM1)/Ltot

            # if j%40==0:
            #     t1 = time.time()-t0
            #     print 'Ltot %.0f at'%(Ltot*1.e6), t1
            #     print '%d out of %d done'%(j+1, len(Ltots))
            #     print 'expected time remaining: %d seconds'%(t1/(j+1)*len(Ltots)-t1)


    LM1smesh , Ltotsmesh = np.meshgrid(LM1s,Ltots)
    dimarray=np.array([1,len(Ltots),len(LM1s),1])
    frac_out_through_M2s = ((Ltotsmesh-LM1smesh)/(Ltotsmesh)).reshape(dimarray)
    out_through_M2s = avg_p_ZPLs*frac_out_through_M2s #use broadcasting

    if len(Ltots)>1:       
        i_max_out_through_M2s = np.argmax(out_through_M2s,axis=1)
        max_out_through_M2s = np.amax(out_through_M2s,axis=1)
        Ltot_max_out_through_M2s = Ltots[i_max_out_through_M2s]
    else:
        i_max_out_through_M2s =0
        max_out_through_M2s =0
        Ltot_max_out_through_M2s=0

    if show_plots:
        fig,ax = plt.subplots()
        ax.plot(dLs, p_cav_length/np.sum(p_cav_length),'o')
        ax.set_xlabel('dL w.r.t. resonance (nm)')
        ax.set_ylabel('probability distribution')
        ax2=ax.twinx()
        ax2.plot(dLs*1.e9,into_ZPLs,'o')
        ax2.set_ylabel('into ZPL')
        plt.show()
        plt.close()

        fig,ax = plt.subplots()
        ax.plot(dLs*1.e9,lambda_ress*1.e9,'o')
        ax.set_xlabel('dL w.r.t. resonance (nm)')
        ax.set_ylabel('resonant wavelength')
        plt.show()
        plt.close()


    if save_data:
        filename = '%s_avgZPLvsvibs_method_%s_R_%.1f_td_%.2f_ta0_%.2f_%s.hdf5'%(tb.get_timestamp_from_now()[8:14],method,R*1.e6,t_d*1.e6,t_a0*1.e6,tag)

        if not os.path.exists(os.path.join(data_folder, filename)):
            print 'creating analysis file'
            mode = 'w'    

        f = h5py.File(os.path.join(data_folder, filename), 'w')
        g = f.require_group('results')

        f['/results/avg_ZPLs'] = avg_p_ZPLs
        # f['/results/out_through_M2s'] = out_through_M2s #do not save, can reprocudec easily
        f['/results/dnu']=dnu_from_losses
        f['/results/lambda_ress'] = lambda_ress
        f['/results/intoZPLs'] = intoZPLs
        f['/results/intoZPL0s'] = intoZPL0s
        f['/results/w0s'] = w0s
        f['/results/mode_volumes'] = mode_volumes

        data_dict = {}
        data_dict['vibs_dL']=vib_dLs
        # data_dict['max_out_through_M2s']=max_out_through_M2s #ot saving, since too large, and can be easily reconstructed! 
        # data_dict['Ltot_max_out_through_M2s']=Ltot_max_out_through_M2s #ot saving, since too large, and can be easily reconstructed!
        data_dict['Rs']=Rs
        data_dict['t_d']=t_d
        data_dict['t_a_g']=t_a_g
        data_dict['t_a']=t_a0
        data_dict['method']=method
        data_dict['Ltots'] = Ltots
        data_dict['LM1s'] = LM1s
        data_dict['dLmax'] = dLmax
        data_dict['nr_pts'] = nr_pts
        data_dict['lambda_ref'] = lambda_ref
        data_dict['res_wl_search_range'] = res_wl_search_range

        data_dict['realistic_mirrors'] = realistic_mirrors
        data_dict['AR_coating'] = AR_coating


        for k in data_dict:
            g.attrs[k] = data_dict[k]
                
        f.close()

    return intoZPLs,lambda_ress,avg_p_ZPLs,max_out_through_M2s,Ltot_max_out_through_M2s


def calculate_avg_into_ZPL_vs_tds(t_ds=np.linspace(1.e-6,15.e-6,15), t_a_gs = np.ones(15)*1.2e-6,R=18.e-6,vib_dL=0.4e-9, dLmax=1.0e-9,LM1=146e-6,Ltot=1000e-6,lambda_ZPL=637.e-9,
            nr_pts=51, AR_coating = False, realistic_mirrors =False, method = 'numeric',save_data=False,tag=''):
    
    t0 = time.time()
    
    lambda_ress = np.zeros((len(t_ds),nr_pts))
    avg_p_ZPLs = np.zeros(len(t_ds))

    for i,t_d in enumerate(t_ds):
        t_a_g = t_a_gs[i]
        a,lambda_ress[i,:],avg_p_ZPLs[i],b,c = calculate_avgintoZPL_vs_vibrations(t_d=t_d,t_a_g=t_a_g,Rs=[R],vib_dLs=[vib_dL], dLmax=dLmax, Ltots = [Ltot],LM1s = [LM1],lambda_ZPL=lambda_ZPL,
            nr_pts = nr_pts, show_plots=False, plot_r=False,AR_coating=AR_coating,realistic_mirrors=realistic_mirrors, method=method,save_data=False,tag=tag)

        t1 = time.time()-t0
        if i%10==0:
            print '%d out of %d done'%(i+1, len(t_ds))
            print 'expected time remaining: %d seconds'%(t1/(i+1)*len(t_ds)-t1)

    if save_data:
        filename = '%s_avgZPLvstds_method_%s_R_%.1f_vibdL_%.2f_tag_%.2f_AR_%s_RealM_%s_%s.hdf5'%(tb.get_timestamp_from_now()[8:14],method,R*1.e6,vib_dL*1.e6,t_a_g*1.e6,str(AR_coating),str(realistic_mirrors),tag)

        if not os.path.exists(os.path.join(data_folder, filename)):
            print 'creating analysis file'
            mode = 'w'    

        f = h5py.File(os.path.join(data_folder, filename), 'w')
        g = f.require_group('results')

        f['/results/avg_ZPLs'] = avg_p_ZPLs
        f['/results/lambda_ress'] = lambda_ress

        data_dict = {}
        data_dict['vibs_dL']=vib_dL
        data_dict['Rs']=R
        data_dict['t_ds']=t_ds
        data_dict['t_a_gs']=t_a_gs
        data_dict['method']=method
        data_dict['Ltots'] = Ltot
        data_dict['LM1s'] = LM1
        data_dict['dLmax'] = dLmax
        data_dict['nr_pts'] = nr_pts
        data_dict['tag']=tag

        data_dict['realistic_mirrors'] = realistic_mirrors
        data_dict['AR_coating'] = AR_coating

        for k in data_dict:
            g.attrs[k] = data_dict[k]
                
        f.close()

    fig,ax = plt.subplots()
    ax.plot(t_ds*1.e6,avg_p_ZPLs)
    ax.set_xlabel('diamond thickness (um)')
    ax.set_ylabel('avg emission into ZPL')
    plt.show()
    plt.close()

    return avg_p_ZPLs

def calculate_Fp_w_vibrations(t_d=4.e-6,t_a_g=1.2e-6,R=18.e-6,vib_dL=0.4e-9,losses=1000e-6,dnu=3.5e9,
            lambda_ZPL=637.e-9,dLmax=1.e-9,nr_pts=51,show_plots=False,
            fast_approx=False,AR_coating=False,plot_r=False,LM1=242e-6): #dL in nm
    """
    Calculates the average emission into the ZPL, when taking vibrations into account.
    Assumes a gaussian distribution of cavity length with FWHM of vib_dL as a result of vibrations
    When using fast_approx this function does the same as calculate_avgintoZPL_w_vibrations with method=numeric.
    """
    s = Cavity(t_d, t_a_g, R,dnu=dnu,lambda_i=lambda_ZPL,realistic_mirrors=True,AR_coating=AR_coating, res_search_range=0.01e-6)
    s.find_res_condition(plot_r=plot_r,nr_pts=2001)
    t_a0 = s.t_a
    s.calculate_w0()
    s.analyse_cavity(plot_Evac=False,nr_points=20001)   
    dnu_from_losses = s.calculate_dnu_from_LM1_LM2(LM1 = LM1,LM2 = losses-LM1)#calculate the linewidth based on the losses!
    s.dnu = dnu_from_losses #need this to get the right into ZPL calculation
    s.calculate_Purcell()
    intoZPL0 = s.calculate_into_ZPL()

    #In the high-finesse limit, the resonance condition does not end up at 637nm. The cavity object below is used to calculate the resonant wavelength.
    sr = Cavity(t_d, t_a0, R,dnu=dnu, res_wl_search_range=1e-10,AR_coating=AR_coating,realistic_mirrors=True)
    sr.find_res_wavelength(plot_r=False,nr_pts=101)

    dLs = np.linspace(-dLmax,dLmax,nr_pts)
    p_cav_length = sim_vib.gaussian(t_a0+dLs,t_a0,vib_dL)

    into_ZPLs = np.zeros(len(dLs))
    lambda_ress = np.zeros(len(dLs))
    Purcells=np.zeros(len(dLs))
    dnus_si=np.zeros(len(dLs))

    t0 = time.time()

    for i,dL in enumerate(dLs):
        si = Cavity(t_d, t_a0+dL, R,dnu=dnu, res_wl_search_range=1e-10,AR_coating=AR_coating,realistic_mirrors=True)
        si.find_res_wavelength(plot_r=False,nr_pts=101)

        lambda_ress[i] = si.lambda_i

        if not fast_approx: #in the fast approximation, do not re-calculate the emission intoZPL on resonance; rather assume the change is negligable for each.
            si.calculate_w0()
            si.analyse_cavity(plot_Evac=False)
            dnu_from_losses_si = si.calculate_dnu_from_LM1_LM2(LM1 = LM1,LM2 = losses-LM1)#(LM1 = losses/2,LM2 = losses/2)#fix base losses for now!!242e-6
            si.dnu = dnu_from_losses_si
            dnus_si[i]=dnu_from_losses_si
            Purcells[i] = si.calculate_Purcell()
            intoZPLonR = si.calculate_into_ZPL()
            into_ZPL = intoZPLonR  *spectral_overlap(sr.lambda_i,si.lambda_i,dnu_from_losses_si)
            into_ZPLs[i] = into_ZPL

        t1 = time.time()-t0
        if i%50==0:
            print 'dL %.2f nm at'%(dL*1.e9), t1
            print '%d out of %d done'%(i+1, len(dLs))
            print 'expected time remaining: %d seconds'%(t1/(i+1)*len(dLs)-t1)

    if fast_approx:
        dnus_si = np.ones(len(into_ZPLs))*dnu_from_losses
        into_ZPLs = intoZPL0  *spectral_overlap(sr.lambda_i,lambda_ress,dnu_from_losses)

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


    return into_ZPLs,Purcells,lambda_ress,avg_p_ZPL,dnus_si

# def calc_avg_ZPL_vs_vibrations(t_d=4.e-6,t_a_g=1.2e-6,R=18.e-6,vib_dLs=np.linspace(0.05,1.,20)*1.e-9,dnu=3.5e9,
#         lambda_ZPL=637.e-9,dLmax=None,nr_pts=31,save_data=False,analytic_version=False,AR_coating=False):
#     """
#     calculate average emission into ZPL into cavity vs FWHM of vibrations
#     """
#     if dLmax==None:
#         dLmaxs = vib_dLs*1.8
#     else:
#         dLmaxs = [dLmax]*len(vib_dLs)

#     avg_ZPLs = np.zeros(len(vib_dLs))
#     t0 = time.time()

#     for i,vib_dL in enumerate(vib_dLs):
#         dLmax = dLmaxs[i]
#         # print vib_dL
#         if analytic_version:
#             a,b,avg_ZPL = calculate_Fp_w_fullyanalytic_vibrations(t_d=t_d,t_a_g=t_a_g,R=R,vib_dL=vib_dL,dnu=dnu,
#                 lambda_ZPL=lambda_ZPL,dLmax=dLmax,nr_pts=nr_pts,show_plots=False,AR_coating=AR_coating)
#         else:
#             a,b,c,avg_ZPL = calculate_Fp_w_vibrations(t_d=t_d,t_a_g=t_a_g,R=R,vib_dL=vib_dL,dnu=dnu,
#                 lambda_ZPL=lambda_ZPL,dLmax=dLmax,nr_pts=nr_pts,fast_approx=False,show_plots=True,AR_coating=AR_coating)
#         # print avg_ZPL
#         avg_ZPLs[i]=avg_ZPL
#         t1 = time.time()-t0
#         print '%d out of %d done'%(i+1, len(vib_dLs))
#         print 'expected time remaining: %d seconds'%(t1/(i+1)*len(vib_dLs)-t1)


#     if save_data:
#         filename = '%s_avgZPLvsvibs_R_%.1f_td_%.2f_tag_%.2f_dnu_%.2f'%(tb.get_timestamp_from_now()[8:14],R*1.e6,t_d*1.e6,t_a_g*1.e6,dnu*1.e-9)
#         data_dict = {}
#         data_dict['vibs_dL']=vib_dLs
#         data_dict['avg_ZPLs']=avg_ZPLs
#         data_dict['R']=R
#         data_dict['t_d']=t_d
#         data_dict['t_a_g']=t_a_g
#         data_dict['dnu']=dnu
#         data_dict['dLmax'] = dLmax
#         data_dict['nr_pts'] = nr_pts
#         save_to_json_file(data_folder,filename,data_dict)


#     fig,ax = plt.subplots()
#     ax.plot(vib_dLs,avg_ZPLs)
#     ax.set_xlabel('FWHM vibrations (nm)')
#     ax.set_ylabel('avg emission into ZPL')
#     plt.show()
#     plt.close()

#     return avg_ZPLs

# def calc_avg_ZPL_vs_finesse(t_d=4.e-6,t_a_g=1.2e-6,R=18.e-6,dnus=np.linspace(1,10,10)*1.e9,vib_dL=0.4,lambda_ZPL=637.e-9,
#         dLmax=1.e-9,nr_pts=51,save_data=False,AR_coating=False):
#     """
#     calculate average emission into ZPL into cavity vs losses
#     """
#     avg_ZPLs = np.zeros(len(dnus))
#     t0 = time.time()
#     for i,dnu in enumerate(dnus):
#         # print vib_dL
#         if analytic_version:
#             a,b,avg_ZPL = calculate_Fp_w_fullyanalytic_vibrations(t_d=t_d,t_a_g=t_a_g,R=R,vib_dL=vib_dL,dnu=dnu,
#                     lambda_ZPL=lambda_ZPL,dLmax=dLmax,nr_pts=nr_pts,fast_approx=True,AR_coating=AR_coating)
#         else:
#             a,b,c,avg_ZPL = calculate_Fp_w_vibrations(t_d=t_d,t_a_g=t_a_g,R=R,vib_dL=vib_dL,dnu=dnu,
#                     lambda_ZPL=lambda_ZPL,dLmax=dLmax,nr_pts=nr_pts,fast_approx=True,AR_coating=AR_coating)
#         # print avg_ZPL
#         t1 = time.time()-t0
#         # print 'at', t1
#         print '%d out of %d done'%(i+1, len(dnus))
#         print 'expected time remaining: %d seconds'%(t1/(i+1)*len(dnus)-t1)

#         avg_ZPLs[i]=avg_ZPL

#     if save_data:
#         filename = '%s_avgZPLvsF_R_%.1f_td_%.2f_tag_%.2f_vibZPL_%.2f'%(tb.get_timestamp_from_now()[8:14],R*1.e6,t_d*1.e6,t_a_g*1.e6,vib_dL*1.e9)
#         data_dict = {}
#         data_dict['vib_dL']=vib_dL
#         data_dict['avg_ZPLs']=avg_ZPLs
#         data_dict['R']=R
#         data_dict['t_d']=t_d
#         data_dict['t_a_g']=t_a_g
#         data_dict['dnus']=dnus
#         data_dict['dLmax'] = dLmax
#         data_dict['nr_pts'] = nr_pts
#         save_to_json_file(data_folder,filename,data_dict)

#     fig,ax = plt.subplots()
#     ax.plot(dnus*1.e-9,avg_ZPLs)
#     ax.set_xlabel('cavity linewidth (GHz)')
#     ax.set_ylabel('avg emission into ZPL')
#     plt.show()
#     plt.close()



#     return avg_ZPLs



# def calc_avg_ZPL_vs_t_d_s(t_d_s=np.linspace(1.e-6,15.e-6,15),t_a_g=1.2e-6,R=18.e-6,vib_dL=0.3e-9,losses=1000e-6,
#         lambda_ZPL=637.e-9,dLmax=0.54e-9,nr_pts=31,save_data=False,mode='',analytic_version=False,AR_coating=False):
#     """
#     calculate average emission into ZPL into cavity vs thickness of diamond
#     """
#     avg_ZPLs = np.zeros(len(t_d_s))
#     t0 = time.time()

#     for i,t_d in enumerate(t_d_s):
#         if i==0: 
#             plot_r=True
#         else:
#             plot_r=False
#         if AR_coating:
#             dnu = (scipy.constants.c/(2*(t_d*n_diamond+t_a_g+lambda_ZPL/4)))/(2*math.pi/losses)#finesse to dnu calculation
#         else:
#             dnu = (scipy.constants.c/(2*(t_d*n_diamond+t_a_g)))/(2*math.pi/losses)#finesse to dnu calculation

#         if dnu<1:
#             nr_points=601
#         elif dnu<2:
#             nr_points=401
#         elif dnu<4:
#             nr_points=201
#         else:
#             nr_points=nr_pts
#         if analytic_version:
#             a,b,avg_ZPL = calculate_Fp_w_fullyanalytic_vibrations(t_d=t_d,t_a_g=t_a_g,R=R,vib_dL=vib_dL,losses=losses,
#                 lambda_ZPL=lambda_ZPL,dLmax=dLmax,nr_pts=nr_pts,show_plots=False,AR_coating=AR_coating)
#         else:
#             a,b,c,avg_ZPL = calculate_Fp_w_vibrations(t_d=t_d,t_a_g=t_a_g,R=R,vib_dL=vib_dL,dnu=dnu,losses=losses,
#                 lambda_ZPL=lambda_ZPL,dLmax=dLmax,nr_pts=nr_points,fast_approx=True,show_plots=False,AR_coating=AR_coating,plot_r=plot_r)
#         # print avg_ZPL
#         avg_ZPLs[i]=avg_ZPL
#         t1 = time.time()-t0
#         if not analytic_version:
#             print '%d out of %d done'%(i+1, len(t_d_s))
#             print 'expected time remaining: %d seconds'%(t1/(i+1)*len(t_d_s)-t1)


#     if save_data:
#         filename = '%s_avgZPLvstds_R_%.1f_tag_%.2f_vibdL_%.2f,loss_%.2f_%s'%(tb.get_timestamp_from_now()[8:14],R*1.e6,t_a_g*1.e6,vib_dL*1.e9,losses*1.e6,mode)
#         data_dict = {}
#         data_dict['vib_dL']=vib_dL
#         data_dict['avg_ZPLs']=avg_ZPLs
#         data_dict['R']=R
#         data_dict['t_d_s']=t_d_s
#         data_dict['t_a_g']=t_a_g
#         data_dict['dnu']=dnu
#         data_dict['dLmax'] = dLmax
#         data_dict['nr_pts'] = nr_pts
#         data_dict['mode'] = mode
#         data_dict['losses'] = losses
#         save_to_json_file(data_folder,filename,data_dict)


#     fig,ax = plt.subplots()
#     ax.plot(t_d_s*1.e6,avg_ZPLs)
#     ax.set_xlabel('diamond thickness (um)')
#     ax.set_ylabel('avg emission into ZPL')
#     plt.show()
#     plt.close()

#     return avg_ZPLs



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
    Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_2, 2.94525*L_2,rho_even)
    Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_1, 0.94808*L_1,rho_odd)
    Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_2, 2.99742*L_2,rho_even)

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





