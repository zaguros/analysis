##############################################
# A.M.J. Zwerver
# a.m.j.zwerver@student.tudelft.nl
#adapted by SvD
##############################################
# Simulation of the resonant modes of a cavity with a diamond membrane

import numpy as np
from scipy import interpolate as intp
import scipy as sc
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import pi,sqrt,sin
import collections
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot as tools_plot
from analysis.lib.tools import toolbox as tb
import os


from pylab import *




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


c = sc.constants.c # speed of light
t_d = 3.e-6
n_air = 1.0
n_diamond = 2.41
n_H = 2.146 #n_H (Ta2O5)
n_L = 1.4585 # n_L
number_of_layers = 2
na = n_air
nb = n_air
lab0 = 637e-9
pi=math.pi

rho_in = (na - n_H)/(na + n_H)
rho = (n_H - n_L)/(n_H + n_L)
rho_out = (n_H - nb)/(n_H + nb) ##n_(i-1)-n_i/(n_(i-1)+n_i)
rho_ad = (n_air-n_diamond)/(n_diamond+n_air)

l_H = lab0/(4*n_H)
l_L = lab0/(4*n_L)

D_ad = (1/(1+rho_ad)) * np.array([[1, rho_ad] , [rho_ad,1]])
D_da = linalg.inv(D_ad)


class electric_field_cavity():

    def __init__(self, t_d, t_a_g, R, **kw):
        self.t_d = t_d
        self.t_a_g = t_a_g
        self.R = R
        self.lambda_cav = kw.pop('lambda_cav', 637.e-9)
        self.n_1 = kw.pop('n_1', 2.15) #(Ta2O5)
        self.n_2 = kw.pop('n_2', 1.46) #SiO2
        self.M = kw.pop('M', 10) #number of mirror layers
        self.type = kw.pop('type','hybrid') #cavity type. can be hybrid or air

    def boundary_matrix(self,n1,n2):
        rho = (n1-n2)/(n1+n2)
        tau = 1+rho
        matrix = (1./tau)*np.array(((1,rho),(rho,1)))
        return matrix

    def propagation_matrix(self,n,t,lambda_i):
        k_c = 2*math.pi*n/lambda_i
        matrix = np.array(((exp(1j*k_c*t),0),(0,exp(-1j*k_c*t))))
        return matrix

    def electric_field_from_transfer_matrix(self,matrix,E):
        out = np.dot(matrix,E)
        return out  

    def ns_in_cavity(self):
        if type == 'hybrid':
            ns_mirror = np.append(np.array((self.n_1,self.n_2)*self.M),np.array((self.n_1)))
            ns_cav = np.array((n_diamond, n_air))
            self.ns = np.concatenate((ns_mirror,ns_cav,ns_mirror))
        elif type == 'air':
            ns_mirror = np.append(np.array((self.n_1,self.n_2)*M),np.array((n_1)))
            ns_cav = np.array([n_air])
            self.ns = np.concatenate((ns_mirror,ns_cav,ns_mirror)) 
        else:
            'specify valid cavity type (hyrbid or air)!'
        return self.ns

    def ts_in_cavity(self):
        if type == 'hybrid':
            ts_mirror = np.append(np.array((self.lambda_cav/(4*self.n_1),self.lambda_cav/(4*self.n_2))*self.M),np.array(self.lambda_cav/(4*self.n_1)))
            ts_cav = np.array((self.t_d,self.t_a))
            self.ts = np.concatenate((ts_mirror,ts_cav,ts_mirror))
        elif type == 'air':
            ts_mirror = np.append(np.array((self.lambda_cav/(4*self.n_1),self.lambda_cav/(4*self.n_2))*self.M),np.array(self.lambda_cav/(4*self.n_1)))
            ts_cav = np.array([self.t_a])
            self.ts = np.concatenate((ts_mirror,ts_cav,ts_mirror))
        return self.ts


    def electric_field_distribution(self,na=n_air,nb=n_air,lambda_i = 637.e-9,nr_points=30001):
        """
        calculates the electric field disrtibution, assuming there is no outgoing field.
        input:
        ts          lengths of all elements. right-most element first
        ns          refractive indices of all elements. right-most element first 
        na          refractive index of left-most element
        nb          refractive index of right-most element
        lambdas     wavelengths at which to evaluate
        """
        self.L=np.sum(self.ts)
        self.zs = np.linspace(0,self.L,nr_points)
        m1 = boundary_matrix(self.ns[0],na)
        mf = boundary_matrix(nb,self.ns[-1])
        self.ns_vs_z = np.zeros(nr_points)
        Efw_vs_z = np.zeros((nr_points),dtype='complex')
        Ebw_vs_z= np.zeros((nr_points),dtype='complex')
        
        t_tot=0

        Er = np.transpose(np.array((1,0),dtype='complex'))

        for i,(t,n) in enumerate(zip(ts,ns)):
            ki =  2*math.pi*n/lambda_i
            if i==0:
                c = np.identity(2)
                m = m1
                E=Er
            else:
                m = boundary_matrix(ns[i],ns[i-1])
            c = np.dot(m,c)
            E = electric_field_from_transfer_matrix(m,E)

            i_ts = np.where(((zs <= L-t_tot) & (zs > L-t_tot-t)))
            self.ns_vs_z[i_ts] = n
            Efw_vs_z[i_ts] = E[0]*exp(1j*(L-t_tot-zs[i_ts])*2*math.pi*n/lambda_i)#E2_fw*exp(1j*(L-t_tot-zs[i_t])*2.*math.pi*n/lambda_i)
            Ebw_vs_z[i_ts] = E[1]*exp(-1j*(L-t_tot-zs[i_ts])*2*math.pi*n/lambda_i)#E2_bw*exp(-1j*(L-t_tot-zs[i_t])*2.*math.pi*n/lambda_i)

            p_i = propagation_matrix(n,t,lambda_i)
            c = np.dot(p_i,c)

            E = electric_field_from_transfer_matrix(p_i,E)


            t_tot=t_tot+t

        c = np.dot(mf,c)
        
        self.Etot_vs_z = Efw_vs_z+Ebw_vs_z #E
        H_vs_z = (Efw_vs_z-Ebw_vs_z)*n #H
        Z_vs_z = Etot_vs_z/ H_vs_z# 

        r = c[1,0]/c[0,0] #reflectivity
        Etot_vs_z=np.abs(Etot_vs_z/c[0,0])
        H_vs_z=np.abs(H_vs_z/c[0,0])

        return zs,ns_vs_z,Etot_vs_z,H_vs_z,Z_vs_z,r





def trans_per_element(matrix):
    t = matrix[0][0]-(matrix[0][1]*matrix[1][0])/matrix[1][1]
    T =np.conjugate(t)*t
    return T

def refl_per_element(matrix):
    t = -matrix[1][0]/matrix[1][1]
    T =np.conjugate(t)*t
    return T





def electric_field_distribution(ts,ns,na,nb,lambda_i = 637.e-9,nr_points=31):
    """
    calculates the electric field disrtibution, assuming there is no outgoing field.
    input:
    Ls          lengths of all elements. right-most element first
    ns          refractive indices of all elements. right-most element first 
    na          refractive index of left-most element
    nb          refractive index of right-most element
    lambdas     wavelengths at which to evaluate
    """
    L=np.sum(ts)
    zs = np.linspace(0,L,nr_points)
    m1 = boundary_matrix(ns[0],na)
    mf = boundary_matrix(nb,ns[-1])
    ns_vs_z = np.zeros(nr_points)
    # E_vs_z = np.zeros((nr_points,2),dtype='complex')
    Efw_vs_z = np.zeros((nr_points),dtype='complex')
    Ebw_vs_z= np.zeros((nr_points),dtype='complex')
    # Etot_vs_z= np.zeros((nr_points),dtype='complex')
    # H_vs_z= np.zeros((nr_points),dtype='complex')
    # Z_vs_z= np.zeros((nr_points),dtype='complex')
    
    t_tot=0

    Er = np.transpose(np.array((1,0),dtype='complex'))

    for i,(t,n) in enumerate(zip(ts,ns)):
        ki =  2*math.pi*n/lambda_i
        if i==0:
            c = np.identity(2)
            m = m1
            E=Er
        else:
            m = boundary_matrix(ns[i],ns[i-1])
        c = np.dot(m,c)
        E = electric_field_from_transfer_matrix(m,E)

        i_ts = np.where(((zs <= L-t_tot) & (zs > L-t_tot-t)))
        ns_vs_z[i_ts] = n
        Efw_vs_z[i_ts] = E[0]*exp(1j*(L-t_tot-zs[i_ts])*2*math.pi*n/lambda_i)#E2_fw*exp(1j*(L-t_tot-zs[i_t])*2.*math.pi*n/lambda_i)
        Ebw_vs_z[i_ts] = E[1]*exp(-1j*(L-t_tot-zs[i_ts])*2*math.pi*n/lambda_i)#E2_bw*exp(-1j*(L-t_tot-zs[i_t])*2.*math.pi*n/lambda_i)

        p_i = propagation_matrix(n,t,lambda_i)
        c = np.dot(p_i,c)

        E = electric_field_from_transfer_matrix(p_i,E)


        t_tot=t_tot+t

    c = np.dot(mf,c)
    
    Etot_vs_z = Efw_vs_z+Ebw_vs_z #E
    H_vs_z = (Efw_vs_z-Ebw_vs_z)*n #H
    Z_vs_z = Etot_vs_z/ H_vs_z# 

    r = c[1,0]/c[0,0] #reflectivity
    Etot_vs_z=np.abs(Etot_vs_z/c[0,0])
    H_vs_z=np.abs(H_vs_z/c[0,0])

    # E_vs_z = E_vs_z/c[0,0]#and multiply by E0+ = 1 

    return zs,ns_vs_z,Etot_vs_z,H_vs_z,Z_vs_z,r



def electric_field_in_hybrid_cavity(t_a,t_d,n_1,n_2,M,lambda_c = 637.e-9,lambda_i=637e-9,nr_points=31):
    """
    calculate electric field distribution in a hybrid diamond-air cavity
    """
    ns = ns_in_hybrid_cavity(n_1,n_2,M)
    ts = ts_in_hybrid_cavity(n_1,n_2,M,lambda_c,t_d,t_a)
    zs,ns_vs_z,Etot_vs_z,H_vs_z,Z_vs_z,r=electric_field_distribution(ts,ns,n_air,n_air,lambda_i=lambda_i,nr_points=nr_points)

    return ts,ns,zs,ns_vs_z,Etot_vs_z,H_vs_z,Z_vs_z,r



def electric_field_in_air_cavity(t_a,n_1,n_2,M,lambda_c = 637.e-9,lambda_i=637e-9,nr_points=31):
    """
    calculate electric field distribution in a cavity
    """
    ns = ns_in_air_cavity(n_1,n_2,M)
    ts = ts_in_air_cavity(n_1,n_2,M,lambda_c,t_a)
    zs,ns_vs_z,Etot_vs_z,H_vs_z,Z_vs_z,r=electric_field_distribution(ts,ns,n_air,n_air,lambda_i=lambda_i,nr_points=nr_points)

    return ts,ns,zs,ns_vs_z,Etot_vs_z,H_vs_z,Z_vs_z,r


def cavity_reflectivity(ts,ns,na,nb,lambda_i = 637.e-9):
    """
    calculates the cavity reflectivity, assuming there is no outgoing field.
    input:
    Ls          lengths of all elements. right-most element first
    ns          refractive indices of all elements. right-most element first 
    na          refractive index of left-most element
    nb          refractive index of right-most element
    lambdas     wavelengths at which to evaluate
    """
    L=np.sum(ts)
    m1 = boundary_matrix(ns[0],na)
    mf = boundary_matrix(nb,ns[-1])

    Er = np.transpose(np.array((1,0)))

    for i,(t,n) in enumerate(zip(ts,ns)):

        ki =  2*math.pi*n/lambda_i
        if i==0:
            c = np.identity(2)
            m = m1
            E=Er
        else:
            m = boundary_matrix(ns[i],ns[i-1])
        c = np.dot(m,c)
 
        pi = propagation_matrix(n,t,lambda_i)
        c = np.dot(pi,c)

    c = np.dot(mf,c)
    r = c[1,0]/c[0,0] #reflectivity

    return r


def max_Evac(mode_volume,n_max,lambda_cav=637.e-9):
    """
    Uses the mode volume to determine the maximum vacuum electric field (hbar omega/(epsilon*V)) (epsilon = n^2*epsilon_0)
    """
    E = np.sqrt(sc.constants.hbar*2*math.pi*(sc.constants.c/lambda_cav)/(2*sc.constants.epsilon_0*mode_volume))/n_max
    return E

def calculate_E_max_in_n(zs,Etot_vs_z,ns_vs_z,lambda_cav=637.e-9,n_z0=n_diamond):
    """
    calculate the maximum electric field density and its location, in the cavity region with refractive index n_z0
    """
    dz = np.abs(zs[0]-zs[1])
    arg_zs_in_cavity = np.where((ns_vs_z<n_z0+0.001)&(ns_vs_z>n_z0-0.001))
    zs_in_cavity = zs[arg_zs_in_cavity]
    arg_shallow_zs_in_cavity = arg_zs_in_cavity[0][-int(lambda_cav/2/n_z0/dz):]
    z0_i = arg_shallow_zs_in_cavity[0] + np.argmax(Etot_vs_z[arg_shallow_zs_in_cavity])
    z0 = zs[z0_i]
    E_z0 = Etot_vs_z[z0_i]
    dz0 = np.abs(zs_in_cavity[-1]-z0)
    
    # print 'index corresponding to z0',z0_i 
    # print 'z0',z0*1.e6
    # print 'dz0 = ',dz0*1.e9,'nm = ',dz0*n_z0/lambda_cav,'x lambda'

    return E_z0,z0,dz0

def calculate_energy_dist_length(zs,Etot_vs_z,ns_vs_z,E_z0,n_z0):
    """
    calculate energy distribution length int((n(z)^2)*(E(z))^2 dz)/(n(z0))^2*E(z0)^2
    """
    dz = np.abs(zs[0]-zs[1])
    length = np.real(np.sum((ns_vs_z**2)*(Etot_vs_z**2))*dz/(n_z0**2*E_z0**2))
    return 2*length #note factor 2 that seems necessary... 
    
def calculate_energy_in_m(zs,Etot_vs_z,ns_vs_z,n_m):
    arg_zs_in_m = np.where((ns_vs_z<n_m+0.01)&(ns_vs_z>n_m-0.01))
    ns_vs_z_in_m = ns_vs_z[arg_zs_in_m]
    Etot_vs_z_in_m = Etot_vs_z[arg_zs_in_m]
    dz = np.abs(zs[0]-zs[1])
    Etot = 1./2*sc.constants.epsilon_0*np.sum((ns_vs_z_in_m**2)*(Etot_vs_z_in_m**2))*dz 
    return Etot #in J/m^2 

def calculate_mode_volume(beam_waist,effective_length):
#     beam_waist = cav_sim.calc_waist()
    mode_volume = np.pi*beam_waist**2/4.*effective_length
    return mode_volume 

def calculate_w0(L,R,lambda_cav):
    """
    calculate the physical beam waist in cavity part with refr index n_z0
    parameters:
    L  optical cavity length
    R  radius of curvature
    """
    return ((lambda_cav/np.pi)**0.5)*(L*(R-L))**(1/4.)


def analyse_cavity(zs,Etot_vs_z,ns_vs_z,lambda_cav,n_z0,R,folder,t_d,t_a,plot_Evac=False):
    E_z0,z0,dz0 = calculate_E_max_in_n(zs,Etot_vs_z,ns_vs_z,lambda_cav=lambda_cav,n_z0=n_z0)
    eff_cav_length = calculate_energy_dist_length(zs,Etot_vs_z,ns_vs_z,E_z0,n_z0)

    beam_waist = calculate_w0(eff_cav_length,R,lambda_cav)
    mode_volume = calculate_mode_volume(beam_waist,eff_cav_length)#t_a+t_d*2.4)
    
    Evac_max = max_Evac(mode_volume,n_z0,lambda_cav)
    Etot_vs_z = Etot_vs_z*Evac_max/E_z0

    Etot_diamond = calculate_energy_in_m(zs,Etot_vs_z,ns_vs_z,n_diamond)
    Etot_air = calculate_energy_in_m(zs,Etot_vs_z,ns_vs_z,n_air)
    

    if plot_Evac:
        print 'energy distribution length = ',eff_cav_length*1.e6, 'um'
        print 'beam waist', beam_waist*1.e6, 'um'
        print 'mode volume', mode_volume*(1.e6)**3, 'um^3'
        print 'max Evac',Evac_max/1000.,'kV/m' 
        today = tb.get_timestamp_from_now()[:8]
        title_string = today+'_td_%.2f_ta_%.2f_lambdacav_%.1f_R_%.1f'%(t_d*1.e6,t_a*1.e6,lambda_cav*1.e9,R*1.e6)

        fig,ax = plt.subplots(figsize=(12,6))
        ax.plot(zs*1.e6,ns_vs_z)
        ax.set_ylabel('n')
        ax.set_xlabel('z (um)')
        ax.set_title(title_string)
        ax2=ax.twinx()
        ax2.plot(zs*1.e6,Etot_vs_z/1000.,'g')
        ax2.plot([z0*1.e6],[Evac_max/1000.],'ok',label='dz = %.1f nm, Emax = %.1f kV/m'%(dz0*1.e9,E_z0))
        ax2.set_ylabel('E (kV/m)')
        ax2.legend()
        plt.savefig(folder+title_string+'.png')

    return zs,ns_vs_z,Etot_vs_z,Evac_max,z0,dz0,Etot_diamond,Etot_air

def analyse_hybrid_cavity(t_a,t_d,n_1,n_2,lambda_cav,R,folder,nr_pts=30001,plot_Evac=False):
    ts,ns,zs,ns_vs_z,Etot_vs_z,H_vs_z,Z_vs_z,r = electric_field_in_hybrid_cavity(t_a,t_d,n_1,n_2,10,lambda_i=lambda_cav,nr_points=nr_pts)
    return analyse_cavity(zs,Etot_vs_z,ns_vs_z,lambda_cav,n_diamond,R,folder,t_d,t_a,plot_Evac=plot_Evac)

def analyse_air_cavity(t_a,n_1,n_2,lambda_cav,R,folder,nr_pts=30001,plot_Evac=False):
    ts,ns,zs,ns_vs_z,Etot_vs_z,H_vs_z,Z_vs_z,r = electric_field_in_air_cavity(t_a,n_1,n_2,10,lambda_i=lambda_cav,nr_points=nr_pts)
    return analyse_cavity(zs,Etot_vs_z,ns_vs_z,lambda_cav,n_air,R,folder,t_d,t_a,plot_Evac=plot_Evac)

def find_res_condition(t_a_g,t_d,lambda_cav,nr_pts=161,plot_r=False):
    ns = ns_in_hybrid_cavity(n_1,n_2,10)

    t_as = np.linspace(t_a_g-100e-9,t_a_g+100e-9,nr_pts)
    r_is = np.zeros(len(t_as))
    for i,t_ai in enumerate(t_as):
        ts = ts_in_hybrid_cavity(n_1,n_2,10,lambda_cav,t_d,t_ai)
        r = cavity_reflectivity2(ts,ns,n_air,n_air)
        r_is[i]=np.abs(r)**2

    t_ai=t_as[np.argmin(r_is)]
    r_i = r_is[np.argmin(r_is)]
    
    if plot_r:
        fig,ax = plt.subplots()
        ax.plot(t_as*1.e6,r_is)
        plt.show()
        plt.close()
    
    t_as = np.linspace(t_ai-1e-9,t_ai+1e-9,4*nr_pts)
    r_iis = np.zeros(len(t_as))
    for i,t_aii in enumerate(t_as):
        ts = ts_in_hybrid_cavity(n_1,n_2,10,lambda_cav,t_d,t_aii)
        r = cavity_reflectivity2(ts,ns,n_air,n_air)
        r_iis[i]=np.abs(r)**2


    t_aii=t_as[np.argmin(r_iis)]
    r_ii = r_iis[np.argmin(r_iis)]
    if plot_r:
        print 'air gap',t_aii*1.e6, ' um; r = ', r_ii

        fig,ax = plt.subplots()
        ax.plot(t_as*1.e6,r_iis)
        plt.show()
        plt.close()
    
    return t_aii,r_ii,r_iis#,t_ai,r_i,r_is#,




def dielmirror_recursive_refl_coeff(Gamma_ip1, k_i, l_i,rho_i):
    Gamma_i = (rho_i+Gamma_ip1*exp(-2j*k_i*l_i) )/(1+rho_i*Gamma_ip1*exp(-2j*k_i*l_i) )
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

# def cavity_reflectivity(t_d,t_a,n_1,n_2,M,lambdas):
#     rho_AD = (n_air - n_diamond)/(n_air+n_diamond)
#     ks_d = 2*math.pi*n_diamond/lambdas
#     ks_a = 2*math.pi*n_air/lambdas

#     Gamma,R = dielmirror_reflectivity(n_1,n_2,M,lambdas,na=n_air,nb=n_diamond)

#     Gamma = dielmirror_recursive_refl_coeff(Gamma, ks_d, t_d, rho_AD)

#     #propagation of gamma through air
#     Gamma = Gamma*np.exp(-2j*ks_a*t_a)

#     #second mirror
#     Gamma,R = dielmirror_reflectivity(n_1,n_2,M,lambdas,na=n_air,nb=n_air)

#     Refl = np.abs(Gamma)**2
#     return Gamma, Refl


def calculate_transfer_matrices_mirrors(wavelength, lambda_design=637e-9,n_air=n_air,n_H=n_H,n_L=n_L,na=n_air,nb=n_air,number_of_layers=number_of_layers):
    k_H = 2.*pi*n_H/wavelength
    k_L = 2.*pi*n_L/wavelength
    t_H = lambda_design/(4*n_H)
    t_L = lambda_design/(4*n_L)
    v_H =np.exp(1j*k_H*t_H)


    rho = (n_H-n_L)/(n_H+n_L)
    rho_in = (na-n_H)/(na+n_H)
    rho_out = (n_H-nb)/(n_H+nb)

    B1 = 1./(1+rho)*np.array(((1,rho),(rho,1)))
    P1 = np.array(((exp(1j*k_H*t_H),0),(0,exp(-1j*k_H*t_H))))
    B2 = 1./(1-rho)*np.array(((1,-rho),(-rho,1)))
    P2 = np.array(((exp(1j*k_L*t_L),0),(0,exp(-1j*k_L*t_L))))

    F = np.dot(B2,P2)
    F = np.dot(P1,F)
    F = np.dot(B1,F)
    # print 'F style 1 \n',F2
    # F = (1./(1.-rho**2))*np.array([[(v_H*v_L)-rho**2*v_H/v_L, (-2j*rho/v_H)*np.sin(k_L*l_L)],[2j*rho*v_H*np.sin(k_L*l_L), 1/(v_H*v_L)-((rho**2)*v_L/v_H)]])   
    # print 'F style 2\n',F

    B_in = 1./(1+rho_in)*np.array(((1,rho_in),(rho_in,1)))
    P_in = np.array(((exp(1j*k_H*t_H),0),(0,exp(-1j*k_H*t_H))))
    F_in = np.dot(B_in,P_in)
    # F_in = (1./(1.+rho_in))*np.array(((v_H, rho_in/v_H),(rho_in*v_H, 1/v_H)))
    F_out = (1./(1.+rho_out))*np.array(((1.,rho_out),(rho_out,1.)))
    # print 'F style 1\n',F_in
    # print 'F style 2 \n',F_in2

    N = number_of_layers
    # add all layers of the mirror. 
    H = F_in
    for i in range(N):
        H = np.dot(H,F)

    M_gma = np.dot(H,F_out)
    M_amg = M_gma#linalg.inv(M_gma)

    return(M_gma,M_amg)

def mirror_reflectivity_from_matrix(lambdas,number_of_layers=10):

    r_coeff = np.zeros(len(lambdas),dtype='complex')
    R = np.zeros(len(lambdas),dtype='complex')

    for i,wavelength in enumerate(lambdas):
        M_gma,M_amg = calculate_transfer_matrices_mirrors(wavelength=wavelength,number_of_layers=number_of_layers)
        # print M_gma
        # print M_gma[0,0]
        # print M_gma[1,0]
        # print M_gma[0,0]/M_gma[0,1]
        r_coeff[i] = M_gma[1,0]/M_gma[0,0]
        # print M_gma[0,1]
        # print r_coeff[i]
        R[i] = np.abs(r_coeff[i])**2
        R[i] = np.conjugate(r_coeff[i])*r_coeff[i]
    return R,r_coeff


def calculate_cavity_transmission(wavelengths, lengths, t_d = 3.e-6, D_ad = D_ad,D_da = D_da,n_diamond=n_diamond):
    Ts=np.zeros((len(wavelengths),len(lengths)))
    # Ts_inv=np.zeros((len(wavelengths),len(lengths)))

    for ii,wavelength in enumerate(wavelengths):
        M_gma, M_amg = calculate_transfer_matrices_mirrors(wavelength)

        T1 =1-0.999865  #tranmission of the mirror on the air side
        T2 =1-0.9999 #ranmission of the mirror on the diamond side

        M_amg = 1j/np.sqrt(np.real(T1))*np.array([[-1,-np.sqrt(1-np.real(T1))],[np.sqrt(1-np.real(T1)),1]]) #tranmission matrix on air side, no losses
        M_gma =  1j/np.sqrt(np.real(T2))*np.array([[1,np.sqrt(1-np.real(T2))],[-np.sqrt(1-np.real(T2)),-1]]) #transmission matrix mirror on diamond side, no losses

        for j,l in enumerate(lengths):

            L_air = np.array([[np.exp(-2j*pi*(l)/wavelength),0],[0,np.exp(2j*pi*(l)/wavelength)]])
            L_diamond = np.array([[np.exp(-2j*pi*n_diamond*t_d/wavelength),0],[0,np.exp(2j*pi*n_diamond*t_d/wavelength)]])

            L_air2 = linalg.inv(L_air)
            L_diamond2 = linalg.inv(L_diamond)
            
            S = np.identity(2)

            for i in [M_amg, L_air, D_da, L_diamond, D_ad, M_gma]:
                S = np.dot(S,i)

            t_small = S[0][0]-(S[0][1]*S[1][0])/S[1][1]
            T = np.conjugate(t_small)*t_small
            Ts[ii,j]= T   
            
            # Sinv = linalg.inv(S)            
            # t_small_inv = Sinv[0][0]-(Sinv[0][1]*Sinv[1][0])/Sinv[1][1]
            # T_inv = np.conjugate(t_small_inv)*t_small_inv
            # Ts_inv[ii,j]= T_inv


    return Ts

def calculate_mirror_reflectivity(wavelengths,number_of_layers=number_of_layers):
    M_gmas = np.zeros((len(wavelengths),2,2))
    M_amgs = np.zeros((len(wavelengths),2,2))

    #print np.shape(M_gmas)
    for i,wavelength in enumerate(wavelengths):
        M_gmas[i], M_amgs[i] = calculate_transfer_matrices_mirrors(wavelength, number_of_layers=number_of_layers)
    print M_gmas
    print M_amgs
    #print M_gmas
    #print M_gmas[:,0,0]
    #E_out = np.array([[1],[0]])
    #E_outs = np.tile(E_out,(len(wavelengths),1,1))
    ##print np.shape(E_outs)
    #print E_outs
    #print E_out
    #print M_gmas[0]
    #print E_outs[0]
    # E_in = np.dot(M_gmas[0],E_outs[0])
    #print E_in
    #M_gmas_d = np.dot(D_ad,M_gmas)
    Ts   = (np.conjugate(M_gmas[:,0,1])/M_gmas[:,0,0])
    Ts2   = (np.conjugate(M_amgs[:,0,1])/M_amgs[:,0,0])

    #t_small_to_a = M_gmas[:,0,0]-(M_gmas[:,0,1]*M_gmas[:,1,0])/M_gmas[:,1,1]
    #t_small_to_a = M_amgs[:,0,0]-(M_amgs[:,0,1]*M_amgs[:,1,0])/M_amgs[:,1,1]
    #t_small_to_d = M_gmas_d[:,0,0]-(M_gmas_d[:,0,1]*M_gmas_d[:,1,0])/M_gmas_d[:,1,1]
    #Ts = np.conjugate(t_small_to_a)*t_small_to_a
    fig,ax = plt.subplots()
    ax.plot(wavelengths,Ts)
    ax.plot(wavelengths,Ts2)
    return Ts

    # print E_in

    #r = E_in[1][0]/E_in[0][0]

def extract_res_frequencies(Ts,wavelength2,length,frequency2,plot_Tss=True, plot_fits=True,sweep_td=False,t_ds=t_d):
    # Xs = []
    # Ys = []
    # LWs = []

    linewidths = np.zeros(np.shape(np.transpose(Ts)))
    Qs = np.zeros(np.shape(np.transpose(Ts)))
    Fs = np.zeros(np.shape(np.transpose(Ts)))
    FSRs = np.zeros(np.shape(np.transpose(Ts)))

    offset = 0.
    amplitude = 1.
    gamma = 0.008 # in THz
    freqfit_old = 0.

    for jj,Tss in enumerate(np.transpose(Ts)):#look at Ts for each length
        length_j = [length[jj]]
        if sweep_td:
            td_j = t_ds[jj]
        else:
            td_j = t_ds
        if plot_Tss:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(frequency2,np.exp(Tss),'o')
            #ax.set_xlim(473.8,474)
        maxima = np.r_[True, Tss[1:] > Tss[:-1]] & np.r_[Tss[:-1] > Tss[1:], True] #returns True for local maximum
        Tss_maxima = Tss*maxima
        Tss_maxima_indices = np.nonzero(Tss_maxima)[0]
        # remove the end-indices from the list of local maxima - they are artifacts:
        Tss_maxima_indices=Tss_maxima_indices[ (Tss_maxima_indices !=0) & (Tss_maxima_indices != (len(frequency2)-1)) ]
        
        for k,index in  enumerate(Tss_maxima_indices):
            frequency2_part = np.linspace(frequency2[index]-0.08,frequency2[index]+0.08,301)
            wavelength2_part = c/frequency2_part*1.e-12
            Tss_part = calculate_cavity_transmission(wavelength2_part, length_j,t_d=td_j)[:,0]

            guess_freq = frequency2[index]
            x0 = guess_freq
            fixed = [offset]

            p0, fitfunc, fitfunc_str = common.fit_lorentz(offset, amplitude, x0, gamma)
            fit_result = fit.fit1d(frequency2_part,Tss_part, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
            
            freqfit = fit_result['params_dict']['x0']
            gammafit = fit_result['params_dict']['gamma']
            gammaerror = fit_result['error_dict']['gamma']
            chisq = fit_result['chisq']

            if plot_fits or chisq > 1e-5:
                ax2 = tools_plot.plot_fit1d(fit_result, np.linspace(frequency2_part[0],frequency2_part[-1],10*len(frequency2_part)), ret='ax', ax=None,color='r',show_guess=True, plot_data=True, label = 'fit', add_txt=False)
                ax2.set_xlabel('frequency (THz)', fontsize=14)
                ax2.set_ylabel("Transmission", fontsize=14)
                ax2.set_title('Cavity length is '+ str(round(length[jj]*1e6,2)) + ' um \n Linewidth is ' + str(round(gammafit*1000,3)) + '$\pm$ '+ str(round(gammaerror*1000,2)) + ' GHz')
                plt.show()
                plt.close()
            if chisq > 1.e-5:
                print 'fit didnt go well: chisq = ', chisq
                continue

            # Xs = np.append(Xs,length_j[0])
            # Ys = np.append(Ys,freqfit)
            # LWs = np.append(LWs,gammafit)

            line_width = int(len(frequency2)*gammafit/(frequency2[0]-frequency2[-1]))+10
            for l in np.arange(2*line_width+1) - line_width:
                if 0 < index+l < len(frequency2):
                    linewidths[jj,index+l] = gammafit*1.e3 #put it at the position of freq_guess, in GHz. Ideally, fix this
                    Qs[jj,index+l] = frequency2[index]/gammafit
                    if freqfit > freqfit_old:
                        FSRs[jj,index+l] = freqfit - freqfit_old
                        Fs[jj,index+l] = (freqfit - freqfit_old)/(gammafit)
            freqfit_old = freqfit

    return np.transpose(linewidths), np.transpose(Qs), np.transpose(Fs), np.transpose(FSRs)#, Xs, Ys, LWs


def extract_res_lengths(Ts,wavelength2,length,plot_Tss=True, plot_fits=True,t_d = t_d):
    # Xs = []
    # Ys = []
    # LWs = []

    linewidths = np.zeros(np.shape(Ts))#linewidths = np.zeros(np.shape(np.transpose(Ts)))
    Fs = np.zeros(np.shape(Ts))#Qs = np.zeros(np.shape(np.transpose(Ts)))
    FSRs = np.zeros(np.shape(Ts))#FSRs = np.zeros(np.shape(np.transpose(Ts)))

    offset = 0.
    amplitude = 1.
    gamma = 0.3e-9 # in m
    lengthfit_old = 0 #dummy number

    for jj,Tss in enumerate(Ts):#look at Ts for each frequency
        wl_j = [wavelength2[jj]]
        
        if plot_Tss:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(length,np.exp(Tss),'o')
            plt.show()
            plt.close()

        maxima = np.r_[True, Tss[1:] > Tss[:-1]] & np.r_[Tss[:-1] > Tss[1:], True] #returns True for local maximum
        Tss_maxima = Tss*maxima
        Tss_maxima_indices = np.nonzero(Tss_maxima)[0]
        # remove the end-indices from the list of local maxima - they are artifacts:
        Tss_maxima_indices=Tss_maxima_indices[ (Tss_maxima_indices !=0) & (Tss_maxima_indices != (len(length)-1)) ]

        for k,index in  enumerate(Tss_maxima_indices[:]):
            length_part = np.linspace(length[index]-0.5e-9,length[index]+0.5e-9,151)
            Tss_part = calculate_cavity_transmission(wl_j, length_part,t_d = t_d)[0,:]
            guess_len = length[index]
            x0 = guess_len
            fixed = [offset]
            p0, fitfunc, fitfunc_str = common.fit_lorentz(offset, amplitude, x0, gamma)

            fit_result = fit.fit1d(length_part,Tss_part, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
            lengthfit = fit_result['params_dict']['x0']
            gammafit = fit_result['params_dict']['gamma']
            gammaerror = fit_result['error_dict']['gamma']
            chisq = fit_result['chisq']

            if plot_fits or chisq > 1.e-5:
                print chisq
                ax2 = tools_plot.plot_fit1d(fit_result, np.linspace(length_part[0],length_part[-1],10*len(length_part)), ret='ax', ax=None,color='r',show_guess=True, plot_data=True, label = 'fit', add_txt=False)
                ax2.set_xlabel('length (um)', fontsize=14)
                ax2.set_ylabel("Transmission", fontsize=14)
                ax2.set_title('Frequency is '+ str(round((c/wl_j[0]*1.e-12-470.4)*1.e3,2)) + ' GHz \n Linewidth in length is ' + str(round(gammafit*1.e9,3)) + '$\pm$ '+ str(round(gammaerror*1.e9,3)) + ' nm')
                plt.show()
                plt.close()
            if chisq > 1.e-5:
                print 'fit didnt go well: chisq = ', chisq
                continue

            #Xs = np.append(Xs,lengthfit)
            #Ys = np.append(Ys,wl_j[0])
            #LWs = np.append(LWs,gammafit)

            line_width = int(len(length)*gammafit/(length[0]-length[-1]))+10
            for l in np.arange(2*line_width+1) - line_width:
                if 0 < index+l < len(length):
                    linewidths[jj,index+l] = gammafit*1.e9 #put it at the position of lengthguess. Ideally, fix this
                    if lengthfit > lengthfit_old:
                        FSRs[jj,index+l] = lengthfit - lengthfit_old
                        Fs[jj,index+l] = (lengthfit - lengthfit_old)/(gammafit)
            lengthfit_old = lengthfit

    return linewidths, FSRs, Fs #, Xs, Ys, LWs#np.transpose(linewidths), np.transpose(Qs), np.transpose(FSRs), Xs, Ys, LWs

def plot_2d(X,Y,Zs,fig_string='',output_linewidth='frequency',xlabel='x',ylabel='y',titlelabel=None,plot_type='lw', exclude_zero = False):
    """
    function that makes a 2d plot, excluding the '0's in the data. This is suitable for plotting properties of resonances.
    PARAMETERS:
    X - X part of meshgrid
    Y - Y part of meshgrid
    Zs - 2d array of data
    fig_string - name for figure saving
    output_linewidth - either 'frequency' or 'length'. default: 'frequency'
    xlabel - string defining the xlabel.
    ylabel - string defining the ylabel
    titlelabel - optional. default: None
    plot_type - the type of data in Zs: either 'lw', 'F', 'Q', 'FSR', 'T'. Default: 'T'
    exclude_zero - if True, excludes the '0's in the data. This is suitable for plotting properties of resonances. Default: 'False'
    """

    if titlelabel == None: #if no titlelabel is specified, base the plot title on the plot type
        if plot_type == 'lw':
            if output_linewidth == 'frequency':
                title_string = 'linewidth (GHz)'
            elif output_linewidth == 'length':
                title_string = 'linewidth (nm)'
        elif plot_type == 'FSR':
            title_string = 'Free Spectral Range'
        elif plot_type == 'F':
            title_string = 'Finesse'
        elif plot_type == 'Q':
            title_string = 'Quality factor'
        elif plot_type == 'T':
            title_string = 'Transmission (log(T))'
        else:
            print 'no good plot_type specified'
            title_string = ''

    else:
        title_string = titlelabel

    print 'plotting '+title_string
    # plot linewidths at resonance vs length lws

    fig = plt.figure()
    ax = fig.add_subplot(111)
    CS = ax.pcolor(X,Y,Zs,vmin = np.min(Zs[nonzero(Zs)]))#,vmax =np.max(FSRs[nonzero(FSRs)]))
    CS.cmap.set_under('w')
    CS.cmap.set_over('k')
    ax.set_xlim(X[0][0],X[-1][-1])
    ax.set_ylim(Y[0][0],Y[-1][-1])
    ax.set_ylabel(ylabel,fontsize=15)
    ax.set_xlabel(xlabel,fontsize=15)
    if plot_type == 'T' or titlelabel == None:
        ax.set_title(title_string,fontsize = 15)
    else:
        ax.set_title(title_string + ' in '+output_linewidth,fontsize=15)
    colorbar = plt.colorbar(CS,ax=ax)

    try:
        fig.savefig(data_folder+'/'+tb.get_timestamp_from_now()[8:]+'__'+plot_type+'__'+fig_string)
    except: 
        print 'Figure not saved'

    plt.show()
    plt.close()


def sweep_t_d():
    """
    Function that varies the diamond thickness, and the frequency of the ingoing light, while keeping the air-part of the cavity length fixed.
    It plots the cavity transmission against these two parameters
    and extracts and plots the linewidth in frequency.
    """

    output_linewidth = 'frequency'

    tdpoints = 5001
    ypoints = 21
    t_ds = np.linspace(1.5e-6,4.5e-6,tdpoints)
    wavelengths = np.linspace(636.e-9,638e-9,ypoints)
    frequencies = c/wavelengths*1.e-12 #in THz
    lengths = np.ones(tdpoints)*1.0e-6

    X,Y = np.meshgrid(t_ds*1.e6,frequencies)
    xlabel = 'diamond length (um) + cavity length {} um'.format(lengths[0]*1.e6)
    ylabel = 'frequency (THz)'

    Ts = np.zeros((len(wavelengths),len(t_ds)))
    for i,t_d in enumerate(t_ds):
        Ts[:,i] = calculate_cavity_transmission(wavelengths, np.ones(1)*lengths[i], t_d = t_d)[:,0]
    
    Ts = np.log(Ts)
    fig_string = 'vs_freq_tds__l_{:.1f}_n1_{:.2f}_n2_{:.2f}_N_{:d}.png'.format(lengths[0]*1.e6,n_H,n_L,number_of_layers)


    plot_2d(X,Y,Ts,fig_string=fig_string,output_linewidth=output_linewidth,xlabel=xlabel,ylabel=ylabel,plot_type='T')

    lws,Qs,FSRs,Fs = extract_res_frequencies(Ts,wavelength2=wavelengths,length=lengths,frequency2=frequencies,sweep_td=True,t_ds=t_ds,plot_Tss=False,plot_fits=False)

    plot_2d(X,Y,lws,fig_string=fig_string,output_linewidth=output_linewidth,xlabel=xlabel,ylabel=ylabel,plot_type='lw', exclude_zero = True)
    #plot_2d(X,Y,Fs,fig_string=fig_string,output_linewidth=output_linewidth,xlabel=xlabel,ylabel=ylabel,plot_type='F', exclude_zero = True)


def fix_f():
    """
    Function that varies the diamond thickness, and the air-part of the cavity length, for a fixed frequency of the ingoing light.
    It plots the cavity transmission against these two parameters
    and extracts and plots the linewidth in frequency.
    """
    output_linewidth = 'frequency'

    tdpoints = 201
    lpoints = 201
    wlpoints = 21
    wlpoint = int((wlpoints-1)/2)
    t_ds = np.linspace(2.5e-6,3.5e-6,tdpoints)
    lengths = np.linspace(1.e-6,2.e-6,lpoints)
    wavelengths = np.linspace(636.e-9,638.e-9,wlpoints)
    frequencies = c/wavelengths*1.e-12 #in THz

    X,Y = np.meshgrid(lengths*1.e6,t_ds*1.e6)
    xlabel = 'air length (um)'
    ylabel = 'diamond length (um)'
    lw_titlelabel = 'Linewidth (GHz) - f = {:.0f} GHz'.format((frequencies[wlpoint]-470.4)*1.e3)
    fig_string = 'vs_l_tds__f_{:.0f}GHz_n1_{:.2f}_n2_{:.2f}_N_{:d}.png'.format((frequencies[wlpoint]-470.4)*1.e3,n_H,n_L,number_of_layers)

    Ts = np.zeros((len(wavelengths),len(lengths),len(t_ds)))
    lws = np.zeros((len(wavelengths),len(lengths),len(t_ds)))
    Qs = np.zeros((len(wavelengths),len(lengths),len(t_ds)))
    FSRs = np.zeros((len(wavelengths),len(lengths),len(t_ds)))
    Fs = np.zeros((len(wavelengths),len(lengths),len(t_ds)))

    for j,t_d in enumerate(t_ds):
        Ts[:,:,j] = calculate_cavity_transmission(wavelengths, lengths, t_d)
    
    Ts = np.log(Ts)
    
    plot_2d(X,Y,Ts[wlpoint,:,:],fig_string=fig_string,output_linewidth=output_linewidth,xlabel=xlabel,ylabel=ylabel,plot_type='T')

    for j,t_d in enumerate(t_ds):
        lws[:,:,j],Qs[:,:,j],FSRs[:,:,j],Fs[:,:,j] = extract_res_frequencies(Ts[:,:,j],wavelength2=wavelengths,length=lengths,frequency2=frequencies,sweep_td=False,t_ds=t_d,plot_Tss=False,plot_fits=False)

    plot_2d(X,Y,lws[wlpoint,:,:],fig_string=fig_string,output_linewidth=output_linewidth,xlabel=xlabel,ylabel=ylabel,titlelabel = lw_titlelabel,plot_type='lw', exclude_zero = True)

def fix_t_d():
    """
    Function that varies the frequency of the ingoing light, and the air-part of the cavity length, for a fixed diamond length.
    It plots the cavity transmission against these two parameters
    and extracts and plots the linewidth in frequency, or in length.
    """
    t_d = 3.5.e-6
    output_linewidth = 'frequency' #set to length or frequency. 

    if output_linewidth == 'length':
        xpoints = 2001#151 #length
        ypoints = 151#701# #wavelength
    elif output_linewidth == 'frequency':
        ypoints = 301 #length
        xpoints = 2001# #wavelength       

    wavelength2 = np.linspace(700.e-9,600.e-9,ypoints)
    frequency2 = c/wavelength2*1.e-12 #in THz
    length = np.linspace(4e-6,5.75e-6,xpoints)


    X,Y = np.meshgrid(length*1.e6,frequency2)    
    xlabel = 'cavity length (um) + {:.1f} um diamond'.format(t_d*1.e6)
    ylabel = 'frequency (THz)'

    fig_string = 'in_{}_t_d_{:.1f}_n1_{:.2f}_n2_{:.2f}_N_{:d}.png'.format(output_linewidth,t_d*1.e6,n_H,n_L,number_of_layers)


    print 'calculating cavity transmission'
    Ts = calculate_cavity_transmission(wavelength2, length, t_d = t_d)
    Ts = np.log(Ts)

    plot_2d(X,Y,Ts,fig_string=fig_string,output_linewidth=output_linewidth,xlabel=xlabel,ylabel=ylabel,plot_type='T')

    print 'calculating resonances'
    if output_linewidth == 'length':
        lws,FSRs,Fs = extract_res_lengths(Ts,wavelength2=wavelength2,length=length,plot_Tss=False,plot_fits=False,t_d=t_d)#,Xs,Ys,LWs
    elif output_linewidth == 'frequency':
        lws,Qs,FSRs,Fs = extract_res_frequencies(Ts,wavelength2=wavelength2,frequency2=frequency2,length=length,t_ds = t_d,sweep_td=False,plot_Tss=False,plot_fits=False) ##,Xs,Ys,LWs

    plot_2d(X,Y,lws,fig_string=fig_string,output_linewidth=output_linewidth,xlabel=xlabel,ylabel=ylabel,plot_type='lw', exclude_zero = True)
    plot_2d(X,Y,Fs,fig_string=fig_string,output_linewidth=output_linewidth,xlabel=xlabel,ylabel=ylabel,plot_type='F', exclude_zero = True)




