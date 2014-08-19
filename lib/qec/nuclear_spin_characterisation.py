#Module that contains functions to simulate nuclear spin responses to decoupling sequences
import numpy as np

#Physical constants needed in  module
gamma_c = 1.0705e3 #g-factor for C13 in Hz/G

def calc_hyperfine_from_tau(tau_k,k,B_field):
    '''
    calculates the parallel component fo the hyperfine interaction in Hz using formula (3) from Taminiau et al. (2012) PRL
    This formula is valid in the approximation of high magnetic field (Omega_Larmor >> Omega_HF)
    ------
    inputs
    ------
    tau:        tau in seconds where the dip is observed
    k:           estimated order of the dip
    B_field:   Magnetic field strength in Gauss
    ------
    returns
    ------
    HF_par:   parallel component of the hyperfine interaction in Hz
    '''

    omega_L = 2*np.pi *gamma_c*B_field #in Radial freq

    HF_par = ((2*k-1)*np.pi/tau_k-2*omega_L)
    HF_par = HF_par/(2*np.pi) #convert back to freq in Hz

    return HF_par

def calc_tau_from_HF(HF_par,k,B_field):
    '''
    calculates the resonance tau  from the hyperfine interaction in Hz  using formula (3) from Taminiau et al. (2012) PRL
    This formula is valid in the approximation of high magnetic field (Omega_Larmor >> Omega_HF)
    ------
    inputs
    ------
    HF_par:  Parallel component of the HF in Hz
    k:           order of the dip
    B_field:   Magnetic field strength in Gauss
    ------
    returns
    ------
    tau_k:     tau in seconds of dip order k
    '''
    HF_par = HF_par *2*np.pi #convert to radial freq
    omega_L = 2*np.pi *gamma_c*B_field #in radial frequency
    tau_k = (2*k-1)*np.pi/(2*omega_L+HF_par)
    return tau_k

def dyn_dec_signal(HFs_par, HFs_orth, B_field, N ,tau):
    '''
    Takes the HF interaction strengths (paralel and orthogonal), the magnetic field strenght
    and an array of times and returns the signal at those times for that specific spin.
M    ------
    inputs
    ------
    HFs_par:        list of parallel component of HF strength in Hz
    HFs_orth:       list of orthogonal component of HF strength in Hz
    B_field:           Magnetic field in Gauss
    N:                   number of pulses
    tau:                time in s
    -------
    returns
    -------
    M:       list of signals of individual simulated spins
    measured signal is M.prod(axis=0)

    '''
    #Model parameters
    omega_larmor = 2*np.pi*gamma_c*B_field #radial frequency
    tau_larmor = 2*np.pi/omega_larmor #time in seconds


    print 'tau larmor = %s' %tau_larmor

    if np.size(tau)!=1:
        M=np.zeros([np.size(HFs_par),np.size(tau)])
    else:
        M = np.zeros([np.size(HFs_par),np.size(N)])
    for i,HF_par in enumerate(HFs_par):
        HF_par = HF_par*2*np.pi #Convert to radial frequency
        HF_orth = HFs_orth[i]*2*np.pi #convert to radial frequency
        omega_tilde = np.sqrt((HF_par+omega_larmor)**2+HF_orth**2)
        alpha = omega_tilde*tau
        beta = omega_larmor*tau
        mx = HF_orth/omega_tilde
        mz = (HF_par+omega_larmor)/omega_tilde
        vec_term = mx**2 *((1-np.cos(alpha))*(1-np.cos(beta)))/(1+np.cos(alpha)*np.cos(beta)-mz*np.sin(alpha)*np.sin(beta))
        angle_term = np.sin(N*np.arccos(np.cos(alpha)*np.cos(beta)-mz*np.sin(alpha)*np.sin(beta))/2)**2

        M[i,:]= 1-(vec_term*angle_term)

    # Signal = M.prod(axis=0)
    # F = ((M+1)/2)
    return M
