#Module that contains functions to simulate nuclear spin responses to decoupling sequences
import numpy as np

def fingerprint():
    pass

def dyn_dec_signal(HFs_par,HFs_orth,B_field,N,tau):
    '''
    Takes the HF interaction strengths (paralel and orthogonal), the magnetic field strenght and an array of times and returns the signal at those times for that specific spin.
    ------
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
    #physical constants
    gamma_c = 1.071e3 #g-factor for C13 in Hz/G
    #Model parameters
    omega_larmor = 2*np.pi*gamma_c*B_field #radial frequency
    tau_larmor = 2*np.pi/omega_larmor #time in seconds


    print 'tau larmor = %s' %tau_larmor
    M=np.zeros([np.size(HFs_par),np.size(tau)])
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
