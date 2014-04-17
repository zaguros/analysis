import numpy as np
from scipy import linalg
from Calc_axis import Calc_axis

def Char_Gate(NV,res ,B_field):
    """
    Characterize the gate, take the NV centre, the resonance paramters and the Bfield as input
    returns the fidelity with which an x-gate can be implemented.
    """
    #physical constants
    gamma_c = 1.071e3 #g-factor for C13 in Hz/G
    #Model parameters
    omega_larmor = 2*np.pi*gamma_c*B_field
    tau_larmor = 2*np.pi/omega_larmor
    tau = res[0]
    n_pulses = int(res[1]*2) #So that we do a pi -pulse

    M = np.zeros(np.shape(NV)[0])
    for idC in range(np.shape(NV)[0]):
        A= 2*np.pi*NV[idC,0]
        B= 2*np.pi*NV[idC,1]  #Converts to radial frequency in Hz/G
        omega_tilde = np.sqrt((A+omega_larmor)**2+B**2)
        alpha = omega_tilde*tau
        beta = omega_larmor*tau
        mx = B/omega_tilde
        mz = (A+omega_larmor)/omega_tilde
        vec_term = mx**2 *((1-np.cos(alpha))*(1-np.cos(beta)))/(1+np.cos(alpha)*np.cos(beta)-mz*np.sin(alpha)*np.sin(beta))
        angle_term = np.sin(n_pulses*np.arccos(np.cos(alpha)*np.cos(beta)-mz*np.sin(alpha)*np.sin(beta))/2)**2
        M[idC] = 1-(vec_term*angle_term)

    Signal = M.prod()
    F = (1-(Signal+1)/2)
    return F
