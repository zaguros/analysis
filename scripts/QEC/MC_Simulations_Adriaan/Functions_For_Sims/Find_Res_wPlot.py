import numpy as np
from scipy import linalg #For matrix exponential

from find_mins import find_mins
from check_pulse_t import check_pulse_t

def Find_Res(HyperfineStr ,B_field  ,tau_max,tau_step , max_gate_time):
    """
    Finds resonant driving conditions and number of pulses required for a hyperfine strength
    input HyperfineStr in real freq
    but number of pulses unreliable because of dependency on timestep
    """
    #physical constants
    gamma_c = 1.071e3 #g-factor for C13 in MHz/G
    #Model parameters
    timesteps = np.arange(tau_step,tau_max,tau_step)
    omega_larmor = 2*np.pi*gamma_c*B_field
    tau_larmor = 2*np.pi/omega_larmor

    A= 2*np.pi*HyperfineStr[0]
    B= 2*np.pi*HyperfineStr[1]  #Converts to radial frequency in Hz
    omega_tilde = np.sqrt((A+omega_larmor)**2+B**2)
    mx = B/omega_tilde
    mz = (A+omega_larmor)/omega_tilde

    ax_prod = np.zeros(np.size(timesteps))
    #Char evolution
    for idt, tau in enumerate(timesteps):
        alpha = omega_tilde*tau
        beta = omega_larmor*tau
        ax_prod[idt] = 1-mx**2 *((1-np.cos(alpha))*(1-np.cos(beta)))/(1+np.cos(alpha)*np.cos(beta)-mz*np.sin(alpha)*np.sin(beta))
    #Find Resonances as minima
    dip_ind = find_mins(ax_prod)
    theta = np.zeros(np.size(dip_ind))
    for idi, ind in enumerate(dip_ind):  #Either save V0's and calc only needed theta or calc all and not save V0
        tau = timesteps[ind]
        alpha = omega_tilde*tau
        beta = omega_larmor*tau
        theta[idi] = np.pi-np.arccos(np.cos(alpha)*np.cos(beta)-mz*np.sin(alpha)*np.sin(beta))#pi - angle because
    n_pulses = np.divide(np.pi,(2*(np.pi-np.abs(np.pi-theta)))).astype(int)
    tlist = timesteps[dip_ind]
    res = check_pulse_t(tlist, n_pulses, max_gate_time)
    return (ax_prod,res)
