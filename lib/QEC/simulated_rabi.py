#Module that contains functions to simulate rabi oscillations
import numpy as np


def Rabi_formula(Omega_R,Delta = 0, t = 0):
    '''
    Rabi formula,
    ------
    inputs
    ------
    Omega_R:              Rabi frequency in Hz
    Delta:                    Detuning in Hz
    t:                           Time in s
    --------
    returns
    --------
    Signal of Rabi oscillation at given input parameters
    '''
    Signal = 1-Omega_R**2/(Delta**2 +Omega_R**2) *np.sin(t/2*np.sqrt(Delta**2+Omega_R**2))**2
    return Signal

def hyperfine_split_rabi(Omega_R,Delta = 0, t=0, init_nucl=[0,1,0],hyperfine_split = [-2.16e6,0,2.16e6], individual_pops = False) :
    '''
    Signal of rabi oscillation when addressing hyperfine split lines
    ------
    inputs
    ------
    Omega_R:              Rabi frequency in Hz
    Delta:                    Detuning in Hz
    t:                           Time in s
    init_nucl:               list of percentages initialised in different HF-split states
    hyperfine_split:     list of splittings relative to center freq in Hz
    individual_pops:    boolean, determines if function returns signal or individual populations
    --------
    returns
    --------
    when individual pops ==False
    Signal of Rabi oscillation at given input parameters
    when individual pops ==True
    list of individual populations at given input parameters
    '''
    P = np.zeros(np.size(init))
    for ind,initial_fraction in enumerate(init):
        P[ind] = initial_fraction*Rabi_formula(Omega_R,Delta+hyperfine_split[ind], t)
    if individual_pops == False:
        signal = np.sum(P)
        return signal
    else:
        return P

def Rabi_oscillation(Omega_R,Delta,t_list):
    '''
    Rabi oscillation ,
    ------
    inputs
    ------
    Omega_R:              Rabi frequency in Hz
    Delta:                    Detuning in Hz
    t_list:                     list of times in s
    --------
    returns
    --------
    rabi_oscillation for times specified. Uses Rabi_formula function
    '''
    S = []
    for t in t_list:
        S.append(Rabi_formula(Omega_R, Delta, t))
    return S

def HF_split_rabi_osc(Omega_R,Delta,t_list, init_nucl=[0,1,0],hyperfine_split = [-2.16e6,0,2.16e6]):
    '''
    Signal of rabi oscillation when addressing hyperfine split lines
    ------
    inputs
    ------
    Omega_R:              Rabi frequency in Hz
    Delta:                    Detuning in Hz
    t:                           Time in s
    init_nucl:               list of percentages initialised in different HF-split states
    hyperfine_split:     list of splittings relative to center freq in Hz
    individual_pops:    boolean, determines if function returns signal or individual populations
    --------
    returns
    --------
    Signal of Rabi oscillation at given input parameters
    uses hyperfine_split_rabi function
    '''
    S = []
    for t in t_list:
        S.append(hyperfine_split_rabi(Omega_R,Delta,t, init_nucl,hyperfine_split))
    return S

