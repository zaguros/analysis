"""
Module that helps reconstructing the density matrix from measurement
data and calculating fidelities.
"""

import numpy as np
import sympy
import error

def dm(statevec):
    '''
    returns the density matrix for the input state vector.
    '''
    return np.outer(statevec, statevec.conj())

def measured_single_qubit_dm(mx, my, mz, u_mx=0., u_my=0., u_mz=0.):
    '''
    returns the density matrix for measurement results on a single qubit
    in the x,y,z bases. we take the measurement result to be the probability
    to measure the qubit state |0>.
    returns the density matrix and the uncertainty thereof.
    '''
    a = mz; u_a = u_mz
    x = mx - .5; u_x = u_mx
    y = .5 - my; u_y = u_my

    return (np.array([[a, x+1j*y],[x-1j*y, 1.-a]]), 
            np.array([[u_a, u_x+1j*u_y],[u_x-1j*u_y, u_a]]))#1-u_a

def fidelity(statevec, dm, u_dm):
    '''
    returns the fidelity (and its uncertainty) of the measured density
    matrix with a given state vector.
    '''
    return (np.dot(statevec.conj(), np.dot(dm, statevec)).real, 
            np.dot(statevec.conj(), np.dot(u_dm, statevec)).real)



