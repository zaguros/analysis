"""
Module that helps reconstructing the density matrix from measurement
data and calculating fidelities.
"""

import numpy as np
import sympy
import error
from analysis.lib.math import error

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

def fidelity(statevec, dm, u_dm,ustatevec=np.array([0,0])):
    '''
    returns the fidelity (and its uncertainty) of the measured density
    matrix with a given state vector.
    '''
   
    f = error.Formula()
    beta,a,x,b,alpha = sympy.symbols('beta,a,x,b,alpha')
    '''
    f.formula = np.real(np.dot(np.array([alpha,beta]).conj(), np.dot(np.array([[a,x+1j*iy],[x-1j*iy,1.-a]]), np.array([alpha,beta]))))
    '''
    f.formula=beta*np.conj(beta)+ a*(alpha*np.conj(alpha)-beta*np.conj(beta))+x*(np.conj(alpha)*beta+alpha*np.conj(beta))+1j*(alpha*np.conj(beta)-np.conj(alpha)*beta)*b
    

    f.values[beta]=statevec[1]
    f.values[a]=float(dm[0,0])
    f.values[x]=float(np.real(dm[0,1]))
    f.values[b]=float(np.imag(dm[0,1]))

    f.uncertainties[beta]=ustatevec[1]
    f.uncertainties[a]=u_dm[0,0]
    f.uncertainties[x]=float(np.real(u_dm[0,1]))
    f.uncertainties[b]=float(np.imag(u_dm[0,1]))
    print float(np.imag(u_dm[0,1]))
    fid,ufid = f.num_eval(alpha,statevec[0],ustatevec[0])

    return (fid,ufid)

