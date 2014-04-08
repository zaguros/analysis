import numpy as np
from scipy import linalg
def Calc_axis(M):
    """
    Calculates the cartesian projection of the rotation axis that the input matrix describes
    Returns a numpy array. Returns nan for the identity matrix
    """
    eff_exp = 1j*linalg.logm(M)
    x =np.real(0.5*(eff_exp[0,1]+eff_exp[1,0]))
    y =np.imag(0.5*(eff_exp[1,0]-eff_exp[0,1]))
    z =np.real(0.5*(eff_exp[0,0]-eff_exp[1,1]))
    n = np.array([x,y,z])
    n = n/linalg.norm(n)
    return n
