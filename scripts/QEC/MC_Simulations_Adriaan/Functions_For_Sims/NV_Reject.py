import numpy as np
def NV_Reject(NV,A_max):
    """
    Compare the size of the hyperfine interactions of an NV array with A_max and return False if it is greater than A_max
    Input is array of shape (:,2) output is boolean
    """
    A = np.sqrt(np.add(np.square(NV[:,0]),np.square(NV[:,1])))
    return np.max(A)>A_max

