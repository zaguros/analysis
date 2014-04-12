import numpy as np
def Find_Carb(NV,A_Min):
    """
    Compares Hyperfine strengths to a threshold and returns indices if they are larger than the threshold
    Inputs a np array of shape [:,2] and a float. Outputs a 1d array of ints
    """
    A = np.sqrt(np.add(np.square(NV[:,0]),np.square(NV[:,1])))
    ind = np.where(A>A_Min)
    return ind[0]
