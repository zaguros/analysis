import numpy as np
def addressable_C(F,F_Min=0.95):
    """
    Returns the number of usable gates per NV
    input is an 1D array containing fidelities and a threshold. output is an integer
    """
    n = (np.array(F)>F_Min).sum()
    return n

