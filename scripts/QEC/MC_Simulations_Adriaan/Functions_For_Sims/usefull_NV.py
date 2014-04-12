import numpy as np
def usefull_NV(number_carbons, Gate_min = 5):
    """
    Judges if a NV has enough addressable gates
    takes a 2x2 array as input and returns a 1D array containing number of usefull_NV's
    """
    u=number_carbons>Gate_min
    n_use= (u.sum(0).astype(float))
    return n_use
