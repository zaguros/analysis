import numpy as np
def find_mins(vector):
    """
    Simple dip finding algorithm, checks if value is smaller than the value before and after
    Returns an np.array of indices
    """
    vect_L = np.r_[0,vector[:-1]]
    vect_R = np.r_[vector[1:],0]
    dips = np.logical_and(vector < vect_L,vector< vect_R)
    ind = np.where(dips)[0]
    return ind
