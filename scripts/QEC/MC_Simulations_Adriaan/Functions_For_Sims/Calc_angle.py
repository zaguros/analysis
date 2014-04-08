import numpy as np
from scipy import linalg
    def Calc_angle (M):
        """
        Returns the radial angle of rotation for a rotation matrix
        NOTE only works for small angles
        """
        theta = 2*np.arccos(np.trace(M)/2)
        return theta
