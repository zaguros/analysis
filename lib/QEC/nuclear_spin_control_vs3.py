''' A module to calculate the 13C nuclear and electron spin dynamics
under dynamical decoupling gates. By THT 2016'''

import numpy as np
from scipy import linalg

# import analysis.lib.QEC.hyperfine_params as hf
from matplotlib import pyplot as plt
import matplotlib.cm as cm

# ### import the hyperfine parameters ###
# import hyperfine_params as hf_params; reload(hf_params)
# hf = hf_params.hyperfine_params

# ### import the experimental values for tau and N ###
# # import measurement.scripts.lt2_scripts.setup.msmt_params as msmt_params; reload(msmt_params)

# ### import the theoretically tuned values for tau and N ###
# import gate_params as gate_params; reload(gate_params)
# mp = gate_params.gp



#### Basic structure ##### 
### Class: Gate
### Class: State


### START new file ###


### Pauli matrixes
identity =  np.array([[1, 0],  [ 0, 1]])
sx       = np.array([[0, 1],  [ 1, 0]])
sy       = np.array([[0, -1j],[1j, 0]])
sz       = np.array([[1, 0],  [0, -1]]) 

### Basic spin rotations

X = linalg.expm(-1j*sx*np.pi) 
X = linalg.expm(sx) 


# X = (-1j*sx*np.pi).expm();   mX = (1j*sx*np.pi).expm()
# Y = (-1j*sy*np.pi).expm();   mY = (1j*sy*np.pi).expm()
# Z = (-1j*sz*np.pi).expm();   mZ = (1j*sz*np.pi).expm()
# x = (-1j*sx*np.pi/2).expm(); mx = (1j*sx*np.pi/2).expm()
# y = (-1j*sy*np.pi/2).expm(); my = (1j*sy*np.pi/2).expm()
# z = (-1j*sz*np.pi/2).expm(); mz = (1j*sz*np.pi/2).expm()
# H = (-1j*(sx+sz)/np.sqrt(2)*np.pi).expm()


print 'YES!'









