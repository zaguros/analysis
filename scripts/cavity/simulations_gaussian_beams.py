import sympy

from sympy.solvers import solve
from sympy import Symbol
import numpy as np

import math

### below functions do not work well when entering floats. Using table generated from mathemtaica instead.
# def solve_hybrid_cavity(La=4,d=4,R=39,lambda_i=637.e-9,n_diamond=2.41):
#     """
#     Solves for gaussian beam parameters in a hybrid cavity. For more details see mathematica file in 'projects/cavities/simulations/gaussian_beam_optics'
#     """
#     s_La = Symbol('s_La',real=True,positive=True)
#     s_d = Symbol('s_d',real=True,positive=True)
#     s_R = Symbol('s_R',real=True,positive=True)
#     s_z0a = Symbol('s_z0a',real=True,positive=True)
#     s_z0d = Symbol('s_z0d',real=True,positive=True)
#     s_delta_za = Symbol('s_delta_za',real=True,positive=True)
#     s_n_diamond = Symbol('s_n_diamond',real=True,positive=True)
#     print 'La',La,'d',d,'R',R
#     s= solve( [ (s_La + s_d- s_delta_za)*(1 + (s_z0a/(s_La + s_d- s_delta_za))**2 ) - s_R, 
#         (s_d- s_delta_za)*(1 + (s_z0a/(s_d- s_delta_za))**2) -  s_d*(1 + (s_z0d/s_d)**2), 
#         (s_z0a)*(1 + ((s_d- s_delta_za)/s_z0a)**2 ) - (s_z0d/s_n_diamond)*(1 + (s_d/s_z0d)**2) ], 
#         [s_z0a, s_z0d, s_delta_za] )[0]
#     print s
#     z0a,z0d,delta_za = s.subs({s_La:La,s_d:d,s_R:R,s_n_diamond:n_diamond})
#     print z0a,z0d,delta_za
#     w0 = math.sqrt(lambda_i*z0d*1.e-6/(n_diamond*math.pi))
#     return z0a*1.e-6,z0d*1.e-6,delta_za*1.e-6, w0


# def solve_hybrid_cavity2(s_La=4,s_d=4,s_R=39,lambda_i=637.e-9,s_n_diamond=2.41):
#     """
#     Solves for gaussian beam parameters in a hybrid cavity. For more details see mathematica file in 'projects/cavities/simulations/gaussian_beam_optics'
#     """

#     s_z0a = Symbol('s_z0a',real=True,positive=True)
#     s_z0d = Symbol('s_z0d',real=True,positive=True)
#     s_delta_za = Symbol('s_delta_za',real=True,positive=True)
#     print 'La',s_La,'d',s_d,'R',s_R
#     s= solve( ( (s_La + s_d- s_delta_za)*(1 + (s_z0a/(s_La + s_d- s_delta_za))**2 ) - s_R, 
#         (s_d- s_delta_za)*(1 + (s_z0a/(s_d- s_delta_za))**2) -  s_d*(1 + (s_z0d/s_d)**2), 
#         (s_z0a)*(1 + ((s_d- s_delta_za)/s_z0a)**2 ) - (s_z0d/s_n_diamond)*(1 + (s_d/s_z0d)**2) ), 
#         s_z0a, s_z0d, s_delta_za, rational=False )
#     print s
#     z0a,z0d,delta_za = s[0]#.subs({s_La:La,s_d:d,s_R:R,s_n_diamond:n_diamond})
#     print z0a,z0d,delta_za
#     w0 = math.sqrt(lambda_i*z0d*1.e-6/(n_diamond*math.pi))
#     return z0a*1.e-6,z0d*1.e-6,delta_za*1.e-6, w0

def get_w0_from_table(La=4.e-6,d=4.e-6,R=37.e-6):
    """
    Get w0 from the tables created using mathemtaica script "gaussian_beams_in_hyrbid_cavity"
    this script assumes lamdba_i = 637.e-9, and n_diamond = 2.4.
    current stepsize in diamond thickness is 0.5 um; stepsize in air is 0.2um. 
    precision error from approximating La, d and R from table is ~0.01 um in waist.
    """
    R_int = int(R*1.e6) #determine which table to use
    data = np.loadtxt(r'K:/ns/qt/Diamond/Projects/Cavities/simulations/gaussian_beam_optics/w0s/w0table_R_%s.dat'%str(R_int))
    ds = data[0,1:]
    Las = data[1:,0]
    w0s = data[1:,1:] 
    La_idx = (np.abs(Las-La*1.e6)).argmin()
    d_idx = (np.abs(ds-d*1.e6)).argmin()
    # print max(ds),max(Las)
    return w0s[La_idx,d_idx]*1.e-6