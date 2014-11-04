"""
This script uses a specific function of nuclear_spin_control to
estimate the outcome of the carbon tomographz along Z for various 
electronic input states.
NK 2014-11-03
"""

import numpy as np
import os, sys
import qutip

sys.path.append("/measuring/")
from analysis.lib.tools import toolbox, plot
# from analysis.lib.tools import toolbox
# from analysis.lib.tools import plot
from analysis.lib.QEC import nuclear_spin_control as nsc
reload(nsc)

#create some useful states.
ket0, bra0, ket1, bra1, rho0, rho1, rhom, ketx,brax,ketmx,bramx,rhox,rhomx,kety,bray,ketmy,bramy,rhoy,rhomy = nsc.basic_spin_states()



nsc.nuclear_TomoZ(carbon_nrs=[1],N=48,Tau=4.994e-6,Bfield=404.69,ms='-1',input_state=qutip.tensor(rho0,rho1))
