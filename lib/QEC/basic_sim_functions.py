import numpy as np
import qutip
import analysis.lib.QEC.hyperfine_params as hf
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import analysis.lib.QEC.basic_sim_functions as bs

### import the hyperfine parameters ###
import hyperfine_params as hf_params; reload(hf_params)
hf = hf_params.hyperfine_params
### import the experimental values for tau and N ###

import measurement.scripts.lt2_scripts.setup.msmt_params as msmt_params; reload(msmt_params)
mp = msmt_params.cfg['samples']['Hans_sil1']


def pauli():
    '''Define pauli spin matrices'''
    identity = qutip.qeye(2)
    sx = qutip.sigmax()/2
    sy = qutip.sigmay()/2
    sz = qutip.sigmaz()/2
    return identity, sx, sy, sz

def basic_spin_rotations():
    ''' define some simple spin rotations'''
    X = (-1j*sx*np.pi).expm();   mX = (1j*sx*np.pi).expm()
    Y = (-1j*sy*np.pi).expm();   mY = (1j*sy*np.pi).expm()
    Z = (-1j*sz*np.pi).expm();   mZ = (1j*sz*np.pi).expm()
    x = (-1j*sx*np.pi/2).expm(); mx = (1j*sx*np.pi/2).expm()
    y = (-1j*sy*np.pi/2).expm(); my = (1j*sy*np.pi/2).expm()
    z = (-1j*sz*np.pi/2).expm(); mz = (1j*sz*np.pi/2).expm()
    return X,Y,Z,x,y,z,mX,mY,mZ,mx,my,mz

def basic_spin_states():
    '''define some basic spin states'''
    ket0 = qutip.basis(2,0)
    bra0 = qutip.basis(2,0).dag()
    ket1 = qutip.basis(2,1)
    bra1 = qutip.basis(2,1).dag()
    rho0 = ket0*bra0
    rho1 = ket1*bra1
    rhom = (rho0+rho1)/2
    ketx = 1/np.sqrt(2)*(qutip.basis(2,0)+qutip.basis(2,1))
    brax = 1/np.sqrt(2)*(qutip.basis(2,0).dag()+qutip.basis(2,1).dag())
    ketmx = 1/np.sqrt(2)*(qutip.basis(2,0)-qutip.basis(2,1))
    bramx = 1/np.sqrt(2)*(qutip.basis(2,0).dag()-qutip.basis(2,1).dag())
    kety = 1/np.sqrt(2)*(qutip.basis(2,0)+1j*qutip.basis(2,1))
    bray = kety.dag()
    ketmy = 1/np.sqrt(2)*(qutip.basis(2,0)-1j*qutip.basis(2,1))
    bramy = ketmy.dag()
    rhox =ketx*brax
    rhomx = ketmx*bramx
    rhoy =kety*bray
    rhomy = ketmy*bramy
    return ket0, bra0, ket1, bra1, rho0, rho1, rhom, ketx,brax,ketmx,bramx,rhox,rhomx,kety,bray,ketmy,bramy,rhoy,rhomy

### create a set of usefull simple states and gates
Id, sx, sy, sz = pauli()                                # Electron spin operators
Id, Ix, Iy, Iz = pauli()                                # Nuclear spin operators
X,Y,Z,x,y,z,mX,mY,mZ,mx,my,mz = basic_spin_rotations()  # Basic gates
ket0, bra0, ket1, bra1, rho0, rho1, rhom, ketx,brax,ketmx,bramx,rhox,rhomx,kety,bray,ketmy,bramy,rhoy,rhomy = basic_spin_states() # Basic states

def any_pure_state(alpha,beta,return_psi = False,return_rho = True):
    '''gives out your psi and if wanted your rho for a state alpha 0 + beta 1 '''
    psi = alpha*ket0+beta*ket1
    rho = psi*psi.dag()
    # print psi
    # print rho
    if return_psi == True:
        return psi
    if return_rho == True:
        return rho

def any_mixed_state(alpha,beta):
    '''gives out a mixture of rho0 and rho1'''
    rho = alpha *rho0+beta*rho1
    return rho

###########################
### Auxilairy functions ###
###########################

def print_matrix(Qobject):
    print np.round(Qobject.full()*100)/100
    print type(np.round(Qobject.full()*100)/100)

def get_C13_hyperfine_params(carbon_nrs, ms = '+1'):
    '''
    load hyperfine paramters for a given list of carbon_nrs
    ms = '+1' or '-1' indicates which electron transition is used
    (we alter the  sign of the parallel component of the hypefine interaction
     for the ms=-1 transition)
    '''
    A_par   = []
    A_perp  = []

    for kk, carbon_nr in enumerate(carbon_nrs):
        perp    =  2*np.pi*hf['C' + str(carbon_nr)]['perp']
        if ms == '+1':
            par     = 2*np.pi*hf['C' + str(carbon_nr)]['par']
        elif ms == '-1':
            par     = -2*np.pi*hf['C' + str(carbon_nr)]['par']

        A_par.append(par)
        A_perp.append(perp)

    return A_par, A_perp

###################################
### Nuclear evolution and gates ###
###################################

def nuclear_rotation_matrix(tau, omega_Larmor, A_par, A_perp):
    ''' Function to calculate a C13 evolution matrix'''

    #Hamiltonian for ms=0 and ms=+/-1
    H0 = omega_Larmor * Iz
    H1 = (A_par+omega_Larmor)*Iz + A_perp*Ix
    #Evolution during tau for ms=0 and ms=+/-1
    expH0 = (-1j*H0*tau).expm();    expH1 = (-1j*H1*tau).expm()
    #Evolution during a decouple unit
    V0 = expH0*expH1*expH1*expH0;   V1 = expH1*expH0*expH0*expH1

    return V0, V1

def nuclear_gate(number_of_pulses, tau, omega_Larmor, A_par, A_perp):
    '''Gives the evolution matrix for number_of_pulses pulses'''

    V0, V1 = nuclear_rotation_matrix(tau, omega_Larmor, A_par, A_perp)

    U0 = V0 ** (number_of_pulses/2)
    U1 = V1 ** (number_of_pulses/2)

    return U0, U1

def c13_gate(carbon_nr, number_of_pulses, tau, B_field=304.22):
    '''calculates the evolution matrices for a single Carbon spin
    For an Ren gate the ideal gate is rounded tot 1/sqrt(2)
    For a phase gate (set phase from None to the required phase) the ideal gate is just calculated using the rotation matrix
    When phase_y is true, the phase gets the opposite sign, this is to simulate what we really do in the experiment - should be improved (JULIA)
    Note: this is now standard on 'True' as this is the only phase gate we used up to now!
    '''    

    #Hamiltonian for ms=0 and ms=+/-1
    omega_Larmor = 2 * np.pi * B_field * 1.07e3
    A_par = 2 * np.pi * hf['C' + str(carbon_nr)]['par']
    A_perp = 2 * np.pi *hf['C' + str(carbon_nr)]['perp']

    U0, U1 = nuclear_gate(number_of_pulses, tau, omega_Larmor, A_par, A_perp)

    ## if REN gate, create nice gate, else, round at least a bit
    if number_of_pulses == mp['C' + str(carbon_nr) + '_Ren_N'][0] and tau == mp['C' + str(carbon_nr) + '_Ren_tau'][0]:
        U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
        U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    else:
        U0id = np.round(U0.full()*100)/100.
        U1id = np.round(U1.full()*100)/100.

    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)

    return U0, U1, U0id, U1id

def nuclear_Ren_matrix(carbon_nr,B_field=304.22):
    ''' Gives out the non-ideal and ideal Ren matrix elements for a single Carbon spin '''

    number_of_pulses = mp['C' + str(carbon_nr) + '_Ren_N'][0]
    tau = mp['C' + str(carbon_nr) + '_Ren_tau'][0]

    U0, U1, U0id, U1id = c13_gate(carbon_nr,number_of_pulses, tau, B_field)

    return U0, U1, U0id, U1id

def Ren_gate(carbon_nr, B_field=304.22):
    '''create a Ren gate for given carbon number, only interacting with the electron spin '''

    U0, U1, U0id, U1id = nuclear_Ren_matrix(carbon_nr, B_field)

    Ren = qutip.tensor(rho0,U0)+qutip.tensor(rho1,U1)
    Ren_id = qutip.tensor(rho0,U0id)+qutip.tensor(rho1,U1id)

    return Ren, Ren_id


def waittime(carbon_nr, B_field=304.22,ms = 'dec'):
    '''
    calculates the matrix elements for a carbon spin during a waittime
    ms is the state of the electron spin during this time, can be: 
    'dec': electron is decoupled
    'ms0': electron in zero
    'ms1': electron in +1
    'msm1': electron in -1
    IN PROGRESS
    '''

    omega_Larmor = 2 * np.pi * B_field * 1.07e3
    A_par = 2 * np.pi * hf['C' + str(carbon_nr)]['par']
    A_perp = 2 * np.pi *hf['C' + str(carbon_nr)]['perp']

    H0 = omega_Larmor * Iz
    H1 = (A_par+omega_Larmor)*Iz + A_perp*Ix
    
    precession_freq = mp['C' + str(carbon_nr) +'_freq']*2*np.pi

    current_phase = total_time*precession_freq
    phase_dif = (phase-current_phase)%(2*np.pi)
    dec_time =  (phase_dif)/precession_freq

    tau = dec_time/4

    U0, U1 = nuclear_gate(2, tau, omega_Larmor, A_par, A_perp)

    if return_gate == True:
        return U0, U1
    if return_tau == True:
        return tau


