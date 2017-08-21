import numpy as np
import qutip
import analysis.lib.QEC.hyperfine_params as hf
from matplotlib import pyplot as plt
import matplotlib.cm as cm

### import the hyperfine parameters ###
import hyperfine_params as hf_params; reload(hf_params)
hf = hf_params.hyperfine_params

def pauli(spin=0.5):
    '''Define pauli spin matrices'''
    
        ## Spin 1/2
    if spin == 0.5:
        identity = qutip.qeye(2)
        sx = qutip.sigmax()/2
        sy = qutip.sigmay()/2
        sz = qutip.sigmaz()/2

        ## Spin 1
    if spin == 1:
        identity = qutip.qeye(3)
        sx = qutip.jmat(1,'x')
        sy = qutip.jmat(1,'y')
        sz = qutip.jmat(1,'z')

    return identity, sx, sy, sz

### Create
Id, sx, sy, sz      = pauli(0.5)                              # Electron spin 1/2 operators
Id1, sx1, sy1, sz1  = pauli(1)                                # Electron spin 1 operators
Id, Ix, Iy, Iz      = pauli(0.5)                              # Nuclear spin 1/2 operators, C13, N15
Id1, Ix1, Iy1, Iz1  = pauli(1)                                # Nuclear spin 1/2 operators, N14

def basic_spin_rotations():
    ''' define some simple spin rotations'''
    
        ## Spin 1/2
    X = (-1j*sx*np.pi).expm();   mX = (1j*sx*np.pi).expm()
    Y = (-1j*sy*np.pi).expm();   mY = (1j*sy*np.pi).expm()
    Z = (-1j*sz*np.pi).expm();   mZ = (1j*sz*np.pi).expm()
    x = (-1j*sx*np.pi/2).expm(); mx = (1j*sx*np.pi/2).expm()
    y = (-1j*sy*np.pi/2).expm(); my = (1j*sy*np.pi/2).expm()
    z = (-1j*sz*np.pi/2).expm(); mz = (1j*sz*np.pi/2).expm()

    return X,Y,Z,x,y,z,mX,mY,mZ,mx,my,mz

def basic_spin_rotations_spin1():
    ''' define some simple spin rotations'''
    
        ### m1 transition pulses
    ssx =   qutip.qutrit_ops()[0]+qutip.qutrit_ops()[4]/2+qutip.qutrit_ops()[4].dag()/2
    ssy =   qutip.qutrit_ops()[0]-qutip.qutrit_ops()[4]/2*1j+qutip.qutrit_ops()[4].dag()/2*1j
    ssz =   qutip.qutrit_ops()[0]+qutip.qutrit_ops()[1]/2-qutip.qutrit_ops()[2].dag()/2

        ## Spin 1/2
    X = (-1j*ssx*np.pi).expm();   mX = (1j*ssx*np.pi).expm()
    Y = (-1j*ssy*np.pi).expm();   mY = (1j*ssy*np.pi).expm()
    Z = (-1j*ssz*np.pi).expm();   mZ = (1j*ssz*np.pi).expm()
    x = (-1j*ssx*np.pi/2).expm(); mx = (1j*ssx*np.pi/2).expm()
    y = (-1j*ssy*np.pi/2).expm(); my = (1j*ssy*np.pi/2).expm()
    z = (-1j*ssz*np.pi/2).expm(); mz = (1j*ssz*np.pi/2).expm()

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

X,Y,Z,x,y,z,mX,mY,mZ,mx,my,mz             = basic_spin_rotations()       # Basic gates
X1,Y1,Z1,x1,y1,z1,mX1,mY1,mZ1,mx1,my1,mz1 = basic_spin_rotations_spin1() # Basic gates, spin 1
ket0, bra0, ket1, bra1, rho0, rho1, rhom, ketx,brax,ketmx,bramx,rhox,rhomx,kety,bray,ketmy,bramy,rhoy,rhomy = basic_spin_states() # Basic states

###########################
### Auxilairy functions ###
###########################

def print_matrix(Qobject):
    print np.round(Qobject.full()*1e9)/1e9
    print np.round(Qobject.full()*100)/100


###################################
### Nuclear evolution and gates ###
###################################

def Hamiltonian():

    ### Magnetic field
    Bz = 400        #Gauss
    Bx   = 0
    By   = 0
    
    ### Constants
    ZFS             = 2.87e9        #Zero field splitting, Hz
    Quad            = 5e6           #Quadrupole splitting, Hz
    A_zz            = 2.18e6        #Hyperfine interaction, Hz
    A_xx            = 0*20e6        #Transversal hyperfine interaction
    A_yy            = A_xx         

    ye = 2.8e6      #Hz/Gauss
    yn = 0.3e3      #Hz/Gauss

    ### Electron    
    H_e   =    ZFS*sz1**2 + Bz*ye*sz1 + Bx*ye*sx1 + By*ye*sy1 

    ### Electron + Nitrogen
    H_en  = (ZFS*qutip.tensor(sz1**2,Id1)   + 
            Bz*ye*qutip.tensor(sz1,Id1)     + 
            Bx*ye*qutip.tensor(sx1,Id1)     + 
            By*ye*qutip.tensor(sy1,Id1)     +
            Quad*qutip.tensor(Id1,Iz1**2)   +
            Bz*yn*qutip.tensor(Id1,Iz1)     +
            Bx*yn*qutip.tensor(Id1,Ix1)     +
            By*yn*qutip.tensor(Id1,Iy1)     +
            A_zz*qutip.tensor(sz1,Iz1)      +
            A_xx*qutip.tensor(sx1,Ix1)      +
            A_yy*qutip.tensor(sy1,Iy1)      )
    ### C13 spin 

    return H_en

def Evolution(H, tau_list, N):
    ''' H is the hamiltonian'''

    ### Initial states
    psi_electron =    (qutip.basis(3,1) + qutip.basis(3,2))*1./np.sqrt(2)
    rho_electron =    psi_electron*psi_electron.dag()

    # print 'electron state is'
    # print rho_electron

    rho_14N     =  qutip.basis(3,1)*qutip.basis(3,1).dag()

    # print '14N state is'
    # print rho_14N

    rho_tot     =   qutip.tensor(rho_electron,rho_14N)

    # print 'tot state is'
    # print rho_tot

    ### Pi-pulses
    pi_pulse_electron = qutip.qutrit_ops()[0]+qutip.qutrit_ops()[4]+qutip.qutrit_ops()[4].dag()
    pi_pulse = qutip.tensor(pi_pulse_electron,Id1)

    # print 'pi pulse is'
    # print pi_pulse

    
    #sequence and RO
    S = np.zeros(len(tau_list))
    for i, tau in enumerate(tau_list):

        ### Free evolution 
        expH = (-1j*H*tau).expm()
   
        U    = expH*pi_pulse*expH*expH*pi_pulse*expH
        U_tot   = U**(N/2)
        
        rho_final = U_tot*rho_tot*U_tot.dag() 

        rho_el_final = rho_final.ptrace(0)                  # Trace out the nuclear spin
        S[i] = qutip.fidelity(rho_electron, rho_el_final)**2


    ## plot ##
    f, ax = plt.subplots(1)
    ax.plot(tau_list, S, 'o-', lw=1)
    ax.set_title('P(ms=0)'); ax.set_xlabel('tau (us)')
    plt.show()
    return S[i]




        # seq3 = electron_x*Ren*electron_y
        # rho_final = seq3*rho_seq2a*seq3.dag()


        # rho_el_final = rho_final.ptrace(0)                  # Trace out the nuclear spin
        # #S[i] = qutip.expect(sz, rho_el_final) + 1./2       # Z measurement two alternative ways
        # S[i] = qutip.fidelity(rho0, rho_el_final)**2



### How do we handle the fact that we do hard pi pulses, while we have new eigenstates? 
### Step 1, calculate new eigenstates (electron, nuclear mixtures)
### Step 2, we do still initialize and measure the "ms=0" eiegenstates (?)
### Step 3, pulses are still applie between the ms=


### How does this work for the C13 spin? 
# Electron remains diagonal because hyperfine Sx is small compared to Sz. 
# For the C13 spin this is not true, because Ix cannot be neglected. Sz*Ix needs be kept. But Sx*Iz, etc, not.
# 










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


