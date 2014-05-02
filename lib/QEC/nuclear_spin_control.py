''' A module to calculate the 13C nuclear and electron spin dynamics
under dynamical decoupling gates. By THT '''

import numpy as np
import qutip
import analysis.lib.QEC.hyperfine_params as hf
from matplotlib import pyplot as plt

### import the hyperfine parameters ###
import hyperfine_params as hf_params; reload(hf)
hf = hf_params.hyperfine_params

### import the experimental values for tau and N ###
import measurement.scripts.lt2_scripts.setup.msmt_params as msmt_params; reload(msmt_params)
mp = msmt_params.cfg['samples']['Hans_sil1']


#######################
### Basic functions ###
#######################

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
    return ket0, bra0, ket1, bra1, rho0, rho1, rhom

### create a set of usefull simple states and gates
Id, sx, sy, sz = pauli()                                # Electron spin operators
Id, Ix, Iy, Iz = pauli()                                # Nuclear spin operators
X,Y,Z,x,y,z,mX,mY,mZ,mx,my,mz = basic_spin_rotations()  # Basic gates
ket0, bra0, ket1, bra1, rho0, rho1, rhom = basic_spin_states() # Basic states

###########################
### Auxilairy functions ###
###########################

def print_matrix(Qobject):
    print np.round(Qobject.full()*100)/100

def find_phase_gate(total_time, carbon_nr, axis_phase):
    '''function to determine the parameters of the preceding DD phase
    gate to set the phase, NOTE, implementation pending on the experimental implementation'''
    pass


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
    '''calculates the evolution matrices for a single
    Carbon spin, electron is always qubit 1'''

    omega_Larmor = 2 * np.pi * B_field * 1.07e3
    A_par = 2 * np.pi * hf['C' + str(carbon_nr)]['par']
    A_perp = 2 * np.pi *hf['C' + str(carbon_nr)]['perp']

    U0, U1 = nuclear_gate(number_of_pulses, tau, omega_Larmor, A_par, A_perp)

    gate = qutip.tensor(rho0,U0) + qutip.tensor(rho1,U1)
    return gate

def Ren_gate(carbon_nr, B_field=304.22, phase=0):
    '''create a Ren gate for given carbon number'''

    number_of_pulses = mp['C' + str(carbon_nr) + '_Ren_N']
    tau = mp['C' + str(carbon_nr) + '_Ren_tau']
    print number_of_pulses
    print tau
    Ren = c13_gate(carbon_nr, number_of_pulses, tau, B_field=304.22)

    return Ren

def xn_gate(carbon_nr, phase):
    pass

def c13_gate_multiqubit(carbon_nr, number_of_pulses, tau, B_field):
    '''calculates the evolution matrices a multiqubit space,
    the electron is always qubit 1'''
    pass

###################
### Experiments ###
###################

def nuclear_rabi_no_init(carbon_nr, tau, nr_of_pulses_list=np.linspace(0,100,51), B_field=304.22):
    '''nuclear Rabi experiment without init
    scheme: x - Ren(N) - x - RO'''

    #initial states
    rho_nuc_i = rhom
    rho_el_i  = rho0
    rho_init  = qutip.tensor(rho_el_i,rho_nuc_i)
    print 'initial state = '
    print_matrix(rho_init)

    #electron gates
    electron_x = qutip.tensor(x,Id)
    electron_mx = qutip.tensor(mx,Id)

    #sequence and RO
    S = np.zeros(len(nr_of_pulses_list))
    for i, N in enumerate(nr_of_pulses_list):
        gate = c13_gate(carbon_nr, N, tau, B_field)         # Define nuclear spin gate

        seq  = electron_x*gate*electron_x                   # Define gate sequence
        rho_final = seq*rho_init*seq.dag()                  # Apply gate sequence

        rho_el_final = rho_final.ptrace(0)                  # Trace out the nuclear spin
        #S[i] = qutip.expect(sz, rho_el_final) + 1./2       # Z measurement two alternative ways
        S[i] = qutip.fidelity(rho0, rho_el_final)**2


    ## plot ##
    f, ax = plt.subplots(1)
    ax.plot(nr_of_pulses_list, S, 'o-', lw=1)
    ax.set_title('P(ms=0)'); ax.set_xlabel('N')
    plt.show()
    return S[i]

def nuclear_ramsey_no_init(carbon_nr, tau_wait, N_wait_list, B_field=304.22):
    '''nuclear Rabi experiment without init
    scheme: x - Ren - DD_wait - Ren - x - RO'''

    #initial states
    rho_nuc_i = rhom
    rho_el_i  = rho0
    rho_init  = qutip.tensor(rho_el_i,rho_nuc_i)
    print 'initial state = '
    print_matrix(rho_init)

    #wait time pulses, in this case taken from the msmt_params
    N   = mp['C' + str(carbon_nr) + '_Ren_N']
    tau = mp['C' + str(carbon_nr) + '_Ren_tau']

    #gates
    electron_x  = qutip.tensor(x,Id)
    electron_mx = qutip.tensor(mx,Id)
    Ren         = c13_gate(carbon_nr, N, tau, B_field)

    #sequence and RO
    S = np.zeros(len(N_wait_list))
    for i, N_wait in enumerate(N_wait_list):

        DD_wait = c13_gate(carbon_nr, N_wait, tau_wait, B_field)         # Define DD waiting gate

        seq  = electron_mx*Ren*DD_wait*Ren*electron_x                   # Define gate sequence
        rho_final = seq*rho_init*seq.dag()                  # Apply gate sequence

        rho_el_final = rho_final.ptrace(0)                  # Trace out the nuclear spin
        #S[i] = qutip.expect(sz, rho_el_final) + 1./2       # Z measurement two alternative ways
        S[i] = qutip.fidelity(rho0, rho_el_final)**2

    ## plot ##
    f, ax = plt.subplots(1)
    ax.plot(N_wait_list*2*tau_wait*1e6, S, 'o-', lw=1)
    ax.set_title('P(ms=0)'); ax.set_xlabel('Evolution_time (us)')
    plt.show()
    return S[i]

def nuclear_init_gate(carbon_nr, init_state):
    '''function that returns a gate sequence for nuclear spin initialization
    seq = y - Ren - x - Rz - Ren '''
    pass

##########################################
### Nuclear evolution characterization ###
##########################################

def calc_operator_rotation_axis_and_angle(operator):
    '''Calculate the angle and axis of rotation of a given
    single qubit unitary operator'''

    #Get eigenstates and eigenvalues and calculate the rotation axis and angle from them
    eig_vals, eig_states = operator.eigenstates()
    angle = -1*np.angle(eig_vals[0])*2
    X_axis_projection  = qutip.expect(qutip.sigmax(),eig_states[0])*np.sign(angle)
    Y_axis_projection  = qutip.expect(qutip.sigmay(),eig_states[0])*np.sign(angle)
    Z_axis_projection  = qutip.expect(qutip.sigmaz(),eig_states[0])*np.sign(angle)

    return np.array([X_axis_projection, Y_axis_projection, Z_axis_projection]), np.abs(angle)

def characterize_c13_DD_unit(carbon_nr, B_field=304.22, tau_list = np.linspace(10,5000,500)):

    A_par = 2*np.pi*hf['C' + str(carbon_nr)]['par']
    A_perp = 2*np.pi*hf['C' + str(carbon_nr)]['perp']

    print A_par/2./np.pi
    print A_perp/2./np.pi

    characterize_DD_unit(A_par,A_perp,B_field=304.22,tau_list=tau_list)

def characterize_DD_unit(A_par = 2*np.pi*100e3, A_perp = 2*np.pi*30e3, B_field = 304.22, tau_list = np.linspace(10,5000,500), N=32):
    '''gives a full characterization of the rotation matrix
    for a single DD unit tau - X - 2tau - X - tau'''

    angle0 = np.zeros(len(tau_list))
    angle1 = np.zeros(len(tau_list))
    X_proj_0 = np.zeros(len(tau_list)); Y_proj_0 = np.zeros(len(tau_list)); Z_proj_0 = np.zeros(len(tau_list))
    X_proj_1 = np.zeros(len(tau_list)); Y_proj_1 = np.zeros(len(tau_list)); Z_proj_1 = np.zeros(len(tau_list))
    innerprod =np.zeros(len(tau_list)); pulses_for_pi2 = np.zeros(len(tau_list)); signal=np.zeros(len(tau_list))

    omega_Larmor = 2 * np.pi * B_field * 1.07e3

    for i in range(len(tau_list)):
        tau = tau_list[i]*1e-9
        if i%10 == 0:
            print str(tau_list[i]) + ' out of ' + str(max(tau_list))
        V0, V1 = nuclear_rotation_matrix(tau, omega_Larmor, A_par, A_perp)
        axes0, angle0[i] = calc_operator_rotation_axis_and_angle(V0)
        axes1, angle1[i] = calc_operator_rotation_axis_and_angle(V1)

        X_proj_0[i] = axes0[0]
        Y_proj_0[i] = axes0[1]
        Z_proj_0[i] = axes0[2]
        X_proj_1[i] = axes1[0]
        Y_proj_1[i] = axes1[1]
        Z_proj_1[i] = axes1[2]

        innerprod[i]    = X_proj_0[i]*X_proj_1[i] + Y_proj_0[i]*Y_proj_1[i] + Z_proj_0[i]*Z_proj_1[i]
        pulses_for_pi2[i]  = (np.pi/2./(np.pi-abs(np.pi-angle0[i].real))).real
        signal[i] = ( V0**(N/2) * (V1**(N/2)).dag() ).tr().real/4+1./2

    #plots
    plt.close('all')

    f, ax = plt.subplots(3,3)
    ax[0,0].plot(tau_list/1e3,X_proj_0, '-', lw=1,label = 'data')
    ax[0,0].set_title('X projection ms=0'); ax[0,0].set_xlabel('tau (us)')

    ax[1,0].plot(tau_list/1e3,X_proj_1, '-', lw=1,label = 'data')
    ax[1,0].set_title('X projection ms=1'); ax[1,0].set_xlabel('tau (us)')

    ax[2,0].plot(tau_list/1e3,(X_proj_0-X_proj_1), '-', lw=1,label = 'data')
    ax[2,0].set_title('X projection ms=0 - X projection ms=1'); ax[2,0].set_xlabel('tau (us)')

    ax[0,1].plot(tau_list/1e3,Y_proj_0, '-', lw=1,label = 'data')
    ax[0,1].set_title('Y projection ms=0'); ax[0,1].set_xlabel('tau (us)')

    ax[1,1].plot(tau_list/1e3,Y_proj_1, '-', lw=1,label = 'data')
    ax[1,1].set_title('Y projection ms=1'); ax[1,1].set_xlabel('tau (us)')

    ax[2,1].plot(tau_list/1e3,(Y_proj_0-Y_proj_1), '-', lw=1,label = 'data')
    ax[2,1].set_title('Y projection ms=0 - Y projection ms=1'); ax[2,1].set_xlabel('tau (us)')

    ax[0,2].plot(tau_list/1e3,Z_proj_0, '-', lw=1,label = 'data')
    ax[0,2].set_title('Z projection ms=0'); ax[0,2].set_xlabel('tau (us)')

    ax[1,2].plot(tau_list/1e3,Z_proj_1, '-', lw=1,label = 'data')
    ax[1,2].set_title('Z projection ms=1'); ax[1,2].set_xlabel('tau (us)')

    ax[2,2].plot(tau_list/1e3,(Z_proj_0-Z_proj_1), '-', lw=1,label = 'data')
    ax[2,2].set_title('Z projection ms=0 - Z projection ms=1'); ax[2,2].set_xlabel('tau (us)')

    f2, ax2 = plt.subplots(4,1)
    ax2[0].plot(tau_list/1e3,innerprod, '-', lw=1,label = 'data')
    ax2[0].set_title('axis innerporduct'); ax2[0].set_xlabel('tau (us)')
    ax2[1].plot(tau_list/1e3,angle0/np.pi, '-', lw=1,label = 'data')
    ax2[1].set_title('rotation_angle'); ax2[1].set_xlabel('tau (us)'); ax2[1].set_ylim(0,2)
    ax2[2].plot(tau_list/1e3,pulses_for_pi2, '-', lw=1,label = 'data')
    ax2[2].set_title('nr of pulses'); ax2[2].set_xlabel('tau (us)')
    ax2[3].plot(tau_list/1e3,signal, '-', lw=1,label = 'data')
    ax2[3].set_title('signal'); ax2[3].set_xlabel('tau (us)')

    plt.show()

print 'succes'
