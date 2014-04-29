''' A module to calculate the 13C nuclear and electron spin dynamics
under dynamical decoupling gates. By THT '''

import numpy as np
import qutip
import analysis.lib.QEC.hyperfine_params as hf
from matplotlib import pyplot as plt

def pauli():
    '''Define pauli spin matrices'''
    identity = qutip.qeye(2)
    sx = qutip.sigmax()/2
    sy = qutip.sigmay()/2
    sz = qutip.sigmaz()/2
    print 'test'
    return identity, sx, sy, sz

def spin_rotations():
    ''' define some simple spin rotations'''
    X = (-1j*sx*np.pi).expm()
    Y = (-1j*sy*np.pi).expm()
    Z = (-1j*sz*np.pi).expm()
    x = (-1j*sx*np.pi/2).expm()
    y = (-1j*sy*np.pi/2).expm()
    z = (-1j*sz*np.pi/2).expm()
    return X,Y,Z,x,y,z

Id, sx, sy, sz = pauli()
Id, Ix, Iy, Iz = pauli()
X,Y,Z,x,y,z = spin_rotations()


def c13_rotation_matrix(tau, omega_Larmor, A_par, A_perp):
    ''' Function to calculate a C13 evolution matrix'''

    #Hamiltonian for ms=0 and ms=+/-1
    H0 = omega_Larmor * Iz
    H1 = (A_par+omega_Larmor)*Iz + A_perp*Ix

    #Evolution during tau for ms=0 and ms=+/-1
    expH0 = (-1j*H0*tau).expm();    expH1 = (-1j*H1*tau).expm()
    #Evolution during a decouple unit
    V0 = expH0*expH1*expH1*expH0;   V1 = expH1*expH0*expH0*expH1
    #Evolution during number_of_pulses/2 units
    #U0 = V0**(number_of_pulses/2);  U1 = V1**(number_of_pulses/2)

    return V0, V1#, U0, U1

def c13_gate(carbon_number, number_of_pulses, tau):
    print 'test'

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

def calc_operators_inner_product(operator1, operator2):
    '''calculate inner product of rotation axis'''
    pass

def characterize_DD_unit(carbon_nr = 1, B_field = 300, tau_list = np.linspace(0,5000,501)):
    '''gives a full characterization of the rotation matrix
    for a single DD unit tau - X - 2tau - X - tau'''

    angle0 = np.zeros(len(tau_list))
    angle1 = np.zeros(len(tau_list))
    X_proj_0 = np.zeros(len(tau_list)); Y_proj_0 = np.zeros(len(tau_list)); Z_proj_0 = np.zeros(len(tau_list))
    X_proj_1 = np.zeros(len(tau_list)); Y_proj_1 = np.zeros(len(tau_list)); Z_proj_1 = np.zeros(len(tau_list))
    innerprod =np.zeros(len(tau_list))

    A_par = 2*np.pi*100e3
    A_perp = 2*np.pi*30e3
    omega_Larmor = 2*np.pi*300 * 1.07e3

    for i in range(len(tau_list)):
        tau = tau_list[i]*1e-9
        if i%10 == 0:
            print str(tau_list[i]) + ' out of ' + str(max(tau_list))
        V0, V1 = c13_rotation_matrix(tau, omega_Larmor, A_par, A_perp)
        axes0, angle0[i] = calc_operator_rotation_axis_and_angle(V0)
        axes1, angle1[i] = calc_operator_rotation_axis_and_angle(V1)

        X_proj_0[i] = axes0[0]
        Y_proj_0[i] = axes0[1]
        Z_proj_0[i] = axes0[2]
        X_proj_1[i] = axes1[0]
        Y_proj_1[i] = axes1[1]
        Z_proj_1[i] = axes1[2]

        innerprod[i] = X_proj_0[i]*X_proj_1[i] + Y_proj_0[i]*Y_proj_1[i] + Z_proj_0[i]*Z_proj_1[i]

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

    f2, ax2 = plt.subplots(2,1)
    ax2[0].plot(tau_list/1e3,innerprod, '-', lw=1,label = 'data')
    ax2[0].set_title('axis innerporduct'); ax2[0].set_xlabel('tau (us)')
    ax2[1].plot(tau_list/1e3,angle0/np.pi, '-', lw=1,label = 'data')
    ax2[1].set_title('rotation_angle'); ax2[1].set_xlabel('tau (us)'); ax2[1].set_ylim(0,2)

    plt.show()


print 'succes'
