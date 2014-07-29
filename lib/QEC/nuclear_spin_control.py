''' A module to calculate the 13C nuclear and electron spin dynamics
under dynamical decoupling gates. By THT '''

import numpy as np
import qutip
import analysis.lib.QEC.hyperfine_params as hf
from matplotlib import pyplot as plt

### import the hyperfine parameters ###
import hyperfine_params as hf_params; reload(hf_params)
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
    ketx = 1/np.sqrt(2)*(qutip.basis(2,0)+qutip.basis(2,1))
    brax = 1/np.sqrt(2)*(qutip.basis(2,0).dag()+qutip.basis(2,1).dag())
    ketmx = 1/np.sqrt(2)*(qutip.basis(2,0)-qutip.basis(2,1))
    bramx = 1/np.sqrt(2)*(qutip.basis(2,0).dag()-qutip.basis(2,1).dag())
    rhox =ketx*brax
    rhomx = ketmx*bramx
    return ket0, bra0, ket1, bra1, rho0, rho1, rhom, ketx,brax,ketmx,bramx,rhox,rhomx

### create a set of usefull simple states and gates
Id, sx, sy, sz = pauli()                                # Electron spin operators
Id, Ix, Iy, Iz = pauli()                                # Nuclear spin operators
X,Y,Z,x,y,z,mX,mY,mZ,mx,my,mz = basic_spin_rotations()  # Basic gates
ket0, bra0, ket1, bra1, rho0, rho1, rhom, ketx,brax,ketmx,bramx,rhox,rhomx = basic_spin_states() # Basic states

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


def nuclear_Ren_matrix(carbon_nr,B_field=304.22):
    ''' difference to Ren_gate is that this gives two matrices, can combine'''

    #Hamiltonian for ms=0 and ms=+/-1
    omega_Larmor = 2 * np.pi * B_field * 1.07e3
    A_par = 2 * np.pi * hf['C' + str(carbon_nr)]['par']
    A_perp = 2 * np.pi *hf['C' + str(carbon_nr)]['perp']
    number_of_pulses = mp['C' + str(carbon_nr) + '_Ren_N'][0]
    tau = mp['C' + str(carbon_nr) + '_Ren_tau'][0]

    U0, U1 = nuclear_gate(number_of_pulses, tau, omega_Larmor, A_par, A_perp)

    return U0, U1


def nuclear_gate(number_of_pulses, tau, omega_Larmor, A_par, A_perp):
    '''Gives the evolution matrix for number_of_pulses pulses'''

    V0, V1 = nuclear_rotation_matrix(tau, omega_Larmor, A_par, A_perp)

    U0 = V0 ** (number_of_pulses/2)
    U1 = V1 ** (number_of_pulses/2)

    return U0, U1

def c13_gate(carbon_nr, number_of_pulses, tau, B_field=304.22, return_indiv = False, return_id = False, phase = None):
    '''calculates the evolution matrices for a single
    Carbon spin, electron is always qubit 1'''

    omega_Larmor = 2 * np.pi * B_field * 1.07e3
    A_par = 2 * np.pi * hf['C' + str(carbon_nr)]['par']
    A_perp = 2 * np.pi *hf['C' + str(carbon_nr)]['perp']

    U0, U1 = nuclear_gate(number_of_pulses, tau, omega_Larmor, A_par, A_perp)

    gate = qutip.tensor(rho0,U0) + qutip.tensor(rho1,U1)
    if return_indiv == False:
        return gate
    elif return_indiv == True and return_id == False:
        return U0, U1
    elif return_indiv == True and return_id == True:
        U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
        U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
        if phase != None:
            # print 'phase'+str(phase/np.pi)
            U0id = (-1j*sz*phase).expm()
            U1id = (-1j*sz*phase).expm()
            # print U0id
        U0id = qutip.Qobj(U0id)
        U1id = qutip.Qobj(U1id)
        return U0, U1, U0id, U1id

def Ren_gate(carbon_nr, B_field=304.22, phase=0):
    '''create a Ren gate for given carbon number'''

    number_of_pulses = mp['C' + str(carbon_nr) + '_Ren_N'][0]
    tau = mp['C' + str(carbon_nr) + '_Ren_tau'][0]
    print number_of_pulses
    print tau
    Ren = c13_gate(carbon_nr, number_of_pulses, tau, B_field=304.22)

    return Ren

def waittime(carbon_nr, time, B_field=304.22,return_indiv = False):
    '''calculates the evolution matrices for a single
    Carbon spin, electron is always qubit 1'''

    omega_Larmor = 2 * np.pi * B_field * 1.07e3
    A_par = 2 * np.pi * hf['C' + str(carbon_nr)]['par']
    A_perp = 2 * np.pi *hf['C' + str(carbon_nr)]['perp']

    H0 = omega_Larmor * Iz
    H1 = (A_par+omega_Larmor)*Iz + A_perp*Ix

    expH0 = (-1j*H0*time).expm();    expH1 = (-1j*H1*time).expm()
    Utot = qutip.tensor(rho0,expH0) + qutip.tensor(rho1,expH1)
    print_matrix(Utot)

    if return_indiv == False:
        return Utot
    if return_indiv == True:
        return expH0, expH1

def phase_gate(carbon_nr, phase, B_field=304.22,total_time = 0,return_gate = False, return_tau = False):
    '''calculates the evolution matrices for a single
    Carbon spin, electron is always qubit 1
    NOTE: now only works for total time = 0'''

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



def c13_gate_multiqubit(carbon_nrs, number_of_pulses, tau, B_field, gate_on_C = [], return_for_one = False, phase = None):
    '''calculates the evolution matrices a multiqubit space,
    the electron is always qubit 1
    To calculate the ideal case, you can give the function which C13 is adressed with gate_on_C
    This number should be the number of the C13 in the list in carbon_nrs. Ie: carbon_nrs = [1,2,4]
    gate on C2 means gate_on_C = [1]
    If return_for_one = True, only the carbon to be addressed is addressed, the others get the Id gate, to avoid phases in the simulation.
    NOTE: it only works for 3 C13 spins now, because of the way the final gate is constructed.
    '''
    U0={}
    U1={}
    U0_id={}
    U1_id={}

    gate_0 = gate_0_id = rho0
    gate_1 = gate_1_id = rho1

    for ii in range(len(carbon_nrs)):

        U0['C_'+str(ii)], U1['C_'+str(ii)], U0_id['C_'+str(ii)], U1_id['C_'+str(ii)] = c13_gate(carbon_nrs[ii], number_of_pulses, tau, B_field, return_indiv = True, return_id = True,  phase = phase)

        # if U0_id['C_'+str(ii)][0,0]== -1j*U0_id['C_'+str(ii)][0,1]:
        #     print 'Qubit '+str(ii)+' has a -/+ Ren gate, switched in simulation'
        #     U0_id['C_'+str(ii)], U1_id['C_'+str(ii)] =U1_id['C_'+str(ii)], U0_id['C_'+str(ii)]
        #     U0['C_'+str(ii)], U1['C_'+str(ii)] =U1['C_'+str(ii)], U0['C_'+str(ii)]

        if ii not in gate_on_C:
            U0_id['C_'+str(ii)] = Id
            U1_id['C_'+str(ii)] = Id
            if return_for_one == True:
                U0['C_'+str(ii)] = Id
                U1['C_'+str(ii)] = Id

        gate_0 = qutip.tensor(gate_0,U0['C_'+str(ii)])
        gate_1 = qutip.tensor(gate_1,U1['C_'+str(ii)])
        gate_0_id = qutip.tensor(gate_0_id,U0_id['C_'+str(ii)])
        gate_1_id = qutip.tensor(gate_1_id,U1_id['C_'+str(ii)])


    gate = gate_0+gate_1
    gate_id = gate_0_id+gate_1_id

    return gate, gate_id

def check_phase_gate():
    for phase in [np.pi]:#[0, np.pi/2, np.pi, np.pi*3/2, np.pi*2]:
        tau = phase_gate(1, phase, B_field=304.22,total_time = 0,return_gate = False, return_tau = True)
        U0, U1, U0id, U1id = c13_gate(1, 2, tau, 304.22, return_indiv = True, return_id = True)
        print phase
        print U0, U1, U0id, U1id

###################
### Pauli Sets ###
###################

def single_qubit_pauli(rho, do_plot = False):
    ii=-0.5
    pauli_set = []
    ii_list = []
    xticks_list = ['I','X','Y','Z']

    for oper in [Id, 2*sx,2*sy,2*sz]:

        pauli_set.append(qutip.expect(oper,rho))
        ii_list.append(ii)
        ii = ii+1


    if do_plot ==True:
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.bar(ii_list, pauli_set, width=1)
        plt.xticks(np.arange(0, 4, 1.0))
        ax.set_xticklabels(xticks_list)
        ax.set_xlim(-0.5,len(x_ticks_list)-0.5)
        ax.set_ylim(-1,1)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
            tick.label.set_rotation('vertical')
    return pauli_set, ii_list, xticks_list

def multi_qubit_pauli(rho,carbon_nrs=[1,1,1],do_plot=False, give_fid = False, alpha=None, beta=None,use_el=False,log_phase_corr=False, title = None):
    ''' This function works to perform two and three-qubit tomography
    it either just takes the expectation values of the given density matrix (when use_el = False)
    or performs a measurement of the expectation values using the electron spin (close to the experiment) NOTE: this does not work yet for two qubits
    '''
    no_qubits = len(carbon_nrs)
    ii=-0.5
    pauli_set = []
    ii_list = []
    xticks_list = ['X','Y','Z']
    final_x_tick_list = []
    oper_list = [2*sx,2*sy,2*sz]
    final_oper_list =[]

    if no_qubits ==2:
        for ff in xticks_list:
            final_x_tick_list.append(ff+'I')
        for ff in xticks_list:
            final_x_tick_list.append('I'+ff)

        for ff in oper_list:
            final_oper_list.append(qutip.tensor(ff,Id))
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(Id,ff))

        for jj in oper_list:
            for kk in oper_list:
                final_oper_list.append(qutip.tensor(jj,kk))

        for jj in xticks_list:
            for kk in xticks_list:
                final_x_tick_list.append(jj+kk)

        for jj in oper_list:
            for kk in oper_list:
                final_oper_list.append(qutip.tensor(jj,kk))


    if no_qubits ==3:
        for ff in xticks_list:
            final_x_tick_list.append(ff+'I'+'I')
        for ff in xticks_list:
            final_x_tick_list.append('I'+ff+'I')
        for ff in xticks_list:
            final_x_tick_list.append('I'+'I'+ff)

        for ff in oper_list:
            final_oper_list.append(qutip.tensor(ff,Id,Id))
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(Id,ff,Id))
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(Id,Id,ff))

        for jj in oper_list:
            for kk in oper_list:
                final_oper_list.append(qutip.tensor(jj,kk,Id))
        for jj in xticks_list:
            for kk in xticks_list:
                final_x_tick_list.append(jj+kk+'I')

        for jj in oper_list:
            for kk in oper_list:
                final_oper_list.append(qutip.tensor(jj,kk,Id))
        for jj in xticks_list:
            for kk in xticks_list:
                final_x_tick_list.append(jj+'I'+kk)

        for jj in oper_list:
            for kk in oper_list:
                final_oper_list.append(qutip.tensor(Id,jj,kk))
        for jj in xticks_list:
            for kk in xticks_list:
                final_x_tick_list.append('I'+jj+kk)


        for jj in oper_list:
            for kk in oper_list:
                for hh in oper_list:
                    final_oper_list.append(qutip.tensor(jj,kk,hh))
        for jj in xticks_list:
            for kk in xticks_list:
                for hh in xticks_list:
                    final_x_tick_list.append(jj+kk+hh)

    if use_el == False:
        for oper in final_oper_list:
            # print qutip.expect(oper,rho)
            pauli_set.append(qutip.expect(oper,rho))
            ii_list.append(ii)
            ii = ii+1


    if use_el == True:
        rho_in = qutip.tensor(rho0,rho)

        xel = qutip.tensor(x,Id,Id,Id)
        mxel = qutip.tensor(mx,Id,Id,Id)
        Xel = qutip.tensor(X,Id,Id,Id)
        yel = qutip.tensor(y,Id,Id,Id)
        myel = qutip.tensor(my,Id,Id,Id)

        tau_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_tau'][0]
        number_of_pulses_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_N'][0]
        Ren_C1, Ren_C1_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C1, tau_Ren_C1, 304.22, gate_on_C = [0], return_for_one = True)

        tau_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_tau'][0]
        number_of_pulses_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_N'][0]
        Ren_C2, Ren_C2_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C2, tau_Ren_C1, 304.22, gate_on_C = [1], return_for_one = True)


        tau_z_C1 = phase_gate(carbon_nrs[0], np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
        Rz_C1, Rz_C1_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C1, 304.22, gate_on_C = [0], return_for_one = True, phase = np.pi/2)


        tau_z_C2 = phase_gate(carbon_nrs[0], np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
        Rz_C2, Rz_C2_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C2, 304.22, gate_on_C = [1], return_for_one = True, phase = np.pi/2)


        tau_z_C1 = phase_gate(carbon_nrs[0], -np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
        Rmz_C1, Rmz_C1_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C1, 304.22, gate_on_C = [0], return_for_one = True, phase = -np.pi/2)
        total_time = 2*tau_Ren_C3*number_of_pulses_Ren_C3+ 2*tau_z_C1*2


        tau_z_C2 = phase_gate(carbon_nrs[0], -np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
        Rmz_C2, Rmz_C2_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C2, 304.22, gate_on_C = [1], return_for_one = True, phase = -np.pi/2)

        if no_qubits == 3:
            tau_Ren_C3 = mp['C' + str(carbon_nrs[2]) + '_Ren_tau'][0]
            number_of_pulses_Ren_C3 = mp['C' + str(carbon_nrs[2]) + '_Ren_N'][0]
            Ren_C3, Ren_C3_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C3, tau_Ren_C3, 304.22, gate_on_C = [2], return_for_one = True)

            tau_z_C3 = phase_gate(carbon_nrs[2], np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False, return_tau = True)
            Rz_C3, Rz_C3_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C3, 304.22, gate_on_C = [2], return_for_one = True, phase = np.pi/2)

            tau_z_C3 = phase_gate(carbon_nrs[2], -np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False, return_tau = True)
            Rmz_C3, Rmz_C3_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C3, 304.22, gate_on_C = [2], return_for_one = True, phase = -np.pi/2)


        Ren_gates = {}
        Ren_gates['C1']= Ren_C1_id
        Ren_gates['C2']= Ren_C2_id
        if no_qubits == 3:
            Ren_gates['C3']= Ren_C3_id

        Rz_gates = {}
        Rz_gates['C1']= Rz_C1_id
        Rz_gates['C2']= Rz_C2_id
        if no_qubits == 3:
            Rz_gates['C3']= Rz_C3_id

        Rmz_gates = {}
        Rmz_gates['C1']= Rmz_C1_id
        Rmz_gates['C2']= Rmz_C2_id
        if no_qubits == 3:
            Rmz_gates['C3']= Rmz_C3_id
    # print Ren_gates['C1']

        for ev in final_x_tick_list:
            ii_list.append(ii)
            ii = ii+1
                # single-qubit-ev
            if ev.count('I') == 2:
                if ev[0]!='I': C_nr = 'C1'
                elif ev[1]!='I': C_nr = 'C2'
                elif ev[2]!='I': C_nr = 'C3'
                seq_elm =myel* Ren_gates[C_nr]*xel
                if 'Z' in ev:
                    seq_elm = myel*Ren_gates[C_nr]*xel*Rz_gates[C_nr]*Ren_gates[C_nr]
                elif 'Y' in ev:
                    seq_elm = myel*Ren_gates[C_nr]*xel*Rz_gates[C_nr]
                    print 'Y'
            # two-qubit-ev
            elif ev.count('I') ==1:
                # print '2'
                if ev[0]=='I':
                    C_nr1 = 'C2'
                    C_nr2 = 'C3'
                elif ev[1]=='I':
                    C_nr1 = 'C1'
                    C_nr2 = 'C3'
                elif ev[2]=='I':
                    C_nr1 = 'C1'
                    C_nr2 = 'C2'

                seq_elm = xel*Ren_gates[C_nr1]*Ren_gates[C_nr2]*xel
                for i in range(3):
                    if ev[i] == 'Z':
                        seq_elm =seq_elm *Rz_gates['C'+str(i+1)] *Ren_gates['C'+str(i+1)]
                    elif ev[i] == 'Y':
                        seq_elm = seq_elm*Rz_gates['C'+str(i+1)]
            # three-qubit-ev
            elif 'I'not in ev:
                # print '3'
                seq_elm = yel*Ren_gates['C1']*Ren_gates['C2']*Ren_gates['C3']*xel
                for i in range(3):
                    if ev[i] == 'Z':
                        if log_phase_corr == True:
                            seq_elm = yel*Ren_gates['C1']*Ren_gates['C2']*Ren_gates['C3']*mxel
                        seq_elm =seq_elm*Rz_gates['C'+str(i+1)]* Ren_gates['C'+str(i+1)]
                    elif ev[i] == 'Y':
                        if log_phase_corr == True:
                            seq_elm = yel*Ren_gates['C1']*Ren_gates['C2']*Ren_gates['C3']*mxel
                        seq_elm =seq_elm*Rmz_gates['C'+str(i+1)]

            rho_el = (seq_elm*rho_in*seq_elm.dag()).ptrace(0)

            expect_value = qutip.expect(2*sz,rho_el)
            pauli_set.append(expect_value)

    if do_plot == True:
        if no_qubits == 3:
            figsize =(24,3.5)
        elif no_qubits == 2:
            figsize =(8,3.5)
        fig = plt.figure(figsize=figsize)
        ax = plt.subplot(111)
        ax.bar(ii_list, pauli_set, width=1)
        plt.xticks(np.arange(0, len(final_oper_list)-1, 1.0))
        ax.set_xticklabels(final_x_tick_list)
        ax.set_xlim(-0.5,len(final_x_tick_list)-0.5)
        # ax.set_ylim(-1,1)
        ax.set_ylim(-1,1)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
            tick.label.set_rotation('vertical')
        if title != None:
            ax.set_title(title)
        plt.show()
        print 'yes'



    if give_fid == True:
        psi_ideal = alpha*qutip.tensor(ketx,ketx,ketx)+beta*qutip.tensor(ketmx,ketmx,ketmx)
        rho_ideal = psi_ideal*psi_ideal.dag()
        Fid = qutip.fidelity(rho,rho_ideal)**2
        return pauli_set, ii_list, final_x_tick_list, Fid
    else:
        return pauli_set, ii_list, final_x_tick_list

###################
### Experiments ###
###################

def nuclear_rabi_no_init(carbon_nrs, tau, nr_of_pulses_list=np.linspace(0,300,76), B_field=304.225):
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

    S_tot = np.ones(len(nr_of_pulses_list))
    for kk, carbon_nr in enumerate(carbon_nrs):
        #sequence and RO
        S = np.zeros(len(nr_of_pulses_list))
        for i, N in enumerate(nr_of_pulses_list):
            gate = c13_gate(carbon_nr, N, tau, B_field)         # Define nuclear spin gate

            seq  = electron_mx*gate*electron_x                   # Define gate sequence
            rho_final = seq*rho_init*seq.dag()                  # Apply gate sequence

            rho_el_final = rho_final.ptrace(0)                  # Trace out the nuclear spin
            S[i] = qutip.expect(sz, rho_el_final)*2            # Z measurement two alternative ways
            #S[i] = qutip.fidelity(rho0, rho_el_final)**2

        ## Cumulative signal ##
        S_tot = S_tot*S


    ## plot ##
    f, ax = plt.subplots(1)
    ax.plot(nr_of_pulses_list, (S_tot+1)/2., 'o-', lw=1)
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
    Ren         = c13_gate(carbon_nr, N[0], tau[0], B_field)

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

def nuclear_ramsey_no_init_no_DD(carbon_nr, tau_wait_list, B_field=304.22):
    '''nuclear Rabi experiment without init or DD
    scheme: y - Ren - x - wait - y - Ren - x - RO'''

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
    electron_y = qutip.tensor(y,Id)
    Ren         = c13_gate(carbon_nr, N[0], tau[0], B_field)


    #sequence and RO
    S = np.zeros(len(tau_wait_list))
    for i, tau_wait in enumerate(tau_wait_list):

        wait_gate = waittime(carbon_nr, tau_wait, B_field)         # Define DD waiting gate

        seq1  = electron_x*Ren*electron_y                   # Define gate sequence
        rho_seq1 = seq1*rho_init*seq1.dag()

        seq2 = wait_gate
        rho_seq2 = seq2*rho_seq1*seq2.dag()
        #electron dephasing
        rho_seq2a =    0.5*qutip.tensor(z,Id)*rho_seq2*qutip.tensor(z,Id).dag() + 0.5*qutip.tensor(mz,Id)*rho_seq2*qutip.tensor(mz,Id).dag()

        seq3 = electron_x*Ren*electron_y
        rho_final = seq3*rho_seq2a*seq3.dag()


        rho_el_final = rho_final.ptrace(0)                  # Trace out the nuclear spin
        #S[i] = qutip.expect(sz, rho_el_final) + 1./2       # Z measurement two alternative ways
        S[i] = qutip.fidelity(rho0, rho_el_final)**2

    ## plot ##
    f, ax = plt.subplots(1)
    ax.plot(tau_wait_list, S, 'o-', lw=1)
    ax.set_title('P(ms=0)'); ax.set_xlabel('Evolution_time (us)')
    plt.show()
    return S[i]


######################
### Initialization ### NOTE: FOR MULTIPLE C13 SPINS THE PHASE GATES SHOULD BE ADDED!!
######################

def nuclear_init_single(carbon_nr,do_plot = False, method = 'SWAP'):
    '''function that returns a density matrix for an initialized C13 spin note: Z-gate not yet calculated in right way
    for method = SWAP    seq = y - Ren - x - Rz - Ren
    for  method = MBI    seq = y - Ren - x          nuclear spin initialized in |y> if electron was in |0>'''

    rho = qutip.tensor(rho0,rhom)

    yel = qutip.tensor(y,Id)
    xel = qutip.tensor(x,Id)

    U0, U1 = nuclear_Ren_matrix(carbon_nr,B_field=304.22)

    U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
    U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)
    Ren = qutip.tensor(rho0,U0)+qutip.tensor(rho1,U1)
    Ren_id = qutip.tensor(rho0,U0id)+qutip.tensor(rho1,U1id)

    Rz = qutip.tensor(Id,z)
    print_matrix(Rz)
    if method == 'SWAP':

        seq = Ren*Rz*xel*Ren*yel
        seq_id = Ren_id*Rz*xel*Ren_id*yel

        rho_final = seq*rho*seq.dag()
        rho_final_id = seq_id*rho*seq_id.dag()

        rho_nucl = rho_final.ptrace(1)
        rho_nucl_id = rho_final_id.ptrace(1)

    elif method == 'MBI':
        seq = xel*Ren*yel
        seq_id = xel*Ren_id*yel

        rho_final = seq*rho*seq.dag()
        rho_final_id = seq_id*rho*seq_id.dag()

        # measure electron to be in 0 and renormalize
        el0 = qutip.tensor(rho0,Id)
        rho_final_id = el0*rho_final_id*el0.dag()
        rho_final = el0*rho_final*el0.dag()

        # renormalize
        norm_id = qutip.fidelity(rho0,rho_final_id.ptrace([0]))
        norm = qutip.fidelity(rho0,rho_final.ptrace([0]))

        rho_nucl = 1/norm**2*rho_final.ptrace([1])
        rho_nucl_id = 1/(norm_id**2)*rho_final_id.ptrace([1])


    if do_plot == True:
        fig = plt.figure()
        ax = plt.subplot(111)

        for rho_n in [rho_nucl_id]:
            pauli_set, ii_list, x_ticks_list = single_qubit_pauli(rho_n)
            ax.bar(ii_list, np.real(pauli_set), width=1,alpha = 0.5)

        plt.xticks(np.arange(0, 4, 1.0))
        ax.set_xticklabels(x_ticks_list)
        ax.set_ylim(-1,1)
        print 'Fidelity to ideal state:'
        print qutip.fidelity(rho_nucl,rho_nucl_id)
        plt.show()


    return rho_nucl,rho_nucl_id


def three_spin_encoding(carbon_nrs = [1,1,1],alpha=1/np.sqrt(2),beta=1/np.sqrt(2),do_plot=True):
    ''' encodes your three chosen C13 spins in the state alpha(xxx)+beta(-x-x-x)
    note: it depends on the direction of the Ren gate for each C13 spin if it will be + or -x, see documentation'''

    # define density matrices
    rho_C1, rho_C1_id = nuclear_init_single(carbon_nrs[0])
    rho_C2, rho_C2_id = nuclear_init_single(carbon_nrs[1])
    rho_C3, rho_C3_id = nuclear_init_single(carbon_nrs[2])

    rho_el = any_pure_state(alpha,1j*beta,return_psi = False,return_rho = True)

    rho = qutip.tensor(rho_el,rho_C1,rho_C2,rho_C3)
    rho_id = qutip.tensor(rho_el,rho_C1_id,rho_C2_id,rho_C3_id)

    rho_C = qutip.tensor(rho_C1,rho_C2,rho_C3)
    rho_C_id = qutip.tensor(rho_C1_id,rho_C2_id,rho_C3_id)


    #define gates
    xel = qutip.tensor(x,Id,Id,Id)
    yel = qutip.tensor(y,Id,Id,Id)


    tau_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_N'][0]
    Ren_C1, Ren_C1_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C1, tau_Ren_C1, 304.22, gate_on_C = [0], return_for_one = True)

    tau_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_N'][0]
    Ren_C2, Ren_C2_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C2, tau_Ren_C1, 304.22, gate_on_C = [1], return_for_one = True)

    tau_Ren_C3 = mp['C' + str(carbon_nrs[2]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C3 = mp['C' + str(carbon_nrs[2]) + '_Ren_N'][0]
    Ren_C3, Ren_C3_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C3, tau_Ren_C3, 304.22, gate_on_C = [2], return_for_one = True)

    tau_z_C1 = phase_gate(carbon_nrs[0], np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
    Rz_C1, Rz_C1_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C1, 304.22, gate_on_C = [0], return_for_one = True, phase = np.pi/2)

    tau_z_C2 = phase_gate(carbon_nrs[0], np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
    Rz_C2, Rz_C2_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C2, 304.22, gate_on_C = [1], return_for_one = True, phase = np.pi/2)


    tau_z_C3 = phase_gate(carbon_nrs[2], np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False, return_tau = True)
    Rz_C3, Rz_C3_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C3, 304.22, gate_on_C = [2], return_for_one = True, phase = np.pi/2)


    # define and apply full sequence
    seq = xel*Rz_C3*Rz_C2*Rz_C1*Ren_C3*Ren_C2*Ren_C1
    seq_id = xel*Rz_C3_id*Rz_C2_id*Rz_C1_id*Ren_C3_id*Ren_C2_id*Ren_C1_id
    # seq_id = xel*Rz_C3_id*Rz_C2*Rz_C1_id*Ren_C3_id*Ren_C2_id*Ren_C1_id
    rho_after = seq*rho*seq.dag()
    rho_after_id = seq_id*rho_id*seq_id.dag()
    #measure electron to be in zero
    el0 = qutip.tensor(rho0,Id,Id,Id)
    rho_final_id = el0*rho_after_id*el0.dag()
    rho_final = el0*rho_after*el0.dag()

    # renormalize
    norm_id = qutip.fidelity(rho0,rho_final_id.ptrace([0]))
    norm = qutip.fidelity(rho0,rho_final.ptrace([0]))

    rho_final = 1/norm**2*rho_final.ptrace([1,2,3])
    rho_final_id = 1/(norm_id**2)*rho_final_id.ptrace([1,2,3])

    if do_plot == True:
        fig = plt.figure(figsize=(24,3.5))
        ax = plt.subplot(111)
        for rho_n in [rho_final_id]:
            pauli_set, ii_list, x_ticks_list = multi_qubit_pauli(rho_n)
            ax.bar(ii_list, np.real(pauli_set), width=1,alpha = 0.5)

        plt.xticks(np.arange(0, len(x_ticks_list), 1.0))
        ax.set_xticklabels(x_ticks_list)
    print 'Fidelity to ideal state:'
    print qutip.fidelity(rho_final,rho_final_id)

    return rho_final, rho_final_id

def check_entangled_state(alpha = 1./np.sqrt(2), beta = 1/np.sqrt(2), state ='+++'):
    ''' Creates a simple density matrix for a state that we wish to encode in.
     alplha xxx+beta -x-x-x and permutations: +-+ etc'''
 # define density matrices

    rho_el = any_pure_state(alpha,1j*beta,return_psi = False,return_rho = True)

    rho = qutip.tensor(rho_el,rho0,rho0,rho0)

    #define gates
    xel = qutip.tensor(x,Id,Id,Id)
    yel = qutip.tensor(y,Id,Id,Id)

    U0C1 = x
    U1C1 = mx
    U0C2 = x
    U1C2 = mx
    U0C3 = x
    U1C3 = mx

    if state[0] == '-':
        U0C1 = mx
        U1C1 = x
    elif state[1] == '-':
        U0C2 = mx
        U1C2 = x
    if state[2] == '-':
        U0C3 = mx
        U1C3 = x

    Ren_C1 = qutip.tensor(rho0,U0C1,Id,Id)+qutip.tensor(rho1,U1C1,Id,Id)
    Ren_C2 = qutip.tensor(rho0,Id,U0C2,Id)+qutip.tensor(rho1,Id,U1C2,Id)
    Ren_C3 = qutip.tensor(rho0,Id,Id,U0C3)+qutip.tensor(rho1,Id,Id,U1C3)

    Rz_C1 = qutip.tensor(Id,z,Id,Id)
    Rz_C2 = qutip.tensor(Id,Id,z,Id)
    Rz_C3 = qutip.tensor(Id,Id,Id,z)


    # define and apply full sequence
    seq = xel*Rz_C3*Rz_C2*Rz_C1*Ren_C3*Ren_C2*Ren_C1
    rho_after = seq*rho*seq.dag()

    #measure electron to be in zero
    el0 = qutip.tensor(rho0,Id,Id,Id)
    rho_final = el0*rho_after*el0.dag()

    norm = qutip.fidelity(rho0,rho_final.ptrace([0]))

    rho_final_1 = 1/norm**2*rho_final.ptrace([1,2,3])

    psi_final = alpha*qutip.tensor(ketx,ketx,ketx)+beta*qutip.tensor(ketmx,ketmx,ketmx)
    rho_final_2 = psi_final*psi_final.dag()

    fig = plt.figure(figsize=(24,3.5))
    ax = plt.subplot(111)

    pauli_set_1, ii_list, x_ticks_list = multi_qubit_pauli(rho_final_1)
    pauli_set_2, ii_list, x_ticks_list = multi_qubit_pauli(rho_final_2)
    ax.bar(ii_list, np.real(pauli_set_1), width=1,alpha = 0.5)
    ax.bar(ii_list, np.real(pauli_set_2), width=1,alpha = 0.5,color = 'red')

    plt.xticks(np.arange(0, len(x_ticks_list), 1.0))
    ax.set_xticklabels(x_ticks_list)
    # ax.set_title('entangled state = 1/sqrt(2)*(' + str(alpha*np.sqrt(2)) +' |xxx> +' + str(beta*np.sqrt(2)) + ' |-x-x-x>)')
    ax.set_title('entangled state = (' + str(alpha) +' |xxx> +' + str(beta) + ' |-x-x-x>)')
    ax.set_xlim(-0.5,len(x_ticks_list)-0.5)
    ax.set_ylim(-1,1)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
        tick.label.set_rotation('vertical')
    return rho_final

def two_spin_encoding(carbon_nrs = [1,1],alpha=1/np.sqrt(2),beta=1/np.sqrt(2)):
    ''' encodes your three chosen C13 spins in the state alpha(xxx)+beta(-x-x-x)
    note: it depends on the direction of the Ren gate for each C13 spin if it will be + or -x, see documentation'''

    # define density matrices
    rho_C1, rho_C1_id = nuclear_init_single(carbon_nrs[0])
    rho_C2, rho_C2_id = nuclear_init_single(carbon_nrs[1])

    psi_el = alpha*ket0+1j*beta*ket1
    rho_el = psi_el*psi_el.dag()

    rho = qutip.tensor(rho_el,rho_C1,rho_C2)
    rho_id = qutip.tensor(rho_el,rho_C1_id,rho_C2_id)

    rho_C = qutip.tensor(rho_C1,rho_C2)
    rho_C_id = qutip.tensor(rho_C1_id,rho_C2_id)


    #define gates
    xel = qutip.tensor(x,Id,Id)
    yel = qutip.tensor(y,Id,Id)

    U0, U1 = nuclear_Ren_matrix(carbon_nrs[0],B_field=304.22)
    U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
    U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)
    Ren_C1 = qutip.tensor(rho0,U0,Id)+qutip.tensor(rho1,U1,Id)
    Ren_C1_id = qutip.tensor(rho0,U0id,Id)+qutip.tensor(rho1,U1id,Id)
    Rz_C1 = qutip.tensor(Id,z,Id)

    U0, U1 = nuclear_Ren_matrix(carbon_nrs[1],B_field=304.22)
    U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
    U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)
    Ren_C2 = qutip.tensor(rho0,Id,U0)+qutip.tensor(rho1,Id,U1)
    Ren_C2_id = qutip.tensor(rho0,Id,U0id)+qutip.tensor(rho1,Id,U1id)
    Rz_C2 = qutip.tensor(Id,Id,z)


    # define and apply full sequence
    seq = xel*Rz_C2*Rz_C1*Ren_C2*Ren_C1
    seq_id = xel*Rz_C2*Rz_C1*Ren_C2_id*Ren_C1_id
    rho_after = seq*rho*seq.dag()
    rho_after_id = seq_id*rho_id*seq_id.dag()

    #measure electron to be in zero
    el0 = qutip.tensor(rho0,Id,Id)
    rho_final_id = el0*rho_after_id*el0.dag()
    rho_final = el0*rho_after*el0.dag()

    norm_id = qutip.fidelity(rho0,rho_final_id.ptrace([0]))
    norm = qutip.fidelity(rho0,rho_final.ptrace([0]))

    rho_final = 1/norm**2*rho_final.ptrace([1,2])
    rho_final_id = 1/norm_id**2*rho_final_id.ptrace([1,2])

    print 'Fidelity to ideal state:'
    print qutip.fidelity(rho,rho_id)

    fig = plt.figure(figsize=(60,4))
    ax = plt.subplot(111)
    for rho_n in [rho_final,rho_final_id]:
        pauli_set, ii_list, x_ticks_list = multi_qubit_pauli(rho_n)
        ax.bar(ii_list, np.real(pauli_set), width=1,alpha = 0.5)
    ax.set_xlim(-0.5,len(x_ticks_list)-0.5)
    ax.set_ylim(-1,1)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
        tick.label.set_rotation('vertical')
    plt.xticks(np.arange(0, len(x_ticks_list), 1.0))
    ax.set_xticklabels(x_ticks_list)
    plt.show()


#######################
### Error detection ###
#######################

def parity_msmt(qubits=[0,1],carbon_nrs = [1,1,1],alpha=1/np.sqrt(2),beta=1/np.sqrt(2),error_phase = 0, error_list = ['Q1','Q2'], do_plot = True, return_fid = False,return_state = False):
    '''implements an error or not on the two qubits that are selected, then parity is measured
    '''

    rho_enc_start, rho_enc_start_id = three_spin_encoding(carbon_nrs=carbon_nrs,alpha=alpha,beta=beta,do_plot=False)
    rho_enc = qutip.tensor(rho0,rho_enc_start)
    rho_enc_id = qutip.tensor(rho0,rho_enc_start_id)

    #define gates
    xel = qutip.tensor(x,Id,Id,Id)
    Xel = qutip.tensor(X,Id,Id,Id)
    yel = qutip.tensor(y,Id,Id,Id)

    tau_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_N'][0]
    Ren_C1, Ren_C1_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C1, tau_Ren_C1, 304.22, gate_on_C = [0], return_for_one = True)

    tau_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_N'][0]
    Ren_C2, Ren_C2_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C2, tau_Ren_C1, 304.22, gate_on_C = [1], return_for_one = True)

    tau_Ren_C3 = mp['C' + str(carbon_nrs[2]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C3 = mp['C' + str(carbon_nrs[2]) + '_Ren_N'][0]
    Ren_C3, Ren_C3_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C3, tau_Ren_C3, 304.22, gate_on_C = [2], return_for_one = True)

    total_time = 2*tau_Ren_C2*number_of_pulses_Ren_C2+2*tau_Ren_C3*number_of_pulses_Ren_C3

    tau_z_C1 = phase_gate(carbon_nrs[0], np.pi/2, B_field=304.22,total_time = total_time ,return_gate = False,return_tau = True)
    Rz_C1, Rz_C1_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C1, 304.22, gate_on_C = [0], return_for_one = True, phase = np.pi/2)
    total_time = 2*tau_Ren_C3*number_of_pulses_Ren_C3+ 2*tau_z_C1*2


    tau_z_C2 = phase_gate(carbon_nrs[0], np.pi/2, B_field=304.22,total_time = total_time ,return_gate = False,return_tau = True)
    Rz_C2, Rz_C2_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C2, 304.22, gate_on_C = [1], return_for_one = True, phase = np.pi/2)

    total_time = 2*tau_z_C1*2+2*tau_z_C2*2

    tau_z_C3 = phase_gate(carbon_nrs[2], np.pi/2, B_field=304.22,total_time = total_time ,return_gate = False, return_tau = True)
    Rz_C3, Rz_C3_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C3, 304.22, gate_on_C = [2], return_for_one = True, phase = np.pi/2)

    # implement error

    for error in error_list:
        if error =='Q1':
            tau_z_C1 = phase_gate(carbon_nrs[qubits[0]], error_phase, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
            Rz_C1, Rz_C1_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C1, 304.22, gate_on_C = [qubits[0]], return_for_one = True, phase = error_phase)
            rho_enc = Rz_C1*rho_enc*Rz_C1.dag()
            rho_enc_id = Rz_C1_id*rho_enc_id*Rz_C1_id.dag()

        if error =='Q2':
            tau_z_C2 = phase_gate(carbon_nrs[qubits[1]], error_phase, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
            Rz_C2, Rz_C2_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C2, 304.22, gate_on_C = [qubits[1]], return_for_one = True, phase = error_phase)
            rho_enc = Rz_C2*rho_enc*Rz_C2.dag()
            rho_enc_id = Rz_C2_id*rho_enc_id*Rz_C2_id.dag()

        if error =='Q3':
            tau_z_C3 = phase_gate(carbon_nrs[2], error_phase, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
            Rz_C3, Rz_C3_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C3, 304.22, gate_on_C = [2], return_for_one = True, phase = error_phase)
            rho_enc = Rz_C3*rho_enc*Rz_C3.dag()
            rho_enc_id = Rz_C3_id*rho_enc_id*Rz_C3_id.dag()

    # detect error
    seq = yel*Ren_C2*Ren_C1*yel
    seq_id = yel*Ren_C2_id*Ren_C1_id*yel
    rho_after =seq*rho_enc*seq.dag()
    rho_after_id =seq_id*rho_enc_id*seq_id.dag()

    # measure electron state
    rho_el_after = rho_after.ptrace(0)
    rho_el_after_id = rho_after_id.ptrace(0)
    rho_el_before = rho_enc_id.ptrace(0)

    print rho_el_after_id

    print 'Fidelity to zero:'
    print qutip.fidelity(rho0,rho_el_after)

    print 'Ideal fidelity to zero:'
    print qutip.fidelity(rho0,rho_el_after_id)
    if return_fid ==True:
        return qutip.fidelity(rho0,rho_el_after), qutip.fidelity(rho0,rho_el_after_id)

    if return_state ==True:
        return rho_after, rho_after_id

    if do_plot ==True:
        fig = plt.figure()
        ax = plt.subplot(111)
        for rho_el in [rho_el_after_id, rho0]:
            pauli_set, ii_list, x_ticks_list = single_qubit_pauli(rho_el, do_plot = False)
            ax.bar(ii_list, np.real(pauli_set), width=1,alpha = 0.5)
        ax.set_xlim(-0.5,len(x_ticks_list)-0.5)
        ax.set_ylim(-1,1)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
            tick.label.set_rotation('vertical')
        plt.xticks(np.arange(0, len(x_ticks_list), 1.0))
        ax.set_xticklabels(x_ticks_list)
        ax.set_title('electron state for encoded state 1/sqrt(2)'+str(alpha*np.sqrt(2))+ '|xxx>'+str(beta*np.sqrt(2))+ ' |-x-x-x>')

        #project electron
        el0 = qutip.tensor(rho_el_after_id,Id,Id,Id)
        rho_final_id = el0*rho_after_id*el0.dag()
        rho_final = el0*rho_after*el0.dag()

        norm_id = qutip.fidelity(rho_el_after_id,rho_final_id.ptrace([0]))
        norm = qutip.fidelity(rho_el_after_id,rho_final.ptrace([0]))

        rho_enc_final = 1/norm**2*rho_final.ptrace([1,2,3])
        rho_enc_final_id = 1/norm_id**2*rho_final_id.ptrace([1,2,3])

        fig = plt.figure(figsize=(25,4))
        ax = plt.subplot(111)
        rho_enc_exp = [rho_enc_start,rho_enc_final]
        colors = ["red", "blue"]
        labels = ['exp before parity msmst','exp after parity msmt']
        for ii in [0,1]:
            pauli_set, ii_list, x_ticks_list = multi_qubit_pauli(rho_enc_exp[ii], do_plot = False)
            ax.bar(ii_list, np.real(pauli_set), color = colors[ii], width=1,alpha = 0.5,label = labels[ii])
        ax.set_xlim(-0.5,len(x_ticks_list)-0.5)
        ax.set_ylim(-1,1)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
            tick.label.set_rotation('vertical')
        plt.xticks(np.arange(0, len(x_ticks_list), 1.0))
        ax.set_xticklabels(x_ticks_list)
        ax.legend()

        fig = plt.figure(figsize=(25,4))
        ax = plt.subplot(111)
        rho_enc =[rho_enc_start_id,rho_enc_final_id]
        colors = ["red", "blue"]
        labels = ['id before parity msmst','id after parity msmt']
        for ii in [0,1]:
            pauli_set, ii_list, x_ticks_list = multi_qubit_pauli(rho_enc[ii], do_plot = False)
            ax.bar(ii_list, np.real(pauli_set), color = colors[ii],width=1,alpha = 0.5,label = labels[ii])
        ax.set_xlim(-0.5,len(x_ticks_list)-0.5)
        ax.set_ylim(-1,1)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
            tick.label.set_rotation('vertical')
        plt.xticks(np.arange(0, len(x_ticks_list), 1.0))
        ax.set_xticklabels(x_ticks_list)
        ax.legend()

def error_curves():
    colors = [ "blue", "red", "green"]
    labels = ['Error on one qubit in parity msmt','Error both qubits in parity msmt','Error on third qubit' ]
    error_types = ['Q1','Q2','Q3']#[['Q1'],['Q1','Q2'],['Q3']]

    alpha_list =[1]# [1,0,1/np.sqrt(2),1/np.sqrt(2),1/np.sqrt(2),1/np.sqrt(2)]
    beta_list = [0]#[0,1,1/np.sqrt(2),-1/np.sqrt(2),1j/np.sqrt(2),-1j/np.sqrt(2)]
    for mm in range(len(alpha_list)):
        fig = plt.figure(figsize=(10,10))
        ax = plt.subplot(111)
        for jj in range(len(error_types)):
            error_phase_list = np.linspace(0,np.pi,10)
            fid = np.zeros(len(error_phase_list))
            fid_id = np.zeros(len(error_phase_list))
            for ii in range(len(error_phase_list)):
                fid[ii], fid_id[ii] = parity_msmt(qubits=[0,1],carbon_nrs = [1,1,1],alpha=alpha_list[mm],beta=beta_list[mm],error_phase = error_phase_list[ii], error_list = error_types[jj], do_plot = False, return_fid = True)
            error_phase_list = error_phase_list/np.pi
            fid_sq = [i ** 2 for i in fid]
            # ax.plot(error_phase_list, fid_sq, color = colors[jj],label = labels[jj])
            fid_id_sq = [i ** 2 for i in fid_id]
            ax.plot(error_phase_list, fid_id_sq, color = colors[jj],label = labels[jj]+' ideal',linestyle = 'dashed')

        ax.set_ylim(-0.1,1.1)
        ax.legend(loc=3)
        ax.set_xlabel('Error_angle (units of pi')
        ax.set_ylabel('Electron probability to measure |0>')
        ax.set_title('state = (' + str(alpha_list[mm]) +' |xxx> + ' + str(beta_list[mm]) + ' |-x-x-x>)')

def full_QEC_experiment(alpha=1/np.sqrt(2),beta=1/np.sqrt(2),carbon_nrs = [4,4,4],points = 5,do_plot=True,reset = True):

    error_angle_list = np.linspace(0,np.pi,points)

    rho_enc_start, rho_enc_start_id = three_spin_encoding(carbon_nrs=carbon_nrs,alpha=alpha,beta=beta,do_plot=False)
    rho_enc = qutip.tensor(rho0,rho_enc_start)
    rho_enc_id = qutip.tensor(rho0,rho_enc_start_id)

    #define gates
    xel = qutip.tensor(x,Id,Id,Id)
    mxel = qutip.tensor(mx,Id,Id,Id)
    Xel = qutip.tensor(X,Id,Id,Id)
    yel = qutip.tensor(y,Id,Id,Id)
    myel = qutip.tensor(my,Id,Id,Id)

    tau_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_N'][0]
    Ren_C1, Ren_C1_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C1, tau_Ren_C1, 304.22, gate_on_C = [0], return_for_one = True)

    tau_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_N'][0]
    Ren_C2, Ren_C2_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C2, tau_Ren_C1, 304.22, gate_on_C = [1], return_for_one = True)

    tau_Ren_C3 = mp['C' + str(carbon_nrs[2]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C3 = mp['C' + str(carbon_nrs[2]) + '_Ren_N'][0]
    Ren_C3, Ren_C3_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C3, tau_Ren_C3, 304.22, gate_on_C = [2], return_for_one = True)

    tau_z_C1 = phase_gate(carbon_nrs[0], np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
    Rz_C1, Rz_C1_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C1, 304.22, gate_on_C = [0], return_for_one = True, phase = np.pi/2)
    total_time = 2*tau_Ren_C3*number_of_pulses_Ren_C3+ 2*tau_z_C1*2


    tau_z_C2 = phase_gate(carbon_nrs[0], np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
    Rz_C2, Rz_C2_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C2, 304.22, gate_on_C = [1], return_for_one = True, phase = np.pi/2)


    tau_z_C3 = phase_gate(carbon_nrs[2], np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False, return_tau = True)
    Rz_C3, Rz_C3_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C3, 304.22, gate_on_C = [2], return_for_one = True, phase = np.pi/2)

    # implement error

    error_probability= np.zeros(len(error_angle_list))
    F_0_simple = np.zeros(len(error_angle_list))
    F_temp_0 = np.zeros(len(error_angle_list))
    F_temp_1 = np.zeros(len(error_angle_list))
    F_temp_00 = np.zeros(len(error_angle_list))
    F_temp_01 = np.zeros(len(error_angle_list))
    F_temp_10 = np.zeros(len(error_angle_list))
    F_temp_11 = np.zeros(len(error_angle_list))

    Fid_nc = np.zeros(len(error_angle_list))
    Fid_c = np.zeros(len(error_angle_list))
    Fid_pc = np.zeros(len(error_angle_list))

    Three_Fid_nc = np.zeros(len(error_angle_list))
    Three_Fid_c = np.zeros(len(error_angle_list))
    Three_Fid_pc = np.zeros(len(error_angle_list))



    for ii in range(len(error_angle_list)):
        error_phase = error_angle_list[ii]
        print 'error_angle '+str(error_phase)
        rho_enc_exp_id = rho_enc_id
        rho_enc_exp = rho_enc

        for error in ['Q1','Q2','Q3']:
            if error =='Q1':
                tau_z_C1 = phase_gate(carbon_nrs[0], error_phase, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
                Rz_C1, Rz_C1_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C1, 304.22, gate_on_C = [0], return_for_one = True, phase = error_phase)
                rho_enc_exp = Rz_C1*rho_enc_exp*Rz_C1.dag()
                rho_enc_exp_id = Rz_C1_id*rho_enc_exp_id*Rz_C1_id.dag()

            if error =='Q2':
                tau_z_C2 = phase_gate(carbon_nrs[1], error_phase, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
                Rz_C2, Rz_C2_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C2, 304.22, gate_on_C = [1], return_for_one = True, phase = error_phase)
                rho_enc_exp = Rz_C2*rho_enc_exp*Rz_C2.dag()
                rho_enc_exp_id = Rz_C2_id*rho_enc_exp_id*Rz_C2_id.dag()

            if error =='Q3':
                tau_z_C3 = phase_gate(carbon_nrs[2], error_phase, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
                Rz_C3, Rz_C3_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C3, 304.22, gate_on_C = [2], return_for_one = True, phase = error_phase)
                rho_enc_exp = Rz_C3*rho_enc_exp*Rz_C3.dag()
                rho_enc_exp_id = Rz_C3_id*rho_enc_exp_id*Rz_C3_id.dag()

        #first parity msmt

        # detect error
        seq = yel*Ren_C2*Ren_C1*yel
        seq_id = yel*Ren_C2_id*Ren_C1_id*yel

        # seq = yel*Ren_C3*Ren_C2*yel
        # seq_id = yel*Ren_C3_id*Ren_C2_id*yel
        rho_after =seq*rho_enc_exp*seq.dag()
        rho_after_id =seq_id*rho_enc_exp_id*seq_id.dag()

        # measure electron state
        rho_el_after = rho_after.ptrace(0)
        rho_el_after_id = rho_after_id.ptrace(0)
        rho_el_before = rho_enc_id.ptrace(0)
        F_0_simple[ii] = qutip.fidelity(rho_el_after_id,rho0)**2

        # print ' after first'
        # print rho_el_after_id

        # measure electron state in 0
        el0 = qutip.tensor(rho0,Id,Id,Id)
        rho_final_0_id = el0*rho_after_id*el0.dag()
        rho_final_0 = el0*rho_after*el0.dag()

        F_temp_0[ii] =  qutip.fidelity(rho_final_0_id.ptrace(0),rho0)**2

            # measure electron state in 1 and flip to 0 if reset = True

        el1 = qutip.tensor(rho1,Id,Id,Id)
        rho_final_id = el1*rho_after_id*el1.dag()
        rho_final = el1*rho_after*el1.dag()

        if reset == True:
            rho_final_1_id = Xel*rho_final_id*Xel.dag()
            rho_final_1 = Xel*rho_final*Xel.dag()
        else:
            rho_final_1_id = rho_final_id
            rho_final_1 = rho_final

        F_temp_1[ii] =  qutip.fidelity(rho_final_1_id.ptrace(0),rho0)**2

        # Second parity measurement after measured 0

            # detect error
        seq = yel*Ren_C3*Ren_C2*yel
        seq_id = yel*Ren_C3_id*Ren_C2_id*yel

        rho_after_0 =seq*rho_final_0*seq.dag()
        rho_after_0_id =seq_id*rho_final_0_id*seq_id.dag()

            # measure electron state
        rho_el_after_0 = rho_after_0.ptrace(0)
        rho_el_after_0_id = rho_after_0_id.ptrace(0)
        # F_00_simple[ii] = qutip.fidelity(rho_el_after_0_id,rho0)**2
        # print ' after 0'
        # print rho_el_after_0_id
            # measure electron state in 0
        rho_final_00_id = el0*rho_after_0_id*el0.dag()
        rho_final_00 = el0*rho_after_0*el0.dag()

        F_temp_00[ii] =  qutip.fidelity(rho_final_00_id.ptrace(0),rho0)**2

            # measure electron state in 1 and reset

        rho_final_id = el1*rho_after_0_id*el1.dag()
        rho_final = el1*rho_after_0_id*el1.dag()

        if reset == True:
            rho_final_01_id = Xel*rho_final_id*Xel.dag()
            rho_final_01 = Xel*rho_final*Xel.dag()
        else:
            rho_final_01_id = rho_final_id
            rho_final_01 = rho_final

        F_temp_01[ii] =  qutip.fidelity(rho_final_01_id.ptrace(0),rho0)**2

        # Second parity measurement after measured 1

            # detect error

        seq = yel*Ren_C3*Ren_C2*yel
        seq_id = yel*Ren_C3_id*Ren_C2_id*yel

        rho_after_1 =seq*rho_final_1*seq.dag()
        rho_after_1_id =seq_id*rho_final_1_id*seq_id.dag()

            # measure electron state
        rho_el_after_1 = rho_after_1.ptrace(0)
        rho_el_after_1_id = rho_after_1_id.ptrace(0)
        # F_01_simple[ii] = qutip.fidelity(rho_el_after_1_id,rho0)**2
        # print ' after 1'
        # print rho_el_after_1_id

            # measure electron state in 0
        rho_final_10_id = el0*rho_after_1_id*el0.dag()
        rho_final_10 = el0*rho_after_1*el0.dag()


        F_temp_10[ii] =  qutip.fidelity(rho_final_10_id.ptrace(0),rho0)**2

            # measure electron state in 1 and reset
        rho_final_id = el1*rho_after_1_id*el1.dag()
        rho_final = el1*rho_after_1_id*el1.dag()

        rho_final_11_id = Xel*rho_final_id*Xel.dag()
        rho_final_11 = Xel*rho_final*Xel.dag()

        F_temp_11[ii] =  qutip.fidelity(rho_final_11_id.ptrace(0),rho0)**2


        tau_z_C1 = phase_gate(carbon_nrs[0], np.pi, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
        RZ_C1, RZ_C1_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C1, 304.22, gate_on_C = [0], return_for_one = True, phase = np.pi)

        tau_z_C2 = phase_gate(carbon_nrs[0], np.pi, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
        RZ_C2, RZ_C2_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C2, 304.22, gate_on_C = [1], return_for_one = True, phase = np.pi)

        tau_z_C3 = phase_gate(carbon_nrs[2], np.pi, B_field=304.22,total_time = 0 ,return_gate = False, return_tau = True)
        RZ_C3, RZ_C3_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C3, 304.22, gate_on_C = [2], return_for_one = True, phase = np.pi)

        # Apply corrections (pr around Z for qubit that shows error, Pi around X (using double Ren gate) on one qubit for 10 and 01 outcome)

        seq_00 = RZ_C2*Ren_C2*Ren_C2
        seq_01 = RZ_C1
        seq_10 = RZ_C3
        seq_11 = Ren_C1*Ren_C1

        if reset == True:
            seq_00_id = qutip.tensor(Id,Id,Id,Id)#RZ_C2_id#*RZ_C1_id
            seq_11_id = RZ_C2_id
            seq_01_id = RZ_C3_id*Ren_C3_id*Ren_C3_id#RZ_C2_id#qutip.tensor(Id,Id,Id,Id)
            seq_10_id =  RZ_C1_id*Ren_C2_id*Ren_C2_id #Ren_C3_id*Ren_C3_id #RZ_C3_id*Ren_C3_id*Ren_C3_id

        else:
            # For y my
            # seq_00_id = RZ_C2_id#RZ_C2_id
            # seq_11_id = RZ_C3_id*Ren_C1_id*Ren_C1_id
            # seq_01_id = RZ_C1_id*Ren_C3_id*Ren_C3_id#RZ_C2_id#qutip.tensor(Id,Id,Id,Id)
            # seq_10_id = qutip.tensor(Id,Id,Id,Id)  #Ren_C3_id*Ren_C3_id #RZ_C3_id*Ren_C3_id*Ren_C3_id

            #for y y
            seq_10_id = RZ_C2_id#*RZ_C1_id
            seq_01_id = Ren_C3_id*Ren_C3_id*RZ_C3_id
            seq_11_id = Ren_C1_id*Ren_C1_id*RZ_C1_id#RZ_C2_id#qutip.tensor(Id,Id,Id,Id)
            seq_00_id = qutip.tensor(Id,Id,Id,Id)  #Ren_C3_id*Ren_C3_id #RZ_C3_id*Ren_C3_id*Ren_C3_id

        rho_final_10_c = seq_10*rho_final_10*seq_10.dag()
        rho_final_10_c_id = seq_10_id*rho_final_10_id*seq_10_id.dag()

        rho_final_01_c = seq_01*rho_final_01*seq_01.dag()
        rho_final_01_c_id = seq_01_id*rho_final_01_id*seq_01_id.dag()

        rho_final_00_c = seq_00*rho_final_00*seq_00.dag()
        rho_final_00_c_id = seq_00_id*rho_final_00_id*seq_00_id.dag()

        rho_final_11_c = seq_11*rho_final_11*seq_11.dag()
        rho_final_11_c_id = seq_11_id*rho_final_11_id*seq_11_id.dag()

        rho_total = rho_final_00+rho_final_01+rho_final_10+rho_final_11
        rho_total_id = rho_final_00_id+rho_final_01_id+rho_final_10_id+rho_final_11_id

        rho_total_c = rho_final_00_c+rho_final_01_c+rho_final_10_c+rho_final_11_c
        rho_total_c_id = rho_final_00_c_id+rho_final_01_c_id+rho_final_10_c_id+rho_final_11_c_id

        # apply only logical phase correction

        # seq_00 = Ren_C2*Ren_C2
        # seq_01 = RZ_C1
        # seq_10 = RZ_C3
        # seq_11 = Ren_C1*Ren_C1
        if reset == True:
            seq_11_id = qutip.tensor(Id,Id,Id,Id)
            seq_10_id = Ren_C1_id*Ren_C1_id
            seq_01_id = Ren_C3_id*Ren_C3_id
            seq_00_id = qutip.tensor(Id,Id,Id,Id)
        else:
            seq_00_id = qutip.tensor(Id,Id,Id,Id)
            seq_01_id = Ren_C1_id*Ren_C1_id
            seq_11_id = Ren_C3_id*Ren_C3_id
            seq_10_id = qutip.tensor(Id,Id,Id,Id)

        rho_final_10_pc = seq_10*rho_final_10*seq_10.dag()
        rho_final_10_pc_id = seq_10_id*rho_final_10_id*seq_10_id.dag()

        rho_final_01_pc = seq_01*rho_final_01*seq_01.dag()
        rho_final_01_pc_id = seq_01_id*rho_final_01_id*seq_01_id.dag()

        rho_final_00_pc = seq_00*rho_final_00*seq_00.dag()
        rho_final_00_pc_id = seq_00_id*rho_final_00_id*seq_00_id.dag()

        rho_final_11_pc = seq_11*rho_final_11*seq_11.dag()
        rho_final_11_pc_id = seq_11_id*rho_final_11_id*seq_11_id.dag()

        rho_total = rho_final_00+rho_final_01+rho_final_10+rho_final_11
        rho_total_id = rho_final_00_id+rho_final_01_id+rho_final_10_id+rho_final_11_id

        rho_total_pc = rho_final_00_pc+rho_final_01_pc+rho_final_10_pc+rho_final_11_pc
        rho_total_pc_id = rho_final_00_pc_id+rho_final_01_pc_id+rho_final_10_pc_id+rho_final_11_pc_id


        # measure fidelities

        error_probability[ii] = np.sin(error_phase/2)**2

        pauli_set_nc, ii_list, x_ticks_list, Three_Fid_nc[ii] = multi_qubit_pauli(rho_total_id.ptrace([1,2,3]),give_fid = True, alpha=alpha, beta=beta)
        pauli_set_c, ii_list, x_ticks_list, Three_Fid_c[ii]= multi_qubit_pauli(rho_total_c_id.ptrace([1,2,3]),give_fid = True, alpha=alpha, beta=beta)
        pauli_set_pc, ii_list, x_ticks_list, Three_Fid_pc[ii]= multi_qubit_pauli(rho_total_pc_id.ptrace([1,2,3]),give_fid = True, alpha=alpha, beta=beta)
        pauli_set_init, ii_list, x_ticks_list = multi_qubit_pauli(rho_enc_start_id,give_fid = False, alpha=alpha, beta=beta)

        # fig = plt.figure()
        # ax = plt.subplot(111)
        # ax.bar(ii_list, np.real(pauli_set_c), width=1,alpha = 0.5,color='blue')
        # # ax.bar(ii_list, np.real(pauli_set_pc), width=1,alpha = 0.5, color = 'red')
        # ax.bar(ii_list, np.real(pauli_set_init), width=1,alpha = 0.5, color = 'green')
        # ax.set_xlim(-0.5,len(x_ticks_list)-0.5)
        # ax.set_ylim(-1,1)
        # ax.set_xticklabels(x_ticks_list)
        # ax.set_title('state = (' + str(alpha) +' |xxx> + ' + str(beta) + ' |-x-x-x>)')

        if alpha == 1 and beta == 0:
            Fid_nc[ii] = 1/2.*(1+pauli_set_nc[0])
            Fid_c[ii] = 1/2.*(1+pauli_set_c[0])
            Fid_pc[ii] = 1/2.*(1+pauli_set_pc[0])
        elif alpha == 0 and beta == 1:
            Fid_nc[ii] = 1/2.*(1-pauli_set_nc[0])
            Fid_c[ii] = 1/2.*(1-pauli_set_c[0])
            Fid_pc[ii] = 1/2.*(1-pauli_set_pc[0])
        elif alpha == 1/np.sqrt(2) and beta == 1/np.sqrt(2):
            Fid_nc[ii] = 1/2.*(1+pauli_set_nc[62])
            Fid_c[ii] = 1/2.*(1+pauli_set_c[62])
            Fid_pc[ii] = 1/2.*(1+pauli_set_pc[62])
        elif alpha == 1/np.sqrt(2) and beta == -1/np.sqrt(2):
            Fid_nc[ii] = 1/2.*(1-pauli_set_nc[62])
            Fid_c[ii] = 1/2.*(1-pauli_set_c[62])
            Fid_pc[ii] = 1/2.*(1-pauli_set_pc[62])
        elif alpha == 1/np.sqrt(2) and beta == 1j/np.sqrt(2):
            Fid_nc[ii] = 1/2.*(1-pauli_set_nc[53])
            Fid_c[ii] = 1/2.*(1-pauli_set_c[53])
            Fid_pc[ii] = 1/2.*(1-pauli_set_pc[53])
        elif alpha == 1/np.sqrt(2) and beta == -1j/np.sqrt(2):
            Fid_nc[ii] = 1/2.*(1+pauli_set_nc[53])
            Fid_c[ii] = 1/2.*(1+pauli_set_c[53])
            Fid_pc[ii] = 1/2.*(1+pauli_set_pc[53])


    if do_plot == True:
        fig = plt.figure(figsize=(10,10))
        ax = plt.subplot(111)
        ax.plot(error_probability,Fid_nc,color='blue',label='no correction', linestyle='solid')
        ax.plot(error_probability,Fid_c,color='red', label = 'correction')
        ax.plot(error_probability,Fid_pc,color='green', label = 'logical phase correction', linestyle='dashed')

        ax.set_ylim(-0.1,1.1)
        ax.legend(loc=3)
        ax.set_xlabel('Error probability')
        ax.set_ylabel('Fidelity initial ideal (decoded) state')
        ax.set_title('Error curve state = (' + str(alpha) +' |xxx> + ' + str(beta) + ' |-x-x-x>)')

        fig = plt.figure(figsize=(10,10))
        ax = plt.subplot(111)
        ax.plot(error_probability,Three_Fid_nc,color='blue',label='no correction', linestyle='solid')
        ax.plot(error_probability,Three_Fid_c,color='red', label = 'correction')
        ax.plot(error_probability,Three_Fid_pc,color='green', label = 'logical phase correction', linestyle='dashed')

        ax.set_ylim(-0.1,1.1)
        ax.legend(loc=3)
        ax.set_xlabel('Error probability')
        ax.set_ylabel('Three-qubit Fidelity initial ideal state')
        ax.set_title('Error curve state = (' + str(alpha) +' |xxx> + ' + str(beta) + ' |-x-x-x>)')

    return error_probability, Fid_nc, Fid_pc, Fid_c

def full_error_curves(points = 15):
    state_list = ['z','mz','x','mx','y','my']
    alpha_list = [1,0,1/np.sqrt(2),1/np.sqrt(2),1/np.sqrt(2),1/np.sqrt(2)]
    beta_list = [0,1,1/np.sqrt(2),-1/np.sqrt(2),1j/np.sqrt(2),-1j/np.sqrt(2)]
    Fid_nc = {}
    Fid_pc = {}
    Fid_c = {}
    for ii in range(len(alpha_list)):
        error_probability, Fid_nc[state_list[ii]] , Fid_pc[state_list[ii]] , Fid_c[state_list[ii]] = full_QEC_experiment(alpha=alpha_list[ii],beta=beta_list[ii],points = points,do_plot=True)

    PF_nc = 1/4.*(Fid_nc['z']+Fid_nc['mz']+Fid_nc['x']+Fid_nc['mx']+Fid_nc['y']+Fid_nc['my'])-1/2.
    PF_pc = 1/4.*(Fid_pc['z']+Fid_pc['mz']+Fid_pc['x']+Fid_pc['mx']+Fid_pc['y']+Fid_pc['my'])-1/2.
    PF_c = 1/4.*(Fid_c['z']+Fid_c['mz']+Fid_c['x']+Fid_c['mx']+Fid_c['y']+Fid_c['my'])-1/2.
    # fig = plt.figure(figsize=(10,10))
    # ax = plt.subplot(111)
    # ax.plot(error_probability,PF_pc,color='blue',label='no correction')
    # ax.plot(error_probability,PF_c,color='red', label = 'correction')

    # ax.set_ylim(-0.1,1.1)
    # ax.legend(loc=3)
    # ax.set_xlabel('Error probability')
    # ax.set_ylabel('Process fidelity initial ideal (decoded) state')

def three_qubit_ev_via_el(rho,carbon_nrs = [1,1,1],ev = 'ZZZ'):

    rho_in = qutip.tensor(rho0,rho)
    print rho_in
    #define gates
    xel = qutip.tensor(x,Id,Id,Id)
    mxel = qutip.tensor(mx,Id,Id,Id)
    Xel = qutip.tensor(X,Id,Id,Id)
    yel = qutip.tensor(y,Id,Id,Id)
    myel = qutip.tensor(my,Id,Id,Id)

    tau_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_N'][0]
    Ren_C1, Ren_C1_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C1, tau_Ren_C1, 304.22, gate_on_C = [0], return_for_one = True)

    tau_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_N'][0]
    Ren_C2, Ren_C2_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C2, tau_Ren_C1, 304.22, gate_on_C = [1], return_for_one = True)

    tau_Ren_C3 = mp['C' + str(carbon_nrs[2]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C3 = mp['C' + str(carbon_nrs[2]) + '_Ren_N'][0]
    Ren_C3, Ren_C3_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C3, tau_Ren_C3, 304.22, gate_on_C = [2], return_for_one = True)

    tau_z_C1 = phase_gate(carbon_nrs[0], np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
    Rz_C1, Rz_C1_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C1, 304.22, gate_on_C = [0], return_for_one = True, phase = np.pi/2)
    total_time = 2*tau_Ren_C3*number_of_pulses_Ren_C3+ 2*tau_z_C1*2


    tau_z_C2 = phase_gate(carbon_nrs[0], np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
    Rz_C2, Rz_C2_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C2, 304.22, gate_on_C = [1], return_for_one = True, phase = np.pi/2)


    tau_z_C3 = phase_gate(carbon_nrs[2], -np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False, return_tau = True)
    Rmz_C3, Rmz_C3_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C3, 304.22, gate_on_C = [2], return_for_one = True, phase = -np.pi/2)

    tau_z_C1 = phase_gate(carbon_nrs[0], -np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
    Rmz_C1, Rmz_C1_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C1, 304.22, gate_on_C = [0], return_for_one = True, phase = -np.pi/2)
    total_time = 2*tau_Ren_C3*number_of_pulses_Ren_C3+ 2*tau_z_C1*2


    tau_z_C2 = phase_gate(carbon_nrs[0], -np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
    Rmz_C2, Rmz_C2_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C2, 304.22, gate_on_C = [1], return_for_one = True, phase = -np.pi/2)


    tau_z_C3 = phase_gate(carbon_nrs[2], -np.pi/2, B_field=304.22,total_time = 0 ,return_gate = False, return_tau = True)
    Rz_C3, Rz_C3_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C3, 304.22, gate_on_C = [2], return_for_one = True, phase = -np.pi/2)

    Ren_gates = {}
    Ren_gates['C1']= Ren_C1_id
    Ren_gates['C2']= Ren_C2_id
    Ren_gates['C3']= Ren_C3_id

    Rz_gates = {}
    Rz_gates['C1']= Rz_C1_id
    Rz_gates['C2']= Rz_C2_id
    Rz_gates['C3']= Rz_C3_id

    Rmz_gates = {}
    Rmz_gates['C1']= Rmz_C1_id
    Rmz_gates['C2']= Rmz_C2_id
    Rmz_gates['C3']= Rmz_C3_id
    # print Ren_gates['C1']

    # single-qubit-ev

    if ev.count('I') == 2:
        if ev[0]!='I': C_nr = 'C1'
        elif ev[1]!='I': C_nr = 'C2'
        elif ev[2]!='I': C_nr = 'C3'
        seq_elm =myel* Ren_gates[C_nr]*xel
        if 'Z' in ev:
            seq_elm = myel*Ren_gates[C_nr]*xel*Rz_gates[C_nr]*Ren_gates[C_nr]
        elif 'Y' in ev:
            seq_elm = myel*Ren_gates[C_nr]*xel*Rz_gates[C_nr]
            print 'Y'

    # two-qubit-ev
    elif ev.count('I') ==1:
        print '2'
        if ev[0]=='I':
            C_nr1 = 'C2'
            C_nr2 = 'C3'
        elif ev[1]=='I':
            C_nr1 = 'C1'
            C_nr2 = 'C3'
        elif ev[2]=='I':
            C_nr1 = 'C1'
            C_nr2 = 'C2'

        seq_elm = xel*Ren_gates[C_nr1]*Ren_gates[C_nr2]*xel
        for i in range(3):
            if ev[i] == 'Z':
                seq_elm =seq_elm *Rz_gates['C'+str(i+1)] *Ren_gates['C'+str(i+1)]
            elif ev[i] == 'Y':
                seq_elm = seq_elm*Rz_gates['C'+str(i+1)]

    # three-qubit-ev
    elif 'I'not in ev:
        print '3'
        seq_elm = myel*Ren_gates['C1']*Ren_gates['C2']*Ren_gates['C3']*xel
        for i in range(3):
            if ev[i] == 'Z':
                seq_elm =seq_elm*Rz_gates['C'+str(i+1)]* Ren_gates['C'+str(i+1)]
            elif ev[i] == 'Y':
                seq_elm =seq_elm*Rz_gates['C'+str(i+1)]

    print seq_elm
    rho_el = (seq_elm*rho_in*seq_elm.dag()).ptrace(0)
    # print rho_el
    #measure electron
    # rho_el = rho_after.ptrace(0)

    expect_value = qutip.expect(2*sz,rho_el)

    return expect_value

def test_pauli_sets():
    alpha_list = [1,1/np.sqrt(2),1/np.sqrt(2)]# [1,0,1/np.sqrt(2),1/np.sqrt(2),1/np.sqrt(2),1/np.sqrt(2)]
    beta_list = [0,1/np.sqrt(2),1j/np.sqrt(2)]#[0,1,1/np.sqrt(2),-1/np.sqrt(2),1j/np.sqrt(2),-1j/np.sqrt(2)]
    for ii in range(len(alpha_list)):
        rho,rho_id = three_spin_encoding(carbon_nrs = [1,1,1],alpha=alpha_list[ii],beta=beta_list[ii],do_plot=False)
        multi_qubit_pauli(rho_id,do_plot=True,use_el=False)
        multi_qubit_pauli(rho_id,do_plot=True,use_el=True)

def two_qb_entanglement_parity(carbon_nrs = [1,4], initial_states = 'ZZ' ,combination = '-X-X'):
    ''' this function performs a parity meassurement between two C13 spins on one of the 4 combinations ( XX, X-X, -XX,-X-X)
    on a certain initial product state (either 'ZZ' or 'XX')
    '''
    if initial_states == 'ZZ':
        method = 'SWAP'
    elif initial_states == 'XX':
        method = 'MBI'

    # initialize Carbon spins
    rho_C1, rho_C1_id = nuclear_init_single(carbon_nrs[0],method = method)
    rho_C2, rho_C2_id = nuclear_init_single(carbon_nrs[1],method = method)

    rho_enc = qutip.tensor(rho0,rho_C1,rho_C2)
    rho_enc_id = qutip.tensor(rho0,rho_C1_id,rho_C2_id)

    #define gates
    xel = qutip.tensor(x,Id,Id)
    Xel = qutip.tensor(X,Id,Id)
    yel = qutip.tensor(y,Id,Id)

    tau_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C1 = mp['C' + str(carbon_nrs[0]) + '_Ren_N'][0]
    Ren_C1, Ren_C1_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C1, tau_Ren_C1, 304.22, gate_on_C = [0], return_for_one = True)

    tau_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_tau'][0]
    number_of_pulses_Ren_C2 = mp['C' + str(carbon_nrs[1]) + '_Ren_N'][0]
    Ren_C2, Ren_C2_id = c13_gate_multiqubit(carbon_nrs, number_of_pulses_Ren_C2, tau_Ren_C1, 304.22, gate_on_C = [1], return_for_one = True)

    tau_z_C1 = phase_gate(carbon_nrs[0], np.pi, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
    RZ_C1, RZ_C1_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C1, 304.22, gate_on_C = [0], return_for_one = True, phase = np.pi)

    tau_z_C2 = phase_gate(carbon_nrs[0], np.pi, B_field=304.22,total_time = 0 ,return_gate = False,return_tau = True)
    RZ_C2, RZ_C2_id = c13_gate_multiqubit(carbon_nrs, 2, tau_z_C2, 304.22, gate_on_C = [1], return_for_one = True, phase = np.pi)

    # perform parity measurement
    if combination == 'XX':
        # seq = yel*Ren_C2*Ren_C1*yel
        seq_id = yel*Ren_C2_id*Ren_C1_id*yel
    elif combination == '-XX':
        # seq = yel*Ren_C2*Ren_C1*yel
        seq_id = yel*RZ_C1_id*Ren_C2_id*Ren_C1_id*yel
    elif combination == 'X-X':
        # seq = yel*Ren_C2*Ren_C1*yel
        seq_id = yel*RZ_C2_id*Ren_C2_id*Ren_C1_id*yel
    elif combination == '-X-X':
        # seq = yel*Ren_C2*Ren_C1*yel
        seq_id = yel*RZ_C1_id*RZ_C2_id*Ren_C2_id*Ren_C1_id*yel

    # rho_after =seq*rho_enc*seq.dag()
    rho_after_id =seq_id*rho_enc_id*seq_id.dag()

    # measure electron in 0 and renormalize
    el0 = qutip.tensor(rho0,Id,Id)
    rho_final_0_id = el0*rho_after_id*el0.dag()
    # rho_final_0 = el0*rho_after*el0.dag()

    norm_id = qutip.fidelity(rho0,rho_final_0_id.ptrace([0]))
    # norm = qutip.fidelity(rho0,rho_final_0.ptrace([0]))

    # rho_final_0 = 1/norm**2*rho_final_0.ptrace([1,2])
    rho_final_0_id = 1/(norm_id**2)*rho_final_0_id.ptrace([1,2])



    # measure electron in 1 and renormalize
    el1 = qutip.tensor(rho1,Id,Id)
    rho_final_1_id = el1*rho_after_id*el1.dag()
    # rho_final_1 = el1*rho_after*el1.dag()

    norm_id = qutip.fidelity(rho1,rho_final_1_id.ptrace([0]))
    # norm = qutip.fidelity(rho1,rho_final_1.ptrace([0]))

    # rho_final_1 = 1/norm**2*rho_final_1.ptrace([1,2])
    rho_final_1_id = 1/(norm_id**2)*rho_final_1_id.ptrace([1,2])


    # multi_qubit_pauli(rho_enc_id.ptrace([1,2]),carbon_nrs=carbon_nrs,do_plot=True, give_fid = False, alpha=None, beta=None,use_el=False,log_phase_corr=False)
    # multi_qubit_pauli(rho_after_id.ptrace([1,2]),carbon_nrs=carbon_nrs,do_plot=True, give_fid = False, alpha=None, beta=None,use_el=False,log_phase_corr=False)

    multi_qubit_pauli(rho_final_0_id,carbon_nrs=carbon_nrs,do_plot=True, give_fid = False, alpha=None, beta=None,use_el=False,log_phase_corr=False, title = 'Carbon nrs '+str(carbon_nrs)+ ', initial state '+initial_states+ ', parity msmt '+combination +', electron in 0')
    # ax.set_title('Carbon nrs '+str(Carbon_nrs)+ ', initial state '+initial_states+ ', parity msmt '+combination +', electron in 0')
    multi_qubit_pauli(rho_final_1_id,carbon_nrs=carbon_nrs,do_plot=True, give_fid = False, alpha=None, beta=None,use_el=False,log_phase_corr=False,title = 'Carbon nrs '+str(carbon_nrs)+ ', initial state '+initial_states+ ', parity msmt '+combination +', electron in 1')
    # ax.set_title('Carbon nrs '+str(Carbon_nrs)+ ', initial state '+initial_states+ ', parity msmt '+combination +', electron in 0')

def test_2qb_parity():
    for combination in ['XX','X-X','-XX','-X-X']:
        two_qb_entanglement_parity(carbon_nrs = [1,4], initial_states = 'ZZ' ,combination = combination)

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

def characterize_c13_DD_unit(carbon_nrs, B_field=304.22, tau_list = np.linspace(10,5000,500)):
    '''
    carbon_nrs is a list of carbons to calcuate for
    '''
    A_par = []
    A_perp = []

    for kk, carbon_nr in enumerate(carbon_nrs):
        par     = 2*np.pi*hf['C' + str(carbon_nr)]['par']
        perp    = 2*np.pi*hf['C' + str(carbon_nr)]['perp']

        A_par.append(par)
        A_perp.append(perp)

        print 'Carbon_nr = ' + str(carbon_nr)
        print 'Parallel hyperfine = '       + str(par/2./np.pi)
        print 'Perpendicular hyperfine = '  + str(perp/2./np.pi)

    characterize_DD_unit(A_par,A_perp,B_field=B_field,tau_list=tau_list)

def characterize_DD_unit(A_par = [2*np.pi*100e3], A_perp = [2*np.pi*30e3], B_field = 304.22, tau_list = np.linspace(10,5000,500), N=32):
    '''gives a full characterization of the rotation matrix
    for a single DD unit tau - X - 2tau - X - tau'''

    angle0 = np.zeros(len(tau_list))
    angle1 = np.zeros(len(tau_list))
    X_proj_0 = np.zeros(len(tau_list)); Y_proj_0 = np.zeros(len(tau_list)); Z_proj_0 = np.zeros(len(tau_list))
    X_proj_1 = np.zeros(len(tau_list)); Y_proj_1 = np.zeros(len(tau_list)); Z_proj_1 = np.zeros(len(tau_list))
    innerprod =np.zeros(len(tau_list)); pulses_for_pi2 = np.zeros(len(tau_list)); signal=np.zeros(len(tau_list))
    print B_field
    omega_Larmor = 2 * np.pi * B_field * 1.07e3

    X_proj_0_all = numpy.zeros((len(A_par),len(tau_list))); Y_proj_0_all = numpy.zeros((len(A_par),len(tau_list))); Z_proj_0_all = numpy.zeros((len(A_par),len(tau_list)))
    X_proj_1_all = numpy.zeros((len(A_par),len(tau_list))); Y_proj_1_all = numpy.zeros((len(A_par),len(tau_list))); Z_proj_1_all = numpy.zeros((len(A_par),len(tau_list)))

    innerprod_all       =   numpy.zeros((len(A_par),len(tau_list)))
    pulses_for_pi2_all  =   numpy.zeros((len(A_par),len(tau_list)))
    signal_all          =   numpy.zeros((len(A_par),len(tau_list)))

    angle0_all          =   numpy.zeros((len(A_par),len(tau_list)))

    for kk in range(len(A_par)):

        for i in range(len(tau_list)):
            tau = tau_list[i]*1e-9
            if i%10 == 0:
                print str(tau_list[i]) + ' out of ' + str(max(tau_list))
            V0, V1 = nuclear_rotation_matrix(tau, omega_Larmor, A_par[kk], A_perp[kk])
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

        X_proj_0_all[kk] =  X_proj_0
        Y_proj_0_all[kk] =  Y_proj_0
        Z_proj_0_all[kk] =  Z_proj_0
        X_proj_1_all[kk] =  X_proj_1
        Y_proj_1_all[kk] =  Y_proj_1
        Z_proj_1_all[kk] =  Z_proj_1

        innerprod_all[kk]       =   innerprod
        pulses_for_pi2_all[kk]  =   pulses_for_pi2
        signal_all[kk]          =   signal
        angle0_all[kk]          =   angle0


    #plots
    plt.close('all')

    f, ax = plt.subplots(3,3)

    for kk in range(len(A_par)):

        ax[0,0].plot(tau_list/1e3,X_proj_0_all[kk], '-', lw=1,label = 'data')
        ax[0,0].set_title('X projection ms=0'); ax[0,0].set_xlabel('tau (us)')

        ax[1,0].plot(tau_list/1e3,X_proj_1_all[kk], '-', lw=1,label = 'data')
        ax[1,0].set_title('X projection ms=1'); ax[1,0].set_xlabel('tau (us)')

        ax[2,0].plot(tau_list/1e3,(X_proj_0_all[kk]-X_proj_1_all[kk]), '-', lw=1,label = 'data')
        ax[2,0].set_title('X projection ms=0 - X projection ms=1'); ax[2,0].set_xlabel('tau (us)')

        ax[0,1].plot(tau_list/1e3,Y_proj_0_all[kk], '-', lw=1,label = 'data')
        ax[0,1].set_title('Y projection ms=0'); ax[0,1].set_xlabel('tau (us)')

        ax[1,1].plot(tau_list/1e3,Y_proj_1_all[kk], '-', lw=1,label = 'data')
        ax[1,1].set_title('Y projection ms=1'); ax[1,1].set_xlabel('tau (us)')

        ax[2,1].plot(tau_list/1e3,(Y_proj_0_all[kk]-Y_proj_1_all[kk]), '-', lw=1,label = 'data')
        ax[2,1].set_title('Y projection ms=0 - Y projection ms=1'); ax[2,1].set_xlabel('tau (us)')

        ax[0,2].plot(tau_list/1e3,Z_proj_0_all[kk], '-', lw=1,label = 'data')
        ax[0,2].set_title('Z projection ms=0'); ax[0,2].set_xlabel('tau (us)')

        ax[1,2].plot(tau_list/1e3,Z_proj_1_all[kk], '-', lw=1,label = 'data')
        ax[1,2].set_title('Z projection ms=1'); ax[1,2].set_xlabel('tau (us)')

        ax[2,2].plot(tau_list/1e3,(Z_proj_0_all[kk]-Z_proj_1_all[kk]), '-', lw=1,label = 'data')
        ax[2,2].set_title('Z projection ms=0 - Z projection ms=1'); ax[2,2].set_xlabel('tau (us)')

    f2, ax2 = plt.subplots(4,1)

    for kk in range(len(A_par)):

        ax2[0].plot(tau_list/1e3,innerprod_all[kk], '-', lw=1,label = 'data')
        ax2[0].set_title('axis innerporduct'); ax2[0].set_xlabel('tau (us)')
        ax2[1].plot(tau_list/1e3,angle0_all[kk]/np.pi, '-', lw=1,label = 'data')
        ax2[1].set_title('rotation_angle'); ax2[1].set_xlabel('tau (us)'); ax2[1].set_ylim(0,2)
        ax2[2].plot(tau_list/1e3,pulses_for_pi2_all[kk], '-', lw=1,label = 'data')
        ax2[2].set_title('nr of pulses'); ax2[2].set_xlabel('tau (us)')
        ax2[3].plot(tau_list/1e3,signal_all[kk], '-', lw=1,label = 'data')
        ax2[3].set_title('signal'); ax2[3].set_xlabel('tau (us)')

    plt.show()

print 'succes'
