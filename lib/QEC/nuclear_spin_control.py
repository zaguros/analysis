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
    return ket0, bra0, ket1, bra1, rho0, rho1, rhom

### create a set of usefull simple states and gates
Id, sx, sy, sz = pauli()                                # Electron spin operators
Id, Ix, Iy, Iz = pauli()                                # Nuclear spin operators
X,Y,Z,x,y,z,mX,mY,mZ,mx,my,mz = basic_spin_rotations()  # Basic gates
ket0, bra0, ket1, bra1, rho0, rho1, rhom = basic_spin_states() # Basic states

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

def waittime(carbon_nr, time, B_field=304.22,return_indiv = False):
    '''calculates the evolution matrices for a single
    Carbon spin, electron is always qubit 1'''

    omega_Larmor = 2 * np.pi * B_field * 1.07e3
    A_par = 2 * np.pi * hf['C' + str(carbon_nr)]['par']
    A_perp = 2 * np.pi *hf['C' + str(carbon_nr)]['perp']

    H0 = omega_Larmor * Iz
    H1 = (A_par+omega_Larmor)*Iz + A_perp*Ix

    expH0 = (-1j*H0*tau).expm();    expH1 = (-1j*H1*tau).expm()
    Utot = qutip.tensor(rho0,expH0) + qutip.tensor(rho1,expH1)
    print_matrix(Utot)

    if return_indiv == False:
        return Utot
    if return_indiv == True:
        return expH0, expH1


def xn_gate(carbon_nr, phase):
    pass

def c13_gate_multiqubit(carbon_nr, number_of_pulses, tau, B_field):
    '''calculates the evolution matrices a multiqubit space,
    the electron is always qubit 1'''
    pass



###################
### Pauli Sets ###
###################

def single_qubit_pauli(rho, do_plot = False):
    ii=-0.5
    pauli_set = []
    ii_list = []
    xticks_list = ['II','X','Y','Z']

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
    return pauli_set, ii_list, xticks_list
def multi_qubit_pauli(rho,do_plot=False):
    
    no_qubits = np.size(qutip.dims(rho))/2
    ii=-0.5
    pauli_set = []
    ii_list = []

    xticks_list = ['X','Y','Z']
    final_x_tick_list = []
    oper_list = [2*sx,2*sy,2*sz]
    final_oper_list =[]


    if no_qubits ==2:
        for ff in xticks_list:
            final_x_tick_list.append('I'+ff)
        for ff in xticks_list:
            final_x_tick_list.append(ff+'I')
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(Id,ff))
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(ff,Id))
        for jj in oper_list:
            for kk in oper_list:
                final_oper_list.append(qutip.tensor(jj,kk))
        for jj in xticks_list:
            for kk in xticks_list:
                final_x_tick_list.append(jj+kk)

    if no_qubits ==3:
        
        for ff in xticks_list:
            final_x_tick_list.append('I'+'I'+ff)
        for ff in xticks_list:
            final_x_tick_list.append('I'+ff+'I')
        for ff in xticks_list:
            final_x_tick_list.append(ff+'I'+'I')
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(Id,Id,ff))
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(Id,ff,Id))
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(ff,Id,Id))
        for jj in oper_list:
            for kk in oper_list:
                for hh in oper_list:
                    final_oper_list.append(qutip.tensor(jj,kk,hh))
        for jj in xticks_list:
            for kk in xticks_list:
                for hh in xticks_list:
                    final_x_tick_list.append(jj+kk+hh)
    
    if no_qubits ==4:
        for ff in xticks_list:
            final_x_tick_list.append('I'+'I'+'I'+ff)
        for ff in xticks_list:
            final_x_tick_list.append('I'+'I'+ff+'I')
        for ff in xticks_list:
            final_x_tick_list.append('I'+ff+'I'+'I')
        for ff in xticks_list:
            final_x_tick_list.append(ff+'I'+'I'+'I')
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(Id,Id,Id,ff))
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(Id,Id,ff,Id))
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(Id,ff,Id,Id))
        for ff in oper_list:
            final_oper_list.append(qutip.tensor(ff,Id,Id,Id))
        for jj in oper_list:
            for kk in oper_list:
                for hh in oper_list:
                    for mm in oper_list:
                        final_oper_list.append(qutip.tensor(jj,kk,hh,mm))
        for jj in xticks_list:
            for kk in xticks_list:
                for hh in xticks_list:
                    for mm in xticks_list:
                        final_x_tick_list.append(jj+kk+hh+mm)
    
    for oper in final_oper_list:                   
        pauli_set.append(qutip.expect(oper,rho))
        ii_list.append(ii)
        ii = ii+1
    

    if do_plot == True:
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.bar(ii_list, pauli_set, width=1)
        plt.xticks(np.arange(0, len(final_oper_list)-1, 1.0))
        ax.set_xticklabels(final_x_tick_list)
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

def nuclear_ramsey_no_init_no_DD(carbon_nr, tau_wait, N_wait_list, B_field=304.22):
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


######################
### Initialization ### NOTE: FOR MULTIPLE C13 SPINS THE PHASE GATES SHOULD BE ADDED!!
######################

def nuclear_init_single(carbon_nr,do_plot = False):
    '''function that returns a gate sequence for a single nuclear spin initialization
    seq = y - Ren - x - Rz - Ren  note: Z-gate not yet calculated in right way'''

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

    seq = Ren*Rz*xel*Ren*yel
    seq_id = Ren_id*Rz*xel*Ren_id*yel

    rho_final = seq*rho*seq.dag()  
    rho_final_id = seq_id*rho*seq_id.dag()  


    rho_nucl = rho_final.ptrace(1)
    rho_nucl_id = rho_final_id.ptrace(1)

    if do_plot == True:
        fig = plt.figure()
        ax = plt.subplot(111)

        for rho_n in [rho_nucl,rho_nucl_id]:
            pauli_set, ii_list, x_ticks_list = single_qubit_pauli(rho_n)
            ax.bar(ii_list, np.real(pauli_set), width=1,alpha = 0.5)

        plt.xticks(np.arange(0, 4, 1.0))
        ax.set_xticklabels(x_ticks_list)

        print 'Fidelity to ideal state:'
        print qutip.fidelity(rho_nucl,rho_nucl_id)
    return rho_nucl, rho_nucl_id


def three_spin_encoding(carbon_nrs = [1,1,1],alpha=1/np.sqrt(2),beta=1/np.sqrt(2),do_plot=True):
    ''' encodes your three chosen C13 spins in the state alpha(xxx)+beta(-x-x-x)
    note: it depends on the direction of the Ren gate for each C13 spin if it will be + or -x, see documentation'''

    # define density matrices
    rho_C1, rho_C1_id = nuclear_init_single(carbon_nrs[0])
    rho_C2, rho_C2_id = nuclear_init_single(carbon_nrs[1])
    rho_C3, rho_C3_id = nuclear_init_single(carbon_nrs[2])

    psi_el = alpha*ket0-1j*beta*ket1
    rho_el = psi_el*psi_el.dag()

    rho = qutip.tensor(rho_el,rho_C1,rho_C2,rho_C3)
    rho_id = qutip.tensor(rho_el,rho_C1_id,rho_C2_id,rho_C3_id)

    rho_C = qutip.tensor(rho_C1,rho_C2,rho_C3)
    rho_C_id = qutip.tensor(rho_C1_id,rho_C2_id,rho_C3_id)


    #define gates
    xel = qutip.tensor(x,Id,Id,Id)
    yel = qutip.tensor(x,Id,Id,Id)

    U0, U1 = nuclear_Ren_matrix(carbon_nrs[0],B_field=304.22)
    U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
    U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)
    Ren_C1 = qutip.tensor(rho0,U0,Id,Id)+qutip.tensor(rho1,U1,Id,Id)
    Ren_C1_id = qutip.tensor(rho0,U0id,Id,Id)+qutip.tensor(rho1,U1id,Id,Id)
    Rz_C1 = qutip.tensor(Id,z,Id,Id)

    U0, U1 = nuclear_Ren_matrix(carbon_nrs[1],B_field=304.22)
    U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
    U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)
    Ren_C2 = qutip.tensor(rho0,Id,U0,Id)+qutip.tensor(rho1,Id,U1,Id)
    Ren_C2_id = qutip.tensor(rho0,Id,U0id,Id)+qutip.tensor(rho1,Id,U1id,Id)
    Rz_C2 = qutip.tensor(Id,Id,z,Id)

    U0, U1 = nuclear_Ren_matrix(carbon_nrs[2],B_field=304.22)
    U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
    U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)
    Ren_C3 = qutip.tensor(rho0,Id,Id,U0)+qutip.tensor(rho1,Id,Id,U1)
    Ren_C3_id = qutip.tensor(rho0,Id,Id,U0id)+qutip.tensor(rho1,Id,Id,U1id)
    Rz_C3 = qutip.tensor(Id,Id,Id,z)

    # define and apply full sequence
    seq = xel*Rz_C3*Rz_C2*Rz_C1*Ren_C3*Ren_C2*Ren_C1
    seq_id = xel*Rz_C3*Rz_C2*Rz_C1*Ren_C3_id*Ren_C2_id*Ren_C1_id
    rho_after = seq*rho*seq.dag()
    rho_after_id = seq_id*rho_id*seq_id.dag()

    #measure electron to be in zero
    el0 = qutip.tensor(rho0,Id,Id,Id)
    rho_final_id = el0*rho_after_id*el0.dag()
    rho_final = el0*rho_after*el0.dag()

    norm_id = qutip.fidelity(rho0,rho_final_id.ptrace([0]))
    norm = qutip.fidelity(rho0,rho_final.ptrace([0]))

    rho_final = 1/norm**2*rho_final.ptrace([1,2,3])
    rho_final_id = 1/norm_id**2*rho_final_id.ptrace([1,2,3])

    if do_plot == True:
        fig = plt.figure()
        ax = plt.subplot(111)
        for rho_n in [rho_final,rho_final_id]:
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

    psi_el = alpha*ket0-1j*beta*ket1
    rho_el = psi_el*psi_el.dag()
    print rho_el
    rho = qutip.tensor(rho_el,rho0,rho0,rho0)

    #define gates
    xel = qutip.tensor(x,Id,Id,Id)
    yel = qutip.tensor(x,Id,Id,Id)

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

    rho_final = 1/norm**2*rho_final.ptrace([1,2,3])


    fig = plt.figure(figsize=(100,5))
    ax = plt.subplot(111)

    pauli_set, ii_list, x_ticks_list = multi_qubit_pauli(rho_final)
    ax.bar(ii_list, np.real(pauli_set), width=1,alpha = 0.5)

    plt.xticks(np.arange(0, len(x_ticks_list), 1.0))
    ax.set_xticklabels(x_ticks_list)
    # ax.set_title('entangled state = 1/sqrt(2)*(' + str(alpha*np.sqrt(2)) +' |xxx>' + str(beta*np.sqrt(2)) + ' |-x-x-x>)')
    ax.set_title('entangled state = (' + str(alpha) +' |xxx>' + str(beta) + ' |-x-x-x>)')
    return rho_final

def two_spin_encoding(carbon_nrs = [1,1],alpha=1/np.sqrt(2),beta=1/np.sqrt(2)):
    ''' encodes your three chosen C13 spins in the state alpha(xxx)+beta(-x-x-x)
    note: it depends on the direction of the Ren gate for each C13 spin if it will be + or -x, see documentation'''

    # define density matrices
    rho_C1, rho_C1_id = nuclear_init_single(carbon_nrs[0])
    rho_C2, rho_C2_id = nuclear_init_single(carbon_nrs[1])

    psi_el = alpha*ket0-1j*beta*ket1
    rho_el = psi_el*psi_el.dag()

    rho = qutip.tensor(rho_el,rho_C1,rho_C2)
    rho_id = qutip.tensor(rho_el,rho_C1_id,rho_C2_id)

    rho_C = qutip.tensor(rho_C1,rho_C2)
    rho_C_id = qutip.tensor(rho_C1_id,rho_C2_id)


    #define gates
    xel = qutip.tensor(x,Id,Id)
    yel = qutip.tensor(x,Id,Id)

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

    fig = plt.figure()
    ax = plt.subplot(111)
    for rho_n in [rho_final,rho_final_id]:
        pauli_set, ii_list, x_ticks_list = multi_qubit_pauli(rho_n)
        ax.bar(ii_list, np.real(pauli_set), width=1,alpha = 0.5)

    plt.xticks(np.arange(0, len(x_ticks_list), 1.0))
    ax.set_xticklabels(x_ticks_list)

#######################
### Error detection ###
#######################

def parity_msmt(qubits=[0,1],carbon_nrs = [1,1,1],alpha=1/np.sqrt(2),beta=1/np.sqrt(2), error_list = ['Q1','Q2']):
    '''implements an error or not on the two qubits that are selected, then parity is measured
    no partial error yet'''

    rho_enc_start, rho_enc_start_id = three_spin_encoding(carbon_nrs=carbon_nrs,alpha=alpha,beta=beta,do_plot=False)
    rho_enc = qutip.tensor(rho0,rho_enc_start)
    rho_enc_id = qutip.tensor(rho0,rho_enc_start_id)

    #define gates
    xel = qutip.tensor(x,Id,Id,Id)
    Xel = qutip.tensor(X,Id,Id,Id)
    yel = qutip.tensor(x,Id,Id,Id)

    U0, U1 = nuclear_Ren_matrix(carbon_nrs[qubits[0]],B_field=304.22)
    U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
    U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)

    C1_0 = [rho0,Id,Id,Id]
    C1_0[qubits[0]+1] = U0
    C1_0_id = [rho0,Id,Id,Id]
    C1_0_id[qubits[0]+1] = U0id
    C1_1 = [rho1,Id,Id,Id]
    C1_1[qubits[0]+1] = U1
    C1_1_id = [rho1,Id,Id,Id]
    C1_1_id[qubits[0]+1] = U1id

    Ren_C1 = qutip.tensor(C1_0)+qutip.tensor(C1_1)
    Ren_C1_id = qutip.tensor(C1_0_id)+qutip.tensor(C1_1_id)
    ZC1 = [Id,Id,Id,Id]
    ZC1[qubits[0]+1] = z
    Rz_C1 = qutip.tensor(ZC1)

    U0, U1 = nuclear_Ren_matrix(carbon_nrs[qubits[1]],B_field=304.22)
    U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
    U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)

    C2_0 = [rho0,Id,Id,Id]
    C2_0[qubits[1]+1] = U0
    C2_0_id = [rho0,Id,Id,Id]
    C2_0_id[qubits[1]+1] = U0id
    C2_1 = [rho1,Id,Id,Id]
    C2_1[qubits[1]+1] = U1
    C2_1_id = [rho1,Id,Id,Id]
    C2_1_id[qubits[1]+1] = U1id

    Ren_C2 = qutip.tensor(C2_0)+qutip.tensor(C2_1)
    Ren_C2_id = qutip.tensor(C2_0_id)+qutip.tensor(C2_1_id)
    ZC2 = [Id,Id,Id,Id]
    ZC2[qubits[1]+1] = z
    Rz_C2 = qutip.tensor(ZC2)

    # implement error

    for error in error_list:
        if error =='Q1':
            rho_enc = (Rz_C1*Rz_C1)*rho_enc*(Rz_C1*Rz_C1).dag()
            rho_enc_id = (Rz_C1*Rz_C1)*rho_enc_id*(Rz_C1*Rz_C1).dag()
        if error =='Q2':
            rho_enc = (Rz_C2*Rz_C2)*rho_enc*(Rz_C2*Rz_C2).dag()
            rho_enc_id = (Rz_C2*Rz_C2)*rho_enc_id*(Rz_C2*Rz_C2).dag()

    # detect error
    seq = yel*Ren_C2*Ren_C1*yel
    seq_id = yel*Ren_C2_id*Ren_C1_id*yel
    rho_after =seq*rho_enc*seq.dag()
    rho_after_id =seq_id*rho_enc_id*seq_id.dag()
    # measure electron state

    rho_el_after = rho_after.ptrace(0)
    rho_el_after_id = rho_after_id.ptrace(0)

    print 'Fidelity to ideal el after state:'
    print qutip.fidelity(rho_el_after,rho_el_after_id)

    print 'Fidelity to opposite el after state:'
    print qutip.fidelity(rho_el_after,X*rho_el_after_id*X.dag())


    fig = plt.figure()
    ax = plt.subplot(111)
    for rho_el in [rho_el_after,rho_el_after_id]:
        pauli_set, ii_list, x_ticks_list = single_qubit_pauli(rho_el, do_plot = False)
        ax.bar(ii_list, np.real(pauli_set), width=1,alpha = 0.5)

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

    fig = plt.figure()
    ax = plt.subplot(111)
    rho_enc_exp = [rho_enc_start,rho_enc_final]
    colors = ["red", "blue"]
    labels = ['exp before parity msmst','exp after parity msmt']
    for ii in [0,1]:
        pauli_set, ii_list, x_ticks_list = multi_qubit_pauli(rho_enc_exp[ii], do_plot = False)
        ax.bar(ii_list, np.real(pauli_set), color = colors[ii], width=1,alpha = 0.5,label = labels[ii])

    plt.xticks(np.arange(0, len(x_ticks_list), 1.0))
    ax.set_xticklabels(x_ticks_list) 
    ax.legend() 

    fig = plt.figure()
    ax = plt.subplot(111)
    rho_enc =[rho_enc_start_id,rho_enc_final_id]
    colors = ["red", "blue"]
    labels = ['id before parity msmst','id after parity msmt']
    for ii in [0,1]:
        pauli_set, ii_list, x_ticks_list = multi_qubit_pauli(rho_enc[ii], do_plot = False)
        ax.bar(ii_list, np.real(pauli_set), color = colors[ii],width=1,alpha = 0.5,label = labels[ii])

    plt.xticks(np.arange(0, len(x_ticks_list), 1.0))
    ax.set_xticklabels(x_ticks_list)  
    ax.legend() 

def three_qubit_msmt_via_el(carbon_nrs = [1,1,1],alpha=1/np.sqrt(2),beta=1/np.sqrt(2),state ='+++' ,meas = 'XXX'):

    # rho_enc_start, rho_enc_start_id = three_spin_encoding(carbon_nrs=carbon_nrs,alpha=alpha,beta=beta,do_plot=False)
    rho_enc_start_id = check_entangled_state(alpha = alpha, beta = beta, state = state) # for now
    rho_enc_id = qutip.tensor(rho0,rho_enc_start_id)

    #define gates
    xel = qutip.tensor(x,Id,Id,Id)
    yel = qutip.tensor(x,Id,Id,Id)

    U0, U1 = nuclear_Ren_matrix(carbon_nrs[0],B_field=304.22)
    U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
    U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)
    Ren_C1 = qutip.tensor(rho0,U0,Id,Id)+qutip.tensor(rho1,U1,Id,Id)
    Ren_C1_id = qutip.tensor(rho0,U0id,Id,Id)+qutip.tensor(rho1,U1id,Id,Id)
    Rz_C1 = qutip.tensor(Id,z,Id,Id)

    U0, U1 = nuclear_Ren_matrix(carbon_nrs[1],B_field=304.22)
    U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
    U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)
    Ren_C2 = qutip.tensor(rho0,Id,U0,Id)+qutip.tensor(rho1,Id,U1,Id)
    Ren_C2_id = qutip.tensor(rho0,Id,U0id,Id)+qutip.tensor(rho1,Id,U1id,Id)
    Rz_C2 = qutip.tensor(Id,Id,z,Id)

    U0, U1 = nuclear_Ren_matrix(carbon_nrs[2],B_field=304.22)
    U0id = np.round(U0.full()*np.sqrt(2))/np.sqrt(2)
    U1id = np.round(U1.full()*np.sqrt(2))/np.sqrt(2)
    U0id = qutip.Qobj(U0id)
    U1id = qutip.Qobj(U1id)
    Ren_C3 = qutip.tensor(rho0,Id,Id,U0)+qutip.tensor(rho1,Id,Id,U1)
    Ren_C3_id = qutip.tensor(rho0,Id,Id,U0id)+qutip.tensor(rho1,Id,Id,U1id)
    Rz_C3 = qutip.tensor(Id,Id,Id,z)
    # project to el
    seq_id = yel*Ren_C2_id*Ren_C1_id*yel
    rho_after_id =seq_id*rho_enc_id*seq_id.dag()

#project electron
    rho_el_after_id = rho_after_id.ptrace(0)
    el0 = qutip.tensor(rho_el_after_id,Id,Id,Id)
    rho_final_id = el0*rho_after_id*el0.dag()

    norm_id = qutip.fidelity(rho_el_after_id,rho_final_id.ptrace([0]))

    rho_enc_final_id = 1/norm_id**2*rho_final_id.ptrace([1,2,3])

    fig = plt.figure()
    ax = plt.subplot(111)
    pauli_set, ii_list, x_ticks_list = single_qubit_pauli(rho_el_after_id, do_plot = False)
    ax.bar(ii_list, np.real(pauli_set), width=1,alpha = 0.5)

    plt.xticks(np.arange(0, len(x_ticks_list), 1.0))
    ax.set_xticklabels(x_ticks_list)
    ax.set_title('electron state')

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

    characterize_DD_unit(A_par,A_perp,B_field=B_field,tau_list=tau_list)

def characterize_DD_unit(A_par = 2*np.pi*100e3, A_perp = 2*np.pi*30e3, B_field = 304.22, tau_list = np.linspace(10,5000,500), N=32):
    '''gives a full characterization of the rotation matrix
    for a single DD unit tau - X - 2tau - X - tau'''

    angle0 = np.zeros(len(tau_list))
    angle1 = np.zeros(len(tau_list))
    X_proj_0 = np.zeros(len(tau_list)); Y_proj_0 = np.zeros(len(tau_list)); Z_proj_0 = np.zeros(len(tau_list))
    X_proj_1 = np.zeros(len(tau_list)); Y_proj_1 = np.zeros(len(tau_list)); Z_proj_1 = np.zeros(len(tau_list))
    innerprod =np.zeros(len(tau_list)); pulses_for_pi2 = np.zeros(len(tau_list)); signal=np.zeros(len(tau_list))
    print B_field
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
