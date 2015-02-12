### Tools for fitting QEC data, THT

import numpy as np
import os
from analysis.lib.tools import toolbox; reload(toolbox)
from analysis.lib.tools import plot; reload(plot)
from analysis.lib.fitting import fit, common, ramsey;reload(common); reload(fit)
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import qutip


def fit_init_fidelity_numerical(g_F1, g_A_par_1, g_A_perp_1,  g_F2, g_A_par_2, g_A_perp_2, g_det, g_t):
    '''Fit function for an electron Ramsey with/without C13 initialization
    g_F1        -  Fidelity of the strongly coupled C13
    g_A_par_1   -  Hyperfine of the strongly coupled C13
    g_A_perp_1  -  Hyperfine of the strongly coupled C13
    g_F2        -  Fidelity of the initialized C13
    g_A_par_2   -  Hyperfine of the initialized C13
    g_A_perp_2  -  Hyperfine of the initialized C13
    g_det       -  Detuning (carrier frequency)
    g_a         -  Offset, maybe fix to 0.5
    g_A         -  Amplitude maybe fix to 0.5
    g_t         -  Decay of Ramsey

    The function should take the complete data set (3 measurements) at once and fit it to a numberical model
    '''

    fitfunc_str = '''numerical solution'''

    ### Parameters
        ### Spin 1
    F1          = fit.Parameter(g_F1, 'F1')
    A_par_1     = fit.Parameter(g_A_par_1, 'A_par_1')
    A_perp_1    = fit.Parameter(g_A_perp_1, 'A_perp_1')

        ### Spin 2
    F2          = fit.Parameter(g_F2, 'F2')
    A_par_2     = fit.Parameter(g_A_par_2, 'A_par_2')
    A_perp_2    = fit.Parameter(g_A_perp_2, 'A_perp_2')

        ### Ramsey
    det = fit.Parameter(g_det, 'det')
    t   = fit.Parameter(g_t, 't')

    p0 = [F1, A_par_1, A_perp_1, F2, A_par_2, A_perp_2, det, t]

    def fitfunc(x):

        L = len(x)/3

        x_noinit   = x[0:L]
        x_up        = x[L:2*L]
        x_down      = x[2*L:3*L]

        ### Initial states
        rho_c1          = F1()*rho0 + (1-F1())*rho1
        rho_c2_up       = F2()*rho0 + (1-F2())*rho1
        rho_c2_down     = (1-F2())*rho0 + F2()*rho1
        rho_c2_noinit   = 0.5*rho0 + 0.5*rho1
        rho_up             = qutip.tensor(rhox,rho_c1,rho_c2_up)
        rho_down           = qutip.tensor(rhox,rho_c1,rho_c2_down)
        rho_noinit         = qutip.tensor(rhox,rho_c1,rho_c2_noinit)

        ### Hamiltonian
        H  =  2*np.pi*(det()*qutip.tensor(szz,Id,Id) + 430e3*(qutip.tensor(Id,sz,Id) + qutip.tensor(Id,Id,sz)) +
              A_par_1()*qutip.tensor(szz,sz,Id) + A_perp_1()*qutip.tensor(szz,sx,Id) +
              A_par_2()*qutip.tensor(szz,Id,sz) + A_perp_2()*qutip.tensor(szz,Id,sx) )

        Fid_up     = np.zeros(len(x_up))
        Fid_down   = np.zeros(len(x_down))
        Fid_noinit = np.zeros(len(x_noinit))

        for ii, tau in enumerate(x_up):
            expH = (-1j*H*tau).expm()

            ### State evolution
            rho_final_up        = expH * rho_up * expH.dag()
            rho_final_down      = expH * rho_down * expH.dag()
            rho_final_noinit    = expH * rho_noinit * expH.dag()

            ### Trace out the electron
            rho_el_final_up     = rho_final_up.ptrace(0)                  # Trace out the nuclear spins
            rho_el_final_down   = rho_final_down.ptrace(0)                  # Trace out the nuclear spins
            rho_el_final_noinit = rho_final_noinit.ptrace(0)                  # Trace out the nuclear spins

            ### Final fidelity
            Fid_up[ii]      = 0.5 + 0.5*(2*qutip.fidelity(rhox, rho_el_final_up)**2 -1)* np.exp(-(tau/t())**2)
            Fid_down[ii]    = 0.5 + 0.5*(2*qutip.fidelity(rhox, rho_el_final_down)**2 -1)* np.exp(-(tau/t())**2)
            Fid_noinit[ii]  = 0.5 + 0.5*(2*qutip.fidelity(rhox, rho_el_final_noinit)**2 -1)* np.exp(-(tau/t())**2)

        Fid = np.r_[Fid_noinit, Fid_up, Fid_down]
        return Fid

    return p0, fitfunc, fitfunc_str

### Basic states and gates
Id = qutip.qeye(2)
sx = qutip.sigmax()/2
sz = qutip.sigmaz()/2
szz = (qutip.sigmaz()-qutip.qeye(2))/2

ket0 = qutip.basis(2,0)
bra0 = qutip.basis(2,0).dag()
ket1 = qutip.basis(2,1)
bra1 = qutip.basis(2,1).dag()
rho0 = ket0*bra0
rho1 = ket1*bra1
ketx = 1/np.sqrt(2)*(qutip.basis(2,0)+qutip.basis(2,1))
brax = 1/np.sqrt(2)*(qutip.basis(2,0).dag()+qutip.basis(2,1).dag())
rhox =ketx*brax

### Fake data
data_F1          = 0.5
data_A_par_1     = 180e3
data_A_perp_1    = 0
data_F2          = 0.9
data_A_par_2     = 21.2e3
data_A_perp_2    = 43e3
data_det         = 500e3
data_t           = 5e-6

p0, fitfunc, fitfunc_str = fit_init_fidelity_numerical(data_F1, data_A_par_1, data_A_perp_1,
                                                       data_F2, data_A_par_2, data_A_perp_2,
                                                       data_det, data_t)
xdata = np.r_[np.linspace(2e-6,8e-6,20),np.linspace(2e-6,8e-6,20),np.linspace(2e-6,8e-6,20)]
ydata = fitfunc(xdata)

### plot the fake data
plt.figure(1)
L = len(xdata)/3
plt.plot(xdata[0:L],ydata[0:L],'bo')
plt.plot(xdata[L:2*L],ydata[L:2*L],'ro')
plt.plot(xdata[2*L:3*L],ydata[2*L:3*L],'go')

### Fitting
guess_F1          = 0.4
guess_A_par_1     = 180e3
guess_A_perp_1    = 0
guess_F2          = 1
guess_A_par_2     = 21.2e3
guess_A_perp_2    = 43e3
guess_det         = 400e3
guess_t           = 7e-6

p0, fitfunc, fitfunc_str = fit_init_fidelity_numerical(guess_F1, guess_A_par_1, guess_A_perp_1,
                                                       guess_F2, guess_A_par_2, guess_A_perp_2,
                                                       guess_det, guess_t)

print 'fitting, this might take a while'
fit_result = fit.fit1d(xdata, ydata, fit_init_fidelity_numerical,
            guess_F1, guess_A_par_1, guess_A_perp_1,
            guess_F2, guess_A_par_2, guess_A_perp_2,
            guess_det, guess_t,
            fixed=[1,2,4,5],
            do_print=True, ret=True)

print fit_result['params_dict']
### plot result

p02, fitfunc2, fitfunc_str2 = fit_init_fidelity_numerical(
                            fit_result['params_dict']['F1'],
                            guess_A_par_1,#fit_result['params_dict']['A_par_1'],
                            guess_A_perp_1,#fit_result['params_dict']['A_perp_1'],
                            fit_result['params_dict']['F2'],
                            guess_A_par_2,#fit_result['params_dict']['A_par_2'],
                            guess_A_perp_2,#fit_result['params_dict']['A_perp_2'],
                            fit_result['params_dict']['det'],
                            fit_result['params_dict']['t'])

x_temp = np.r_[np.linspace(xdata[0],xdata[-1],50), np.linspace(xdata[0],xdata[-1],50), np.linspace(xdata[0],xdata[-1],50)]
y_temp = fitfunc2(x_temp)
L = len(x_temp)/3

plt.plot(x_temp[0:L], y_temp[0:L],'b')
plt.plot(x_temp[L:2*L], y_temp[L:2*L],'r')
plt.plot(x_temp[2*L:3*L], y_temp[2*L:3*L],'g')

plt.show()



### Fitfunction

# def fit_QEC(g_O, g_A, g_p):
#     '''Fit function for QEC process fidelity data
#     g_O -  Offset, given by the fidelity of the state that is insensitive to errors
#     g_A -  Amplitude, g_iven by the fidelity of the states that are sensitive
#     g_p -  Avegage probabililty to correct single qubit errors
#     '''

#     fitfunc_str = '''test'''

#     O   = fit.Parameter(g_O , 'O')
#     A   = fit.Parameter(g_A, 'A')
#     p   = fit.Parameter(g_p, 'p')

#     p0 = [O, A, p]

#     def fitfunc(x):
#         '''test'''
#         return (O() + A()*(  1-3*x+3*x**2-2*x**3 + 3*(2*p()-1)*(x-3*x**2+2*x**3)))

#     return p0, fitfunc, fitfunc_str

# ##################################
# ### Make fake data and plot it ###
# ##################################
# data_O = 0.4
# data_A = 0.4
# data_p = 0.85

# p0, fitfunc, fitfunc_str = fit_QEC(data_O, data_A, data_p)

# xdata = np.linspace(0,1,20)
# ydata = fitfunc(xdata)

# plt.figure()
# plt.plot(xdata,ydata,'bo')

# ###########
# ### Fit ###
# ###########

# guess_O = 0.5
# guess_A = 0.5
# guess_p = 1

# p0, fitfunc, fitfunc_str = fit_QEC(guess_O, guess_A, guess_p)

# fit_result = fit.fit1d(xdata, ydata, fit_QEC,
#         guess_O, guess_A, guess_p,
#         fixed=[],
#         do_print=True, ret=True)

# ### show guess ###
# if 1:
#     x_temp = np.linspace(xdata[0],xdata[-1],200)
#     y_temp = fitfunc(x_temp)

#     plt.plot(x_temp, y_temp, 'b', lw=2)

# p02, fitfunc2, fitfunc_str2 = fit_QEC(fit_result['params_dict']['O'], fit_result['params_dict']['A'], fit_result['params_dict']['p'])

# plt.plot(x_temp, fitfunc2(x_temp),'r')


# plt.show()
