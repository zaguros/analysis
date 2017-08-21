import numpy as np
import os
from analysis.lib.tools import toolbox, plot;
from analysis.lib.m2.ssro import mbi, sequence, pqsequence
from matplotlib import pyplot as plt
# import matplotlib as mpl
from analysis.lib.fitting import fit, common
import copy as cp

reload(fit)
reload(mbi)
reload(common)
reload(toolbox)

CR_after_check = False  # global variable that let's us post select whether or not the NV was ionized


def carbon_phase_noise_decay_curve(A0, N, omega, sigma_t, phi):
    return A0*np.cos(phi)*np.exp(-0.5*(sigma_t**2)*(omega**2)*N)


def fit_carbon_phase_noise_decay_curve(g_a, g_A, g_x0, g_T, g_n, g_f, g_phi):
    fitfunc_str = 'carbon phase noise decay curve'

    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    x0 = fit.Parameter(g_x0, 'x0')
    T = fit.Parameter(g_T, 'T')
    n = fit.Parameter(g_n, 'n')
    f = fit.Parameter(g_f, 'f')
    phi = fit.Parameter(g_phi, 'phi')

    A0_A = fit.Parameter(g_A0_A, 'A0_A')


    p0 = [a, A, x0, f, phi]
    def fitfunc(x):
        return a() + A() * np.cos(2*np.pi*( f()*x + phi()/360.))

    return p0, fitfunc, fitfunc_str


def plot_noisey_carbon_correlations(
        Ns, sigma_t,
        A0_A, A0_B,
        d_omega_A, d_omega_B,
        detuning_A, detuning_B,
        phi0_A=0.0, phi0_B=0.0,
        **kwargs
    ):

    phi_A = (phi0_A + Ns*detuning_A) / 180.0 * np.pi
    phi_B = (phi0_B + Ns*detuning_B) / 180.0 * np.pi

    p_term = carbon_phase_noise_decay_curve(A0_A, Ns, d_omega_A + d_omega_B, sigma_t, phi_A+phi_B)
    m_term = carbon_phase_noise_decay_curve(A0_B, Ns, d_omega_A - d_omega_B, sigma_t, phi_A-phi_B)

    # plt.plot(Ns, p_term, label="+ term")
    # plt.plot(Ns, m_term, label="- term")

    plt.plot(Ns, 0.5*p_term + 0.5*m_term, **kwargs)