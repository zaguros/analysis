import numpy as np

# own modules
import fit
#=======
#import analysis.lib.fitting.fit as fit


### common fitfunctions
def fit_cos(g_f, g_a, g_A, g_phi, *arg):
    fitfunc_str = 'A * cos(2pi * (f*x + phi/360) ) + a'

    f = fit.Parameter(g_f, 'f')
    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    phi = fit.Parameter(g_phi, 'phi')

    p0 = [f, a, A,phi] #Note: If you do not want to use a fit argument set fixed when using in fit1d 

    def fitfunc(x):
        return a() + A() * np.cos(2*np.pi*( f()*x + phi()/360.))

    return p0, fitfunc, fitfunc_str


def fit_decaying_cos(g_f, g_a, g_A, g_phi,g_t, *arg):
    fitfunc_str = 'A *exp(-x/t) cos(2pi * (f*x + phi/360) ) + a'

    f = fit.Parameter(g_f, 'f')
    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    phi = fit.Parameter(g_phi, 'phi')
    t   = fit.Parameter(g_t, 't')
    print 'guessed frequency is '+str(g_f)
    p0 = [f, a, A,phi,t]

    def fitfunc(x):
        return a() + A()*np.exp(-x/t()) * np.cos(2*np.pi*( f()*x + phi()/360.))

    return p0, fitfunc, fitfunc_str

def fit_double_decaying_cos(g_f1, g_A1, g_phi1, g_t1, g_f2, g_A2, g_phi2, g_t2, *arg):
    ''' quite a specific function, for electron nuclear control, maybe place somewhere else '''
    fitfunc_str = '''(A1 *exp(-x/t1) cos(2pi * (f1*x + phi1/360) ) + a1)*
                     (A2 *exp(-x/t2) cos(2pi * (f2*x + phi2/360) ) + a2)/2+1/2'''

    f1 = fit.Parameter(g_f1, 'f1')
    #a1 = fit.Parameter(g_a1, 'a1')
    A1 = fit.Parameter(g_A1, 'A1')
    phi1 = fit.Parameter(g_phi1, 'phi1')
    t1   = fit.Parameter(g_t1, 't1')

    f2 = fit.Parameter(g_f2, 'f2')
    #a2 = fit.Parameter(g_a2, 'a2')
    A2 = fit.Parameter(g_A2, 'A2')
    phi2 = fit.Parameter(g_phi2, 'phi2')
    t2   = fit.Parameter(g_t2, 't2')

    #p0 = [f1, a1, A1, phi1, t1, f2, a2, A2, phi2, t2]
    p0 = [f1, A1, phi1, t1, f2, A2, phi2, t2]

    def fitfunc(x):
        return (1 - A1() + A1()*np.exp(-x/t1()) * np.cos(2*np.pi*( f1()*x + phi1()/360.)))*(1-A2() + A2()*np.exp(-x/t2()) * np.cos(2*np.pi*( f2()*x + phi2()/360.)))/2+0.5

    return p0, fitfunc, fitfunc_str

def fit_exp_decay_with_offset(g_a, g_A, g_tau, *arg):
    """
    fitfunction for an exponential decay,
        y(x) = A * exp(-x/tau) + a

    Initial guesses (in this order):
        g_a : offset
        g_A : initial Amplitude
        g_tau : decay constant

    """
    fitfunc_str = 'A * exp(-x/tau) + a'

    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    tau = fit.Parameter(g_tau, 'tau')
    p0 = [a, A, tau]

    def fitfunc(x):
        return a() + A() * np.exp(-x/tau())

    return p0, fitfunc, fitfunc_str

def fit_double_exp_decay_with_offset(g_a, g_A, g_tau, g_A2, g_tau2, *arg):
    """
    fitfunction for an exponential decay,
        y(x) = A * exp(-x/tau)+ A2 * exp(-x/tau2) + a

    Initial guesses (in this order):
        g_a : offset
        g_A : initial Amplitude
        g_tau : decay constant
        g_A2 : initial Amplitude 2
        g_tau2 : decay constant 2
    """
    fitfunc_str = 'A * exp(-x/tau)+ A2 * exp(-x/tau2) + a'

    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    tau = fit.Parameter(g_tau, 'tau')
    A2 = fit.Parameter(g_A2, 'A2')
    tau2 = fit.Parameter(g_tau2, 'tau2')
    p0 = [a, A, tau, A2, tau2]

    def fitfunc(x):
        return a() + A() * np.exp(-x/tau()) + A2() * np.exp(-x/tau2())

    return p0, fitfunc, fitfunc_str

def fit_exp_decay_shifted_with_offset(g_a, g_A, g_tau, g_x0, *arg):
    """
    fitfunction for an exponential decay,
        y(x) = A * exp(-(x-x0)/tau) + a

    Initial guesses (in this order):
        g_a : offset
        g_A : initial Amplitude
        g_tau : decay constant
        g_x0 : x offset
    """

    fitfunc_str = 'A * exp(-(x-x0)/tau) + a'

    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    tau = fit.Parameter(g_tau, 'tau')

    x0 = fit.Parameter(g_x0, 'x0')
    p0 = [a, A, tau, x0]

    def fitfunc(x):
        return a() + A() * np.exp(-(x-x0())/tau())

    return p0, fitfunc, fitfunc_str

def fit_saturation(g_A, g_xsat, *arg):
    """
    fitfunction for a saturation (e.g., the NV PL)
        y(x) = A * x / (x + xsat)

    I.g.:
        g_A : maximum signal (at x=infinity)
        g_xsat : saturation point
    """

    fitfunc_str = 'A * x / (x + x_sat)'

    A = fit.Parameter(g_A, 'A')
    xsat = fit.Parameter(g_xsat, 'xsat')
    p0 = [A, xsat]

    def fitfunc(x):
        return A() * x / (x + xsat())

    return p0, fitfunc, fitfunc_str


def fit_saturation_with_offset_linslope(g_a, g_b, g_A, g_xsat, *arg):
    """
    fitfunction for a saturation (e.g., the NV PL)
        y(x) = a + b*x + A * x / (x + xsat)

    I.g.:
        g_a : offset
        g_b : linear slope
        g_A : maximum signal (at x=infinity)
        g_xsat : saturation point
    """

    fitfunc_str = 'a + b*x + A * x / (x + x_sat)'

    a = fit.Parameter(g_a, 'a')
    b = fit.Parameter(g_b, 'b')
    A = fit.Parameter(g_A, 'A')
    xsat = fit.Parameter(g_xsat, 'xsat')
    p0 = [a, b, A, xsat]

    def fitfunc(x):
        return a() + b()*x + A() * x / (x + xsat())

    return p0, fitfunc, fitfunc_str

def fit_poly(indices, *arg):
    fitfunc_str = 'sum_n ( a[n] * x**n )'

    idx = 0
    p0 = []
    for i,a in enumerate(indices):
        p0.append(fit.Parameter(a, 'a%d'%i))

    def fitfunc(x):
        val = 0
        for i in range(len(indices)):
            val += p0[i]() * x**i
        return val

    return p0, fitfunc, fitfunc_str

def fit_parabole(g_o, g_A, g_c, *arg):
    fitfunc_str = 'o + A * (x-c)**2'

    o = fit.Parameter(g_o, 'o')
    A = fit.Parameter(g_A, 'A')
    c = fit.Parameter(g_c, 'c')
    p0 = [o, A, c]

    def fitfunc(x):
        return o() + A() * (x-c())**2

    return p0, fitfunc, fitfunc_str


def fit_AOM_powerdependence(g_a, g_xc, g_k, *arg):
    fitfunc_str = 'a * exp(-exp(-k*(x-xc)))'

    a = fit.Parameter(g_a, 'a')
    xc = fit.Parameter(g_xc, 'xc')
    k = fit.Parameter(g_k, 'k')

    p0 = [a, xc, k]

    def fitfunc(x):
        return a() * np.exp(-np.exp(-k()*(x-xc())))

    return p0, fitfunc, fitfunc_str

def fit_gauss(g_a, g_A, g_x0, g_sigma):
### i think there should be a factor 2 infront of the sigma
    fitfunc_str = 'a + A * exp(-(x-x0)**2/(2*sigma**2))'

    a = fit.Parameter(g_a, 'a')
    x0 = fit.Parameter(g_x0, 'x0')
    A = fit.Parameter(g_A, 'A')
    sigma = fit.Parameter(g_sigma, 'sigma')

    p0 = [a, x0, A, sigma]

    def fitfunc(x):
        return a() + A() * np.exp(-(x-x0())**2/(2*sigma()**2))
    return p0, fitfunc, fitfunc_str

def fit_general_exponential(g_a, g_A, g_x0, g_T, g_n):
    fitfunc_str = 'a + A * exp(-((x-x0)/T )**n)'

    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    x0 = fit.Parameter(g_x0, 'x0')
    T = fit.Parameter(g_T, 'T')
    n = fit.Parameter(g_n, 'n')

    p0 = [a, A, x0, T, n]

    def fitfunc(x):
        return a() + A() * np.exp(-(x-x0())**n()/(T()**n()))
    return p0, fitfunc, fitfunc_str

def fit_2gauss(g_a1, g_A1, g_x01, g_sigma1, g_A2, g_x02, g_sigma2):
### i think there should be a factor 2 infront of the sigma
    fitfunc_str = 'a1 + A1 * exp(-(x-x01)**2/(2*sigma1**2)) +\
            A2 * exp(-(x-x02)**2/(2*sigma2**2))'

    a1 = fit.Parameter(g_a1, 'a1')
    x01 = fit.Parameter(g_x01, 'x01')
    A1 = fit.Parameter(g_A1, 'A1')
    sigma1 = fit.Parameter(g_sigma1, 'sigma1')

    x02 = fit.Parameter(g_x02, 'x02')
    A2 = fit.Parameter(g_A2, 'A2')
    sigma2 = fit.Parameter(g_sigma2, 'sigma2')


    p0 = [a1, x01, A1, sigma1, x02, A2, sigma2]

    def fitfunc(x):
        return a1()+A1()*np.exp(-(x-x01())**2/(2*sigma1()**2))+\
                A2()*np.exp(-(x-x02())**2/(2*sigma2()**2))
    return p0, fitfunc, fitfunc_str

def fit_lorentz(g_a, g_A, g_x0, g_gamma):
    fitfunc_str = 'a + 2*A/np.pi*gamma/(4*(x-x0)**2+gamma**2)'

    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    x0 = fit.Parameter(g_x0, 'x0')
    gamma = fit.Parameter(g_gamma, 'gamma')

    p0 = [a, A, x0, gamma]

    def fitfunc(x):
        return a() + 2*A()/np.pi*gamma()/(4*(x-x0())**2+gamma()**2)

    return p0, fitfunc, fitfunc_str

def fit_2lorentz(g_a1, g_A1, g_x01, g_gamma1, g_A2, g_x02, g_gamma2):
    fitfunc_str = 'a1 + 2*A1/np.pi*gamma1/(4*(x-x01)**2+gamma1**2) \
            + 2*A2/np.pi*gamma2/(4*(x-x02)**2+gamma2**2)'

    a1 = fit.Parameter(g_a1, 'a1')
    A1 = fit.Parameter(g_A1, 'A1')
    x01 = fit.Parameter(g_x01, 'x01')
    gamma1 = fit.Parameter(g_gamma1, 'gamma1')

    A2 = fit.Parameter(g_A2, 'A2')
    x02 = fit.Parameter(g_x02, 'x02')
    gamma2 = fit.Parameter(g_gamma2, 'gamma2')

    p0 = [a1, A1, x01, gamma1, A2, x02, gamma2]

    def fitfunc(x):
        return a1()+2*A1()/np.pi*gamma1()/(4*(x-x01())**2+gamma1()**2)+\
                2*A2()/np.pi*gamma2()/(4*(x-x02())**2+gamma2()**2)

    return p0, fitfunc, fitfunc_str


def fit_line(g_a, g_b, *arg):
    """
    fitfunction for a line
        y(x) = a + b*x

    I.g.:
        g_a : offset
        g_b : linear slope
    """

    fitfunc_str = 'a + b*x'

    a = fit.Parameter(g_a, 'a')
    b = fit.Parameter(g_b, 'b')
    #xsat = fit.Parameter(g_xsat, 'xsat')
    p0 = [a, b]

    def fitfunc(x):
        return a() + b()*x

    return p0, fitfunc, fitfunc_str


def fit_general_exponential_dec_cos(g_a, g_A, g_x0, g_T, g_n,g_f,g_phi):
    # NOTE: Order of arguments has changed to remain consistent with fitting params order 
    # NOTE: removed g_x0=0 as a default argument. This should be handed explicitly to the function to prevent confusion
    # Fits with a general exponential modulated by a cosine
    fitfunc_str = 'a + A * exp(-((x-x0)/T )**n*cos(2pi *(f*x+phi/360) )'

    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    x0 = fit.Parameter(g_x0, 'x0')
    T = fit.Parameter(g_T, 'T')
    n = fit.Parameter(g_n, 'n')
    f = fit.Parameter(g_f, 'f')
    phi = fit.Parameter(g_phi, 'phi')


    p0 = [a, A, x0, T, n,f,phi]
    def fitfunc(x):
        return a() + A() * np.exp(-((x-x0())/T())**n())*np.cos(2*np.pi*( f()*x + phi()/360.))
    return p0, fitfunc, fitfunc_str

