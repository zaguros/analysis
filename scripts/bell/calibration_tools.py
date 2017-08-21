import numpy as np
import matplotlib.pyplot as plt

from analysis.lib.fitting import fit
from analysis.lib.tools import plot
from analysis.lib.tools import toolbox

from analysis.lib.m2.ssro import sequence
reload(sequence)
from analysis.lib.m2.ssro import ssro
reload(ssro)

def fit_rabi(folder=None, ax = None, f_guess=0.9, A_guess=1,fit_phi = False,fit_k = False):
    """
    fit (1 - A) + A * exp(-kx) * cos(2pi f (x - x0)) from Sequence in folder
    """
    a, ax, x, y = plot_result(folder,ax)

    f = fit.Parameter(f_guess, 'f')
    A = fit.Parameter(A_guess, 'A')
    x0 = fit.Parameter(0, 'x0')
    k = fit.Parameter(0, 'k')
    p0 = [f, A]
    if fit_phi:
        p0.append(x0)
    if fit_k:
        p0.append(k)
    fitfunc_str = '(1 - A) + A * exp(-kx) * cos(2pi f (x - x0))'

    def fitfunc(x) : 
        return (1.-A()) + A() * np.exp(-k()*x) * \
            np.cos(2*np.pi*(f()*(x - x0())))

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

def fit_linear(folder=None, ax = None, a_guess=1., b_guess=0.):
    """
    fit a x + b from Sequence in folder
    """
    a, ax, x, y = plot_result(folder,ax)

    a = fit.Parameter(a_guess, 'a')
    b = fit.Parameter(b_guess, 'b')
    fitfunc_str = 'a x + b'

    def fitfunc(x):
        return a() * x + b()

    fit_result = fit.fit1d(x,y, None, p0=[a, b], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

def fit_gaussian(folder=None, ax = None, x0_guess=0., a_guess=0.5, c_guess=15,**kw):
    """
    fit o + a * exp( - ( (x-x0)/c )^2) from Sequence in folder
    """
    a, ax, x, y = plot_result(folder,ax)

    x0 = fit.Parameter(x0_guess, 'x0')
    a = fit.Parameter(a_guess, 'a')
    o = fit.Parameter(0.5, 'o')
    c = fit.Parameter(c_guess, 'c')
    fixed = kw.pop('fixed',[])
    fitfunc_str = 'o + a * exp( - ( (x-x0)/c )^2) '

    def fitfunc(x):
        return o() + a() * np.exp( -((x-x0())/ c())**2)

    fit_result = fit.fit1d(x,y, None, p0=[o,x0,a,c], fixed = fixed, fitfunc=fitfunc,
            fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
            plot_data=False,**kw)

    return fit_result

def fit_parabolic(folder=None, ax = None,x0_guess=0., of_guess=0.5, a_guess=1.,**kw):
    """
    fit (1-of) + a (x-x0)**2 from Sequence in folder
    """
    do_print      = kw.pop('do_print', False)
    fixed = kw.pop('fixed', [])
    a0, ax, x, y = plot_result(folder,ax,**kw)
    x0 = fit.Parameter(x0_guess, 'x0')
    of = fit.Parameter(of_guess, 'of')
    a = fit.Parameter(a_guess, 'a')

    x_min = kw.pop("x_min", None)
    x_max = kw.pop("x_max", None)

    fit_x = x
    fit_y = y

    if x_min is not None:
        mask = fit_x > x_min
        fit_x = fit_x[mask]
        fit_y = fit_y[mask]

    if x_max is not None:
        mask = fit_x < x_max
        fit_x = fit_x[mask]
        fit_y = fit_y[mask]

    fitfunc_str = '(1-of) + a (x-x0)**2'
    def fitfunc_parabolic(x):
        return (1.-of()) + a() * (x-x0())**2
    
    fit_result = fit.fit1d(fit_x, fit_y, None, p0=[of, a, x0], fitfunc=fitfunc_parabolic,
        fitfunc_str=fitfunc_str, do_print=do_print, ret=True)
    fit_result['u_y'] = a0.u_p0
    plot.plot_fit1d(fit_result, np.linspace(np.min(fit_x), np.max(fit_x), 201), ax=ax,
        plot_data=False, **kw)
        
        
    return fit_result

def plot_result(folder=None, ax = None, ret=True, save_name = 'ssro',**kw):
    if folder==None:
        folder= tb.latest_data('CORPSEPiCalibration')
    if ax == None:
        fig, ax = plt.subplots(1,1, figsize=(4.5,4))
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name=save_name)
    a.get_electron_ROC()
    print len(a.p0),len(a.sweep_pts)
    a.plot_result_vs_sweepparam(ax=ax, name=save_name)
    # print a.normalized_ssro
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    print 'min x = ', x[np.argmin(y)]
    print 'min y =', np.amin(y)
    ax.set_title(a.timestamp)
    if ret:
    	return a, ax, x, y