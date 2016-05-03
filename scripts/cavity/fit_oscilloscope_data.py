### fit oscilloscope data


import numpy as np 
from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot

def get_linewidth(data,EOM_freq, g_a1 = 0.5, g_A1 = 0.04 , g_x01 = 0.0, g_gamma1 = 0.1, g_dx = 0.35, g_A2 = 0.02, g_gamma2 = 0.1, plot_data = True):
    """
    load data, fit 3 lorentzians, get linewidth and the error
    input:
    data  = numpy ndarray with the data
    plot_data = False/True, default True
    g_a1 = offset guess parameter for fit
    g_A1 = guess area of peak 1 for fit
    g_x01 = guess x value of the fit 
    g_gamma1 = guess linewidth of middle peak
    g_dx = guess separation between peaks 
    g_A2 = guess area of peak 2 and 3 for fit
    g_gamma 2 = guess linewidth of peak 2 and 3 for fit
    output:
    lw = linewidth of middle peak
    u_lw = error of the linewidth
    """
    x = data[:,0]
    y = data[:,1]
    fixed=[]
    print x,y
    p0, fitfunc, fitfunc_str = common.fit_3lorentz_symmetric(g_a1, g_A1, g_x01, g_gamma1, g_dx, g_A2, g_gamma2)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)

    x01 = fit_result['params_dict']['x01']
    dx = fit_result['params_dict']['dx']
    gamma1 = fit_result['params_dict']['gamma1']
    gamma2 = fit_result['params_dict']['gamma2']
    u_gamma1 = fit_result['error_dict']['gamma1']
    scaling = EOM_freq/dx #scale the known EOM freq with the separation here.
    lw = gamma1*scaling #scale the linewidth to get linewidht in frequency
    u_lw = u_gamma1*scaling

    fig,ax = plt.subplots(figsize=(8,4))
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],len(x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)
    ax.set_xlabel("Frequency (GHz)", fontsize = 14)
    ax.set_ylabel("Intensity (a.u.)", fontsize = 14)
    ax.set_xlim(x[0],x[-1])
    plt.show()
    return lw, u_lw