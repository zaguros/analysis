### Tools for fitting QEC data, THT

import numpy as np
import os
from analysis.lib.tools import toolbox; reload(toolbox)
from analysis.lib.tools import plot; reload(plot)
from analysis.lib.fitting import fit, common, ramsey;reload(common); reload(fit)
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt


### Fitfunction

def fit_QEC(g_O, g_A, g_p):
    '''Fit function for QEC process fidelity data
    g_O -  Offset, given by the fidelity of the state that is insensitive to errors
    g_A -  Amplitude, g_iven by the fidelity of the states that are sensitive
    g_p -  Avegage probabililty to correct single qubit errors
    '''

    fitfunc_str = '''test'''

    O   = fit.Parameter(g_O , 'O')
    A   = fit.Parameter(g_A, 'A')
    p   = fit.Parameter(g_p, 'p')

    p0 = [O, A, p]

    def fitfunc(x):
        '''test'''
        return (O() + A()*(  1-3*x+3*x**2-2*x**3 + 3*(2*p()-1)*(x-3*x**2+2*x**3)))

    return p0, fitfunc, fitfunc_str

##################################
### Make fake data and plot it ###
##################################
data_O = 0.4
data_A = 0.4
data_p = 0.85

p0, fitfunc, fitfunc_str = fit_QEC(data_O, data_A, data_p)

xdata = np.linspace(0,1,20)
ydata = fitfunc(xdata)

plt.figure()
plt.plot(xdata,ydata,'bo')

###########
### Fit ###
###########

guess_O = 0.5
guess_A = 0.5
guess_p = 1

p0, fitfunc, fitfunc_str = fit_QEC(guess_O, guess_A, guess_p)

fit_result = fit.fit1d(xdata, ydata, fit_QEC,
        guess_O, guess_A, guess_p,
        fixed=[],
        do_print=True, ret=True)

### show guess ###
if 1:
    x_temp = np.linspace(xdata[0],xdata[-1],200)
    y_temp = fitfunc(x_temp)

    plt.plot(x_temp, y_temp, 'b', lw=2)

p02, fitfunc2, fitfunc_str2 = fit_QEC(fit_result['params_dict']['O'], fit_result['params_dict']['A'], fit_result['params_dict']['p'])

plt.plot(x_temp, fitfunc2(x_temp),'r')


plt.show()
