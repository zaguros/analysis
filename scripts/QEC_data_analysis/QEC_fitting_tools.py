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
    g_A -  Amplitude, given by the fidelity of the states that are sensitive  
    g_p -  Avegage probabililty to correct single qubit errors
    '''
 
    fitfunc_str = '''  '''

    O   = fit.Parameter(g_O , 'O')
    A   = fit.Parameter(g_A, 'A')
    p   = fit.Parameter(g_p, 'p')

    p0 = [O, A, p]

    def fitfunc(x):
    	return (O() + A()*(  1-3*x+3*x**2-2*x**3 + 3*(2*p()-1)*(x-3*x**2+2*x**3)))
    
	return p0, fitfunc, fitfunc_str

### Make fake data and plot it
data_O = 0.5 
data_A = 0.5
data_p = 1

p0_data, fitfunc_data, fitfunc_str = fit_QEC(data_O, data_A, data_p)

xdata = linspace(0,1,20)
ydata = fitfunc(xdata)

fig = a_up.default_fig(figsize=(10,5))
ax  = a_up.default_ax(fig)
ax.set_xlim(xdata[0]-0.05,xdata[-1]+0.05)

ax.plot(xdata,   ydata, markersize = 6, label = 'data') 


# ### Fit: inital guesses and fitfunc

# guess_O = 0.5 
# guess_A = 0.5
# guess_p = 1

# p0, fitfunc, fitfunc_str = fit_QEC(guess_O, guess_A, guess_p)

# fit_result = fit.fit1d(x_all, y_all, fit_initialization_fidelity_all,
#         guess_F1, guess_A1, guess_F2, guess_A2,
#         guess_f, guess_a, guess_A, guess_tau,
#          fixed=[3],
#         do_print=True, ret=True)




# ### show guess ###
# if 1:
#     x_temp = linspace(xdata[0],xdata[-1],200)
#     y_temp = fitfunc(x_temp)

#     ax.plot(x_temp, guess_curves_no_init, 'b', lw=2)

# ### Fitting

