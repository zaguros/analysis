import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt

from scipy.stats.stats import pearsonr
reload(common)


N = [4,8,16,32,64]#,128,256]
N = [1,2,4,8,16,32,64,128,256,512,1024]
T_list=[1.32,2.2,3.58,5.9,11.6]#,16.26,21.23]
T_list=[]
#################### plot and fit scaling with N
fig = plt.figure(3,figsize=(8,6))
ax2 = fig.add_subplot(111)
ax2.plot(N,T_list,'bo')




    
B_guess = T_list[0]
n_guess = 0.78         

#A = fit.Parameter(A_guess, 'A')
B = fit.Parameter(B_guess, 'B')
n=  fit.Parameter(n_guess, 'n')

def fitfunc(x):
            return  B()*x**n()

print 'running fit'
fit_result = fit.fit1d(N, T_list, None, p0 = [ B, n],
                fitfunc = fitfunc, do_print=True, ret=True, fixed=[])

#A0 = fit_result['params_dict']['A']
B0 = fit_result['params_dict']['B']
n0 = fit_result['params_dict']['n']

    # plotting the fitted function
plot.plot_fit1d(fit_result, np.linspace(min(N), max(N), 1000), ax=ax2, plot_data=True)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('N')
ax2.set_ylabel('coherence time (ms)')


################ Hans data

N = [4,8,16,32,64]#,128,256]
#N = [1,2,4,8,16]#,32,64]
T_list=[2.4,4.5,5.3,9.6,19]#,16.26,21.23]
#################### plot and fit scaling with N
fig = plt.figure(4,figsize=(8,6))
ax3 = fig.add_subplot(111)
ax3.plot(N,T_list,'bo')




    
B_guess = T_list[0]
n_guess = 0.78         

#A = fit.Parameter(A_guess, 'A')
B = fit.Parameter(B_guess, 'B')
n=  fit.Parameter(n_guess, 'n')

def fitfunc(x):
            return  B()*x**n()

print 'running fit'
fit_result = fit.fit1d(N, T_list, None, p0 = [ B, n],
                fitfunc = fitfunc, do_print=True, ret=True, fixed=[])

#A0 = fit_result['params_dict']['A']
B0 = fit_result['params_dict']['B']
n0 = fit_result['params_dict']['n']

    # plotting the fitted function
plot.plot_fit1d(fit_result, np.linspace(min(N), max(N), 1000), ax=ax3, plot_data=True)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel('N')
ax3.set_ylabel('coherence time (ms)')
