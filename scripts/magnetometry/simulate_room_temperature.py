
import numpy as np
import pylab as plt
from matplotlib import rc

N = np.arange (50000)+1

a0 = 0.031*N
a1 = 0.021*N
C = 1./(1+(2*(a0+a1))/((a0-a1)**2))**0.5

matplotlib.rc('xtick', labelsize=15) 
matplotlib.rc('ytick', labelsize=15) 
plt.figure()
plt.loglog (N, 1-C, 'RoyalBlue', linewidth=3)
plt.xlabel ('number of msmnt read-outs', fontsize=15)
plt.ylabel ('read-out error', fontsize=15)
plt.show()

f = plt.figure()
plt.semilogx (N[200:], C[200:], 'RoyalBlue', linewidth=3)
#plt.hlines(y=0.75, xmin = 100, xmax=1e5, color = 'crimson')
#plt.hlines(y=0.88, xmin = 100, xmax=1e5, color = 'crimson')
plt.xlabel ('number of msmnt read-outs', fontsize=15)
plt.ylabel ('read-out fidelity', fontsize=15)
f.savefig ('d:/measuring/roomT_fidelity.svg')
plt.show()