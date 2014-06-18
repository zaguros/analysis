import numpy as np
import pylab as plt

t_list = np.linspace(0.0,300.0,1000) #Time in ns

Omega_R = 2*np.pi/(2*120.0)
Delta = 2*np.pi * 2.16e6*1e-9

print Omega_R
print Delta

Omega_prnt = Omega_R/(2*np.pi)*1e3
Delta_prnt = Delta/(2*np.pi)*1e3
prefactor = Omega_R**2/(Delta**2 +Omega_R**2)
print 'prefactor = %.2f' %prefactor


P0 = 1-Omega_R**2/(Omega_R**2) *np.sin(t_list/2*np.sqrt(Omega_R**2))**2
P1 = 1-Omega_R**2/(Delta**2 +Omega_R**2) *np.sin(t_list/2*np.sqrt(Delta**2+Omega_R**2))**2

Sum_osc = (P0+2*P1)/3.0

plt.plot(t_list,Sum_osc)
plt.xlabel('time ns')
plt.ylabel('F(|0>)')
plt.title(r'Simulated Electron Rabi: $\omega_R$ =%.1f  MHz, $\Delta$ = %.1f MHz' %(Omega_prnt, Delta_prnt))
plt.grid(True)

plt.show()
