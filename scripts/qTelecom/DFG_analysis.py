
import numpy as np
import pylab as plt
from analysis.lib import fitting
from analysis.lib.fitting import fit,esr, common

def fit_sweep_pump_power(x, y, L=0.042):
	guess_b = 1
	guess_a = max(y)

	a = fit.Parameter(guess_a, 'a')
	b = fit.Parameter(guess_b, 'b')
	
	p0 = [a, b]
	fitfunc_str = ''

	def fitfunc(x):
		return a()*np.sin(L*(b()*x)**0.5)**2


	fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=[],
        	do_print=True, ret=True)
	a_fit = fit_result['params_dict']['a']
	b_fit = fit_result['params_dict']['b']
	print 'a= ',a_fit
	print 'b=',b_fit
	return a_fit, b_fit



pump = np.array ([50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 820, 900, 1000])
signal = np.array ([1.47, 2.79, 4.06, 5.183, 6.3, 7.26, 8.02, 8.87, 9.6, 10.2, 10.72, 11.4, 11.8, 12.17, 12.57, 12.85, 13.1, 13.5])
ind = np.argsort(pump)
pump = pump[ind]*0.001
eta = signal[ind]*1587./(400*637) # red input power : 400 uW


plt.plot (pump, eta*100, 'r')
plt.plot (pump, eta*100, 'ob')
plt.xlabel ('pump [W]', fontsize = 15)
plt.ylabel ('Conversion efficiency [%]', fontsize=15)
plt.show()

a, b = fit_sweep_pump_power (x = pump, y = eta)
p = np.linspace (pump[0], pump[-1], 10000)
y = a*np.sin(0.042*(b*p)**0.5)**2
plt.plot (pump, eta*100, 'ob')
plt.plot (p, y*100, 'r')
plt.xlabel ('pump [W]', fontsize = 15)
plt.ylabel ('Conversion efficiency [%]', fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.title("2015-08-06\nFirst DFG signal\nSignal without fiber", fontsize = 15)
plt.show()
