
''' Measurement set to calibrate the evolution frequencies for the different carbon spins 
plus the phases they pick up during carbon gates'''

import numpy as np
from analysis.scripts.QEC import carbon_ramsey_analysis as cr 
reload(cr)

A_list = np.zeros(15)
u_A_list = np.zeros(15)
# fit
n=0
timestamp = '20150108_145530'
while n<15:
	

	timestamp, folder = toolbox.latest_data('phase_C5',older_than = timestamp, return_timestamp = True)
	print timestamp

	ssro_timestamp, ssro_calib_folder = toolbox.latest_data('SSRO',older_than = timestamp, return_timestamp=True)

	A, u_A = cr.Carbon_Ramsey(timestamp=None, 
	                       offset = 0.5, amplitude = 0.5, x0=0, decay_constant = 1e5, exponent = 2, 
	                       frequency = 1/360., phase =0, 
	                       plot_fit = True, show_guess = False,fixed = [2,3,4,5],
			            return_phase = False,
			            return_freq = False,
			            return_amp = True,
			            return_results = True,
					title = 'phase_C5')
	A_list[n] = np.abs(A)
	u_A_list[n] = u_A

	n=n+1

print A_list
print u_A_list
print len(A_list)
print 2*np.sum(A_list)/len(A_list)
print 2*np.sum(u_A_list**2)/len(u_A_list)
