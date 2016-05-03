### Calculating cavity characteristics with both FSR and linewidth

import numpy as np
from matplotlib import pyplot as plt
# Inserting FSR from spectrometer analysis script
# Inserting linewidth from laser scan analysis script

#Calculating the Finesse of the cavity with an assumption for the linewidth

v1 = 673.e-9 # wavelength position one
v2 = 690.e-9 #wavelength position two
dv = 2.09*1.e9 # linewidth in Hz
error_dv = 0.017*1.e9 #error linewidth in Hz
#FSR in wavelength and frequency
FSR = v2 - v1
FSR_freq = (3.e8/v1) - (3.e8/v2)

# Finesse calculation
F = FSR_freq/dv
F_error = FSR_freq/error_dv
print 'The Finesse is', round(F,0)

#Cavity length calculation
L = 3.e8/(2*FSR_freq) # Calculating the length of the cavity in micrometer
print 'The Cavity Length is', round(L*1.e6,2), 'um.'

#Quality factor calculation
Q = (3.e8/v1)/dv
print 'The Q is', round(Q,0)


# X = np.linspace(-2,10,13)
# print len(X)
# Y = [6.76,6.85,6.87,7.01,7.19,7.17,7.28,7.44,7.64,7.56,7.68,7.87,7.83]
# print len(Y)
# plt.plot(X,Y, 'ro')
# m, b = np.polyfit(X, Y, 1)
# plt.plot(X, m*X + b, '-')
# print m
# plt.xlabel("Voltage (V)", fontsize = 14)
# plt.ylabel("Cavity length (um)", fontsize = 14)

# plt.show()

# #Calculating the effective radius of curvature

# df_trans=fabs(peak_freq_1[-1]-peak_freq_1[-2]) # value from data: 2 nm (not conclusive yet!)
# T=(df_trans/FSR_freq)*pi
# ROC =l*(1/(1-(cos(T)**2)))
# print 'The ROC is', round(ROC*1.e6,2),'um.' 

# #Calculating the beam waist and mode volume 

# w_0= sqrt((peak_WL[-2]*1.e-9)/pi)*(l*(ROC-l))**(1/4)
# V_0 = (pi*(w_0**2)*l)/4

# print 'The beam waist is', round(w_0*1.e6,2)
# print 'The mode volume is', round(V_0*1.e12,2)

# Calculating the fine piezo calibration




