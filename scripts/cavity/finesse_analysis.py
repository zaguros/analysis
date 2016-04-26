### Calculating cavity characteristics with both FSR and linewidth


# Inserting FSR from spectrometer analysis script
# Inserting linewidth from laser scan analysis script

#Calculating the Finesse of the cavity with an assumption for the linewidth

v1 = 620.e-9 # wavelength position one
v2 = 640.e-9 #wavelength position two
dv = 3*1.e9 # linewidth in Hz

#FSR in wavelength and frequency
FSR = v2 - v1
FSR_freq = (3.e8/v1) - (3.e8/v2)

# Finesse calculation
F = FSR_freq/dv
print 'The Finesse is', round(F,0)

#Cavity length calculation
L = 3.e8/(2*2.4*FSR_freq) # Calculating the length of the cavity in micrometer
print 'The Cavity Length is', round(L*1.e6,2), 'um.'

#Quality factor calculation
Q = (3.e8/v1)/dv
print 'The Q is', round(Q,0)

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




