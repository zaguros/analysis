import numpy as np

################################################
######### SAMPLE INPUT PARAMETERS ##############
################################################

N_larger = 820.     #number of readout events with more than 0 photons
N_equal = 180.      #number of readout events with 0 photons

F0 = 0.82885
F1 = 0.98575

F0_err = 0.01
F1_err = 0.01

def correct(N_larger, N_equal, F0, F1, F0_err, F1_err, verbose = True):
    """
    Returns the read-out corrected amplitudes for given input parameters.
    Input: 
        N_larger        number of readout events with more than 0 photons
        N_equal         number of readout events with 0 photons
        F0              Fidelity of ms = 0 readout obtained by a SSRO measurement
        F1              Fidelity of ms = 1 readout 
        F0_err          Uncertainty in F0 (1 sigma)
        F1_err          Uncertainty in F1 (1 sigma)

    Output is a dictionary with elements:
        r0              corrected ms = 0 amplitude
        sigma_r0        standard deviation
        r1              corrected ms = 1 amplitude
        sigma_r1        standard deviation
    """
    
    #########################################
    ####### CALCULATION OF THE READOUT ERRORS
    #########################################

    m0 = N_larger/float(N_larger + N_equal)
    m1 = N_equal/float(N_larger + N_equal)

    sigma_m0 = np.sqrt((N_equal/float(N_larger+N_equal)**2)**2 * N_larger \
            + (N_larger/float(N_larger+N_equal)**2)**2 * N_equal)
    sigma_m1 =  np.sqrt((N_larger/float(N_equal+N_larger)**2)**2 * N_equal \
            + (N_equal/float(N_equal+N_larger)**2)**2 * N_larger)

    if verbose:
        print 'Analyzing read-out errors for state psi = a |0> + b |1>\n'
        print 'Measured quantities with statistical errors (1 sigma):'
        print 'a = %5.3f(%5.3f) and b = %5.3f(%5.3f)\n'%(m0,sigma_m0, m1, sigma_m1)

    sigma_F0 = F0_err
    sigma_F1 = F1_err

    norm = F0+F1-1

    r0 = 1/norm*(m0*F1+m1*(F1-1))

    r1 = 1/norm*(m0*(F0-1)+F0*m1)

    sigma_r0 = np.sqrt(((m0*F0+m1*(F1-1))/norm**2)**2 * sigma_F0**2 \
            + ((m0*(F0-1)+m1*F0)/norm**2)**2 *sigma_F1**2 \
            + (F1/norm)**2 *sigma_m0**2 \
            + ((F1-1)/norm)**2 *sigma_m1**2)

    sigma_r1 = np.sqrt(((m1-F1*(m0+m1))/norm**2)**2 * sigma_F0**2 \
            + ((m0*(F0-1)+m1*F0)/norm**2)**2 * sigma_F1**2 \
            + ((F0-1)/norm)**2 * sigma_m0**2 \
            + (F0/norm)**2 * sigma_m1**2)

    if verbose:
        print 'After correction for readout (1 sigma errors):'
        print 'a = %5.3f(%5.3f) and b = %5.3f(%5.3f)\n'\
                %(r0,sigma_r0, r1, sigma_r1)

    return {'r0' : r0, 'sigma_r0' : sigma_r0, 'r1' : r1, 'sigma_r1' : sigma_r1}



########################################
#### ENTANGLEMENT EVENT CALCULATION ####
########################################


N_equal_min = 0
N_equal_max = 50

N_larger_min = 50
N_larger_max = 100

b = np.zeros([N_equal_max-N_equal_min+1,N_equal_max-N_equal_min+1])
c = np.zeros([N_larger_max-N_larger_min+1,N_larger_max-N_larger_min+1])
for N_larger in range(N_larger_min,N_larger_max+1):
    for N_equal in range(N_equal_min,N_equal_max+1):
        a = correct(N_larger, N_equal, F0, F1, F0_err, F1_err, verbose = False)
        b[N_larger-N_larger_min, N_equal-N_equal_min] = a['sigma_r0']
        c[N_larger-N_larger_min, N_equal-N_equal_min] = a['sigma_r1']

N_larger = range(N_larger_min,N_larger_max+1)
N_equal = range(N_equal_min,N_equal_max+1)

X,Y = meshgrid(N_equal, N_larger)

plt.figure()
plt.pcolor(X, Y, b)
plt.colorbar()
tick_xlocs = range(N_equal_min, N_equal_max, 5)
tick_ylocs = range(N_larger_min, N_larger_max, 5)

#string the shit
xticklabels = yticklabels = list()
for k in tick_xlocs:
    xticklabels.append(str(k))
for k in tick_ylocs:
    yticklabels.append(str(k))
    
#plt.xticks(tick_xlocs, xticklabels)
#plt.yticks(tick_ylocs, yticklabels)

plt.xlabel('$N_{=0}$')
plt.ylabel('$N_{>0}$')
plt.title('Final error of the read-out corrected state')
plt.xlim([N_equal_min,N_equal_max])
plt.ylim([N_larger_min,N_larger_max])



