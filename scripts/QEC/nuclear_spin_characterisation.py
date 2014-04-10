#Simulates the response of a single carbon atom to a number of decoupling pulses
#Script by MAR
import numpy as np
import matplotlib.pyplot as plt
import h5py

################
#input arguments ##
################
def main(labpc = False):

    ##############
    ## exp params ##
    ##############
    N = 16 #Integer
    tau=np.linspace(0, 15e-6, 15000) #seconds
    B_Field = 300 #Gauss
    pulseF = .86 #correction factor on final signal for plotting (pulse fidelity )



    ####################
    ## Hyperfine strenghts ##
    ####################

    HF_par = [30e3,-23e3, 178e3]
    HF_orth =[100e3,38e3,80e3]

    fingerprint_signal = Fingerprint(HF_par,HF_orth,B_Field,N,tau)

    #########################
    ##### Importing data #######
    #########################
    if labpc == True:
        pass
        #Script that uses toolbox functions
    else:
        filepaths = []
        filepath = '/Users/Adriaan/Documents/Python_Programs/Data/fingerprint_data/TIMESTAMP_DecouplingSequence_Fingerprint_Hans_sil1_N1024/TIMESTAMP_DecouplingSequence_Fingerprint_Hans_sil1_N1024.hdf5'
        # timestamps = ['223450', '223801', '224119', '224446', '231158','225215', '225614', '230023', '231645', '232118', '232603', '233057', '233605']
        timestamps = ['223450', '223801', '224119', '224446', '231158','225215', '225614', '230023', '231645', '232118', '232603', '233057', '233605', '234654', '111324', '111725', '112126', '112529', '112930', '113614', '114015', '114416', '114818', '115220', '115622', '120024', '120426','120825']#, '130753']


        for i, timestamp in enumerate(timestamps):
            filepaths.append(filepath.replace("TIMESTAMP",timestamp))
        data_tau= []
        signal = []
        for filep in filepaths:
            d_tau, sig =  load_dd_fingerprintdata(filep)
            data_tau.extend(d_tau)
            signal.extend(sig)


    ##########
    #Plotting
    ##########
    plt.figure()
    plt.xlabel('tau(us)')
    plt.ylabel('Signal')
    # plt.plot(tau*1e6,fingerprint_signal[0,:]*pulseF,'r',label='test_spin1')
    # plt.plot(tau*1e6,fingerprint_signal[1,:]*pulseF,'g',label = 'test_spin2')
    plt.plot(tau*1e6,fingerprint_signal[2,:]*pulseF,'m',label = 'test_spin3')
    # plt.plot(tau*1e6,fingerprint_signal.prod(axis=0)*pulseF,'c--',label ='Combined signal of test spins',linewidth=.5)
    plt.plot(data_tau,signal,'b.--',linewidth = 0.5,label='measured_data')
    plt.legend(loc=4)
    plt.show()




def Fingerprint(HFs_par,HFs_orth,B_field,N,tau,input_in_radial_freq = False):
    '''
    Takes the HF interaction strengths (paralel and orthogonal), the magnetic field strenght and an array of times and returns the signal at those times for that specific spin.
    '''
    #physical constants
    gamma_c = 1.071e3 #g-factor for C13 in Hz/G
    #Model parameters
    omega_larmor = 2*np.pi*gamma_c*B_field #radial frequency
    tau_larmor = 2*np.pi/omega_larmor #time in seconds


    print 'tau larmor = %s' %tau_larmor
    M=np.zeros([np.size(HFs_par),np.size(tau)])
    for i,HF_par in enumerate(HFs_par):
        if input_in_radial_freq == False:
            HF_par = HF_par*2*np.pi #Convert to radial frequency
            HF_orth = HFs_orth[i]*2*np.pi #convert to radial frequency
        omega_tilde = np.sqrt((HF_par+omega_larmor)**2+HF_orth**2)
        alpha = omega_tilde*tau
        beta = omega_larmor*tau
        mx = HF_orth/omega_tilde
        mz = (HF_par+omega_larmor)/omega_tilde
        vec_term = mx**2 *((1-np.cos(alpha))*(1-np.cos(beta)))/(1+np.cos(alpha)*np.cos(beta)-mz*np.sin(alpha)*np.sin(beta))
        angle_term = np.sin(N*np.arccos(np.cos(alpha)*np.cos(beta)-mz*np.sin(alpha)*np.sin(beta))/2)**2

        M[i,:]= 1-(vec_term*angle_term)

    Signal = M.prod(axis=0)
    F = ((M+1)/2)
    return F

def load_dd_fingerprintdata(h5filepath='/Users/Adriaan/Documents/Python_Programs/Data/fingerprint_data/233605_DecouplingSequence_Fingerprint_Hans_sil1_N1024/233605_DecouplingSequence_Fingerprint_Hans_sil1_N1024.hdf5'):
    '''
    Loads fingerprint data when not on a lab pc
    '''
    f = h5py.File(h5filepath,'r')
    name = f.keys()[0]
    g = f[name]
    adwingrpname = g.keys()[0]
    adwingrp = g[adwingrpname]
    sweep_length = adwingrp.attrs['sweep_length']
    reps_per_ROsequence = adwingrp.attrs['reps_per_ROsequence']

    ssro_results = np.reshape((adwingrp['ssro_results'].value),(reps_per_ROsequence,sweep_length))
    sweep_pts = adwingrp.attrs['sweep_pts'] #in us
    RO_data = np.sum(ssro_results,axis=0)/float(reps_per_ROsequence)

    xdata = sweep_pts
    ydata = RO_data
    return xdata,ydata

main()
# load_dd_fingerprintdata()
