#Simulates the response of a single carbon atom to a number of decoupling pulses
#Script by MAR
import numpy as np
import matplotlib.pyplot as plt
import h5py
import analysis.lib.qec.nuclear_spin_characterisation as SC

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
    HF_par = [30e3,-23e3, -274e3]
    HF_orth =[100e3,38e3,15e3]

    M = SC.dyn_dec_signal(HF_par,HF_orth,B_Field,N,tau)
    fingerprint_signal = ((M+1)/2)
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

    plt.plot(tau*1e6,fingerprint_signal[0,:]*pulseF,'r',label='test_spin1')
    plt.plot(tau*1e6,fingerprint_signal[1,:]*pulseF,'g',label = 'test_spin2')
    # plt.plot(tau*1e6,fingerprint_signal[2,:]*pulseF,'m',label = 'test_spin3')
    # plt.plot(tau*1e6,fingerprint_signal.prod(axis=0)*pulseF,'c--',label ='Combined signal of test spins',linewidth=.5)
    plt.plot(data_tau,signal,'b.--',linewidth = 0.5,label='measured_data')
    plt.legend(loc=4)
    plt.show()


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
