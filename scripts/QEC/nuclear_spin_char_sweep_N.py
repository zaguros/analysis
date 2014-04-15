#Simulates the response of a single carbon atom to a number of decoupling pulses
#Script by MAR
import numpy as np
import matplotlib.pyplot as plt
import h5py
import analysis.lib.qec.nuclear_spin_characterisation as SC
reload(SC)

################
#input arguments ##
################

def main(labpc = False):


    ##############
    ## exp params ##
    ##############
    N = np.array(range(0,124,4)) #Integer
    # tau= 6.522e-6 #seconds
    tau = 6.622e-6
    B_Field = 304.12 #Gauss
    pulseF = .86 #correction factor on final signal for plotting (pulse fidelity )


    ####################
    ## Hyperfine strenghts ##
    # ####################
    # HF_par = [28e3,-62.5e3, 27e3]
    # HF_orth =[90e3,130e3,30e3]

    HF_par = [28e3,-62.5e3, 27e3]
    HF_orth =[26e3,130e3,30e3]

    Mn = SC.dyn_dec_signal(HF_par,HF_orth,B_Field,N,tau)
    sweep_N_signal = ((Mn+1)/2)
    #########################
    ##### Importing data #######
    #########################
    if labpc == True:
        pass
        #Script that uses toolbox functions
    else:
        filepaths = []
        # datestamp = '20140410'
        # timestamps = ['165509']
        datestamp = '20140408'
        # timestamps = ['174302']
        timestamps = ['164854']

        datadir = '/Users/Adriaan/Documents/teamdiamond/data/'
        filepath = datadir+'/'+datestamp+'/TIMESTAMP_DecouplingSequence_Hans_sil1_C2_sweep_N/TIMESTAMP_DecouplingSequence_Hans_sil1_C2_sweep_N.hdf5'





        for i, timestamp in enumerate(timestamps):
            filepaths.append(filepath.replace("TIMESTAMP",timestamp))


        data_N= []
        signal = []
        for filep in filepaths:
            d_N, sig =  load_dd_fingerprintdata(filep)
            data_N.extend(d_N)
            signal.extend(sig)


    ##########
    #Plotting
    ##########
    plt.figure()
    plt.xlabel('N')
    plt.ylabel('Signal')

    plt.plot(N,sweep_N_signal[0,:]*pulseF,'r.-',label='test_spin1')
    plt.plot(N,sweep_N_signal[1,:]*pulseF,'g.-',label='test_spin2')
    plt.plot(data_N,signal,'b.-',linewidth = 0.5,label='measured_data')
    plt.legend(loc=4)
    plt.show()


def load_dd_fingerprintdata(h5filepath):
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
# estimate_HFs()
# load_dd_fingerprintdata()
