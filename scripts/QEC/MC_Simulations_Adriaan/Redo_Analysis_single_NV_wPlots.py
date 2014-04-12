import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt

import os   #needed for file and folder management
import datetime #Needed to put timestamps on filenames
import sys
sys.path.insert(0, 'Functions_For_Sims') #Required because functions in different subfolder
import time #Used for timing purposes

#Import custom functions

from NV_Reject import NV_Reject
from Find_Carb import Find_Carb
from addressable_C import addressable_C
from Find_Res_wPlot import Find_Res
from Char_Gate import Char_Gate
from FingerPrint import Fingerprint



def Redo_analysis (NV_List, idNV =0, idC = 0, idR = 0, B_Field = 600):
    """
    Plot the ax_prod that finds the resonances and the fingerprint for a specifc Carbon of one of the analyzed NV-Centres at a specified B-field
    Inputs
    """

    #Phyiscal Model paramters
    Carbon_Conc = 0.011 #1,1% C13
    #Simulation Parameters
    N_NV = 200#Number of NV Centres to generate
    Gridsize = 25
    A_Max = 180e3 #Max Hyperfine interaction strength in Hz
    A_Min = 10e3 # Min Hyperfine interaction strength !Allowed to be lower at lower Conc
    tau_max= 6e-6 # max resonance time to look for in s !Allowed to be longer if lower Conc
    tau_step = 1e-9 #timestep in s , was 2ns nu 1ns beste wat we experimenteel kunnen doen
    max_gate_time = 2e-3 #max gate time in s was 500e-6 kan hoger van tim (dubbel)!Allowed to be higher lower Conc
    F_Min = 0.85 # Minimum Gate Fidelity
    Gate_Min = 5 # Minimum number of gates

    NV = NV_List[idNV]
    times = np.arange(tau_step,tau_max,tau_step)
    axprod = np.zeros(np.size(times))
    print 'looking at NV centre #' +str(idNV)
    if NV_Reject(NV,A_Max):
        print 'center rejected because strongly coupled'
        return
    else:
        print 'All carbons have coupling < ' +str(A_Max*1e-3) +'kHz'
        indices = Find_Carb(NV,A_Min)
        print 'found ' +str(np.size(indices)) +' carbons that have coupling stronger than ' +str(A_Min *1e-3) +'kHz'
        print [indices]

        ind = indices [idC] #select the carbon atom specified by the function
        print 'looking at Carbon ' +str(ind)
        HyperfineStr = NV[ind,:]
        axprod,res = Find_Res(HyperfineStr,B_Field,tau_max,tau_step,max_gate_time) #Find the first resonances
        res = res[np.where(res[:,1]!=0)] #Reject False positives
        print 'found ' +str(np.shape(res)[0]) +' resonances'


        F ,F_of_Res= Fingerprint(NV,res[idR] ,B_Field,times)
        print 'Looking at resonance #' +str(idR) +'that has a Fidelity of ' +str(F_of_Res) + 'for ' +str(int(res[idR,1])) +'pulses'
    plt.figure()
    plt.ylabel(r'Product of $\hat{n_0} \hat{n_1}$')
    plt.xlabel(r'$\tau$')
    plt.vlines(res[:,0], -1, 1, colors='r',linewidth =2 )
    plt.plot(times,axprod)
    plt.ylim([-1.2,1.1])

    plt.figure()
    plt.ylabel(r'Fidelity of $\pi$ gate')
    plt.xlabel(r'$\tau$')
    plt.plot(times,F)
    plt.vlines(res[:,0], 0, 1, colors='r')
    plt.vlines(res[idR,0],0,1, colors = 'g', linewidth = 2)
    plt.ylim([0,1.1])


    plt.show()
    return


