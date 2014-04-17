import numpy as np
import pylab as plt
import h5py

def Rabi_evolution(transition_driven='msm1'):
    #######################
    #### Model Parameters  ####
    #######################

    #Initial nitrogen-state Population
    a1 = .95#1/3.0#2/3.0  #part in the central dip
    a2 = 1-a1 #2/3.0  #part in the two side dips


    A1 = 1#.9
    A2 = 0.65
    A3 = 0.0
    #initial guess electron-state population
    eP1 =  1# 0.892/A1#Part of pupulation initialised in ms0
    eP2 =0#(1-eP1)*0.99 #Part of pupulation initialised in ms-1
    eP3 = 1-eP1-eP2 #Part of pupulation initialised in ms+1

    print 'C ms0 = %.2f' % eP1
    print 'C ms-1= %.2f' % eP2
    print 'C ms+1 = %.2f' % eP3

    #Check if populations add up to 1
    if eP1+eP2+eP3 !=1:
        print "caution! sum of populations != 1"
    if A1 > 1:
        print 'caution probability of collecting photon from ms=0 >1'

    t_list = np.linspace(0.0,3000.0,1000) #Time in ns
    Omega_R = 2*np.pi/(2*100.0)
    Delta = 2*np.pi * 2.16e6*1e-9

    Omega_prnt = Omega_R/(2*np.pi)*1e3
    Delta_prnt = Delta/(2*np.pi)*1e3
    prefactor = Omega_R**2/(Delta**2 +Omega_R**2)
    print'Omega_R = %.2f' %Omega_prnt
    print 'prefactor = %.2f' %prefactor


    P0 = 1-Omega_R**2/(Omega_R**2) *np.sin(t_list/2*np.sqrt(Omega_R**2))**2
    PD = 1-Omega_R**2/(Delta**2 +Omega_R**2) *np.sin(t_list/2*np.sqrt(Delta**2+Omega_R**2))**2

    EnvelopeFunction = 1- np.sin(t_list/2*((np.sqrt(Delta**2+Omega_R**2))-Omega_R))**2/2.0
    PT = a1*P0 +a2*PD

    if transition_driven =='msm1':
        P_ms0 = eP1*PT+eP2*(1-PT)
        P_msm1 = eP2*PT +eP1*(1-PT)
        P_msp1 = eP3+t_list*0
    elif transition_driven =='msp1':
        P_ms0 = eP1*PT+eP3*(1-PT)
        P_msm1 = eP2+t_list*0
        P_msp1 = eP3*PT +eP1*(1-PT)

    Sim_Osc = A1*P_ms0+A2*P_msm1+A3*P_msp1
    #Normalised with factor

    ###########################
    #### plotting of the data #######
    ###########################
    fig, ax1  = plt.subplots(nrows=1)

    ax1.plot(t_list*1e-3,Sim_Osc,label = 'simulated oscillation')
    ax1.plot(t_list*1e-3,EnvelopeFunction, label = 'envelope')

    # ax1.plot(sweep_pts,normalized_RO_data, 'ro',label ='measured data')

    ax1.set_xlabel('time us')
    ax1.set_ylabel(r'$\geq 1$ photon detected (ms0)')
    # ax1.set_ylim(0,1)
    ax1.grid(True)

    # print normalized_RO_data[0]

# Rabi_evolution(transition_driven='msm1')
Rabi_evolution(transition_driven='msp1')


plt.show()
