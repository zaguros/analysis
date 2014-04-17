import numpy as np
import pylab as plt
import h5py

def Rabi_evolution(transition_driven='msm1'):
    if transition_driven =='msm1':
        h5filepath='/Users/Adriaan/Documents/teamdiamond/data_for_analysis/20140313/193451_ElectronRabi_Hans_sil1_Rabi-1/193451_ElectronRabi_Hans_sil1_Rabi-1.hdf5'
        # h5filepath='/Users/Adriaan/Documents/teamdiamond/data_for_analysis/20140312/172721_ElectronRabi_Hans_sil1_Rabi-1/172721_ElectronRabi_Hans_sil1_Rabi-1.hdf5'
        # h5filepath = '/Users/Adriaan/Documents/teamdiamond/data_for_analysis/20140319/104721_ElectronRabi_Hans_sil1_Rabi-1/104721_ElectronRabi_Hans_sil1_Rabi-1.hdf5'
    else:
        # h5filepath='/Users/Adriaan/Documents/teamdiamond/data_for_analysis/20140313/193621_ElectronRabi_Hans_sil1_Rabi+1/193621_ElectronRabi_Hans_sil1_Rabi+1.hdf5'
        # h5filepath='/Users/Adriaan/Documents/teamdiamond/data_for_analysis/20140312/172801_ElectronRabi_Hans_sil1_Rabi+1/172801_ElectronRabi_Hans_sil1_Rabi+1.hdf5'
        # h5filepath = '/Users/Adriaan/Documents/teamdiamond/data_for_analysis/20140319/104750_ElectronRabi_Hans_sil1_Rabi+1/104750_ElectronRabi_Hans_sil1_Rabi+1.hdf5'
        h5filepath ='/Users/Adriaan/Documents/teamdiamond/data_for_analysis/202845_PulsarMBIElectronRabi_hans1_finding_fast_rabi_frequency/202845_PulsarMBIElectronRabi_hans1_finding_fast_rabi_frequency.hdf5'
    #######################
    #### Model Parameters  ####
    #######################

    #Initial nitrogen-state Population
    a1 = 3/3.0  #part in the central dip
    a2 = 0/3.0  #part in the two side dips


    A1 = .9
    A2 = 0.65
    A3 = 0.0

    #initial guess electron-state population
    eP1 =  1 # 0.892/A1#Part of pupulation initialised in ms0
    eP2 =0#(1-eP1)*0.99 #Part of pupulation initialised in ms-1
    eP3 = 1-eP1-eP2 #Part of pupulation initialised in ms+1



    print 'C ms0 = %.2f' % eP1
    print 'C ms-1= %.2f' % eP2
    print 'C ms+1 = %.2f' % eP3

    print 'A1 =%.2f' %A1
    print 'A2 =%.2f' %A2
    print 'A3 =%.2f' %A3

    #Check if populations add up to 1
    if eP1+eP2+eP3 !=1:
        print "caution! sum of populations != 1"
    if A1 > 1:
        print 'caution probability of collecting photon from ms=0 >1'

    t_list = np.linspace(0.0,300.0,1000) #Time in ns

    #Rabi Frequency dependent on transition being driven
    if transition_driven == 'msm1':
        Omega_R = 2*np.pi/(105.0)
    elif transition_driven == 'msp1':
        Omega_R = 2*np.pi/(2*100.0)

    Delta = 2*np.pi * 2.16e6*1e-9

    Omega_prnt = Omega_R/(2*np.pi)*1e3
    Delta_prnt = Delta/(2*np.pi)*1e3
    prefactor = Omega_R**2/(Delta**2 +Omega_R**2)
    print'Omega_R = %.2f' %Omega_prnt
    print 'prefactor = %.2f' %prefactor


    P0 = 1-Omega_R**2/(Omega_R**2) *np.sin(t_list/2*np.sqrt(Omega_R**2))**2
    PD = 1-Omega_R**2/(Delta**2 +Omega_R**2) *np.sin(t_list/2*np.sqrt(Delta**2+Omega_R**2))**2

    PT = a1*P0 +a2*PD

    if transition_driven =='msm1':
        P_ms0 = eP1*PT+eP2*(1-PT)
        P_msm1 = eP2*PT +eP1*(1-PT)
        P_msp1 = eP3+t_list*0
    elif transition_driven =='msp1':
        P_ms0 = eP1*PT+eP3*(1-PT)
        P_msm1 = eP2+t_list*0
        P_msp1 = eP3*PT +eP1*(1-PT)

    #Sim_Osc = P_ms0 +eP2 +eP3 #Adds initial populations of eP2 and eP3 to signal so that it starts at 1
    Sim_Osc = A1*P_ms0+A2*P_msm1+A3*P_msp1
    #Normalised with factor

    ###########################
    ##### Importing the data #######
    ###########################
    f = h5py.File(h5filepath,'r')
    name = f.keys()[0]
    g = f[name]
    adwingrpname = g.keys()[1]
    adwingrp = g[adwingrpname]

    reps = adwingrp['completed_reps'].value
    sweep_pts = adwingrp.attrs['sweep_pts'] #in ns
    RO_data = adwingrp['RO_data'].value
    normalized_RO_data = RO_data/(float(reps/len(sweep_pts)))


    ###########################
    #### plotting of the data #######
    ###########################
    fig, ax1  = plt.subplots(nrows=1)
    # if transition_driven == 'msm1':
    #     ax1.set_title(r'Simulated Rabi driving $\mathrm{ms}_0 \leftrightarrow \mathrm{ms}_{-1}$ with: $\omega_R$ =%.1f  MHz, $\Delta$ = %.1f MHz' %(Omega_prnt, Delta_prnt))
    # else:
    #     ax1.set_title(r'Simulated Rabi driving $\mathrm{ms}_0 \leftrightarrow \mathrm{ms}_{+1}$ with: $\omega_R$ =%.1f  MHz, $\Delta$ = %.1f MHz' %(Omega_prnt, Delta_prnt))

    # ax0.plot(t_list,P_ms0, label ='ms0')
    # ax0.plot(t_list,P_msm1,label='ms-1')
    # ax0.plot(t_list,P_msp1,label='ms+1')
    # ax0.set_xlabel('time ns')
    # ax0.set_ylabel('population in msX state')
    # ax0.set_ylim(-.01,1)
    # ax0.grid(True)

    # ax0.legend(['ms0','ms-1','ms+1'])

    ax1.plot(t_list,Sim_Osc,label = 'simulated oscillation')
    ax1.plot(sweep_pts,normalized_RO_data, 'ro',label ='measured data')

    ax1.set_xlabel('time ns')
    ax1.set_ylabel(r'$\geq 1$ photon detected (ms0)')
    ax1.legend(['simulated oscillation','measured data (no RO corr)'],loc=0)
    ax1.set_ylim(0,1)
    ax1.grid(True)

    print normalized_RO_data[0]

Rabi_evolution(transition_driven='msm1')
Rabi_evolution(transition_driven='msp1')


plt.show()
