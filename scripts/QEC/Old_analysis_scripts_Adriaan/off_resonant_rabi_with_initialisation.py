import numpy as np
import pylab as plt
def Rabi_evolution(transition_driven='msp1'):
    #Initial nitrogen-state Population
    a1 = 1/3.0  #part in the central dip
    a2 = 2/3.0  #part in the two side dips

    #initial electron-state population
    eP1 = .9#Part of pupulation initialised in ms0
    eP2 = .04 #Part of pupulation initialised in ms-1
    eP3 = 1-eP1-eP2 #Part of pupulation initialised in ms+1
    #Check if populations add up to 1
    if eP1+eP2+eP3 !=1:
        print "caution! sum of populations != 1"

    t_list = np.linspace(0.0,300.0,1000) #Time in ns

    #Rabi Frequency dependent on transition being driven
    if transition_driven == 'msm1':
        Omega_R = 2*np.pi/(100.0)
    elif transition_driven == 'msp1':
        Omega_R = 2*np.pi/(2*120.0)

    Delta = 2*np.pi * 2.16e6*1e-9

    print Omega_R
    print Delta

    Omega_prnt = Omega_R/(2*np.pi)*1e3
    Delta_prnt = Delta/(2*np.pi)*1e3
    prefactor = Omega_R**2/(Delta**2 +Omega_R**2)
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

    Normalised_Osc = P_ms0 +eP2 +eP3
    fig, (ax0, ax1)  = plt.subplots(nrows=2)

    if transition_driven == 'msm1':
        ax0.set_title(r'Simulated Rabi driving $\mathrm{ms}_0 \leftrightarrow \mathrm{ms}_{-1}$ with: $\omega_R$ =%.1f  MHz, $\Delta$ = %.1f MHz' %(Omega_prnt, Delta_prnt))
    else:
        ax0.set_title(r'Simulated Rabi driving $\mathrm{ms}_0 \leftrightarrow \mathrm{ms}_{+1}$ with: $\omega_R$ =%.1f  MHz, $\Delta$ = %.1f MHz' %(Omega_prnt, Delta_prnt))
    ax0.plot(t_list,P_ms0, label ='ms0')
    ax0.plot(t_list,P_msm1,label='ms-1')
    ax0.plot(t_list,P_msp1,label='ms+1')
    ax0.set_xlabel('time ns')
    ax0.set_ylabel('population in msX state')
    ax0.set_ylim(-.01,1)
    ax0.grid(True)
    # handles, labels = ax0.get_legend_handles_labels()
    # print(handles,labels)
    ax0.legend(['ms0','ms-1','ms+1'])

    ax1.plot(t_list,Normalised_Osc)
    ax1.set_xlabel('time ns')
    ax1.set_ylabel('normalised population in ms0')
    ax1.set_ylim(0,1)
    ax1.grid(True)

Rabi_evolution(transition_driven='msm1')
Rabi_evolution(transition_driven='msp1')


plt.show()
