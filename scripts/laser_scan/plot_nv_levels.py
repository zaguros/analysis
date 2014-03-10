from analysis.lib.nv import nvlevels
import numpy as np
import pylab
from matplotlib import pyplot as plt
from analysis.lib.tools import plot

def get_NV_states(B_field = [0., 0., 0.], Ex = 0., ZFS = 2.8769):

    es_spectrum =np.sort(nvlevels.get_ES(
                                                        E_field=[Ex,0.,0.], 
                                                        B_field=[B_field[0],B_field[1],B_field[2]],
                                                        Ee0=-1.94,
                                                        transitions=False,
                                                        )[0])

    print 'B_field = ', str(B_field)
    print 'Strain is = ', str(Ex*2) 
    print ''

    ms0_energy = 0
    msm_energy = ((ZFS - B_field[2]*2.8e-3)**2 + (B_field[0]*2.8e-3)**2 + (B_field[1]*2.8e-3)**2)**0.5 
    msp_energy = ((ZFS + B_field[2]*2.8e-3)**2 + (B_field[0]*2.8e-3)**2 + (B_field[1]*2.8e-3)**2)**0.5
    gs_spectrum = [ms0_energy, msm_energy, msp_energy]

    print 'ground state spectrum is ' + str(gs_spectrum) 
    print ''    
    print 'excited state spectrum is ' + str(es_spectrum)
    print ''

    return gs_spectrum, es_spectrum


def plot_excited_state_vs_B(Ex = 0, b_range = np.linspace(0,1000,100), b_direction = 'Z'):
    ''' Plots the excited state levels as function of magnetic field over the b_range (in Gauss) and 
    with orientation b_direction'''
    spectrum=np.zeros((6,))
    for B in b_range:
        if b_direction == 'Z':
            spectrum=np.vstack((spectrum,np.sort(nvlevels.get_ES(
                                                            E_field=[Ex,0.,0.], 
                                                            B_field=[0.,0.,B],
                                                            Ee0=-1.94,
                                                            transitions=False,
                                                            )[0])))
        elif b_direction == 'X':
            spectrum=np.vstack((spectrum,np.sort(nvlevels.get_ES(
                                                            E_field=[Ex,0.,0.], 
                                                            B_field=[B,0.,0.],
                                                            Ee0=-1.94,
                                                            transitions=False,
                                                            )[0])))
        elif b_direction == 'Y':
            spectrum=np.vstack((spectrum,np.sort(nvlevels.get_ES(
                                                            E_field=[Ex,0.,0.], 
                                                            B_field=[0.,B,0.],
                                                            Ee0=-1.94,
                                                            transitions=False,
                                                            )[0])))
        else: 
            print 'b_direction must be "X", "Y" or "Z"'
    
    spectrum=spectrum[1:]

    pylab.figure() 
    for i in range(6):
        pylab.plot(b_range,spectrum[:,i])

def plot_excited_state_vs_strain(B = [0., 0., 0.], Ex_range = np.linspace(0,20,100)/2):
    ''' B_field in gauss and as a vector [B_x, B_y, B_z], 
    Ex is approx strain_splitting divided by 2
     '''
    spectrum=np.zeros((6,))
    for E_x in Ex_range:
        spectrum=np.vstack((spectrum,np.sort(nvlevels.get_ES(
                                                            E_field=[E_x,0.,0.], 
                                                            B_field=[B[0],B[1],B[2]],
                                                            Ee0=-1.94,
                                                            transitions=False,
                                                            )[0])))
    spectrum=spectrum[1:]
    pylab.figure() 
    for i in range(6):
        pylab.plot(Ex_range,spectrum[:,i]) 


def get_ESR_transitions(B_field = [0., 0., 0.], Ex = 0., ZFS = 2.8769):

    gs_spectrum, es_spectrum = get_NV_states(B_field = B_field, Ex = Ex, ZFS = ZFS)

    gs_ESR_msm = gs_spectrum[1] 
    gs_ESR_msp = gs_spectrum[2]

    es_ESR_Ex_Ep1 = abs(es_spectrum[0]-es_spectrum[3])
    es_ESR_Ex_Ep2 = abs(es_spectrum[1]-es_spectrum[3])

    es_ESR_Ey_Ep1 = abs(es_spectrum[0]-es_spectrum[2])
    es_ESR_Ey_Ep2 = abs(es_spectrum[1]-es_spectrum[2])

    es_ESR_Ex_A1 = abs(es_spectrum[4]-es_spectrum[3])
    es_ESR_Ex_A2 = abs(es_spectrum[5]-es_spectrum[3])

    es_ESR_Ey_A1 = abs(es_spectrum[4]-es_spectrum[2])
    es_ESR_Ey_A2 = abs(es_spectrum[5]-es_spectrum[2])

    es_ESR_Ex_Ey = abs(es_spectrum[2]-es_spectrum[3])

    print 'es_ESR_Ex_Ep1 = ' + str(es_ESR_Ex_Ep1)
    print 'es_ESR_Ex_Ep2 = ' + str(es_ESR_Ex_Ep2)
    print ''   
    print 'es_ESR_Ey_Ep1 = ' + str(es_ESR_Ey_Ep1)
    print 'es_ESR_Ey_Ep2 = ' + str(es_ESR_Ey_Ep2)
    print ''   
    print 'es_ESR_Ex_A1 = ' + str(es_ESR_Ex_A1)
    print 'es_ESR_Ex_A2 = ' + str(es_ESR_Ex_A2)
    print ''   
    print 'es_ESR_Ey_A1 = ' + str(es_ESR_Ey_A1)
    print 'es_ESR_Ey_A2 = ' + str(es_ESR_Ey_A2)
    print ''   
    print 'es_ESR_Ex_Ey = ' + str(es_ESR_Ex_Ey)


def get_optical_transitions(B_field = [0., 0., 0.], Ex = 0., ZFS = 2.8769):
    pass


'''
###Optical Spectrum###
## For now manusally check the ordering
Ex_experimental = 68.7 - es_spectrum[3]


trans_Ex =  es_spectrum[3] +Ex_experimental
trans_Ey =  es_spectrum[2] +Ex_experimental

trans_msp_Ep1 =  es_spectrum[0] - msp1 +Ex_experimental
trans_msp_Ep2 =  es_spectrum[1] - msp1 +Ex_experimental

trans_msm_Ep1 =  es_spectrum[0] - msm1 +Ex_experimental
trans_msm_Ep2 =  es_spectrum[1] - msm1 +Ex_experimental

trans_msp_A1 =  es_spectrum[4] - msp1 +Ex_experimental
trans_msp_A2 =  es_spectrum[5] - msp1 +Ex_experimental

trans_msm_A1 =  es_spectrum[4] - msm1 +Ex_experimental
trans_msm_A2 =  es_spectrum[5] - msm1 +Ex_experimental


print 'Optical transitions'
print 'Ex: ' + str(trans_Ex)
print 'Ey: ' + str(trans_Ey)

print 'trans_msp_Ep1: ' + str(trans_msp_Ep1)
print 'trans_msp_Ep2: ' + str(trans_msp_Ep2)

print 'trans_msm_Ep1: ' + str(trans_msm_Ep1)
print 'trans_msm_Ep2: ' + str(trans_msm_Ep2)

print 'trans_msp_A1: ' + str(trans_msp_A1)
print 'trans_msp_A2: ' + str(trans_msp_A2)

print 'trans_msm_A1: ' + str(trans_msm_A1)
print 'trans_msm_A2: ' + str(trans_msm_A2)
'''