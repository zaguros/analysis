from analysis.lib.nv import nvlevels
import numpy as np

##################
### ESR Spectra###
##################

Bz = 209.6 #B_field in Gauss
Ex = 2.5/2 #Ex is approximately half the strain splitting

es_spectrum =np.sort(nvlevels.get_ES(
                                                        E_field=[Ex,0.,0.], 
                                                        B_field=[0.,0.,Bz],
                                                        Ee0=-1.94,
                                                        transitions=False,
                                                        )[0])

ms0 = 0
msm1 = 2.8769 - Bz*2.8e-3
msp1 = 2.8769 + Bz*2.8e-3
gs_spectrum = [ms0, msm1, msp1]

print 'B_field= ', str(Bz)
print 'Strain is= ', str(Ex*2)

print 'excited state spectrum is ' + str(es_spectrum)
print 'excited state spectrum is ' + str(gs_spectrum)

### ESR transitions need to manually check the ordering ###
gs_ESR_msm = msm1
gs_ESR_msp = msp1

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

print 'es_ESR_Ey_Ep1 = ' + str(es_ESR_Ey_Ep1)
print 'es_ESR_Ey_Ep2 = ' + str(es_ESR_Ey_Ep2)

print 'es_ESR_Ex_A1 = ' + str(es_ESR_Ex_A1)
print 'es_ESR_Ex_A2 = ' + str(es_ESR_Ex_A2)

print 'es_ESR_Ey_A1 = ' + str(es_ESR_Ey_A1)
print 'es_ESR_Ey_A2 = ' + str(es_ESR_Ey_A2)

print 'es_ESR_Ex_Ey = ' + str(es_ESR_Ex_Ey)


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


################
###Sweeping B###
################

b_range=np.linspace(0,1000,100) #gauss

spectrum=np.zeros((6,))
for B_z in b_range:
	spectrum=np.vstack((spectrum,np.sort(nvlevels.get_ES(
														E_field=[Ex,0.,0.], 
														B_field=[0.,0.,B_z],
														Ee0=-1.94,
														transitions=False,
														)[0])))
spectrum=spectrum[1:]

pylab.figure() 
for i in range(6):
	pylab.plot(b_range,spectrum[:,i])

#################
###Sweeping E ###
#################
Ex_range = np.linspace(0,20,100)/2 #Ex is approx strain_splitting divided by 2

spectrum=np.zeros((6,))
for E_x in Ex_range:
    spectrum=np.vstack((spectrum,np.sort(nvlevels.get_ES(
                                                        E_field=[E_x,0.,0.], 
                                                        B_field=[0.,0.,Bz],
                                                        Ee0=-1.94,
                                                        transitions=False,
                                                        )[0])))
spectrum=spectrum[1:]

pylab.figure() 
for i in range(6):
    pylab.plot(Ex_range,spectrum[:,i])

