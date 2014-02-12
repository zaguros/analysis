from analysis.lib.nv import nvlevels
import numpy as np
b_range=np.linspace(0,1000,100) #gauss
Ex=8/2. #Ex is approx strain_splitting divided by 2
spectrum=np.zeros((6,))
for Bz in b_range:
	spectrum=np.vstack((spectrum,np.sort(nvlevels.get_ES(
														E_field=[Ex,0.,0.], 
														B_field=[0.,0.,Bz],
														Ee0=-1.94,
														transitions=False,
														)[0])))
spectrum=spectrum[1:]

for i in range(6):
	plot(b_range,spectrum[:,i])