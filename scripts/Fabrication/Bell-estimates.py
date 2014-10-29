
from numpy import *
import pylab as plt
from matplotlib import rc

#import plots


distance=1.7 #km
dcounts_shot = 1.2e-6	
attanuation = 	8	#db per km
other_losses =	1	#db
detection_W = 	30	#ns
tailcounts=	14	#E-4
LDE_fidelity=	0.85	
rep_rate=	50	#khz
pts=100
tailcounts=np.linspace(0,20,pts)
fid=[]
time=[]
for i in np.arange(pts):
	P_NV_shot = tailcounts[i]*0.0001*10**(-(attanuation/2*distance+other_losses)/(10))
	P_NV_NV = 0.5*P_NV_shot**2
	P_NV_DC=2*P_NV_shot*dcounts_shot
	P_DC_DC=dcounts_shot**2
	F=(P_NV_NV*LDE_fidelity+P_NV_DC*1/4+P_DC_DC*1/4)/(P_NV_NV+P_NV_DC+P_DC_DC)
	T=1/(0.5*rep_rate*P_NV_NV*1000*60)
	fid.append(F)
	time.append(T)
'''
plt.plot(tailcounts,time,color='Crimson')
plt.plot(tailcounts,fid,color='RoyalBlue')
plt.xlabel ('Photon detection probability [10^-4]', fontsize = 12)
plt.ylabel ('Bell state Fidelity', fontsize = 12)   
#plt.ylim([0.75,0.85])
#plt.ylim ([0, 1])
'''
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
line1, = ax1.plot(tailcounts,fid,color='RoyalBlue',linewidth=2)
line2, = ax2.plot(tailcounts,time,color='Crimson',linewidth=2)
ax1.set_ylim([0.75,0.85])
ax2.set_yscale('log')
ax2.set_ylim([10,100])
ax1.set_ylabel ('Bell state Fidelity', fontsize = 16)  
ax2.set_ylabel ('Time per entanglement event [min]', fontsize = 16) 
ax1.set_xlabel ('Photon detection probability [$10^{-4}$]', fontsize = 16)
ax1.yaxis.label.set_color('RoyalBlue')
ax2.yaxis.label.set_color('Crimson')
ax2.set_yticks([10,50,100])
for label in ax2.get_yticklabels():
	label.set_fontsize(14)
for label in ax1.get_yticklabels():
	label.set_fontsize(14)	
for label in ax1.get_xticklabels():
	label.set_fontsize(14)

plt.show()




