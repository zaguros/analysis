from numpy import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.stats import norm
from analysis.lib.nv import NV_generator as Gen
reload(Gen)




def nr_spins_in_radius(Carb_Conc = 0.011, N = 20, radius = 2, hyperfine_difference = 1e3):
	
	Ap_sorted, Ao_sorted,r_sorted = Gen.Generate_NV(Carb_Conc=Carb_Conc,N_NV=1,N=N,do_sphere = True)
	print max(r_sorted)
	for jj, r in enumerate(r_sorted):
		if r>radius*1e-9:
			break
	print 'Number of C13 spins within '+str(radius)+' nm = ' + str(jj)

	m = 0
	for Apar in Ap_sorted[0:jj]: # only take C13 spins within 2 nm from NV
		Atemp = abs(Ap_sorted[0:jj]-Apar) #take difference parallel hyperfine coupling wrt to other C13 spins
		n = 0
		for x in Atemp:
			if x < hyperfine_difference:
				n = n+1
		if n==1: # only difference with itself should be 0
			m = m+1
			
	print 'Number of C13 spins with parallel hyperfine coupling of more than ' +str(hyperfine_difference*1e-3) + ' kHz difference to other spins in radius: '+str(m)

def plot_r_vs_number(Carb_Conc = 0.011, N = 25):
	
	Ap, Ao,r = Gen.Generate_NV(Carb_Conc=Carb_Conc,N_NV=1,N=N,do_sphere=True)

	r_sorted = sort(r)[0:len(r)/2.]
	fig,ax = plt.subplots()
	ii = range(len(r_sorted))
	ax.plot(ii,r_sorted)
	ax.set_xlabel('Carbon number')
	ax.set_ylabel('Distance to NV (nm)')
	plt.show()

	# fig,ax = plt.subplots()
	# ii = range(len(Ap_sorted))
	# ax.plot(ii[0:jj],Ap_sorted[0:jj])
	# plt.show()