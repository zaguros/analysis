from analysis.lib.nv import nvlevels
reload(nvlevels)
import numpy as np
from matplotlib import pyplot as plt

save_path = 'H:/My Documents/Cavities/Simulations/nv_levels/'

def plot_ES_b_dependence(strain_splitting=2.):
	b_range=np.linspace(0,1000,100) #gauss
	Ex=strain_splitting/2. #Ex is approx strain_splitting divided by 2
	spectrum=np.zeros((6,))
	#plt.close()
	for Bz in b_range:
		spectrum=np.vstack((spectrum,np.sort(nvlevels.get_ES(
															E_field=[Ex,0.,0.], 
															B_field=[0.,0.,Bz],
															Ee0=-1.94,
															transitions=False,
															)[0])))
	spectrum=spectrum[1:]

	for i in range(6):
		plt.plot(b_range,spectrum[:,i])
		plt.savefig(save_path+'ES_B_dependence_Ex_'+str(Ex)+'.png')

def plot_GS_b_dependence(strain_splitting=2.):
	b_range=np.linspace(0,1000,100) #gauss
	Ex=strain_splitting/2. #Ex is approx strain_splitting divided by 2
	spectrum=np.zeros((3,))
	plt.close()
	for Bz in b_range:
		spectrum=np.vstack((spectrum,np.sort(nvlevels.get_GS(
															E_field=[Ex,0.,0.], 
															B_field=[0.,0.,Bz],
															)[0])))
	spectrum=spectrum[1:]

	for i in range(3):
		plt.plot(b_range,spectrum[:,i])
		plt.savefig(save_path+'GS_B_dependence.png')

def plot_ES_e_dependence(Bz=20.):
	e_range=np.linspace(0,10,100) #gauss
	spectrum=np.zeros((6,))
	# plt.close()
	plt.figure()
	for Ex in e_range:
		spectrum=np.vstack((spectrum,np.sort(nvlevels.get_ES(
															E_field=[Ex,0.,0.], 
															B_field=[0.,0.,Bz],
															Ee0=-1.94,
															transitions=False,
															)[0])))
	spectrum=spectrum[1:]


	for i in range(6):
		plt.plot(e_range,spectrum[:,i])
		plt.title('Excited state levels E dependence, Bz = '+str(Bz))
		plt.xlabel('Strain E (GHz)')
		plt.ylabel('relative energy (GHz)')
		plt.savefig(save_path+'ES_E_dependence_Bz_'+str(Bz)+'.png')

def plot_transitions_b_dependence(strain_splitting=4.11):
	b_range=np.linspace(250,500,50) #gauss
	Ex=strain_splitting/2. #Ex is approx strain_splitting divided by 2
	no_transitions=8
	spectrum=np.zeros((no_transitions,))
	# plt.close()
	# plt.figure()
	for Bz in b_range:
		spectrum=np.vstack((spectrum,nvlevels.get_optical_transitions(show_E_transitions=True,
																show_A_transitions=True,
																show_FB_E_transitions=False,	 
                            									show_FB_A_transitions=False, 
                            									show_E_prime_flip_transitions=False,
															E_field=[Ex,0.,0.], 
															B_field=[0,0.,Bz],
															Ee0=-1.94,
															)))
	spectrum=spectrum[1:]

	for i in range(no_transitions):
		plt.title('Transitions B dependence, Strain = '+str(Ex))
		plt.xlabel('Magnetic field Bz (G)')
		plt.ylabel('relative energy difference (GHz)')
		# if i < 4:
		style = 3
		# else: style =0.5
		plt.plot(b_range,spectrum[:,i],linewidth = style)
	plt.axvline(x=400, ymin=-10, ymax = 8, linewidth=2, color='k')
	plt.savefig(save_path+'Transitions_B_dependence_strain_'+str(Ex)+'.png')

def plot_transitions_bx_dependence(Bz=300,strain_splitting=2.5):
	bx_range=np.linspace(0,500,100) #gauss
	Ex=strain_splitting/2. #Ex is approx strain_splitting divided by 2
	no_transitions=18
	spectrum=np.zeros((no_transitions,))
	plt.close()
	for Bx in bx_range:
		spectrum=np.vstack((spectrum,nvlevels.get_optical_transitions(
															E_field=[Ex,0.,0.], 
															B_field=[Bx,0.,Bz],
															Ee0=-1.94,
															show_FB_E_transitions=True
															)))
	spectrum=spectrum[1:]

	for i in range(no_transitions):

		plt.plot(bx_range,spectrum[:,i])
		plt.savefig(save_path+'Transitions_Bx_dependence_Ex_'+str(Ex)+'_Bz_'+str(Bz)+'.png')

def plot_transitions_E_dependence(Bz=300.):
	e_range=np.linspace(0,10,100) #gauss
	no_transitions=18

	spectrum=np.zeros((no_transitions,))
	plt.figure()
	for Ex in e_range:
		spectrum=np.vstack((spectrum,nvlevels.get_optical_transitions(
															show_E_transitions=True,
															show_A_transitions=True,
															show_FB_E_transitions=True, 
                        									show_FB_A_transitions=True, 
                        									show_E_prime_flip_transitions=True,
															E_field=[Ex,0.,0.], 
															B_field=[0,0.,Bz],
															Ee0=-1.94,
															)))
	spectrum=spectrum[1:]
	# print spectrum
	# print len(spectrum)
	for i in range(no_transitions):
		if i < 8:
			style = 2
		else: style =0.5
		plt.plot(e_range,spectrum[:,i],linewidth = style)
		plt.title('Transitions E dependence, Bz = '+str(Bz))
		plt.xlabel('Strain E (GHz)')
		plt.ylabel('relative energy difference (GHz)')
		plt.savefig(save_path+'Transitions_E_dependence_Bz_'+str(Bz)+'.png')

def plot_transitions_in_laserscan_Ex_Ey(Ex,Ey,height=1,B_field = [0.,0.,300.]):
	x2,y2 = nvlevels.get_transitions_ExEy_plottable(Ex,Ey,height,B_field=B_field)
	plt.plot(x2,y2)

def plot_p_ms0_vs_E(Bz = 0):
	e_range=np.linspace(0,20,100) #GHz
	p_Ex_ms0 = []
	p_Ey_ms0 = []
	p_Eprimexy_ms0 = []
	for Ex in e_range:
		p_Ex_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  3))
		if nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2) > nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0):
			p_Ey_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2))
			p_Eprimexy_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0))
		else:
			p_Eprimexy_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2))
			p_Ey_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0))
	fig=plt.figure()
	ax1 = fig.add_subplot(111)
	lns1 = ax1.plot(e_range,p_Ex_ms0,color = 'red',label = 'Ex')
	ax2 = ax1.twinx()
	lns2 = ax2.plot(e_range,p_Ey_ms0,color='blue',label = 'Ey')
	plt.title('ms=0 probability')#, Bz = '+str(Bz))
	plt.xlabel('Strain E (GHz)')
	ax1.set_ylabel('p(ms=0) for Ex', color='r')
	ax2.set_ylabel('p(ms=0) for Ey', color='b')
	lns = lns1+lns2
	labs = [l.get_label() for l in lns]
	ax1.legend(lns, labs, loc=3)
	# ax1.axis([min(e_range),max(e_range),0.99,1])
	for tl in ax2.get_yticklabels():
	    tl.set_color('b')
	for tl in ax1.get_yticklabels():
	    tl.set_color('r')


def plot_p_ms0_vs_E_incl_mixing(Bz = 400):
	e_range=np.linspace(0,20,100) #GHz
	color_Ex = [(1,0,0),(0.8,0,0),(0.6,0,0),(0.4,0,0),(0.2,0,0)]
	color_Ey = [(0,0,1),(0,0,0.8),(0,0,0.6),(0,0,0.4),(0,0,0.2)]
	T_list = [2,3,4,5,6]
	fig=plt.figure()
	ax1 = fig.add_subplot(111)	
	ax2 = ax1.twinx()
	labs = []
	for i,T in enumerate(T_list):	
		p_Ex_ms0 = []
		p_Ey_ms0 = []
		p_Eprimexy_ms0 = []
		p_ExMix = []
		p_EyMix = []
		for Ex in e_range:
			p_Ex_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  3))
			if nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2) > nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0):
				p_Ey_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2))
				p_Eprimexy_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0))
			else:
				p_Eprimexy_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2))
				p_Ey_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0))
		

		pmix = nvlevels.mixing_probability(T)
		for ii in range(len(p_Ex_ms0)):
			p_ExMix.append(p_Ex_ms0[ii]*(1-pmix)+p_Ey_ms0[ii]*pmix)
			p_EyMix.append(p_Ey_ms0[ii]*(1-pmix)+p_Ex_ms0[ii]*pmix)


	# fig=plt.figure()
	# ax1 = fig.add_subplot(111)
		lns1 = ax1.plot(e_range,p_ExMix,color = color_Ex[i],label = 'Ex')
		lns2 = ax2.plot(e_range,p_EyMix,color=color_Ey[i],label = 'Ey')
	plt.title('ms=0 probability')#, Bz = '+str(Bz))
	plt.xlabel('Strain E (GHz)')
	ax1.set_ylabel('p(ms=0) for Ex', color='r')
	ax2.set_ylabel('p(ms=0) for Ey', color='b')
		# lns = lns1+lns2
		# labs = [l.get_label() for l in lns]
		# ax1.legend(lns, labs, loc=3)
		# ax1.axis([min(e_range),max(e_range),0.99,1])
	for tl in ax2.get_yticklabels():
	    tl.set_color('b')
	for tl in ax1.get_yticklabels():
	    tl.set_color('r')

def plot_p_ms0_vs_Bx_incl_mixing(Bx_range= np.linspace(0,100,50)):
	Ex = 2.155
	Bz=402.7
	color_Ex = [(1,0,0),(0.8,0,0),(0.6,0,0),(0.4,0,0),(0.2,0,0)]
	color_Ey = [(0,0,1),(0,0,0.8),(0,0,0.6),(0,0,0.4),(0,0,0.2)]
	T_list = [2,3,4,5,6]
	fig=plt.figure()
	ax1 = fig.add_subplot(111)	
	ax2 = ax1.twinx()
	labs = []
	for i,T in enumerate(T_list):	
		p_Ex_ms0 = []
		p_Ey_ms0 = []
		p_Eprimexy_ms0 = []
		p_ExMix = []
		p_EyMix = []
		for Bx in Bx_range:
			p_Ex_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  3,Bx=Bx))
			if nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2,Bx=Bx) > nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1,Bx=Bx)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0,Bx=Bx):
				p_Ey_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2,Bx=Bx))
				p_Eprimexy_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1,Bx=Bx)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0,Bx=Bx))
			else:
				p_Eprimexy_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2,Bx=Bx))
				p_Ey_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1,Bx=Bx)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0,Bx=Bx))
		

		pmix = nvlevels.mixing_probability(T)
		for ii in range(len(p_Ex_ms0)):
			p_ExMix.append(p_Ex_ms0[ii]*(1-pmix)+p_Ey_ms0[ii]*pmix)
			p_EyMix.append(p_Ey_ms0[ii]*(1-pmix)+p_Ex_ms0[ii]*pmix)

		lns1 = ax1.plot(Bx_range,p_ExMix,color = color_Ex[i],label = 'Ex')
		lns2 = ax2.plot(Bx_range,p_EyMix,color=color_Ey[i],label = 'Ey')
	plt.title('ms=0 probability')#, Bz = '+str(Bz))
	plt.xlabel('Bx (G)')
	ax1.set_ylabel('p(ms=0) for Ex', color='r')
	ax2.set_ylabel('p(ms=0) for Ey', color='b')

def plot_p_ms0_vs_Bz_incl_mixing(Bz_range= np.linspace(0,800,400)):
	Ex = 2.155
	# Bz=402.7
	color_Ex = [(1,0,0)]#,(0.8,0,0),(0.6,0,0),(0.4,0,0),(0.2,0,0)]
	color_Ey = [(0,0,1)]#,(0,0,0.8),(0,0,0.6),(0,0,0.4),(0,0,0.2)]
	T_list = [4]#[2,3,4,5,6]
	fig=plt.figure()
	ax1 = fig.add_subplot(111)	
	ax2 = ax1.twinx()
	labs = []
	for i,T in enumerate(T_list):	
		p_Ex_ms0 = []
		p_Ey_ms0 = []
		p_Eprimexy_ms0 = []
		p_ExMix = []
		p_EyMix = []
		for Bz in Bz_range:
			p_Ex_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  3))
			if nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2) > nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0):
				p_Ey_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2))
				p_Eprimexy_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0))
			else:
				p_Eprimexy_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  2))
				p_Ey_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  1)+nvlevels.get_ms0_fraction_incl_B(2*Ex,Bz,  0))
		

		pmix = nvlevels.mixing_probability(T)
		for ii in range(len(p_Ex_ms0)):
			p_ExMix.append(p_Ex_ms0[ii]*(1-pmix)+p_Ey_ms0[ii]*pmix)
			p_EyMix.append(p_Ey_ms0[ii]*(1-pmix)+p_Ex_ms0[ii]*pmix)

		lns1 = ax1.plot(Bz_range,p_ExMix,color = color_Ex[i],label = 'Ex')
		lns2 = ax2.plot(Bz_range,p_EyMix,color=color_Ey[i],label = 'Ey')
	plt.title('ms=0 probability')#, Bz = '+str(Bz))
	plt.xlabel('Bz (G)')
	ax1.set_ylabel('p(ms=0) for Ex', color='r')
	ax2.set_ylabel('p(ms=0) for Ey', color='b')


def plot_p_ms0_vs_E_vs_B():
	e_range=np.linspace(0,10,100) #GHz
	b_list = [0,300,400,500,600]
	color_Ex = [(1,0,0),(0.8,0,0),(0.6,0,0),(0.4,0,0),(0.2,0,0)]
	color_Ey = [(0,0,1),(0,0,0.8),(0,0,0.6),(0,0,0.4),(0,0,0.2)]
	fig=plt.figure()
	ax1 = fig.add_subplot(111)	
	ax2 = ax1.twinx()
	labs = []
	for i,Bz in enumerate(b_list):
		p_Ex_ms0 = []
		p_Ey_ms0 = []
		p_Eprimexy_ms0 = []
		for Ex in e_range:
			p_Ex_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex, Bz, 3))
			if nvlevels.get_ms0_fraction_incl_B(2*Ex, Bz, 2) > nvlevels.get_ms0_fraction_incl_B(2*Ex, Bz, 1)+nvlevels.get_ms0_fraction_incl_B(2*Ex, Bz, 0):
				p_Ey_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex, Bz, 2))
				p_Eprimexy_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex, Bz, 1)+nvlevels.get_ms0_fraction_incl_B(2*Ex, Bz, 0))
			else:
				p_Eprimexy_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex, Bz, 2))
				p_Ey_ms0.append(nvlevels.get_ms0_fraction_incl_B(2*Ex, Bz, 1)+nvlevels.get_ms0_fraction_incl_B(2*Ex, Bz, 0))

		lns1 = ax1.plot(e_range,p_Ex_ms0,color = color_Ex[i],label = 'Ex')
		lns2 = ax2.plot(e_range,p_Ey_ms0,color=color_Ey[i],label = 'Ey')
	plt.title('ms=0 probability')#, Bz = '+str(Bz))
	plt.xlabel('Strain E (GHz)')
	ax1.set_ylabel('p(ms=0) for Ex', color='r')
	ax2.set_ylabel('p(ms=0) for Ey', color='b')
	# 	lns = lns1+lns2
	# 	labs.append = [l.get_label() for l in lns]
	# ax1.legend(lns, labs, loc=3)
	ax1.axis([min(e_range),max(e_range),0.99,1])
	for tl in ax2.get_yticklabels():
	    tl.set_color('b')
	for tl in ax1.get_yticklabels():
	    tl.set_color('r')
	# color_Ex = [(1,0,0),(1,0.1,0),(1,0.2,0),(1,0.3,0),(1,0.4,0)]
	# color_Ey = [(0,0,1),(0,0.1,1),(0,0.2,1),(0,0.3,1),(0,0.4,1)]