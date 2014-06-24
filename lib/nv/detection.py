
def get_Ey_SSRO_stats_from_saturation(strain_splitting, saturation_counts):
	dead_time = 50e-9  #due to combined apd/counting lifetime
	NV_effective_lifetime = 12e-9 + 1/100*150e-9 #ms0 cycling + singlet cycling

	real_countrates=saturation_counts*1/(1 - (dead_time*saturation_counts))

	eta=real_countrates*NV_effective_lifetime
	return get_RO_stats_from_detection_efficiency(strain_splitting, eta)

def get_Ey_SSRO_stats_from_detection_efficiency(strain_splitting, eta):
	
	pflip = get_pflip(strain_splitting)
	F0 = eta/(pflip+eta-pflip*eta) # = 1-Sum[(1-eta)^i *pflip*(1-pflip)^(i-1) ,{i,1,Infinity}]; 
								   # see also page 135 of bas' master thesis
	navg = eta/pflip
	return F0, eta

def get_pflip(strain_splitting):
	p0 = 0.012 #spin flip probability at zero strain
	S0 = 15.   #GHz, the spoint where Ey becomes degenerate with E12
	pflip = p0*(S0/(S0 - strain_splitting)) #phenomological model for the spin flip probability 
											#TODO: verify with nvlevels eigenstate model
	return pflip