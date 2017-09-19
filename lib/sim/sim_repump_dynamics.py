import analysis.lib.sim.electron_nuclear_sim as en_sim
reload(en_sim)
from scipy.stats import norm
%matplotlib inline


def simulate_repump(time_duration = 50,dt=0.1,repetitions = 6000):
		
	time_duration = 1000
	dt = 2
	repetitions = 4000

	total_points = np.ceil(time_duration/np.float(dt) + 1).astype(int)
	points = np.arange(total_points).astype(int)
	t_array = dt * points

	num_states = 6
	pops = np.zeros([repetitions,total_points,num_states])

	init_state = 0
	rateMat = repumping_rateMat()

	pops = np.zeros([repetitions,total_points])
	dtRateMat = dt*rateMat

	if np.max(dtRateMat) > 0.8:
	    print 'dt not small enough! Probably gonna mess up'
	    print 'dt should be no bigger than ', 0.8/np.max(rateMat)
	st = time.time()
	for jj in range(repetitions):
	    rands = np.random.rand(total_points,num_states)
	    evrands = np.random.rand(total_points)
	    current_state = init_state
	    for ii in points:
	        decays = rands[ii] < dtRateMat[current_state]
	        events = np.where(decays)[0]
	        n_events = np.size(events)
	        if n_events == 1:
	            current_state = events[0]
	        elif n_events > 1:
	            current_state = events[int(evrands[ii] * n_events)]
	        pops[jj,ii] = current_state

	print time.time() - st

	populations = np.zeros([total_points,num_states])

	for ind in range(num_states):
	    populations[:,ind] = np.mean(pops==ind,axis =0)
	populations[:,2] = 1-populations[:,2]
	plt.plot(t_array,populations)
	plt.ylim([0,1])
	plt.show()



def simulate_double_pulsing(time_duration = 50,pulse_center = 5,dt = 0.1,FWHMs = np.array([0.5,1.6,2.3]),repetitions = 6000):
	
	mus = FWHMs/2.35482004503

	initial_excited_state_pop = 1.0
	NV_lifetime = 12.0
	decay_prob_per_step = dt/NV_lifetime

	init_state = en_sim.ket0.full()
	zero_state = en_sim.ket0.full()
	total_points = np.ceil(time_duration/np.float(dt) + 1).astype(int)
	points = np.arange(total_points).astype(int)
	t_array = dt * points
	excited_pops = np.zeros([np.shape(mus)[0],total_points])
	emitted_photons = np.zeros([np.shape(mus)[0],total_points])
	double_pulse_probs = np.zeros(np.shape(mus)[0])
	p_photon_during_pulse = np.zeros(np.shape(mus)[0])
	p_photon_after_pulse = np.zeros(np.shape(mus)[0])

	for kk,mu in enumerate(mus):
	    pulse_cut = pulse_center + 3*mu
	    pulse_cut_ind = np.argmin(np.abs(t_array - pulse_cut))

	    unitary = make_optical_pulse_unitary(t_array,pulse_center,mu,np.arcsin(np.sqrt(initial_excited_state_pop)))
	    
	    excited_pop = np.zeros([repetitions,total_points])
	    emitted_photon = np.zeros([repetitions,total_points])
	    for jj in range(repetitions):
	        decay_chance_for_step = np.random.rand(total_points)
	        bright_state_chance_for_step = np.random.rand(total_points)
	        state = init_state
	        for ii in points:
	            state =  np.dot(unitary[ii], state)
	            excited_pop[jj,ii] = np.abs(state[1])**2
	            if decay_chance_for_step[ii] < decay_prob_per_step:
	                state = zero_state
	                if bright_state_chance_for_step[ii] < excited_pop[jj,ii]:
	                    emitted_photon[jj,ii] = 1

	    excited_pops[kk] = np.mean(excited_pop,axis=0)
	    emitted_photons[kk] = np.mean(emitted_photon,axis=0)
	    photon_during_pulse = np.sum(emitted_photon[:,0:pulse_cut_ind],axis=1).astype(int)
	    photon_after_pulse = np.sum(emitted_photon[:,pulse_cut_ind:],axis=1).astype(int)
	    
	    p_photon_during_pulse[kk] = np.mean(photon_during_pulse)
	    p_photon_after_pulse[kk] = np.mean(photon_after_pulse)
	    double_pulse_probs[kk] = np.mean(np.logical_and(photon_during_pulse,photon_after_pulse))

	cond_prob_double_pulse = double_pulse_probs/p_photon_after_pulse


def make_optical_pulse_unitary(t_array,pulse_center,mu,pulse_theta):
	pulse_EField_shape = np.sqrt(norm.pdf(t_array, pulse_center, mu))
    # Need to correct power to ensure full rabi, could also do this easily for arbitrary pulse shapes
    EFieldFactor = 2*pulse_theta/(2**(0.75) * (np.pi)**(0.25) * np.sqrt(mu))
    instaneous_rabi_freq = pulse_EField_shape*EFieldFactor
    unitary = np.zeros([total_points,2,2],dtype=complex)
    for ii in points:
        unitary[ii] = en_sim.spin_y_rotation(instaneous_rabi_freq[ii] * dt).full()
    return unitary
