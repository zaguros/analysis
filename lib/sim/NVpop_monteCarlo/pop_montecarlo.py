def monteCarlo(init_state,rateMat,time_duration,dt,repetitions):

    total_points = np.ceil(time_duration/np.float(dt) + 1).astype(int)
    points = np.arange(total_points).astype(int)
    t_array = np.float(dt) * points
    pops_expanded = np.zeros([repetitions,total_points,num_states])
    pops = np.zeros([repetitions,total_points])
    dtRateMat = np.float(dt)*rateMat

    if np.max(dtRateMat) > 0.8:
        print 'dt not small enough! Probably gonna mess up'
        print 'dt should be no bigger than ', 0.8/np.max(rateMat)

    for jj in range(repetitions):
        rands = np.random.rand(total_points,num_states)
        evrands = np.random.rand(total_points)
        current_state = init_state
        rates = dtRateMat[current_state]
        for ii, rand in enumerate(rands):
            pops[jj,ii] = current_state
            decays = rand < rates
            events = np.where(decays)[0]
            n_events = np.size(events)
            if n_events == 1:
                current_state = events[0]
                rates = dtRateMat[current_state]
            elif n_events > 1:
                current_state = events[int(evrands[ii] * n_events)]
                rates = dtRateMat[current_state]

    for ind in range(num_states):
        pops_expanded[:,:,ind] = pops==ind

    return t_array,pops_expanded


def monteCarloNewSchool(int init_state,np.ndarray[double, ndim=2,negative_indices=False,mode='c'] rateMat,double time_duration,double dt,int repetitions):
    
    num_states = np.shape(rateMat)[0]
    total_points = np.ceil(time_duration/np.float(dt) + 1).astype(int)
    
    pops = np.empty([repetitions,total_points], dtype=np.int)
    rates = np.empty(num_states)
    dtRateMat = dt* rateMat
            
    state_array = np.arange(num_states)
    
    t_array = dt * np.arange(total_points).astype(int)
    
    pops_expanded = np.empty([repetitions,total_points,num_states])
    
    if np.max(dtRateMat) > 0.8:
        print 'dt not small enough! Probably gonna mess up'
        print 'dt should be no bigger than ', 0.8/np.max(rateMat)

    for jj in range(repetitions):
        current_state = init_state
        for x in range(num_states):
            rates[x] = dtRateMat[init_state,x]
        
        for ii in range(total_points):
            
            pops[jj,ii] = current_state
            
            for x in range(num_states): # Fisher Yates Shuffle
                y = int(rand()*randFactorN) # random int from 0 to num_states
                t_state = state_array[y]
                state_array[y] = state_array[x]
                state_array[x] = t_state
 
            for x in range(num_states): # Iterate through states testing for decay
                t_state = state_array[x]
                if rates[t_state] > (rand()*randFactor):
                    current_state = t_state
                    for x in range(num_states):
                        rates[x] = dtRateMat[t_state,x]
                    break
           
    for jj in range(repetitions):
        for ii in range(total_points):
            for ind in range(num_states):
                pops_expanded[jj,ii,ind] = pops[jj,ii]==ind

    return t_array,pops_expanded
