import cython
import numpy as np
cimport numpy as np
from libc.stdlib cimport rand, RAND_MAX

@cython.boundscheck(False)

def monteCarlo(np.ndarray[int, ndim=1,negative_indices=False,mode='c'] init_state,np.ndarray[double, ndim=2,negative_indices=False,mode='c'] rateMat,double time_duration,double dt,int repetitions):
    
    cdef int init_state_per_rep
    if np.size(init_state) == 1:
        init_state_per_rep = 0
    elif np.size(init_state) == repetitions:
        init_state_per_rep = 1
    else:
        print "Incorrect init_state size!"

    cdef int num_states = np.shape(rateMat)[0]
    cdef int total_points = np.ceil(time_duration/np.float(dt) + 1).astype(int)
    cdef int current_state
    
    cdef extern from "stdlib.h":
        double drand48()

    cdef np.ndarray[double, ndim=1,negative_indices=False,mode='c'] t_array = dt * np.arange(total_points).astype(int)
    cdef np.ndarray[np.int_t, ndim=2,negative_indices=False,mode='c'] pops = np.empty([repetitions,total_points], dtype=np.int)
    
    cdef np.ndarray[double, ndim=1,negative_indices=False,mode='c'] cumBranching = np.empty(num_states)
    cdef np.ndarray[double, ndim=1,negative_indices=False,mode='c'] probs = dt*np.sum(rateMat,axis=1).astype(float)
    probs[probs == 0.0] = -1.0 # To avoid tedious divide by zero errors
    cdef np.ndarray[double, ndim=2,negative_indices=False,mode='c'] branchingMat = np.cumsum(dt*rateMat,axis=1)/np.tile(probs,[num_states,1]).T

    cdef np.ndarray[double, ndim=3,negative_indices=False,mode='c'] pops_expanded = np.empty([repetitions,total_points,num_states])
    
    cdef double prob, branching_rand

    if np.max(probs) > 0.8:
        print 'dt not small enough! Probably gonna mess up'
        print 'dt should be no bigger than ', 0.8/np.max(probs)

    for jj in range(repetitions):
        if init_state_per_rep == 1:
            current_state = init_state[jj]
        else:
            current_state = init_state[0]

        prob = probs[current_state]
        for x in range(num_states): # Stupid way to slice a matrix.
            cumBranching[x] = branchingMat[current_state,x]

                
        for ii in range(total_points):
            
            if prob == -1.0: # Pretty boring if never going to decay from this state
                pops[jj,ii:] = current_state
                break
            else:
                pops[jj,ii] = current_state
                if drand48() < prob:
                    branching_rand = drand48()
                    for y in range(num_states): # Iterate through states testing for decay
                        if cumBranching[y] > branching_rand:
                            current_state = y
                            prob = probs[current_state]
                            for x in range(num_states):
                                cumBranching[x] = branchingMat[current_state,x]
                            break
           
    for jj in range(repetitions):
        for ii in range(total_points):
            for ind in range(num_states):
                pops_expanded[jj,ii,ind] = pops[jj,ii]==ind

    return t_array,pops_expanded,pops

def monteCarlo_old(int init_state,np.ndarray[double, ndim=2,negative_indices=False,mode='c'] rateMat,double time_duration,double dt,int repetitions):
    
    cdef int num_states = np.shape(rateMat)[0]
    cdef int total_points = np.ceil(time_duration/np.float(dt) + 1).astype(int)
    cdef int current_state, y, t_state
    
    cdef np.ndarray[np.int_t, ndim=2,negative_indices=False,mode='c'] pops = np.empty([repetitions,total_points], dtype=np.int)
    cdef np.ndarray[double, ndim=1,negative_indices=False,mode='c'] rates = np.empty(num_states)
    cdef np.ndarray[double, ndim=2,negative_indices=False,mode='c'] dtRateMat = dt* rateMat
            
    cdef np.ndarray[np.int_t, ndim=1,negative_indices=False,mode='c'] state_array = np.arange(num_states)
    
    cdef np.ndarray[double, ndim=1,negative_indices=False,mode='c'] t_array = dt * np.arange(total_points).astype(int)
    
    cdef np.ndarray[double, ndim=3,negative_indices=False,mode='c'] pops_expanded = np.empty([repetitions,total_points,num_states])
    
    cdef double randFactorN = 1.0/(RAND_MAX*np.double(num_states))
    cdef double randFactor = 1.0/RAND_MAX
    
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
