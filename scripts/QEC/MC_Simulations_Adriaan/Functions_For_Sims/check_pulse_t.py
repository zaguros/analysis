import numpy as np
def check_pulse_t(tlist, n_pulses, max_gate_time):
    """
    compares time it applies to do a gate to a certain max value and passes on a 2D array
    containing the corresponding resonance times and number of pules
    """
    gate_time = 4*np.multiply(tlist,n_pulses) #n_pulses is pulse blokken
    ind = np.where(gate_time<max_gate_time)

    resonances = np.column_stack((tlist[ind],n_pulses[ind] ))
    return resonances
