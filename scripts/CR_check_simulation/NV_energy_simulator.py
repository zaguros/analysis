import numpy as np
import math
import matplotlib as plt
import random
import time
import types
import os
import msvcrt

class NV_energy_simulator(Instrument):
    def __init__(self, name, plot_name=''):
        Instrument.__init__(self, name)
        # Gaussian related
        peak_width_gate = 1.0*10**6;
        peak_width_newfocus = 1.0*10**6; #GHz
        peak_width_yellow = 0.5*10**6;
        peak_counts_cr_check = 35.0;
        peak_counts_repump = 25;
        resonance_start_gate = 57.0*10**6;
        resonance_start_newfocus = 45.0*10**6; #GHz
        resonance_start_yellow = 17.0*10**6;

        # Time related
        sleep_time = 0.05

        # Jump and drift related
        jump_chance = 0.005*sleep_time
        drift_gate = 0*sleep_time
        drift_yellow = 0*sleep_time

        # init of vars
        resonance_gate = resonance_start_gate;
        resonance_newfocus = resonance_start_newfocus;
        resonance_yellow = resonance_start_yellow;
        p_charge_repump = 1

    def signal(self):
        # system simulated as we do a CR check
        repump_duration = 300.0 * 10**(-6)        # s, 300.0 * 10**(-6) SSRO value
        CR_duration = 50.0 * 10**(-6)             # s, 50.0 * 10**(-6) SSRO value
        experiment_duration = 0.001               # s, 
        count_duration = 0.1                      # s, 0.1 SSRO value
        repump_threshold = 5
        experiment_threshold = 20
        ionisation_constant = 0.03
        mode = 1
        do_CR_mod = 0
        state = 1
        peak_counts_cr_check = 35.0
        peak counts_repump = 25

        while True:
            if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): 
            print 'Quit by user'
            return False; 

            self.NV_control()

            counts_per_repump_array = []
            counts_per_cr_check_array = []
            time_start = 0
            t = 0

            while (t - time_start < count_duration):
                if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): 
                print 'Quit by user'
                return False;
                if mode == 1: # Do charge repump

                    self.NV_control()
                    if (state == 0):      
                        counts_per_repump_array.append( np.random.poisson(self.lorentzian_yellow_normalized*gate_counts_repump) + 0.01*random.random() )
                    else:
                        counts_per_repump_array.append( 0.01*random.random() )

                    if (self.p_charge_repump > random.random()):
                        state = -1

                    t += repump_duration
                    if do_CR_mod=1:
                        mode = 2
                    else:
                        mode = 3

                if mode == 2: # Do CR mod

                if mode == 3: # Do CR check

                    self.NV_control()
                    if (state == -1):
                        counts_per_cr_check_array.append( np.random.poisson(peak_counts_cr_check*self.lorentzian_gate_normalized*self.p_spin_repump) + 0.01*random.random() )
                    else:
                        counts_per_cr_check_array.append( 0.01*random.random() )

                    if (self.p_ionisation < random.random()):
                        state = 0

                    t += CR_duration
                    if (counts_per_cr_check_array[-1] < repump_threshold):
                        mode = 1
                    elif (counts_per_cr_check_array[-1] < experiment_threshold):
                        mode = 3
                    else:
                        t += experiment_duration
                        mode = 3


            counts_per_repump = np.sum(counts_per_repump_array)
            counts_per_cr_check = np.sum(counts_per_cr_check_array)

    def NV_control(self):

        get yellow
        get gate 
        get newfocus 

        self.gaussian_newfocus = (1/(peak_width_newfocus*math.sqrt(2*np.pi))) * math.exp( - (newfocus-resonance_newfocus)**2 /(2.0*(peak_width_newfocus**2)) )
        self.gaussian_newfocus_max = (1/(peak_width_newfocus*math.sqrt(2*np.pi)))
        self.gaussian_newfocus_normalized = self.gaussian_newfocus/self.gaussian_newfocus_max
        if (self.gaussian_newfocus_normalized > 0.5):
            self.p_spin_repump = 1
        else:
            self.p_spin_repump = 2*self.gaussian_newfocus_normalized

        self.lorentzian_gate = (1/np.pi)*( (peak_width_gate/2.0) / ((gate-resonance_gate)**2 + (peak_width_gate/2.0)**2) )
        self.lorentzian_gate_max = (1/np.pi)*( 1 / (peak_width_gate/2.0) )
        self.lorentzian_gate_normalized = self.lorentzian_gate/self.lorentzian_gate_max

        self.p_ionisation = ionisation_constant*lorentzian_gate_normalized*p_spin_repump
   

        self.lorentzian_yellow = (1/np.pi)*( (peak_width_yellow/2.0) / ((yellow-resonance_yellow)**2 + (peak_width_yellow/2.0)**2) )
        self.lorentzian_yellow_max = (1/np.pi)*( 1 / (peak_width_yellow/2.0) )
        self.lorentzian_yellow_normalized = self.lorentzian_yellow/self.lorentzian_yellow_max
        if (self.lorentzian_yellow_normalized > 0.5):
            self.p_charge_repump = 1
        else:
            self.p_charge_repump = 2*self.lorentzian_yellow_normalized

        #  # jumps 
        # p_jump = random.random()
        # if p_jump<jump_chance:
        #     resonance_gate += 200.0*random.uniform(-1,1) 
        #     p_yellow = random.random()
        #     if p_yellow < 0.5:
        #         resonance_yellow += 0.2 * (10**6) *random.uniform(-1,1)  
        #     else:
        #         if random.random() > 0.5: 
        #             p_yellow_resonance = random.uniform(0.8,1) 
        #         else:
        #             p_yellow_resonance = random.uniform(-0.8,-1)
        #         resonance_yellow += 0.5 * (10**6) *p_yellow_resonance

        # # and drifts
        # resonance_gate += drift_gate
        # resonance_yellow += drift_yellow




