import numpy as np
import math
from matplotlib import pyplot as plt
import random
import time
import types
import os
import msvcrt
from scipy.stats import cauchy

class NV_energy_simulator(object):

    def __init__(self):
        # Gaussian related
        self.peak_width_gate = 0.05*10**6;
        self.peak_width_newfocus = 0.3*10**6; #GHz
        self.peak_width_yellow = 0.05*10**6;
        self.peak_counts_cr_check = 0.13;
        self.peak_counts_repump = 0.75;
        self.dark_counts = 0.0002
        resonance_start_gate = 57.0*10**6;
        resonance_start_newfocus = 45.0*10**6; #GHz
        resonance_start_yellow = 17.0*10**6;
        self.ionisation_constant = 0.3

        # Time related
        sleep_time = 1

        # Jump and drift related
        self.jump_chance = 0.005*sleep_time
        self.drift_gate = 0*sleep_time
        self.drift_yellow = 0*sleep_time

        # init of vars
        self.resonance_gate = resonance_start_gate;
        self.resonance_newfocus = resonance_start_newfocus;
        self.resonance_yellow = resonance_start_yellow;
        self.p_charge_repump = 1

    def signal(self):
        # system simulated as we do a CR check
        repump_duration = 0.0003                    # s, 300.0 * 10**(-6) SSRO value
        CR_duration = 0.00005                       # s, 50.0 * 10**(-6) SSRO value
        experiment_duration = 0.001                 # s, 
        count_duration = 0.1                        # s, 0.1 SSRO value
        repump_threshold = 5
        experiment_threshold = 20
        mode = 1
        do_CR_mod = 0
        state = 1
        step_increment = 0
        time_array = []
        counts_per_repump_graph = []
        counts_per_cr_check_graph = []

        plt.ion()
        plt.show()

        while True:
            if (msvcrt.kbhit() and (msvcrt.getch() == 'q')) or (step_increment>=1000): 
                print 'Quit by user'
                return False; 

            self.NV_control()

            counts_per_repump_array = []
            counts_per_cr_check_array = []
            time_start = 0
            t = 0

            while (t - time_start < count_duration):
                if (msvcrt.kbhit() and (msvcrt.getch() == 'q')): 
                    print 'stopped'
                    return False;
                if mode == 1: # Do charge repump

                    self.NV_control()
                    if (state == 0):      
                        counts_per_repump_array.append( np.random.exponential(self.lorentzian_yellow_normalized*self.peak_counts_repump) + self.dark_counts*random.random() )
                    else:
                        counts_per_repump_array.append( self.dark_counts*random.random() )

                    if (self.p_charge_repump > random.random()):
                        state = -1

                    t += repump_duration
                    if (do_CR_mod == 1):
                        mode = 2
                    else:
                        mode = 3

                if mode == 2: # Do CR mod
                    print 'CR mod, apparently....'
                if mode == 3: # Do CR check

                    self.NV_control()
                    if (state == -1):
                        counts_per_cr_check_array.append( np.random.exponential(self.peak_counts_cr_check*self.lorentzian_gate_normalized*self.p_spin_repump) + self.dark_counts*random.random() )
                    else:
                        counts_per_cr_check_array.append( self.dark_counts*random.random() )

                    if (self.p_ionisation > random.random()):
                        state = 0

                    t += CR_duration
                    if (counts_per_cr_check_array[-1] < repump_threshold):
                        mode = 1
                    elif (counts_per_cr_check_array[-1] < experiment_threshold):
                        mode = 3
                    else:
                        t += experiment_duration
                        mode = 3
            step_increment += 1.0
            time_array.append(0.1*step_increment)
            counts_per_repump_graph.append(np.sum(counts_per_repump_array))
            counts_per_cr_check_graph.append(np.sum(counts_per_cr_check_array))

            plt.subplot(211)
            plt.plot(time_array,counts_per_cr_check_graph)
            plt.subplot(212)
            plt.plot(time_array,counts_per_repump_graph, 'r')
            # plt.show()


    def NV_control(self):
        gate = (57.0-0.0)*10**6;
        newfocus = (45.0-0.0)*10**6; #GHz
        yellow = (17.0-0.0)*10**6;   #GHz

        self.gaussian_newfocus = (1/(self.peak_width_newfocus*math.sqrt(2.0*np.pi))) * math.exp( - (newfocus-self.resonance_newfocus)**2 /(2.0*(self.peak_width_newfocus**2)) )
        self.gaussian_newfocus_max = (1/(self.peak_width_newfocus*math.sqrt(2.0*np.pi)))
        self.gaussian_newfocus_normalized = self.gaussian_newfocus/self.gaussian_newfocus_max
        if (self.gaussian_newfocus_normalized > 0.5):
            self.p_spin_repump = 1
        else:
            self.p_spin_repump = 2*self.gaussian_newfocus_normalized

        self.lorentzian_gate = (1/np.pi)*( (self.peak_width_gate/2.0) / ((gate-self.resonance_gate)**2 + (self.peak_width_gate/2.0)**2) )
        self.lorentzian_gate_max = (1/np.pi)*( 1 / (self.peak_width_gate/2.0) )
        self.lorentzian_gate_normalized = self.lorentzian_gate/self.lorentzian_gate_max

        self.p_ionisation = self.ionisation_constant*self.lorentzian_gate_normalized*self.p_spin_repump


        self.lorentzian_yellow = (1/np.pi)*( (self.peak_width_yellow/2.0) / ((yellow-self.resonance_yellow)**2 + (self.peak_width_yellow/2.0)**2) )
        self.lorentzian_yellow_max = (1/np.pi)*( 1 / (self.peak_width_yellow/2.0) )
        self.lorentzian_yellow_normalized = self.lorentzian_yellow/self.lorentzian_yellow_max
        if (self.lorentzian_yellow_normalized > 1):
            self.p_charge_repump = 1
        else:
            self.p_charge_repump = self.lorentzian_yellow_normalized

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

        # and drifts
        self.resonance_gate += self.drift_gate
        self.resonance_yellow += self.drift_yellow














