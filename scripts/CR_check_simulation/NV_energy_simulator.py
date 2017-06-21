import numpy as np
import math
import matplotlib
matplotlib.use('GTKAgg')
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from IPython import display
import random
import time
import types
import os
import msvcrt
from scipy.stats import cauchy

class NV_energy_simulator(object):

    def __init__(self):
        # Gaussian related
        self.peak_width_gate = 0.10*10**6;
        self.peak_width_newfocus = 0.25*10**6; #GHz
        self.peak_width_yellow = 0.12*10**6;
        self.peak_counts_cr_check = 0.30;
        self.peak_counts_repump = 0.35;
        self.dark_counts = 0.0005
        resonance_start_gate = 57.0*10**6;
        resonance_start_newfocus = 45.0*10**6; #GHz
        resonance_start_yellow = 17.0*10**6;
        self.ionisation_constant = 0.3

        # Time related
        sleep_time = 0.1

        # Jump and drift related
        self.do_jumps = True
        self.do_drift = False
        self.jump_chance = 0.00015*sleep_time
        self.drift_gate = 0.5*sleep_time
        self.drift_yellow = 0.5*sleep_time

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
        repump_threshold = 0.01
        experiment_threshold = 0.05
        plot_length = 1000.0
        do_plot = 1.0

        # init
        mode = 1
        do_CR_mod = 0
        state = 1
        step_increment = 0
        time_array = np.arange(-100,0,0.1)
        counts_per_cr_check_graph = np.zeros(plot_length)
        counts_per_repump_graph = np.zeros(plot_length)
        resonance_gate_graph = np.zeros(plot_length)
        resonance_yellow_graph = np.zeros(plot_length)
        x = 0
        y = 0

        # create fig
        f, (ax1,ax2,ax3,ax4) = plt.subplots(4, 1, sharex='col', sharey='row', figsize=(12,10))
        plt.ion()
        display.clear_output(wait=True)
        display.display(plt.gcf())

        while True:
            if (msvcrt.kbhit() and (msvcrt.getch() == 'q')) or (step_increment>=1000):
                print 'Quit by user'
                return False; 

            # init before each round
            counts_per_repump_array = []
            counts_per_cr_check_array = []
            time_start = 0
            t = 0
            cur_time = time.time()

            # calculation
            while (t - time_start < count_duration):
                if (msvcrt.kbhit() and (msvcrt.getch() == 'q')):
                    print 'stopped'
                    return False;
                if mode == 1: # Do charge repump

                    self.NV_control()
                    if (state == 0):      
                        counts_per_repump_array.append( np.random.exponential(self.lorentzian_yellow_normalized*self.peak_counts_repump*self.p_charge_repump/(self.p_spin_repump**0.5)) + self.dark_counts*random.random() )
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
                        # maybe want to change this calculation into something different JM
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

            # rotate array
            counts_per_cr_check_graph = np.roll(counts_per_cr_check_graph,-1)
            counts_per_repump_graph = np.roll(counts_per_repump_graph,-1)
            resonance_yellow_graph = np.roll(resonance_yellow_graph,-1)
            resonance_gate_graph = np.roll(resonance_gate_graph,-1)
            time_array = np.roll(time_array, -1)

            # Fill array
            step_increment += 1.0
            time_array[-1] = 0.1*step_increment
            counts_per_cr_check_graph[-1] = (np.sum(counts_per_cr_check_array))
            counts_per_repump_graph[-1] = (np.sum(counts_per_repump_array))
            resonance_gate_graph[-1] = self.resonance_gate
            resonance_yellow_graph[-1] = self.resonance_yellow

            # print data
            if (0.1*step_increment > do_plot):
                ax1.cla()
                ax2.cla()
                ax3.cla()
                ax4.cla()            
                ax1.set_title('Counts per CR check')
                ax1.set_ylim([0,45])
                ax1.plot(time_array,counts_per_cr_check_graph, 'b')
                ax2.set_title('Counts per repump')
                ax2.set_ylim([0,35])
                ax2.plot(time_array,counts_per_repump_graph, 'r')
                ax3.set_title('gate resonance value')
                ax3.set_ylim([55.5*10**6,58.5*10**6])
                ax3.plot(time_array,resonance_gate_graph, 'b')
                ax4.set_title('yellow resonance value')
                ax4.set_ylim([15.5*10**6,18.5*10**6])
                ax4.plot(time_array,resonance_yellow_graph, 'b')
                plt.tight_layout()
                display.clear_output(wait=True)
                display.display(plt.gcf())
                do_plot += 1.0
                # print do_plot
            
            if time.time() < (cur_time+0.1):
                time.sleep((0.1 - (time.time()-cur_time)))
                # qt.msleep((0.1 - (time.time()-cur_time)))



    def NV_control(self):
        # Control values
        gate = (57.0-0.0)*10**6;
        newfocus = (45.0-0.0)*10**6; #GHz
        yellow = (17.0-0.0)*10**6;   #GHz
        # next step -> creating virtual instruments (PID's and auto optimizer) JM
        # PID's; yellowfrq, gate, E_primer
        # babysitter/auto_optimizer

        # Create values
        self.gaussian_newfocus = (1/(self.peak_width_newfocus*math.sqrt(2.0*np.pi))) * math.exp( - (newfocus-self.resonance_newfocus)**2 /(2.0*(self.peak_width_newfocus**2)) )
        self.gaussian_newfocus_max = (1/(self.peak_width_newfocus*math.sqrt(2.0*np.pi)))
        self.gaussian_newfocus_normalized = self.gaussian_newfocus/self.gaussian_newfocus_max

        self.lorentzian_gate = (1/np.pi)*( (self.peak_width_gate/2.0) / ((gate-self.resonance_gate)**2 + (self.peak_width_gate/2.0)**2) )
        self.lorentzian_gate_max = (1/np.pi)*( 1 / (self.peak_width_gate/2.0) )
        self.lorentzian_gate_normalized = self.lorentzian_gate/self.lorentzian_gate_max

        self.lorentzian_yellow = (1/np.pi)*( (self.peak_width_yellow/2.0) / ((yellow-self.resonance_yellow)**2 + (self.peak_width_yellow/2.0)**2) )
        self.lorentzian_yellow_max = (1/np.pi)*( 1 / (self.peak_width_yellow/2.0) )
        self.lorentzian_yellow_normalized = self.lorentzian_yellow/self.lorentzian_yellow_max

        self.p_spin_repump = self.gaussian_newfocus_normalized*self.gaussian_newfocus_normalized
        self.p_ionisation = self.ionisation_constant*self.lorentzian_gate_normalized*self.p_spin_repump
        self.p_charge_repump = 2.0*self.lorentzian_yellow_normalized*self.lorentzian_yellow_normalized

         # Jumps 
        if (self.do_jumps):
            p_jump = random.random()
            if p_jump<self.jump_chance:
                self.resonance_gate += 150000.0*random.uniform(-1,1) 
                if random.random() < 0.3:
                    self.resonance_yellow += 0.15 * (10**6) *random.uniform(-1,1)  
                else:
                    p_chance = random.random()
                    if p_chance > 0.95: 
                        p_yellow_resonance = random.uniform(0.8,1) 
                    elif p_chance < 0.05:
                        p_yellow_resonance = random.uniform(-0.8,-1)
                    else:
                        p_yellow_resonance = 0
                    self.resonance_yellow += 0.45 * (10**6) *p_yellow_resonance

        # Drifts
        if (self.do_drift):
            self.resonance_gate += self.drift_gate
            self.resonance_yellow += self.drift_yellow