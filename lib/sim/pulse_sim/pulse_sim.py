"""
Module to parse pulse sequences that have been dumped on the way to the AWG, and extract the useful information
Peter H, 11/2016
"""

### imports
import sys, os
import numpy as np
import h5py
from matplotlib import pyplot as plt
from tabulate import tabulate
import pickle as pkl
import h5py

def get_AWG_seq_file(basedir,name):
    name = os.path.splitext(name)[0]

    with open(os.path.join(basedir,name+'.pickle'), 'rb') as f:  # Python 3: open(..., 'rb')
        combined_seq,combined_list_of_elements = pkl.load(f)
        f.close()
    return combined_seq,combined_list_of_elements

def get_pulse_elem_for_seq_element(seq_elem,combined_list_of_elements):
    return next((x for x in combined_list_of_elements if x.name == seq_elem['wfname']), None) # Pull out appropriate pulse sequence

def get_pulse_elem_by_name(name,combined_list_of_elements):
    return next((x for x in combined_list_of_elements if x.name == name), None) # Pull out appropriate pulse sequence

def get_seq_element_pulses(elem, elem_start = 0, verbose = False):
    pulses = np.array([]).reshape(0,5)
    for key,val in elem.pulses.iteritems():
        if verbose:
            print key
            print val.name
            
        has_amplitude = (val.amplitude > 0) if hasattr(val,'amplitude') else False
        has_amplitude = (has_amplitude and val.env_amplitude > 0) if hasattr(val,'env_amplitude') else has_amplitude
        
        if ('MW_Imod' in val.channels or 'MW_Qmod' in val.channels) and has_amplitude:
            if val.pi2_pulse == True:
                rotation = np.pi/2
            else:
                rotation = np.pi
            special = 0
            pulse_mid_point = elem_start + (val.effective_start()+val.effective_stop())/2 
            pulse_length = (val.effective_stop()-val.effective_start())
            entry = np.array((pulse_mid_point,pulse_length,\
                    np.pi * val.phase/180, rotation,special)).reshape(-1,1).T
            pulses = np.vstack([pulses,entry])
        
        elif 'AOM_Newfocus' in val.channels and has_amplitude:
            rotation = 0
            phase = 0
            special = 1
            pulse_mid_point = elem_start + (val.effective_start()+val.effective_stop())/2 
            pulse_length = (val.effective_stop()-val.effective_start())
            entry = np.array((pulse_mid_point,pulse_length,\
                    phase, rotation,special)).reshape(-1,1).T
            pulses = np.vstack([pulses,entry])
       
        elif 'adwin_sync' in val.channels and has_amplitude:
            rotation = 0
            phase = 0
            special = 2
            pulse_mid_point = elem_start + (val.effective_start()+val.effective_stop())/2 
            pulse_length = (val.effective_stop()-val.effective_start())
            entry = np.array((pulse_mid_point,pulse_length,\
                    phase, rotation,special)).reshape(-1,1).T
            pulses = np.vstack([pulses,entry])
            
    pulses = pulses[pulses[:,0].argsort()]   # Sort by pulse start
    return pulses

def group_seq_elems(combined_seq,combined_list_of_elements):
    jumps_or_go_tos = []
    sequence_overview = []
    pulse_array = []
    current_group_pulses = np.array([]).reshape(0,5)
    time = 0
    first_group_elem_name = None 
    trigger_wait = False
    
    for seq_elem in combined_seq.elements:
        
        if seq_elem['goto_target'] != None: # End of group!
                jumps_or_go_tos.append(seq_elem['goto_target'])
            
        if seq_elem['jump_target'] != None: # End of group!
                jumps_or_go_tos.append(seq_elem['jump_target'])
    
    for seq_elem in combined_seq.elements:
        
        elem = get_pulse_elem_for_seq_element(seq_elem,combined_list_of_elements)

        current_elem_pulses = get_seq_element_pulses(elem,verbose = False)
        reps = int(seq_elem['repetitions'])
        
        if seq_elem['goto_target'] != None or seq_elem['jump_target'] != None:
            reps = 1 # Dont repeat an element if you can jump out of it.
            
        pulses = np.array([]).reshape(0,5)
            
        for x in range(reps):
            if np.size(current_elem_pulses):
                temp_pulses = np.copy(current_elem_pulses)
                temp_pulses[:,0] += time
                pulses = np.vstack([pulses,temp_pulses])
            time += elem.length()
            
        if first_group_elem_name == None:
            first_group_elem_name = seq_elem['name']
        
        start_of_group = False
        end_of_group = False
        
        if seq_elem['trigger_wait'] != False or seq_elem['name'] in jumps_or_go_tos: # Start of new group!
            
            if np.size(current_group_pulses):
                sequence_overview.append({'name' : first_group_elem_name,\
                               'trigger_wait' : int(trigger_wait),\
                               'goto_target' : str(seq_elem['goto_target']),\
                               'jump_target' : str(seq_elem['jump_target']), 'final_time' : float(time)}) 
                pulse_array.append(current_group_pulses)
            first_group_elem_name = seq_elem['name']
            trigger_wait = seq_elem['trigger_wait']
            current_group_pulses = pulses
            
            start_of_group = True
         
        if seq_elem['goto_target'] != None or seq_elem['jump_target'] != None:
            end_of_group = True
            
            if not(start_of_group): current_group_pulses = np.vstack([current_group_pulses,pulses])
            sequence_overview.append({'name' : first_group_elem_name,\
                               'trigger_wait' : int(trigger_wait),\
                               'goto_target' : str(seq_elem['goto_target']),\
                               'jump_target' : str(seq_elem['jump_target']), 'final_time' : float(time)}) 
            pulse_array.append(current_group_pulses) 
            time = 0
            first_group_elem_name = None
            trigger_wait = False
            current_group_pulses = np.array([]).reshape(0,5)
        elif end_of_group == False and start_of_group == False:
            current_group_pulses = np.vstack([current_group_pulses,pulses])
            
    
    if np.size(current_group_pulses):
                sequence_overview.append({'name' : first_group_elem_name,\
                               'trigger_wait' : int(trigger_wait),\
                               'goto_target' : str(seq_elem['goto_target']),\
                               'jump_target' : str(seq_elem['jump_target']), 'final_time' : float(time)}) 
                pulse_array.append(current_group_pulses)

    return [sequence_overview,pulse_array]

def save_grouped_pulses(basedir,name,grouped_seq):
    
    name = os.path.splitext(name)[0]
    
    with h5py.File(os.path.join(basedir,name+'.h5'), 'w') as hf:
        
        for key in grouped_seq[0][0].keys():
            key_data = np.array([ item[key] for item in grouped_seq[0]])
            hf.create_dataset('seq_overview/'+key, data= key_data)
        
        for i,pulses in enumerate(grouped_seq[1]):
            hf.create_dataset('pulses/'+str(i), data=pulses)

def load_grouped_pulses(basedir,name):
    
    name = os.path.splitext(name)[0]
    
    with h5py.File(os.path.join(basedir,name+'.h5'), 'r') as hf:
        
        num_pulses = np.max(np.array(hf['pulses'].keys()).astype(int))
        overview_keys = [str(j) for j in (hf['seq_overview'].keys())]
        overview_data = {}
        for key in overview_keys:
            overview_data[key] = np.array(hf.get('seq_overview/'+key))
        
        grouped_seq = []
        for i in range(num_pulses):
            overview = {}
            for key in overview_keys:
                overview[key] = overview_data[key][i]
            pulses = np.array(hf.get('pulses/'+str(i)))
            grouped_seq.append([overview,pulses])
    
    return grouped_seq

def batch_process_folder(basedir, verbose = False):
    for file in os.listdir(basedir):
        if file.endswith(".pickle"):
            if verbose:
                print file
            combined_seq,combined_list_of_elements = get_AWG_seq_file(basedir,file)
            grouped_seq = group_seq_elems(combined_seq,combined_list_of_elements)
            save_grouped_pulses(basedir,file,grouped_seq)

def print_group_summaries(grouped_seq):
    print 'Extracted pulse groups'
    tab_data = []
    for i, group in enumerate(grouped_seq[0]):
        tab_data.append([i, group['name'][:30], group['trigger_wait'],\
                         str(group['goto_target'])[:30], str(group['jump_target'])[:30]])
    print tabulate(tab_data, headers=['#', 'First elem', 'Wait', 'Go To', 'Jump To'])

        
def draw_seq_element_pulses(pulses, start = 0, elem_length = 0, time_res = 1e-9):
    
    total_duration =  elem_length if elem_length > 0 else (pulses[-1][0]+pulses[-1][1]/2 - start)
    plot_duration = 1.1*total_duration
    start_steps = int(start/time_res)
    time_steps = int(plot_duration/time_res)
    time = np.arange(time_steps)*time_res*1e6
    I_channel = np.zeros(time_steps)
    Q_channel = np.zeros(time_steps)
    repump_channel = np.zeros(time_steps)
    for pulse in pulses:
        pulse_start = int((pulse[0]-pulse[1]/2)/time_res) - start_steps
        pulse_stop = int((pulse[0]+pulse[1]/2)/time_res) - start_steps
        
        if pulse[4] == 0:
            I_channel[pulse_start:pulse_stop] = np.sin(pulse[2]) * pulse[3]/np.pi
            Q_channel[pulse_start:pulse_stop] = np.cos(pulse[2]) * pulse[3]/np.pi
        else:
            repump_channel[pulse_start:pulse_stop] = 1
            
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.plot(time,I_channel, label='I mod')
    plt.plot(time,Q_channel, label='Q mod')
    plt.plot(time,repump_channel, label='repump')
    plt.plot([total_duration*1e6,total_duration*1e6],[-1,1],'--k')

    plt.ylim([-1.1,1.1])
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Pulse amplitude')
    plt.show()
    plt.close()