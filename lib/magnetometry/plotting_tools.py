import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, cm
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
reload(magnetometry)
load_data=False

def analyze_saved_simulations (timestamp,error_bars=False):
    mgnt_exp = magnetometry.AdaptiveMagnetometry(N=14, tau0=20e-9)
    mgnt_exp.error_bars=error_bars
    mgnt_exp.load_analysis (timestamp=timestamp)
    mgnt_exp.calculate_scaling()
    return mgnt_exp

def add_scaling_plot(timestamp, ax, exp_data, label, marker_settings, color):
    #adds a scaling plot to axis 'ax', loading from analyzed data with a given 'timestamp'
    #exp_data=boolena, if 'True' then data is plotted with markers and errorbars are calculated, 
    #otherwise it is considered a simulation, and plotted like a line
    #label, string for legend
    data_file = analyze_saved_simulations (timestamp=timestamp, error_bars=exp_data)

    ax.plot (data_file.total_time, data_file.sensitivity, marker_settings,color=color, label=label)
    if exp_data: 
        ax.fill_between (data_file.total_time, data_file.sensitivity-data_file.err_sensitivity, data_file.sensitivity+data_file.err_sensitivity, color=color, alpha=0.2)
    return ax

def compare_2plots(timestamp1, timestamp2, title):
    f = plt.figure(figsize=(8,6))
    ax = f.add_subplot(1,1,1)
    ax = add_scaling_plot (timestamp = timestamp1, exp_data=True, ax=ax, label = 'adapt', marker_settings='o', color='b')
    ax = add_scaling_plot (timestamp = timestamp2, exp_data=True, ax=ax, label = 'non adapt', marker_settings='^', color='r')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel ('total ramsey time [$\mu$s]', fontsize=15)
    ax.set_ylabel ('$V_{H}T$', fontsize=15)
    plt.title (title, fontsize=20)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(15)
    plt.legend()
    plt.show()

def compare_multiple_plots(timestamps, labels, title, colours = None):
    #MAX 10 plots!!!! Then no more markers!
    f = plt.figure(figsize=(8,6))
    ax = f.add_subplot(1,1,1)
    markers = ['o', '^', 's','v', '>', '<',  '*', 'D', '+','|']
    ccc = np.linspace(0,1,len(timestamps))
    for i,k in enumerate(timestamps):
        if colours==None:
            c = cm.Set1(ccc[i])
        else:
            c = colours[i]
        ax = add_scaling_plot (timestamp = k, exp_data=True, ax=ax, label = labels[i], marker_settings=markers[i], color=c)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel ('total ramsey time [$\mu$s]', fontsize=15)
    ax.set_ylabel ('$V_{H}T$', fontsize=15)
    plt.title (title, fontsize=20)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(15)
    plt.legend()
    plt.show()