import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, cm
from analysis.lib.fitting import fit
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
reload(magnetometry)
load_data=False

def analyze_saved_simulations (timestamp,error_bars=False):
    mgnt_exp = magnetometry.AdaptiveMagnetometry(N=14, tau0=20e-9)
    mgnt_exp.error_bars=error_bars
    mgnt_exp.load_analysis (timestamp=timestamp)
    mgnt_exp.calculate_scaling()
    return mgnt_exp

def add_scaling_plot(timestamp, ax, exp_data, label, marker_settings, color, return_all = False):
    #adds a scaling plot to axis 'ax', loading from analyzed data with a given 'timestamp'
    #exp_data=boolena, if 'True' then data is plotted with markers and errorbars are calculated, 
    #otherwise it is considered a simulation, and plotted like a line
    #label, string for legend
    data_file = analyze_saved_simulations (timestamp=timestamp, error_bars=exp_data)

    ax.plot (data_file.total_time, data_file.sensitivity, marker_settings,color=color, label=label)
    if exp_data: 
        ax.fill_between (data_file.total_time, data_file.sensitivity-data_file.err_sensitivity, data_file.sensitivity+data_file.err_sensitivity, color=color, alpha=0.2)
    if return_all:
        return ax, data_file.total_time, data_file.sensitivity, data_file.err_sensitivity, data_file.G, data_file.F, data_file.fid0
    else:
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


def generate_data_dict (timestamps, save_name = None):
    data_dict = {}
    for i,k in enumerate(timestamps):
        data_file  = analyze_saved_simulations (timestamp=k,error_bars=True)
        data_dict [k] = {'total_time': data_file.total_time, 'sensitivity':data_file.sensitivity, 'err_sensitivity':data_file.err_sensitivity, 'F':data_file.F, 'G':data_file.G, 'fid0':data_file.fid0}
    
    return data_dict


def compare_best_sensitivities (data_dict_array, legend_array, title, colours=None):
    #compare_multiple_plots must be run before, to have a data_dict

    f = plt.figure(figsize=(8,6))
    ax = f.add_subplot(1,1,1)
    markers = ['o', '^', 's','v', '>', '<',  '*', 'D', '+','|']
    ccc = np.linspace(0,1,len(data_dict_array))

    idx = 0
    plt.figure()    
    for data_dict in data_dict_array:
        list_F = []
        list_sens = []
        list_err_sens = []

        for i,k in enumerate(data_dict):
            list_F.append(data_dict[k]['F'])
            sensitivity = data_dict[k]['sensitivity']
            #print sensitivity
            min_sens = min(sensitivity)
            list_sens.append(min_sens)
            min_idxs = [iii for iii, val in enumerate(sensitivity) if val == min_sens]  
            m = min_idxs[0]
            list_err_sens.append(data_dict[k]['err_sensitivity'][m])
        
        F = np.array(list_F)
        ind = np.argsort(F)
        F=F[ind]
        sens = np.array(list_sens)[ind]
        err_sens = np.array(list_err_sens)[ind]

        if colours==None:
            c = cm.Set1(ccc[idx])
        else:
            c = colours[idx]

        ax.plot (F, sens, color=c, label = legend_array[idx])
        ax.errorbar (F, sens,  yerr = err_sens, fmt=markers[idx], color=c)
        idx = idx + 1

    ax.set_yscale ('log')
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(15)

    #plt.fill_between (F, sens-err_sens, sens+err_sens, color='RoyalBlue', alpha=0.2)
    ax.set_xlabel ('F', fontsize=15)
    ax.set_ylabel ('best sensitivity [VH*T]', fontsize=15)
    ax.set_title (title, fontsize=15)
    ax.legend()
    plt.show()



def compare_scaling_fits (data_dict_array, legend_array, title='', colours=None, first_N = 3, last_N=10):
    #compare_multiple_plots must be run before, to have a data_dict

    f = plt.figure(figsize=(8,6))
    ax = f.add_subplot(1,1,1)
    markers = ['o', '^', 's','v', '>', '<',  '*', 'D', '+','|']
    ccc = np.linspace(0,1,len(data_dict_array))

    idx = 0
    plt.figure()    
    for data_dict in data_dict_array:
        list_F = []
        list_scaling = []
        list_err_scaling = []

        for i,k in enumerate(data_dict):
            list_F.append(data_dict[k]['F'])
            x0 = np.log10(data_dict[k]['total_time'][first_N:last_N]/20e-9)
            y0 = np.log10(data_dict[k]['sensitivity'][first_N:last_N]/20e-9)
            err_array = np.log10(data_dict[k]['err_sensitivity'][first_N:last_N]/20e-9)

            guess_b = -1
            guess_a = 0*y0[-1]+1
            
            a = fit.Parameter(guess_a, 'a')
            b = fit.Parameter(guess_b, 'b')
            
            p0 = [a, b]
            fitfunc_str = ''

            def fitfunc(x):
                return a()+b()*x


            fit_result = fit.fit1d(x0,y0, None, p0=p0, fitfunc=fitfunc, fixed=[],
                    do_print=False, ret=True, err_y = err_array)
            b_fit = fit_result['params_dict']['b']
            b_err = fit_result['error_dict']['b']
            list_scaling.append(b_fit)
            list_err_scaling.append(b_fit)

        
        F = np.array(list_F)
        ind = np.argsort(F)
        F=F[ind]
        scaling_coeff = np.array(list_scaling)[ind]
        err_scaling_coeff = np.array(list_err_scaling)[ind]

        if colours==None:
            c = cm.Set1(ccc[idx])
        else:
            c = colours[idx]

        ax.plot (F, scaling_coeff, color=c, label = legend_array[idx])
        ax.plot (F, scaling_coeff, markers[idx], color=c)
        idx = idx + 1

    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(15)

    #plt.fill_between (F, sens-err_sens, sens+err_sens, color='RoyalBlue', alpha=0.2)
    ax.set_ylim([-1.1, 1.1])
    ax.set_xlabel ('F', fontsize=15)
    ax.set_ylabel ('scaling coefficient', fontsize=15)
    ax.set_title (title, fontsize=15)
    ax.legend()
    plt.show()
