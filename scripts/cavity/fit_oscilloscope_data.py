### fit oscilloscope data


import numpy as np 
from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
import pandas as pd


def load_data(filepath):
    """
    Loading csv -files into pandas dataframe and convert to array
    Input: 
    filepath  - the filepath of csv file
    Output: 
    data  -  numpy ndarray with (time, intensity)
    """
    data = pd.read_csv(filepath, skiprows=16, names = ["X","Y","none"],usecols=[0,1]) #creating a dataframe in pandas and importing the data
    data = data.as_matrix()
    x = data[:,0]*1.e5 # x-axis does not matter, factor for clarity. 
    y = data[:,1]
    return data, x, y



def get_linewidth(data,EOM_freq, g_a1 = 0.5, g_A1 = 0.04 , g_x01 = 0.0, g_gamma1 = 0.1, g_dx = 0.35, g_A2 = 0.02, g_gamma2 = 0.1, plot_data = True):
    """
    fit 3 lorentzians, get linewidth and the error
    input:
    data  = numpy ndarray with the data
    plot_data = False/True, default True
    g_a1 = offset guess parameter for fit
    g_A1 = guess area of peak 1 for fit
    g_x01 = guess x value of the fit 
    g_gamma1 = guess linewidth of middle peak
    g_dx = guess separation between peaks 
    g_A2 = guess area of peak 2 and 3 for fit
    g_gamma2 = guess linewidth of peak 2 and 3 for fit
    output:
    lw = linewidth of middle peak
    u_lw = error of the linewidth
    """

    n_xticks = 10
    fixed=[]
    
    p0, fitfunc, fitfunc_str = common.fit_3lorentz_symmetric(g_a1, g_A1, g_x01, g_gamma1, g_dx, g_A2, g_gamma2)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)

    x01 = fit_result['params_dict']['x01']
    dx = fit_result['params_dict']['dx']
    gamma1 = fit_result['params_dict']['gamma1']
    gamma2 = fit_result['params_dict']['gamma2']
    u_gamma1 = fit_result['error_dict']['gamma1']
    scaling = EOM_freq/dx #scale the known EOM freq with the separation here.
    lw = gamma1*scaling #scale the linewidth to get linewidth in frequency
    u_lw = u_gamma1*scaling
    return lw, u_lw, scaling, fit_result

def plot_data_and_fit(data, EOM_freq):

    fig,ax = plt.subplots(figsize=(8,4))
    lw, u_lw, scaling, fit_result = get_linewidth(data, EOM_freq)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],len(x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)
    ax = set_axes(ax)
    plt.show()

    # ax.set_xlabel("Frequency (GHz)", fontsize = 14)
    # ax.set_ylabel("Intensity (a.u.)", fontsize = 14)
    # ax.set_xlim(x[0],x[-1])
    # #rescaling x-axis in GHz
    # X_min_freq = ax.get_xlim()[0]*scaling
    # X_max_freq = ax.get_xlim()[-1]*scaling
    # xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],n_xticks) 
    # xticklabels = np.linspace(X_min_freq,X_max_freq,n_xticks)
    # xticklabels_round=[]
    # for j in xticklabels:
    #   round_ = round(j,0)
    #   xticklabels_round = np.append(xticklabels_round,round_)
    # ax.set_xticklabels(xticklabels_round)
    # ax.set_xticks(xticks)

def set_axes(ax):
    ax.set_xlabel("Frequency (GHz)", fontsize = 14)
    ax.set_ylabel("Intensity (a.u.)", fontsize = 14)

    ax.grid(False)
    ax.set_axis_bgcolor('white')
    #ax.set_xlim(x[0],x[-1])

    ax.tick_params(which = 'both', direction = 'out')
    xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],int((ax.get_xlim()[0]-ax.get_xlim()[-1])/2+1))
    xticklabels = np.linspace(ax.get_xlim()[0]*scaling,ax.get_xlim()[-1]*scaling,int((ax.get_xlim()[0]-ax.get_xlim()[-1])/2+1))
    xticklabels = np.round(xticklabels,1)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=0)

    return ax