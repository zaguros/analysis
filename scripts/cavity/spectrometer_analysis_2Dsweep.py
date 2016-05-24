# analyse 2D data of spectrometer 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import matplotlib.image as mpimg
import scipy 
import time

from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
import analysis.scripts.cavity.spectrometer_analysis as sa
reload(sa)



# parameters to vary per measurement Note: you might have to change the vmin and vmax of the colorbar inside the script! 
V_min = -2
# V_min = 4
V_max = 10
n_diamond = 2.4 #refractive index diamond
c = 3.e8 #speed of light


def get_data(data_dir):
    wavelengths,filenumbers,intensities = sa.load_data_from_folder(data_dir)
    print data_dir
    return wavelengths,filenumbers,intensities 

def plot_data(data_dir,wavelengths,intensities,vmax = None):

    fig,ax = plt.subplots()
    ax=sns.heatmap(intensities, vmax = vmax,cmap='YlGnBu',ax=ax)

    ax = set_axes(ax,wavelengths)

    try: 
        print os.path.join(data_dir, '2D_plot.jpg')
        fig.savefig(os.path.join(data_dir, '2D_plot.jpg'))
    except:
        print('could not save figure')


def set_axes(ax,wavelengths):
    ax.set_xlabel("Voltage (V)", fontsize = 14)
    ax.set_ylabel("Frequency (THz)", fontsize = 14)

    ax.grid(False)
    ax.set_axis_bgcolor('white')


    ax.tick_params(which = 'both', direction = 'out')
    xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],int((V_max-V_min)/2+1))
    xticklabels = np.linspace(V_min,V_max,int((V_max-V_min)/2+1))
    xticklabels = np.round(xticklabels,1)

    yticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],7)
    ytickslabels = np.linspace(wavelengths[-1],wavelengths[0],7)#in THz
    ytickslabels =np.round(ytickslabels,0).astype(int)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=0)

    ax.set_yticks(yticks)
    ax.set_yticklabels(ytickslabels)

    return ax


def plot_from_2D_data(data_dir,vmax = None):
    print data_dir
    wavelengths,filenumbers,intensities = get_data(data_dir)
    plot_data(data_dir,wavelengths,intensities,vmax = vmax)    

def peaks_from_2D_data(data_dir,order_peak_detection=200):
    wavelengths,filenumbers,intensities = get_data(data_dir)

    x=np.array([])
    u_x = np.array([])
    y=np.array([])

    print 'getting peak locations'
    for i,intensity in enumerate(np.transpose(intensities)):
        print i,'/',len(np.transpose(intensities))
        indices,peak_wavelengths,peak_intensity = sa.approximate_peak_location(wavelengths,intensity,order_peak_detection=order_peak_detection)
        x0s,u_x0s = sa.fit_peak(wavelengths,intensity,indices,peak_wavelengths,peak_intensity)
        x = np.append(x,x0s)
        u_x = np.append(u_x,u_x0s)
        for j in np.arange(len(x0s)):
            y = np.append(y,i)

    # make a scatter plot of the peaks
    fig,ax = plt.subplots(figsize =(6,4))
    ax.scatter(y,x)
    ax.errorbar(y, x, ls='none',marker=None,yerr= u_x)

    ax.set_xlim((filenumbers[0],filenumbers[-1]))
    ax.set_ylim((wavelengths[-1], wavelengths[0]))

    ax=set_axes(ax,wavelengths)

    return fig,ax

def overlap_peaks_and_modes(diamond_thickness=4.e-6,cavity_length = 1.e-6,conversion_factor = -150.e-9,order_peak_detection=70):
    fig,ax = peaks_from_2D_data(order_peak_detection=order_peak_detection)
    ax = plot_diamond_modes(diamond_thickness=diamond_thickness,ax = ax)
    ax = plot_air_modes(cavity_length=cavity_length,ax=ax, conversion_factor=conversion_factor)

    title ='d={}um_L={}um_cf={}nmpV'.format(str(diamond_thickness*1e6),str(cavity_length*1.e6),str(conversion_factor*1.e9))

    ax.text(ax.get_xlim()[0] + (ax.get_xlim()[-1]-ax.get_xlim()[0])/4,ax.get_ylim()[0],title)

    #add an axis at the top with the cavity length 
    ax2 = ax.twiny()
    xticks = np.linspace(ax2.get_xlim()[0],ax2.get_xlim()[-1],int((V_max-V_min)/2+1))
    xticklabels2 =np.linspace(cavity_length*1.e6,cavity_length*1.e6+(conversion_factor*(V_max-V_min)*1.e6),int((V_max-V_min)/2+1))
    xticklabels2 = np.round(xticklabels2,2)
    print xticklabels2
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticklabels2,rotation=0)
    ax2.set_xlabel('cavity length (um)',fontsize = 14)

    try: 
        print os.path.join(data_dir, 'overlap_peaks_and_modes{}.png'.format(title))
        fig.savefig(os.path.join(data_dir, 'overlap_peaks_and_modes{}.png'.format(title)))
    except:
        print('could not save figure')


    plt.close(fig)

def pure_diamond_modes(diamond_thickness=4.e-6):
    max_nr_modes = 40
    nu_diamond = np.zeros(max_nr_modes)
    for N in np.arange(max_nr_modes):
        nu_diamond[N] = (N * c / (2 * n_diamond *diamond_thickness))/1.e12 # in THz

    return nu_diamond

def plot_diamond_modes(diamond_thickness=4.e-6,ax = None):
    return_fig = False
    if ax == None:
        return_fig = True
        fig,ax = plt.subplots()

    nu_diamond = pure_diamond_modes(diamond_thickness)

    for N,nu in enumerate(nu_diamond):
        ax.plot([-1000,1000],[nu,nu], lw=2)
        ax.text(ax.get_xlim()[-1],nu, 'N={}'.format(N))

    if return_fig:
        return fig,ax

    return ax

def pure_air_modes(cavity_length=1.e-6,conversion_factor = -150.e-9):
    delta_V = V_max - V_min
    delta_L = delta_V*(conversion_factor) # in m
    print 'delta L',delta_L

    max_nr_modes = 20
    nu_air = np.zeros((max_nr_modes,3))
    for N in np.arange(max_nr_modes):
        nu_air[N,0] = (N * c / (2 * cavity_length))/1.e12 # in THz
        nu_air[N,1] = (N * c / (2 * (cavity_length+delta_L)))/1.e12 # in THz

    return nu_air

def plot_air_modes(cavity_length=1.e-6,ax = None,conversion_factor = -150.e-9):
    return_fig = False
    if ax == None:
        return_fig = True
        fig,ax = plt.subplots()

    nu_air = pure_air_modes(cavity_length=cavity_length,conversion_factor=conversion_factor)

    ax.get_ylim()
    ax.get_xlim()

    for N,nu in enumerate(nu_air):
        ax.plot(ax.get_xlim(),[nu[0],nu[1]], lw=2)
        ax.text(ax.get_xlim()[0],nu[0], 'N={}'.format(N))

    if return_fig:
        return fig,ax

    return ax

