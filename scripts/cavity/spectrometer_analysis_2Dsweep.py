# analyse 2D data of spectrometer 
import pandas as pd
import numpy as np
import math
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
n_diamond = 2.419 #refractive index diamond
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

def shifted_plot_from_2D_data(data_dir):
    wavelengths,filenumbers,intensities = get_data(data_dir)
 
    fig,ax = plt.subplots(figsize=(6,4))
    ax=set_axes(ax,wavelengths)

    for i,intensity in enumerate(np.transpose(intensities)):
        ax.plot(wavelengths,intensity-i*15000)
    
    plt.show()

    try: 
        print os.path.join(data_dir, 'shifted_plots.png')
        fig.savefig(os.path.join(data_dir, 'shifted_plots.png'))
    except:
        print('could not save figure')

    plt.close()

def peaks_from_1D_data(wavelengths,intensity,**kw):
    '''
    function that finds the peaks in 1D data with 
    x = wavelengths,
    y = intensity,
    output parameters:
    x0s - array of xs corresponding to peak locations
    u_x0s - uncertainty in xs from lorentzian fit
    '''
    plot_fit = kw.pop('plot_fit',False)
    plot_peak_locations = kw.pop('plot_peak_locations',False)

    indices,peak_wavelengths,peak_intensity = sa.approximate_peak_location(wavelengths,intensity,**kw)
    x0s,u_x0s,success = sa.fit_peak(wavelengths,intensity,indices,peak_wavelengths,peak_intensity,plot_fit = plot_fit)
    peak_intensity_x0s = peak_intensity[np.where(peak_intensity*success >0)]
    # print len(peak_intensity),len(peak_intensity_x0s)

    if plot_peak_locations:
        fig, ax = plt.subplots()
        ax.plot(wavelengths,intensity)
        # ax.plot(peak_wavelengths,peak_intensity,'o', mfc=None, mec='r', mew=2, ms=8)
        ax.plot(x0s,peak_intensity_x0s,'+', mfc=None, mec='r', mew=2, ms=8)
        plt.show()
        plt.close()

    return x0s,u_x0s    

def peaks_from_2D_data(data_dir,**kw):
    '''
    function that finds the peaks in 2D data in data_dir.
    Output: a 2d scatter plot showing the peak locations.
    TODO: save the obtained peak locations
    '''
    save_fig = kw.pop('save_fig',True)

    wavelengths,filenumbers,intensities = get_data(data_dir)

    x=np.array([])
    u_x = np.array([])
    y=np.array([])

    print 'getting peak locations'
    for i,intensity in enumerate(np.transpose(intensities)):
        print i,'/',len(np.transpose(intensities))
        # if i>0:
        #     break
        x0s,u_x0s = peaks_from_1D_data(wavelengths,intensity,**kw)
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

    if save_fig:
        try: 
            print os.path.join(data_dir, 'peaks.png')
            fig.savefig(os.path.join(data_dir, 'peaks.png'))
        except:
            print('could not save figure')


    return fig,ax

def overlap_peaks_and_modes(data_dir,diamond_thickness=4.e-6,cavity_length = 5.e-6,
        conversion_factor = 307.e-9,nr_points=31, **kw):
    '''
    function that plots the fitted peak locations in 2D data in data_dir, 
    and overlaps it with the analytically derived diamond and air modes.
    Input parameters:
    data_dir - directory containing 2D data
    diamond_thickness - diamond thickness used to obtain analytic result for resonance frequency
    cavity_length - diamond thickness used to obtain analytic result for resonance frequency
    conversion_factor - the piezo conversion factor. at RT:307 nm/V,. at LT:
    nr_points - the number of points used for plotting analytic results of resonances
    **keywords for peak-finding:
    plot_fit - whether to plot the fit of each resonance found in the data

    '''
    fig,ax = peaks_from_2D_data(data_dir,**kw)
    # ax = plot_diamond_modes(diamond_thickness=diamond_thickness,ax = ax)
    # ax = plot_air_modes(cavity_length=cavity_length,diamond_thickness=diamond_thickness,
    #     ax=ax, conversion_factor=conversion_factor,nr_points=nr_points)
    ax = plot_diamond_air_modes(cavity_length=cavity_length,diamond_thickness=diamond_thickness,
        ax=ax,conversion_factor=conversion_factor,nr_points=nr_points)


    title ='d={}um_L={}um_cf={}nmpV'.format(str(diamond_thickness*1e6),str(cavity_length*1.e6),str(conversion_factor*1.e9))

    ax.text(ax.get_xlim()[0] + (ax.get_xlim()[-1]-ax.get_xlim()[0])/4,ax.get_ylim()[0],title)

    #add an axis at the top with the cavity length 
    ax2 = ax.twiny()
    xticks = np.linspace(ax2.get_xlim()[0],ax2.get_xlim()[-1],int((V_max-V_min)/2+1))
    xticklabels2 =np.linspace(cavity_length*1.e6,cavity_length*1.e6+(conversion_factor*(V_max-V_min)*1.e6),int((V_max-V_min)/2+1))
    xticklabels2 = np.round(xticklabels2,2)

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
    max_nr_modes = 100
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
        ax.plot(ax.get_xlim(),[nu,nu], lw=2)
        # ax.text(ax.get_xlim()[-1],nu, 'N={}'.format(N))

    if return_fig:
        return fig,ax

    return ax

def pure_air_modes(cavity_length=1.e-6,conversion_factor = -150.e-9,nr_points=31):
    delta_V = V_max - V_min
    delta_L = delta_V*(conversion_factor) # in m
    Ls = np.linspace(cavity_length,cavity_length+delta_L,nr_points)

    max_nr_modes = 180
    nu_air = np.zeros((max_nr_modes,nr_points))
    for N in np.arange(max_nr_modes):
        for i,L in enumerate(Ls):
            nu_air[N,i] = (N * c / (2 * L))/1.e12 # in THz

    return nu_air

def plot_air_modes(cavity_length=1.e-6,diamond_thickness=0.e-6,ax = None,conversion_factor = -150.e-9,nr_points=31):
    return_fig = False
    if ax == None:
        return_fig = True
        fig,ax = plt.subplots()

    nu_air = pure_air_modes(cavity_length=cavity_length-diamond_thickness,conversion_factor=conversion_factor,nr_points=nr_points)
    xs = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],nr_points)

    for N,nu in enumerate(nu_air):
        ax.plot(xs,nu, lw=2)
        # ax.text(ax.get_xlim()[0],nu[0], 'N={}'.format(N))

    if return_fig:
        return fig,ax

    return ax

def diamond_air_mode_freq(N=1,cavity_length=1.e-6, diamond_thickness=4.e-6):
    Ltot = cavity_length+(n_diamond-1)*diamond_thickness
    Lred = cavity_length-(n_diamond+1)*diamond_thickness
    nu = c / (2*math.pi*Ltot) * \
        (math.pi*N - (-1)**N * math.asin( (n_diamond-1)/(n_diamond+1) * \
        math.sin( (N*math.pi*Lred/Ltot))))
    return nu

def diamond_air_modes(cavity_length = 1.e-6, diamond_thickness = 4.e-6, conversion_factor = 100e-9,nr_points=31):
    delta_V = V_max - V_min
    delta_L = delta_V*(conversion_factor) # in m
    print delta_L
    Ls = np.linspace(cavity_length,cavity_length+delta_L,nr_points)

    max_nr_modes = 180
    nu_diamond_air = np.zeros((max_nr_modes,nr_points))

    for N in np.arange(max_nr_modes):
        for i,L in enumerate(Ls):
            nu_diamond_air[N,i] = diamond_air_mode_freq(N=N,cavity_length=L,diamond_thickness=diamond_thickness)/1.e12 # in THz

    return nu_diamond_air

def plot_diamond_air_modes(cavity_length=1.e-6,diamond_thickness=4.e-6,ax = None,conversion_factor = -150.e-9,nr_points=31):
    return_fig = False
    if ax == None:
        return_fig = True
        fig,ax = plt.subplots()

    nu_diamond_air = diamond_air_modes(cavity_length=cavity_length,diamond_thickness=diamond_thickness,conversion_factor=conversion_factor,nr_points=nr_points)
    xs = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],nr_points)


    for N,nu in enumerate(nu_diamond_air):
        ax.plot(xs,nu, lw=2)
        # if (nu[0]<ax.get_ylim()[-1]) and (nu[0]>ax.get_ylim()[0]):
        #     ax.text(ax.get_xlim()[0],nu[0], 'N={}'.format(N))

    if return_fig:
        return fig,ax

    return ax

