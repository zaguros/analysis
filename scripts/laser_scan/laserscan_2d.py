'''
Updated script to analyze linescan data using our more modern scripts
2016 PH
'''

import numpy as np
import os
from analysis.lib.tools import toolbox, plot; 
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
from analysis.lib.fitting import fit, common
import copy as cp
from scipy.interpolate import interp1d

def plot_freq_against_index(contains = '',plot_range = None, interp_res = 0.01, **kw):

    title = kw.pop('title',None)

    folder = toolbox.latest_data(contains)
    filename = toolbox.measurement_filename(folder, ext = 'dat')

    V,f,c,i,st = np.loadtxt(filename,unpack=True)

    Y=np.unique(i)
    if plot_range == None:
        X = np.arange(np.min(f),np.max(f),interp_res)
    else:
        X = np.arange(plot_range[0],plot_range[1],interp_res)
    Z=np.zeros((len(X), len(Y)))

    for j,y in enumerate(Y):
        fltr=i==y #(dd[:,2]<(y+0.05)) & (dd[:,2]>(j-0.05)) #
        xx=f[np.where(fltr)]
        zz=c[np.where(fltr)]
        interp_func=interp1d(xx,zz, bounds_error=False, fill_value=0.)
        Z[:,j]=np.nan_to_num(interp_func(X))
    
    fig = plt.figure()
    ax=plt.subplot(111)

    ax.imshow(Z, aspect='auto', origin='lower',vmin=0,vmax = np.max(Z),extent=[min(Y),max(Y),X[0],X[-1]], cmap='binary')
    ax.set_xlabel('Scan nr')
    ax.set_ylabel('Laser frequency [GHz]')

    if title != None:
        plt.title(title + ' '+os.path.split(filename)[1])
    else:
        plt.title(os.path.split(filename)[1])

    save_and_close_plot(folder,name = 'Laserscan_2d')

def plot_overlaid(contains = '',show_indiv = False,fit = False, compensate_drift = False, plot_range = None, interp_res = 0.005, **kw):

    title = kw.pop('title',None)

    folder = toolbox.latest_data(contains)

    filename = toolbox.measurement_filename(folder, ext = 'dat')

    V,f,c,i,st = np.loadtxt(filename,unpack=True)

    if not show_indiv:
        if plot_range == None:
            if compensate_drift:
                X = np.arange(np.min(f)-np.mean(f),np.max(f)-np.mean(f),interp_res)
            else:
                X = np.arange(np.min(f),np.max(f),interp_res)
                
        else:
            X = np.arange(plot_range[0],plot_range[1],interp_res)
        Z = np.zeros(len(X))

    Y=np.unique(i)

    fig = plt.figure()
    ax=plt.subplot(111)


    for j,y in enumerate(Y):
        fltr=i==y #(dd[:,2]<(y+0.05)) & (dd[:,2]>(j-0.05)) #
        xx=f[np.where(fltr)]
        zz=c[np.where(fltr)]
        if compensate_drift:
            avg_pos = np.sum(xx * zz)/np.sum(zz)
            xx = xx - avg_pos

        if show_indiv:
            ax.plot(xx,zz)
        else:
            interp_func=interp1d(xx,zz, bounds_error=False, fill_value=0.)
            Z = Z + np.nan_to_num(interp_func(X))
    
    if not show_indiv:
        ax.plot(X,Z)

        if fit:
            fitfunc = fit_laserscan(X,Z)
            ax.plot(X,fitfunc['fitfunc'](X))

    if plot_range != None:
        ax.set_xlim(plot_range)

    ax.set_ylabel('Counts')
    ax.set_xlabel('Laser frequency [GHz]')

    if title != None:
        plt.title(title + ' '+os.path.split(filename)[1])
    else:
        plt.title(os.path.split(filename)[1])

    #save_and_close_plot(folder,name = 'Laserscan_2d')

def save_and_close_plot(f,name = 'Results'):
    
    #plt.savefig(os.path.join(f,name + '.pdf'),format='pdf')
    #plt.savefig(os.path.join(f,name + '.png'),format='png')
    plt.show()
    plt.close('all')

def fit_laserscan(x,y):

    p0, fitfunc, fitfunc_str= common.fit_lorentz(0,np.max(y),np.sum(x * y)/np.sum(y),0.2)

    return fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc,
            fitfunc_str=fitfunc_str, do_print=True, ret=True)

def get_tstamp_from_folder(folder):
    return folder[18:18+15]
