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

def plot_single(contains = '',ind = 0, fit = False, plot_range = None, **kw):

    title = kw.pop('title',None)

    folder = toolbox.latest_data(contains)

    filename = toolbox.measurement_filename(folder, ext = 'dat')

    V,f,c,i,st = np.loadtxt(filename,unpack=True)


    fig = plt.figure()
    ax=plt.subplot(111)

    fltr=i==ind #(dd[:,2]<(y+0.05)) & (dd[:,2]>(j-0.05)) #
    xx=f[np.where(fltr)]
    zz=c[np.where(fltr)]
   
    ax.plot(xx,zz)

    if fit:
        X = np.arange(np.min(xx),np.max(xx),0.005)
        fitfunc = fit_laserscan(xx,zz)
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

def spectral_diffusion(contains = '',print_avg_counts = False, count_filter = None, plot_range = None,**kw):

    title = kw.pop('title',None)

    folder = toolbox.latest_data(contains)

    filename = toolbox.measurement_filename(folder, ext = 'dat')

    V,f,c,i,st = np.loadtxt(filename,unpack=True)

    Y=np.unique(i)

    centres = np.array([])
    errors = np.array([])
    diff = np.array([])
    overlaps = np.array([])

    avg_counts = 0.

    for j,y in enumerate(Y):
        fltr=i==y #(dd[:,2]<(y+0.05)) & (dd[:,2]>(j-0.05)) #
        xx=f[np.where(fltr)]
        zz=c[np.where(fltr)]

        fitfunc = fit_laserscan(xx,zz,do_print=False)

        if count_filter == None or np.sum(zz) > count_filter:

            centres = np.append(centres,fitfunc['params'][2])
            errors =  np.append(errors,fitfunc['params'][3])

            if centres.size != 1:
                diff = np.append(diff, centres[-1] - centres[-2])
                overlaps =  np.append(overlaps,overlap(diff[-1],errors[-1]))

        avg_counts = (avg_counts * j + np.sum(zz))/(j+1)

        #print j,fitfunc['params'][3]
    ys = range(len(centres))

    fig = plt.figure()
    ax=plt.subplot(111)
    ax.errorbar(ys,centres,yerr=errors)
    ax.set_ylabel('Centre')
    ax.set_xlabel('Scan ind')

    if print_avg_counts:
        print 'Avg counts ', avg_counts

    if plot_range != None:
        ax.set_ylim(plot_range)

    if title != None:
        plt.title(title + ' '+os.path.split(filename)[1])
    else:
        plt.title(os.path.split(filename)[1])

    fig = plt.figure()
    ax=plt.subplot(111)
    ax.plot(ys[:-1],1e3*diff)
    ax.set_ylabel('Difference (MHz)')
    ax.set_xlabel('Scan ind')

    fig = plt.figure()
    ax=plt.subplot(111)
    ax.plot(ys[:-1],overlaps)
    ax.set_ylabel('Overlap with previous')
    ax.set_xlabel('Scan ind')

    print "Mean overlap ", np.mean(overlaps)
    

    #save_and_close_plot(folder,name = 'Laserscan_2d')

def save_and_close_plot(f,name = 'Results'):
    
    #plt.savefig(os.path.join(f,name + '.pdf'),format='pdf')
    #plt.savefig(os.path.join(f,name + '.png'),format='png')
    plt.show()
    plt.close('all')

def fit_laserscan(x,y,do_print = True):

    p0, fitfunc, fitfunc_str= common.fit_lorentz(0,np.max(y),np.sum(x * y)/np.sum(y),0.2)

    return fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc,
            fitfunc_str=fitfunc_str, do_print=do_print, ret=True)

def overlap(x_diff,gamma):
        # As laid out in PhysRevLett.100.133601
        lorentz = lambda x,x0,gamma : 2/np.pi*gamma/(4*(x-x0)**2+gamma**2)

        gamma = abs(gamma)

        X = np.linspace(-5*gamma,5*gamma+abs(x_diff),np.round(10*gamma/0.01))

        e_field_0 = np.sqrt(lorentz(X,0,gamma))/np.sqrt(np.sum(lorentz(X,0,gamma)))
        e_field_diff = np.sqrt(lorentz(X,x_diff,gamma))/np.sqrt(np.sum(lorentz(X,x_diff,gamma)))

        return np.sum(e_field_0*e_field_diff)**2

def get_tstamp_from_folder(folder):
    return folder[18:18+15]
