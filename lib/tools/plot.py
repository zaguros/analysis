import os, sys, time
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from analysis.lib.fitting import fit

# TODO implement formatting options
def plot_fit1d(res, fit_xvals, ax=None, ret=None, **kw):
    '''
    function to plot a fitresult as returned by analysis.lib.fitting.fit.fit1d().
    '''

    # know keywords:
    print_info = kw.pop('print_info', True)
    info_xy = kw.pop('info_xy', 'auto')
    plot_data = kw.pop('plot_data', True)

    if res == None:
        return False
    
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if plot_data:
        
        if 'yerr' in res.keys():
            ax.errorbar(res['x'],res['y'], fmt='o', yerr=res['yerr'])
        else:    
            ax.plot(res['x'], res['y'], 'o')
    
    ax.plot(fit_xvals, res['fitfunc'](fit_xvals), '-', lw=2)

    if print_info:
        params_str = res['fitfunc_str'] + '\n' + fit.str_fit_params(res)
        
        if info_xy == 'auto':
            info_x = ax.get_xlim()[0] + (ax.get_xlim()[1]-ax.get_xlim()[0])*0.02
            info_y = ax.get_ylim()[0] + (ax.get_ylim()[1]-ax.get_ylim()[0])*0.02
        else:
            info_x = info_xy[0]
            info_y = info_xy[1]

        plt.text(info_x, info_y, params_str, size='x-small',
                color='k', ha='left', va='bottom', 
                bbox=dict(facecolor='white', alpha=0.5))
    
    if ret == None:
        return
    if ret == 'ax':
        return ax

    return

#### Density matrix plotting
def plot_dm(dm, **kw):
    dmr = dm.real
    dmi = dm.imag

    fig = plt.figure(figsize=(12,4))
    axr = fig.add_subplot(121, projection='3d')
    axi = fig.add_subplot(122, projection='3d')

    lx= len(dm[0])
    ly= len(dm[:,0])
    xticks = [(r'$| %d \rangle$' % i) for i in range(lx)]
    yticks = [(r'$\langle %d |$' % i) for i in range(ly)]

    plot_mat(dmr, axr, title=r'$\Re(\rho)$', 
            xticks=xticks, yticks=yticks, **kw)
    plot_mat(dmi, axi, title=r'$\Im(\rho)$', 
            xticks=xticks, yticks=yticks, **kw)

def plot_mat(mat, ax=None, **kw):
    '''
    plot a matrix as 3d bar chart.
    '''
    xlabel = kw.pop('xlabel', '')
    ylabel = kw.pop('ylabel', '')
    xticks = kw.pop('xticks', [])
    yticks = kw.pop('yticks', [])
    title = kw.pop('title', '')
    zlabel = kw.pop('zlabel', '')
    ret = kw.pop('ret', None)

    if ax==None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    
    lx= len(mat[0])            # Work out matrix dimensions
    ly= len(mat[:,0])
    xpos = np.arange(0,lx,1)    # Set up a mesh of positions
    ypos = np.arange(0,ly,1)
    xpos, ypos = np.meshgrid(xpos+0.25, ypos+0.25)

    xpos = xpos.flatten()   # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx*ly)

    dx = 0.75 * np.ones_like(zpos)
    dy = dx.copy()
    dz = mat.flatten()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, alpha=0.5, color='RoyalBlue', 
            zsort='average')

    if xticks == []:
        xticks = [str(i) for i in np.arange(lx)]
    if yticks == []:
        yticks = [str(i) for i in np.arange(ly)]
    
    ax.w_xaxis.set_ticks(np.arange(lx) + 0.5)
    ax.w_yaxis.set_ticks(np.arange(ly) + 0.5)
    ax.w_xaxis.set_ticklabels(xticks)
    ax.w_yaxis.set_ticklabels(yticks)
    ax.set_zlim3d([-1,1])
    ax.set_zlabel(zlabel)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if ret == 'ax':
        return ax

    return




