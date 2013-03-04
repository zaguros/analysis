import os, sys, time
import matplotlib
from matplotlib import pyplot as plt

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
            ax.errorbar(res['x'],res['y'],fmt='o',yerr=res['yerr'])
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


