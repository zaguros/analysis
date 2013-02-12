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
    info_xy = kw.pop('info_xy', (0.4, 0.85))
    plot_data = kw.pop('plot_data', False)

    if res == None:
        return False
    
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if plot_data:
        ax.plot(res['x'], res['y'], 'o', mfc='w', mew=1)
    
    ax.plot(fit_xvals, res['fitfunc'](fit_xvals), '-', lw=2)

    if print_info:
        params_str = res['fitfunc_str'] + '\n' + fit.str_fit_params(res)
        plt.figtext(info_xy[0], info_xy[1], params_str, size='x-small',
                color='k', ha='left', va='top', 
                bbox=dict(facecolor='white', alpha=0.5))
    
    if ret == None:
        return
    if ret == 'ax':
        return ax

    return


