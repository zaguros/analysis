
import sys
import os
sys.path.append("d:/measuring")
import numpy as np
import h5py

from matplotlib import pyplot as plt
import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(mbi)
from analysis.lib.tools import toolbox as tb
from analysis.lib.m2 import m2
from analysis.lib.m2.ssro import ssro, mbi, sequence, pqsequence
from analysis.lib.nv import nvlevels
from analysis.lib.lde import tail_cts_per_shot_v4 as tail
from analysis.lib.pq import pq_tools, pq_plots
from analysis.lib.math import readout_correction as roc
from analysis.lib.math import error
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
reload(m2)
reload(tb)
reload(ssro)
reload(mbi)
reload(sequence)
reload(pqsequence)
reload(tail)
reload(pq_tools)
reload(pq_plots)



def fit_ramsey(title=None, timestamp = None, ssro_calib_timestamp = None, ax = None,do_fit = False, Rmbi_guess=0.9, theta_guess= 1.5705, phi_guess=1.5705, show_guess = False):


    """
    fit Rmbi**2Cos(theta)**2 + Rmbi**2Sin(theta)**2(cos(x+phi)+sin(x+phi)
    """
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data(title)

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')

    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'
        print ssro_calib_folder

    fig, ax = plt.subplots(1,1, figsize=(4.5,4))
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()    

    x = a.sweep_pts.reshape(-1)[:]

    y= ((a.p0.reshape(-1))-0.5)*2
    x = np.linspace(0,720,len(y))

    y_err = a.u_p0.reshape(-1)[:]


    Rmbi = fit.Parameter(Rmbi_guess, 'Rmbi')
    theta = fit.Parameter(theta_guess, 'theta')
    phi = fit.Parameter(phi_guess, 'phi')
    
    p0 = [Rmbi,theta, phi]

    fitfunc_str = 'Rmbi^2 * Cos(theta)^2 + Rmbi^2 * Sin(theta)^2 * (Cos(x+phi)+Sin(x+phi)'

    ax.errorbar(x,y,yerr=y_err)

    def fitfunc(x) : 
        # return (Rmbi()*np.cos(theta()+np.pi*x/180.)+0*phi())
        return (Rmbi()**2)*(np.cos(theta()*np.pi/180.)**2 + np.sin(theta()*np.pi/180.)**2 * (np.cos(x*np.pi/180.+np.pi*phi()/180.)+np.sin(x*np.pi/180.+np.pi*phi()/180.)))

    if show_guess:
        ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)

    if do_fit == True:
        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc,
            fitfunc_str=fitfunc_str, do_print=True, ret=True,fixed = [])

        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
            plot_data=False, print_info=False)

    ax.set_title(folder)
    ax.set_ylabel('Contrast')
    ax.set_xlabel('Phase')

    ax.set_ylim([-1.1,1.1])
    # print 'fitted minimal fidelity {0:.3f} +- {1:.3f}'.format((1-fit['params'][0]),fit['error'][0])
    # print 'Fitted minimum at {0:.3f} +- {1:.3f}'.format((fit['params'][2]),fit['error'][2])