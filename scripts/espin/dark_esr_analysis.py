import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot
reload(plot)
### settings
timestamp =None#'20150409204752' #YYYYmmddHHMMSS

guess_offset = 1
guess_x0 = 2.816e9#2.926

guess_splitN = 2.19e-3
guess_splitC = .150e-3 #12.78
guess_splitC2 = .0250e-3 #12.78
guess_width = .5e-3
guess_splitB = 0*0.0008

# guess_splitC = .8e-3 #12.78
#guess_width = 0.2e-3
guess_sigma = guess_width
guess_amplitude = 0.3

# try fitting
#guess_offset = 1.0
guess_A_min1 = 0.0
guess_A_plus1 = 0.0
guess_A_0 = 0.3
#guess_x0 = 3.730
#guess_sigma = 0.435e-3
guess_Nsplit = guess_splitN

def analyze_dark_esr(folder,center_guess = False, ax=None, ret=None,min_dip_depth = 0.9 , **kw):

    if ax == None:
        fig, ax = plt.subplots(1,1)
    ssro_calib_folder = toolbox.latest_data(contains='AdwinSSRO_SSROCalibration')
    print ssro_calib_folder
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)

    x = a.sweep_pts # convert to MHz
    y = a.p0.reshape(-1)[:]
    # ax.plot(x,y)
    a.plot_result_vs_sweepparam(ret=ret, name='ssro', ax=ax)
    # ax.set_ylim(0,1.05)
    # ax.set_ylim([0., 0.10])

    if center_guess == True:
        #guess_x0 = float(raw_input('Center guess?'))
        guess_x0 = 1
    else:
        j=0
        print min_dip_depth
        #print y[21]
        while y[j]>min_dip_depth and j < len(y)-2:  #y[j]>0.93*y[j+1]: # such that we account for noise
            k = j
            j += 1
        #j = len(y)-2
        if k > len(y)-5:
            print 'Could not find dip'
            return
        else:
            print 'k'+str(k)
            print len(y)
            guess_x0 = x[k]+ guess_splitB #convert to GHz and go to middle dip
            guess_ctr=guess_x0 + guess_Nsplit
            print 'guess_ctr= '+str(guess_x0)

    ### fitfunction
    A_min1 = fit.Parameter(guess_A_min1, 'A_min1')
    A_plus1 = fit.Parameter(guess_A_plus1, 'A_plus1')
    A_0 = fit.Parameter(guess_A_0, 'A_0')
    o = fit.Parameter(guess_offset, 'o')
    x0 = fit.Parameter(guess_x0, 'x0')
    sigma = fit.Parameter(guess_sigma, 'sigma')
    Nsplit = fit.Parameter(guess_Nsplit, 'Nsplit')
    def fitfunc(x):
        # return o() - A_min1()*np.exp(-((x-(x0()-splitting-Nsplit()))/sigma())**2) \
        #         - A_min1()*np.exp(-((x-(x0()+splitting-Nsplit()))/sigma())**2) \
        #         - A_plus1()*np.exp(-((x-(x0()-splitting+Nsplit()))/sigma())**2) \
        #         - A_plus1()*np.exp(-((x-(x0()+splitting+Nsplit()))/sigma())**2) \
        #         - A_0()*np.exp(-((x-(x0()+Nsplit()))/sigma())**2) \
        #         - A_0()*np.exp(-((x-(x0()-Nsplit()))/sigma())**2)
        return o() - np.abs(A_min1())*np.exp(-((x-(x0()-Nsplit()))/sigma())**2) \
                - np.abs(A_plus1())*np.exp(-((x-(x0()+Nsplit()))/sigma())**2) \
                - np.abs(A_0())*np.exp(-((x-x0())/sigma())**2) \
    
    try:
        #fit_result = fit.fit1d(x, y, None, p0 = [A_min1, A_plus1, A_0, sigma, o, x0, Nsplit],
        #fitfunc = fitfunc, do_print=True, ret=True, fixed=[0,1,6])
        
        #fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
        #         guess_amplitude, guess_width, guess_ctr,
        #         # (2, guess_splitN),
                  #(2, guess_splitC),
                  #(2, guess_splitC2),
                 # (2, guess_splitB),
                 #(3, guess_splitN),

        #         do_print=True, ret=True, fixed=[0,1,2,3])
        
        plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), ax=ax, plot_data=False, **kw)
        Norm=(fit_result['params'][0]+fit_result['params'][1]+fit_result['params'][2])
        Population_left=fit_result['params'][0]/Norm
        Population_middle=fit_result['params'][2]/Norm
        Population_right=fit_result['params'][1]/Norm
        print '############################'
        print 'Population left ' , Population_left
        print 'Population middle ' , Population_middle
        print 'Population right ' , Population_right
        print '#############################'
    except Exception:

        if center_guess ==True:
            guess_ctr = float(raw_input('Center guess?'))
        print guess_ctr
        fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
                guess_amplitude, guess_width, guess_ctr,
                # (2, guess_splitN),
                # (2, guess_splitC),
                # (2, guess_splitC2),
                #(2, guess_splitB),
                (3, guess_splitN),
                do_print=True, ret=True, fixed=[0])
        plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), ax=ax, plot_data=False, **kw)
        

    print' Minimum is ', np.min(y)
    print 'at', x[np.argmin(y)] , 'GHz'
    
    ax.set_xlabel('MW frq (GHz)')
    ax.set_ylabel(r'fidelity wrt. $|0\rangle$')
    ax.set_title(a.timestamp+'\n'+a.measurementstring)
    ax.set_ylim([0.5,1.05])
    plt.savefig(os.path.join(folder, 'darkesr_analysis.png'),
            format='png')
    #ret='f0'
    if ret == 'f0':
        #print fit_result[]
        f0 = fit_result['params_dict']['x0']
        u_f0 = fit_result['error_dict']['x0']

        ax.text(f0, 0.4, '$f_0$ = ({:.3f} +/- {:.3f})'.format(
            (f0-2.8)*1e3, u_f0*1e3), ha='center')

        return (f0-2.8)*1e3, u_f0*1e3
    return fit_result
### script
if __name__ == '__main__':
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
        if folder==None:
            folder = toolbox.latest_data(timestamp) 
    else:
        folder = toolbox.latest_data('DarkESR')
    print folder
    fit_result=analyze_dark_esr(folder, center_guess= False)








