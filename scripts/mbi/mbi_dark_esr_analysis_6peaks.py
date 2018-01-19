import os, sys
import numpy as np
import h5py
import logging
import sympy

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, esr
from analysis.lib.tools import plot
from analysis.lib.math import error

### settings
timestamp = '20171204/222919' #None#'20151102_125238'#'20150613_093116' # '20150324_213946'#  '120943'#'171251'#None #'190948' #

guess_offset = 1
guess_A_min1 = 0.1
guess_A_plus1 = 0.1
guess_A_0 = 0.6
guess_x0 = 1746.5
guess_sigma = .2
guess_Nsplit = 2.179
guess_Csplit = 0.190

### set the ssro_calib folder

ssro_folder = toolbox.latest_data('SSROCalibration', folder = r'/Volumes/staff-bulk-tnw-ns/qt/Taminiaulab/Manual Backup 20171204/data/')


### fitfunction
A_min1_pC13 = fit.Parameter(guess_A_min1, 'A_min1_pC13')
A_min1_mC13 = fit.Parameter(guess_A_min1, 'A_min1_mC13')
A_plus1_pC13 = fit.Parameter(guess_A_plus1, 'A_plus1_pC13')
A_plus1_mC13 = fit.Parameter(guess_A_plus1, 'A_plus1_mC13')
A_0_pC13 = fit.Parameter(guess_A_0, 'A_0_pC13')
A_0_mC13 = fit.Parameter(guess_A_0, 'A_0_mC13')
o = fit.Parameter(guess_offset, 'o')
x0 = fit.Parameter(guess_x0, 'x0')
sigma = fit.Parameter(guess_sigma, 'sigma')
Nsplit = fit.Parameter(guess_Nsplit, 'Nsplit')
Csplit = fit.Parameter(guess_Csplit, 'Csplit')

def fitfunc(x):
    return o() - A_min1_mC13()*np.exp(-((x-(x0()-Nsplit()-Csplit()))/sigma())**2) \
            - A_min1_pC13()*np.exp(-((x-(x0()-Nsplit()+Csplit()))/sigma())**2) \
            - A_plus1_mC13()*np.exp(-((x-(x0()+Nsplit()-Csplit()))/sigma())**2) \
            - A_plus1_pC13()*np.exp(-((x-(x0()+Nsplit()+Csplit()))/sigma())**2) \
            - A_0_pC13()*np.exp(-((x-x0()+Csplit())/sigma())**2) \
            - A_0_mC13()*np.exp(-((x-x0()-Csplit())/sigma())**2) \
### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp, folder = r'/Volumes/staff-bulk-tnw-ns/qt/Taminiaulab/Manual Backup 20171204/data/')
    print folder

else:
    #folder = toolbox.latest_data('PostInitDarkESR')
    folder = toolbox.latest_data('DESR')


a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')

a.get_electron_ROC(ssro_calib_folder = ssro_folder)
ax = a.plot_results_vs_sweepparam(ret='ax',name='ssro')
x = a.sweep_pts
y = a.p0.reshape(-1)
wrong_population = None

if (len(y)==5):
    off = (y[0]+y[4])*0.5
    a_p1 = off-y[1]
    a_0 = off-y[2]
    a_m1 = off-y[3]
    max_pop = min(y[1], y[2], y[3])
    norm = a_p1+a_0+a_m1
    Population_left=a_p1/norm
    Population_middle=a_0/norm
    Population_right=a_m1/norm
    wrong_population = max_pop

else:
    # try fitting

    fit_result = fit.fit1d(x, y, None, p0 = [A_min1_pC13,A_min1_mC13, A_plus1_pC13, A_plus1_mC13, A_0_pC13, A_0_mC13, o, x0, sigma, Nsplit, Csplit],
            fitfunc = fitfunc, do_print=True, ret=True, fixed=[], VERBOSE = True)



    plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), plot_data=False, ax=ax, do_print = True)

    # print fit_result
    
    Norm=(fit_result['params'][0]+fit_result['params'][1]+fit_result['params'][2]+fit_result['params'][3]+fit_result['params'][4]+fit_result['params'][5])
    Norm_error=np.sqrt((fit_result['error'][0]**2+fit_result['error'][1]**2+fit_result['error'][2]**2+fit_result['error'][3]**2+fit_result['error'][4]**2+fit_result['error'][5]**2))

    Population_left=(fit_result['params'][0]+fit_result['params'][1])/Norm
    Population_left_error = np.sqrt((1/Norm)**2*fit_result['error'][0]**2+(fit_result['params'][0]/Norm**2)**2*Norm_error**2)
    Population_middle=(fit_result['params'][4]+fit_result['params'][5])/Norm
    Population_right=(fit_result['params'][2]+fit_result['params'][3])/Norm
    
    ax.set_ylim(-0.05,1.1)


    plt.savefig(os.path.join(folder, 'mbi_darkesr_analysis.pdf'),
            format='pdf')
    plt.savefig(os.path.join(folder, 'mbi_darkesr_analysis.png'),
            format='png')

    '''
    pol = error.Formula()
    a0, am1, ap1 = sympy.symbols('a0, am1, ap1')
    pol.formula = am1 / (a0 + ap1 + am1)
    pol.values[a0] = A_0()
    pol.values[am1] = A_min1()
    pol.values[ap1] = A_plus1()
    pol.uncertainties[a0] = fit_result['error_dict']['A_0']
    pol.uncertainties[am1] = fit_result['error_dict']['A_min1']
    pol.uncertainties[ap1] = fit_result['error_dict']['A_plus1']

    print 'Spin polarization = %.3f +/- %.3f' \
            % (float(pol.value()), float(pol.uncertainty()))
    '''

print '############################'
print 'Population left ' , Population_left
print 'population left error ',  Population_left_error
print 'Population middle ' , Population_middle
print 'Population right ' , Population_right
print '#############################'
if (not(wrong_population)==None):
    print 'Non-initialized population: ', wrong_population
    print 'Electron initialization:  ', off


