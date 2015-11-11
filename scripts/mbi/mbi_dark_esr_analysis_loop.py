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
timestamp = '20150614_001010' ### one can use this timestamp as starting point.
analysis_length = 39
average_data = True
analysis_type = 'double_lorentz' #### possibilities 'peseudo-voigt' 'gauss' 'lorentz' 'double_gauss', 'double_lorentz'

show_guess = False

folder_list = []

if timestamp == None:
    timestamp = '20800101_101010' ## some future date... only younger data is considered.

for i in range(analysis_length):
    timestamp, folder = toolbox.latest_data(contains='DESR',
                                        return_timestamp =True,
                                        older_than=timestamp,
                                        raise_exc=False)
    folder_list.append(folder)

# print folder_list

guess_offset = 1.00
guess_A_min1 = 0.07
guess_A_plus1 = 0.02
guess_A_0 = 0.02
guess_x0 = 1746.668427
guess_sigma = 0.250/2.
guess_Nsplit = 2.195435#2.182 #

### needed for a linearized pseudo-voigt profile
guess_eta = 0.5


### fitfunction
A_min1 = fit.Parameter(guess_A_min1, 'A_min1')
A_plus1 = fit.Parameter(guess_A_plus1, 'A_plus1')
A_0 = fit.Parameter(guess_A_0, 'A_0')
o = fit.Parameter(guess_offset, 'o')
x0 = fit.Parameter(guess_x0, 'x0')
sigma = fit.Parameter(guess_sigma, 'sigma')
Nsplit = fit.Parameter(guess_Nsplit, 'Nsplit')


if analysis_type == 'pseudo-voigt':
    eta = fit.Parameter(guess_eta, 'eta')



    def fitfunc(x):
        s = 0.727981*sigma() ### conversion factor gaussian/lorentzian

        lorentzm1 = .5*s/((x-x0()+Nsplit())**2+(0.5*s)**2)/np.pi
        lorentz0 = .5*s/((x-x0())**2+(0.5*s)**2)/np.pi
        lorentzp1 = .5*s/((x-x0()-Nsplit())**2+(0.5*s)**2)/np.pi

        F1 =  - 1.826*A_min1()*lorentzm1 - 1.826*A_plus1()*lorentzp1 - 1.826*A_0()*lorentz0

        gaussm1 = - A_min1()*np.exp(-((x-(x0()-Nsplit()))/sigma())**2)
        gauss0 = - A_plus1()*np.exp(-((x-(x0()+Nsplit()))/sigma())**2)
        gaussp1 = - A_0()*np.exp(-((x-x0())/sigma())**2)

        F2 = gaussm1 + gauss0 + gaussp1

        return o() + eta()*F1 + (1-eta())*F2
    
    

elif analysis_type == 'lorentz':

    def fitfunc(x):
        lorentzm1 = .5*sigma()/((x-x0()+Nsplit())**2+(0.5*sigma())**2)/np.pi
        lorentz0 = .5*sigma()/((x-x0())**2+(0.5*sigma())**2)/np.pi
        lorentzp1 = .5*sigma()/((x-x0()-Nsplit())**2+(0.5*sigma())**2)/np.pi

        return o() - A_min1()*lorentzm1 - A_plus1()*lorentzp1 - A_0()*lorentz0

elif analysis_type == 'double_gauss':

    guess_Csplit = 0.917e-3
    Csplit = fit.Parameter(guess_Csplit, 'Csplit')
    A_min12 = fit.Parameter(guess_A_min1, 'A_min12')
    A_plus12 = fit.Parameter(guess_A_plus1, 'A_plus12')
    A_02 = fit.Parameter(guess_A_0, 'A_02')

    def fitfunc(x):
        
        gaussm1 = - A_min1()*np.exp(-((x-(x0()-Nsplit()-Csplit()))/sigma())**2) - A_min12()*np.exp(-((x-(x0()-Nsplit()+Csplit()))/sigma())**2)
        gauss0 = - A_plus1()*np.exp(-((x-(x0()+Nsplit()-Csplit()))/sigma())**2)- A_plus12()*np.exp(-((x-(x0()+Nsplit()+Csplit()))/sigma())**2)
        gaussp1 = - A_0()*np.exp(-((x-x0()-Csplit())/sigma())**2)- A_02()*np.exp(-((x-(x0()+Csplit()))/sigma())**2)

        return o() + gaussm1 + gauss0 + gaussp1


elif analysis_type == 'double_lorentz':

    guess_Csplit = 180e-3
    Csplit = fit.Parameter(guess_Csplit, 'Csplit')
    A_min12 = fit.Parameter(guess_A_min1, 'A_min12')
    A_plus12 = fit.Parameter(guess_A_plus1, 'A_plus12')
    A_02 = fit.Parameter(guess_A_0, 'A_02')

    def fitfunc(x):
        s = sigma()
        lorentzm1 = -  A_min1()*.5*s/((x-x0()+Nsplit()-Csplit())**2+(0.5*s)**2)/np.pi - A_min12()*.5*s/((x-x0()+Nsplit()+Csplit())**2+(0.5*s)**2)/np.pi

        lorentzp1 = -  A_plus1()*.5*s/((x-x0()-Nsplit()-Csplit())**2+(0.5*s)**2)/np.pi - A_plus12()*.5*s/((x-x0()-Nsplit()+Csplit())**2+(0.5*s)**2)/np.pi

        lorentz0 = -  A_0()*.5*s/((x-x0()-Csplit())**2+(0.5*s)**2)/np.pi - A_02()*.5*s/((x-x0()+Csplit())**2+(0.5*s)**2)/np.pi
        
        return o() + lorentzm1 + lorentzp1 + lorentz0

else: ### gaussian analysis
    def fitfunc(x):

        csplit = 1.07e-6

        gaussm1 = - A_min1()*(np.exp(-((x-(x0()-Nsplit()) + csplit)/sigma())**2) + np.exp(-((x-(x0()-Nsplit()) - csplit)/sigma())**2))
        gauss0 =  - A_0()*(np.exp(-((x-x0() + csplit)/sigma())**2) + np.exp(-((x-x0() - csplit)/sigma())**2))
        gaussp1 = - A_plus1()*(np.exp(-((x-(x0()+Nsplit()+csplit))/sigma())**2) + np.exp(-((x-(x0()+Nsplit()-csplit))/sigma())**2))
        return o() + gaussm1 + gauss0 + gaussp1
    



###################################################

# get data

###################################################

cum_pts = 0
cum_sweep_pts       = np.empty(0)
cum_p0              = np.empty(0)
cum_u_p0              = np.empty(0)

cum_normalized_ssro = np.empty(0)

for kk,folder in enumerate(folder_list):

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    cum_pts += a.pts

    if kk == 0:
        cum_sweep_pts = a.sweep_pts
        if average_data:
            cum_p0 = a.p0/float(analysis_length)
            cum_u_p0 = a.u_p0**2
        else:
            cum_p0 = a.p0
            cum_u_p0 = a.u_p0

        # reps_per_datapoint = a.reps
    else:
        if average_data:
            cum_p0 = a.p0/float(analysis_length) + cum_p0
            cum_u_p0 = a.u_p0**2+ cum_u_p0
        else:
            cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
            cum_p0 = np.concatenate((cum_p0, a.p0))
            cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))
        
if average_data:
    cum_u_p0 = np.sqrt(cum_u_p0)/float(analysis_length)
#sorting_order=cum_sweep_pts.argsort()
#cum_sweep_pts.sort()
#cum_p0=cum_p0[sorting_order]
#cum_u_p0=cum_u_p0[sorting_order]

a.pts   = cum_pts
a.sweep_pts = cum_sweep_pts
a.p0    = cum_p0
a.u_p0  = cum_u_p0


#########################################################

# fit and plot

#########################################################


ax = a.plot_results_vs_sweepparam(ret='ax',name='ssro',fmt='-')
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
    if show_guess:
        ax.plot(np.linspace(x[0],x[-1],501), fitfunc(np.linspace(x[0],x[-1],501)), '--', lw=2)

    if analysis_type == 'pseudo-voigt':
        fixed = [3,4,6]
        fit_result = fit.fit1d(x, y, None, p0 = [A_min1, A_plus1, A_0, o, x0, sigma, Nsplit,eta],
                fitfunc = fitfunc, do_print=True, ret=True, fixed=fixed)


    elif analysis_type == 'lorentz':
        fixed = [3,4,6]
        fit_result = fit.fit1d(x, y, None, p0 = [A_min1, A_plus1, A_0, o, x0, sigma, Nsplit],
                fitfunc = fitfunc, do_print=True, ret=True, fixed=fixed)

    elif 'double' in analysis_type:
        fixed = [6]
        fit_result = fit.fit1d(x, y, None, p0 = [A_min1, A_min12,A_plus1,A_plus12, A_0,A_02, o, x0, sigma, Nsplit, Csplit],
                fitfunc = fitfunc, do_print=True, ret=True, fixed=fixed)

    else:
        fixed = [4,6]
        fit_result = fit.fit1d(x, y, None, p0 = [A_min1, A_plus1, A_0, o, x0, sigma, Nsplit],
                fitfunc = fitfunc, do_print=True, ret=True, fixed=fixed)

    


    plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), plot_data=False, ax=ax)

    if not 'double' in analysis_type:
        Norm=(fit_result['params_dict']['A_min1']+fit_result['params_dict']['A_0']+fit_result['params_dict']['A_plus1'])
        Norm_error=np.sqrt((fit_result['error_dict']['A_min1']**2+fit_result['error_dict']['A_0']**2+fit_result['error_dict']['A_plus1']**2))

        print fit_result['params_dict']
        Population_left=fit_result['params_dict']['A_min1']/Norm
        Population_left_error = np.sqrt((1/Norm)**2*fit_result['error_dict']['A_min1']**2+(fit_result['params_dict']['A_min1']/Norm**2)**2*Norm_error**2)
        Population_middle=fit_result['params_dict']['A_0']/Norm
        Population_right=fit_result['params_dict']['A_plus1']/Norm

    else:
        Norm=(fit_result['params_dict']['A_min1']+fit_result['params_dict']['A_min12']+fit_result['params_dict']['A_0']+fit_result['params_dict']['A_02']+fit_result['params_dict']['A_plus1']+fit_result['params_dict']['A_plus12'])
        Norm_error=np.sqrt((fit_result['error_dict']['A_min1']**2+fit_result['error_dict']['A_min12']**2+fit_result['error_dict']['A_0']**2+fit_result['error_dict']['A_02']**2+fit_result['error_dict']['A_plus1']**2+fit_result['error_dict']['A_plus12']**2))

        print fit_result['params_dict']
        Population_left=(fit_result['params_dict']['A_min1']+fit_result['params_dict']['A_min12'])/Norm
        Population_left_error = np.sqrt((1/Norm)**2*(fit_result['error_dict']['A_min1']+fit_result['error_dict']['A_min12'])**2+((fit_result['params_dict']['A_min1']+fit_result['params_dict']['A_min12'])/Norm**2)**2*Norm_error**2)
        Population_middle=(fit_result['params_dict']['A_0']+fit_result['params_dict']['A_02'])/Norm
        Population_right=(fit_result['params_dict']['A_plus1']+fit_result['params_dict']['A_plus12'])/Norm

    ax.set_ylim(0.6,1.1)
    # ax.set_xlim(1.743e3,1.749e3)

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


