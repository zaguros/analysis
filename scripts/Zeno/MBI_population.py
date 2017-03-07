"""
This script loops over a series of DarkESR measurements (with initialized nitrogen)
and analyses the population of the nitrogen after MBI via a double lorentzian fit per nitrogen state (strongly coupled carbon)
NK 2015
"""

import os, sys
import numpy as np
import h5py
import logging
import sympy
import matplotlib as mpl
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, esr
from analysis.lib.tools import plot
from analysis.lib.math import error




### plot settings:

figscaling = 0.66*0.72
mpl.rcParams['font.sans-serif'] = 'Lucida Grande'
mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['mathtext.rm'] = 'sans'
# mpl.rcParams['mathtext.default'] = 'rm'
# mpl.rcParams['lines.solid_joinstyle'] = 
# mpl.rcParams['errorbar.capsize'] = 3
# font = {'family' : 'sans-serif',
#     'weight' : 'normal',
#     'size'   : 30}
GlobalErrorCapsize = 0
fit_lw = 2 *figscaling
bar_width = 0.5

## writes a bunch of parameters to matplotlib config.
## input: effective figure scaling factor
mpl.rcParams['axes.linewidth'] = 1.8*figscaling
mpl.rcParams['xtick.major.width'] = 1.8*figscaling
mpl.rcParams['ytick.major.width'] = 1.8*figscaling
mpl.rcParams['font.size'] = (22-8*(figscaling-0.66*0.72)/(1-0.66*0.72))*figscaling
mpl.rcParams['axes.titlesize'] = mpl.rcParams['font.size']
mpl.rcParams['legend.fontsize'] = mpl.rcParams['font.size']
mpl.rcParams['legend.labelspacing'] = 0.5*figscaling
mpl.rcParams['legend.columnspacing'] = 1.5*figscaling
mpl.rcParams['legend.handletextpad'] = 0.3*figscaling
mpl.rcParams['legend.handlelength'] = 1.*figscaling
mpl.rcParams['legend.borderpad'] = 0.2+0.2*(figscaling-0.66*0.72)/(1-0.66*0.72)
mpl.rcParams['lines.markersize'] = 7*figscaling #was 6
mpl.rcParams['figure.figsize'] = (6*figscaling,4*figscaling)
mpl.rcParams['lines.markeredgewidth'] = 0.3/figscaling

### settings
timestamp = '20150613_093117'#'20150609_231010' ### one can use this timestamp as starting point.
ssro_calib_timestamp = '20150610_174113'

save_folder = r'D:/measuring/data/Zeno_results'

analysis_length = 39 ### how many separate measurements did you take?
average_data = True
analysis_type = 'double_lorentz' 

show_guess = False


ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'

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
guess_A_plus1 = 0.002
guess_A_0 = 0.002
guess_x0 = 1746.663387
guess_sigma = 0.250/2.
guess_Nsplit = 2.195433#2.182 #
guess_Csplit = 91.681e-3



### fitfunction
   
def Npop_double_lorentz_fit(g_Csplit,g_Nsplit,g_offset,g_x0,g_sigma,g_A_min1,g_A_0,g_A_plus1):

    fitfunc_str = "Fit 6 lorentzians to MBI data"

    Csplit = fit.Parameter(g_Csplit, 'Csplit')
    A_min12 = fit.Parameter(g_A_min1, 'A_min12')
    A_plus12 = fit.Parameter(g_A_plus1, 'A_plus12')
    A_02 = fit.Parameter(g_A_0, 'A_02')
    A_min1 = fit.Parameter(g_A_min1, 'A_min1')
    A_plus1 = fit.Parameter(g_A_plus1, 'A_plus1')
    A_0 = fit.Parameter(g_A_0, 'A_0')
    o = fit.Parameter(g_offset, 'o')
    x0 = fit.Parameter(g_x0, 'x0')
    sigma = fit.Parameter(g_sigma, 'sigma')
    Nsplit = fit.Parameter(g_Nsplit, 'Nsplit')


    p0 = [Csplit,Nsplit,o,x0,sigma,A_min1,A_min12,A_0,A_02,A_plus1,A_plus12]

    def fitfunc(x):
        s = sigma()
        lorentzm1 = -  A_min1()*.5*s/((x-x0()+Nsplit()-Csplit())**2+(0.5*s)**2)/np.pi - A_min12()*.5*s/((x-x0()+Nsplit()+Csplit())**2+(0.5*s)**2)/np.pi

        lorentzp1 = -  A_plus1()*.5*s/((x-x0()-Nsplit()-Csplit())**2+(0.5*s)**2)/np.pi - A_plus12()*.5*s/((x-x0()-Nsplit()+Csplit())**2+(0.5*s)**2)/np.pi

        lorentz0 = -  A_0()*.5*s/((x-x0()-Csplit())**2+(0.5*s)**2)/np.pi - A_02()*.5*s/((x-x0()+Csplit())**2+(0.5*s)**2)/np.pi
        
        return o() + lorentzm1 + lorentzp1 + lorentz0
    return p0, fitfunc, fitfunc_str

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
    a.get_electron_ROC(ssro_calib_folder)
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


# ax = a.plot_results_vs_sweepparam(ret='ax',name='ssro',fmt='-')
x = a.sweep_pts
y = a.p0.reshape(-1)
wrong_population = None


fig = plt.figure(figsize=[7/1.5,5/1.5])
ax = plt.subplot()

### manipulate the x-axis and switch to a detuning axis in units of MHz:
x = x - 1746.668427

plt.errorbar(x,y,a.u_p0,fmt =None)
### start fitting

guess_offset = 1.00
guess_A_min1 = 0.07
guess_A_plus1 = 0.002
guess_A_0 = 0.002
guess_x0 = 0.0
guess_sigma = 0.250/2.
guess_Nsplit = 2.195435#2.182 #
guess_Csplit = 91.681e-3

p0, fitfunc,fitfunc_str = Npop_double_lorentz_fit(guess_Csplit,guess_Nsplit,guess_offset,guess_x0,guess_sigma,guess_A_min1,guess_A_0,guess_A_plus1)
# try fitting
if show_guess:
    ax.plot(np.linspace(x[0],x[-1],501), fitfunc(np.linspace(x[0],x[-1],501)), '--', lw=1)



fixed = [2,3]
fit_result = fit.fit1d(x, y, None, p0 = p0,
        fitfunc = fitfunc, do_print=True, ret=True, fixed=fixed)



plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), plot_data=False, ax=ax,add_txt = False,lw=1)




Norm=(fit_result['params_dict']['A_min1']+fit_result['params_dict']['A_min12']+fit_result['params_dict']['A_0']+fit_result['params_dict']['A_02']+fit_result['params_dict']['A_plus1']+fit_result['params_dict']['A_plus12'])
Norm_error=np.sqrt((fit_result['error_dict']['A_min1']**2+fit_result['error_dict']['A_min12']**2+fit_result['error_dict']['A_0']**2+fit_result['error_dict']['A_02']**2+fit_result['error_dict']['A_plus1']**2+fit_result['error_dict']['A_plus12']**2))

print fit_result['params_dict']
Population_left=(fit_result['params_dict']['A_min1']+fit_result['params_dict']['A_min12'])/Norm
Population_left_error = np.sqrt((1/Norm)**2*(fit_result['error_dict']['A_min1']+fit_result['error_dict']['A_min12'])**2+((fit_result['params_dict']['A_min1']+fit_result['params_dict']['A_min12'])/Norm**2)**2*Norm_error**2)
Population_middle=(fit_result['params_dict']['A_0']+fit_result['params_dict']['A_02'])/Norm
Population_middle_error = np.sqrt((1/Norm)**2*(fit_result['error_dict']['A_0']+fit_result['error_dict']['A_02'])**2+((fit_result['params_dict']['A_0']+fit_result['params_dict']['A_02'])/Norm**2)**2*Norm_error**2)

Population_right=(fit_result['params_dict']['A_plus1']+fit_result['params_dict']['A_plus12'])/Norm
Population_right_error = np.sqrt((1/Norm)**2*(fit_result['error_dict']['A_plus1']+fit_result['error_dict']['A_plus12'])**2+((fit_result['params_dict']['A_plus1']+fit_result['params_dict']['A_plus12'])/Norm**2)**2*Norm_error**2)

ax.set_ylim(0.58,1.02)
ax.set_xlim(-3,3)
ax.set_yticks([0.6,0.8,1.0])
# ax.set_xticks([-2,0,2])
ax.set_xlabel('Detuning (MHz)')
ax.set_ylabel(r'p($m_s = 0$)')
# ax.set_xlim(1.743e3,1.749e3)

plt.savefig(os.path.join(save_folder, 'mbi_Nitrogen_population.pdf'),
        format='pdf')
plt.savefig(os.path.join(save_folder, 'mbi_Nitrogen_population.png'),
        format='png')

print "plot is saved in",save_folder

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
print 'Population middle error' , Population_middle_error
print 'Population right ' , Population_right
print 'Population right error ' , Population_right_error
print '#############################'
if (not(wrong_population)==None):
    print 'Non-initialized population: ', wrong_population
    print 'Electron initialization:  ', off


