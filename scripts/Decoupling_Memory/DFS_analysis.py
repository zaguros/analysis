'''Written by MAB 10-3-15 for a general coherence msmt with a "free" exponent'''

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
from math import floor, log10



if os.name == 'posix':
        DBdir = r'/Users/'+os.getlogin()+r'/Dropbox/QEC LT/Decoupling memory/General_Data/'
else:
    DBdir = r'D:/Dropbox/QEC LT/Decoupling memory/General_Data/'

reload(common)
reload(mbi)
from scipy import optimize

C36_DFS={
'clas |xx>' : 'From20150504143809_To20150504152559_NuclearDD_111_1_sil18_sweep_evolution_time_auto_C3&6_TomoXX_ROpositive_el1_N01_part.txt',
'Logic +X' : 'From20150502002510_To20150502015926_NuclearDD_111_1_sil18_sweep_evolution_time_auto_C3&6_LogicpX_TomoXX_ROpositive_el1_N01_part.txt',
'Logic -X' : 'From20150502021521_To20150502034937_NuclearDD_111_1_sil18_sweep_evolution_time_auto_C3&6_LogicmX_TomoXX_ROpositive_el1_N01_part.txt',}

C36_DFS={
# r'$|xx\rangle$,        ' : 'From20150504143809_To20150504152559_NuclearDD_111_1_sil18_sweep_evolution_time_auto_C3&6_TomoXX_ROpositive_el1_N01_part.txt',
r'$|00\rangle+|11\rangle$' : 'From20150502002510_To20150502015926_NuclearDD_111_1_sil18_sweep_evolution_time_auto_C3&6_LogicpX_TomoXX_ROpositive_el1_N01_part.txt',
r'$|01\rangle+|10\rangle$' : 'From20150502021521_To20150502034937_NuclearDD_111_1_sil18_sweep_evolution_time_auto_C3&6_LogicmX_TomoXX_ROpositive_el1_N01_part.txt'}

C26_DFS={
 'clas |xx>' : 'From20150504163050_To20150504171904_NuclearDD_111_1_sil18_sweep_evolution_time_auto_C2&6_TomoXX_ROpositive_el1_N01_part.txt',
 'Logic +X' : 'From20150504043413_To20150504052522_NuclearDD_111_1_sil18_sweep_evolution_time_auto_C2&6_LogicpX_TomoXX_ROpositive_el1_N01_part.txt',
 'Logic -X' : 'From20150504053437_To20150504062550_NuclearDD_111_1_sil18_sweep_evolution_time_auto_C2&6_LogicmX_TomoXX_ROpositive_el1_N01_part.txt',}

C26_DFS={
 'Logic +X' : 'From20150504043413_To20150504052522_NuclearDD_111_1_sil18_sweep_evolution_time_auto_C2&6_LogicpX_TomoXX_ROpositive_el1_N01_part.txt',
 'Logic -X' : 'From20150504053437_To20150504062550_NuclearDD_111_1_sil18_sweep_evolution_time_auto_C2&6_LogicmX_TomoXX_ROpositive_el1_N01_part.txt',}
# Elec_Norbert={
# '128'   : 'FirstTS_20150429_110632_LastTS_20150430_104617_Electron_DD_N128_XY8_Pts72_Reps250.txt'}
def un2str(x, xe, precision=1):
    """pretty print nominal value and uncertainty

    x  - nominal value
    xe - uncertainty
    precision - number of significant digits in uncertainty

    returns shortest string representation of `x +- xe` either as
        x.xx(ee)e+xx
    or as
        xxx.xx(ee)"""
    # base 10 exponents
    x_exp = int(floor(log10(x)))
    xe_exp = int(floor(log10(xe)))

    # uncertainty
    un_exp = xe_exp-precision+1
    un_int = round(xe*10**(-un_exp))

    # nominal value
    no_exp = un_exp
    no_int = round(x*10**(-no_exp))

    # format - nom(unc)exp
    fieldw = x_exp - no_exp
    fmt = '%%.%df' % fieldw
    result1 = (fmt + '(%.0f)e%d') % (no_int*10**(-fieldw), un_int, x_exp)

    # format - nom(unc)
    fieldw = max(0, -no_exp)
    fmt = '%%.%df' % fieldw
    result2 = (fmt + '(%.0f)') % (no_int*10**no_exp, un_int*10**max(0, un_exp))

    # return shortest representation
    if len(result2) <= len(result1):
        return result2
    else:
        return result1

def fit_general_exponential(g_a, g_A, g_x0, g_T, g_n):
    fitfunc_str = 'a + A * exp(-((x-x0)/T )**n)'

    a = fit.Parameter(g_a, 'a')
    A = fit.Parameter(g_A, 'A')
    x0 = fit.Parameter(g_x0, 'x0')
    T = fit.Parameter(g_T, 'T')
    n = fit.Parameter(g_n, 'n')

    p0 = [a, A, x0, T, n]

    def fitfunc(x):
        return a() + A() * np.exp(-(x-x0())**n()/(T()**n()))
    return p0, fitfunc, fitfunc_str

def DFS_comparison(msmts_dict,
    offset = 0.5, 
    x0 = 0,  
    amplitude = 0.4,  
    decay_constant = 0.426546, 
    exponent = 1.75, fixed = [0,2,4]):

    plt.close('all')

    # msmts_list = ['pX,mX,XX']
    
    color = ['r','b','g','m','k']
    x = []
    y = []
    y_err = []
    fig = plt.figure(figsize=(7,5))
    # mpl.rcParams['axes.linewidth'] = 2
    # max_x = 8
    max_x = 1.2
    min_x = 0.
    # fig3 = plt.figure(figsize=(10,8))
    Amps = []
    Amps_err = []
    Ts = []
    Ts_err = []
    exps = []
    exps_err = []
    ax = fig.add_subplot(111)
    # ax3 = fig3.add_subplot(111)
    # p0 = [amplitude, decay_constant, exponent]
    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
             x0, decay_constant,exponent)
    for aaa, dictms in enumerate(msmts_dict):
        #msmts_list = list(dictms.keys())
        msmts_list =[r'$|00\rangle+|11\rangle$' , r'$|01\rangle+|10\rangle$']
        # msmts_list =[r'$|00\rangle+|11\rangle$']
        for ii,msmt in enumerate(msmts_list):
            array = np.loadtxt(DBdir+dictms[msmt],skiprows=1)
            x.append(array[:,0])
            y.append(array[:,1])
            y_err.append(array[:,2])
            # print x
            # print y
            # print y_err
            # optimize.leastsq(fit_func,x,y,p0)
            p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
                 x0, decay_constant,exponent)
            fit_result = fit.fit1d(x[-1],y[-1], None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
            
            Ts.append(fit_result['params_dict']['T'])
            Ts_err.append(fit_result['error_dict']['T'])
            Amps.append(fit_result['params_dict']['A'])
            Amps_err.append(fit_result['error_dict']['A'])
            if 4 in fixed:
                # ax.errorbar(x,y,fmt='o',yerr=y_err, label='N=' + str(N)+' T=' + '%.2f' % fit_result['params_dict']['T'] + '+-' + '%.2f' % fit_result['error_dict']['T']
                # +', A=' '%.2f' % fit_result['params_dict']['A'] + '+-' + '%.2f' % fit_result['error_dict']['A'],color=color[ii])
                ax.errorbar(x[-1],y[-1],fmt='o',yerr=y_err[-1], label=msmt,color=color[ii])
            else:
                exps.append(fit_result['params_dict']['n'])
                exps_err.append(fit_result['error_dict']['n'])
                ax.errorbar(x[-1],y[-1],fmt='o',yerr=y_err[-1], label=msmt,color=color[ii])
                # ax.errorbar(x[-1],y[-1],fmt='o',yerr=y_err[-1], label=msmt + ',  T=' + '%.2f' % fit_result['params_dict']['T'] + '+-' + '%.2f' % fit_result['error_dict']['T']
                # +', A=' '%.2f' % fit_result['params_dict']['A'] + '+-' + '%.2f' % fit_result['error_dict']['A'] +', n=' '%.2f' % fit_result['params_dict']['n'] + '+-' + '%.2f' % fit_result['error_dict']['n'],color=color[ii])
            if aaa == 0:
                ls = '-'
            else:
                ls = '--'
            plot.plot_fit1d(fit_result, np.linspace(-0.1,max_x,1001), ax=ax, plot_data=False,print_info=False,linestyle=ls,color=color[ii])
            
    # ax.set_xscale('log')
    ax.hlines([0.5],min_x,max_x,linestyles='dotted', linewidth = 2)
    ax.hlines([0.],min_x,max_x,linestyles='dotted', linewidth = 2)
    # plt.text(1, 0.8, r'$\left< XX \right>$', fontsize=30)
    ax.hlines([1.],min_x,max_x,linestyles='dotted', linewidth = 2)
    ax.set_xlim(min_x,max_x)
    ax.set_ylim(0.45,0.85)
    plt.xticks([0.0,0.5,1.0])
    plt.yticks([0.5,0.6,0.7,0.8],['0.0','0.2','0.4','0.6'])
    ax.set_xlabel('Free evolution time (s)',fontsize = 15)
    ax.set_ylabel(r'$\langle$XX$\rangle$',fontsize = 20)
    ax.tick_params(axis='x', which='major', labelsize=15)
    ax.tick_params(axis='y', which='major', labelsize=15)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2.)
    plt.legend(loc = 'upper right',fontsize = 20,numpoints = 1,frameon = False,columnspacing=0.5,handletextpad=0.0)

    # plt.savefig(r'/Users/'+os.getlogin()+r'/Dropbox/DFS36XX.pdf', bbox_inches='tight')
    #mpl.rcParams['axes.linewidth'] = 2
    plt.show()

DFS_comparison([C36_DFS])



# DD_scaling_mult_exp(C1_oldXY_msmts)
