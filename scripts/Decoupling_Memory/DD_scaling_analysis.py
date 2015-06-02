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
        DBdir = r'/Users/'+os.getlogin()+r'/Dropbox/QEC LT/Decoupling memory/Electron_DD_Data/'
else:
    DBdir = r'D:/Dropbox/QEC LT/Decoupling memory/Electron_DD_Data/'

reload(common)
reload(mbi)
from scipy import optimize

C2_XYauto_msmts={
'1' : '20150414191406_C2_N1_posneg_Pts8_Reps280.txt',
'4' : '20150414204529_C2_N4_posneg_Pts8_Reps280.txt',
'8' : '20150414230710_C2_N8_posneg_Pts8_Reps280.txt',
'16': '20150415025358_C2_N16_posneg_Pts8_Reps280.txt',
'32': '20150415074830_C2_N32_posneg_Pts8_Reps280.txt'}

C1_RFamp_msmts={
'1' : '20150415222144_C1_N1_posneg_Pts8_Reps220.txt',
'4' : '20150415233508_C1_N4_posneg_Pts8_Reps220.txt',
'8' : '20150416012948_C1_N8_posneg_Pts8_Reps220.txt',
'16': '20150416043428_C1_N16_posneg_Pts8_Reps220.txt',
'32': '20150416043428_C1_N16_posneg_Pts8_Reps220.txt'}

C1_oldXY_msmts={
'1' : '20150410234153_C1_N1_positive_Pts8_Reps300.txt',
'4' : '20150411004456_C1_N4_positive_Pts8_Reps300.txt',
'8' : '20150411023445_C1_N8_posneg_Pts8_Reps300.txt',
'16': '20150411054458_C1_N16_posneg_Pts8_Reps300.txt',
'32': '20150411104646_C1_N32_posneg_Pts8_Reps300.txt'}

C1_XYauto_msmts={
'1' : '20150416214508_C1_N1_posneg_Pts8_Reps250.txt',
'4' : '20150416230733_C1_N4_posneg_Pts8_Reps250.txt',
'8' : '20150417011548_C1_N8_posneg_Pts8_Reps250.txt',
'16': '20150417044207_C1_N16_posneg_Pts8_Reps250.txt',
'32': '20150417090916_C1_N32_posneg_Pts8_Reps250.txt'}

Elec_Norbert={
'128'   : 'FirstTS_20150429_110632_LastTS_20150430_104617_Electron_DD_N128_XY8_Pts72_Reps250.txt',
'256'   : 'FirstTS_20150429_113126_LastTS_20150430_100149_Electron_DD_N256_XY8_Pts55_Reps250.txt',
'512'   : 'FirstTS_20150429_115420_LastTS_20150430_101511_Electron_DD_N512_XY8_Pts45_Reps250.txt',
'1024'  : 'FirstTS_20150429_120802_LastTS_20150430_103340_Electron_DD_N1024_XY8_Pts35_Reps250.txt',
'2048'  : 'FirstTS_20150430113907_LastTS_20150430104832_Electron_DD_N2048_XY8_Pts30_Reps250.txt'}

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

def DD_scaling_elec(msmts,
    offset = 0.5, 
    x0 = 0,  
    dip_std = 1.,
    amplitude = 0.5,  
    decay_constant = 70, 
    correct_dips = True,
    exponent = 3., fixed = [0,2,4]):

    plt.close('all')

    Nlist = [1,4,8,16,32]
    Nlist = [128,256,512,1024,2048]
    # Nlist = [2048]
    color = ['r','g','b','m','k']
    
    
    fig = plt.figure(figsize=(10,8))
    mpl.rcParams['axes.linewidth'] = 2
    max_x = 8
    max_x = 1000
    min_x = 1
    # fig3 = plt.figure(figsize=(10,8))
    Amps_Tot = []
    Amps_err_Tot = []
    Ts_Tot = []
    Ts_err_Tot = []
    exps_Tot = []
    exps_err_Tot = []


    ax = fig.add_subplot('111')
    # ax3 = fig3.add_subplot(111)
    # p0 = [amplitude, decay_constant, exponent]
    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
             x0, decay_constant,exponent)
    x = []
    y = []
    y_err = []
    exponents = [3.]

    for ii,N in enumerate(Nlist):
        
        Amps = []
        Amps_err = []
        Ts = []
        Ts_err = []
        exps = []
        exps_err = []
        
        array = np.loadtxt(DBdir+msmts[str(N)],skiprows=1)
        
        
        sorting_order = array[:,0].argsort()
        array[:,0].sort()
        array[:,1]=array[:,1][sorting_order]
        array[:,2]=array[:,2][sorting_order]

        print array[:,1]
        print array[:,0]
        if correct_dips:
            x_old = array[:,0]
            y_old = array[:,1]
            y_err_old = array[:,2]
            x_new=[]
            y_new=[]
            y_err_new=[]
            for jj in range(len(x_old)-1):
                if all([ (y_old[jj] + dip_std*y_err_old[jj]) > value for value in y_old[jj+1::]]):
                    x_new.append(x_old[jj])
                    y_new.append(y_old[jj])
                    y_err_new.append(y_err_old[jj])
                else:
                    print y_old[jj], y_old[jj+1]
            print len(x_new), len(y_new), len(y_err_new)
            print x_new
            print y_new
            print y_err_new
            x.append(np.array(x_new))
            y.append(np.array(y_new))
            y_err.append(np.array(y_err_new))
        else:
            x.append(array[:,0])
            y.append(array[:,1])
            y_err.append(array[:,2])
        # if 'correct'

        for exp_nr, exponent in enumerate(exponents):
            print 'NUMBER', exp_nr
            if exp_nr == len(exponents)-1 and False:
                fixed = [0, 2]
            else:
                fixed = [0, 1, 2,4]
            
        
        # print x
        # print y
        # print y_err
        # optimize.leastsq(fit_func,x,y,p0)
            p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
                 x0, decay_constant,exponent)
            fit_result = fit.fit1d(x[-1],y[-1], None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
            
            Ts.append(fit_result['params_dict']['T'])
            Ts_err.append(fit_result['error_dict']['T'])
            
            # Amps.append(fit_result['params_dict']['A'])
            # Amps_err.append(fit_result['error_dict']['A'])
            

            if 4 in fixed:
                # ax.errorbar(x,y,fmt='o',yerr=y_err, label='N=' + str(N)+' T=' + '%.2f' % fit_result['params_dict']['T'] + '+-' + '%.2f' % fit_result['error_dict']['T']
                # +', A=' '%.2f' % fit_result['params_dict']['A'] + '+-' + '%.2f' % fit_result['error_dict']['A'],color=color[ii])
                if N > 1000:
                    ax.errorbar(x[-1],y[-1],fmt='o',yerr=y_err[-1], label='N=' + str(N)+ ', T=' + un2str(fit_result['params_dict']['T'], fit_result['error_dict']['T']) ,color=color[ii])
                    # ax.errorbar(x[-1],y[-1],fmt='o',yerr=y_err[-1], label='N=' + str(N)+ ', T=' + un2str(fit_result['params_dict']['T'], fit_result['error_dict']['T']) + ', A=' + un2str(fit_result['params_dict']['A'], fit_result['error_dict']['A']),color=color[ii])
                else:
                    # ax.errorbar(x[-1],y[-1],fmt='o',yerr=y_err[-1], label='N=' + str(N)+ ',   T=' + un2str(fit_result['params_dict']['T'], fit_result['error_dict']['T']) + ', A=' + un2str(fit_result['params_dict']['A'], fit_result['error_dict']['A']),color=color[ii])
                    ax.errorbar(x[-1],y[-1],fmt='o',yerr=y_err[-1], label='N=' + str(N)+ ',   T=' + un2str(fit_result['params_dict']['T'], fit_result['error_dict']['T']),color=color[ii])
                print 'HELLO MOTO'
            else:
                exps.append(fit_result['params_dict']['n'])
                exps_err.append(fit_result['error_dict']['n'])
                # ax.errorbar(x[-1],y[-1],fmt='o',yerr=y_err[-1], label='N=' + str(N)+', T=' + un2str(fit_result['params_dict']['T'], fit_result['error_dict']['T'])
                # + ', A=' + un2str(fit_result['params_dict']['A'], fit_result['error_dict']['A']) +', n=' + un2str(fit_result['params_dict']['n'], fit_result['error_dict']['n']),color=color[ii])
                ax.errorbar(x[-1],y[-1],fmt='o',yerr=y_err[-1], label='N=' + str(N)+', T=' + un2str(fit_result['params_dict']['T'], fit_result['error_dict']['T'])
                 +', n=' + un2str(fit_result['params_dict']['n'], fit_result['error_dict']['n']),color=color[ii])
            plot.plot_fit1d(fit_result, np.linspace(min_x,max_x,1001), ax=ax, plot_data=False,print_info=False,color=color[ii])
        
        Amps_Tot.append(Amps)
        Amps_err_Tot.append(Amps_err)
        Ts_Tot.append(Ts)
        Ts_err_Tot.append(Ts_err)
        exps_Tot.append(exps)
        exps_err_Tot.append(exps_err)

    print Amps_Tot
    Amps_Tot = list(map(list, zip(*Amps_Tot)))
    Amps_err_Tot = list(map(list, zip(*Amps_err_Tot)))
    Ts_Tot = list(map(list, zip(*Ts_Tot)))
    Ts_err_Tot = list(map(list, zip(*Ts_err_Tot)))
    exps_Tot = list(map(list, zip(*exps_Tot)))
    exps_err_Tot = list(map(list, zip(*exps_err_Tot)))
    ax.set_xscale('log')
    ax.hlines([0.5],min_x,max_x,linestyles='dotted', linewidth = 2)

    ax.hlines([1.],min_x,max_x,linestyles='dotted', linewidth = 2)
    ax.set_xlim(min_x,max_x)
    ax.set_xlabel('Free evolution time (ms)',fontsize = 25)
    ax.set_ylabel('Fidelity',fontsize = 20)
    ax.tick_params(axis='x', which='major', labelsize=20)
    ax.tick_params(axis='y', which='major', labelsize=20)
    mpl.rcParams['axes.linewidth'] = 2
    plt.legend(loc = 'lower left',fontsize = 16)

    # plt.show()
    logx = np.log10(Nlist)
    def fit_func(x, a, b):
            return a*x + b

    def fit_func(x, b):
            return (2./3.)*x + b

    for exp_nr, exponent in enumerate(exponents):
        if exp_nr == 0:
            fig3 = plt.figure(figsize=(10,8))
            ax3 = fig3.add_subplot(111)

            
    
        print Ts_Tot[exp_nr]
        # sdfklaj
        logy = np.log10(Ts_Tot[exp_nr])

        print 'logx', logx
        print 'logy', logy

        params, covs = optimize.curve_fit(fit_func, logx, logy)

        # print 'Fitting gives'
        # print 'gamma =', params[0]
        # print 'T_2 =  ', 10**(params[1])
        # gamma = params[0]
        # T_20 = 10**(params[1])
        # gamma_err = np.sqrt(covs[0][0])
        # T_20_err = np.sqrt(covs[1][1])*T_20

        print 'Fitting gives'
        print 'T_2 =  ', 10**(params[0])
        T_20 = 10**(params[0])
        T_20_err = np.sqrt(covs[0][0])*T_20

        # print logx
        # print logy
        print covs
        print 
        print params

        print covs[0][0]
        # print 10**(covs[0][1])
        # for a in range(len(Ts)):
        #     yerrp.append( abs(np.log(Ts[a]+Ts_err[a])-np.log(Ts[a])))
        #     yerrm.append( abs(np.log(Ts[a]-Ts_err[a])-np.log(Ts[a])))
        if exp_nr == 0:
            fig2 = plt.figure(figsize=(8,6))
            ax2 = fig2.add_subplot(111)
        ax2.errorbar(Nlist,Ts_Tot[exp_nr], yerr=Ts_err_Tot[exp_nr], fmt='o',color=color[exp_nr])
        if exp_nr == 0:
            ax2.set_xscale('log')
            ax2.set_yscale('log')
            ax2.set_xlabel('Number of Pulses',fontsize = 20)
            ax2.set_ylabel('Coherence time (ms)',fontsize = 20)

        xforplot = np.linspace(1,ax2.get_xlim()[1],501)
        if exp_nr == len(exponents)-1 and False:
            ax2.plot(xforplot,10**(params[1])*xforplot**(params[0]),color=color[exp_nr],
                label='exp= Free, gamma= %.3f +- %.3f, T_20 = %.2f +- %.2f' % (gamma, gamma_err, T_20, T_20_err))
            
            # ax2.plot(xforplot,10**(params[1])*xforplot**(params[0]), color =color[exp_nr] ,
            #     label= 'exp=' + '%.2f' % exponent + ', gamma=' + '%.3f' % params[0] + '+-' + '%.4f' % covs[0][0])
        else:
            # ax2.plot(xforplot,10**(params[1])*xforplot**(params[0]),color=color[exp_nr],
            #     label='exp= %.2f, gamma= %.3f +- %.3f, T_20 = %.2f +- %.2f' % (exponent,gamma, gamma_err, T_20, T_20_err))
            ax2.plot(xforplot,10**(params[0])*xforplot**(2./3.),color=color[exp_nr],
                label='exp= %.2f, gamma=2/3, T_20 = %.2f +- %.2f' % (exponent, T_20, T_20_err))
        # 'T_2=' + '%.2f' % (10**(params[0][1])) + '+-' + '%.2f' % (10**(params[1][1])) 
    plt.legend(loc='lower left', prop={'size':13})
    ax2.set_xlim(10.**2,10.**4)
    plt.show()
DD_scaling_elec(Elec_Norbert)

if os.name == 'posix':
        DBdir = r'/Users/'+os.getlogin()+r'/Dropbox/QEC LT/Decoupling memory/XYdata/'
else:
    DBdir = r'D:/Dropbox/QEC LT/Decoupling memory/XYdata/'

def DD_scaling(msmts,
    offset = 0.5, 
    x0 = 0,  
    amplitude = 0.2,  
    decay_constant = 1.2, 
    exponent = 1.7, fixed = [0,2,4]):

    plt.close('all')

    Nlist = [1,4,8,16,32]
    # Nlist = [128,256,512,1024,2048]
    color = ['r','g','b','m','k']
    x = np.zeros((8,len(Nlist)))
    y = np.zeros((8,len(Nlist)))
    y_err = np.zeros((8,len(Nlist)))
    fig = plt.figure(figsize=(10,8))
    mpl.rcParams['axes.linewidth'] = 2

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
    for ii,N in enumerate(Nlist):
        array = np.loadtxt(DBdir+msmts[str(N)],skiprows=1)
        x[:,ii]=array[:,0]
        y[:,ii]=array[:,1]
        y_err[:,ii] = array[:,2]

        # optimize.leastsq(fit_func,x,y,p0)
        p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
             x0, decay_constant,exponent)
        fit_result = fit.fit1d(x[:,ii],y[:,ii], None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
        
        Ts.append(fit_result['params_dict']['T'])
        Ts_err.append(fit_result['error_dict']['T'])
        Amps.append(fit_result['params_dict']['A'])
        Amps_err.append(fit_result['error_dict']['A'])
        if 4 in fixed:
            ax.errorbar(x[:,ii],y[:,ii],fmt='o',yerr=y_err[:,ii], label='N=' + str(N)+' T=' + '%.2f' % fit_result['params_dict']['T'] + '+-' + '%.2f' % fit_result['error_dict']['T']
            +', A=' '%.2f' % fit_result['params_dict']['A'] + '+-' + '%.2f' % fit_result['error_dict']['A'],color=color[ii])
        else:
            exps.append(fit_result['params_dict']['n'])
            exps_err.append(fit_result['error_dict']['n'])
            ax.errorbar(x[:,ii],y[:,ii],fmt='o',yerr=y_err[:,ii], label='N=' + str(N)+' T=' + '%.2f' % fit_result['params_dict']['T'] + '+-' + '%.2f' % fit_result['error_dict']['T']
            +', A=' '%.2f' % fit_result['params_dict']['A'] + '+-' + '%.2f' % fit_result['error_dict']['A'] +', n=' '%.2f' % fit_result['params_dict']['n'] + '+-' + '%.2f' % fit_result['error_dict']['n'],color=color[ii])
        plot.plot_fit1d(fit_result, np.linspace(0,8,1001), ax=ax, plot_data=False,print_info=False,color=color[ii])
        
    
    ax.hlines([0.5],0,8,linestyles='dotted', linewidth = 2)
    ax.set_xlim(0,8)
    ax.set_xlabel('Free evolution time (s)',fontsize = 25)
    ax.set_ylabel('Fidelity',fontsize = 20)
    ax.tick_params(axis='x', which='major', labelsize=20)
    ax.tick_params(axis='y', which='major', labelsize=20)
    mpl.rcParams['axes.linewidth'] = 2
    plt.legend(prop={'size':15})


    fig3 = plt.figure(figsize=(10,8))
    ax3 = fig3.add_subplot(111)
    y_norm = np.zeros(np.shape(y))
    y_err_norm = np.zeros(np.shape(y))
    xforplot = np.linspace(0,8,1001)
    for ii,N in enumerate(Nlist):
        y_norm[:,ii] = (y[:,ii]-0.5)/Amps[ii]
        y_err_norm[:,ii] = (y_err[:,ii])/Amps[ii]
        # y_norm[:,ii] = y[:,ii]
        # y_err_norm[:,ii] = y_err[:,ii]
        if 4 in fixed:
            # ax3.errorbar(x[:,ii],y_norm[:,ii],fmt='o',yerr=y_err_norm[:,ii], label='N=' + str(N)+' T=' + '%.2f' % Ts[ii] + '+-' + '%.2f' % Ts_err[ii]
            #     +', A=' '%.2f' % Amps[ii] + '+-' + '%.2f' % Amps_err[ii],color=color[ii])
            if N>9:
                ax3.errorbar(x[:,ii],y_norm[:,ii],fmt='o',yerr=y_err_norm[:,ii], label='N=' + str(N)+', T=' + un2str(Ts[ii], Ts_err[ii]), color=color[ii])
                ax3.plot(xforplot, np.exp(-xforplot**exponent/(Ts[ii]**exponent)),color=color[ii])
            else:
                ax3.errorbar(x[:,ii],y_norm[:,ii],fmt='o',yerr=y_err_norm[:,ii], label='N=' + str(N)+',   T=' + un2str(Ts[ii], Ts_err[ii]), color=color[ii])
                ax3.plot(xforplot, np.exp(-xforplot**exponent/(Ts[ii]**exponent)),color=color[ii]) 
            # ax3.plot(xforplot, 0.5+Amps[ii]*np.exp(-xforplot**exponent/(Ts[ii]**exponent)),color=color[ii])

        else:
            ax3.errorbar(x[:,ii],y_norm[:,ii],fmt='o',yerr=y_err_norm[:,ii], label='N=' + str(N)+' T=' + '%.2f' % Ts[ii] + '+-' + '%.2f' % Ts_err[ii]
                +', A=' '%.2f' % Amps[ii] + '+-' + '%.2f' % Amps_err[ii] + ', exp=' '%.2f' % exps[ii] + '+-' + '%.2f' % exps_err[ii],color=color[ii])
            ax3.plot(xforplot, np.exp(-xforplot**exps[ii]/(Ts[ii]**exps[ii])),color=color[ii])
            # ax3.plot(xforplot, 0.5+Amps[ii]*np.exp(-xforplot**exps[ii]/(Ts[ii]**exps[ii])),color=color[ii])

    logx = np.log10(Nlist)
    print 'y', Ts
    logy = np.log10(Ts)
    print 'logx', logx
    print 'logy', logy

    plt.legend(loc = 'center left', fontsize=18,frameon=False)
    ax3.set_xscale('log')
    ax3.set_xlim(np.min(x),8)
    ax3.hlines([1.],np.min(x),8,linestyles='dotted', linewidth = 2)
    ax3.hlines([0.],np.min(x),8,linestyles='dotted', linewidth = 2)

    ax3.set_xlabel('Free evolution time (s)',fontsize = 22)
    ax3.set_ylabel('Normalized Signal',fontsize = 22)
    ax3.tick_params(axis='x', which='major', labelsize=22)
    ax3.tick_params(axis='y', which='major', labelsize=22)

    # ax3.legend(loc='lower left')

    def fit_func(x, a, b):
        return a*x + b
    params, covs = optimize.curve_fit(fit_func, logx, logy)
    print 'Fitting gives'
    print 'gamma =', params[0]
    print 'T_2 =  ', 10**(params[1])
    gamma = params[0]
    T_20 = 10**(params[1])
    gamma_err = np.sqrt(covs[0][0])
    T_20_err = np.sqrt(covs[1][1])*T_20
    # print logx
    # print logy
    print covs
    print 
    print params

    print covs[0][0]
    # print 10**(covs[0][1])
    # for a in range(len(Ts)):
    #     yerrp.append( abs(np.log(Ts[a]+Ts_err[a])-np.log(Ts[a])))
    #     yerrm.append( abs(np.log(Ts[a]-Ts_err[a])-np.log(Ts[a])))
    
    fig2 = plt.figure(figsize=(8,6))
    ax2 = fig2.add_subplot(111)
    ax2.errorbar(Nlist,Ts, yerr=Ts_err, fmt='o',color='r')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Number of Pulses',fontsize = 20)
    ax2.set_ylabel('Coherence time (s)',fontsize = 20)

    xforplot = np.linspace(1,ax2.get_xlim()[1],501)
    #http://wiki.scipy.org/Cookbook/FittingData#head-5eba0779a34c07f5a596bbcf99dbc7886eac18e5
    ax2.plot(xforplot,10**(params[1])*xforplot**(params[0]),color='r',
        label='gamma= %.2f +- %.4f, T_20 = %.2f +- %.4f' % (gamma, gamma_err, T_20, T_20_err))
    # 'T_2=' + '%.2f' % (10**(params[0][1])) + '+-' + '%.2f' % (10**(params[1][1])) 
    ax2.legend()
    info_x = ax2.get_xlim()[0] + (ax2.get_xlim()[-1]-ax2.get_xlim()[0])*0.02
    info_y = ax2.get_ylim()[0] + (ax2.get_ylim()[-1]-ax2.get_ylim()[0])*0.02
    # ax2.text(info_x, info_y, '$T_2$'=, size='x-small',
    #                 color='k', ha='left', va='bottom',
    #                 bbox=dict(facecolor='white', alpha=0.5))
    
    plt.show()

# DD_scaling(C2_XYauto_msmts)

def DD_scaling_mult_exp(msmtss,
    offset = 0.5, 
    x0 = 0,  
    amplitude = 0.4,  
    decay_constant = 1.5, 
    exponent = 1.75, fixed = [0,2, 4]):
    plt.close('all')

    Nlist = [1,4,8,16,32]
    color = ['r','g','b','m','k']
    x = np.zeros((8,len(Nlist)))
    y = np.zeros((8,len(Nlist)))
    y_err = np.zeros((8,len(Nlist)))
    fig = plt.figure(figsize=(10,8))
    # fig3 = plt.figure(figsize=(10,8))
    
    ax = fig.add_subplot(111)
    exponents = [1.25, 1.5, 1.75, 1.5]
    # ax3 = fig3.add_subplot(111)
    # p0 = [amplitude, decay_constant, exponent]
    
    for iiii, msmts in enumerate(msmtss):
        x = np.zeros((8,len(Nlist)))
        y = np.zeros((8,len(Nlist)))
        y_err = np.zeros((8,len(Nlist)))
        exponents = [1.25, 1.5, 1.75, 1.5]
        for exp_nr, exponent in enumerate(exponents):
            print 'NUMBER', exp_nr
            if exp_nr == len(exponents)-1:
                fixed = [0, 2]
            else:
                fixed = [0, 2, 4]
            Amps = []
            Amps_err = []
            Ts = []
            Ts_err = []
            exps = []
            exps_err = []

            for ii,N in enumerate(Nlist):
                # print 'jaaaaaa', msmts
                array = np.loadtxt(DBdir+msmts[str(N)],skiprows=1)
                x[:,ii]=array[:,0]
                y[:,ii]=array[:,1]
                y_err[:,ii] = array[:,2]

                # optimize.leastsq(fit_func,x[:,ii],y[:,ii],p0)
                p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
                     x0, decay_constant,exponent)
                fit_result = fit.fit1d(x[:,ii],y[:,ii], None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
                
                Ts.append(fit_result['params_dict']['T'])
                Ts_err.append(fit_result['error_dict']['T'])
                Amps.append(fit_result['params_dict']['A'])
                Amps_err.append(fit_result['error_dict']['A'])
                if 4 in fixed:
                    ax.errorbar(x[:,ii],y[:,ii],fmt='o',yerr=y_err[:,ii], label='N=' + str(N)+' T=' + '%.2f' % fit_result['params_dict']['T'] + '+-' + '%.2f' % fit_result['error_dict']['T']
                    +', A=' '%.2f' % fit_result['params_dict']['A'] + '+-' + '%.2f' % fit_result['error_dict']['A'],color=color[ii])
                else:
                    exps.append(fit_result['params_dict']['n'])
                    exps_err.append(fit_result['error_dict']['n'])
                    ax.errorbar(x[:,ii],y[:,ii],fmt='o',yerr=y_err[:,ii], label='N=' + str(N)+' T=' + '%.2f' % fit_result['params_dict']['T'] + '+-' + '%.2f' % fit_result['error_dict']['T']
                    +', A=' '%.2f' % fit_result['params_dict']['A'] + '+-' + '%.2f' % fit_result['error_dict']['A'] +', n=' '%.2f' % fit_result['params_dict']['n'] + '+-' + '%.2f' % fit_result['error_dict']['n'],color=color[ii])
                plot.plot_fit1d(fit_result, np.linspace(0,8,1001), ax=ax, plot_data=False,print_info=False,color=color[ii])
                
            if exp_nr == len(exponents)-1 and iiii==0:   
                ax.set_xlim(0,8)
                ax.set_xlabel('Free evolution time (s)',fontsize = 20)
                ax.set_ylabel('Fidelity',fontsize = 20)
                plt.legend()

            if exp_nr == 0 and iiii==0:
                fig3 = plt.figure(figsize=(10,8))
                ax3 = fig3.add_subplot(111)
            y_norm = np.zeros(np.shape(y))
            y_err_norm = np.zeros(np.shape(y))
            xforplot = np.linspace(0,8,1001)
            for ii,N in enumerate(Nlist):
                y_norm[:,ii] = (y[:,ii]-0.5)/Amps[ii]
                y_err_norm[:,ii] = (y_err[:,ii])/Amps[ii]
                if 4 in fixed:
                    ax3.errorbar(x[:,ii],y_norm[:,ii],fmt='o',yerr=y_err_norm[:,ii], label='N=' + str(N)+' T=' + '%.2f' % Ts[ii] + '+-' + '%.2f' % Ts_err[ii]
                        +', A=' '%.2f' % Amps[ii] + '+-' + '%.2f' % Amps_err[ii],color=color[ii])
                    ax3.plot(xforplot, np.exp(-xforplot**exponent/(Ts[ii]**exponent)),color=color[ii])
                else:
                    ax3.errorbar(x[:,ii],y_norm[:,ii],fmt='o',yerr=y_err_norm[:,ii], label='N=' + str(N)+' T=' + '%.2f' % Ts[ii] + '+-' + '%.2f' % Ts_err[ii]
                        +', A=' '%.2f' % Amps[ii] + '+-' + '%.2f' % Amps_err[ii] + ', exp=' '%.2f' % exps[ii] + '+-' + '%.2f' % exps_err[ii],color=color[ii])
                    ax3.plot(xforplot, np.exp(-xforplot**exps[ii]/(Ts[ii]**exps[ii])),color=color[ii])

            logx = np.log10(Nlist)
            print 'y', Ts
            logy = np.log10(Ts)
            print 'logx', logx
            print 'logy', logy
            if exp_nr == 0 and iiii==0:
                ax3.set_xscale('log')
                ax3.set_xlim(np.min(x),8)
                ax3.set_xlabel('Free evolution time (s)',fontsize = 20)
                ax3.set_ylabel('Normalized signal',fontsize = 20)
                ax3.legend(loc='lower left')

            def fit_func(x, a, b):
                return a*x + b
            params, covs = optimize.curve_fit(fit_func, logx, logy)
            print 'Fitting gives'
            print 'gamma =', params[0]
            print 'T_2 =  ', 10**(params[1])
            gamma = params[0]
            T_20 = 10**(params[1])
            gamma_err = np.sqrt(covs[0][0])
            T_20_err = np.sqrt(covs[1][1])*T_20
            # print logx
            # print logy
            print covs
            print 
            print params

            print covs[0][0]
            # print 10**(covs[0][1])
            # for a in range(len(Ts)):
            #     yerrp.append( abs(np.log(Ts[a]+Ts_err[a])-np.log(Ts[a])))
            #     yerrm.append( abs(np.log(Ts[a]-Ts_err[a])-np.log(Ts[a])))
            if exp_nr == 0 and iiii==0:
                fig2 = plt.figure(figsize=(8,6))
                ax2 = fig2.add_subplot(111)
            ax2.errorbar(Nlist,Ts, yerr=Ts_err, fmt='o',color=color[exp_nr])
            if exp_nr == 0 and iiii==0:
                ax2.set_xscale('log')
                ax2.set_yscale('log')
                ax2.set_xlabel('Number of Pulses',fontsize = 20)
                ax2.set_ylabel('Coherence time (s)',fontsize = 20)

            xforplot = np.linspace(1,ax2.get_xlim()[1],501)
            if 4 in fixed:
                ax2.plot(xforplot,10**(params[1])*xforplot**(params[0]),color=color[exp_nr],
                    label='exp= %.2f, gamma= %.3f +- %.3f, T_20 = %.2f +- %.2f' % (exponent, gamma, gamma_err, T_20, T_20_err))
                # ax2.plot(xforplot,10**(params[1])*xforplot**(params[0]), color =color[exp_nr] ,
                #     label= 'exp=' + '%.2f' % exponent + ', gamma=' + '%.3f' % params[0] + '+-' + '%.4f' % covs[0][0])
            else:
                ax2.plot(xforplot,10**(params[1])*xforplot**(params[0]),color=color[exp_nr],
                    label='exp= Free, gamma= %.3f +- %.3f, T_20 = %.2f +- %.2f' % (gamma, gamma_err, T_20, T_20_err))
                # ax2.plot(xforplot,10**(params[1])*xforplot**(params[0]),color = color[exp_nr]  ,
                #     label= 'exp= free, gamma=' + '%.3f' % params[0] + '+-' + '%.4f' % covs[0][0])
            # 'T_2=' + '%.2f' % (10**(params[0][1])) + '+-' + '%.2f' % (10**(params[1][1])) 
        


    info_x = ax2.get_xlim()[0] + (ax2.get_xlim()[-1]-ax2.get_xlim()[0])*0.02
    info_y = ax2.get_ylim()[0] + (ax2.get_ylim()[-1]-ax2.get_ylim()[0])*0.02
    # ax2.text(info_x, info_y, '$T_2$'=, size='x-small',
    #                 color='k', ha='left', va='bottom',
    #                 bbox=dict(facecolor='white', alpha=0.5))
    ax2.legend(loc='lower left', prop={'size':13})
    plt.show()
# DD_scaling(C1_XYauto_msmts)

# DD_scaling_mult_exp([C2_XYauto_msmts,C1_XYauto_msmts])

# DD_scaling_mult_exp(C1_oldXY_msmts)
