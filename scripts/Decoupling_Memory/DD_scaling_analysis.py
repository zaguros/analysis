'''Written by MAB 10-3-15 for a general coherence msmt with a "free" exponent'''

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt

reload(common)
reload(mbi)
from scipy import optimize

C2_msmts={
'1' : '20150414191406_C2_N1_posneg_Pts8_Reps280.txt',
'4' : '20150414204529_C2_N4_posneg_Pts8_Reps280.txt',
'8' : '20150414230710_C2_N8_posneg_Pts8_Reps280.txt',
'16': '20150415025358_C2_N16_posneg_Pts8_Reps280.txt',
'32': '20150415074830_C2_N32_posneg_Pts8_Reps280.txt'}

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
        return a + A * np.exp(-(x-x0)**n/(T**n))
    return p0, fitfunc, fitfunc_str




def DD_scaling(msmts,
    offset = 0.5, 
    x0 = 0,  
    amplitude = 0.4,  
    decay_constant = 2, 
    exponent = 1.5, fixed = [0,2,4]):
    Head_folder = 'D:/Dropbox/QEC LT/Decoupling memory/XYdata/'

    Nlist = [1,4,8,16,32]
    color = ['r','g','b','m','k']
    x = np.zeros((8,len(Nlist)))
    y = np.zeros((8,len(Nlist)))
    y_err = np.zeros((8,len(Nlist)))
    fig = plt.figure(figsize=(10,8))
    Amps = []
    Amps_err = []
    Ts = []
    Ts_err = []
    ax = fig.add_subplot(111)
    # p0 = [amplitude, decay_constant, exponent]
    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
             x0, decay_constant,exponent)
    for ii,N in enumerate(Nlist):
        array = np.loadtxt('D:/Dropbox/QEC LT/Decoupling memory/XYdata/'+C2_msmts[str(N)],skiprows=1)
        x[:,ii]=array[:,0]
        y[:,ii]=array[:,1]
        y_err[:,ii] = array[:,2]

        # optimize.leastsq(fit_func,x[:,ii],y[:,ii],p0)
        p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
             x0, decay_constant,exponent)
        ax.errorbar(x[:,ii],y[:,ii],fmt='o',yerr=y_err[:,ii],label=str(N),color=color[ii])
        fit_result = fit.fit1d(x[:,ii],y[:,ii], None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
        Ts.append(fit_result['params_dict']['T'])
        Ts_err.append(fit_result['error_dict']['T'])
        Amps.append(fit_result['params_dict']['A'])
        Amps_err.append(fit_result['error_dict']['A'])

        plot.plot_fit1d(fit_result, np.linspace(0.,8,1001), ax=ax, plot_data=False,color=color[ii])
        # print fit_result
        #1.2 0.44047845 
        #1.8 0.44188627 
    # print Ts
    ax.set_xlim(0,8)
    ax.hlines([0.5],0,8,linestyles='dotted',linewidth = 2)
    ax.set_xlabel('Free evolution time (s)',fontsize = 25)
    ax.set_ylabel('Fidelity',fontsize = 25)
    plt.legend()


    logx = np.log10(Nlist)
    logy = np.log10(Ts)

    def fit_func(x, a, b):
        return a*x + b
    params = optimize.curve_fit(fit_func, logx, logy)
    print params[0]

    yerrp = []
    yerrm = []
    # print logx
    # print logy

    for a in range(len(Ts)):
        yerrp.append( abs(np.log(Ts[a]+Ts_err[a])-np.log(Ts[a])))
        yerrm.append( abs(np.log(Ts[a]-Ts_err[a])-np.log(Ts[a])))
    fig2 = plt.figure(figsize=(8,6))
    ax2 = fig2.add_subplot(111)
    ax2.errorbar(logx,logy,fmt='o',yerr=[yerrm, yerrp])
    xforplot = np.linspace(min(logx),max(logx)*1.1,101)
    ax2.plot(xforplot,fit_func(xforplot,params[0][0],params[0][1]))

    plt.show()


DD_scaling(C2_msmts)


def Carbon_T_mult(timestamp=None, older_than =None, posneg = True, folder_name = 'Hahn', measurement_name = 'adwindata', ssro_calib_timestamp =None,
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            decay_constant = 200, 
            exponent = 2, 

            partstr = 'part', plot_fit = True, do_print = True, fixed = [0,2], show_guess = False):
    ''' 
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''
    figsize=(6,4.7)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        timestamp, folder = toolbox.latest_data(folder_name, older_than=older_than, return_timestamp=True)
    try:
        Number_of_pulses = int(folder[folder.rfind('_N')+2:folder.rfind('_N')+4].replace('_',''))
        Adressed_carbon = int(folder[folder.rfind('_C')+2:folder.rfind('_C')+3])
        print 'N=', Number_of_pulses, '   C=', Adressed_carbon
        DD_msmt = True
    except:
        DD_msmt = False
    
    print 'Timestamp', timestamp
    print 'Folder', folder

    if partstr in folder:
        numberstart = folder.find(partstr)+len(partstr)
        numberofparts = int(folder[numberstart:len(folder)])
        basis_str_pos = folder[folder.rfind('\\')+7:numberstart]
        if posneg:
            posneg_str = 'posneg'
            if 'positive' in basis_str_pos:
                basis_str_neg = basis_str_pos.replace('positive', 'negative')
            else:
                basis_str_neg = basis_str_pos
                basis_str_pos = basis_str_neg.replace('negative', 'positive')
        else:
            if 'positive' in basis_str_pos:
                posneg_str = 'positive'
            else:
                posneg_str = 'negative'
    else:
        numberofparts = 1
        if 'positive' in folder:
            posneg_str = 'positive'
            basis_str_pos = folder[folder.rfind('\\')+7:len(folder)]
            basis_str_neg = basis_str_pos.replace('positive', 'negative')
        else:
            posneg_str = 'negative'
            basis_str_neg = folder[folder.rfind('\\')+7:len(folder)]
            basis_str_pos = basis_str_neg.replace('negative', 'positive')


        # basis_str_pos , basis_str_neg = basis_str_neg , basis_str_pos



    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print ssro_calib_folder

    cum_pts = 0
    print numberofparts
    if posneg:
        
        for kk in range(numberofparts):
            if partstr in folder:
                folder_pos = toolbox.latest_data(basis_str_pos+str(kk+1), older_than = older_than)
                folder_neg = toolbox.latest_data(basis_str_neg+str(kk+1), older_than = older_than)
            else:
                folder_pos = toolbox.latest_data(basis_str_pos, older_than = older_than)
                folder_neg = toolbox.latest_data(basis_str_neg, older_than = older_than)
            a = mbi.MBIAnalysis(folder_pos)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder)
            cum_pts += a.pts

            b = mbi.MBIAnalysis(folder_neg)
            b.get_sweep_pts()
            b.get_readout_results(name='adwindata')
            b.get_electron_ROC(ssro_calib_folder)
            print a.p0, b.p0
            if kk == 0:
                cum_sweep_pts = a.sweep_pts
                cum_p0 = (a.p0+(1-b.p0))/2.
                cum_u_p0 = np.sqrt(a.u_p0**2+b.u_p0**2)/2
                reps_per_datapoint = a.reps
            else:
                cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
                cum_p0 = np.concatenate((cum_p0, (a.p0+(1-b.p0))/2))
                cum_u_p0 = np.concatenate((cum_u_p0, np.sqrt(a.u_p0**2+b.u_p0**2)/2))

        a.pts   = cum_pts
        a.sweep_pts = cum_sweep_pts
        a.p0    = cum_p0
        a.u_p0  = cum_u_p0

    else:
        for kk in range(numberofparts):
            folder = toolbox.latest_data(basis_str_pos+str(kk+1), older_than = older_than)
            a = mbi.MBIAnalysis(folder)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC(ssro_calib_folder)
            reps_per_datapoint = a.reps
            cum_pts += a.pts

            if kk == 0:
                cum_sweep_pts = a.sweep_pts
                cum_p0 = a.p0
                cum_u_p0 = a.u_p0
            else:
                cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
                cum_p0 = np.concatenate((cum_p0, a.p0))
                cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

        a.pts   = cum_pts
        a.sweep_pts = cum_sweep_pts
        a.p0    = cum_p0
        a.u_p0  = cum_u_p0


    sorting_order=a.sweep_pts.argsort()
    a.sweep_pts.sort()
    a.p0=a.p0[sorting_order]
    print a.p0
    print a.u_p0
    a.u_p0=a.u_p0[sorting_order]

    ax=a.plot_results_vs_sweepparam(ret='ax',ax=None, 
                                    figsize=figsize, 
                                    ylim=(0.0,1.05)
                                    )
    fit_results = []
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    # print np.shape(a.sweep_pts[:]),np.shape(a.p0[:,0]),np.shape(a.u_p0[:])
    # print np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose()
    ax.plot(x,y)
    if DD_msmt:
        savestr = timestamp + '_C' + str(Adressed_carbon) + '_N' + str(Number_of_pulses) + '_' + posneg_str + '_Pts' + str(cum_pts) + '_Reps' + str(reps_per_datapoint) + '.txt'    
        save_folder_str = 'D:/Dropbox/QEC LT/Decoupling memory/XYdata/' + savestr
        np.savetxt(save_folder_str, np.vstack((a.sweep_pts[:],a.p0[:,0],a.u_p0[:,0])).transpose(),header=savestr)

    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
             x0, decay_constant,exponent)

         #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    ## plot data and fit as function of total time
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    fit_results.append(fit_result)

    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')

    return fit_results


def Carbon_T(timestamp=None, folder_name = 'Hahn', measurement_name = 'adwindata', ssro_calib_timestamp =None,
            offset = 0.5, 
            x0 = 0,  
            amplitude = 0.5,  
            decay_constant = 200, 
            exponent = 2, 
            plot_fit = True, do_print = True, fixed = [2], show_guess = False):
    ''' 
    Inputs:
    timestamp: in format yyyymmdd_hhmmss or hhmmss or None.
    measurement_name: list of measurement names
    List of parameters (order important for 'fixed') 
    offset, amplitude, decay_constant,exponent,frequency ,phase 
    '''
    figsize=(6,4.7)

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data(folder_name)

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Hans_sil1'
        print ssro_calib_folder

    print folder
    fit_results = []
    
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')

    a.get_electron_ROC(ssro_calib_folder)
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]

    ax=a.plot_results_vs_sweepparam(ret='ax',ax=None, 
                                    figsize=figsize, 
                                    ylim=(0.0,1.0)
                                    )

    
    
    ax.plot(x,y)

    p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, 
             x0, decay_constant,exponent)

         #plot the initial guess
    if show_guess:
        ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

    ## plot data and fit as function of total time
    if plot_fit == True:
        plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],1001), ax=ax, plot_data=False)

    fit_results.append(fit_result)

    plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    format='png')

    return fit_results

def fit1d(x, y, fitmethod, *arg, **kw):
    """
    example: from analysis.lib.fitting import fit,common
             x=np.array([0,1,2,3,4])
             y=np.array([2,12,22,32,42])
             fit_result=fit.fit1d(x,y,common.fit_line,2,8,ret=True,
                    fixed=[0],do_print=True)
             
    
    """
    # process known kws
    do_print = kw.pop('do_print', False)
    ret = kw.pop('ret', False)
    fixed = kw.pop('fixed', [])

    # err_y = kw.pop ('err_y', None)
    #if False :
    #    if (len(err_y) != len(y)):
    #       print 'Data and error arrays have non-matching lengths!'
    #       err_y = None

    # if (len(err_y) != len(y)):
    #   print 'Data and error arrays have non-matching lengths!'
    #   err_y = None

    # use the standardized fitmethod: any arg is treated as initial guess
    if fitmethod != None:
        p0, fitfunc, fitfunc_str = fitmethod(*arg)
    else:
        p0 = kw.pop('p0')
        fitfunc = kw.pop('fitfunc')
        fitfunc_str = kw.pop('fitfunc_str', '')        
    
    # general ability to fix parameters
    fixedp = []
    for i,p in enumerate(p0):
        if i in fixed:
            fixedp.append(p)
    for p in fixedp:
        p0.remove(p)
   
    # convenient fitting method with parameters; see scipy cookbook for details
    def f(params):
        i = 0
        for p in p0:
            p.set(params[i])
            i += 1

        # if (err_y != None):
        #   return ((y-fitfunc(x))/(err_y))
        # else:
        return y - fitfunc(x)

    if x is None: x = arange(y.shape[0])
    p = [param() for param in p0]
    
    # do the fit and process
    p1, cov, info, mesg, success = optimize.leastsq(f, p, full_output=True, maxfev=len(x)*100)
    if not success or cov == None: # FIXME: find a better solution!!!
        success = False
        print 'ERROR: Fit did not converge !'
        print 'reason: ',mesg
        #return False
    result = result_dict(p1, cov, info, mesg, success, x, y, p0, 
            fitfunc, fitfunc_str)

    # package the result neatly
    
    #print 'info',info
    #print 'p1',p1
    #print 'cov',cov

    #print 'dof',dof
    if do_print and success:
        print_fit_result(result)

    if ret:
        return result       

    return