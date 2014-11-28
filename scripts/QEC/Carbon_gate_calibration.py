import numpy as np
import os,sys

sys.path.append("/measuring/")
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(plot)


def Carbon_control_sweep_N_zoom(timestamp=None, measurement_name = ['adwindata'], 
            A = [0.5, 0.5],
            fitfunc_type = 'single', plot_fit = True, do_print = False, show_guess = True):
    ''' Function to analyze data for optimization of the number of pulses for a controlled C13 gate. 
    '''
    
    if timestamp != None:
        folder1 = toolbox.data_from_time(timestamp)
    else:
        folder1 = toolbox.latest_data('Decoupling')

    folder_list=[folder1]


    for k in ['pi','x','-x','y','-y']:
        folder1 = toolbox.latest_data(contains ='sweep_N_'+ k +'_tau', older_than=timestamp,return_timestamp = False)
        folder_list=folder_list+[folder1]




    fit_results = []
    for i,folder in enumerate(folder_list):
        for k in range(0,len(measurement_name)):
            a = mbi.MBIAnalysis(folder)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC()
            if i==0:
                ax = a.plot_results_vs_sweepparam(ret='ax')
                ax.set_ylim(-0.05,1.05)

            x = a.sweep_pts.reshape(-1)[:]
            y = a.p0.reshape(-1)[:]
            y_u = a.u_p0.reshape(-1)[:]


            p0, fitfunc, fitfunc_str = common.fit_poly(A)
            
            if show_guess:
                ax.plot(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(x[0],x[-1],201)), ':', lw=2)
            
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[3],err_y=y_u)
            if plot_fit == True:
                plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax, plot_data=True,print_info=False)

            fit_results.append(fit_result)
            print folder


    #Let's calculate crossing points of the various RO types. x vs -x, y vs -y and pi vs no pulse.
    # the coefficient a0 gives the offfset and a1 gives the slope. we want a0+a1*x=a0'+a1'*x
    # solving for x gives: x=(a0-a0')/(a1'-a1)
    comparison_list=['pi vs no pulse', 'x vs -x', 'y vs -y']
    for i,f in enumerate([[0,1],[2,3],[4,5]]):
        fitres1=fit_results[f[0]]['params_dict']
        fitres2=fit_results[f[1]]['params_dict']
        err1=fit_results[f[0]]['error']
        err2=fit_results[f[1]]['error']
        intersection = (fitres1['a0']-fitres2['a0'])/(fitres2['a1']-fitres1['a1'])

        #calculate the error on the intersection via gaussian propagation of uncertainty.
        err_a01=(1./(fitres2['a1']-fitres1['a1']))*err1[0]
        err_a11=((fitres1['a0']-fitres2['a0'])*(1/(fitres2['a1']-fitres1['a1']))**2)*err1[1]
        err_a02=(-1./(fitres2['a1']-fitres1['a1']))*err2[0]
        err_a12=(-(fitres1['a0']-fitres2['a0'])*(1/(fitres2['a1']-fitres1['a1']))**2)*err2[1]

        intersection_error = np.sqrt(err_a01**2+err_a11**2+err_a02**2+err_a12**2)

        print 'The point of intersection for the curves ' + comparison_list[i] + ' is: '
        print str(intersection) + ' +- ' + str(intersection_error)


    plt.savefig(os.path.join(folder_list[0], 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder_list[0], 'analyzed_result.png'),
    format='png')

    # diff = np.abs(y - 0.5)
    # print diff
    # print 'Optimum number of pulses N = ' + str(x[np.argmin(diff)])
    # print 'with y-0.5 = ' + str(y[np.argmin(diff)]-0.5) + ' +/- ' + str(y_u[np.argmin(diff)])

            # freq = fit_results[0]['params_dict']['f1']
            # period = 1/freq 
            # print 'Period is %s pulses ' %(period)
            # # N_pi = round(period*.5/2)*2.0 
            # N_pi2 = round(period/2*.25)*2.0
            # # print 'Pi pulse: %s pulses' %N_pi
            # print 'Pi2 pulse: %s pulses' %N_pi2
    return fit_results

Carbon_control_sweep_N_zoom(timestamp='20141128_101807')

