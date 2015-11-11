import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot
from analysis.lib.math import error
reload(plot)
# fit_startup = False
def simple_plot_contrast(timestamps = [None,None], tag = '', measurement_name = ['adwindata'],folder_name ='Tomo',
        ssro_calib_timestamp =None, save = True,
        do_plot = True, return_data = False,
        guess_frq = 1/4000.,
        guess_amp = 1.5,
        guess_k = 0.,
        guess_phi = 0.,
        guess_o = 0.3):

    ### SSRO calibration
    if ssro_calib_timestamp == None: 
        ssro_calib_folder = toolbox.latest_data('SSRO')
    else:
        ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil8'

    if timestamps[0] == None: 
        folder_a = toolbox.latest_data(contains='positive' + tag)
        folder_b = toolbox.latest_data(contains='negative' + tag)
    elif len(timestamps)==1:        
        folder_b = toolbox.data_from_time(timestamps[0])      
        print folder_b
        folder_a = toolbox.latest_data(contains = 'pos', older_than = timestamps[0])   
        print folder_a
    else:
        folder_a = toolbox.data_from_time(timestamps[0])      
        folder_b = toolbox.data_from_time(timestamps[1])           
    
    a = mbi.MBIAnalysis(folder_a)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder)
    y_a= ((a.p0.reshape(-1)[:])-0.5)*2
    y_err_a = 2*a.u_p0.reshape(-1)[:] 


    b = mbi.MBIAnalysis(folder_b)
    b.get_sweep_pts()
    b.get_readout_results(name='adwindata')
    b.get_electron_ROC(ssro_calib_folder)
    y_b= ((b.p0.reshape(-1)[:])-0.5)*2
    y_err_b = 2*b.u_p0.reshape(-1)[:] 

    x = a.sweep_pts.reshape(-1)[:]
    # x = range(len(y_a)) 
    # ax = a.plot_results_vs_sweepparam(ret='ax', name = 'adwindata' )


    
    ### Combine data
    y = (y_a - y_b)/2.
    y_err =  1./2*(y_err_a**2 + y_err_b**2)**0.5 

    if x[-1] < x[0]:
        x = x[::-1]
        y = y[::-1]

    fig,ax = plt.subplots()

    o = fit.Parameter(guess_o, 'o')
    f = fit.Parameter(guess_frq, 'f')
    A = fit.Parameter(guess_amp, 'A')
    phi = fit.Parameter(guess_phi, 'phi')
    k = fit.Parameter(guess_k, 'k')
    p0 = [A,o,f,phi]
    fitfunc_str = ''

    ax.errorbar(x,y,yerr = y_err, marker = 'o',ls = '')

    def fitfunc(x):
        return (o()-A()) + A() * np.exp(-(k()*x)**2) * np.cos(2*np.pi*(f()*x - phi()))

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, fixed=[],
            do_print=True, ret=True)

    ax.set_title(folder_a)

    
    

    y_fit = fit_result['fitfunc'](np.arange(0,60,1))

    for i, y in enumerate(abs(y_fit)):
        if y == min(abs(y_fit)):
            x_opt = i


    ax.text(5,0.9,'pulses for pi/2: ' +str(x_opt))

    plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax,
           plot_data=False)

    ax.set_ylim([-1,1])

    plt.savefig(os.path.join(folder_a, 'mbi_erabi_analysis.pdf'),
            format='pdf')
    plt.savefig(os.path.join(folder_a, 'mbi_erabi_analysis.png'),
            format='png')




# print 'The pulse length shift is:' + str(t_shift)

# else:
    # fit_result = fit.fit1d(a.sweep_pts[:], a.p0.reshape(-1)[:], rabi.fit_rabi_fixed_upper,
        # guess_frq, guess_amp, guess_phi, guess_k, fixed=[],
        # do_print=True, ret=True)
    # plot.plot_fit1d(fit_result, np.linspace(0,a.sweep_pts[-1], 201), ax=ax,
        # plot_data=False)

    # plt.savefig(os.path.join(folder, 'electronrabi_analysis.pdf'),
        # format='pdf')



### FFT
# p0_fft = np.fft.fft(a.p0.reshape(-1), n=32)
# frq = np.fft.fftfreq(32, d=(a.sweep_pts[1]-a.sweep_pts[0])/1e3)
#
# fig = a.default_fig()
# ax = a.default_ax(fig)
# ax.plot(frq, p0_fft, 'o-')
