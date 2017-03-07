import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt

reload(common)

def electron_DD_analysis(timestamp=None, measurement_name = ['adwindata'], offset = 0.5, 
                        amplitude = 0.5, position =0, T2 = 800, 
                        power=2, plot_fit = True, do_print = False, 
                        show_guess = False):
    ''' Function to analyze simple decoupling measurements. Loads the results and fits them to a simple exponential.
    Inputs:
    timestamp: list of timestamps in format in format yyyymmdd_hhmmss or hhmmss or None. (Added: List functionality. NK 20150320)
    measurement_name: list of measurement names
    Based on electron_T1_anal,
    '''

    if timestamp != None:
        if type(timestamp)==str:
            folder_list=[toolbox.data_from_time(timestamp)]
        else:
            folder_list = []
            for t in timestamp:
                folder_list.append(toolbox.data_from_time(t))
    else:
        folder_list = [toolbox.latest_data('Decoupling')]

    fit_results = []
    for ii,f in enumerate(folder_list):
        for k in range(0,len(measurement_name)):
            a = mbi.MBIAnalysis(f)
            a.get_sweep_pts()
            a.get_readout_results(name=measurement_name[k])
            a.get_electron_ROC()
            if ii==0:
                ax = a.plot_results_vs_sweepparam(ret='ax')
            else:
                ax.errorbar(a.sweep_pts.reshape(-1)[:],a.p0.reshape(-1)[:],yerr=a.u_p0.reshape(-1)[:],fmt='o')
            x = a.sweep_pts.reshape(-1)[:]
            y = a.p0.reshape(-1)[:]

            p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)

            #plot the initial guess
            if show_guess:
                ax.plot(np.linspace(0,x[-1],201), fitfunc(np.linspace(0,x[-1],201)), ':', lw=2)

            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=[0,2])

            ## plot data and fit as function of total time
            if plot_fit == True:
                plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)

            fit_results.append(fit_result)

    plt.savefig(os.path.join(folder_list[0], 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder_list[0], 'analyzed_result.png'),
    format='png')

            ## plot data and fit as function of tau
            #plt.figure()
            #plt.plot(np.linspace(0,x[-1],201), fit_result['fitfunc'](np.linspace(0,x[-1],201)), ':', lw=2)
            #plt.plot(x, y, '.', lw=2)




    return fit_results





'''
# fit_startup = False

timestamp = '20140403184323' #'20130907183620' # None
guess_frq = 1./200.
guess_amp = 0.5
guess_k = 0.
guess_phi = 0.
guess_o = 1.

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
ax = a.plot_results_vs_sweepparam(ret='ax', )

t_shift = 0
# t_start=0
# t_0=10

x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]
'''
''' add the fit function:
'''
'''
o = fit.Parameter(guess_o, 'o')
f = fit.Parameter(guess_frq, 'f')
A = fit.Parameter(guess_amp, 'A')
phi = fit.Parameter(guess_phi, 'phi')
k = fit.Parameter(guess_k, 'k')
p0 = [f, A]
fitfunc_str = ''

def fitfunc(x) :
    return (o()-A()) + A() * exp(-k()*x) * cos(2*pi*(f()*x - phi()))

fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)

plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax,
        plot_data=False)

plt.savefig(os.path.join(folder, 'mbi_erabi_analysis.pdf'),
        format='pdf')
plt.savefig(os.path.join(folder, 'mbi_erabi_analysis.png'),
        format='png')
'''



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
