import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt


N_HF_splt = 2.17e-3
detuning = 1.0e-6

timestamp = '20141016151548' #'20140715095713' #None #
guess_f1 = 1.10e-3#detuning - N_HF_splt
guess_A1 = 0.15
guess_phi1 = 0.
guess_f2 = 0.91e-3#detuning
guess_A2 = 0.15
guess_phi2 = 0.
guess_f3 = detuning + N_HF_splt
guess_A3 = -0.15
guess_phi3 = 0.


guess_tau = 4560
guess_a = 0.5

c_green = (9/255.,232/255.,94/255.)
c_grey = (64/255.,78/255.,77/255.)#(240/255.,242/255.,166/255.)
c_blue = (68/255.,204/255.,255/255.)
c_red = (150/255.,52/255.,132/255.)
c_orange = (242/255.,129/255.,35/255.)
c_orange_2 = (242/255.,129/255.,35/255.)

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('ElectronRamsey')
print folder
a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC(r'D:\measuring\data\20141016\150451_AdwinSSRO_SSROCalibration_111_1_sil18')
x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]
y_err = a.u_p0.reshape(-1)[:]
# ax = a.plot_results_vs_sweepparam(ret='ax', name = 'adwindata')
p0, fitfunc_0, fitfunc_str_0 = ramsey.fit_ramsey_gaussian_decay(guess_tau, guess_a, (guess_f1, guess_A1, guess_phi1), 
	(guess_f2,guess_A2, guess_phi2))
fit_xvals_0=np.linspace(0,a.sweep_pts[-1],1000)
# ax.plot(fit_xvals_0,fitfunc_0(fit_xvals_0), 'r--', lw=1)
# #fit_xvals=np.linspace(res['x'][0],res['x'][-1],fit_num_points)
fit_result = fit.fit1d(a.sweep_pts, a.p0.reshape(-1), ramsey.fit_ramsey_gaussian_decay,
        guess_tau, guess_a, (guess_f1, guess_A1, guess_phi1), (guess_f2,guess_A2, guess_phi2), fixed=[4,7], #,4,7,10],
        do_print=True, ret=True)
# #fit_result = False
# if fit_result != False :
# 	plot.plot_fit1d(fit_result, np.linspace(0,a.sweep_pts[-1],201), ax=ax,
#         plot_data=False)

fig, ax = plt.subplots(figsize = (10,5))
# ax.plot(np.linspace(x[0],x[-1],201), fitfunc_0(np.linspace(x[0],x[-1],201)), ':', lw=2)




res = fit_result
fit_xvals=np.linspace(res['x'][0],res['x'][-1],1001)
ax.plot(fit_xvals*1e-3, res['fitfunc'](fit_xvals), linestyle = '-',color = c_red, linewidth = 2 )
x = x*1e-3
errlines = ax.errorbar(x,y,yerr = y_err, color = c_red,ls = '', marker = 'o',markersize = 7,markeredgecolor = c_red,capsize=5)
# set(errlines, linewidth=6)
xticks = np.arange(0,x[-1]+1,2)
yticks = np.linspace(0,1,3)
ax.set_ylim([0,1])
ax.set_xlim([0-0.2,x[-1]+0.2])
plt.xticks(xticks)
plt.yticks(yticks)
ax.set_yticks(np.arange(-1,1.1,0.2), minor = True)
ax.set_xticks(np.arange(0,x[-1]+1,1), minor = True)
ax.set_xlabel('Free evolution time (ms)',fontsize = 30)
ax.set_ylabel('State fidelity',fontsize = 30)
# ax.hlines([0.5],x[0]-1,x[-1]+1,linestyles='dotted',linewidth = 1.5)
ax.set_ylim([0,1])
ax.tick_params(axis='x', which='major', labelsize=30)
ax.tick_params(axis='y', which='major', labelsize=30)
mpl.rcParams['axes.linewidth'] = 1
ax.tick_params('both', length=4, width=1, which='major')
folder = r'D:\measuring\data\QEC_data\figs\final figures'
plt.savefig(os.path.join(folder, 'eramsey.pdf'),format='pdf',bbox_inches='tight')
plt.savefig(os.path.join(folder, 'eramsey.png'),format='png',bbox_inches='tight')




### FFT
# p0_fft = np.fft.fft(a.p0.reshape(-1), n=32)
# frq = np.fft.fftfreq(32, d=(a.sweep_pts[1]-a.sweep_pts[0])/1e3)
#
# fig = a.default_fig()
# ax = a.default_ax(fig)
# ax.plot(frq, p0_fft, 'o-')
