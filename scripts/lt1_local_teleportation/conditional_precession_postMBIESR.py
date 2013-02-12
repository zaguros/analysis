import numpy as np
from matplotlib import pyplot as plt
from analysis.lib.m2.ssro import mbi
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot

timestamp = None # '20130207202558'

if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('PostInitDarkESR')

fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)

a = mbi.PostInitDarkESRAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results()

ax.errorbar(a.sweep_pts, a.normalized_ssro[:,0], fmt='o',
        yerr=a.u_normalized_ssro[:,0])


xleft,yleft = a.sweep_pts[:15], a.normalized_ssro[:15,0]
xright,yright = a.sweep_pts[15:], a.normalized_ssro[15:,0]

fitleft = fit.fit1d(xleft, yleft, common.fit_gauss, 
        0.81, -0.2, 2856.6e6, 0.5e6,
        do_print=True, ret=True)
if fitleft != None and fitleft != False:
    plot.plot_fit1d(fitleft, xleft, ax=ax, print_info=False)
    leftdip = fitleft['params_dict']['a'] + fitleft['params_dict']['A']
else:
    leftdip = min(yleft)

fitright = fit.fit1d(xleft, yleft, common.fit_gauss, 
        0.81, -0.2, 2856.6e6 + 12.8e6, 0.5e6,
        do_print=True, ret=True)
if fitright != None and fitright != False:
    plot.plot_fit1d(fitright, xright, ax=ax, print_info=False)
    rightdip = fitright['params_dict']['a'] + fitright['params_dict']['A']
else:
    rightdip = min(yright)


plt.text(2.852e9, 0.55, 'left: %.2f, right: %.2f, ratio(right): %.2f' \
        % (leftdip, rightdip, rightdip/(leftdip+rightdip)), 
        ha='left', va='bottom')

ax.set_xlim(2.85e9, 2.875e9)
ax.set_ylim(0.5, 0.9)
ax.set_title(a.timestamp+'\n'+a.measurementstring)
ax.set_xlabel('MW frequency (Hz)')
ax.set_ylabel(r'uncorrected fidelity $F(|0\rangle)$')


