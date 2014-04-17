'''Script to analyze the dynamical decoupling data
by THT'''

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
from matplotlib.widgets import Cursor

reload(common)

## Data location ##
measurement_name = ['adwindata']

timestamp = ['20140407_223450', '20140407_223801', '20140407_224119',
'20140407_224446', '20140407_231158','20140407_225215', '20140407_225614',
 '20140407_230023', '20140407_231645', '20140407_232118', '20140407_232603',
 '20140407_233057', '20140407_233605', '20140407_234654', '20140408_111324',
 '20140408_111725', '20140408_112126', '20140408_112529', '20140408_112930',
 '20140408_113614', '20140408_114015', '20140408_114416', '20140408_114818',
 '20140408_115220', '20140408_115622', '20140408_120024', '20140408_120426',
 '20140408_120825', '20140408_130753']#,

cum_pts = 0
cum_normalized_ssro = np.empty(0)

for kk in range(len(timestamp)):
    folder = toolbox.data_from_time(timestamp[kk])

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
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

#ax = a.plot_results_vs_sweepparam(ret='ax')

fig = a.default_fig(figsize=(18,4))
ax = a.default_ax(fig)
ax.set_xlim(4,6)
ax.plot(a.sweep_pts, a.p0, '.-b', lw=1)


timestamp = ['20140411_101732','20140411_102646','20140411_104011',
'20140411_105345'] #with 32 pulses
cum_pts = 0
cum_normalized_ssro = np.empty(0)

for kk in range(len(timestamp)):
    folder = toolbox.data_from_time(timestamp[kk])

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
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


ax.plot(a.sweep_pts, a.p0, '.-r', lw=1)


plt.savefig(os.path.join(folder, 'combined_result.pdf'),
format='pdf')
plt.savefig(os.path.join(folder, 'combined_result.png'),
format='png')

zoom = 13.0
fig = a.default_fig(figsize=(8,4))
ax = a.default_ax(fig)
ax.set_xlim(12,14)
ax.plot(a.sweep_pts, a.p0, '.-b', lw=1)








