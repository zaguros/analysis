'''Script to analyze the dynamical decoupling data
by THT'''

import numpy as np
import os
from measurement.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
from matplotlib.widgets import Cursor

reload(common)

## Data location ##
measurement_name = ['adwindata']

timestamp = ['20140407_223450', '20140407_223801', '20140407_224119', '20140407_224446', '20140407_231158','20140407_225215', '20140407_225614', '20140407_230023', '20140407_231645', '20140407_232118', '20140407_232603', '20140407_233057', '20140407_233605', '20140407_234654', '20140408_111324', '20140408_111725', '20140408_112126', '20140408_112529', '20140408_112930', '20140408_113614', '20140408_114015', '20140408_114416', '20140408_114818', '20140408_115220', '20140408_115622', '20140408_120024', '20140408_120426','20140408_120825', '20140408_130753']#,
timestamp2 = ['20140411_101732','20140411_102646', '20140411_103753']

cum_pts_a = 0
cum_normalized_ssro_a = np.empty(0)

cum_pts_c = 0
cum_normalized_ssro_c = np.empty(0)

for kk in range(len(timestamp)):
    folder = toolbox.data_from_time(timestamp[kk])

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    cum_pts_a += a.pts

    if kk == 0:
        cum_sweep_pts_a = a.sweep_pts
        cum_p0_a = a.p0
        cum_u_p0_a = a.u_p0
    else:
        cum_sweep_pts_a = np.concatenate((cum_sweep_pts_a, a.sweep_pts))
        cum_p0_a = np.concatenate((cum_p0_a, a.p0))
        cum_u_p0_a = np.concatenate((cum_u_p0_a, a.u_p0))

a.pts_a   = cum_pts_a
a.sweep_pts_a = cum_sweep_pts_a
a.p0_a    = cum_p0_a
a.u_p0_a  = cum_u_p0_a

fig = a.default_fig(figsize=(18,4))
ax = a.default_ax(fig)
ax.plot(a.sweep_pts_a, a.p0_a, '.-b', lw=1)

for cc in range(len(timestamp2)):
    folder = toolbox.data_from_time(timestamp[cc])

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    cum_pts_c += a.pts

    if cc == 0:
        cum_sweep_pts_c = a.sweep_pts
        cum_p0_c = a.p0
        cum_u_p0_c = a.u_p0
    else:
        cum_sweep_pts_c = np.concatenate((cum_sweep_pts_c, a.sweep_pts))
        cum_p0_c = np.concatenate((cum_p0_c, a.p0))
        cum_u_p0_c = np.concatenate((cum_u_p0_c, a.u_p0))

a.pts_c   = cum_pts_c
a.sweep_pts_c = cum_sweep_pts_c
a.p0_c    = cum_p0_c
a.u_p0_c  = cum_u_p0_c

ax.plot(a.sweep_pts_c, a.p0_c, '.-r', lw=1)

cursor = Cursor(ax, useblit=True, color='gray', linewidth =1)

print folder
plt.savefig(os.path.join(folder, 'combined_result.pdf'),
format='pdf')
plt.savefig(os.path.join(folder, 'combined_result.png'),
format='png')

zoom = 13.0
fig = a.default_fig(figsize=(8,4))
ax = a.default_ax(fig)
ax.set_xlim(zoom-1000e-3,zoom+1000e-3)
ax.plot(a.sweep_pts, a.p0, '.-b', lw=1)








