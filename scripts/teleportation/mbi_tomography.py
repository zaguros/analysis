import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot
from analysis.lib.math import error


timestamp = None

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('MBIMeasurement')

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results()
a.get_N_ROC(0.99, 0.02, 0.94, 0.01, 0.96, 0.01)

a.save()    
a.plot_results_vs_sweepparam()

x11=a.p0.reshape(-1)[0]
x12=0.5-a.p0.reshape(-1)[1]
y12=a.p0.reshape(-1)[2] -0.5
u_x11=a.u_p0.reshape(-1)[0]
u_x12=a.u_p0.reshape(-1)[1]
u_y12=a.u_p0.reshape(-1)[2]
densitymatrix=array([[x11,x12+y12*1j],[x12-y12*1j,1-x11]])
plot.plot_dm(densitymatrix)
densitymatrix2=np.matrix([[x11,x12+y12*1j],[x12-y12*1j,1-x11]])
print 'Density Matrix:'
print densitymatrix2
print 'Error in [x11,x12,y12]:'
print [u_x11,u_x12,u_y12]


a.finish()



# zip(a.sweep_pts. a.p0)

# ax = a.plot_results_vs_sweepparam(ret='ax', )

# fit_result = fit.fit1d(a.sweep_pts, a.p0.reshape(-1), rabi.fit_rabi_fixed_upper,
#         guess_frq, guess_amp, guess_phi, guess_k, fixed=[2,3],
#         do_print=True, ret=True)
# plot.plot_fit1d(fit_result, np.linspace(0,a.sweep_pts[-1],201), ax=ax,
#         plot_data=False)
# 
# plt.savefig(os.path.join(folder, 'electronrabi_analysis.pdf'),
#         format='pdf')


