import os, sys
import numpy as np
import h5py
import logging
import sympy

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, esr
from analysis.lib.tools import plot
from analysis.lib.math import error

### settings
timestamp = None


### fitfunction

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)

else:
    #folder = toolbox.latest_data('PostInitDarkESR')
    folder = toolbox.latest_data()

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')

a.get_electron_ROC()
ax = a.plot_results_vs_sweepparam(ret='ax',name='adwindata')
x = a.sweep_pts
y = a.p0.reshape(-1)
y_u = a.u_p0.reshape(-1)

average 	= np.average(y) 
average_u 	= np.sqrt(np.sum(y_u**2))/np.size(y_u) 

print 'Average = %s +/- %s' %(average,  average_u)

ax.set_ylim(-0.05,1.05)

plt.savefig(os.path.join(folder, 'mbi_darkesr_analysis.pdf'),
        format='pdf')
plt.savefig(os.path.join(folder, 'mbi_darkesr_analysis.png'),
        format='png')
