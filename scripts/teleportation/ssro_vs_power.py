import os, sys
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib.fitting import fit, common
from measurement.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.m2.ssro import ssro


timestamp = None
analyze_runs = True


if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('AdwinSSRO_RO_saturation_power')

if analyze_runs:
    a = ssro.SSROAnalysis(folder)

    for g in a.g.items():
        gn = g[0]
        a.get_run(gn)
        a.cpsh_hist(a.ro_counts, a.reps, name=gn)
        a.readout_relaxation(a.ro_time, a.ro_counts, a.reps, a.binsize, name=gn)
        a.spinpumping(a.sp_time, a.sp_counts, a.reps, a.binsize, name=gn)
        a.charge_hist(a.cr_counts, name=gn)
        a.fidelity(a.ro_counts, a.reps, a.binsize, 0, name=gn)
        
    plt.close('all')
    a.finish()

