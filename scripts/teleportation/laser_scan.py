import os
import numpy as np
from matplotlib import pyplot as plt
from measurement.lib.tools import toolbox
from analysis.lib.nv import nvlevels

# settings
timestamp = None


### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('LaserFrequencyScan')

# get data
fn = os.path.join(folder, os.path.split(folder)[1]+'.dat')
d = np.loadtxt(fn)

frq = d[:,1]
cts = d[:,2]

fig = plt.Figure()

