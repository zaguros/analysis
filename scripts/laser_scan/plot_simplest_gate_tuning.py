import os
import numpy as np
from matplotlib import pyplot as plt
from measurement.lib.tools import toolbox
from analysis.lib.nv import nvlevels

title= 'Pippin Sil5 Voltage on 22 MW stripline'

#gate_V = [0,100*90./2000.,500*90./2000.]
#E_y = [47.9,47.3,45.84]
#E_x = [54.6,55.0,56.7]

#gate_V = [0,-100,-250, -500, -700,-900, -1000]*90./2000.
#E_y = [47.11, 47.58, 48.3,  49.8,  50.74, 52.09,52.48]
#E_x = [53.53, 53.18, 53.27, 52.82, 52.97, 52.97, 53.91]

#gate_V = np.array([-700,-800,-900, -1000, -1100,-1200, -1300, -1500, -1700,-2000])
#gate_V= gate_V*90./2000.
#E_y =  	 [51.08, 51.53, 51.96,  52.33,  52.61, 53.03, 53.36, 53.58, 53.5,54.84]
#E_x =    [53.75, 53.7,  53.73,  53.79,  54.15, 54.57, 55.04, 56.01,57.04,59.09]

gate_V = np.array([0.0,9.0,18.0,27,45.0,63,90])
E_y =  	 [53.75,54.5,55.5,56.66,58.9,59.6,57.24]
E_x =    [64.46,63.82,63.3,63.0,63.54,64.65,68.32]



fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(gate_V, E_x, 'bo', ms = 8, label = '$E_x$')
ax.plot(gate_V, E_y, 'ro', ms = 8, label = '$E_y$')

ax.set_xlabel('gate voltage [V]')
ax.set_ylabel('frequency [GHz]')
ax.set_title(title)


ax.legend()