import os
import numpy as np
from matplotlib import pyplot as plt
from measurement.lib.tools import toolbox
from analysis.lib.nv import nvlevels



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

#gate_V = np.array([0.0,9.0,18.0,27,45.0,63,90])
#E_y =  	 [53.75,54.5,55.5,56.66,58.9,59.6,57.24]
#E_x =    [64.46,63.82,63.3,63.0,63.54,64.65,68.32]

#gate_V = np.array([0.0,-9.0,-18.0,-27,-36.0,-45.0,-90])
#E_y =  	 [53.42,54.80,56.32,57.81,58.81,59.43,54.44]
#E_x =    [64.37,63.63,62.62,62.94,63.19,64.18,72.09]


#gate_V = np.array([0.0,9.0,18.0])
#E_y =  	 [56.23,55.79,54.79]
#E_x =    [68.01,68.51,69.4]

#gate_V = np.array([0.0,-9.0,-18.0,-36,-45,-54,-63,-72,-81,-90])
#E_y =  	 [56.36,57.41,58.29,60.01,60.82,61.54,62.67,63.53,63.72,62.99]
#E_x =    [68.20,67.34,66.51,64.77,64.15,63.53,65.64,64.75,64.64,65.53]

#gate_V = np.array([-54,-36,-18,0,18])
#E_y =  	 [55.27,55.90,56.36,56.24,55.73]
#E_x =    [68.51,68.14,67.94,68.09,68.60]
#
#gate_V = np.array([0,9,27,45,54,63,72])
#E_y =  	 [56.05,56.74,57.96,57.45,56.79,58.82,57.76]
#E_x =    [67.76,67.07,65.24,65.16,65.54,66.62,67.31]

#gate_V = np.array([0,9,13.5,18,27])
#E_y =  	 [72.14,73.36,73.37,73.35,72.71]
#E_x =    [74.61,73.71,73.91,73.95,74.59]

#gate_V = np.array(
#		 [0,	-18,	-36,	-54,	-72,	-90])
#gate_V = np.array(
#		 [0,	18,		36])
#E_y =  	 [45.1,	45.2,	44.5]
#E_x =    [53.2,	53.5,	53.9]


gate_V = np.array(
		 [-18,-36,-54,-72,-90])
E_y =  	 [67.22,68.42,69.79,71.09,71.23]
E_x =    [73.2,72.85,72.77,73.01,74.42]
#gate_V = np.array([ -12,	-9,	-6,	0,	0,	0,	5,	9,	18])
#E_y =  	 [74.37, 74.59, 74.41, 74.35, 74.51, 74.45, 74.51, 74.71, 75.06]
#E_x =    [77.19, 77.10,  -1, -1 ,-1 , -1 , 76.5, 76.85, 77.93]
#E_p = 	 [68.51, 68.68, 68.35, 68.14, 68.33, 68.33, -1, 68.62, 69.28]

#Yellow_tuning_on_Ey = [27.66, 27.73, 26.56, 25.72,26.10,26.04, 26.36, 27.62, 30.64]
#Yellow_tuning_on_Ex = [27.82, 28.06, -1,-1,-1,-1, 26.98,28.04, 30.80 ]


title= 'the111no1 Sil8 gate 21 during SSRO'
fig = plt.figure()
ax = fig.add_subplot(111)
#ax.set_ylim([73,8])

ax.plot(gate_V, E_x, 'bo', ms = 8, label = '$E_x$')
ax.plot(gate_V, E_y, 'ro', ms = 8, label = '$E_y$')
#ax.plot(gate_V, E_p, 'go', ms = 8, label = '$E\'$')


ax.set_xlabel('gate voltage [uW]')
ax.set_ylabel('frequency [GHz]')
ax.set_title(title)

ax.legend(loc=0)


#fig2 = plt.figure()
#ax2 = fig.add_subplot(111)
#ax2.set_ylim([25,32])

#ax2.plot(gate_V, Yellow_tuning_on_Ex, 'bo', ms = 8, label = 'tun.  on $E_x$')
#ax2.plot(gate_V, Yellow_tuning_on_Ey, 'ro', ms = 8, label = 'tun. on $E_y$')

#ax.set_xlabel('gate voltage [uW]')
#ax.set_ylabel('frequency [GHz]')
#ax.set_title(title)


#ax2.legend(loc=0)