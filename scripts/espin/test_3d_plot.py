import os
import numpy as np
from matplotlib import pyplot as plt
from measurement.lib.tools import toolbox
from analysis.lib.nv import nvlevels
from mpl_toolkits.mplot3d import Axes3D




#gate_V = np.array(
#		 [0,9,27,36,54,31.5,27,27,31.5, 22.5, 27, 54,27])
#E_y =  	 [52.4,52.67,52.91,54.28,52.07,54.89,55.13, 53.3,52.92,54.16, 53.47,51.89,55.19]
#E_x =    [53.57,53.52,56.17,58.55,60.47,57.96,57.36,56.35,56.94,55.87,56.41,59.84,57.35]

target_pt = [20,45]
target_Ey = [56.0,57.3] 
target_Ex = [57.9,59.2] 


D1_voltage = [-1200,-1200,-1200,-1200,-1239,-1239,-1400,-1400,-1400,-1400,-1400,-1400]
stripline_voltage = [350,370,330,171,171,400,400,350,300,250,200,150]

E_y = [56.29,56.06,55.39,55.76,55.65,54.30,54.28,54.62,55.10,55.45,55.71,55.95] 
E_x = [57.61,57.70,57.19,56.38,56.24,57.34,58,57.73,57.57,57.28,57.13,57.14]

E_p = [49.6,49.77,49.03,49.01,48.91,48.59,48.78,48.83,49.15,49.23,49.33,49.33]


title= 'Pippin test'
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
#ax.set_ylim([73,8])

ax.scatter(stripline_voltage,D1_voltage,E_y,zdir=u'z',s=9,c=u'b')
ax.scatter(stripline_voltage,D1_voltage,E_x,zdir=u'z',s=9,c=u'r')
ax.scatter(stripline_voltage,D1_voltage,E_p,zdir=u'z',s=9,c=u'g')

#ax.plot(gate_V, E_x, 'bo', ms = 8, label = '$E_x$')
#ax.plot(gate_V, E_y, 'ro', ms = 8, label = '$E_y$')
#ax.plot(target_pt, target_Ey, 'mo', ms = 10, label = '$target$')
#ax.plot(target_pt, target_Ex, 'go', ms = 10, label = '$target$')


ax.set_xlabel('gate voltage [V]')
ax.set_ylabel('frequency [GHz]')
ax.set_title(title)

ax.legend(loc=0)