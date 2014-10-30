import os
import numpy as np
from matplotlib import pyplot as plt
from measurement.lib.tools import toolbox
from analysis.lib.nv import nvlevels

# settings
timestamp = None

Ex = 86.72
Ey = 82.41
Bz=402.69
# lines = nvlevels.get_ES_ExEy(Ex,Ey)
lines = np.sort(nvlevels.get_optical_transitions(E_field=[(Ex-Ey)/2.,0.,0.],B_field = [0.,0.,Bz],Ee0=-1.94,
					show_E_transitions=True,show_A_transitions=True,show_FB_E_transitions=False, 
                            show_FB_A_transitions=False, show_E_prime_flip_transitions=False))+Ex-1.192
Eflip_lines = np.sort(nvlevels.get_optical_transitions(E_field=[(Ex-Ey)/2.,0.,0.],B_field = [0.,0.,Bz],Ee0=-1.94,
					show_E_transitions=False,show_A_transitions=False,show_FB_E_transitions=False, 
                            show_FB_A_transitions=False, show_E_prime_flip_transitions=True))+Ex-1.192
FB_lines = np.sort(nvlevels.get_optical_transitions(E_field=[(Ex-Ey)/2.,0.,0.],B_field = [0.,0.,Bz],Ee0=-1.94,
					show_E_transitions=False,show_A_transitions=False,show_FB_E_transitions=False, 
                            show_FB_A_transitions=True, show_E_prime_flip_transitions=False))+Ex-1.192
### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('LaserFrequencyScan')

# get data
fn = os.path.join(folder, os.path.split(folder)[1]+'.dat')
d = np.loadtxt(fn)

_folder, _file = os.path.split(fn)
_folder, _ = os.path.split(_folder)
_folder, _date = os.path.split(_folder)
title = _date+'/'+_file[:-4]

frq = d[:,1]
cts = d[:,2]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(frq, cts, 'k-', lw=2)

ax.set_xlabel('frequency (GHz) - 470.4 THz')
ax.set_ylabel('counts (Hz)')
ax.set_title(title)

y0,y1 = ax.get_ylim()
ax.vlines(lines, 0, y1, colors='r', lw=1)
ax.vlines(Eflip_lines, 0, y1, colors='b', lw=1)
ax.vlines(FB_lines, 0, y1, colors='g', lw=1)

# ax.set_ylim(0,y1)
ax.set_xlim(min(frq), max(frq))
ax.set_ylim(0, 3000)
ax.text(lines[0], 2000, '$ms_p - E_y$', ha='center', va='top',
        color='r', backgroundcolor='w',rotation='vertical')
ax.text(lines[1], 2000, '$ms_m - E_x$', ha='center', va='top',
        color='r', backgroundcolor='w',rotation='vertical')
ax.text(lines[2], 2000, '$ms_0 - E_y$', ha='center', va='top',
        color='r', backgroundcolor='w',rotation='vertical')
ax.text(lines[4], 2000, '$ms_0 - E_x$', ha='center', va='top',
        color='r', backgroundcolor='w',rotation='vertical')
ax.text(lines[3], 2000, '$ms_p - A_1$', ha='center', va='top',
        color='r', backgroundcolor='w',rotation='vertical')
ax.text(lines[5], 2000, '$ms_m - A_1$', ha='center', va='top',
        color='r', backgroundcolor='w',rotation='vertical')
ax.text(lines[6], 2000, '$ms_p - A_2$', ha='center', va='top',
        color='r', backgroundcolor='w',rotation='vertical')
ax.text(lines[7], 2000, '$ms_m - A_2$', ha='center', va='top',
        color='r', backgroundcolor='w',rotation='vertical')




