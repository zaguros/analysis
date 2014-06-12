import os
import numpy as np
from matplotlib import pyplot as plt
from measurement.lib.tools import toolbox
from analysis.lib.nv import nvlevels
from analysis.lib.m2.ssro import ssro

reps = 5000

#gate_V = np.array([0,-100,-250, -500, -700,-900])*90./2000.
#E_y = np.array([47.11, 47.58, 48.3,  49.8,  50.74, 52.09])
#E_x = np.array([53.53, 53.18, 53.27, 52.82, 52.97, 52.97])

#ssro_data_E_y = np.array(['20140513102849', '20140513105022', '20140513105541', '20140513110828', '20140513111136', '20140513114023'])
gate_V = np.array([-700,-800,-900, -1000, -1100,-1200, -1300, -1500, -1700,-2000])
gate_V= gate_V*90./2000.
E_y =  	 np.array([51.08, 51.53, 51.96,  52.33,  52.61, 53.03, 53.36, 53.58, 53.5,54.84])
E_x =    np.array([53.75, 53.7,  53.73,  53.79,  54.15, 54.57, 55.04, 56.01,57.04,59.09])


ssro_data_E_y = np.array(['20140513124215', '20140513125308', '20140513125914', '20140513130340', '20140513130712', '20140513131304','20140513131824', '20140513132605','20140513134100','20140513134217'])

cpsh = []
fid = []
fid_err = []

for i,d in enumerate(ssro_data_E_y):
	folder = toolbox.data_from_time(d)
	a = ssro.SSROAnalysis(folder)
	a.get_run('ms0')
	cpsh_run = np.sum(a.ro_counts[:,0:-1], axis=1)
	cpsh.append(sum(cpsh_run)/float(reps))
	#f = h5py.File(folder, 'r')
	#cts_m0 = 
    #times = f['fidelity/time'].value


fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(gate_V, cpsh, 'bo', ms = 8)
ax.set_xlabel('strain splitting [GHz]')
ax.set_ylabel('counts per shot $E_y$')