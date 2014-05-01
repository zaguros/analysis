import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, ramsey
from analysis.lib.tools import plot
from analysis.lib.math import error

reload(toolbox)

timestamp = '192057'#None#'155611'#'154018'#'184802'#None#'20130530183305'#None #
t_ramsey=1/2.19e6
gamma=2.8e6
bs=[8]
dtime=20e-6
### script
if timestamp != None:
    folder = toolbox.latest_data(timestamp)
else:
    folder = toolbox.latest_data()

a = sequence.SequenceAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='ssro')
fig = a.default_fig(figsize=(24,4))
ax = a.default_ax(fig)
for j in arange(len(bs)):
    binsize=bs[j]

    a.sweep_pts=np.ones(a.reps/binsize)
    db=np.ones(a.reps/binsize)
    a.normalized_ssro=np.ones(a.reps/binsize)
    a.u_normalized_ssro=np.ones(a.reps/binsize)

    print binsize
    for i in arange(a.reps/binsize):

        a.normalized_ssro[i]=sum(a.ssro_results[i*binsize:(i+1)*binsize])/float(binsize)
        a.u_normalized_ssro[i]=(a.normalized_ssro[i]*(1.-a.normalized_ssro[i])/binsize)**0.5
        a.sweep_pts[i]=i*dtime*binsize
        a.get_electron_ROC()

        db[i]=1000*((math.acos((2*a.normalized_ssro[i])-1)/np.pi)-1)/(t_ramsey*gamma*2)
    a.sweep_name = 'Time (s)'

    a.plot_result_vs_sweepparam(ret='ax', name='ssro',ax=ax)
fig = a.default_fig(figsize=(6,4))
plt.clf()
plt.hist(db,bins=25,weights=np.ones(len(db))/len(db))
plt.xlabel('delta B [mG]')
plt.xlim([-250,0])
plt.ylabel('P')
plt.savefig(os.path.join(folder, 'electronramsey_hist.png'),
        format='png')
'''
ax.legend([bs[0],bs[1],bs[2]])

plt.savefig(os.path.join(folder, 'electronramsey_analysis.pdf'),
        format='pdf')
plt.savefig(os.path.join(folder, 'electronramsey_analysis.png'),
        format='png')
'''
### FFT
# p0_fft = np.fft.fft(a.p0.reshape(-1), n=32)
# frq = np.fft.fftfreq(32, d=(a.sweep_pts[1]-a.sweep_pts[0])/1e3)
#
# fig = a.default_fig()
# ax = a.default_ax(fig)
# ax.plot(frq, p0_fft, 'o-')
