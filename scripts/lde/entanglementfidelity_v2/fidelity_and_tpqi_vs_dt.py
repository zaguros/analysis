import os, sys, time
import pickle
import pprint
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import rcParams

from analysis.lib.fitting import fit
from analysis import config

import MLE

config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'
config.datadir = r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde'
savedir = os.path.join(config.outputdir,
                time.strftime('%Y%m%d')+'-ldefidelity')

rcParams['mathtext.default'] = 'regular'
rcParams['legend.numpoints'] = 1

### set stuff

tpqi_srcdir = '20121121-tpqi'
fidelity_srcdir = '20121121-ldefidelity'

state = 'psi2'
ZZidx = 0
XXidx = 1
XmXidx = 2

win = 140
ch0start = 641
ch1start = 670

### load the tpqi visibility data
### TODO


### load the fidelity data (best guess, i.e., sqrt-term=0)
fid = {}
f = np.load(os.path.join(config.outputdir,fidelity_srcdir,
    'fidelities_bestguess_dtslices.npz'))
for k in f.keys():
    fid[k] = f[k]
f.close()

corr = {}
f = np.load(os.path.join(config.outputdir,fidelity_srcdir,
    'correlations_bestguess_dtslices.npz'))
for k in f.keys():
    corr[k] = f[k]
f.close()

# get the desired raw correlations
_c = corr['raw'+state+'correlations']
winidx = argmin(abs(corr['winvals']-win))
ch0idx = argmin(abs(corr['ch0starts']-ch0start))
ch1idx = argmin(abs(corr['ch1starts']-ch1start))
c = _c[:,winidx,ch0idx,ch1idx,...]
dtvals = corr['dtvals']
edges = np.append(0,dtvals)
meandtvals = edges[:-1] + (edges[1:]-edges[:-1])/2.

# do an MLE for each dt bin
F = np.array([])
u_F = np.array([])
for i,dt in enumerate(dtvals):
    emily = MLE.Emily()
    emily.prob_pts = 101
    emily.ab_pts = 101
    emily.state = state
    emily.N_ZZ = c[i,...,ZZidx,:]
    emily.N_XX = c[i,...,XXidx,:]
    emily.N_XmX = c[i,...,XmXidx,:]
    emily.fac_exact = 0
    emily.likelihood()
    emily.analysis(savefn='MLE_win%d_ch0%d_ch1%d_dt%d' \
            % (win,ch0start,ch1start,dt))
    F = np.append(F, emily.F)
    u_F = np.append(u_F, emily.u_F)
    plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(np.append(0,dtvals), np.append(F[0],F), 'k',
        drawstyle='steps')
ax.errorbar(meandtvals, F, fmt='o', mec='k', mfc='w',
        yerr=u_F, ecolor='k')

ax.set_xlabel('dt')
ax.set_ylabel('Fidelity')

fig.savefig(os.path.join(savedir, 'F_vs_dtslice_win%d_ch0%d_ch1%d.png' \
        % (win,ch0start,ch1start)))


    










