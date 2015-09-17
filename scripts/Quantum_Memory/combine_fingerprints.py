import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot
from analysis.lib.math import error



def get_data(name_contains):
    folder = toolbox.latest_data(name_contains)
    print folder
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='measurement0')
    a.get_electron_ROC()
    #ax = a.plot_results_vs_sweepparam(ret='ax', )

    t_shift = 0
    # t_start=0
    # t_0=10

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    return x,y,folder

ids=np.arange(6)+1
all_x=[]
all_y=[]
for nr in ids:
    x,y,folder= get_data('16'+str(nr))
    
    all_x+=list(x)
    all_y+=list(y)

f = plt.figure(figsize=(25,6))
ax = f.add_subplot(1,1,1)
norm_y=all_y/mean(all_y[0:5])
plt.plot(all_x,norm_y,'-o',color='RoyalBlue')
plt.ylabel(r'Signal normalized to 1',fontsize=25)
plt.xlabel('tau (us)',fontsize=25)
plt.tick_params(axis='x', labelsize=25)
plt.tick_params(axis='y', labelsize=25)
plt.ylim([0.5,1.05])
#plt.xlim([10,15])
#print len(all_x)
miny=np.where(norm_y==min(norm_y))
print all_x[miny[0]]
    #plt.savefig(os.path.join(folder, 'Combined_fingerprint.pdf'),
    #        format='pdf')
    #plt.savefig(os.path.join(folder, 'Combined_fingerprint.png'),
    #        format='png')

