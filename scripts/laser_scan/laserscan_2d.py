'''
Updated script to analyze linescan data using our more modern scripts
2016 PH
'''

import numpy as np
import os
from analysis.lib.tools import toolbox, plot; 
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
from analysis.lib.fitting import fit, common
import copy as cp
from scipy.interpolate import interp1d

def plot_freq_against_index(contains = '', interp_res = 0.01, **kw):

	title = kw.pop('title',None)

	folder = toolbox.latest_data(contains)

	filename = toolbox.measurement_filename(folder, ext = 'dat')

	V,f,c,i,st = np.loadtxt(filename,unpack=True)

	Y=np.unique(i)
	X = np.arange(np.min(f),np.max(f),interp_res)
	Z=np.zeros((len(X), len(Y)))

	for j,y in enumerate(Y):
	    fltr=i==y #(dd[:,2]<(y+0.05)) & (dd[:,2]>(j-0.05)) #
	    xx=f[np.where(fltr)]
	    zz=c[np.where(fltr)]
	    interp_func=interp1d(xx,zz, bounds_error=False, fill_value=0.)
	    Z[:,j]=np.nan_to_num(interp_func(X))
	
	fig = plt.figure()
	ax=plt.subplot(111)

	ax.imshow(Z, aspect='auto', origin='lower',vmin=0,vmax = np.max(Z),extent=[min(Y),max(Y),X[0],X[-1]], cmap='binary')
	ax.set_xlabel('Scan nr')
	ax.set_ylabel('Laser frequency [GHz]')

	if title != None:
		plt.title(title + ' '+os.path.split(filename)[1])
	else:
		plt.title(os.path.split(filename)[1])

	save_and_close_plot(folder,name = 'Laserscan_2d')

def save_and_close_plot(f,name = 'Results'):
    plt.savefig(os.path.join(f,name + '.pdf'),format='pdf')
    plt.savefig(os.path.join(f,name + '.png'),format='png')
    plt.show()
    plt.close('all')


def get_tstamp_from_folder(folder):
    return folder[18:18+15]
