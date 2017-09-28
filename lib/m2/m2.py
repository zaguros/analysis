### imports
import sys, os, time
import logging
import numpy as np
import h5py
from matplotlib import pyplot as plt
from analysis.lib.tools import toolbox


class M2Analysis:

    plot_format = 'png'

    def __init__(self, folder, **kw):
        #print 'analyzing in', folder
        self.folder = folder
        self.load_hdf5data(**kw)
        

    def load_hdf5data(self, **kw):
        self.h5filepath = toolbox.measurement_filename(self.folder)
        h5mode=kw.pop('hdf5_mode', 'r')
        self.f = h5py.File(self.h5filepath,h5mode)
        for k in self.f.keys():
            if type(self.f[k])==h5py.Group and k in os.path.split(self.h5filepath)[1]: # PH added this check because sometimes additional files added to hdf5 file (06/16)
                self.name = k
                self.g = self.f[self.name]      
                break

        self.measurementstring = os.path.split(self.folder)[1]
        self.timestamp = os.path.split(os.path.split(self.folder)[0])[1] \
                + '/' + self.measurementstring[:6]
        self.measurementstring = self.measurementstring[7:]        
        self.default_plot_title = self.timestamp+'\n'+self.measurementstring      

    def finish(self):
        self.f.close()

    def adwingrp(self, name=''):

        if name != '':
            adwingrpname = name
        else:
            if len(self.g.keys()) == 1:
                adwingrpname = self.g.keys()[0]
            else:
                logging.error("More than one measurement. Please give a name")
                return False
        return self.g[adwingrpname]

    def analysis_h5data(self, name='analysis'):
        if not os.path.exists(os.path.join(self.folder, name+'.hdf5')):
            mode = 'w'    
        else:
            mode = 'r+'            
            
        return h5py.File(os.path.join(self.folder, name+'.hdf5'), mode)

    def default_fig(self, **kw):
        figsize = kw.pop('figsize', (4,4))

        return plt.figure(figsize=figsize, **kw)

    def default_ax(self, fig=None, *arg, **kw):
        if fig == None:
            fig = self.default_fig(*arg, **kw)

        ax = fig.add_subplot(111)
        ax.set_title(self.timestamp+'\n'+self.measurementstring)

        return ax

    def save_fig_incremental_filename(self,fig,savename):
            i=1
            sfn=os.path.join(self.folder, savename+'-{}.png'.format(i))
            while(os.path.exists(sfn)):
                sfn=os.path.splitext(sfn)[0][:-2]+'-{}.png'.format(i)
                i=i+1
            fig.savefig(sfn,format='png')




