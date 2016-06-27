from analysis.lib.m2 import m2
reload(m2)
import numpy as np
from analysis.lib.tools import toolbox as tb
from matplotlib import pyplot as plt
import os

class DisplayScan(m2.M2Analysis):

    def get_data(self):
        self.zfocus = self.f['instrument_settings']['master_of_space'].attrs['z']

        self.xvalues = self.f['x'].value
        self.yvalues = self.f['y'].value
        self.countrates = self.f['countrate'].value

    def plot_data(self,title,save=True):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        colorname='afmhot'
        colors=ax.pcolor(self.xvalues,self.yvalues,self.countrates, vmin=np.amin(self.countrates),vmax=np.amax(self.countrates),cmap=colorname)
        ax.set_xlim([np.amin(self.xvalues),np.amax(self.xvalues)])
        ax.set_ylim([np.amin(self.yvalues),np.amax(self.yvalues)])
        ax.set_title('z = '+str(self.zfocus)+'    '+title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax3 = fig.add_subplot(1,1,1)
        # ax3.axis('off')
        fig.colorbar(colors,ax=ax)
        plt.show()

        if save:
            try:
                fig.savefig(
                    os.path.join(self.folder,colorname+'_scan2d.png'))
            except:
                print 'Figure has not been saved.'


def display_scan(older_than = None, nr_plots =1):
    """function that calls DisplayScan
    PARAMETERS:
    older_than - timestamp that identifies the first data (default = None; gives latest)
    nr_plot - give the total number of plots you want to show. default = 1."""

    for j in np.arange(nr_plots):
        timestamp,folder = tb.latest_data(contains='scan2d',older_than = older_than, return_timestamp=True)
        older_than=timestamp

        title = folder
        a = DisplayScan(folder)
        a.get_data()
        a.plot_data(title)
