
import matplotlib
#mpl.rcParams['text.usetex']=True
#mpl.rcParams['text.latex.unicode']=True
execfile('D:/measuring/analysis/scripts/setup_analysis.py')
from matplotlib import rc, cm
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
from analysis.lib.magnetometry import plotting_tools as pt

reload(magnetometry)
reload(pt)
load_data=False


def compare_cappellaro_G3 ():
    tag = 'simulated_adaptive_magnetometry'
    t_stamps = []
    labels = []
    for i in np.arange(3):
        t_stamps.append(tag+'_N=10G=3F='+str(i))
        labels.append('G=3, F='+str(i))

    colours = ['b', 'r', 'k']
    pt.compare_multiple_plots (timestamps=t_stamps, labels=labels, title = 'compare protocols', colours=colours)


compare_cappellaro_G3()