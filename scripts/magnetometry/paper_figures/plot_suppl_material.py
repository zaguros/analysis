
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


def compare_cappellaro_G3_fid087 ():
    tag = 'simulated_adaptive_magnetometry'
    t_stamps = ['20141215_152517', '20141215_152820', '20141215_153251', '20141215_153913', '20141215_154728', '20141215_155704']
    labels = ['F=0', 'F=1', 'F=2', 'F=3', 'F=4', 'F=5']
    colours = ['k', 'r', 'b', 'g', 'RoyalBlue', 'c']
    pt.compare_multiple_plots (timestamps=t_stamps, labels=labels, title = 'adaptive phase update (G=3) - F0 = 0.87', colours=colours)

def compare_nonadaptive_G3_fid087 ():
    t_stamps = ['20141215_152604', '20141215_152939', '20141215_153454', '20141215_154152', '20141215_155054', '20141215_160056']
    labels = ['F=0', 'F=1', 'F=2', 'F=3', 'F=4', 'F=5']
    colours = ['k', 'r', 'b', 'g', 'RoyalBlue', 'c']
    pt.compare_multiple_plots (timestamps=t_stamps, labels=labels, title = 'non adaptive (G=3) - F0 = 0.87', colours=colours)

t_stamps = ['20141215_152517', '20141215_152820', '20141215_153251', '20141215_153913', '20141215_154728', '20141215_155704']
dict_G3_adptv = pt.generate_data_dict(timestamps = t_stamps)
t_stamps = ['20141215_152604', '20141215_152939', '20141215_153454', '20141215_154152', '20141215_155054', '20141215_160056']
dict_G3_nonAdptv = pt.generate_data_dict(timestamps = t_stamps)
pt.compare_best_sensitivities ([dict_G3_adptv, dict_G3_nonAdptv], title = 'G=3 -- fid_0 = 0.87', legend_array = ['adptv', 'non adptv'])
pt.compare_scaling_fits ([dict_G3_adptv, dict_G3_nonAdptv], title = 'G=3 -- fid_0 = 0.87', legend_array = ['adptv', 'non adptv'])



#compare_cappellaro_G3_fid087()
#compare_nonadaptive_G3_fid087()
#20141215_152517_simulated_adaptive_magnetometry_N=10G=3F=0_fid0=0.87