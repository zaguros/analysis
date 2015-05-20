#######################################
###   PLOT Fig 4                    ###
#######################################

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

def compare_G3_variableF_imperfRO():
    pt.compare_2plots (timestamp1= '20141215_152517', timestamp2='20141215_152604', title = 'G=3, F=0, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_152820', timestamp2='20141215_152939', title = 'G=3, F=1, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_153251', timestamp2='20141215_134529', title = 'G=3, F=2, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_134649', timestamp2='20141215_134741', title = 'G=3, F=3, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_134924', timestamp2='20141215_135030', title = 'G=3, F=4, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_135229', timestamp2='20141215_135352', title = 'G=3, F=5, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_135556', timestamp2='20141215_135723', title = 'G=3, F=6, RO_fid = 0.87')


def compare_G5_variableF_imperfRO():
    pt.compare_2plots (timestamp1= '20141215_134142', timestamp2='20141215_134205', title = 'G=5, F=0, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_134300', timestamp2='20141215_134334', title = 'G=5, F=1, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_134443', timestamp2='20141215_134529', title = 'G=5, F=2, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_134649', timestamp2='20141215_134741', title = 'G=5, F=3, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_134924', timestamp2='20141215_135030', title = 'G=5, F=4, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_135229', timestamp2='20141215_135352', title = 'G=5, F=5, RO_fid = 0.87')
    pt.compare_2plots (timestamp1= '20141215_135556', timestamp2='20141215_135723', title = 'G=5, F=6, RO_fid = 0.87')

def compare_G5_variableF_perfRO():
    pt.compare_2plots (timestamp1= '20141215_140947', timestamp2='20141215_141011', title = 'G=5, F=0, RO_fid = 1.00')
    pt.compare_2plots (timestamp1= '20141215_141102', timestamp2='20141215_141139', title = 'G=5, F=1, RO_fid = 1.00')
    pt.compare_2plots (timestamp1= '20141215_141245', timestamp2='20141215_141333', title = 'G=5, F=2, RO_fid = 1.00')
    pt.compare_2plots (timestamp1= '20141215_141503', timestamp2='20141215_141601', title = 'G=5, F=3, RO_fid = 1.00')
    pt.compare_2plots (timestamp1= '20141215_141745', timestamp2='20141215_141851', title = 'G=5, F=4, RO_fid = 1.00')
    pt.compare_2plots (timestamp1= '20141215_142057', timestamp2='20141215_142214', title = 'G=5, F=5, RO_fid = 1.00')

def compare_cappellaro_varG ():
    t_stamps = ['20141215_143011','20141215_143043','20141215_143125','20141215_143221','20141215_143328','20141215_143447']
    labels = ['G=0', 'G=1', 'G=2', 'G=3', 'G=4', 'G=5']
    pt.compare_multiple_plots (timestamps=t_stamps, labels=labels, title = 'Cappellaro protocol (F=0, RO_fid=1)')

def compare_swarm_opt ():
    t_stamps = ['20141215_161225','20141218_155251','20141216_035209']
    labels = ['G=5, F=2 capp', 'G=5, F=2 swarm', 'G=5, F=7 non_adptv']
    colours = ['b', 'r', 'k']
    pt.compare_multiple_plots (timestamps=t_stamps, labels=labels, title = 'compare protocols', colours=colours)



#compare_G5_variableF_imperfRO()
#compare_G5_variableF_perfRO()
#compare_cappellaro_varG()
compare_swarm_opt()
