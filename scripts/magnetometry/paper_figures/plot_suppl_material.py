
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

def compare_capp_modCapp_supplInfo ():
    t_stamps = ['20150121_134642', '20150113_104934']
    labels = ['adptv phase update', 'modified adptv phase update']
    colours = ['RoyalBlue', 'crimson']
    pt.compare_multiple_plots (timestamps=t_stamps, labels=labels, title = 'G=3, F=0', colours=colours, do_save=True)
    t_stamps = ['20150121_172948', '20150113_164200']
    labels = ['adptv phase update', 'modified adptv phase update']
    colours = ['RoyalBlue', 'crimson']
    pt.compare_multiple_plots (timestamps=t_stamps, labels=labels, title = 'G=3, F=5', colours=colours, do_save=True)


def compare_cappellaro_G2_fid1():
	t_stamps = []
	for i in [0,1,2,3,4,5]:
		t_stamps.append('modCapp_N=10G=2F='+str(i)+'_fid0=1.0')
	labels = ['F=0', 'F=1', 'F=2', 'F=3', 'F=4', 'F=5']
	colours = ['k', 'r', 'b', 'g', 'RoyalBlue', 'c']
	pt.compare_multiple_plots (timestamps=t_stamps, labels=labels, title = 'adaptive (G=2) - F0 = 1', colours=colours)

def return_t_stamps(protocol, N, G, fid, name):
	t_stamps = []
	for i in [0,1,2, 3, 4, 5, 6, 7, 8]:
		t_stamps.append('_'+protocol+'_N='+str(N)+'G='+str(G)+'F='+str(i)+'_fid0='+fid+name)
	return t_stamps	

def compare_protocols (pr1, pr2, N, G, fid, name = ''):

	t_stamps = return_t_stamps (protocol = pr1, N=N, G=G, fid=fid, name=name)
	dict_G2_adptv = pt.generate_data_dict(timestamps=t_stamps)
	plot_colors =  ['RoyalBlue', 'crimson', 'forestgreen', 'deepskyblue', 'gold','darkslategray', 'cornsilk', 'darkkhaki', 'darkviolet', 'limegreen']
	pt.compare_scalings (data_dict = dict_G2_adptv, title = pr1+' (G='+str(G)+') -- fid0 = '+fid+name, do_save=True, add_HL_plot = False, colours = plot_colors)
	pt.compare_variance_with_overhead (data_dict=dict_G2_adptv, title = pr1+' (G='+str(G)+') -- fid0 = '+fid+name, do_save = True, overhead = 300e-6, colours = plot_colors)
	t_stamps = return_t_stamps (protocol = pr2, N=N, G=G, fid=fid, name=name)
	dict_G2_nonadptv = pt.generate_data_dict(timestamps=t_stamps)
	pt.compare_scalings (data_dict = dict_G2_nonadptv, title = pr2+' (G='+str(G)+') -- fid0 = '+name, do_save=True, add_HL_plot = False, colours = plot_colors)
	pt.compare_variance_with_overhead (data_dict=dict_G2_nonadptv, title = pr2+' (G='+str(G)+') -- fid0 = '+fid+name, do_save = True, overhead = 300e-6, colours = plot_colors)
	pt.compare_best_sensitivities ([dict_G2_adptv, dict_G2_nonadptv], title =  'G='+str(G)+' -- fid0 = '+fid+name, legend_array = [pr1, pr2], do_save=True, colours = plot_colors)
	pt.compare_scaling_fits ([dict_G2_adptv, dict_G2_nonadptv], title = 'G='+str(G)+' -- fid0 = '+fid+name, legend_array = [pr1, pr2], do_save=True, colours = plot_colors)

#compare_capp_modCapp_supplInfo ()
pr1 = 'modCapp'
pr2 = 'nnAdptv'
#compare_protocols (pr1=pr1, pr2=pr2, N=10, G=5, fid='1.0', name='_incl_T2')
compare_protocols (pr1=pr1, pr2=pr2, N=10, G=5, fid='0.75', name='_incl_T2')
#compare_protocols (pr1=pr1, pr2=pr2, N=10, G=3, fid='1.0', name='_incl_T2')
#compare_protocols (pr1=pr1, pr2=pr2, N=10, G=3, fid='0.75', name='_incl_T2')

#compare_protocols (N=10, G=2, fid='0.75')
#compare_protocols (N=10, G=3, fid='0.75')
#compare_protocols (N=10, G=4, fid='0.75')
#compare_protocols (N=10, G=5, fid='0.75')
#compare_protocols (N=10, G=10, fid='0.75')
#compare_protocols (N=10, G=10, fid='1.0')
compare_protocols (N=10, G=2, fid='0.87')
compare_protocols (N=10, G=3, fid='0.87')
compare_protocols (N=10, G=4, fid='0.87')
compare_protocols (N=10, G=5, fid='0.87')
compare_protocols (N=10, G=10, fid='0.87')

'''
G=3
t_stamps = ['20141215_152517', '20141215_152820', '20141215_153251', '20141215_153913']
dict_G3_adptv = pt.generate_data_dict(timestamps = t_stamps)
pt.compare_scalings (data_dict = dict_G3_adptv, title = 'Cappellaro (G='+str(G)+') -- fid0 = '+fid+name, do_save=True, add_HL_plot = True)
pt.compare_variance_with_overhead (data_dict=dict_G3_adptv, title = 'Cappellaro (G='+str(G)+') -- fid0 = '+fid+name, do_save = False, overhead = 300e-6)
pt.compare_sensitivity_repRate (data_dict_array = [dict_G3_adptv, dict_G3_adptv], legend_array=['test1', 'test2'], title='test', colours=None, do_save = True, overhead = 200e-6)
'''
#compare_scalings (data_dict_array, title, colours=None)
#compare_cappellaro_G2_fid1()
#compare_cappellaro_G3_fid087()
#compare_nonadaptive_G3_fid087()
#20141215_152517_simulated_adaptive_magnetometry_N=10G=3F=0_fid0=0.87