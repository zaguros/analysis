
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

def plot_RT_contrast ():
	R = np.arange (50000)+1
	a0 = 0.031*R
	a1 = 0.021*R
	C = 1./(1+(2*(a0+a1))/((a0-a1)**2))**0.5

	matplotlib.rc ('xtick', labelsize=20)
	matplotlib.rc ('ytick', labelsize=20)
	f1 = plt.figure (figsize=(10,6))
	plt.semilogx (R, C, linewidth =3)
	plt.xlabel ('read-out repetitions', fontsize=20)
	plt.ylabel ('contrast C', fontsize=20)
	plt.savefig ('D:/measuring/fig_S5.svg')
	plt.show()

def compare_cappellaro_G2_fid1():
	t_stamps = []
	for i in [0,1,2,3,4,5]:
		t_stamps.append('modCapp_N=10G=2F='+str(i)+'_fid0=1.0')
	labels = ['F=0', 'F=1', 'F=2', 'F=3', 'F=4', 'F=5']
	colours = ['k', 'r', 'b', 'g', 'RoyalBlue', 'c']
	pt.compare_multiple_plots (timestamps=t_stamps, labels=labels, title = 'adaptive (G=2) - F0 = 1', colours=colours)

def return_t_stamps(protocol, N, G, fid, name, F_array=None):
	t_stamps = []
	if not(F_array):
		F_array = [0,1,2, 3, 4, 5]

	for i in F_array:
		t_stamps.append('_'+protocol+'_N='+str(N)+'G='+str(G)+'F='+str(i)+'_fid0='+fid+name)
	return t_stamps	

def compare_protocols_old (pr1, pr2, N, G, fid, name = ''):

	t_stamps = return_t_stamps (protocol = pr1, N=N, G=G, fid=fid, name=name, F_array = [3,4,5])
	dict_G2_adptv = pt.generate_data_dict(timestamps=t_stamps)
	plot_colors =  ['RoyalBlue', 'crimson', 'forestgreen', 'deepskyblue', 'gold','darkslategray', 'cornsilk', 'darkkhaki', 'darkviolet', 'limegreen']
	pt.compare_scalings (data_dict = dict_G2_adptv, title = pr1+' (G='+str(G)+') -- fid0 = '+fid+name, do_save=True, add_HL_plot = False, colours = plot_colors)
	pt.compare_variance_with_overhead (data_dict=dict_G2_adptv, title = pr1+' (G='+str(G)+') -- fid0 = '+fid+name, do_save = True, overhead = 300e-6, colours = plot_colors)
	t_stamps = return_t_stamps (protocol = pr2, N=N, G=G, fid=fid, name=name, F_array = [6,7,8])
	dict_G2_nonadptv = pt.generate_data_dict(timestamps=t_stamps)
	pt.compare_scalings (data_dict = dict_G2_nonadptv, title = pr2+' (G='+str(G)+') -- fid0 = '+name, do_save=True, add_HL_plot = False, colours = plot_colors)
	pt.compare_variance_with_overhead (data_dict=dict_G2_nonadptv, title = pr2+' (G='+str(G)+') -- fid0 = '+fid+name, do_save = True, overhead = 300e-6, colours = plot_colors)
	pt.compare_best_sensitivities ([dict_G2_adptv, dict_G2_nonadptv], title =  'G='+str(G)+' -- fid0 = '+fid+name, legend_array = [pr1, pr2], do_save=True, colours = plot_colors)
	pt.compare_scaling_fits ([dict_G2_adptv, dict_G2_nonadptv], title = 'G='+str(G)+' -- fid0 = '+fid+name, legend_array = [pr1, pr2], do_save=True, colours = plot_colors)

def compare_protocols (protocol_array, N, G, fid, name = '', F_array=None, add_HL_plot=False, y_lim=None):

	plot_colors =  ['RoyalBlue', 'crimson', 'forestgreen', 'deepskyblue', 'gold','darkslategray', 'cornsilk', 'darkkhaki', 'darkviolet', 'limegreen']
	dict_array = []
	if not(F_array):
		F_array = [0,1,2,3,4,5]

	for j,pr in enumerate(protocol_array):

		#if (pr=='nnAdptv'):
		#	ff = [0,1,2,3,4,5,6,7]
		#	print 'Non adaptive!!'
		#else:
		#	ff = [0,1,2,3,4,5]
		t_stamps = return_t_stamps (protocol = pr, N=N, G=G, fid=fid[j], name=name[j], F_array = F_array)
		print pr, t_stamps
		dict_gen= pt.generate_data_dict(timestamps=t_stamps)
		dict_array.append(dict_gen)
		pt.compare_scalings (data_dict = dict_gen, title = pr+' (G='+str(G)+') -- fid0 = '+fid[j]+name[j], do_save=True, add_HL_plot = add_HL_plot, colours = plot_colors, y_lim=y_lim)
		pt.compare_variance_with_overhead (data_dict=dict_gen, title = pr+' (G='+str(G)+') -- fid0 = '+fid[j]+name[j], do_save = False, overhead = 300e-6, colours = plot_colors)

	sens_dict = pt.compare_best_sensitivities (dict_array, title =  'G='+str(G)+' -- fid0 = '+fid[j]+name[j], legend_array = protocol_array, do_save=True, colours = plot_colors)
	pt.compare_scaling_fits (dict_array, title = 'G='+str(G)+' -- fid0 = '+fid[j]+name[j], legend_array = protocol_array, do_save=True, colours = plot_colors)
	return sens_dict

def compare_room_temperature ():
	pr1 = 'nnAdptv'
	pr2 = 'swarmOpt'
	sens_dict = compare_protocols (protocol_array = [pr1, pr2, pr1, pr2], N=10, G=5, fid=['0.88', '0.88', '1.0', '1.0'], name='_noT2')
	nn88 = sens_dict['0']
	sw88 = sens_dict['1']
	nn100 = sens_dict['2']
	sw100 = sens_dict['3']

	sw88_B = 1e12*sw88['sensitivity']/((2*np.pi*28*20)**2)
	nn88_B = 1e12*nn88['sensitivity']/((2*np.pi*28*20)**2)
	sw100_B = 1e12*sw100['sensitivity']/((2*np.pi*28*20)**2)
	nn100_B = 1e12*nn100['sensitivity']/((2*np.pi*28*20)**2)

	matplotlib.rc ('xtick', labelsize=25)
	matplotlib.rc ('ytick', labelsize=25)
	f1 = plt.figure(figsize = (10,7))
	plt.semilogy (sw88['F'], 1e3*sw88_B**0.5, 'RoyalBlue', linewidth =3, label = 'sw88')
	plt.semilogy (sw88['F'], 1e3*sw88_B**0.5, marker='^', color='RoyalBlue', markersize=12, label = 'a')
	plt.semilogy (nn88['F'], 1e3*nn88_B**0.5,'RoyalBlue', linestyle='--', linewidth =3, label = 'b')
	plt.semilogy (sw88['F'], 1e3*nn88_B**0.5, marker='v',color='RoyalBlue', markersize=12, label = 'c')
	plt.semilogy (sw88['F'], 1e3*(3600*sw88_B)**0.5, 'crimson', linewidth =3, label = 'd')
	plt.semilogy (sw88['F'], 1e3*(3600*sw88_B)**0.5, marker='h',color='crimson', markersize=10, label = 'e')
	plt.semilogy (nn88['F'], 1e3*(3600*nn88_B)**0.5,'crimson', linestyle='--', linewidth =3, label = 'r')
	plt.semilogy (nn88['F'], 1e3*(3600*nn88_B)**0.5, marker='*',color='crimson', markersize=15, label = 'f')
	plt.semilogy (sw100['F'], 1e3*(50000*sw100_B)**0.5, 'forestgreen', linewidth =3, label = 'g')
	plt.semilogy (sw100['F'], 1e3*(50000*sw100_B)**0.5, marker='o', color='forestgreen', markersize=10, label = 'h')
	plt.semilogy (nn100['F'], 1e3*(50000*nn100_B)**0.5,'forestgreen', linestyle='--', linewidth =3, label = 'i')
	plt.semilogy (nn100['F'], 1e3*(50000*nn100_B)**0.5,marker='s',color='forestgreen', markersize=10, label = 'l')
	plt.ylabel ('minimum sensitivity', fontsize=20)
	plt.xlabel ('F', fontsize=25)
	plt.legend()
	f1.savefig (r"D:\measuring\compare_sens_roomTemperature_nT.svg")
	f1.savefig (r"D:\measuring\compare_sens_roomTemperature_nT.pdf")
	plt.show()

plot_RT_contrast()

#compare_room_temperature()
#def plot_fig_S1 ():

#pr1 = 'nnAdptv'
#pr2 = 'swarmOpt'
#sens_dict = compare_protocols (protocol_array = [pr1, pr1], N=10, G=5, fid=['0.88', '0.88'], name=['_noT2_symmRO', '_noT2'])

#pr1 = 'modCapp'
#pr2 = 'nnAdptv'
#compare_protocols (protocol_array = [pr1, pr2], N=10, G=5, fid=['1.0', '1.0'], name='_incl_T2', add_HL_plot=True, y_lim=[1e-10, 1e-5])
#compare_protocols (protocol_array = [pr1, pr2], N=10, G=5, fid=['0.75', '0.75'], name='_incl_T2', add_HL_plot=True, y_lim=[1e-10, 1e-5])

#compare_capp_modCapp_supplInfo ()
#pr1 = 'modCapp'
#pr2 = 'nnAdptv'
#pr3 = 'swarmOpt'
#compare_protocols (pr1=pr1, pr2=pr2, N=10, G=5, fid='1.0', name='_incl_T2')
#compare_protocols (protocol_array = [pr1, pr2], N=10, G=5, fid=['0.65', '0.65'], name='')
#compare_protocols (pr1=pr1, pr2=pr2, N=10, G=3, fid='1.0', name='_incl_T2')
#compare_protocols (pr1=pr1, pr2=pr2, N=10, G=3, fid='0.75', name='_incl_T2')

#compare_protocols (N=10, G=2, fid='0.75')
#compare_protocols (N=10, G=3, fid='0.75')
#compare_protocols (N=10, G=4, fid='0.75')
#compare_protocols (N=10, G=5, fid='0.75')
#compare_protocols (N=10, G=10, fid='0.75')
#compare_protocols (N=10, G=10, fid='1.0')
#compare_protocols (N=10, G=2, fid='0.87')
#compare_protocols (N=10, G=3, fid='0.87')
#compare_protocols (N=10, G=4, fid='0.87')
#compare_protocols (N=10, G=5, fid='0.87')
#compare_protocols (N=10, G=10, fid='0.87')

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