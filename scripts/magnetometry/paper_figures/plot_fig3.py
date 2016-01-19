from analysis.lib.math import statistics as stat
from matplotlib import rc

rc('text', usetex=False)
execfile('D:\\machielblok/Desktop/PhD/qtlab/analysis/scripts/setup_analysis.py')
from matplotlib import pyplot as plt
reload(plt)
from analysis.scripts.magnetometry import adaptive_magnetometry_simulations as asim
reload(asim)
# These are the "Tableau 20" colors as RGB.  
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  


# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.  
for i in range(len(tableau20)):  
    r, g, b = tableau20[i]  
    tableau20[i] = (r / 255., g / 255., b / 255.)  
markers=['v','^']
colors=[6,0]
mgnt_G5F7_N14_data=asim.analyze_saved_simulations('20150325_142959',error_bars=True)
mgnt_G5F7_N14_sim=asim.analyze_saved_simulations('20141119_094103')

fig = plt.figure(figsize=(3.5,2))
p = fig.add_subplot(111)
d=mgnt_G5F7_N14_data
sim=mgnt_G5F7_N14_sim

def calc_error(result_dict):

    
    err_VH=[]
    start=0
    stop=0
    for dp in np.arange(7*4): # now very specific on nr of periods and points per period. make more general!

        nr_discarded_CR=result_dict['nr_discarded_elements'][dp]   
        nr_success_CR=101-nr_discarded_CR
        start=stop
        stop=start+nr_success_CR

        bs = stat.BootStrap (n_boots = 1000)
        estim_phases = np.array(result_dict['estimated_phase_values'][start:stop])
        bs.set_y (y=estim_phases)
        bs.run_bootstrap_holevo()
        err_VH.append(bs.errH_bootstrap)
    return err_VH
j=0
for i in [2,4]:
    B=d.results_dict[str(i)]['B_field']*1e-6
    VH=d.results_dict[str(i)]['msqe']
    err_VH=calc_error(d.results_dict[str(i)])
    B_sim=sim.results_dict[str(i)]['B_field']*1e-6
    VH_sim=sim.results_dict[str(i)]['msqe']
    
    lbl='N = '+str(i)

    
    plt.plot(B_sim,VH_sim,'-',markersize=12,color=tableau20[colors[j]+1],linewidth=2)
    (_, caps, _) = plt.errorbar(B,VH,yerr=err_VH,fmt=markers[j],elinewidth=1,label=lbl,color=tableau20[colors[j]])

    for cap in caps:
        cap.set_markeredgewidth(1)
    j+=1

plt.xlabel( '$\delta f_b$ (MHz)')
plt.ylabel('$V_{H}$ (rad$^2$)')
plt.tick_params(axis='x')
plt.tick_params(axis='y')
plt.yscale('log')
plt.xlim([-25,25])
#plt.ylim([1e-5,1e4])
plt.legend(prop={'size':18},loc=4)
fig.savefig(r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\analyzed data\Bfield_dep_G5F7_data_and_sim.pdf', bbox_inches='tight')