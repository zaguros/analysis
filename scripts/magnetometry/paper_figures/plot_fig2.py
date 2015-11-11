execfile('D:\\machielblok/Desktop/PhD/qtlab/analysis/scripts/setup_analysis.py')
from analysis.lib.tools import toolbox as tb
from analysis.scripts.magnetometry import adaptive_magnetometry_analysis as ama
reload(ama)
import h5py
#ama.analyze_single_instance(label='det=2.0MHz_N = 1_M=(5, 7)',compare_to_simulations=False)
#ama.analyze_single_instance(label='det=2.0MHz_N = 3_M=(5, 7)',compare_to_simulations=False)
#ama.analyze_single_instance(label='det=2.0MHz_N = 5_M=(5, 7)',compare_to_simulations=False)
#ama.analyze_single_instance(label='det=2.0MHz_N = 7_M=(5, 7)',compare_to_simulations=False)
#ama.analyze_single_instance(label='det=2.0MHz_N = 9_M=(5, 7)',compare_to_simulations=False)
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
  
for i in range(len(tableau20)):  
    r, g, b = tableau20[i]  
    tableau20[i] = (r / 255., g / 255., b / 255.) 
fig = plt.figure(figsize=(7,3))    
folder=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\analyzed data'
fnames=['20141124_171302analysis_adaptive_magnetometry_single_inst_N=1G=5F=7_fid0=0.87.hdf5',
        '20141124_171303analysis_adaptive_magnetometry_single_inst_N=3G=5F=7_fid0=0.87.hdf5',
        '20141124_171305analysis_adaptive_magnetometry_single_inst_N=5G=5F=7_fid0=0.87.hdf5',
        '20141124_171311analysis_adaptive_magnetometry_single_inst_N=7G=5F=7_fid0=0.87.hdf5',
        '20141124_171322analysis_adaptive_magnetometry_single_inst_N=9G=5F=7_fid0=0.87.hdf5'
        ]
fi=0
plt.xlim([0,4])
plt.ylim([10e-5,0.1])
for i,N in enumerate([1,3,5,7,9]):
    l='det=2.0MHz_N = '+str(N)+'_M=(5, 7)'
    #folder=tb.latest_data(contains=l)
    f=fnames[fi]
    fi+=1
    file_name=os.path.join(folder,f)
    print file_name
    
    f = h5py.File(file_name,'r')
    #print f.attrs.keys()
    label='N = ' +str(N)
    data_grp = f['data']
    #print len(data_grp['beta_exp'])
    #print data_grp['beta_exp'].value*1e-6,data_grp['p_exp'].value
    plt.plot(data_grp['beta_exp'].value*1e-6,data_grp['p_exp'].value,label=label,linewidth=1.5,color=tableau20[i*2])
    plt.fill_between(data_grp['beta_exp'].value*1e-6,1e-50,data_grp['p_exp'].value,color='grey',alpha=0.1+0.01*i**3)
    #plt.yscale('log')
    #print np.sum(data_grp.attrs['p_exp'])
    f.close()

plt.xlabel('Magnetic field (MHz)')
plt.ylabel('Probability density')

plt.yscale('log')
plt.yticks([1e0,1e-10,1e-20,1e-30,1e-40,1e-50])
plt.xlim([-25,25])
plt.ylim([10e-50,1])
#plt.legend(prop={'size':12})

#plt.yticks([0.1,0.05,0])
#plt.xlim([0,4])
#plt.ylim([0,0.1])

#plt.gca().yaxis.set_major_locator(ticker.LogLocator(numticks=4))
#plt.gca().yaxis.set_major_locator(MaxNLocator(4))

#plt.gca().yaxis.set_major_locator(matplotlib.ticker.LogLocator(numticks=6))
#plt.xlim([1,3])
#fig.tight_layout()
plt.show()
filename=''
folder=r'M:\tnw\ns\qt\Diamond\Projects\Magnetometry with adaptive measurements\Data\analyzed data'
filename=r'det2MHz_Probability_distr.pdf'
sf=os.path.join(folder,filename)
fig.savefig(sf, bbox_inches='tight',bbox='tight')