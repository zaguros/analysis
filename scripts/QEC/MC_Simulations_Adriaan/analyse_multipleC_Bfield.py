import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os   #needed for file and folder management
import datetime #Needed to put timestamps on filenames
import sys
sys.path.insert(0, 'Functions_For_Sims') #Required because functions in different subfolder
import time #Used for timing purposes

filenames =      ['simulation_data_3month_report/NV_C13_0.011_Raw_Data_20140130_1543plot_data',
    'simulation_data_3month_report/NV_C13_0.003_Raw_Data_20140130_2353plot_data',
    'simulation_data_3month_report/NV_C13_0.0011_Raw_Data_20140201_1409plot_data']
t_start = time.time()

B_Fields =np.zeros(np.size(filenames),dtype=object)
Navg_B =np.zeros(np.size(filenames),dtype=object)
Navg_B_Err=np.zeros(np.size(filenames),dtype=object)

for idfn, filename in enumerate(filenames) :
    print 'loading data from file: ' +str(filename)
    file = open(filename,'r')
    B_Fields[idfn] = pickle.load(file)
    Navg_B[idfn] =pickle.load(file)
    Navg_B_Err[idfn] = pickle.load(file)

mu_ls = [1.1,0.3,0.11]
F_Min_lab = 0.90
T_Max_lab = [1.4,4.8,14]
B_lab = [700,400,150]






##########
figsize= (2.5,2)
fontsize = 8
linewidth = .75
markersize =2
plt.rc('font', size=8)
fig,ax = plt.subplots(figsize=figsize)

# fig,ax = plt.subplots(figsize=figsize)
plt.title(' ' )
plt.ylabel(r'Avg. No. Address. C-13' )
plt.xlabel('Magnetic field (G)')
for idfn, filename in enumerate(filenames):
    labstr = r'$\mu$ = ' +str(mu_ls[idfn]) + r'% , '+r' $T_{Max}$ = '+str(T_Max_lab[idfn])+'ms , ' +r'$F_{Min}$= ' +str(F_Min_lab)
    plt.plot(B_Fields[idfn],Navg_B[idfn],label =(labstr),marker= 'o' ,markersize=2.5)
plt.xlim([B_Fields[0][0]-100,B_Fields[0][-1]+100])
# plt.legend(loc= 'lower right')

# plt.show()


folder_a = '/Users/Adriaan/Documents'
savename = 'Simulations_avgN_vs_Bfield'
fig.savefig(os.path.join(folder_a, savename+'.pdf'),
                format='pdf',bbox_inches='tight')
print' Figure saved in %s' %folder_a
