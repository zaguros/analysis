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

fig,ax = plt.subplots(1)
plt.title(r'Average number of addressable $C^{13}$ in weakly coupled NV centres' )
plt.ylabel(r'$\bar{N}$' )
plt.xlabel('B-Field [Gauss]')
for idfn, filename in enumerate(filenames):
    labstr = r'$\mu$ = ' +str(mu_ls[idfn]) + r'% , '+r' $T_{Max}$ = '+str(T_Max_lab[idfn])+'ms , ' +r'$F_{Min}$= ' +str(F_Min_lab)
    plt.errorbar(B_Fields[idfn],Navg_B[idfn],yerr = Navg_B_Err[idfn],label =(labstr) )
plt.xlim([B_Fields[0][0]-100,B_Fields[0][-1]+100])
plt.legend(loc= 'lower right')

plt.show()