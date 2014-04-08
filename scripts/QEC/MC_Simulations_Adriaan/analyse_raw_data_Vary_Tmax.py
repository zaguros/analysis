import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os   #needed for file and folder management
import datetime #Needed to put timestamps on filenames
import sys
sys.path.insert(0, 'Functions_For_Sims') #Required because functions in different subfolder
import time #Used for timing purposes

t_start = time.time()
try:
    Raw_Data
except NameError:
    filename = 'simulation_data_3month_report/NV_C13_0.011_Raw_Data_20140130_1543'
    #filename = 'simulation_data_werkbespreking/NV_C13_0.003_Raw_Data_20131203_2309'

    print 'loading data from file: ' +str(filename)
    file = open(filename,'r')
    B_Fields = pickle.load(file)
    NV_List = pickle.load(file)
    Raw_Data = pickle.load(file)
    print 'Loading data took' +str(time.time()-t_start)
else:
    print 'data already loaded, reset environment when analyzing new data set'


#Load data


#Things we can tune in post processing of data
max_gate_time_ls =np.arange(0.05e-3,2.0e-3,.005e-3) # Minimum Gate Fidelity 2e-3 #max gate time in s was 500e-6 kan hoger van tim (dubbel)!Allowed to be higher lower Conc
F_Min = 0.90
B_Field = 700
idb =  np.where(B_Fields == B_Field)[0][0]

N = np.shape(Raw_Data)[0]
#data sturctures required for plotting
N_addressable_C = np.zeros([N,np.size(max_gate_time_ls)])

#We need the same loop structure as before to do the post selection this time

#Loop over NV's
for n in np.arange(N):
    #loop over BFields
    for idGt, max_gate_time in enumerate(max_gate_time_ls):
    #Loop over Carbons
        Raw_Data_NV_B = Raw_Data[n,idb]
        if Raw_Data_NV_B!=None: #Check if not empty
            indices = Raw_Data_NV_B[:,0]
            F = np.zeros(np.size(indices))
            for idi, ind  in enumerate(indices):
                res = Raw_Data_NV_B[:,1][idi]
                #Check for total gate time of resonances and select highest F
                if np.all(res == 0):
                    F[idi] = 0
                elif np.shape(res)==(1,3):
                    F[idi] =res[0,0]* (4*res[0,1]*res[0,2]<max_gate_time)
                else:
                    F[idi] =np.max(res[:,0]*(4*res[:,1]*res[:,2] < max_gate_time))
            N_addressable_C[n,idGt] = np.sum(F > F_Min)

N_rejected =np.sum( np.sum(N_addressable_C,axis=1)==0)
x = np.delete(N_addressable_C, np.where(np.sum(N_addressable_C,axis=1) ==0) , axis =0)
print str(N_rejected) +' of ' +str(N) +'NV centres rejected'

#Loop over resonances

Navg_T= np.mean(x,axis=0)
Navg_T_Err =np.std(x,axis=0)/np.sqrt(np.shape(x)[0])

fig,ax = plt.subplots(1)
plt.title(r'Average number of addressable $C^{13}$ in weakly coupled NV centres' )
plt.ylabel(r'$\bar{N}$' )
plt.xlabel('$T_{Max}$ [ms]')
plt.errorbar(max_gate_time_ls*1e3,Navg_T,yerr = Navg_T_Err)
plt.xlim([max_gate_time_ls[0]*1e3,max_gate_time_ls[-1]*1e3])
textstr = '$\mu = 1.1\% $ \n $B_z$ = 700 Gauss \n  $F_{Min}$ = 0.90'
props = dict(boxstyle='round',facecolor = 'wheat', alpha=0)
ax.text(0.7, 0.5, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)
plt.plot([1.68, 1.68], [0, 6], 'r--', lw=0.5)
plt.plot([1.4, 1.4], [0, 6], 'r-', lw=0.5)
plt.plot([1.12, 1.12], [0, 6], 'r--', lw=0.5)


plt.show()


