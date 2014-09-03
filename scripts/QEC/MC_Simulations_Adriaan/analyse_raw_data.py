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
    #filename = 'simulation_data_3month_report/NV_C13_0.003_Raw_Data_20140130_2353'
    #filename = 'simulation_data_3month_report/NV_C13_0.0011_Raw_Data_20140201_1409'
    print 'loading data from file: ' +str(filename)
    file = open(filename,'rb')
    B_Fields = pickle.load(file)
    NV_List = pickle.load(file)
    Raw_Data = pickle.load(file)
    print 'Loading data took' +str(time.time()-t_start)
else:
    print 'data already loaded, reset environment when analyzing new data set'


#Load data


#Things we can tune in post processing of data
max_gate_time = 1.4e-3
#max_gate_time = 4.7e-3
#max_gate_time = 14e-3 #max gate time in s was 500e-6 kan hoger van tim (dubbel)!Allowed to be higher lower Conc

F_Min = 0.90 # Minimum Gate Fidelity


N = np.shape(Raw_Data)[0]
#data sturctures required for plotting
N_addressable_C = np.zeros([N,np.size(B_Fields)])

#We need the same loop structure as before to do the post selection this time

#Loop over NV's
for n in np.arange(N):
    #loop over BFields
    for idb, B_Field in enumerate(B_Fields):
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
            N_addressable_C[n,idb] = np.sum(F > F_Min)

N_rejected =np.sum( np.sum(N_addressable_C,axis=1)==0)
x = np.delete(N_addressable_C, np.where(np.sum(N_addressable_C,axis=1) ==0) , axis =0)
print str(N_rejected) +' of ' +str(N) +'NV centres rejected'

#Loop over resonances
Navg_B = np.mean(x,axis=0)
Navg_B_Err = np.std(x,axis=0)/np.sqrt(np.shape(x)[0])

# fig,ax = plt.subplots(1)
# plt.title(r'Average number of addressable $C^{13}$ in weakly coupled NV centres' )
# plt.ylabel(r'$\bar{N}$' )
# plt.xlabel('B-Field [Gauss]')
# plt.errorbar(B_Fields,Navg_B,yerr = Navg_B_Err)
# plt.xlim([B_Fields[0]-100,B_Fields[-1]+100])
# textstr = '$\mu = 1.1\% $ \n $F_{Min}$ = 0.90  \n  $T_{max}$ = 1.4 ms'
# props = dict(boxstyle='round',facecolor = 'wheat', alpha=0)
# ax.text(0.7, 0.5, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)


# ## 3D Bar Histogram
nbins = 14
hist_data=np.zeros((np.shape(B_Fields)[0],nbins))
for idb, B_Field in enumerate(B_Fields):
    hist, bin_edges=  np.histogram(x[:,idb],bins = nbins, range = (-.5,nbins-.5),density = True)
    hist_data[idb] = hist
data =np.transpose(hist_data)#np.array([
column_names = B_Fields
row_names =np.arange(0,15,1)


# fig = plt.figure()
# ax = Axes3D(fig)
# lx= len(data[0])            # Work out matrix dimensions
# ly= len(data[:,0])
# xpos = np.arange(0,lx,1)    # Set up a mesh of positions
# ypos = np.arange(0,ly,1)
# xpos, ypos = np.meshgrid(xpos+0.25, ypos+0.25)
# xpos = xpos.flatten()   # Convert positions to 1D array
# ypos = ypos.flatten()
# zpos = np.zeros(lx*ly)
# dx = 0.5 * np.ones_like(zpos)
# dy = dx.copy()
# dz = data.flatten()

# ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color='b',alpha =0.5)
# ticksx = np.arange(0.5, np.size(B_Fields), 1)
# plt.xticks(ticksx, column_names)
# ticksy = np.arange(0.5, nbins, 1)
# plt.yticks(ticksy, row_names)
# ax.set_xlabel('B-Field [Gauss]')
# ax.set_ylabel('Number of Gates $N$')
# ax.set_zlabel('P $(n=N)$')
# plt.title(r'Addressable $C^{13}$ in weakly coupled NV centres' )


figsize= (2.5,2)
fontsize = 8
linewidth = .75
markersize =2
plt.rc('font', size=8)

## Histogram
fig,ax = plt.subplots(figsize=figsize)
plt.bar(bin_edges[:-1],hist_data[6])  #hist_data[nB] selects the B_field
# textstr = '$\mu = 1.1\% $ \n $ F_{min} $= 0.90\n $B_z$ = 700 Gauss \n  $T_{max}$ = 1.4ms'
# props = dict(boxstyle='round',facecolor = 'wheat', alpha=0)
# ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)
ax.set_xlabel(r'No. Address. C-13' )
ax.set_ylabel('Probability')

ax.set_xlim(-0.5,12.5)
#plt.title(r'Addressable $C^{13}$ in weakly coupled NV centres' )

# # ## Cumulative Histogram
# fig,ax = plt.subplots(1)
# cumhist = np.cumsum(hist_data[6][::-1])[::-1]
# plt.bar(bin_edges[:-1],cumhist)
# textstr = '$\mu = 1.1\% $ \n $ F_{min} $= 0.90\n $B_z$ = 700 Gauss \n  $T_{max}$ = 1.4ms'
# props = dict(boxstyle='round',facecolor = 'wheat', alpha=0)
# ax.text(0.7, 0.95, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)
# ax.set_xlabel('$n$')
# ax.set_ylabel('P $(n\geq N)$')
# #plt.title(r' $C^{13}$ in weakly coupled NV centres' )

#execfile('analyse_raw_data_Vary_Fidel.py')
#execfile('analyse_raw_data_Vary_Tmax.py')


# plotfilename = filename+ 'plot_data'
# file= open( plotfilename, 'w') #Shelving would be cleaner but so far this works fine
# pickle.dump(B_Fields,file)
# pickle.dump(Navg_B,file)
# pickle.dump(Navg_B_Err,file)
# #pickle.dump(F_Min_ls,file)
# #pickle.dump(Navg_F,file)
# #pickle.dump(Navg_F_Err,file)
# #pickle.dump(max_gate_time_ls,file)
# #pickle.dump(Navg_T,file)
#pickle.dump(Navg_T_Err,file)
#file.close()

# plt.show()



folder_a = '/Users/Adriaan/Documents'
savename = 'Simulations_Histogram_vs_Bfield'
fig.savefig(os.path.join(folder_a, savename+'.pdf'),
                format='pdf',bbox_inches='tight')
print' Figure saved in %s' %folder_a

