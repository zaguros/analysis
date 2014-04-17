#Master script Carbon Addresability Simulation
#Import packages
import numpy as np  #Numerical module of python, used to handle array type manipulations
import matplotlib.pyplot as plt  #used for plotting
import cPickle as pickle #used to store variables in a nice way
import os   #needed for file and folder management
import datetime #Needed to put timestamps on filenames
import sys
sys.path.insert(0, 'Functions_For_Sims') #Required because functions in different subfolder
import time #Used for timing purposes





t_start = time.time()

#Import custom functions
from Generate_NV import Generate_NV
gam_el = 1.760859 *10**11 #Gyromagnetic ratio rad s-1 T-1
gam_n = 67.262 *10**6 #rad s-1 T-1

#Carbon_Conc = 0.011 #1,1% C13
Carbon_Concs = [0.011,0.003,0.0011]
#Simulation Parameters
N_NV = 1000#Number of NV Centres to generate
Gridsizes = [35, 35, 35]
avgT2star_el = np.zeros(3)
avgT2star_carb = np.zeros(3)
T2star_carb_lst = []
for idmu, Carbon_Conc in enumerate(Carbon_Concs):
    Gridsize = Gridsizes[idmu]
    NV_List = Generate_NV(Carbon_Conc,N_NV,Gridsize)
    b_el = np.zeros(N_NV)
    b_Carb = np.zeros(N_NV)
    T2star_el = np.zeros(N_NV)
    T2star_carb = np.zeros(N_NV)

    for n,NV in enumerate(NV_List):
        b_el[n] = 0.5*np.sqrt(np.sum(np.dot(NV[:,0],NV[:,0])))*2*np.pi
        T2star_el[n] = np.sqrt(2)/b_el[n]
        #b_Carb[n] = gam_n/gam_el*0.5*np.sqrt(np.sum(np.dot(NV[:,0],NV[:,0])+np.sum(np.dot(NV[:,1],NV[:,1]))))*2*np.pi
        #orthogonal and parralel
        b_Carb[n] = gam_n/gam_el*0.5*np.sqrt(np.dot(NV[:,0],NV[:,0])+np.dot(NV[:,1],NV[:,1]))*2*np.pi
        #only parallel
        #b_Carb[n] = gam_n/gam_el*0.5*np.sqrt(np.dot(NV[:,0],NV[:,0]))*2*np.pi
        T2star_carb[n] =np.sqrt(2)/b_Carb[n]
    avgT2star_el[idmu] = np.mean(T2star_el)
    avgT2star_carb[idmu]=np.mean(T2star_carb)
    print 'average electron T2 star is : '
    print avgT2star_el[idmu]
    print 'average carbon-13 T2 star is : '
    print avgT2star_carb[idmu]
    print ''
    T2star_carb_lst.append(T2star_carb)


print 'simulation took ' +str(time.time()- t_start) +r'\n'

plt.subplots()
plt.title(r'Carbon $T_2^*$ for $\mu = 1.1\%$ ' )
plt.ylabel(r'normalized bin counts' )
plt.xlabel('$T_{2}^*$ [s]')
n, bins, patches = plt.hist(T2star_carb_lst[0], 100, normed=1,stacked = True,cumulative = False)
textstr = '$\mu = 0.011\% $ '
props = dict(boxstyle='round',facecolor = 'wheat', alpha=0)
    #ax.text(0.7, 0.5, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)

plt.subplots()
#plt.figure()
plt.title(r'Carbon $T_2^*$ for $\mu = 0.3\%$ ' )
plt.ylabel(r'normalized bin counts' )
plt.xlabel('$T_{2}^*$ [s]')
n, bins, patches = plt.hist(T2star_carb_lst[1], 100, normed=1,stacked = True,cumulative = False)
textstr = '$\mu = 0.011\% $ '
props = dict(boxstyle='round',facecolor = 'wheat', alpha=0)

plt.subplots()
#fig3,ax = plt.figure()
plt.title(r'Carbon $T_2^*$ for $\mu = 0.11\%$ ' )
plt.ylabel(r'normalized bin counts' )
plt.xlabel('$T_{2}^*$ [s]')
n, bins, patches = plt.hist(T2star_carb_lst[2], 100, normed=1,stacked = True,cumulative = False)
textstr = '$\mu = 0.011\% $ '
props = dict(boxstyle='round',facecolor = 'wheat', alpha=0)


## Inverse cumulative
plt.subplots()
plt.title(r'Carbon $T_2^*$ for $\mu = 1.1\%$ ' )
plt.ylabel(r'complementary cumulative distribution' )
plt.xlabel('$T_{2}^*$ [s]')
n, bins, patches = plt.hist(T2star_carb_lst[0], 100, normed=1,stacked = True,cumulative = -1)
textstr = '$\mu = 0.011\% $ '
props = dict(boxstyle='round',facecolor = 'wheat', alpha=0)
    #ax.text(0.7, 0.5, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)

plt.subplots()
#plt.figure()
plt.title(r'Carbon $T_2^*$ for $\mu = 0.3\%$ ' )
plt.ylabel(r'complementary cumulative distribution' )
plt.xlabel('$T_{2}^*$ [s]')
n, bins, patches = plt.hist(T2star_carb_lst[1], 100, normed=1,stacked = True,cumulative = -1)
textstr = '$\mu = 0.011\% $ '
props = dict(boxstyle='round',facecolor = 'wheat', alpha=0)

plt.subplots()
#fig3,ax = plt.figure()
plt.title(r'Carbon $T_2^*$ for $\mu = 0.11\%$ ' )
plt.ylabel(r'complementary cumulative distribution' )
plt.xlabel('$T_{2}^*$ [s]')
n, bins, patches = plt.hist(T2star_carb_lst[2], 100, normed=1,stacked = True,cumulative = -1)
textstr = '$\mu = 0.011\% $ '
props = dict(boxstyle='round',facecolor = 'wheat', alpha=0)


plt.show()

