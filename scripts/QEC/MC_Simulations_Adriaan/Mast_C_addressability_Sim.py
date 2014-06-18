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

#Import custom functions
from Generate_NV import Generate_NV
from NV_Reject import NV_Reject
from Find_Carb import Find_Carb
from addressable_C import addressable_C
from usefull_NV import usefull_NV
from Find_Res import Find_Res
from Char_Gate import Char_Gate

#Folder and file to save data to
dir = 'simulation_data'
if not os.path.isdir(dir): os.makedirs(dir)
fname ='NV_C011_F085_Amax_180k_tstep2ns_' + datetime.datetime.now().strftime("%Y%m%d_%H%M")
filename = os.path.join(dir,fname)

#Phyiscal Model paramters
Carbon_Conc = 0.011 #1,1% C13
B_Fields = np.arange(100,1800,100) #B_Field in Gauss
#Simulation Parameters
N_NV = 200#Number of NV Centres to generate
Gridsize = 25
A_Max = 180e3 #Max Hyperfine interaction strength in Hz
A_Min = 10e3 # Min Hyperfine interaction strength !Allowed to be lower at lower Conc
tau_max= 6e-6 # max resonance time to look for in s !Allowed to be longer if lower Conc
tau_step = 1e-9 #timestep in s , was 2ns nu 1ns beste wat we experimenteel kunnen doen
max_gate_time = 2e-3 #max gate time in s was 500e-6 kan hoger van tim (dubbel)!Allowed to be higher lower Conc
F_Min = 0.85 # Minimum Gate Fidelity
Gate_Min = 5 # Minimum number of gates

#Saving settings to file
Names  = np.array(['Settings of simulation' ,'Sim started', 'Carbon_Conc', r'N_NV    ','Gridsize',r'A_Max    ',r'A_Min    ',\
    'tau_max','tau_step','max_gate_time',r'F_min    ','Gate_Min'])
Floats  = np.array(['         ', datetime.datetime.now().strftime("%d-%m-%Y %H:%M"),Carbon_Conc, N_NV,Gridsize,A_Max,A_Min,\
    tau_max,tau_step,max_gate_time,F_Min,Gate_Min])
DAT =  np.column_stack((Names, Floats))
np.savetxt(filename + '.set', DAT ,delimiter=' \t ',fmt="%s")


#Variables to initialise
F=[]
Fidel =[]
Raw_Fidelities = np.zeros([N_NV,np.size(B_Fields)])
N_addressable_C = np.zeros([N_NV,np.size(B_Fields)])


#Start of the actual program
t_start = time.time()
NV_List = Generate_NV(Carbon_Conc,N_NV,Gridsize)
print 'NV_Centres Generated in ' + str(time.time()-t_start) +'s'

for n in range(N_NV):
    NV = NV_List[n]
    print 'addressing NV-center # ' +str(n)
    if NV_Reject(NV,A_Max):
        print 'center rejected'
    else:
        print 'center trough first test'
        t_find_C = time.time()
        indices = Find_Carb(NV,A_Min)
        #print 'found carbons in ' +str(time.time()-t_find_C)
        for idb, B_Field in enumerate(B_Fields): #loop over Field strengths
            t_Bstart = time.time()
            F_of_Carb = np.zeros(np.shape(indices)[0])
            for idi, ind  in enumerate(indices): #loop over individual carbons
                HyperfineStr = NV[ind,:]
                res = Find_Res(HyperfineStr,B_Field,tau_max,tau_step,max_gate_time) #Find the first resonances
                if np.size(res)!=0:
                    res = res[np.where(res[:,1]!=0)] #Reject False positives
                    if np.size(res)!=0: #second check for empy matrix after rejection of false positives
                        F_of_Res = np.zeros(np.shape(res)[0])
                        for k in range(np.shape(res)[0]):
                             F_of_Res[k] = Char_Gate(NV,res[k],B_Field)
                        F_of_Carb[idi]=np.max(F_of_Res)
                    else: F_of_Carb[idi]=0
                else: F_of_Carb[idi] = 0
            N_addressable_C[n,idb] = addressable_C(F_of_Carb,F_Min)
            print 'characterizing resonances for B-field took ' +str(time.time()-t_Bstart)
P_usefull = usefull_NV(N_addressable_C,Gate_Min)/N_NV

#Saving the data
file= open(filename , 'w') #Shelving would be cleaner but so far this works fine
pickle.dump(N_NV,file)
pickle.dump(B_Fields,file)
pickle.dump(P_usefull,file)
pickle.dump(N_addressable_C,file)
pickle.dump(NV_List,file)
file.close()
print 'simulation took ' +str(time.time()-t_start)
##plotting
N_rejected =np.sum( np.sum(N_addressable_C,axis=1)==0)
x = np.delete(N_addressable_C, np.where(np.sum(N_addressable_C,axis=1) ==0) , axis =0)
print str(N_rejected) +' of ' +str(N_NV) +'NV centres rejected'



plt.figure()
plt.ylabel('Probability that NV-Center is usefull')
plt.xlabel('B-Field [Gauss]')
plt.errorbar(B_Fields,P_usefull,yerr = 1/np.sqrt(N_NV))
plt.ylim([0,1])
plt.xlim([B_Fields[0]-100,B_Fields[-1]+100])
plt.figure()
plt.ylabel('Average number of addressable C')
plt.xlabel('B-Field [Gauss]')
plt.errorbar(B_Fields,np.mean(N_addressable_C,axis=0),yerr = np.std(N_addressable_C,axis=0))
plt.xlim([B_Fields[0]-100,B_Fields[-1]+100])

plt.figure()
plt.ylabel('Average number of addressable C in weakly coupled NVs' )
plt.xlabel('B-Field [Gauss]')
plt.errorbar(B_Fields,np.mean(x,axis=0),yerr = np.std(x,axis=0))
plt.xlim([B_Fields[0]-100,B_Fields[-1]+100])
plt.show()
