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
fname ='NV_C011_Raw_Data_ExtraB_Field' + datetime.datetime.now().strftime("%Y%m%d_%H%M")
filename = os.path.join(dir,fname)

#Phyiscal Model paramters
Carbon_Conc = 0.011 #1,1% C13
B_Fields = np.arange(1750,1800,150) #B_Field in Gauss
#Simulation Parameters
#N_NV = 1000#Number of NV Centres to generate
#Gridsize = 25
A_Max = 200e3 #Max Hyperfine interaction strength in Hz
A_Min = 1e3 # Min Hyperfine interaction strength !Allowed to be lower at lower Conc
tau_max= 1e-4 # max resonance time to look for in s !Allowed to be longer if lower Conc
tau_step = 1e-9 #timestep in s , was 2ns nu 1ns beste wat we experimenteel kunnen doen
max_gate_time = 20e-3# choose high because possible to correct in post proc

#Loading Old Dataset
loadname ='simulation_data_werkbespreking/NV_C13_0.0011_Raw_Data_20131204_1522'
print 'loading data from file: ' +str(loadname)
file = open(loadname,'r')
B_Fields = pickle.load(file)
NV_List = pickle.load(file)
Raw_Data = pickle.load(file)

N_NV =np.shape(NV_List)[0]
B_Fields = np.arange(1750,1800,150) #B_Field in Gauss

#Saving settings to file
Names  = np.array(['Settings of simulation' ,'Sim started', 'Carbon_Conc', r'N_NV    ',r'A_Max    ',r'A_Min    ',\
    'tau_max','tau_step','max_gate_time'])
Floats  = np.array(['         ', datetime.datetime.now().strftime("%d-%m-%Y %H:%M"),Carbon_Conc, N_NV,A_Max,A_Min,\
    tau_max,tau_step,max_gate_time])
DAT =  np.column_stack((Names, Floats))
np.savetxt(filename + '.set', DAT ,delimiter=' \t ',fmt="%s")


#Variables to initialise
Raw_Data = np.empty((N_NV,np.size(B_Fields)),dtype=np.object)

# #Start of the actual program
# t_start = time.time()
# NV_List = Generate_NV(Carbon_Conc,N_NV,Gridsize)
# print 'NV_Centres Generated in ' + str(time.time()-t_start) +'s'

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
            Raw_Carbons = np.empty((np.shape(indices)[0],2),dtype = np.object)
            Raw_Carbons[:,0] = indices
            for idi, ind  in enumerate(indices): #loop over individual carbons
                HyperfineStr = NV[ind,:]
                res = Find_Res(HyperfineStr,B_Field,tau_max,tau_step,max_gate_time) #Find the first resonances
                if np.size(res)!=0:
                    res = res[np.where(res[:,1]!=0)] #Reject False positives
                    if np.size(res)!=0: #second check for empy matrix after rejection of false positives
                        F_of_Res = np.zeros(np.shape(res)[0])
                        for k in range(np.shape(res)[0]):
                             F_of_Res[k] = Char_Gate(NV,res[k],B_Field)
                    else:
                        F_of_Res =[]
                        res = [0,0]
                else:
                        F_of_Res =[]
                        res = [0,0]
                if np.size(F_of_Res)==0:
                    Raw_Res = np.hstack((0,res))
                else:
                    Raw_Res = np.hstack((F_of_Res.reshape(np.shape(F_of_Res)[0],1),res))
                Raw_Carbons[idi,1] = Raw_Res
            Raw_Data[n,idb]=Raw_Carbons
            print 'characterizing resonances for B-field took ' +str(time.time()-t_Bstart)
        print 'Characterizing this NV center took' +str(time.time()-t_find_C)
#Saving the data
file= open(filename , 'w') #Shelving would be cleaner but so far this works fine
pickle.dump(B_Fields,file)
pickle.dump(NV_List,file)
pickle.dump(Raw_Data,file)
file.close()
print 'simulation took ' +str(time.time()- t_find_C) +r'/n'


