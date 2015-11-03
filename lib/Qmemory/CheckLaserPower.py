"""
Contains a function to check the maximum laser power after AOM calibration
Loops back in time over several recent AOM calibrations. Plots a comparison of the maximum powers.

NK 2015
"""
import datetime
import numpy as np
import os,h5py
from analysis.lib.tools import toolbox; reload(toolbox)
from matplotlib import pyplot as plt
from analysis.lib.tools import plot

def get_tstamp_from_folder(folder):
    return folder[18:18+15]

def dat_file_from_folder(folder,laser,controller):
    return folder[-(len(laser)+len(controller)+1+6+17):]+'.dat'

def convert_tstamp_to_minutes(tstamp):
    ### there are probably better ways to do this....

    d = int(tstamp[6:8])*60*24
    h = int(tstamp[8:10])*60
    mi = int(tstamp[10:12])
    s = float(tstamp[12:14])/60.
    # print y,m,d,h,mi,s
    return d+h+mi+s


def CheckPower(laser,controller,nr_of_files = 5,older_than = None):
    
    x,y = [],[]

    ### data acquisition

    for i in range(nr_of_files):
        f = toolbox.latest_data(contains = laser+'_'+controller,older_than = older_than)
        older_than = get_tstamp_from_folder(f)

        ## append time stamps to x
        x.append(convert_tstamp_to_minutes(older_than[:8]+older_than[9:]))
        
        if i == 0:
            reference = x[0]

        with open(f+'\\'+dat_file_from_folder(f,laser,controller),'r') as data:
            max_power = 0
            i = 0
            for line in data:
                if i < 15:
                    pass
                else:
                    if float(line[19:-19]) > max_power:
                        max_power = float(line[19:-19])

                i +=1


        ### single out the save power in the AOM calibration file    
        y.append(max_power)


    ### plotting
    x = np.array(x)-x[0]

    fig = plt.figure()
    ax = plt.subplot()
    
    plt.errorbar(x,y,fmt='o-')

    plt.xlabel('time (min)')
    plt.ylabel('max power (uW)')

    plt.show()
    plt.close('all')

