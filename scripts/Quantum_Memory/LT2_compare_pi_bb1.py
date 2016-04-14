"""
Script that compares carbon dephasing for BB1 pulses and regular pi pulses (both hermite)
"""
 
import numpy as np
import os,h5py
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
import analysis.lib.Qmemory.CarbonDephasing as CD
 
older_than_tst_2C_bb1 = '20150826_062655' ### unfortunately it cannot be the same as ...1C_BB1
older_than_tst_2C_pi = '20150827_120000'
older_than_tst_1C_bb1 = '20150826_090000'
older_than_tst_1C_pi = '20150827_120000'
 
c_idents = ['12','13','15','16','23','25','26','35','36','1','2','3','5','6']


for c in c_idents:

### choose correct timestamps

    if len(c) == 1:
        older_than_bb1 = older_than_tst_1C_bb1
        older_than_pi = older_than_tst_1C_pi
        fig = plt.figure()
        ax =  plt.subplot()
        xBB,yBB,yBB_u,folderBB = CD.Sweep_repetitions(carbon = c,older_than = older_than_bb1,plot_result = False)
        xpi,ypi,ypi_u,folderpi = CD.Sweep_repetitions(carbon = c,older_than = older_than_pi,plot_result = False,
                                                            folder_name = 'Memory_NoOf_Repetitions_regularPi_')

        plt.errorbar(xBB,yBB,yBB_u,fmt = 'o',label = 'hermite BB1')
        plt.errorbar(xpi,ypi,ypi_u,fmt = 'o',label = 'hermite pi')
        plt.title('C'+str(c)+' ' + CD.get_tstamp_from_folder(folder))
        plt.xlabel('Number of repetitions')
        plt.ylabel('Bloch vector length')
        plt.legend()
        plt.savefig(os.path.join(folder,'BB1_vs_pi.png'),format='png')
        plt.show()


    if len(c) == 2:

        # select timestamp limit

        older_than_bb1 = older_than_tst_2C_bb1
        older_than_pi = older_than_tst_2C_pi


        # get data
        xXBB,yXBB,yXBB_u,folderBB = CD.Sweep_repetitions(carbon = c,logicstate='X',
                                                            older_than = older_than_bb1,plot_result = False)
        xXpi,yXpi,yXpi_u,folderpi = CD.Sweep_repetitions(carbon = c,logicstate='X',
                                                            older_than = older_than_pi,plot_result = False,
                                                            folder_name = 'Memory_NoOf_Repetitions_regularPi_')

        xmXBB,ymXBB,ymXBB_u,folderBB = CD.Sweep_repetitions(carbon = c,logicstate='mX',
                                                            older_than = older_than_bb1,plot_result = False)
        xmXpi,ymXpi,ymXpi_u,folderpi = CD.Sweep_repetitions(carbon = c,logicstate='mX',
                                                            older_than = older_than_pi,plot_result = False,
                                                            folder_name = 'Memory_NoOf_Repetitions_regularPi_')
        

        #start plotting
        fig = plt.figure()
        ax =  plt.subplot()
        plt.errorbar(xXBB,yXBB,yXBB_u,fmt = 'o',label = '+X: hermite BB1')
        plt.errorbar(xXpi,yXpi,yXpi_u,fmt = 'o',label = '+X: hermite pi')
        plt.errorbar(xmXBB,ymXBB,ymXBB_u,fmt = 'o',label = '-X: hermite BB1')
        plt.errorbar(xmXpi,ymXpi,ymXpi_u,fmt = 'o',label = '-X: hermite pi')

        plt.title('C'+str(c)+' ' + CD.get_tstamp_from_folder(folder))
        plt.xlabel('Number of repetitions')
        plt.ylabel('Bloch vector length')
        plt.legend()
        plt.savefig(os.path.join(folder,'BB1_vs_pi.png'),format='png')
        plt.show()

plt.close('all')

