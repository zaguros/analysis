'''Script to analyze the dynamical decoupling data for eigenstate
by MA''' 

import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt

from scipy.stats.stats import pearsonr
reload(common)


### Inputs ###

## Data location ##
measurement_name = ['adwindata']
timestamp1 = []
timestamp2 = []
timestamp3 = []
timestamp4 = []
timestamp5 = []
timestamp6 = []
timestamp7 = []
timestamp8 = []
timestamp9 = []
timestamp10 = []
timestamp11 = []
timestamp12 = []
timestamp13 = []
timestamp14 = []
timestamp15 = []
timestamp16 = []
timestamp17 = []
timestamp18 = []
timestamp19 = []
timestamp20 = []
timestamp21 = []
timestamp22 = []
timestamp23 = []
timestamp24 = []
timestamp25 = []
timestamp26 = []
timestamp27 = []
timestamp28 = []
timestamp29 = []
timestamp29_2 = [] 
timestamp29_3 = []
timestamp29_4 = []

timestamp30 = []
timestamp31 = []
timestamp32 = []
timestamp33 = []
timestamp34 = []
timestamp35 = []
timestamp36 = []
timestamp37 = []
timestamp38 = []
timestamp39 = []
timestamp40 = []
timestamp41 = []
timestamp42 = []
timestamp43 = []
timestamp44 = []

timestamp64 = []
timestamp128=[]
timestamp256=[]
timestamp512=[]
timestamp1024=[]
timestamp16_m7=[]
timestamp16_p4=[]
timestamp16_m2=[]
temperature_list = []
absolute_time_list = []
temp_stamp1 = []
temp_stamp2 = []
temp_stamp3 = []
temp_stamp11 = []
temp_stamp12 = []
x_list = []
y_list = []
T_list = []
# ##############################################


##################################################### Data for N=2048 #########
new_tsmp = '20170216_004200' ## newer than
old_tsmp = '20170218_135900' ## older than

search_string = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY8sweep_tau_N_2048'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp1.append(old_tsmp)
    
timestamp1 = timestamp1[::-1]



#########################################################
new_tsmp = '20170218_130000' ## newer than
old_tsmp = '20170219_043000' ## older than

search_string = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY8sweep_tau_N_2048'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp2.append(old_tsmp)
    
timestamp2 = timestamp2[::-1]

# ###########################################
new_tsmp = '20170219_120000' ## newer than
old_tsmp = '20170220_043000' ## older than

search_string = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY8sweep_tau_N_2048'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp3.append(old_tsmp)
    
timestamp3 = timestamp3[::-1]

# Data for N= 3072 #######################################################################

# ################################# here data were flipped

new_tsmp = '20170220_131000' ## newer than
old_tsmp = '20170220_143100' ## older than

search_string = '_DecouplingSequence_111_1'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)

    temp_stamp1.append(old_tsmp)
temp_stamp1 =temp_stamp1[::-1]
#timestamp4.append(temp_stamp)

new_tsmp = '20170221_020200' ## newer than
old_tsmp = '20170221_024500' ## older than

search_string = '_DecouplingSequence_111_1'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)

    temp_stamp2.append(old_tsmp)
temp_stamp2 =temp_stamp2[::-1]
#timestamp4.append(temp_stamp)

new_tsmp = '20170220_143100' ## newer than
old_tsmp = '20170220_180000' ## older than


search_string = '_DecouplingSequence_111_1'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)

    temp_stamp3.append(old_tsmp)
temp_stamp3 =temp_stamp3[::-1]
timestamp4=temp_stamp1+temp_stamp2+temp_stamp3
    
#timestamp4 = timestamp4[::-1]

# ################################# ############### ################

new_tsmp = '20170221_044000' ## newer than
old_tsmp = '20170221_200800' ## older than

search_string = '_DecouplingSequence_111_1_'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp5.append(old_tsmp)
    
timestamp5 = timestamp5[::-1]

# #################################

new_tsmp = '20170221_202500' ## newer than
old_tsmp = '20170222_150900' ## older than

search_string = '_DecouplingSequence_111_1'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp6.append(old_tsmp)
    
timestamp6 = timestamp6[::-1]

# ################################# N= 6144 ###############################

new_tsmp = '20170222_154000' ## newer than
old_tsmp = '20170222_203500' ## older than

search_string = '_DecouplingSequence_111_'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp7.append(old_tsmp)
    
timestamp7 = timestamp7[::-1]
#################################

# #################################

new_tsmp = '20170222_203800' ## newer than
old_tsmp = '20170223_012700' ## older than

search_string = '_DecouplingSequence_111'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp8.append(old_tsmp)
    
timestamp8 = timestamp8[::-1]

# #################################

new_tsmp = '20170223_012700' ## newer than
old_tsmp = '20170223_130000' ## older than

search_string = '_DecouplingSequence_111_'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp9.append(old_tsmp)
    
timestamp9 = timestamp9[::-1]


new_tsmp = '20170223_133000' ## newer than
old_tsmp = '20170223_190000' ## older than

search_string = '_DecouplingSequence_111_1_'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp13.append(old_tsmp)
    
timestamp13 = timestamp13[::-1]

# #################################

new_tsmp = '20170223_220000' ## newer than
old_tsmp = '20170224_090000' ## older than

search_string = 'DecouplingSequence_111_1_'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp14.append(old_tsmp)
    
timestamp14 = timestamp14[::-1]

# #################################

new_tsmp = '20170224_100000' ## newer than
old_tsmp = '20170224_180000' ## older than

search_string = 'DecouplingSequence_111'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp15.append(old_tsmp)
    
timestamp15 = timestamp15[::-1]

# #################################

new_tsmp = '20170228_150000' ## newer than
old_tsmp = '20170228_180000' ## older than

search_string = 'DecouplingSequence_111'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)

    temp_stamp11.append(old_tsmp)
temp_stamp11 =temp_stamp11[::-1]

new_tsmp = '20170228_095300' ## newer than
old_tsmp = '20170228_150000' ## older than

search_string = 'DecouplingSequence_111'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)

    temp_stamp12.append(old_tsmp)
temp_stamp12 =temp_stamp12[::-1]

     
timestamp16 = temp_stamp11+temp_stamp12


# ################################# N= 10240 ###############################

new_tsmp = '20170224_183000' ## newer than
old_tsmp = '20170225_013000' ## older than

search_string = '_DecouplingSequence_111_'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp20.append(old_tsmp)
    
timestamp20 = timestamp20[::-1]
#################################

new_tsmp = '20170225_013000' ## newer than
old_tsmp = '20170225_113000' ## older than

search_string = 'DecouplingSequence_111_1_'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp21.append(old_tsmp)
    
timestamp21 = timestamp21[::-1]

# #################################

new_tsmp = '20170225_113000' ## newer than
old_tsmp = '20170226_033000' ## older than

search_string = 'DecouplingSequence_111_1_'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp22.append(old_tsmp)
    
timestamp22 = timestamp22[::-1]

# #################################

new_tsmp = '20170226_033000' ## newer than
old_tsmp = '20170226_110000' ## older than

search_string = '_DecouplingSequence_111'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp23.append(old_tsmp)
    
timestamp23 = timestamp23[::-1]



# #################################

new_tsmp = '20170226_110000' ## newer than
old_tsmp = '20170226_213000' ## older than

search_string = '_DecouplingSequence_111_1_sil18_S'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp24.append(old_tsmp)
    
timestamp24 = timestamp24[::-1]

# #################################

# #################################

new_tsmp = '20170226_213000' ## newer than
old_tsmp = '20170227_083000' ## older than

search_string = '_DecouplingSequence_111_1_sil18_S'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp25.append(old_tsmp)
    
timestamp25 = timestamp25[::-1]

# #################################

new_tsmp = '20170227_182300' ## newer than
old_tsmp = '20170228_023100' ## older than

search_string = '_DecouplingSequence_111_1_sil18_S'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp26.append(old_tsmp)
    
timestamp26 = timestamp26[::-1]


# ################################# tau + 16 ns

new_tsmp = '20170320_225000' ## newer than
old_tsmp = '20170321_231500' ## older than

search_string = '_DecouplingSequence_111_1_sil18_S'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp27.append(old_tsmp)
    
timestamp27 = timestamp27[::-1]


# ################################# tau + 24 ns

new_tsmp = '20170322_072000' ## newer than
old_tsmp = '20170322_151500' ## older than

search_string = '_DecouplingSequence_111_1_sil18_S'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp28.append(old_tsmp)
    
timestamp28 = timestamp28[::-1]

# ################################# tau + 32 ns

new_tsmp = '20170322_213800' ## newer than
old_tsmp = '20170323_051500' ## older than

search_string = '_DecouplingSequence_111_1_sil18_S'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp29.append(old_tsmp)
    
timestamp29 = timestamp29[::-1]

# ################################# tau -16 ns

new_tsmp = '20170323_080000' ## newer than
old_tsmp = '20170323_161500' ## older than

search_string = '_DecouplingSequence_111_1_sil18_S'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp29_2.append(old_tsmp)
    
timestamp29_2 = timestamp29_2[::-1]

# ################################# tau -24 ns

new_tsmp = '20170326_192800' ## newer than
old_tsmp = '20170327_031500' ## older than

search_string = '_DecouplingSequence_111_1_sil18_S'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp29_3.append(old_tsmp)
    
timestamp29_3 = timestamp29_3[::-1]

# ################################# tau -2 ns

new_tsmp = '20170327_060000' ## newer than
old_tsmp = '20170327_141500' ## older than

search_string = '_DecouplingSequence_111_1_sil18_S'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp29_4.append(old_tsmp)
    
timestamp29_4 = timestamp29_4[::-1]
# ################################# #############  ######################

#### #####                        N=1024        #########
 
new_tsmp = '20170217_201500' ## newer than
old_tsmp = '20170218_000200' ## older than

search_string = '_DecouplingSequence_111_1'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp30.append(old_tsmp)
    
timestamp30 = timestamp30[::-1]

new_tsmp = '20170301_150000' ## newer than
old_tsmp = '20170301_200000' ## older than

search_string = '_DecouplingSequence_111_1'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp31.append(old_tsmp)
    
timestamp31 = timestamp31[::-1]


new_tsmp = '20170219_152000' ## newer than
old_tsmp = '20170219_191000' ## older than

search_string = '_DecouplingSequence_111_1'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp32.append(old_tsmp)
    
timestamp32 = timestamp32[::-1]

# ################################# #############  ######################

#### #####                        N=512        #########


new_tsmp = '20170217_174300' ## newer than
old_tsmp = '20170217_201500' ## older than

search_string = '_DecouplingSequence_111_1'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp33.append(old_tsmp)
    
timestamp33 = timestamp33[::-1]

new_tsmp = '20170218_120000' ## newer than
old_tsmp = '20170218_143000' ## older than

search_string = '_DecouplingSequence_111_1'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp34.append(old_tsmp)
    
timestamp34 = timestamp34[::-1]


new_tsmp = '20170219_122000' ## newer than
old_tsmp = '20170219_152000' ## older than

search_string = '_DecouplingSequence_111_1'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp35.append(old_tsmp)
    
timestamp35 = timestamp35[::-1]


################################ N=256 ########

new_tsmp = '20170216_160000' ## newer than
old_tsmp = '20170216_232000' ## older than

search_string = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY8sweep_tau_N_256'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp36.append(old_tsmp)
    
timestamp36 = timestamp36[::-1]

################################ N=128 ########

new_tsmp = '20170216_160000' ## newer than
old_tsmp = '20170216_232000' ## older than

search_string = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY8sweep_tau_N_128'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp37.append(old_tsmp)
    
timestamp37 = timestamp37[::-1]

################################ N=64 ########

new_tsmp = '20170216_160000' ## newer than
old_tsmp = '20170216_232000' ## older than

search_string = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY8sweep_tau_N_64'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp38.append(old_tsmp)
    
timestamp38 = timestamp38[::-1]

################################ N=32 ########

new_tsmp = '20170216_160000' ## newer than
old_tsmp = '20170216_232000' ## older than

search_string = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY8sweep_tau_N_32'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp39.append(old_tsmp)
    
timestamp39 = timestamp39[::-1]

################################ N=16 ########

new_tsmp = '20170309_111400' ## newer than
old_tsmp = '20170309_120000' ## older than

search_string = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY8sweep_tau_N_16'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp40.append(old_tsmp)
    
timestamp40 = timestamp40[::-1]

################################ N=8 ########

new_tsmp = '20170216_160900' ## newer than
old_tsmp = '20170216_232000' ## older than

search_string = 'DecouplingSequence_111_1_sil18_RepT_ShutterNO_XY8sweep_tau_N_8'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp41.append(old_tsmp)
    
timestamp41 = timestamp41[::-1]

################################ N=4 ########

new_tsmp = '20170320_144200' ## newer than
old_tsmp = '20170320_144900' ## older than

search_string = 'DecouplingSequence_111_1_'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp43.append(old_tsmp)
    
timestamp43 = timestamp43[::-1]

################################ N=2 ########

new_tsmp = '20170320_153200' ## newer than
old_tsmp = '20170320_154200' ## older than

search_string = 'DecouplingSequence_111_1_'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp44.append(old_tsmp)
    
timestamp44 = timestamp44[::-1]

################################ N=1 spin echo ########


new_tsmp = '20170308_005000' ## newer than
old_tsmp = '20170308_012000' ## older than

search_string = 'DecouplingSequence_111_1_sil18_'
while toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False) != False:
    old_tsmp, folder = toolbox.latest_data(contains=search_string,
                                        return_timestamp =True,
                                        older_than=old_tsmp,
                                        newer_than=new_tsmp,
                                        raise_exc=False)


    timestamp42.append(old_tsmp)
    
timestamp42 = timestamp42[::-1]

########################################################

# timestamp_list=[timestamp1,timestamp2,timestamp3,timestamp4,timestamp5,timestamp6,timestamp7,timestamp8,timestamp9,timestamp10,timestamp11,timestamp12,timestamp13,timestamp14]#,timestamp14,timestamp15,timestamp16,timestamp17]
timestamp_list_2048=[timestamp1, timestamp2, timestamp3]
timestamp_list_3072=[timestamp4, timestamp5, timestamp6]
timestamp_list_6144=[timestamp7, timestamp8, timestamp9, timestamp13, timestamp14, timestamp15, timestamp16]
timestamp_list_1024=[timestamp30, timestamp31, timestamp32]
timestamp_list_512=[timestamp33, timestamp34, timestamp35]
timestamp_list_10240=[timestamp20,timestamp21,timestamp22,timestamp23,timestamp24,timestamp25,timestamp26,timestamp27,timestamp28,timestamp29,timestamp29_2,timestamp29_3,timestamp29_4]
timestamp_list_8_256=[timestamp42,timestamp44,timestamp43,timestamp41, timestamp40, timestamp39,timestamp38,timestamp37,timestamp36]
timestamp_list= timestamp_list_8_256
combined_timestamp_list=[timestamp_list_512,timestamp_list_1024,timestamp_list_2048,timestamp_list_3072,timestamp_list_6144,timestamp_list_10240]


#labels = ['M*tau_larmor','M*tau_larmor+4 ns','M*tau_larmor-4ns','M*tau_larmor+2ns', 'M*tau_larmor+8ns','M*tau_larmor-8ns','M*tau_larmor+19ns']
#labels = ['M*tau_larmor','M*tau_larmor+4 ns','M*tau_larmor+8ns','M*tau_larmor-4ns', 'M*tau_larmor-8ns','M*tau_larmor+2ns','M*tau_larmor+19ns']
labels = ['1','2','4','8','16','32','64','128','256','512','1024','2048','3072','6144','10240']
N = [1,2,4,8,16,32,64,128,256,512,1024,2048,3072,6144,10240]
color_list =  ['black','darkred','r','orangered','gold','darkgreen','g','b','darkblue','indigo','violet','c','black','darkred','r','orangered','gold','darkgreen','g','b','darkblue','indigo','violet','c']
fmt_label = ['o','s','^','o','s','^','o','s','^','o','s','^','o','s','^','o','s','^','o','s','^','o','s','^']


#print 'time stamp is    ' +str(timestamp4)

offset      = 0.5
amplitude   = 0.45
position    = 0
T2          = 1
power       = 1.5

## other settings ##
plot_fit    = False
show_guess  = False

fit_results = []

cum_pts = 0
cum_sweep_pts       = np.empty(0)
cum_p0              = np.empty(0)
cu                  =np.empty(0) 
cu_u                =np.empty(0) 
cum_u_p0            = np.empty(0)

cum_normalized_ssro = np.empty(0)
e_list=[]



for ii, timestamp in enumerate(timestamp_list):
    #print timestamp

    for kk in range(len(timestamp)):

        folder = toolbox.data_from_time(timestamp[kk])

        a = mbi.MBIAnalysis(folder)
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        a.get_electron_ROC()
        cum_pts += a.pts
        temperature = (a.g.attrs['temp']-100)/0.385
        temperature_list.append(temperature)


        if kk == 0:
            cum_sweep_pts = a.sweep_pts
            cum_p0 = a.p0
            cum_u_p0 = a.u_p0
        else:
            cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
            cum_p0 = np.concatenate((cum_p0, a.p0))
            cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

    
    a.pts   = cum_pts
    a.sweep_pts = cum_sweep_pts
    a.p0    = cum_p0
    a.u_p0  = cum_u_p0

    x_list.append(list(cum_sweep_pts))
    y_list.append(list(cum_p0))
    e_list.append(list(cum_u_p0))

   
 
for jj in range(0, len(combined_timestamp_list)):
    timestamp_list= combined_timestamp_list[jj]

    for ii, timestamp in enumerate(timestamp_list):
        #print timestamp

        for kk in range(len(timestamp)):

            folder = toolbox.data_from_time(timestamp[kk])

            a = mbi.MBIAnalysis(folder)
            a.get_sweep_pts()
            a.get_readout_results(name='adwindata')
            a.get_electron_ROC()
            cum_pts += a.pts
            temperature = (a.g.attrs['temp']-100)/0.385
            temperature_list.append(temperature)


            if kk == 0:
                cum_sweep_pts = a.sweep_pts
                cum_p0 = a.p0
                cum_u_p0 = a.u_p0
            else:
                cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
                cum_p0 = np.concatenate((cum_p0, a.p0))
                cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))


        a.pts   = cum_pts
        a.sweep_pts = cum_sweep_pts
        a.p0    = cum_p0
        a.u_p0  = cum_u_p0

        if ii == 0:
            cum_s = cum_sweep_pts
            cum_e = cum_u_p0
            cum1 =  cum_p0
        elif ii == 1:
            cum2 = cum_p0
        elif ii== 2:
            cum3 = cum_p0
        elif ii == 3:
            cum4 = cum_p0
        elif ii == 4:
            cum5 = cum_p0
        elif ii == 5:
            cum6 = cum_p0
        elif ii == 6:
            cum7 = cum_p0
        elif ii == 7:
            cum8 = cum_p0
        elif ii == 8:
            cum9 = cum_p0
        elif ii == 9:
            cum10 = cum_p0
        elif ii == 10:
            cum11 = cum_p0
        elif ii ==11:
            cum12 = cum_p0
        elif ii ==12:
            cum13 = cum_p0

    max_list = len(cum1)*[0]

    for kk in range(0,len(cum1)):
        if jj in range(0,4):
            max_list[kk]=max(cum1[kk],cum2[kk],cum3[kk])
        elif jj in range(4,5):
            max_list[kk]=max(cum1[kk],cum2[kk],cum3[kk],cum4[kk],cum5[kk],cum6[kk],cum7[kk])
        elif jj in range(5,6):
            max_list[kk]=max(cum1[kk],cum2[kk],cum3[kk],cum4[kk],cum5[kk],cum6[kk],cum7[kk],cum8[kk],cum9[kk],cum10[kk],cum11[kk],cum12[kk],cum13[kk])
            
    x_list.append(list(cum_s))  
    y_list.append(max_list) 
    e_list.append(list(cum_e)) 


    cum_pts = 0
    cum_sweep_pts       = np.empty(0)
    cum_p0              = np.empty(0)
    cum_u_p0            = np.empty(0)
    cum_s=np.empty(0)
    cum1=np.empty(0)
    cum2=np.empty(0)
    cum3=np.empty(0)
    cum4=np.empty(0)
    cum5=np.empty(0)
    cum6=np.empty(0)
    cum7=np.empty(0)


np.save('final_data_x',np.array(x_list))
np.save('final_data_y',np.array(y_list))
       
np.savez('final_data',x_list=x_list, y_list=y_list,e_list=e_list )


######################################### plotting data ################

#ax1.set_color_cycle(color_list)



def plot_final_data(filename='final_data'):
    d = np.load('final_data.npz')
    x_list=d['x_list']
    y_list=d['y_list']
    e_list=d['e_list']

    fig = plt.figure(2,figsize=(30,10))
    ax1 = fig.add_subplot(111)
    N = [1,2,4,8,16,32,64,128,256,512,1024,2048,3072,6144,10240]
    
    for jj in range(0,len(x_list)):

        
        ##############  Fitting the data #########
        x= 1000*np.array(x_list[jj]) /(2*N[jj])
        y= np.array(y_list[jj])
        e= np.array(e_list[jj])

        x = x.reshape(-1)[:]
        y = y.reshape(-1)[:]
        e = e.reshape(-1)[:]


        #########################################################
        T2= 1
        
        if jj ==3 :
            x=np.delete(x,[28,29,30,31,32,33,34,43,44,45,46,47,48,49,50,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130])
            y=np.delete(y,[28,29,30,31,32,33,34,43,44,45,46,47,48,49,50,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130])
        if jj ==4 :
            x=np.delete(x,[28,29,30,31,32,33,34,43,44,45,46,47,48,49,50,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127])
            y=np.delete(y,[28,29,30,31,32,33,34,43,44,45,46,47,48,49,50,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127])

        elif jj in range(5,len(x_list)-2):
            x=np.delete(x,[25,28,29,30,31,32,33,34,43,44,45,46,47,48,49,50,69,70,71,72,73,74,75,76])
            y=np.delete(y,[25,28,29,30,31,32,33,34,43,44,45,46,47,48,49,50,69,70,71,72,73,74,75,76])
            T2=10

        elif jj in range(len(x_list)-3,len(x_list)-2):
            x=np.delete(x,[26,27,28,29,30,31,32,41,42,43,44,45,46,47,48,67,68,69,70,71,72,73,74])
            y=np.delete(y,[26,27,28,29,30,31,32,41,42,43,44,45,46,47,48,67,68,69,70,71,72,73,74])
            T2=10

        elif jj in range(len(x_list)-2,len(x_list)-1):
            x=np.delete(x,[3,7,8,9,11,12,13,15,26,27,28,30,31])
            y=np.delete(y,[3,7,8,9,11,12,13,15,26,27,28,30,31])
            T2=10    

        elif jj in range(len(x_list)-1,len(x_list)):
            x=np.delete(x,[3,7,8,9])#,10,12,13,14,26,27,28,30,31])
            y=np.delete(y,[3,7,8,9])#,10,12,13,14,26,27,28,30,31])
            T2=10
        
        
        
        ax1.errorbar(x.flatten(),y.flatten(),yerr=e[jj],fmt=fmt_label[jj],label=labels[jj],color=color_list[jj])

        
        p0, fitfunc, fitfunc_str = common.fit_general_exponential(offset, amplitude, position, T2, power)

        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=[0,2])

        # plot data and fit as function of total time
        plot_fit= True

        if plot_fit == True:
            plot.plot_fit1d(fit_result, np.linspace(0,np.amax(x)+400.,1001), ax=ax1,add_txt = False, plot_data=False,color=color_list[jj])

        fit_results.append(fit_result)
        T_list.append(fit_result['params_dict']['T'])
        
        
    ##################################
    ax1.set_xlabel('Total evolution time (ms)')
    ax1.set_xscale('log')
    ax1.set_xlim(1.e-2,1.e5)
    ax1.grid()
    ax1.legend()

    plt.savefig(os.path.join(folder, 'selected_max_vs_tot_time_log.pdf'),
    format='pdf')
    plt.savefig(os.path.join(folder, 'selected_max_vs_tot_time_log.png'),
    format='png')

    

    del N[0:2]
    del T_list[0:2]
    del N[-1]
    del T_list[-1]

    N = list(np.array(N)/1)
    #################### plot and fit scaling with N
    fig = plt.figure(3,figsize=(8,6))
    ax2 = fig.add_subplot(111)
    ax2.plot(N,T_list,'bo')




        
    B_guess = T_list[0]/4 
    n_guess = 0.78         

    #A = fit.Parameter(A_guess, 'A')
    B = fit.Parameter(B_guess, 'B')
    n=  fit.Parameter(n_guess, 'n')

    def fitfunc(x):
                return  B()*x**n()

    print 'running fit'
    fit_result = fit.fit1d(N, T_list, None, p0 = [ B, n],
                    fitfunc = fitfunc, do_print=False, ret=True, fixed=[])

    #A0 = fit_result['params_dict']['A']
    B0 = fit_result['params_dict']['B']
    n0 = fit_result['params_dict']['n']

        # plotting the fitted function
    plot.plot_fit1d(fit_result, np.linspace(min(N), max(N), 1000), ax=ax2, plot_data=True)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('N')
    ax2.set_ylabel('coherence time (ms)')





