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

temp_stamp1 = []
temp_stamp2 = []
temp_stamp3 = []
temp_stamp4 = []
temp_stamp5 = []
temp_stamp6 = []
temp_stamp7 = []
# ################################# here data were flipped

new_tsmp = '20170124_145900' ## newer than
old_tsmp = '20170124_171700' ## older than

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

new_tsmp = '20170129_144400' ## newer than
old_tsmp = '20170129_156000' ## older than

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
#timestamp4.append(temp_stamp) ## newer than
new_tsmp = '20170124_171800' ## newer than
old_tsmp = '20170126_072700' ## older than

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



new_tsmp = '20170130_120400' ## newer than
old_tsmp = '20170130_160700' ## older than


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

    temp_stamp4.append(old_tsmp)
temp_stamp4 =temp_stamp4[::-1]

new_tsmp = '20170126_133000' ## older than
old_tsmp = '20170129_125900' ## older than


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

    temp_stamp5.append(old_tsmp)
temp_stamp5 =temp_stamp5[::-1]


new_tsmp = '20170505_185100' ## older than
old_tsmp = '20170510_094000' ## older than


search_string = '_DecouplingSequence_111_1_sil18_tau_'
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

    temp_stamp6.append(old_tsmp)
temp_stamp6 =temp_stamp6[::-1]


new_tsmp = '20170525_205100' ## older than
old_tsmp = '20170529_104000' ## older than


search_string = '_DecouplingSequence_111_1_sil18_tau_'
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

    temp_stamp7.append(old_tsmp)
temp_stamp7 =temp_stamp7[::-1]

timestamp = temp_stamp1+temp_stamp2+temp_stamp3+temp_stamp4+temp_stamp5+temp_stamp6+temp_stamp7
    
#timestamp4 = timestamp4[::-1]

# ################################# ############### ################

for kk in range(len(timestamp)):
    folder = toolbox.data_from_time(timestamp[kk])

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    cum_pts += a.pts

    if kk == 0:
        cum_sweep_pts = a.sweep_pts
        cum_p0 = a.p0
        cum_u_p0 = a.u_p0
        #cum_tau_list = a.tau_list
    #elif kk in [5,10,15,21,26,31]: 
    else:
        cum_sweep_pts = np.concatenate((cum_sweep_pts, a.sweep_pts))
        cum_p0 = np.concatenate((cum_p0, a.p0))
        cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))
        #cum_tau_list = np.concatenate((cum_tau_list, a.tau_list))

a.pts   = cum_pts 
a.sweep_pts = cum_sweep_pts 
a.p0    = cum_p0
a.u_p0  = cum_u_p0


ax = a.plot_results_vs_sweepparam(ret='ax',fmt='o-',figsize=(120,10))
np.save('sweep_tau_32_cum_P',cum_p0)
np.save('sweep_tau_32_cum_sweep_pts',cum_sweep_pts)

