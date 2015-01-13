
from analysis.lib.fitting import dark_esr_auto_analysis; reload(dark_esr_auto_analysis)
import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib.tools import plot
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit,esr
ZFS = 2.877623

f0m = []; u_f0m = [];
f0p = []; u_f0p = [];
f_centre = []
f_diff = []

timestamp = '20141220_140000'
k_list = np.linspace(0,20,21)

for k in range(len(k_list)):

	timestamp, folder = toolbox.latest_data('PulsarDarkESR_magnet_msm1', older_than = timestamp, return_timestamp = True)
	print folder
	f0m_temp, u_f0m_temp = dark_esr_auto_analysis.analyze_dark_esr_double(timestamp = timestamp)

	timestamp, folder = toolbox.latest_data('PulsarDarkESR_magnet_msp1', older_than = timestamp, return_timestamp = True)
	f0p_temp, u_f0p_temp = dark_esr_auto_analysis.analyze_dark_esr_double(timestamp = timestamp)

	f_centre_temp    = (f0m_temp+f0p_temp)/2
	f_diff_temp = (f_centre_temp-ZFS)*1e6

	print timestamp

	f0m.append(f0m_temp)
	u_f0m.append(u_f0m_temp)
	f0p.append(f0p_temp)
	u_f0p.append(u_f0p_temp)
	f_centre.append(f_centre_temp)
	f_diff.append(f_diff_temp)

fig, ax = plt.subplots(1,1)
ax.errorbar(k_list, f0m, u_f0m,'k',label = 'f0m')
ax.errorbar(k_list, f0p, u_f0p,'b',label = 'f0p')
ax.legend()