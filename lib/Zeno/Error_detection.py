"""
Analysis for error detection in the Zeno setting

NK 2015
"""
import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from analysis.lib.QEC import ConditionalParity as CP
from analysis.lib.fitting import fit, common, ramsey;reload(common); reload(fit)
import matplotlib.cm as cm
import matplotlib as mpl; reload(mpl)

reload (CP)
import h5py
import csv

def load_Err_det_data(folder, ssro_calib_folder,
						post_select 			= True,
						nr_of_msmts 			= 1):