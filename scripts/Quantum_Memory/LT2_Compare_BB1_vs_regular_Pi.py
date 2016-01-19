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
	if len(c_idents) == 1:
		older_than_bb1 = older_than_tst_1C_bb1
		older_than_pi = older_than_tst_1C_pi

	

