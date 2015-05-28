'''script to analyze multiple ssro runs '''

import os, sys
if os.name == 'posix':
    sys.path.append("/Users/"+os.getlogin()+"/Documents/teamdiamond/")
else:
    sys.path.append("/measuring/")
from analysis.scripts.QEC import Two_Qubit_Tomography; reload(Two_Qubit_Tomography)

from analysis.lib.tools import toolbox


timestamp_Tomo = '20150504_183000'
while 1:
	timestamp_Tomo, folder = toolbox.latest_data('TomoCheck', older_than = timestamp_Tomo, return_timestamp = True)
	print timestamp_Tomo
	timestamp_SSRO, folder = toolbox.latest_data('AdwinSSRO', older_than = timestamp_Tomo, newer_than = '20150315_000000', return_timestamp = True)
	print timestamp_SSRO
	Two_Qubit_Tomography.BarPlotTomo(timestamp=timestamp_Tomo)