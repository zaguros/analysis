'''script to analyze multiple ssro runs '''

import os, sys
if os.name == 'posix':
    sys.path.append("/Users/"+os.getlogin()+"/Documents/teamdiamond/")
else:
    sys.path.append("/measuring/")
from analysis.lib.m2.ssro import ssro
reload(ssro)
from analysis.lib.tools import toolbox



timestamp = '20150505_081846'

while 1:
	timestamp, folder = toolbox.latest_data('AdwinSSRO', older_than = timestamp, return_timestamp = True)

	print folder
	ssro.ssrocalib(folder = folder) 

# folder_list = [ r'D:\measuring\data\20141219\184417_AdwinSSRO_SSROCalibration_111_1_sil18']
# readout_time = 114

# # folder_list = [
# # r'D:\measuring\data\20150124\102326_AdwinSSRO_SSROCalibration_111_1_sil18',
# # r'D:\measuring\data\20150126\233614_AdwinSSRO_SSROCalibration_111_1_sil18',
# # r'D:\measuring\data\20150127\145611_AdwinSSRO_SSROCalibration_111_1_sil18',
# # r'D:\measuring\data\20150201\183742_AdwinSSRO_SSROCalibration_111_1_sil18',
# # r'D:\measuring\data\20150202\012642_AdwinSSRO_SSROCalibration_111_1_sil18',
# # r'D:\measuring\data\20150204\225213_AdwinSSRO_SSROCalibration_111_1_sil18',
# # r'D:\measuring\data\20150205\110722_AdwinSSRO_SSROCalibration_111_1_sil18',
# # r'D:\measuring\data\20150206\041814_AdwinSSRO_SSROCalibration_111_1_sil18'
# # ]

# f0_tot = 0
# u_f0_tot = 0
# f1_tot = 0
# u_f1_tot = 0

# if 1:
# 	for folder in folder_list:
# 		print 'ok'
# 		f0_temp, u_f0_temp,f1_temp, u_f1_temp = ssro.get_SSRO_calibration(folder, readout_time)
		
# 		f0_tot += f0_temp
# 		u_f0_tot += u_f0_temp**2
# 		f1_tot += f1_temp
# 		u_f1_tot += u_f1_temp**2

# 		print f0_temp
# 		print f1_temp



# 	f0_tot = f0_tot/len(folder_list)
# 	u_f0_tot = u_f0_tot**0.5/len(folder_list)
# 	f1_tot = f1_tot/len(folder_list)
# 	u_f1_tot = (u_f1_tot)**0.5/len(folder_list)

# 	print f0_tot , u_f0_tot,f1_tot , u_f1_tot, 1/2.*(f0_tot+f1_tot),1/2.*(u_f0_tot**2+u_f1_tot**2)**0.5
# <<<<<<< Updated upstream
# =======
# >>>>>>> origin/master
# >>>>>>> Stashed changes
