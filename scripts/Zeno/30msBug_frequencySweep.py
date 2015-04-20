import numpy as np
import os,re,sys
import h5py
sys.path.append(r'D:/measuring')
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
reload(common)
reload(toolbox)

"""
Analyze the zeno data
NK 2014
"""

def analyze_dip(older_than='20150318_025003',newer_than='20150317_152433',ssro_timestamp=None,tag='_ZZC1f0'):

	"""
	Analyzes the test states. In a time window which is specified by older_than and newer_than in timestamp format.
	Analyzes a maximum of N states.
	"""

	search_string='__ROBasis'+tag

	x_arr=[]
	y_arr=[]
	y_u_arr=[]
	fig=plt.figure()
	ax=plt.subplot

	GoOn = True
	jj = 0
	while GoOn:
		#if the desired data set exists, then read the measured values.
		if toolbox.latest_data(contains=search_string,
										return_timestamp =True,
										older_than=older_than,
										newer_than=newer_than,
										raise_exc=False) != False:
			older_than,folder=toolbox.latest_data(contains=search_string,
										return_timestamp =True,
										older_than=older_than,
										newer_than=newer_than,
										raise_exc=False)

			evotime,c1ms0,y,y_err= Zeno_get_2Q_values(older_than,ssro_calib_timestamp=ssro_timestamp)
			x_arr.append(c1ms0)
			if jj == 0:
				### initialize the results and uncertainty arrays
				y_arr = [[y[i]] for i in range(len(evotime))]
				y_u_arr = [[y_err[i]] for i in range(len(evotime))]
				### append results
			else:
				for kk in range(len(evotime)):
					y_arr[kk].append(y[kk])
					y_u_arr[kk].append(y_err[kk])
			jj=1
		else: GoOn = False

	#condense the results for positive and negative read_out
	for i in range(len(x_arr)/2+1):
		if i==0:
			pass
		else:
			for kk in range(len(evotime)):
				y_arr[kk][i-1]=y_arr[kk][2*i-1]/2.-y_arr[kk][2*i-2]/2.
				y_u_arr[kk][i-1]=np.sqrt(y_u_arr[kk][2*i-1]**2+y_u_arr[kk][2*i-2]**2)/2.
	

	for kk in range(len(evotime)):
		y_arr[kk]=y_arr[kk][:len(x_arr)/2]
		y_u_arr[kk] = y_u_arr[kk][:len(x_arr)/2]
	for i in range(len(x_arr)):
		if i%2==0:
			x_arr[i/2]=x_arr[i]
	x_arr=x_arr[:len(x_arr)/2]
	x_arr=[x-431991. for x in x_arr]


	### sort by fictual detuning
	### numpy has a very nice function argsort for this.

	# print np.argsort(x_arr)
	x_arr=np.array(x_arr)
	for kk in range(len(evotime)):
		y_arr[kk] = np.array(y_arr[kk])[np.argsort(x_arr)]
		y_u_arr[kk] = np.array(y_u_arr[kk])[np.argsort(x_arr)]

	x_arr=x_arr[np.argsort(x_arr)]

	for kk in range(len(evotime)):
		plt.errorbar(x_arr,y_arr[kk],y_u_arr[kk],marker='o',label=str(evotime[kk]*1e3)+' us')
	plt.xlabel('Delta f_C1_ms0')
	plt.ylabel('ZZ contrast')
	plt.title('To dip or not to dip: '+folder[17:33])
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	print 'Plot is saved in'
	print folder[17:33]
	plt.savefig(os.path.join(folder,'Zeno_Dip.pdf'),format='pdf')
	plt.savefig(os.path.join(folder,'Zeno_Dip.png'),format='png')
	plt.show()
	plt.close('all')



def Zeno_get_2Q_values(timestamp=None, folder=None,folder_name='Zeno',
						measurement_name = ['adwindata'], 
						ssro_calib_timestamp =None):
	"""
	Returns the relevant RO values for a given timestamp.
	"""

	if timestamp == None and folder==None:
	    timestamp, folder   = toolbox.latest_data(folder_name,return_timestamp =True)
	elif timestamp ==None and folder!=None: 
		pass
	else:
	    folder = toolbox.data_from_time(timestamp)

	if folder != None and timestamp == None:
		d,t = toolbox.get_date_time_string_from_folder(folder)
		timestamp = toolbox.timestamp_from_datetime(t)
	### SSRO calibration
	if ssro_calib_timestamp == None: 
	    ssro_calib_folder = toolbox.latest_data('SSRO',older_than=timestamp)



	else:
	    ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
	    ssro_calib_folder = toolbox.datadir + '/'+ssro_dstmp+'/'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'


	if 'Evotime' in folder:
		c1ms0 = float(folder[114:]) ##assign the frequency from the folder name
	else:
		c1ms0 = float(folder[-8:]) 
	a = mbi.MBIAnalysis(folder)
	a.get_sweep_pts()
	a.get_readout_results(name='adwindata')
	a.get_electron_ROC(ssro_calib_folder)

	x_labels = a.sweep_pts.reshape(-1)
	y= ((a.p0.reshape(-1))-0.5)*2
	y_err = 2*a.u_p0.reshape(-1)

	return x_labels,c1ms0,y,y_err


### original data for the original timing where the dip occured.
# analyze_dip()

### dip for other timings.
# analyze_dip(older_than='20150318_130003',newer_than='20150318_090433',tag='_ZZnewEvotimesC1f0')

### zoom in on the dip.
# analyze_dip(older_than='20150318_113918',newer_than='20150318_100852',tag='_ZZnewEvotimesC1f0')

### reproduce dip for other timings.
#analyze_dip(older_than='20150319_142924',newer_than='20150319_115523')


### for a single zeno measurement and varying times
# analyze_dip(older_than='20150325_235959',newer_than='20150325_191334')

### for a freshly calibrated system. 
## partial data set for 0.0 us
# analyze_dip(older_than='20150326_153739',newer_than='20150326_151247')

## data set for 3 points.
# analyze_dip(older_than='20150326_190007',newer_than='20150326_153903')

### reproduction of the same data set with current calibrations (both sets are taken within the same filling.)
# analyze_dip(older_than='20150327_144603',newer_than='20150327_115035')


## 3 measurements, freshly calibrated.
# analyze_dip(older_than='20150327_162706',newer_than='20150327_144603')

### 3 measurements, freshly calibrated.

analyze_dip(older_than='20150329_173845',newer_than='20150329_160058')

print 'classical zz correlations'
### classically initialized states --> dip is found!!
# analyze_dip(older_than='20150318_152530',newer_than='20150318_141720')

### reproduce the results from above with a higher resolution.
# analyze_dip(older_than='20150318_172439',newer_than='20150318_162218')

### reproduce the dip at 0 detuning with phase gates around the larmor revival.

# analyze_dip(older_than='20150323_153439',newer_than='20150323_144313')

### classically initialized states, different timings --> we don't find a dip
# analyze_dip(older_than='20150318_161722',newer_than='20150318_154202')

### try to reproduce...
# analyze_dip(older_than='20150319_154856',newer_than='20150319_145019')


### with long mw switch risetime

#analyze_dip(older_than='20150319_183913',newer_than='20150319_182934')

### with phase gates around tau_larmor

# analyze_dip(older_than='20150323_140733',newer_than='20150323_135150')

##########################################################################
### for a single zeno measurement

# analyze_dip(older_than='20150325_191334',newer_than='20150325_172735')

### repeat the above measurement for a freshly calibrated system (parameters changed). 
analyze_dip(older_than='20150326_150115',newer_than='20150326_142500')


