from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot
from analysis.lib.tools import toolbox
from matplotlib import pyplot as plt
import numpy as np
import datetime
import time
import calendar
reload(plot)
import os


#practive data: D:\measuring\data\20160803\101635_temp_monitor
def get_temperature(older_than, newer_than):

	x = np.array([])
	y = np.array([])

	while toolbox.latest_data(contains='temp_monitor', older_than=older_than, newer_than=newer_than, raise_exc = False)!=False:
			## Find the data folder
			timestamp, folder = toolbox.latest_data(contains='temp_monitor', older_than=older_than, newer_than=newer_than,return_timestamp = True)
			

			## convert timestamp to Epoch time (absolute time)
			Time_temp1 = datetime.datetime(int(timestamp[:4]), int(timestamp[4:6]), int(timestamp[6:8]), 
										   int(timestamp[8:10]), int(timestamp[10:12]), int(timestamp[12:14]), 0)
			Time_temp1 = datetime.datetime.timetuple(Time_temp1)
			Epoch_timestamp = calendar.timegm(Time_temp1) # in seconds

			## get the data 
			filename 		= toolbox.measurement_filename(folder, ext='dat')
			data 			= np.loadtxt(filename, skiprows=16)
			
			# convert time data to Epoch time
			time_data = data[:,0]*3600 + Epoch_timestamp
			temperature_data = (data[:,1]-100)/0.385

			# Filter out fluke measurements
			indices_flukes      = [i for i,j in enumerate(temperature_data) if j<10]
			time_data 			= np.delete(time_data, indices_flukes)
			temperature_data 	= np.delete(temperature_data, indices_flukes)

			x = np.concatenate((x,time_data)) 				# Time after beginning of measurement [sec]
			y = np.concatenate((y,temperature_data))	
			
			older_than = str(int(timestamp)-1)
	
	## Ordering in time: 
	sort_indices = np.argsort(x)
	x  = np.array(x)[sort_indices]
	y  = np.array(y)[sort_indices]

	print 'folder = ' + folder

	return x,y, folder

def plot_temperature(From, To):

	## Get a significantly larger part of the data
	x,y, folder = get_temperature('20160900_000000','20160700_000000')
	# x,y, folder = get_temperature(To, From)

	##convert timestamp to Epoch time (absolute time)
	Time_temp1 = datetime.datetime(int(From[:4]), int(From[4:6]), int(From[6:8]), 
								   int(From[9:11]), int(From[11:13]), int(From[13:15]), 0)
	Time_temp2 = datetime.datetime(int(To[:4]), int(To[4:6]), int(To[6:8]), 
								   int(To[9:11]), int(To[11:13]), int(To[13:15]), 0)
		
	print Time_temp1
	Time_temp1 = datetime.datetime.timetuple(Time_temp1)
	Time_temp2 = datetime.datetime.timetuple(Time_temp2)
	
	print Time_temp1
	Epoch_From = (calendar.timegm(Time_temp1))/3600. # in hours
	Epoch_To   = (calendar.timegm(Time_temp2))/3600. # in hours

	## convert to hours
	x = x/3600.

	## Find the starting and end points
	index_start = (np.abs(x-Epoch_From)).argmin()
	index_end   = (np.abs(x-Epoch_To)).argmin()

	## Figures:

		### FIGURE temperature plot
	fig = plt.figure(1,figsize=(20,8))

	plt.plot(x[index_start:index_end] - x[index_start] , y[index_start:index_end] ,'.')
	plt.xlim((Epoch_From- x[index_start],Epoch_To- x[index_start]))
	plt.xlabel('time (hours) from ' + From)
	plt.title(folder)

	plt.savefig(os.path.join(folder, 'Temperature vs time.png'), format='png')

		### FIGURE gradient plot
	fig = plt.figure(2,figsize=(20,8))

	xgrad 	 = np.gradient(x)
	gradient = np.gradient(y,xgrad)

	plt.plot(x[index_start:index_end] - x[index_start] , gradient[index_start:index_end] ,'.')
	plt.plot( [0,x[index_end] - x[index_start]] , [0.1, 0.1] ,'r-')
	plt.plot( [0,x[index_end] - x[index_start]] , [-0.1, -0.1] ,'r-')
	
	plt.xlim((Epoch_From- x[index_start],Epoch_To- x[index_start]))
	plt.xlabel('time (hours) from ' + From)
	plt.ylabel('gradient [Kelvin/hour]')
	plt.title(folder)

	plt.savefig(os.path.join(folder, 'Temperature gradient.png'), format='png')

	# 	### FIGURE histograms
	# fig = plt.figure(3,figsize=(20,8))

	# plt.plot(x[index_start:index_end] - x[index_start] , gradient[index_start:index_end] ,'.')
	# plt.plot( [0,x[index_end] - x[index_start]] , [0.1, 0.1] ,'r-')
	# plt.plot( [0,x[index_end] - x[index_start]] , [-0.1, -0.1] ,'r-')
	
	# plt.xlim((Epoch_From- x[index_start],Epoch_To- x[index_start]))
	# plt.xlabel('time (hours) from ' + From)
	# plt.ylabel('gradient [Kelvin/hour]')
	# plt.title(folder)

	# plt.savefig(os.path.join(folder, 'Temperature histogram.png'),
	#          format='png')







