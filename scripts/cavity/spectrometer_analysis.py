### spectrometer analysis 

import numpy as np 
import csv
import pandas as pd
from scipy.signal import argrelextrema

def load_data(filepath):
	"""
	Loading csv -files into pandas dataframe and convert to array
	Input: 
	filepath - the filepath of csv filepath
	Output:
	data_df - data as dataframe
	data - numpy ndarray with (Wavelength,Columns,Intensity)
	"""
	data_df = pd.read_csv(filepath, usecols=[2,4,5])
	data = data_df.as_matrix()
	return data_df, data



def peak_det(data, order):
	"""
	Detecting the peaks in the spectrometer FSR data
	Input:
	data - numpy ndarray with (Wavelength,Columns,Intensity)
	order - order of the peak detection method
	Output:
	indices - indices of the peaks
	peaks_WL - indices converted to wavelength
	peaks_f - indices converted to frequency 
	"""
	
	indices = argrelextrema(data[:,2], np.greater, order=order) # the number 100 is somewhat arbitrary, but seems to work. 

	peak_WL= [] #connect indices with the values for the wavelength, creating arrays
	peak_I=[]
	peak_freq=[]

	for i in indices[0]:
		peak_WL = np.append(peak_WL,data[i,0])
		peak_I = np.append(peak_I,data[i,2])
		peak_f = 3.e8/(data[i,0]*1.e-9)
		peak_freq = np.append(peak_freq,peak_f)

	return indices, peak_WL, peak_I, peak_freq


