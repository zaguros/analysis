import os, sys
import numpy as np
import h5py
import logging
import csv
import pylab as plt

class CVSTrace():

	def __init__(self):
		self.trace = None
		self.x_axis = None

	def load_trace(self, filename, delimiter = ',', quotechar='|', scope_type = None):

		with open(filename, 'rb') as csvfile:
			row_list = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)
			ind = 0
			x = []
			y = []
			for row in row_list:
				ind = ind + 1
				try:
					if (scope_type=='yoko'):
						y.append (float(row[0]))
						print row[0], row[1]
					else:
						x.append(float(row[0]))
						y.append(float(row[1]))
				except:
					pass
					#print 'row ', ind, ' is not a number'

		self.trace = np.asarray(y)
		self.x_axis = np.asarray (x)

	def load_trace_yoko (self, filename, delimiter = ',', quotechar='|'):

		with open(filename, 'rb') as csvfile:
			row_list = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)
			ind = 0
			x = []
			y = []
			for row in row_list:
				try:
					y.append( float(row[1]))
				except:
					print 'no value'
				if (ind<15):
					print ind, row
				if (ind==8):
					self.h_res = float(row[1])
				if (ind==7):
					self.sampling_rate = float(row[1])
				ind = ind+1
		self.trace = np.asarray(y)
		self.x_axis = np.ones(len(self.trace))*self.h_res

	def plot_trace (self):
		plt.figure(figsize = (10,6))
		plt.plot (self.x_axis, self.trace, '.b')
		plt.plot (self.x_axis, self.trace, 'r')
		plt.show()

