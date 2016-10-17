import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc, cm


class compare_functions ():

	def __init__(self):
		self.data = {}
		self.counter = 0 
		self.xlabel = ''
		self.ylabel = ''
		self.title = ''
		self.log_plot = False
		
	def add (self, x, y, label):
		self.counter = self.counter+1
		self.data['x_'+str(self.counter)]=x
		self.data['y_'+str(self.counter)]=y
		self.data['l_'+str(self.counter)]=label
	
	def plot (self, numbers = None,y_log=False):
			
		if (numbers==None):
			numbers = np.arange(self.counter)+1
		fig = plt.figure(figsize=(12,8))
		p = fig.add_subplot(1,1,1)

		colors = cm.jet(np.linspace(0., 0.9, len(numbers)))
		ind = 0 
		for i in numbers:	
			if self.log_plot:	
				plt.semilogy (self.data['x_'+str(i)], self.data['y_'+str(i)], label = self.data['l_'+str(i)], color = colors[ind]) 
			else:
				plt.plot (self.data['x_'+str(i)], self.data['y_'+str(i)], ':', color = colors[ind]) 
				p.plot (self.data['x_'+str(i)], self.data['y_'+str(i)], 'o', color = colors[ind], label = self.data['l_'+str(i)]) 
			ind = ind + 1
		if y_log:
			plt.yscale('log')
		#x0 = self.data['x_1']
		#y0 = self.data['y_1']
		#y = y0[0]/(x0/x0[0])

		#plt.plot (x0, y, ':k')	
		p.tick_params(axis='both', which='major', labelsize=15)
		plt.xlabel (self.xlabel, fontsize = 20)
		plt.ylabel (self.ylabel, fontsize = 20)
		plt.title (self.title)
		plt.legend()
		plt.show()
		return fig

