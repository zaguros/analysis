
from matplotlib import pyplot as plt

class compare_functions ():

	def __init__(self):
		self.data = {}
		self.counter = 0 
		self.xlabel = ''
		self.ylabel = ''
		self.title = ''
		self.log_plot = False
		
	def add (self, x, y, legend):
		self.counter = self.counter+1
		self.data['x_'+str(self.counter)]=x
		self.data['y_'+str(self.counter)]=y
		self.data['l_'+str(self.counter)]=legend
	
	def plot (self, numbers = None):
			
		if (numbers==None):
			numbers = np.arange(self.counter)+1
		
		colors = cm.gist_heat(np.linspace(0., 0.8, len(numbers)))

		ind = 0 
		for i in numbers:	
			if self.log_plot:	
				plt.loglog (self.data['x_'+str(i)], self.data['y_'+str(i)], label = self.data['l_'+str(i)], color = colors[ind]) 
			else:
				plt.plot (self.data['x_'+str(i)], self.data['y_'+str(i)], label = self.data['l_'+str(i)], color = colors[ind]) 
			ind = ind + 1

		x0 = self.data['x_1']
		y0 = self.data['y_1']
		y = y0[0]/(x0/x0[0])
		
		plt.plot (x0, y, ':k')	
		plt.xlabel (self.xlabel)
		plt.ylabel (self.ylabel)
		plt.title (self.title)
		plt.legend()
		plt.show()

