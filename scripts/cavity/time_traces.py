from analysis.lib.tools.oscilloscope_cvs_trace import CVSTrace as Trace
import numpy as np
import pylab as plt
from matplotlib import rc, cm

matplotlib.rc ('xtick', labelsize=15)
matplotlib.rc ('ytick', labelsize=15)

idx = 4
fname = r"D:\Research\cavity\room_temperature\NewFile"+str(idx)+".csv"

t = Trace()
t.load_trace(filename=fname)
x = t.x_axis
y = t.trace

plt.figure (figsize=(17, 5))
plt.plot (x*1000, y, '.b')
plt.plot (x*1000,y,'r')
plt.xlabel ('time [ms]', fontsize=15)
plt.ylabel ('photodiode signal', fontsize=15)
#plt.xlim ([48, 65])
plt.show()

