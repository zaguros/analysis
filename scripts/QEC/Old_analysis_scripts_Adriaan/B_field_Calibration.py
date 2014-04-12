import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## Import data
df = pd.read_excel('magetCalibr.xlsx','Sheet1')

diameter = df['dimensions(mm)'][0]
height = df['dimensions(mm)'][1]
Br =  df['Brem'][0]
#http://www.supermagnete.de/eng/data_table.php
R = diameter /2
D = height
z = df.ix[df['x(mm)']==0]['z(mm)'].values
B_Calc = Br/2 *((D+z)/(np.sqrt(R**2+(D+z)**2))-z/(np.sqrt(R**2 +z**2)))

B_fields = df.ix[df['x(mm)']==0]['B(G)'].values
uB_fields = .2 *B_fields#data[u'uB (G)'].values
label = df['MagnetLabel'][0]

offset = 4.7


df = pd.read_excel('magetCalibr.xlsx','Sheet2')

diameter = df['dimensions(mm)'][0]
height = df['dimensions(mm)'][1]
Br =  df['Brem'][0]
#http://www.supermagnete.de/eng/data_table.php
R = diameter /2
D = height
z2 = df.ix[df['x(mm)']==0]['z(mm)'].values
B_Calc2 = Br/2 *((D+z2)/(np.sqrt(R**2+(D+z2)**2))-z2/(np.sqrt(R**2 +z2**2)))
label2 = df['MagnetLabel'][0]


B_fields2 = df.ix[df['x(mm)']==0]['B(G)'].values
uB_fields2 = .2 *B_fields2#data[u'uB (G)'].values



plt.plot(z,B_Calc,label = str('magnet ' + label + ' calculated'))
plt.errorbar(z+offset,B_fields,uB_fields,label = str('magnet ' + label +' measured' ))
plt.plot(z2,B_Calc2, label =str('magnet ' + label2+ ' calculated'))
plt.errorbar(z2+offset,B_fields2,uB_fields2,label=str('magnet ' +label2+ 'measured'))
plt.title(r'Magnetic Field Calibrations: '+ 'offset in measurements = ' + str(offset)+'mm')
plt.ylabel(r'B [Gauss]' )
plt.xlabel('z [mm]')
plt.xlim([z[0]-1,z[-1]+4])
plt.legend(loc= 'upper right')


plt.show()
