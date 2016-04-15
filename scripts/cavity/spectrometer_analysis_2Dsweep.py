#!/usr/bin/python

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob


indir="K:/ns\qt\Diamond\Projects\Cavities\data/20160330 Spectrometer scripts and data/20160331 Sweep2Dplot reversed" 
outdir="C:\Users\lcoenen\Dropbox\Afstuderen Diamond\Plaatjes"
#os.chdir(indir)

n= 121 #number of files

# load the files in the dataframe, only pick out Intensity and Column. All files in folder. 
dataframes = [pd.read_csv(filename, usecols = [4,5]) for filename in glob.glob(indir + "/*.csv")]
# load the file in the dataframe with a certain name in a certain folder
#dataframes = [pd.read_csv(os.path.join(indir,"2016-03-30 17_47_28 Zoki_pos7_Q500_it20000_600_700 %s.csv") % i, usecols=[4,5]) for i in xrange(1,2*n,2)]

# Group the Intensity data by column and compute the mean
mean_val=[]
for i in dataframes:
	merged=i.groupby('Column')
	mean=merged.agg([np.mean])
	mean_val.append(mean)

# Concatenate in one dataframe
mean_val=pd.concat(mean_val,axis=1)
mean_val = mean_val.rename(columns={"Column":'Frequency'}) 

# Plot using seaborn function

#mean_val_array=mean_val.values
#ax = plt.subplots()
#plt.xlabel("Voltage (V)", fontsize = 18, fontweight='bold')
#plt.ylabel("Frequency THz", fontsize = 18, fontweight='bold')
ax=sns.heatmap(mean_val, xticklabels=False, yticklabels=False, vmin=0, vmax=300, cmap='YlGnBu')
#ax.xaxis.set_ticks_position('top')  # put column labels at the top
#ax.yaxis.set_ylim(-.5,1.5)
#ax.grid(b=True,which='major')
ax.set_title('2D plot cavity with diamond', fontsize = 16, fontweight='bold')
plt.xlabel("Voltage (V) (4.2 - 4.4 V)", fontsize = 14, fontweight='bold')
plt.ylabel("Frequency (THz)", fontsize = 14, fontweight = 'bold')
#xticks=np.arange(-2, 10, 1.0)
#yticks=np.arange(600,700,50)
#ax.set_yticks(yticks)
# ax.set_xticks(xticks)
# ax.set_xticklabels(np.arange(-2, 10, 1))
# ax.set_ytickslabels(yticks)
#ax.set_xticklabels(xticks,rotation=90)
#ax.set_yticklabels(yticks,rotation=90)

#ax.set_xticklabels(cols, rotation=45, ha='right')
plt.show()
#plt.savefig(os.path.join(outdir, "2Dplot.png"), format="png")





















