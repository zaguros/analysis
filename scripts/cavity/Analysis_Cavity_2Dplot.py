#!/usr/bin/python

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

indir="C:\Users\lcoenen\Dropbox\Afstuderen Diamond\Data\Measurements on diamond/" 
#indir="C:\Users\lcoenen\Dropbox\Afstuderen Diamond\Data\Measurement on diamond/" 
outdir="C:\Users\lcoenen\Dropbox\Afstuderen Diamond\Data\Measurements on diamond/" 
#os.chdir(indir)

n= 12 #number of files

# load the files in the dataframe, only pick out Intensity and Column
dataframes = [pd.read_csv(os.path.join(indir,"Zoki_pos7_Q500_it20000_  %s.csv") % i, usecols=[4,5]) for i in xrange(1,n)]

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

#plt.figure()
#plt.xlabel(fontsize = 18, fontweight='bold')
#plt.ylabel(fontsize = 18, fontweight='bold')
ax=sns.heatmap(mean_val, xticklabels=True, yticklabels=True, vmin=0, vmax=3000, cmap='YlGnBu')
#ax.xaxis.set_ticks_position('top')  # put column labels at the top
#ax.yaxis.set_ylim(-.5,1.5)
#ax.grid(b=True,which='major')
ax.set_title('2D plot cavity with diamond')
#ax.set_xticks(pos + 0.4)
yticks=np.arange(600,700,50)
xticks=np.arange(-2,10,1)
ax.set_yticks(yticks)
ax.set_xticks(xticks)
ax.set_xticklabels(xticks,rotation=90)
ax.set_yticklabels(yticks,rotation=90)


#ax.set_xticklabels(cols, rotation=45, ha='right')
plt.show()
plt.savefig(os.path.join(outdir, "2Dplot.eps"), format="eps")





















