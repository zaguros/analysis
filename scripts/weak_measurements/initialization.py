import os
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D


fres_desr,data_desr,udata_desr = sc.plot_dark_esr(sc.get_latest_data('dark',date='20130311'),d='20130311')
fres_zero,data_zero,udata_zero = sc.plot_dark_esr(sc.get_latest_data('init_',date='20130311'),d='20130311')
fres_min1,data_min1,udata_min1 = sc.plot_dark_esr(sc.get_latest_data('init_',date='20130305'),d='20130305')

figure1=plt.figure(1)
figure1.clf()
ax=figure1.add_subplot(111)
ax.errorbar(fres_desr['x'],data_desr,yerr=udata_desr,color='LightGrey',fmt='o',label='No initialization')
xfit=np.linspace(fres_desr['x'].min(),fres_desr['x'].max(),201)
ax.plot(xfit,fres_desr['fitfunc'](xfit),color='LightGrey',ls='-')

ax.errorbar(fres_zero['x'],data_zero,yerr=udata_zero,color='RoyalBlue',fmt='o',label='mI=0')
xfit=np.linspace(fres_zero['x'].min(),fres_zero['x'].max(),201)
ax.plot(xfit,fres_zero['fitfunc'](xfit),color='RoyalBlue',ls='-')

ax.errorbar(fres_min1['x'],data_min1,yerr=udata_desr,color='Crimson',fmt='o',label='mI=-1')
xfit=np.linspace(fres_min1['x'].min(),fres_min1['x'].max(),201)
ax.plot(xfit,fres_min1['fitfunc'](xfit),color='Crimson',ls='-')

ax.set_xlabel('MW freq (GHz)')
ax.set_ylabel('P(ms=0)')
ax.set_ylim([0,1.2])
ax.set_xlim([2.825,2.8335])
ax.set_yticks([0,0.5,1])
ax.set_xticks([2.8245,2.828,2.830,2.832])
ax.xaxis.set_major_formatter(FormatStrFormatter('%0.3f'))
ax.set_title('Dark ESR')
ax.legend(loc=1,prop={'size':10})
