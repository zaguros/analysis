import os
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common,ramsey
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc
from matplotlib import rc

#datafolders=['1154','1252','1258','1303','1310','1347','1351','1316','1326', '1453']
#RO_time=[0,1,2,3,4,5,6,7,9, 11]
#
#
#datafolders=['164825','165126']
#date='20130111'

#datafolders=['111709','112023']
#date='20130115'

datafolders=['183118','183446']
date='20130116'
label=['mI=0','mI=+1']
meas_folder=r'D:\measuring\data'

datafolders=['190058','190909']
date='20130126'
label=['mI=-1','mI=0']
meas_folder=r'D:\measuring\data'

dp=os.path.join(meas_folder, date)
amp=[]
phase=[]
plt.close('all')
plt.figure(75)
j=0
for i in datafolders:    
    tau_guess = 3000
    offset_guess=0.5   
    freq=1e-3
    amp=0.5
    result=sc.plot_data_MBI(sc.get_latest_data(string=i,datapath=dp),fid=(0.8125,0.991900))
    print result['yerr']
    plt.figure(75)
    fit_result = fit.fit1d(result['x'], result['y'], ramsey.fit_ramsey_gaussian_decay, 
                tau_guess, offset_guess, (freq,amp,90),
                do_print = True , ret = True)
    #plot.plot_fit1d(fit_result,np.linspace(result['x'].min(),result['x'].max(),201))
    #measstrent=(90*result['x']/229 +0)/90
    #t=result['x']
    #result['x']=measstrent
    #result['x']=result['x']/191
    plt.errorbar(result['x']/1000.,result['y'],fmt='o',yerr=result['yerr'],label=label[j])
    fit_x=np.linspace(result['x'].min(),result['x'].max(),201)
    plt.errorbar(fit_x/1000.,fit_result['fitfunc'](fit_x),fmt='-r')
    #yerr_weak=np.sqrt(2*(1+np.cos(result['y']))/500.)/sin(result['x']*90*pi/180)
    #plt.errorbar(result['x'],2*(result['y']-0.5)/np.sin(result['x']*90*pi/180),fmt='o',yerr=yerr_weak,label=label[j])
    j=j+1
plt.xlabel (' free evolution time [us]', fontsize = 16)
plt.ylabel ('P (ms=0)', fontsize = 16)   
plt.ylim ([0, 1])
plt.xlim ([0, 3])
plt.legend(loc=1)
plt.show()


