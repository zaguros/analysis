
from numpy import *
import pylab as plt
import numpy as np
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc

#datafolders=['1154','1252','1258','1303','1310','1347','1351','1316','1326', '1453']
#RO_time=[0,1,2,3,4,5,6,7,9, 11]
SIL9=[61.8,66.6,r'D:\measuring\data\20130627\122204_LaserFrequencyScan_red_scan_coarse_gv_0.0\122204_LaserFrequencyScan_red_scan_coarse_gv_0.0.npz',9]
SIL10=[56.9,57.3,r'D:\measuring\data\20130701\140903_LaserFrequencyScan_SIL10_red_mw_gv_0.0\140903_LaserFrequencyScan_SIL10_red_mw_gv_0.0.npz',10]
SIL11=[53.62,55.61,r'D:\measuring\data\20130702\093927_LaserFrequencyScan_SIL11_green_gv_0.0\093927_LaserFrequencyScan_SIL11_green_gv_0.0.npz',11]
SIL18=[55.6,59.5,r'D:\measuring\data\20130702\102752_LaserFrequencyScan_SIL18_green_gv_0.0\102752_LaserFrequencyScan_SIL18_green_gv_0.0.npz',18]
SIL15=[59.6,63.2,r'D:\measuring\data\20130702\112428_LaserFrequencyScan_SIL15_Green_gv_0.0\112428_LaserFrequencyScan_SIL15_Green_gv_0.0.npz',15]
SIL12=[0,0,r'D:\measuring\data\20130702\175735_LaserFrequencyScan_SIL12_green_gv_0.0\175735_LaserFrequencyScan_SIL12_green_gv_0.0.npz',12]
SIL8=[55.5,60.4,r'D:\measuring\data\20130703\132950_LaserFrequencyScan_SIL8_red_gv_0.0\132950_LaserFrequencyScan_SIL8_red_gv_0.0.npz',8]
SIL7=[52.7,59.5,r'D:\measuring\data\20130702\122103_LaserFrequencyScan_SIL7_Green_gv_0.0\122103_LaserFrequencyScan_SIL7_Green_gv_0.0.npz',7]
SIL5=[40.8,50.3,r'D:\measuring\data\20130702\173219_LaserFrequencyScan_SIL5_green_gv_0.0\173219_LaserFrequencyScan_SIL5_green_gv_0.0.npz',5]
SIL1 = [71,78.4,r'D:\measuring\data\20130702\171520_LaserFrequencyScan_SIL1_red_gv_0.0\171520_LaserFrequencyScan_SIL1_red_gv_0.0.npz',1]
datafolders=[SIL1,SIL5,SIL7,SIL8,SIL9,SIL10,SIL11,SIL15,SIL18]
plt.figure(9)
j=0
n=38.5
for i in datafolders:
    
    result = np.load(i[2])
    r=result['data']
    plt.plot(r[:,1],r[:,2]/float(max(r[:,2]))+j)
    if (i[0]!=0):
        plt.text(i[0]-1,j+0.5,'Ey',fontsize=8)
    if (i[1]!=0):
        plt.text(i[1]+0.3,j+0.5,'Ex',fontsize=8)
    name='SIL'+str(i[3])
    plt.text(n,j+0.1,name,fontsize=10)
    j=j+1
    result.close()
plt.xlim([38,80])
plt.ylim([-0.1,j])
plt.xlabel ('Relative Frequency [GHz]')   
plt.ylabel ('Counts (Normalized)')   
'''
datafoldersCond=['1502','4us','6us','8us', '10us', '12us']
RO_timeCond=[2,4,6,8, 10, 12]
ampCond=[]
phaseCond=[]
for i in datafoldersCond:
    result= sc.plot_rabi(sc.get_latest_data(i,datapath=r'D:\measuring\data\20121219'))
    #qt.sleep(4)
    ampCond.append(abs(2*result[0]['params'][1]))
    phaseCond.append(result[0]['params'][3])


result=fit.fit1d(np.array(RO_power), np.array(amp/amp[0]), common.fit_exp_decay_with_offset, 
            0, 1, 5,
            do_plot = True, do_print = True, newfig = False,ret=True)

#resultcond=fit.fit1d(np.array(RO_timeCond), np.array(ampCond/amp[0]), common.fit_exp_decay_with_offset, 
#            0, 1, 5,
#            do_plot = True, do_print = True, newfig = False,ret=True)
print 'normal RO'
print result[0]['params']
#print 'cond RO'
#print resultcond[0]['params']


plt.plot(RO_power,amp/amp[0],'bo')
plt.plot(result[0]['fitx'],result[0]['fity'],'b-')
#plt.plot(resultcond[0]['fitx'],resultcond[0]['fity'],'r-')
#plt.plot (RO_timeCond, ampCond/amp[0], 'r<')
plt.xlabel ('RO time [us]', fontsize = 16)
plt.ylabel ('Contrast', fontsize = 16)   
plt.ylim ([0, 1])
plt.show()

plt.figure(11)
plt.plot(RO_power,mod(phase,180),'bo')
#plt.plot (RO_timeCond, mod(phaseCond, 360), 'r<')
plt.xlabel ('RO Power [nW]', fontsize = 16)
plt.ylabel ('Phase [degree]', fontsize = 16)   
plt.ylim ([0, 360])
plt.show()
'''

