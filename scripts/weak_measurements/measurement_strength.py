import os
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D

#datafolders=['1154','1252','1258','1303','1310','1347','1351','1316','1326', '1453']
#RO_time=[0,1,2,3,4,5,6,7,9, 11]
#
#
#datafolders=['164825','165126']
#date='20130111'

#datafolders=['111709','112023']
#date='20130115'
def make_rho(z,x):

    rho=np.array([[x,z],[1-z,x]])
    return rho
def make_hist(data,title=''):

    column_names = ['','']
    row_names = ['','']

    fig = plt.figure()
    ax = Axes3D(fig)

    lx= len(data[0])            # Work out matrix dimensions
    ly= len(data[:,0])
    xpos = np.arange(0,lx,1)    # Set up a mesh of positions
    ypos = np.arange(0,ly,1)
    xpos, ypos = np.meshgrid(xpos+0.25, ypos+0.25)

    xpos = xpos.flatten()   # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx*ly)

    dx = 0.75 * np.ones_like(zpos)
    dy = dx.copy()
    dz = data.flatten()

    ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color='RoyalBlue',alpha=0.75)

#sh()
    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.set_zlim3d([0,1])
    ax.set_zlabel('')
    ax.set_title(title)
    plt.show()
def calc_meas_strength(x,t_max=210.,theta_min=4.):
    measstren=(90*x/t_max +theta_min)/90
    return measstren
'''
corr=(.728,1-.088)
datafolders=['175338','175003']
date='20130225'
label=['mI=0','mI=+1']
meas_folder=r'D:\measuring\data'

dp=os.path.join(meas_folder, date)
dphase=[]
plt.close('all')
plt.figure(75)
j=0
for i in datafolders:    
    
    result=sc.plot_data_MBI(sc.get_latest_data(string=i,datapath=dp),fid=corr)
    plt.figure(75)
    measstrent=(90*result['x']/214 +2)/90
    t=result['x']
    result['x']=measstrent
    #result['x']=result['x']/191
    plt.errorbar(result['x'],2*(result['y']-0.5),fmt='o',yerr=result['yerr'],label=label[j])
    #plt.errorbar(result['x'],result['y'],fmt='o',yerr=result['yerr'],label=label[j])
    yerr_weak=np.sqrt(2*(1+np.cos(result['y']))/500.)/sin(result['x']*90*pi/180)
    #plt.errorbar(result['x'],2*(result['y']-0.5)/np.sin(result['x']*90*pi/180),fmt='o',yerr=yerr_weak,label=label[j])
    j=j+1
plt.plot(np.linspace(0,1,100),np.exp(-(np.linspace(0,0,100)/1300)**2)*np.sin(np.linspace(0,1,100)*np.pi/2),'r-')
plt.plot(np.linspace(0,1,100),np.exp(-(np.linspace(0,0,100)/1300)**2)*np.sin(np.pi+np.linspace(0,1,100)*np.pi/2),'r-')

plt.xlabel (' Measurement Strength (AU)', fontsize = 16)
plt.ylabel ('< $\sigma_{Z}$ > $_{nitrogen}$', fontsize = 16)   
plt.ylabel ('| $\uparrow$ > $_{nitrogen}$', fontsize = 16)   
plt.ylim ([-1, 1])
plt.xlim ([-0.01, 1])
plt.legend(loc=2)
plt.show()

datafolders=['_25ns','_125ns','_225ns']
date='20130118'
dp=os.path.join(meas_folder, date)
label=['mI=0','mI=+1']
meas_folder=r'D:\measuring\data'
'''

date='20130227'
zfolder='223541'
yfolder='232520'
xfolder='002126'
#dp=os.path.join(meas_folder, date)
result_zmeas=sc.analyse_plot_results_vs_sweepparam(zfolder,yname='P(mI=0)',Nuclcor=True,dataname='Spin_RO',title='strong_meas_res_Z',d=date)
result_zcond=sc.analyse_weakcond_vs_sweepparam(zfolder,yname='P(mI=0)',Nuclcor=True,dataname='cond_Spin_RO',title='',d=date)
date='20130304'
result_ymeas=sc.analyse_plot_results_vs_sweepparam(yfolder,yname='P(mI=0)',Nuclcor=True,dataname='Spin_RO',title='strong_meas_res_y',d=date)
result_ycond=sc.analyse_weakcond_vs_sweepparam(yfolder,yname='P(mI=0)',Nuclcor=True,dataname='cond_Spin_RO',title='',d=date)
date='20130305'
result_xmeas=sc.analyse_plot_results_vs_sweepparam(xfolder,yname='P(mI=0)',Nuclcor=True,dataname='Spin_RO',title='strong_meas_res_x',d=date)
result_xcond=sc.analyse_weakcond_vs_sweepparam(xfolder,yname='P(mI=0)',Nuclcor=True,dataname='cond_Spin_RO',title='',d=date)

zmeas_x=calc_meas_strength(result_zmeas['x'],210,4)
xmeas_x=calc_meas_strength(result_xmeas['x'],210,4)
def plot_backaction():
    figure42=plt.figure(42)
    #plt.figure(42)
    print 'plotting backaction'
    plt.clf()
    plt.errorbar(zmeas_x,result_zcond['y_cond'],fmt='o', yerr=result_zcond['uy_cond'],label='rho 11')
    plt.errorbar(zmeas_x,1-result_zcond['y_cond'],fmt='o', yerr=result_zcond['uy_cond'],label='rho 00')
    plt.errorbar(xmeas_x,result_xcond['y_cond']-0.5,fmt='o',yerr=result_xcond['uy_cond'],label='<Sx>')
    plt.errorbar(xmeas_x,result_ycond['y_cond']-.5,fmt='o',yerr=result_ycond['uy_cond'],label='<Sy>')
    plt.plot(np.linspace(0,1,100),(np.cos(np.linspace(0,1,100)*np.pi/2.)/2.),'r-')
    plt.plot(np.linspace(0,1,100),(np.sin(np.pi+np.linspace(0,1,100)*np.pi/2.)/2.+0.5),'b-')
    plt.plot(np.linspace(0,1,100),(np.sin(np.linspace(0,1,100)*np.pi/2.)/2.+0.5),'g-')
    plt.xlabel (' Measurement Strength (a.u.)', fontsize = 16)
    plt.ylabel ('Density Matrix Element', fontsize = 16)
    plt.title('Backaction conditioned on weak measurement',fontsize=16)
    plt.ylim ([-0.1, 1])
    plt.xlim ([0, 1])
    plt.legend(loc=2)
    plt.show()

    figure43=plt.figure(43)
    #plt.figure(42)
    print 'plotting backaction'
    plt.clf()
    plt.errorbar(zmeas_x,result_zmeas['y'],fmt='o', yerr=result_zmeas['uy'],label='rho 11')
    plt.errorbar(zmeas_x,1-result_zmeas['y'],fmt='o', yerr=result_zmeas['uy'],label='rho 00')
    plt.errorbar(xmeas_x,result_xmeas['y']-0.5,fmt='o',yerr=result_xmeas['uy'],label='<Sx>')
    plt.errorbar(xmeas_x,result_ymeas['y']-0.5,fmt='o',yerr=result_ymeas['uy'],label='<Sy>')
    plt.plot(np.linspace(0,1,100),(np.cos(np.linspace(0,1,100)*np.pi/2.)/2.),'r-')
    plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100)+0.5,'b-')
    plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100)+0.5,'g-')
    plt.xlabel (' Measurement Strength (a.u.)', fontsize = 16)
    plt.ylabel ('Density Matrix Element', fontsize = 16)
    plt.title('Backaction (unconditioned)',fontsize=16)
    plt.ylim ([-.1, 1])
    plt.xlim ([0, 1])
    plt.legend(loc=2)
    plt.show()
plot_backaction()    
'''
result_25ns=sc.plot_rabi(sc.get_latest_data(datafolders[0],datapath=dp),fid=corr)
result_125ns=sc.plot_rabi(sc.get_latest_data(datafolders[1],datapath=dp),fid=corr)
result_225ns=sc.plot_rabi(sc.get_latest_data(datafolders[2],datapath=dp),fid=corr)

date='20130116'
dp=os.path.join(meas_folder, date)
result_zmeas=sc.plot_data_MBI(sc.get_latest_data('180000',datapath=dp),fid=corr)
rho_25 = make_rho(result_zmeas['y'][3],abs(result_25ns[0]['params'][1]))
rho_125 =make_rho(result_zmeas['y'][11],abs(result_125ns[0]['params'][1]))
rho_225 =make_rho(result_zmeas['y'][-1],abs(result_225ns[0]['params'][1]))

rho_25_cond = make_rho(result_zmeas['y_cond'][3],abs(result_25ns[0]['params'][1]))
rho_125_cond =make_rho(result_zmeas['y_cond'][11],abs(result_125ns[0]['params'][1]))
rho_225_cond =make_rho(result_zmeas['y_cond'][-1],abs(result_225ns[0]['params'][1]))

make_hist(rho_25_cond,title='tau = 25 ns')
make_hist(rho_125_cond,title='tau = 125 ns')
make_hist(rho_225_cond,title='tau = 225 ns')

make_hist(rho_25,title='tau = 25 ns')
make_hist(rho_125,title='tau = 125 ns')
make_hist(rho_225,title='tau = 225 ns')
'''
