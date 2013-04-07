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
def calc_meas_strength(x,t_zero,t_star):
    measstren=theta(x,t_zero,t_star)/90.
    return measstren

def theta(tau,t_zero,t_star):
    return 90-2*np.arccos(sqrt(S(tau,t_zero,t_star)))*180./np.pi

def S(tau,t_zero,t_star):
    return np.exp(-(tau/t_star)**2)*np.cos(np.pi/4-(tau+t_zero)*np.pi*.002185/2.)**2

#tau 50 ns, 4us RO 
date='20130401'
foldernamez_4 = '141914'
foldernamex_4 = '143456'

#tau 50 ns, 2 us RO
date='20130401'
foldernamez_2 = '140024'
foldernamex_2 = '133605'

#tau 50 ns, 6 us RO
date='20130401'
foldernamez_6 = '153222'
foldernamex_6 = '150233'

#tau 50 ns, 8 us RO
date='20130401'
foldernamez_8 = '180434'

#tau 50 ns, 15 us RO
date='20130401'
foldernamez_15 = '124303'

date2='20130404'
succesprob_fn=[foldernamez_2,foldernamez_4,foldernamez_6,foldernamez_8,foldernamez_15]
succesprob_fn2=['111620','113803','114334','110834','114951','115605','120436','120919','121458','122045','131359']
RO_times=[2,4,6,8,15,10,20,30,6,40,50,1,2,4,8,50]

def targetstate_analysis(foldernamez,foldernamex, filename='Spin_RO',date=''):
    zdata_norm,zdata_corr=sc.plot_feedback(foldernamez, filename='Spin_RO',d=date)
    xdata_norm,xdata_corr=sc.plot_feedback(foldernamex, filename='Spin_RO',d=date)
    SNfit=sc.fit_sin(xdata_norm['sweep_par'],xdata_corr['FinalRO_SN'],xdata_corr['uFinalRO_SN'])
    FSfit=sc.fit_sin(xdata_norm['sweep_par'],xdata_corr['FinalRO_FS'],xdata_corr['uFinalRO_FS'])
    allfit=sc.fit_sin(xdata_norm['sweep_par'],xdata_corr['FinalRO_All'],xdata_corr['uFinalRO_All'])
    
    xsucces= (xdata_norm['SN'][2]*abs(SNfit['params'][1])+xdata_norm['FS'][2]*abs(FSfit['params'][1]))/(xdata_norm['SN'][2]+xdata_norm['FS'][2])

    ysucces= (xdata_norm['SN'][2]*abs(0.5-SNfit['params'][2])+xdata_norm['FS'][2]*abs(0.5-FSfit['params'][2]))/(xdata_norm['SN'][2]+xdata_norm['FS'][2])


    Sz=[zdata_corr['FinalRO_All'][2],data_corr['FinalRO_Succes'][2],data_corr['FinalRO_SN'][2]]
    Sx=[abs(allfit['params'][1]),xsucces,abs(SNfit['params'][1])]
    Sy=[abs(0.5-allfit['params'][2]),ysucces,abs(0.5-SNfit['params'][2])]
    
    res_all=[Sz[0],Sx[0],Sy[0]]
    res_feedback=[Sz[1],Sx[1],Sy[1]]
    res_heralded=[Sz[2],Sx[2],Sy[2]]
    
    meas_strength = calc_meas_strength(50,12,1400)
    
    ideal=[np.sin(np.pi+meas_strength*np.pi/2.)/2.+0.5,np.cos(meas_strength*np.pi/2.)/2.,0]
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    width=0.2
    x=np.arange(3)

    barall=ax.bar(x, res_all, width, color='RoyalBlue')
    barfee=ax.bar(x+width, res_feedback, width, color='Crimson')
    barher=ax.bar(x+2*width, res_heralded, width, color='DarkGreen')
    baride=ax.bar(x+3*width, ideal, width, color='White')


    ax.set_title('Final State (tau=50ns,RO 6 us)')
    ax.set_xticks(x+width*2)
    ax.set_xticklabels( ('Sz','Sx', 'Sy') )
    ax.legend( (barall[0], barfee[0],barher[0],baride[0]), ('Total', 'Feedback','Heralded','Target') )
    plt.show()

    return res_all,res_feedback,res_heralded

def succesprob(folders,folderstwo,RO_times,date,date2):
    s_heralded=[]
    s_feedback=[]
    for fn in folders:
        zdata_norm,zdata_corr=sc.plot_feedback(fn, filename='Spin_RO',d=date)
        s_heralded.append(zdata_norm['SN'][2])
        s_feedback.append(zdata_norm['SN'][2]+(1-zdata_norm['SN'][2])*zdata_norm['FS'][2])
    for fn in folderstwo:
        zdata_norm,zdata_corr=sc.plot_feedback(fn, filename='Spin_RO',d=date2)
        s_heralded.append(zdata_norm['SN'][1])
        s_feedback.append(zdata_norm['SN'][1]+(1-zdata_norm['SN'][1])*zdata_norm['FS'][1])
    fig=plt.figure()
    ax=fig.add_subplot(111)

    ax.plot(RO_times,s_feedback,'o',color='Crimson',label='Feedback')
    ax.plot(RO_times,s_heralded,'o',color='RoyalBlue',label='Heralded')
    ax.set_ylabel('P(Succes)')
    ax.set_xlabel('RO duration (us)')
    ax.set_title('Probability to reach Target state')
    ax.set_yticks([0,0.25,0.5,0.75,1])
    ax.legend(loc=2)

total,feedback,heralded=targetstate_analysis(foldernamez_2,foldernamex_2)    
#succesprob(succesprob_fn,succesprob_fn2,RO_times,date,date2)
