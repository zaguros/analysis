
import numpy as np
import pylab as plt
from analysis.lib.spin import spin_control as sc
import os

def get_latest_data(string = 'ADwin_SSRO', datapath = '',date=''):
    meas_folder = r'D:\measuring\data'
    if date=='':
        currdate = time.strftime('%Y%m%d')
    else:
        currdate=date
    
    if datapath == '':
        df = os.path.join(meas_folder, currdate)
    else:
        df = datapath
    
    right_dirs = list()

    if os.path.isdir(df):
        for k in os.listdir(df):
            if string in k:
                right_dirs.append(k)
        
        if len(right_dirs) > 0:
            latest_dir = os.path.join(df,right_dirs[len(right_dirs)-1])
        else:
            print 'No measurements containing %s in %s'%(string, df)
            latest_dir=''
        
        print '\nAnalyzing data in %s'%latest_dir

    else:
        print 'Folder %s does not exist'%df
        latest_dir = False

    return latest_dir

# 2 msmnts uncollapse very short
d='20130404'
foldernames=['192943','201204','212222','221836','233234','005033','015405','030325','035504']

# 1 msmnt
d='20130405'
foldernames=['141755','144145','150814','153230','155539','161924','164248','171312']
for n in foldernames:
    datapath=sc.get_latest_data(n,date=d)
    
    files = os.listdir(datapath)
    f='Spin_RO'
    for k in files:
        if f in k:
            spin_ro_file = k

    data = np.load(datapath+'\\'+spin_ro_file)
    if (n==foldernames[0]):
        phase=np.zeros(len(data['FS']))
        rep=np.zeros(len(data['FS']))
    reps_array=data['SN']+data['FF']+data['FS']
    reps=reps_array[0]
    data_norm={}
    data_norm['sweep_par']=data['sweep_par']
    data_norm['sweep_par_name']=data['sweep_par_name']
    data_norm['SN']=data['SN']/(reps+0.)
    data_norm['FS']=data['FS']/(reps-data['SN']+0.)
    data_norm['FF']=data['FF']/(reps-data['SN']+0.)
    data_norm['uSN']=(data_norm['SN']*(1-data_norm['SN'])/reps)**.5
    data_norm['uFS']=(data_norm['FS']*(1-data_norm['FS'])/(reps-data['SN']+0.))**.5
    data_norm['uFF']=(data_norm['FF']*(1-data_norm['FF'])/(reps-data['SN']+0.))**.5

    data_norm['FinalRO_SN']=data['FinalRO_SN']/(data['SN']+0.)
    data_norm['FinalRO_FS']=data['FinalRO_FS']/(data['FS']+0.)
    data_norm['FinalRO_FF']=data['FinalRO_FF']/(data['FF']+0.)
    data_norm['FinalRO_Succes']=(data['FinalRO_SN']+data['FinalRO_FS'])/(data['SN']+data['FS']+0.)
    data_norm['FinalRO_All']=(data['FinalRO_SN']+data['FinalRO_FS']+data['FinalRO_FF'])/(reps+0.)

    print data['FS']
    for i in np.arange(len(data['SN'])):
        phase[i]+=data['FinalRO_SN'][i]
        rep[i]+=data['SN'][i]
    print rep
    phasenorm=phase/(rep+0.)
    x=data['sweep_par']
    data.close()

phasecor,uphasecor=sc.get_nuclear_ROC(phasenorm,rep,sc.get_latest_data('SSRO',date=d))
figure1=plt.figure()
ax=figure1.add_subplot(111)
ax.errorbar(x,phasenorm,fmt='o',yerr=1/sqrt(rep),color='Crimson')
ax.errorbar(x,phasecor,fmt='o',yerr=uphasecor,color='RoyalBlue')
sc.fit_sin(x,phasecor,uphasecor)

