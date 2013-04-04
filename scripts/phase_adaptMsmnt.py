
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




for n in arange (9):
    datapath=get_latest_data(foldername,date=d)
    
    files = os.listdir(datapath)
    f=filename
    for k in files:
        if f in k:
            spin_ro_file = k

    data = np.load(datapath+'\\'+spin_ro_file)
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




