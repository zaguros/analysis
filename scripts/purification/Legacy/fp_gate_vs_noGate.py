import numpy as np
'''
Script to analyze the dynamical decoupling data
Here we compare a fingerprint with and without active gate tuning of the NV optical transitions.
'''
import numpy as np
import os
# import analysis.lib.QEC.nuclear_spin_characterisation as SC #used for simulating FP response of spins
import matplotlib.cm as cm
from matplotlib import pyplot as plt
import nuclear_spin_char as SC;  reload(SC)


import fingerprint_funcs as fp_funcs; reload(fp_funcs)

def fingerprint(xlim =None):
    
    ###################
    ## Data location ##
    ###################
    
    timestamp ='20160112_140443'
    timestampGate ='20160112_143418'
    ssro_calib_folder = 'd:\\measuring\\data\\20160107\\172632_AdwinSSRO_SSROCalibration_Pippin_SIL1'
    a, folder = fp_funcs.load_mult_dat(timestamp,
            number_of_msmts = 25,
            x_axis_step     = 0.1,
            x_axis_start    = 3.5,
            x_axis_pts_per_msmnt= 51,
            ssro_calib_folder=ssro_calib_folder)

    b, folderGate = fp_funcs.load_mult_dat(timestampGate,
            number_of_msmts = 25,
            x_axis_step     = 0.1,
            x_axis_start    = 3.5,
            x_axis_pts_per_msmnt= 51,
            ssro_calib_folder=ssro_calib_folder)

 
    ###############
    ## Plotting ###
    ###############

    fig = a.default_fig(figsize=(10,5))
    ax = a.default_ax(fig)
    if xlim == None:
        ax.set_xlim(3.5,6)
    else:
        ax.set_xlim(xlim)
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 0.5))

    ax.set_ylim(-0.05,1.05)
   
    ax.plot(a.sweep_pts, a.p0, '.-k', lw=0.4,label = 'No gate tuning') #N = 16
    ax.plot(b.sweep_pts, b.p0, '.-b', lw=0.4,label = 'Gate tuning') #N = 16    

    plt.legend(loc=4)

    print folder
    plt.savefig(os.path.join(folder, 'fingerprint.pdf'),
        format='pdf')
    plt.savefig(os.path.join(folder, 'fingerprint.png'),
        format='png')

