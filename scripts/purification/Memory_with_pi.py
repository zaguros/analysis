"""
This script analyses quantum memory measurements that were performed on LT4
Nov 2017 NK
"""


import numpy as np
import os
from analysis.lib.tools import toolbox, plot;
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
# import matplotlib as mpl
from analysis.lib.fitting import fit, common
import copy as cp
import analysis.lib.purification.purify_ssro as pu_ssro;reload(pu_ssro)

reload(fit);reload(mbi);reload(common);reload(toolbox)


#########################
""" global parameters """
#########################

carbon = 3 # can be 3 or 6
analyzed_powers = [6.0,4.0,2.0,1.0,0.5,0.1] ### in uW
get_z = True

base_f2 = None#r"M:\tnw\ns\qt\Diamond\Eigenpapers\17_QMemories2\Data\NoMWErrorsSweepTime"
kws = {'folder': base_f2,'ret':True,'do_plot':False,'do_fit':True,'show_guess':False,
       'do_print':False,'fixed':[0,2,5,6],'return_all':True}
       #### fit function is designed to fit A*Exp(-(t/T)^n)
if carbon == 3:
	older_than = '20171105_215538'
	newer_than = '20171104_122414'
	kws.update({'older_than':older_than})
	kws.update({'newer_than':newer_than})
elif carbon == 6:
	older_than = 0 # not done yet
	newer_than = 0 # not done yet
	kws.update({'older_than':older_than})
	kws.update({'newer_than':newer_than})

#########################
""" helper functions  """
#########################

### plotting and stuff
def plot_fit_vals(fit_dict,vals):
    fig,ax = plt.subplots()
    for val in vals:
        for key in fit_dict.keys():
            val_res,val_res_u,val_x = [],[],[] ## init empty
            for f in fit_dict[key]:
                if f['error_dict'][val] > 1000:
                    print 'had to exclude one result. Errorbar VERY LARGE: ',f['a'].folder
                else:
                    val_x.append(f['a'].g.attrs['LDE_decouple_time']/2.256e-6)
                    val_res.append(f['params_dict'][val]);val_res_u.append(f['error_dict'][val])
            ax.errorbar(val_x,val_res,val_res_u,label = key+' uW')
    plt.xlabel(r'decoupling time ($\tau_L$)')
    plt.ylabel(val)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()
    plt.close('all')

def plot_all_raw(fits):
    fig,ax = plt.subplots()
    for fit in fits:
        fit_x = np.linspace(fit['x'][0],fit['x'][-1],100)
        order = np.round(fit['a'].g.attrs['LDE_decouple_time']/2.256e-6,2)
        ax.errorbar(fit['x'],fit['y'],fit['y_u'],fmt='o',label=r't = %s * $\tau_L$' % order,zorder=5)
        ax.plot(fit_x,fit['fitfunc'](fit_x),color='black',zorder=0)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.xlabel(fit['a'].g.attrs['sweep_name'])

#########################
""" Script execution """
#########################

fit_res_dict = {} ### will host all fit results to play around after finishing the script
for p in analyzed_powers:
	search_string = 'C%s_%suW_X' % (carbon,p)
	kws.update({'contains':search_string})
	fits = pu_ssro.number_of_repetitions(**kws)
	fit_res_dict.update({str(p): fits})
	print str(p)+' uW done!'
plot_fit_vals(fit_res_dict,['T'])
plot_fit_vals(fit_res_dict,['n'])

