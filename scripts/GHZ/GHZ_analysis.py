import numpy as np
import os
import pickle
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from analysis.lib.m2.ssro import ssro
from analysis.lib.QEC import ConditionalParity as CP; reload(CP)
from analysis.lib.fitting import fit, common, ramsey;reload(common); reload(fit)
import matplotlib.cm as cm
import matplotlib as mpl; reload(mpl)
from pylab import *

reload (CP)
import h5py
import csv


RO_correction = True

from matplotlib import pyplot as plt
script_name = 'GHZ_analysis.py'


def get_data(folder,folder_timestamp,RO_correct=True,ssro_calib_timestamp=None,tomo_bases='XXX'):
    a = CP.ConditionalParityAnalysis(folder)
    #a.get_sweep_pts()
    a.get_readout_results(name='adwindata', post_select_GHZ = True)
    orientation_d = a.orientations[3]
    if orientation_d == 'negative':
        multiply_by = -1
    else:
        multiply_by = 1
    #print(a.orientations,multiply_by)

    if RO_correct == True:
        if ssro_calib_timestamp == None: 
            ssro_calib_folder = toolbox.latest_data('SSRO',older_than=folder_timestamp)
            analysis_file=os.path.join(ssro_calib_folder,'analysis.hdf5')
            if not os.path.exists(analysis_file):
                ssro.ssrocalib(ssro_calib_folder)
        else:
            ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
            ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'#GHZ_'+tomo_bases
        
        #print ssro_calib_folder
        a.get_electron_ROC(ssro_calib_folder = ssro_calib_folder, post_select_GHZ = True)


        #take the [0,0] part of the arrays right now; no sweep_pts
        p0_000,u_p0_000 = multiply_by*(a.p0_000[0,0]-0.5)*2,(2*a.u_p0_000[0,0]) 
        p0_001,u_p0_001 = multiply_by*(a.p0_001[0,0]-0.5)*2,(2*a.u_p0_001[0,0])
        p0_010,u_p0_010 = multiply_by*(a.p0_010[0,0]-0.5)*2,(2*a.u_p0_010[0,0])
        p0_011,u_p0_011 = multiply_by*(a.p0_011[0,0]-0.5)*2,(2*a.u_p0_011[0,0])
        p0_100,u_p0_100 = multiply_by*(a.p0_100[0,0]-0.5)*2,(2*a.u_p0_100[0,0])
        p0_101,u_p0_101 = multiply_by*(a.p0_101[0,0]-0.5)*2,(2*a.u_p0_101[0,0])
        p0_110,u_p0_110 = multiply_by*(a.p0_110[0,0]-0.5)*2,(2*a.u_p0_110[0,0])
        p0_111,u_p0_111 = multiply_by*(a.p0_111[0,0]-0.5)*2,(2*a.u_p0_111[0,0])

    else:
        p0_000,u_p0_000 = multiply_by*(a.normalized_ssro_000[0,0]-0.5)*2,(2*a.u_normalized_ssro_000[0,0])
        p0_001,u_p0_001 = multiply_by*(a.normalized_ssro_001[0,0]-0.5)*2,(2*a.u_normalized_ssro_001[0,0])
        p0_010,u_p0_010 = multiply_by*(a.normalized_ssro_010[0,0]-0.5)*2,(2*a.u_normalized_ssro_010[0,0])
        p0_011,u_p0_011 = multiply_by*(a.normalized_ssro_011[0,0]-0.5)*2,(2*a.u_normalized_ssro_011[0,0])
        p0_100,u_p0_100 = multiply_by*(a.normalized_ssro_100[0,0]-0.5)*2,(2*a.u_normalized_ssro_100[0,0])
        p0_101,u_p0_101 = multiply_by*(a.normalized_ssro_101[0,0]-0.5)*2,(2*a.u_normalized_ssro_101[0,0])
        p0_110,u_p0_110 = multiply_by*(a.normalized_ssro_110[0,0]-0.5)*2,(2*a.u_normalized_ssro_110[0,0])
        p0_111,u_p0_111 = multiply_by*(a.normalized_ssro_111[0,0]-0.5)*2,(2*a.u_normalized_ssro_111[0,0])

    p = (a.p000[0,0],a.p001[0,0],a.p010[0,0],a.p011[0,0],a.p100[0,0],a.p101[0,0],a.p110[0,0],a.p111[0,0])
    y = (p0_000,p0_001,p0_010,p0_011,p0_100,p0_101,p0_110,p0_111)
    y_avg = np.mean(y, axis=0)
    #print y_avg
    y_err = (u_p0_000,u_p0_001,u_p0_010,u_p0_011,u_p0_100,u_p0_101,u_p0_110,u_p0_111)

    return(p,y,y_err)

def get_debug_data(folder,folder_timestamp,RO_correct=True,ssro_calib_timestamp=None,tomo_bases="XXX"):
    a = CP.ConditionalParityAnalysis(folder)
    #a.get_sweep_pts()
    a.get_readout_results(name='adwindata', post_select = True, orientation_correct=True)
    orientation_b = a.orientations[1]
    if orientation_b == 'negative':
        multiply_by = -1
    else:
        multiply_by = 1
    #print(a.orientations,multiply_by)

    if RO_correct == True:
        if ssro_calib_timestamp == None: 
            ssro_calib_folder = toolbox.latest_data('SSRO',older_than=folder_timestamp)
            if not os.path.exists(os.path.join(ssro_calib_folder,'analysis.hdf5')):
                ssro.ssrocalib(ssro_calib_folder)
        else:
            ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
            ssro_calib_folder = toolbox.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_111_1_sil18'#GHZ_'+tomo_bases

        #print ssro_calib_folder
        a.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)

        #take the [0,0] part of the arrays right now; no sweep_pts
        p0_0,u_p0_0 = multiply_by*(a.p0_0[0,0]-0.5)*2,(2*a.u_p0_0[0,0]) 
        p0_1,u_p0_1 = multiply_by*(a.p0_1[0,0]-0.5)*2,(2*a.u_p0_1[0,0])

    else:
        p0_0,u_p0_0 = multiply_by*(a.normalized_ssro_0[0,0]-0.5)*2,(2*a.u_normalized_ssro_0[0,0])
        p0_1,u_p0_1 = multiply_by*(a.normalized_ssro_1[0,0]-0.5)*2,(2*a.u_normalized_ssro_1[0,0])
  
    p = (a.p0,a.p1)
    y = (p0_0,p0_1)
    y_avg = np.mean(y, axis=0)
    #print y_avg
    y_err = (u_p0_0,u_p0_1)

    return(p,y,y_err)

def do_plot(folder,timestamp,name,p,y,y_err,x_labels,x, savc=True, return_ax=False, show_plot=True,tomo_bases='XXX',title='',width=None,ylabel=''):
    plt.close('all')
    if width==None:
        fig,ax = plt.subplots() 
    else:
        fig,ax = plt.subplots(figsize=(width,4))
        
    rects = ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
    ax.set_xticks(x)
    # ax.title = timestamp
    ax.set_xticklabels(x_labels)
    ax.set_ylim(-1.1,1.1)
    ax.set_title(str(folder)+'/'+str(timestamp))
    # ax.grid()
    ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    if ylabel=='':
        ax.set_ylabel('<+'+tomo_bases+'>')
    else:
        ax.set_ylabel(ylabel)

    avg_fid = sum(y)/len(y)
    u_avg_fid = sum(y_err)/len(y_err)/sqrt(len(y_err))
 
    # print values on bar plot
    def autolabel(rects):
        for ii,rect in enumerate(rects):
            height = rect.get_height()
            plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2))+'('+str(int(round(y_err[ii]*100)))+')',
                ha='center', va='bottom')
            plt.text(rect.get_x()+rect.get_width()/2., -0.2, str(int(round(p[ii]*100)))+'%',
                ha='center', va='bottom')
            plt.text(0.1,-0.5,'average = '+str(round(avg_fid,2))+'('+str(int(round(u_avg_fid*100)))+')',ha='left', va='bottom')
            plt.text(0.5,0.8,title,ha='left', va='bottom')

    autolabel(rects)
    
    if save and ax != None:
        try:
            fig.savefig(
                os.path.join(folder,name+'.png'))
        except:
            print 'Figure has not been saved.'

    if show_plot:
        plt.show()

    if return_ax:
        return(ax)

def plot_sweep_orientations_GHZ_data(plot=True, save = True, plot_single=True, tomo_bases = 'XXX',extra_tag='',orientations=4,ssro_calib_timestamp=None,return_data=False):
    
    if orientations==4:
        data_labels=[
            'pppp',
            'pppn',
            'ppnp',
            'ppnn',
            'pnpp',
            'pnpn',
            'pnnp',
            'pnnn',
            'nppp',
            'nppn',
            'npnp',
            'npnn',
            'nnpp',
            'nnpn',
            'nnnp',
            'nnnn'
            ] 
    elif orientations==2:
        data_labels = ['pp','pn','np','nn']

    elif orientations==1:
        data_labels = ['pppp','pppn']

    #print data_labels

    p = np.zeros([len(data_labels),8])
    y = np.zeros([len(data_labels),8])
    y_err = np.zeros([len(data_labels),8])

    x_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')
    print 'looking for filename containing '+extra_tag+tomo_bases+'_'

    for k,label in enumerate(data_labels):
        timestamp,folder = toolbox.latest_data(contains=extra_tag+tomo_bases+'_'+label,return_timestamp=True)
        p[k],y[k],y_err[k] = get_data(folder,timestamp,RO_correct=True,ssro_calib_timestamp=ssro_calib_timestamp,tomo_bases=tomo_bases)
        if plot_single:
            do_plot(folder,timestamp,'GHZ_results',p[k],y[k],y_err[k],x_labels,range(len(y[k])),show_plot=False,tomo_bases=tomo_bases)

    print folder
  
    p_avg = np.mean(p, axis = 0)
    y_avg = np.mean(y, axis = 0)
    y_err_avg = np.mean(y_err, axis = 0)
    x = range(len(y_avg)) 

    do_plot(folder,timestamp,'GHZ_results_all_orientations',p_avg,y_avg,y_err_avg,x_labels,x,tomo_bases=tomo_bases,show_plot=plot)

    if return_data == True:
        return(folder,timestamp,p_avg,y_avg,y_err_avg)

def plot_singlerun_GHZ_data(timestamp=None, plot=True, save = True):
    if timestamp == None:
        timestamp, folder   = toolbox.latest_data(folder_name,return_timestamp =True)
    else: 
        folder = toolbox.data_from_time(timestamp) 

    p,y,y_err = get_data(folder,timestamp,RO_correct=True)

    x_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')
    x = range(len(y)) 

    do_plot(folder,timestamp,'GHZ_results',p,y,y_err,x_labels,x)

def plot_debug_GHZ_data(plot=True, save = True,tomo_bases="ZZI",plot_single=True):
    p = np.zeros([4,2])
    y = np.zeros([4,2])
    y_err = np.zeros([4,2])

    x_labels = ('1', '-1')
    data_labels = ['pp','pn','np','nn']

    for k,label in enumerate(data_labels):
        timestamp,folder = toolbox.latest_data(contains='tomo'+tomo_bases+'_'+label,return_timestamp=True)
        print folder
        p[k],y[k],y_err[k] = get_debug_data(folder,timestamp,RO_correct=True,tomo_bases=tomo_bases)
        if plot_single:
            do_plot(folder,timestamp,'GHZ_debug_results',p[k],y[k],y_err[k],x_labels,range(len(y[k])),show_plot=False,tomo_bases=tomo_bases)

    print(p,y,y_err)

    p_avg = np.mean(p, axis = 0)
    y_avg = np.mean(y, axis = 0)
    y_err_avg = np.mean(y_err, axis = 0)
    x = range(len(y_avg)) 

    do_plot(folder,timestamp,'GHZ_debug_results_all_orientations',p_avg,y_avg,y_err_avg,x_labels,x,tomo_bases=tomo_bases)



def plot_sweepphases_GHZ_data(timestamps=[None,None], plot=True, save = True,plot_single=True):
    points = 10
    phases = np.linspace(-360,0,points)

    p = np.zeros([10,8])
    y = np.zeros([10,8])
    y_err = np.zeros([10,8])
    x_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')

    if timestamps[0]==None:
        for kk,phase in enumerate(phases):
            timestamp,folder = toolbox.latest_data(contains='phases_'+str(phase),return_timestamp=True)
            p[kk],y[kk],y_err[kk] = get_data(folder,timestamp)
            if plot_single:
                do_plot(folder,timestamp,'GHZ_results',p[kk],y[kk],y_err[kk],x_labels,range(len(y[kk])),show_plot=False)
    
    print(p)
    p_000=p[:,0]
    print(p_000)
    y_000=y[:,0]
    y_err_000=y_err[:,0]

    fig,ax = plt.subplots() 
    ax.plot(phases,y_000,color = 'k' )#,show_plot=False
    ax.set_ylim(-1.1,1.1)
    ax.set_title(str(folder)+'/'+str(timestamp))
    #ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
    ax.set_ylabel('+- <xxx>')  

    if save and ax != None:
        try:
            fig.savefig(
                os.path.join(folder,'GHZ_vs_phase.png'))
        except:
            print 'Figure has not been saved.'

    plt.show()

def plot_full_tomo(bases= ["XXX"], plot=True):
    
    p = np.zeros([len(bases),8])
    y = np.zeros([len(bases),8])
    y_err = np.zeros([len(bases),8])
    outcome_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')
    basis_labels=bases

    for k,tomo_basis in enumerate(bases):
        folder,timestamp,p[k],y[k],y_err[k]=plot_sweep_orientations_GHZ_data(plot=False,tomo_bases=tomo_basis,show_plot=False,return_data=True)

    for l,outcomes in enumerate(outcome_labels):
        do_plot(folder,timestamp,'tomography'+outcomes,p[:,l],y[:,l],y_err[:,l],basis_labels,range(len(bases)),show_plot=plot,title=outcomes,width=len(basis_labels)*0.7,ylabel='exp. value')

def plot_nonzero_tomo(bases= ["XXX","XYY","YXY","YYX","ZZI","ZIZ","IZZ"], plot=True,plot_single=True):
    
    p = np.zeros([8,len(bases)])
    y = np.zeros([8,len(bases)])
    y_err = np.zeros([8,len(bases)])
    outcome_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')
    basis_labels=bases

    for k,tomo_basis in enumerate(bases):
        folder,timestamp,p[:,k],y[:,k],y_err[:,k]=plot_sweep_orientations_GHZ_data(plot=plot_single,tomo_bases=tomo_basis,extra_tag='_tomo_',return_data=True)

    p_avg = np.sum(p,axis=0)
    y_avg = np.sum(y*p,axis=0)
    y_err_avg = np.sum(y_err*p,axis=0)
    x = range(len(basis_labels))

    do_plot(folder,timestamp,'nonzero_tomo',p_avg,y_avg,y_err_avg,basis_labels,x,show_plot=plot,width=len(basis_labels)*1.2,ylabel='exp. value')


def plot_debug_tomo(mmts= ['IIX_IIX_IIX','IIX_IIX_XXX'], tomo = 'ZZI', init='noinit',plot=True):
    
    p = np.zeros([8,len(mmts)])
    y = np.zeros([8,len(mmts)])
    y_err = np.zeros([8,len(mmts)])
    tag_labels=mmts

    for k,mmts in enumerate(mmts):
        folder,timestamp,p[:,k],y[:,k],y_err[:,k]=plot_sweep_orientations_GHZ_data(plot=False,tomo_bases=tomo,extra_tag=init+'_mmts_'+mmts+'_tomo',orientations=1,return_data=True)
    
    p_avg = np.sum(p,axis=0)
    y_avg = np.sum(y*p,axis=0)
    y_err_avg = np.sum(y_err*p,axis=0)
    x = range(len(tag_labels)) 

    do_plot(folder,timestamp,'sweep_mmts_tomo_'+tomo,p_avg,y_avg,y_err_avg,tag_labels,x,show_plot=plot,width=len(tag_labels)*1.5,tomo_bases=tomo)



