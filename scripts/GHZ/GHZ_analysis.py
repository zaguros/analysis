import numpy as np
import os
import pickle
from analysis.lib.tools import toolbox;
from analysis.lib.m2.ssro import mbi
from analysis.lib.m2.ssro import ssro
from analysis.lib.QEC import ConditionalParity as CP;reload(CP)
from analysis.lib.fitting import fit, common, ramsey
import matplotlib.cm as cm
import matplotlib as mpl; reload(mpl)
from pylab import *


import h5py
import csv
import itertools

RO_correction = True

from matplotlib import pyplot as plt
script_name = 'GHZ_analysis.py'


def get_data(folder,folder_timestamp,RO_correct=True,ssro_calib_timestamp=None,tomo_bases='XXX',get_title=False, return_orientations=False, **kw):
    evaluate_both = kw.pop('evaluate_both',False)

    a = CP.ConditionalParityAnalysis(folder)
    #a.get_sweep_pts()
    a.get_readout_results(name='adwindata', post_select_GHZ = True)
    if a.orientations[3] == 'negative':
        multiply_by = -1
    else:
        multiply_by = 1

    print a.orientations

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
        #take the [0,0] part of the arrays; there are no sweep_pts
        p0_000,u_p0_000 = multiply_by*(a.p0_000[0,0]-0.5)*2,(2*a.u_p0_000[0,0]) 
        p0_001,u_p0_001 = multiply_by*(a.p0_001[0,0]-0.5)*2,(2*a.u_p0_001[0,0])
        p0_010,u_p0_010 = multiply_by*(a.p0_010[0,0]-0.5)*2,(2*a.u_p0_010[0,0])
        p0_011,u_p0_011 = multiply_by*(a.p0_011[0,0]-0.5)*2,(2*a.u_p0_011[0,0])
        p0_100,u_p0_100 = multiply_by*(a.p0_100[0,0]-0.5)*2,(2*a.u_p0_100[0,0])
        p0_101,u_p0_101 = multiply_by*(a.p0_101[0,0]-0.5)*2,(2*a.u_p0_101[0,0])
        p0_110,u_p0_110 = multiply_by*(a.p0_110[0,0]-0.5)*2,(2*a.u_p0_110[0,0])
        p0_111,u_p0_111 = multiply_by*(a.p0_111[0,0]-0.5)*2,(2*a.u_p0_111[0,0])
        print "probablity with ROC", a.p0_000[0,0]
        print "exp  with ROC",p0_000, u_p0_000

    if evaluate_both:
    #else:
        p0_000,u_p0_000 = multiply_by*(a.normalized_ssro_000[0,0]-0.5)*2,(2*a.u_normalized_ssro_000[0,0])
        p0_001,u_p0_001 = multiply_by*(a.normalized_ssro_001[0,0]-0.5)*2,(2*a.u_normalized_ssro_001[0,0])
        p0_010,u_p0_010 = multiply_by*(a.normalized_ssro_010[0,0]-0.5)*2,(2*a.u_normalized_ssro_010[0,0])
        p0_011,u_p0_011 = multiply_by*(a.normalized_ssro_011[0,0]-0.5)*2,(2*a.u_normalized_ssro_011[0,0])
        p0_100,u_p0_100 = multiply_by*(a.normalized_ssro_100[0,0]-0.5)*2,(2*a.u_normalized_ssro_100[0,0])
        p0_101,u_p0_101 = multiply_by*(a.normalized_ssro_101[0,0]-0.5)*2,(2*a.u_normalized_ssro_101[0,0])
        p0_110,u_p0_110 = multiply_by*(a.normalized_ssro_110[0,0]-0.5)*2,(2*a.u_normalized_ssro_110[0,0])
        p0_111,u_p0_111 = multiply_by*(a.normalized_ssro_111[0,0]-0.5)*2,(2*a.u_normalized_ssro_111[0,0])
        print "prob without ROC", a.normalized_ssro_000[0,0]
        print "exp without ROC",p0_000, u_p0_000
        print "WARNING WARNING REMOVE THIS STATEMENT. REMOVE EVALUATE_BOTH OPTION"

    if get_title:
        a_list_name = "".join(aa for aa in a.a_list)
        b_list_name = "".join(aa for aa in a.b_list)
        c_list_name = "".join(aa for aa in a.c_list)
        d_list_name = "".join(aa for aa in a.d_list)
        title = a_list_name+' '+b_list_name+' '+c_list_name+' '+d_list_name
    else:
        title=''

    p = (a.p000[0,0],a.p001[0,0],a.p010[0,0],a.p011[0,0],a.p100[0,0],a.p101[0,0],a.p110[0,0],a.p111[0,0])
    y = (p0_000,p0_001,p0_010,p0_011,p0_100,p0_101,p0_110,p0_111)
    y_err = (u_p0_000,u_p0_001,u_p0_010,u_p0_011,u_p0_100,u_p0_101,u_p0_110,u_p0_111)

    if return_orientations:
        return(p,y,y_err,title,a.orientations)

    return(p,y,y_err,title)


def get_3mmt_data(folder,folder_timestamp,RO_correct=True,ssro_calib_timestamp=None,tomo_bases='XXX',get_title=False,return_orientations=False):
    a = CP.ConditionalParityAnalysis(folder)
    #a.get_sweep_pts()
    a.get_readout_results(name='adwindata', post_select_QEC = True,orientation_correct=True)
    orientations = a.orientations
    orientation_c = a.orientations[2]
    if orientation_c == 'negative':
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
        a.get_electron_ROC(ssro_calib_folder = ssro_calib_folder, post_select_QEC = True)

        #take the [0,0] part of the arrays right now; no sweep_pts
        p0_00,u_p0_00 = multiply_by*(a.p0_00[0,0]-0.5)*2,(2*a.u_p0_00[0,0]) 
        p0_01,u_p0_01 = multiply_by*(a.p0_01[0,0]-0.5)*2,(2*a.u_p0_01[0,0])
        p0_10,u_p0_10 = multiply_by*(a.p0_10[0,0]-0.5)*2,(2*a.u_p0_10[0,0])
        p0_11,u_p0_11 = multiply_by*(a.p0_11[0,0]-0.5)*2,(2*a.u_p0_11[0,0])

    else:
        p0_00,u_p0_00 = multiply_by*(a.normalized_ssro_00[0,0]-0.5)*2,(2*a.u_normalized_ssro_00[0,0])
        p0_01,u_p0_01 = multiply_by*(a.normalized_ssro_01[0,0]-0.5)*2,(2*a.u_normalized_ssro_01[0,0])
        p0_10,u_p0_10 = multiply_by*(a.normalized_ssro_10[0,0]-0.5)*2,(2*a.u_normalized_ssro_10[0,0])
        p0_11,u_p0_11 = multiply_by*(a.normalized_ssro_11[0,0]-0.5)*2,(2*a.u_normalized_ssro_11[0,0])

    if get_title:
        a_list_name = "".join(aa for aa in a.a_list)
        b_list_name = "".join(aa for aa in a.b_list)
        c_list_name = "".join(aa for aa in a.c_list)
        title = a_list_name+' '+b_list_name+' '+c_list_name
    else:
        title=''

    p = (a.p00[0,0],a.p01[0,0],a.p10[0,0],a.p11[0,0])
    y = (p0_00,p0_01,p0_10,p0_11)
    y_avg = np.mean(y, axis=0)
    #print y_avg
    y_err = (u_p0_00,u_p0_01,u_p0_10,u_p0_11)

    if return_orientations:
        return(p,y,y_err,title,orientations)

    return(p,y,y_err,title)

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

        print ssro_calib_folder
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

def do_plot(folder,timestamp,name,p,y,y_err,x_labels,x, savc=True, return_ax=False, show_plot=True,
            tomo_bases='XXX',title='',width=None,ylabel='',print_avg=False,add_title=True,hbar=False):
    plt.close('all')

    if hbar:
        if width==None:
            fig,ax = plt.subplots() 
        else:
            fig,ax = plt.subplots(figsize=(5,0.3*width))
            
        rects = ax.barh(x,y,height = 0.3,xerr=y_err,align ='center',color = 'b',ecolor = 'k' )
        ax.set_yticks(x)
        # ax.title = timestamp
        ax.set_yticklabels(x_labels)
        ax.set_xlim(-1.1,1.1)
        if add_title:
            ax.set_title(str(folder)+'/'+str(timestamp))
        else:
            ax.set_title('')
        # ax.grid()
        ax.vlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        if ylabel=='':
            ax.set_xlabel('<+'+tomo_bases+'>')
        else:
            ax.set_xlabel(ylabel)

        avg_fid = np.mean(y,axis=0)
        u_avg_fid = np.divide(np.sqrt(np.sum((y_err**2),axis=0)),len(y_err))

        #avg_fid = sum(y)/len(y)
        #u_avg_fid = sum(y_err)/len(y_err)/sqrt(len(y_err))
      
        # print values on bar plot
        def autolabel(rects):
            for ii,rect in enumerate(rects):
                # height = rect.get_height()
                # plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2))+'('+str(int(round(y_err[ii]*100)))+')',
                #     ha='center', va='bottom')
                # plt.text(rect.get_x()+rect.get_width()/2., -0.2, str(int(round(p[ii]*100)))+'%',
                #     ha='center', va='bottom')
                r_width = rect.get_width()
                plt.text( 2*r_width, rect.get_y(),str(round(y[ii],2))+'('+str(int(round(y_err[ii]*100)))+')',
                    ha='center', va='bottom')
                #plt.text(-0.2, rect.get_y()+rect.get_height()/2., str(int(round(p[ii]*100)))+'%',
                #    ha='center', va='bottom') 
    else:
        if width==None:
            fig,ax = plt.subplots() 
        else:
            fig,ax = plt.subplots(figsize=(width,4))
            
        rects = ax.bar(x,y,yerr=y_err,align ='center',ecolor = 'k' )
        ax.set_xticks(x)
        # ax.title = timestamp
        ax.set_xticklabels(x_labels)
        ax.set_ylim(-1.1,1.1)
        if add_title:
            ax.set_title(str(folder)+'/'+str(timestamp))
        else:
            ax.set_title('')
        # ax.grid()
        ax.hlines([-1,0,1],x[0]-1,x[-1]+1,linestyles='dotted')
        if ylabel=='':
            ax.set_ylabel('<+'+tomo_bases+'>')
        else:
            ax.set_ylabel(ylabel)

        avg_fid = np.mean(y,axis=0)
        u_avg_fid = np.divide(np.sqrt(np.sum((y_err**2),axis=0)),len(y_err))

        #avg_fid = sum(y)/len(y)
        #u_avg_fid = sum(y_err)/len(y_err)/sqrt(len(y_err))
     
        # print values on bar plot
        def autolabel(rects):
            for ii,rect in enumerate(rects):
                height = rect.get_height()
                plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, str(round(y[ii],2))+'('+str(int(round(y_err[ii]*100)))+')',
                    ha='center', va='bottom')
                plt.text(rect.get_x()+rect.get_width()/2., -0.2, str(int(round(p[ii]*100)))+'%',
                    ha='center', va='bottom')

    autolabel(rects)

    if print_avg:
        plt.text(0.1,-0.5,'average = '+str(round(avg_fid,2))+'('+str(int(round(u_avg_fid*100)))+')',ha='left', va='bottom')
    plt.text(-0.5,0.75,title,ha='left', va='bottom')
    
    if save and ax != None:
        try:
            fig.savefig(
                os.path.join(folder,name+'.png'))
        except:
            print 'Figure has not been saved.'
        # try:
        #     fig.savefig(
        #         os.path.join(folder,name+'.eps'))
        # except:
        #     print 'eps Figure has not been saved.'


    if show_plot:
        plt.show()

    if return_ax:
        return(ax)

def plot_sweep_orientations_GHZ_data(plot=True, save = True, plot_single=True, tomo_bases = 'XXX',
                                    extra_tag='',ssro_calib_timestamp=None,return_data=False,older_than=None,get_title=False, postselect_000 = False,**kw):
    
    nr_orientations = kw.pop('orientations',4)
    RO_correct = kw.pop('RO_correct',True)

    if nr_orientations==4:
        data_labels=['pppp','pppn','ppnp','ppnn','pnpp','pnpn','pnnp','pnnn','nppp','nppn','npnp','npnn','nnpp','nnpn','nnnp','nnnn'] 

    elif nr_orientations==3:
        data_labels = ['ppp','ppn','pnp','pnn','npp','npn','nnp','nnn',] 

    elif nr_orientations==2:
        data_labels = ['pp','pn','np','nn']

    elif nr_orientations==1:
        data_labels = ['pppp','nnnn']

    if nr_orientations != 3:
        p = np.zeros([len(data_labels),8])
        y = np.zeros([len(data_labels),8])
        y_err = np.zeros([len(data_labels) ,8])
        ROms0_p = np.zeros([len(data_labels)])
        ROms0_y = np.zeros([len(data_labels)])
        ROms0_y_err = np.zeros([len(data_labels)])

        x_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')

    elif nr_orientations ==3:
        p = np.zeros([len(data_labels),4])
        y = np.zeros([len(data_labels),4])
        y_err = np.zeros([len(data_labels) ,4])       
        x_labels = ('1,1','1,-1','-1,1','-1,-1')

    print 'looking for filename containing '+extra_tag+tomo_bases+'_'

    for k,label in enumerate(data_labels):
        timestamp,folder = toolbox.latest_data(contains=extra_tag+tomo_bases+'_'+label,older_than=older_than,return_timestamp=True)

        if nr_orientations!=3:
            p[k],y[k],y_err[k],title = get_data(folder,timestamp,RO_correct=RO_correct,ssro_calib_timestamp=ssro_calib_timestamp,tomo_bases=tomo_bases,get_title=get_title,**kw)
        elif nr_orientations ==3:
            p[k],y[k],y_err[k],title,orientations_name = get_3mmt_data(folder,timestamp,RO_correct=True,
                ssro_calib_timestamp=ssro_calib_timestamp,tomo_bases=tomo_bases,get_title=get_title, return_orientations=True)
        # print int(k/2)
        #print p[k,int(k/2)]
        #print y[k,int(k/2)]

        ROms0_p[k] = p[k,int(k/2)]
        ROms0_y[k] = y[k,int(k/2)]
        ROms0_y_err[k] = y_err[k,int(k/2)]        


        if plot_single:
            do_plot(folder,timestamp,'GHZ_results_',p[k],y[k],y_err[k],x_labels,range(len(y[k])),show_plot=False,tomo_bases=tomo_bases,title=title)

        
        try:
            oldest_timestamp
        except NameError:
            oldest_timestamp=timestamp
        oldest_timestamp=min(oldest_timestamp,timestamp)


    p_avg = np.mean(p, axis = 0)
    # y_avg = np.mean(y, axis = 0)
    # y_err_avg = np.mean(y_err, axis = 0)
    #y_avg = np.divide(np.sum(p*y,axis=0),np.sum(p,axis=0))
    #y_err_avg = np.divide(np.sqrt(np.sum((y_err*p)**2, axis=0)),np.sum(p,axis=0))
    y_avg = np.mean(y,axis=0)
    y_err_avg = np.divide(np.sqrt(np.sum((y_err**2),axis=0)),len(data_labels))
    x = range(len(y_avg)) 


    if postselect_000:
        for i in np.arange(8):
            
            p_avg[i] = (p[i*2,i]+p[i*2+1,i])/2
            y_avg[i] = (y[i*2,i]+y[i*2+1,i])/2
            y_err_avg[i] = np.sqrt ( y_err[i*2,i]**2 + y_err[i*2+1,i]**2 ) / 2
        plotname = 'GHZ_results_postselected'
    else:
        plotname = 'GHZ_results_all_orientations_'

    print y_avg
    print y_err_avg

    do_plot(folder,timestamp,plotname,p_avg,y_avg,y_err_avg,x_labels,x,tomo_bases=tomo_bases,show_plot=plot,title=title)


    if return_data == True:
        return(folder,oldest_timestamp,p_avg,y_avg,y_err_avg,title)

def plot_singlerun_GHZ_data(timestamp=None, plot=True, save = True):
    if timestamp == None:
        timestamp, folder   = toolbox.latest_data(folder_name,return_timestamp =True)
    else: 
        folder = toolbox.data_from_time(timestamp) 

    p,y,y_err = get_data(folder,timestamp,RO_correct=True)

    x_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')
    x = range(len(y)) 

    do_plot(folder,timestamp,'GHZ_results',p,y,y_err,x_labels,x)

def plot_debug_GHZ_data(plot=True, save = True,tomo_bases="ZZI",plot_single=True,extra_tag='',return_data=False):
    p = np.zeros([4,2])
    y = np.zeros([4,2])
    y_err = np.zeros([4,2])

    x_labels = ('1', '-1')
    data_labels = ['pp','pn','np','nn']

    for k,label in enumerate(data_labels):
        timestamp,folder = toolbox.latest_data(contains=extra_tag+'tomo'+tomo_bases+'_'+label,return_timestamp=True)
        print folder
        p[k],y[k],y_err[k] = get_debug_data(folder,timestamp,RO_correct=True,tomo_bases=tomo_bases)
        if plot_single:
            do_plot(folder,timestamp,'GHZ_debug_results',p[k],y[k],y_err[k],x_labels,range(len(y[k])),show_plot=False,tomo_bases=tomo_bases)
        try:
            oldest_timestamp
        except NameError:
            oldest_timestamp=timestamp
        oldest_timestamp=min(oldest_timestamp,timestamp)

    p_avg = np.mean(p, axis = 0)
    y_avg = np.mean(y, axis = 0)
    #y_err_avg = np.mean(y_err, axis = 0)
    y_err_avg = np.divide(np.sqrt(np.sum((y_err**2),axis=0)),len(data_labels))

    x = range(len(y_avg)) 

    do_plot(folder,timestamp,'GHZ_debug_results_all_orientations',p_avg,y_avg,y_err_avg,x_labels,x,tomo_bases=tomo_bases)

    if return_data == True:
        return(folder,oldest_timestamp,p_avg,y_avg,y_err_avg,title)

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

def plot_full_tomo(bases= ["XXX"], plot=True,extra_tag='',older_than=None):
    
    p = np.zeros([len(bases),8])
    y = np.zeros([len(bases),8])
    y_err = np.zeros([len(bases),8])
    outcome_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')
    basis_labels=bases

    for k,tomo_basis in enumerate(bases):
        folder,timestamp,p[k],y[k],y_err[k],title=plot_sweep_orientations_GHZ_data(plot=True,tomo_bases=tomo_basis,extra_tag=extra_tag,plot_single=True,return_data=True,older_than=older_than)

    for l,outcomes in enumerate(outcome_labels):
        do_plot(folder,timestamp,'tomography'+outcomes,p[:,l],y[:,l],y_err[:,l],basis_labels,range(len(bases)),show_plot=plot,title=outcomes,width=len(basis_labels)*0.7,ylabel='exp. value')

def plot_nonzero_tomo(bases= ["XXX","XYY","YXY","YYX","ZZI","ZIZ","IZZ","ZXI","XZI"], plot=True,plot_single=True,postselect=False,postselect_000=True):
    
    p = np.zeros([8,len(bases)])
    y = np.zeros([8,len(bases)])
    y_err = np.zeros([8,len(bases)])
    outcome_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')
    basis_labels=bases

    for k,tomo_basis in enumerate(bases):
        folder,timestamp,p[:,k],y[:,k],y_err[:,k],title=plot_sweep_orientations_GHZ_data(plot=plot_single,tomo_bases=tomo_basis,extra_tag='unbranched_tomo_',
            return_data=True,postselect_000=postselect_000)
        if postselect:
            if tomo_basis == 'XXX':
                ps = [0,3,5,6]                
            elif tomo_basis == 'XYY':
                ps = [0,1,2,3]
            elif tomo_basis == 'YXY':
                ps = [0,1,4,5]
            elif tomo_basis == 'YYX':
                ps = [0,2,4,6]
            elif tomo_basis == 'ZZI':
                ps = [2,3,4,5]
            elif tomo_basis == 'ZIZ':
                ps = [1,3,4,6]
            elif tomo_basis == 'IZZ':
                ps = [1,2,5,6]
            elif tomo_basis == 'ZXI':
                ps = [2,3,4,5]
            elif tomo_basis == 'XZI':
                ps = [2,3,4,5]
            for i in ps:
                y[i,k]=-1*y[i,k]
          

    p_avg = np.sum(p,axis=0)

    y_avg = np.mean(y,axis=0)
    y_err_avg = np.divide(np.sqrt(np.sum((y_err**2),axis=0)),len(basis_labels))    

    # y_avg = np.divide(np.sum(y*p,axis=0),np.sum(p,axis=0))
    # y_err_avg = np.divide(np.sqrt(np.sum((y_err*p)**2, axis=0)),np.sum(p,axis=0))
    x = range(len(basis_labels))

    if postselect_000:
        plotname = 'nonzero_tomo_postselected'
    else:
        plotname = 'nonzero_tomo'

    fidelity_GHZ = 1./8.+np.sum(np.abs(y_avg))/8.
    print 'fidelity of the GHZ state is: {}'.format(fidelity_GHZ) 

    do_plot(folder,timestamp,plotname,p_avg,y_avg,y_err_avg,basis_labels,x,show_plot=plot,width=len(basis_labels)*1.2,ylabel='exp. value')

def plot_permuted_old(bases= ["XXX","XYY","YXY","YYX"], plot=True,plot_single=False,postselect=False,number_permuted=6):
    
    p = np.zeros([8,len(bases)*number_permuted])
    y = np.zeros([8,len(bases)*number_permuted])
    y_err = np.zeros([8,len(bases)*number_permuted])
    outcome_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')
    basis_labels=[]
    jj=0

    for k,tomo_basis in enumerate(bases):
        older_than='20150730_000000'
        for ii in np.arange(number_permuted):
            folder,timestamp,p[:,jj],y[:,jj],y_err[:,jj],title=plot_sweep_orientations_GHZ_data(plot=plot_single,tomo_bases=tomo_basis,extra_tag='unbranched_tomo_',return_data=True,get_title=True,older_than=older_than)
            older_than=timestamp

            print title
            #title2 = title[0:3]+'\n'+title[4:7]+'\n'+title[8:11]+'\n'+title[12:15]
            basis_labels = np.append(basis_labels,title)
            if postselect:
                ps = [1,2,4,7]                 
                for ll in ps:
                    y[ll,jj]=-1*y[ll,jj]
            jj=jj+1

    p_avg = np.sum(p,axis=0)
    y_avg = np.mean(y,axis=0)
    y_err_avg = np.divide(np.sum((y_err),axis=0),len(basis_labels))    

    # y_avg = np.divide(np.sum(y*p,axis=0),np.sum(p,axis=0))
    # y_err_avg = np.divide(np.sqrt(np.sum((y_err*p)**2, axis=0)),np.sum(p,axis=0))
    average_y_err = np.divide(np.sqrt(np.sum((y_err_avg)**2, axis=0)),jj)

    x = range(jj)

    print 'average = '+str(round(mean(y_avg),3))+'('+str(int(round(average_y_err*1000,0)))+')'

    do_plot(folder,timestamp,'permutations',p_avg,y_avg,y_err_avg,basis_labels,x,show_plot=plot,width=len(basis_labels)*0.8,ylabel='exp. value',add_title=False,hbar=True)

def plot_permuted(bases= ["XXX","XII","IXI","IIX"], plot=True,plot_single=True,postselect=False,number_permuted=6):
    
    permuted_bases=list(itertools.permutations(bases))
    p = np.zeros([8,len(permuted_bases)])
    y = np.zeros([8,len(permuted_bases)])
    y_err = np.zeros([8,len(permuted_bases)])
    outcome_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')
    basis_labels=[]
    jj=0

    for k,bases in enumerate(permuted_bases):
        print bases
        tomo_basis=bases[3]
        folder,timestamp,p[:,jj],y[:,jj],y_err[:,jj],title=plot_sweep_orientations_GHZ_data(plot=plot_single,tomo_bases=tomo_basis,extra_tag='unbranched_tomo_'+bases[0]+'_'+bases[1]+'_'+bases[2]+'_',return_data=True,get_title=True)
        print title
        #title2 = title[0:3]+'\n'+title[4:7]+'\n'+title[8:11]+'\n'+title[12:15]
        basis_labels = np.append(basis_labels,title)
        if postselect:
            ps = [1,2,4,7]                 
            for ll in ps:
                y[ll,jj]=-1*y[ll,jj]
        jj=jj+1

    p_avg = np.sum(p,axis=0)
    y_avg = np.mean(y,axis=0)
    y_err_avg = np.divide(np.sqrt(np.sum((y_err**2),axis=0)),len(data_labels))    

    # y_avg = np.divide(np.sum(y*p,axis=0),np.sum(p,axis=0))
    # y_err_avg = np.divide(np.sqrt(np.sum((y_err*p)**2, axis=0)),np.sum(p,axis=0))
    average_y_err = np.divide(np.sqrt(np.sum((y_err_avg)**2, axis=0)),len(permuted_bases))
    x = range(jj)

    print 'average = '+str(round(mean(y_avg),3))+'('+str(int(round(average_y_err*1000,0)))+')'

    do_plot(folder,timestamp,'permutations',p_avg,y_avg,y_err_avg,basis_labels,x,show_plot=plot,width=len(basis_labels)*0.8,ylabel='exp. value',add_title=False,hbar=True)

def plot_permuted_carbons(carbons= ["1","2","5"], plot=True,plot_single=True,postselect=False,number_permuted=6):
    
    permuted_carbons=list(itertools.permutations(carbons))
    p = np.zeros([8,2*len(permuted_carbons)])
    y = np.zeros([8,2*len(permuted_carbons)])
    y_err = np.zeros([8,2*len(permuted_carbons)])
    outcome_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')
    basis_labels=[]
    jj=0

    for i, ro_carbons in enumerate([['1','5','2'],['5','1','2']]):
        for k,inv_carbons in enumerate(permuted_carbons):
            print ro_carbons,inv_carbons

            tomo_basis=inv_carbons[2]
            extra_tag='unbranched_carbon_list_ro_'+ro_carbons[0]+ro_carbons[1]+ro_carbons[2]+'_inv_'+inv_carbons[0]+inv_carbons[1]
            print extra_tag
            folder,timestamp,p[:,jj],y[:,jj],y_err[:,jj],title=plot_sweep_orientations_GHZ_data(plot=plot_single,tomo_bases=tomo_basis,extra_tag=extra_tag,return_data=True,get_title=True)
            print title
            title = 'ro '+ro_carbons[0]+ro_carbons[1]+ro_carbons[2]+', inv '+inv_carbons[0]+inv_carbons[1]+inv_carbons[2]
            #title2 = title[0:3]+'\n'+title[4:7]+'\n'+title[8:11]+'\n'+title[12:15]
            basis_labels = np.append(basis_labels,title)
            if postselect:
                ps = [1,2,4,7]                 
                for ll in ps:
                    y[ll,jj]=-1*y[ll,jj]
            jj=jj+1

    p_avg = np.sum(p,axis=0)
    y_avg = np.mean(y,axis=0)
    y_err_avg = np.divide(np.sqrt(np.sum((y_err**2),axis=0)),len(data_labels))    

    # y_avg = np.divide(np.sum(y*p,axis=0),np.sum(p,axis=0))
    # y_err_avg = np.divide(np.sqrt(np.sum((y_err*p)**2, axis=0)),np.sum(p,axis=0))
    average_y_err = np.divide(np.sqrt(np.sum((y_err_avg)**2, axis=0)),len(permuted_carbons))
    x = range(jj)

    print 'average = '+str(round(mean(y_avg),3))+'('+str(int(round(average_y_err*1000,0)))+')'

    do_plot(folder,timestamp,'permutations',p_avg,y_avg,y_err_avg,basis_labels,x,show_plot=plot,width=len(basis_labels)*0.8,ylabel='exp. value',add_title=False,hbar=True)



def plot_contextuality(bases= ["XYY_YXY_YYX_XXX","XII_IYI_IIY_XYY","YII_IXI_IIY_YXY","YII_IYI_IIX_YYX","XII_IXI_IIX_XXX"], plot=True,plot_single=True,postselect=True,**kw):
    
    p = np.zeros([8,len(bases)])
    y = np.zeros([8,len(bases)])
    y_err = np.zeros([8,len(bases)])
    outcome_labels = ('1,1,1', '1,1,-1' , '1,-1,1', '1,-1,-1','-1,1,1', '-1,1,-1' , '-1,-1,1', '-1,-1,-1')
    basis_labels=[]
    jj=0

    for k,tomo_basis in enumerate(bases):
        folder,timestamp,p[:,k],y[:,k],y_err[:,k],title=plot_sweep_orientations_GHZ_data(plot=plot_single,tomo_bases=tomo_basis,extra_tag='_branched_contextuality__',return_data=True,get_title=True,**kw)
        print title
        title2 = title[0:3]+'\n'+title[4:7]+'\n'+title[8:11]+'\n'+title[12:15]
        basis_labels = np.append(basis_labels,title2)
        if postselect:
            ps = [1,2,4,7]                 
            for ll in ps:
                y[ll,k]=-1*y[ll,k]

    p_avg = np.sum(p,axis=0)
    #y_avg = np.sum(y*p,axis=0)
    #y_err_avg = np.sum(y_err*p,axis=0)
    
    y_avg = np.mean(y,axis=0)
    y_err_avg = np.divide(np.sqrt(np.sum((y_err**2),axis=0)),len(data_labels))    

    # y_avg = np.divide(np.sum(y*p,axis=0),np.sum(p,axis=0))
    # y_err_avg = np.divide(np.sqrt(np.sum((y_err*p)**2, axis=0)),np.sum(p,axis=0))
    x = range(len(basis_labels))

    expectation=str(round(-y_avg[0]+y_avg[1]+y_avg[2]+y_avg[3]+y_avg[4],2))
    expectation_err=str(int(round((y_err_avg[0]+y_err_avg[1]+y_err_avg[2]+y_err_avg[3]+y_err_avg[4])/5*100,0)))
    NQI_result='<x1 y2 y2 xyy> + <y1 x2 y3 yxy> + <y1 y2 x3 yyx> + \n<x1 x2 x3 xxx> - <xyy yxy yyx xxx> = '+expectation+'('+expectation_err+')'

    do_plot(folder,timestamp,'nonzero_tomo',p_avg,y_avg,y_err_avg,basis_labels,x,show_plot=plot,width=len(basis_labels)*1.3,ylabel='exp. value',title=NQI_result,
            print_avg=False,add_title=False)


def plot_3mmt_contextuality(bases= ["YI_IX_YX","XI_IY_XY","YY_XX_ZZ","YI_IY_YY","XI_IX_XX","YX_XY_ZZ"], plot=True,plot_single=True,postselect=False):
    
    p = np.zeros([4,len(bases)])
    y = np.zeros([4,len(bases)])
    y_err = np.zeros([4,len(bases)])
    outcome_labels = ('1,1', '1,-1' , '-1,1', '-1,-1')
    basis_labels=[]
    jj=0

    for k,tomo_basis in enumerate(bases):
        folder,timestamp,p[:,k],y[:,k],y_err[:,k],title=plot_sweep_orientations_GHZ_data(plot=plot_single,tomo_bases=tomo_basis,extra_tag='_unbranched_composite_pi__',return_data=True,get_title=True,orientations=3)
        print title
        title2 = title[0:2]+'\n'+title[3:5]+'\n'+title[6:8]
        basis_labels = np.append(basis_labels,title2)
        if postselect:
            ps = [1,2]                 
            for ll in ps:
                y[ll,k]=-1*y[ll,k]

    p_avg = np.sum(p,axis=0)
    y_avg = np.mean(y,axis=0)
    y_err_avg = np.divide(np.sqrt(np.sum((y_err**2),axis=0)),len(data_labels))      
    #y_avg = np.sum(y*p,axis=0)
    # y_err_avg = np.sum(y_err*p,axis=0)
    # y_avg = np.divide(np.sum(y*p,axis=0),np.sum(p,axis=0))
    # y_err_avg = np.divide(np.sqrt(np.sum((y_err*p)**2, axis=0)),np.sum(p,axis=0))
    x = range(len(basis_labels))

    expectation=str(round(y_avg[0]+y_avg[1]-y_avg[2]+y_avg[3]+y_avg[4]+y_avg[5],2))
    expectation_err=str(int(round((y_err_avg[0]+y_err_avg[1]+y_err_avg[2]+y_err_avg[3]+y_err_avg[4]+y_err_avg[5])/6*100,0)))
    NQI_result='<y1 x2 yx> + <x1 y2 xy> - <yy xx zz> + \n<y1 y2 yy> +<x1 x2 xx> + <x1 x2 xx>  = '+expectation+'('+expectation_err+')'

    do_plot(folder,timestamp,'nonzero_tomo',p_avg,y_avg,y_err_avg,basis_labels,x,show_plot=plot,width=len(basis_labels)*1.3,ylabel='exp. value',title=NQI_result,
            print_avg=False,add_title=False)



def plot_debug_tomo(mmts= ['IIX_IIX_IIX','IIX_IIX_XXX'], tomo = 'ZZI', init='noinit',plot=True):
    
    p = np.zeros([8,len(mmts)])
    y = np.zeros([8,len(mmts)])
    y_err = np.zeros([8,len(mmts)])
    tag_labels=mmts

    for k,mmts in enumerate(mmts):
        folder,timestamp,p[:,k],y[:,k],y_err[:,k],title=plot_sweep_orientations_GHZ_data(plot=False,tomo_bases=tomo,extra_tag=init+'_mmts_'+mmts+'_tomo',orientations=2,return_data=True)
    
    p_avg = np.sum(p,axis=0)
    y_avg = np.mean(y,axis=0)
    y_err_avg = np.divide(np.sqrt(np.sum((y_err**2),axis=0)),len(data_labels))    
    # y_avg = np.sum(y*p,axis=0)
    # y_err_avg = np.sum(y_err*p,axis=0)
    x = range(len(tag_labels)) 

    do_plot(folder,timestamp,'sweep_mmts_tomo_'+tomo,p_avg,y_avg,y_err_avg,tag_labels,x,show_plot=plot,width=len(tag_labels)*1.5,tomo_bases=tomo)




def plot_debug_2mmt_tomo(bases= ['X','Y'], extra_tag='',plot=True):
    
    p = np.zeros([2,len(bases)])
    y = np.zeros([2,len(bases)])
    y_err = np.zeros([2,len(bases)])
    tag_labels=bases

    for k,basis in enumerate(bases):
        folder,timestamp,p[:,k],y[:,k],y_err[:,k],title=plot_debug_GHZ_data(plot=False,tomo_bases=basis,extra_tag=extra_tag+'_',return_data=True)
    
    p_avg = np.sum(p,axis=0)
    y_avg = np.mean(y,axis=0)
    y_err_avg = np.divide(np.sqrt(np.sum((y_err**2),axis=0)),len(tag_labels))    
    # y_avg = np.sum(y*p,axis=0)
    # y_err_avg = np.sum(y_err*p,axis=0)
    x = range(len(tag_labels)) 

    do_plot(folder,timestamp,'2mmts_tomo_',p_avg,y_avg,y_err_avg,tag_labels,x,show_plot=plot,width=len(tag_labels)*1.5,tomo_bases='ev')

