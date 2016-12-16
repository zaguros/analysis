"""
This class is used to import data from the oscilloscope in the CSV format 
Specifically used to fit cavity resonances and extract the linewidth  
SvD  9-2016
"""


import numpy as np 
import scipy.constants
import scipy.interpolate
from matplotlib import pyplot as plt
from math import pi, sqrt, sin, cos, fabs
import csv
from itertools import *
from scipy import signal 
import itertools
from scipy.signal import argrelextrema
import pandas as pd
import os
from analysis.scripts.cavity.peakdetect import peakdet
from analysis.lib import fitting
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot

import analysis.scripts.cavity.functions_for_processing_linewidth_data as funcs

reload(common)


### This is the only part where values should be inserted.#### 
n_xticks = 5
n_yticks = 8
 

class oscilloscope_analysis():
    def __init__(self,indir,filename):
        self.indir=indir
        self.filename = filename

    def get_data(self,**kw):
        self.use_timetrace = kw.pop('use_timetrace',True)
        x_min = kw.pop('x_min',0)
        x_max = kw.pop('x_max',-1)
        nr_channels = kw.pop('nr_channels',3) #startign from 0
        data_col = kw.pop('data_col',2) #startign from 0

        # data = pd.read_csv(os.path.join(indir,filename+'.csv'), skiprows=16, names = ["X","Y"],usecols=[0,1]) #creating a dataframe in pandas and importing the data

        # data = pd.read_csv(os.path.join(indir,filename+'.csv'), skiprows=16, names = ["X","None","Y"],usecols=[0,1,2]) #creating a dataframe in pandas and importing the data
        if nr_channels == 3:
            data = pd.read_csv(os.path.join(self.indir,self.filename+'.csv'), skiprows=16, names = ["None","mod","2","3"],usecols=[0,1,2,3]) #creating a dataframe in pandas and importing the data
        elif nr_channels == 4:
            data = pd.read_csv(os.path.join(self.indir,self.filename+'.csv'), skiprows=16, names = ["None","mod","2","3","4"],usecols=[0,1,2,3,4]) #creating a dataframe in pandas and importing the data
        #print data

        #Could be used for x/y ticks and labels
        data_name = str(data_col)

        max_I = data[data_name].max()
        min_I = data[data_name].min()
        y = np.asarray(data[data_name])#[3000:]

        if self.use_timetrace:
            max_X = data['X'].max()
            min_X = data['X'].min()
            x = 1.e3*np.asarray(data['X'])#[3000:] # multiplied by factor 1.e3  to get ms  
        else:
            x = np.arange(len(y))

        mod = np.asarray(data['mod'])
        self.mod = mod[x_min:x_max]
        self.x=x[x_min:x_max]
        self.y=y[x_min:x_max]

        return self.x,self.y


    def plot_data(self,**kw):
        plot_mod = kw.pop('plot_mod',False)
        figsize = kw.pop('figsize',(6,4))
        fmt = kw.pop('fmt','o')
        fig,ax = plt.subplots(figsize=figsize)
        ax.plot(self.x,self.y,fmt)
        if plot_mod:
            scale_y =max(self.y)/ max(self.mod)
            ax.plot(self.x,self.mod*scale_y)
        ax.set_title(self.indir+'/'+self.filename)
        if self.use_timetrace:
            ax.set_xlabel("Time (ms)", fontsize = 14)
        else:
            ax.set_xlabel("datapoints", fontsize = 14)
        plt.savefig(os.path.join(self.indir,self.filename+'.png'))


    def plot_and_fit_single_peak(self,g_a1, g_A1, g_x01, g_gamma1,fixed=[],**kw):


        p0, fitfunc, fitfunc_str = common.fit_lorentz(g_a1, g_A1, g_x01, g_gamma1)
        fit_result = fit.fit1d(self.x,self.y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)
        fig,ax = plt.subplots(figsize=(8,4))
        plot.plot_fit1d(fit_result, np.linspace(self.x[0],self.x[-1],10*len(self.x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)
        ax.set_xlim(self.x[0],self.x[-1])
        if self.use_timetrace:
            ax.set_xlabel("Time (ms)", fontsize = 14)
        else:
            ax.set_xlabel("datapoints", fontsize = 14)
        ax.set_ylabel("Intensity (a.u.)", fontsize = 14)

        linewidth = fit_result['params_dict']['gamma1']
        u_linewidth = fit_result['error_dict']['gamma1']

        plt.savefig(os.path.join(self.indir,self.filename+'_fit.png'))

        plt.show()

        return linewidth, u_linewidth

    def plot_and_fit_with_EOM(self,EOM_freq,**kw):
        """
        input: EOM frequency in GHz
        """
        show_fit = kw.pop('show_fit',False)
        print_fit = kw.pop('print_fit',False)
        plot_name = kw.pop('plot_name','')
        save_fit = kw.pop('save_fit',True)
        g_a1 = kw.pop('g_a1',min(self.y))

        g_x01 = kw.pop('g_x01', self.x[np.argmax(self.y)])
        g_gamma1 = kw.pop('g_gamma1',(self.x[-1]-self.x[0])/20.)
        g_dx = kw.pop('g_dx',(self.x[-1]-self.x[0])/10.)
        g_A1 = kw.pop('g_A1',max(self.y)*g_gamma1)  
        g_A2 = kw.pop('g_A2',g_A1/1.)  
        g_A3 = kw.pop('g_A3',g_A1/1.)
        g_A4 = kw.pop('g_A4',g_A1/8.)  
        g_A5 = kw.pop('g_A5',g_A1/8.)

        fixed = kw.pop('fixed',[])

        # p0, fitfunc, fitfunc_str = common.fit_3lorentz_symmetric_asym_A(g_a1, g_A1, g_x01, g_gamma1, g_dx, g_A2,g_A3)
        # p0, fitfunc, fitfunc_str = common.fit_5lorentz_symmetric_asym_A(g_a1, g_A1, g_x01, g_gamma1, g_dx, g_A2,g_A3,g_A4,g_A5)
        
        p0, fitfunc, fitfunc_str = common.fit_3lorentz_symmetric(g_a1, g_A1, g_x01, g_gamma1, g_dx, g_A2)

        self.fit_result = fit.fit1d(self.x,self.y, None, p0=p0, fitfunc=fitfunc, do_print=print_fit, ret=True, fixed=fixed)

        if 'x01' not in self.fit_result['params_dict']:
            return 0,0,0,0,0,0,10000,g_x01,0 #return dummy values to prevent crashing

        x01 = self.fit_result['params_dict']['x01']
        dx = self.fit_result['params_dict']['dx']
        gamma1 = self.fit_result['params_dict']['gamma1']
        A1 = self.fit_result['params_dict']['A1']
        A2 = self.fit_result['params_dict']['A2']
        # gamma2 = fit_result['params_dict']['gamma2']
        u_gamma1 = self.fit_result['error_dict']['gamma1']
        u_A1 = self.fit_result['error_dict']['A1']
        chi_sq=self.fit_result['chisq']
        # u_gamma2 = fit_result['error_dict']['gamma2']

        scaling = EOM_freq/dx #scale the known EOM freq with the separation here.
        linewidth = gamma1*scaling #scale the linewidth to get linewidht in frequency
        u_linewidth = u_gamma1*scaling
        linewidth_string = 'gamma = '+str(round(linewidth,2))+'+-'+str(round(u_linewidth,3))+'GHz'
        print linewidth_string

        #Plotting

        fig,ax = plt.subplots(figsize=(8,4))
        plot.plot_fit1d(self.fit_result, np.linspace(self.x[0],self.x[-1],10*len(self.x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)
        if self.use_timetrace:
            ax.set_xlabel("Time (ms)", fontsize = 14)
        else:
            ax.set_xlabel("datapoints", fontsize = 14)        
        ax.set_ylabel("Intensity (a.u.)", fontsize = 14)
        ax.set_xlim(self.x[0],self.x[-1])
        xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],n_xticks)
        #rescaling for x-axis in GHz
        X_min_freq = ax.get_xlim()[0]*scaling
        X_max_freq = ax.get_xlim()[-1]*scaling

        ax.set_title(self.indir+'/'+self.filename+'\n'+linewidth_string)
        if save_fit:
            plt.savefig(os.path.join(self.indir,self.filename+'_fit_'+plot_name+'.png'))

        if show_fit:
            plt.show()
        plt.close()

        return linewidth, u_linewidth, A1, u_A1,gamma1,u_gamma1,chi_sq,x01,A2


    def plot_and_fit_birefringence_with_EOM(self,EOM_freq,**kw):
        """
        input: EOM frequency in GHz
        """
        show_fit = kw.pop('show_fit',False)
        g_a1 = kw.pop('g_a1',min(self.y))

        g_x01 = kw.pop('g_X01', self.x[np.argmax(self.y)])
        #g_x01 = kw.pop('g_X01', self.x[np.argmin(self.y)])
        g_gamma1 = kw.pop('g_gamma1',(self.x[-1]-self.x[0])/20.)
        #g_gamma2 = kw.pop('g_gamma2',(self.x[-1]-self.x[0])/20.)
        g_dx = kw.pop('g_dx',(self.x[-1]-self.x[0])/10.)
        g_dx2 = kw.pop('g_dx2',(self.x[-1]-self.x[0])/10.)
        g_A1 = kw.pop('g_A1',max(self.y)*g_gamma1)  
        g_A2 = kw.pop('g_A2',g_A1/8.)  
        #g_A3 = kw.pop('g_A3',g_A1*0.6)
        #g_A4 = kw.pop('g_A4',g_A1/8.)  
        #g_A5 = kw.pop('g_A5',g_A1/8.)

        fixed = kw.pop('fixed',[])

        #print g_a1, g_A1, g_x01, g_gamma1, g_dx, g_A2,g_A3

        p0, fitfunc, fitfunc_str = common.fit_4lorentz_symmetric_sym_A(g_a1, g_A1, g_x01, g_gamma1, g_dx, g_dx2, g_A2)
        # p0, fitfunc, fitfunc_str = common.fit_5lorentz_symmetric_asym_A(g_a1, g_A1, g_x01, g_gamma1, g_dx, g_A2,g_A3,g_A4,g_A5)

    #        p0, fitfunc, fitfunc_str = common.fit_3lorentz_symmetric(g_a1, g_A1, g_x01, g_gamma1, g_dx, g_A2)
        fit_result = fit.fit1d(self.x,self.y, None, p0=p0, fitfunc=fitfunc, do_print=show_fit, ret=True, fixed=fixed)

        #x01 = fit_result['params_dict']['x01']
        dx = fit_result['params_dict']['dx']
        #dx2 = fit_result['params_dict']['dx2']
        gamma1 = fit_result['params_dict']['gamma1']
        #gamma2 = fit_result['params_dict']['gamma2']
        u_gamma1 = fit_result['error_dict']['gamma1']
        #u_gamma2 = fit_result['error_dict']['gamma2']

        scaling = EOM_freq/dx #scale the known EOM freq with the separation here.
        linewidth = gamma1*scaling #scale the linewidth to get linewidht in frequency
        u_linewidth = u_gamma1*scaling
        linewidth_string = 'gamma = '+str(round(linewidth,2))+'+-'+str(round(u_linewidth,3))+'GHz'
        print linewidth_string

        #Plotting

        fig,ax = plt.subplots(figsize=(8,4))
        plot.plot_fit1d(fit_result, np.linspace(self.x[0],self.x[-1],10*len(self.x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)
        if self.use_timetrace:
            ax.set_xlabel("Time (ms)", fontsize = 14)
        else:
            ax.set_xlabel("datapoints", fontsize = 14)        
        ax.set_ylabel("Intensity (a.u.)", fontsize = 14)
        ax.set_xlim(self.x[0],self.x[-1])
        xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],n_xticks)
        #rescaling for x-axis in GHz
        X_min_freq = ax.get_xlim()[0]*scaling
        X_max_freq = ax.get_xlim()[-1]*scaling

        ax.set_title(self.indir+'/'+self.filename+'\n'+linewidth_string)
        plt.savefig(os.path.join(self.indir,self.filename+'_fit.png'))

        if show_fit:
            plt.show()
        plt.close()

        return linewidth, u_linewidth



            #xticklabels = np.linspace(X_min_freq,X_max_freq,n_xticks)

            # xticklabels_round=[]
            # for j in xticklabels:
            #   round_ = round(j,0)
            #   xticklabels_round = np.append(xticklabels_round,round_)

            # ax.set_xticks(xticks)
            # ax.set_xticklabels(xticklabels_round)




        ## you can use this for setting axes ###


        # freq_range = ((max_X1-min_X1)*EOM_freq)/h
        # print (max_X1-min_X1), EOM_freq, h

        # xticks = np.linspace(min_X1,max_X1,n_xticks)
        # xticklabels = np.linspace(min_X1,max_X1,n_xticks)

        # yticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],n_yticks)
        # ytickslabels = np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],n_yticks)

        # xticklabels_round=[]
        # for j in xticklabels:
        #   round_ = round(j,0)
        #   xticklabels_round = np.append(xticklabels_round,round_)

        # ytickslabels_round=[]
        # for i in ytickslabels:
        #   round_ = round(i,0)
        #   ytickslabels_round = np.append(ytickslabels_round,round_)

        # ax.set_yticks(yticks)
        # ax.set_yticklabels(ytickslabels_round)



############
def fit_all_resonances_in_sweep(indir,filename,EOM_freq,nr_lws,x0,windowsize,x_odd_offset,g_gamma,
                                g_dx,g_A,max_chisq,test=False,total_length=125000,data_col=3,max_plots=10,min_plots=0,return_data=False):
    '''
    function that fits the resonances in an oscilloscope trace.
    for particularly noisy/'vibrating' data, it uses a first fit to put the resonance to the middle of the windows
    TODO: improve -> instead of fitting twice, use the highest point in the first window for the middle of the second window 
                    -> have a noisy / clean data option for LT / RT measurements
    '''
    oa_single = oscilloscope_analysis(indir=indir,filename=filename)

    max_u_lw=0.1
    min_A2=0.1
    lws = np.array([])#lws_old## 
    As = np.array([])#lws_old## 
    gammas = np.array([])#lws_old## 
    chisqs= np.array([])
    x01s = np.array([])
    u_lws =np.array([])
    nr_lws=nr_lws
    x0=x0
    success =np.zeros(nr_lws)
    windowsize=windowsize


    
    for i in np.arange(nr_lws):
        if test:
            if i>max_plots-1:
                break
            if i<min_plots:
                continue
        xi = int(x0+(total_length*i/nr_lws))
        print xi



        if i%2==1:
            xi=xi+x_odd_offset
        if (xi-windowsize)>0:
            x_min = xi-windowsize
        else: 
            x_min = 0
        if (xi+windowsize)<total_length-1:
            x_max = xi+windowsize
        else:
            x_max = -1


        oa_single.get_data(use_timetrace=False,x_min = x_min, x_max=x_max,data_col=data_col)    
        oa_single.plot_data()

        lw,u_lw,A1,u_A1,gamma1,u_gamma1,chisq,x01,A2 = \
            oa_single.plot_and_fit_with_EOM(EOM_freq=EOM_freq,
               g_gamma1=g_gamma, g_dx=g_dx,g_A2=g_A,g_A3=g_A,show_fit=test,print_fit=test,save_fit=False)

        if chisq>0.4:
            print 20*'*'
            print 'lw number',int(i+1) 
            print chisq,lw
            print 'fit already not good enough! '
            print 20*'*'
            lws = np.append(lws,lw)
            As = np.append(As,A1)
            gammas = np.append(gammas,gamma1)
            chisqs = np.append(chisqs,chisq)
            x01s = np.append(x01s,x01)
            u_lws = np.append(u_lws,u_lw)

            continue

        if (x01-windowsize)>0:
            x_min = x01-windowsize    
        else: 
            x_min = 0
        if (x01+windowsize)<total_length-1:
            x_max = x01+windowsize
        else:
            x_max = None
        oa_single.get_data(use_timetrace=False,x_min = x_min, x_max=x_max,data_col=data_col)
    #     oa_single.plot_data()    
        lw,u_lw,A1,u_A1,gamma1,u_gamma1,chisq,x02,A2 = \
            oa_single.plot_and_fit_with_EOM(EOM_freq=EOM_freq,
               g_gamma1=g_gamma, g_dx=g_dx,g_A2=g_A,g_A3=g_A,show_fit=True,print_fit=False,plot_name=str(int(round(x01,0))))

        u_lws = np.append(u_lws,u_lw)
        lws = np.append(lws,lw)
        As = np.append(As,A1)
        gammas = np.append(gammas,gamma1)
        chisqs = np.append(chisqs,chisq)
        x01s = np.append(x01s,x01)
        if chisq>max_chisq:
            print 20*'*'
            print 'lw number',int(i+1) 
            print chisq,lw,'+-',u_lw
            print 'fit not good enough! '
            print 20*'*'
            continue
        if u_lw > max_u_lw:
            print 20*'*'
            print 'lw number',int(i+1) 
            print chisq,lw,'+-',u_lw
            print 'fit not good enough! uncertainty in linewidth too large '
            print 20*'*'
            continue
        if A2 < min_A2:
            print 20*'*'
            print 'lw number',int(i+1) 
            print chisq,lw,'+-',u_lw
            print A2
            print 'fit not good enough! A2 is less than 0'
            print 20*'*'
            continue      

        print 20*'*'
        print chisq,lw,'+-',u_lw,'gamma',gamma1
        print 'fit good!!! '
        print 20*'*'
        success[i] = 1

    good_x01s = x01s[np.where(success>0)]
    good_lws = lws[np.where(success>0)]
    good_chisqs = chisqs[np.where(success>0)]
    good_u_lws = u_lws[np.where(success>0)]

    good_fits = {}
    good_fits['x0s'] =good_x01s
    good_fits['lws'] = good_lws
    good_fits['chisqs'] = good_chisqs
    good_fits['u_lws']=good_u_lws
    good_fits['max_chisq']=max_chisq
    good_fits['max_u_lw']=max_u_lw
    good_fits['min_A2']=min_A2

    print '%d good fits out of %d tries'%(len(good_lws),nr_lws)
    print 'the good linewidth fits from %s :'%filename
    print 'average lw = %.2f +- %.2f'%(np.average(good_lws),np.std(good_lws))
    print 'lws = ',good_lws
    print 'u_lws = ',good_u_lws
    print 'chisq = ',good_chisqs
    print 'x0s = ',good_x01s

    funcs.save_to_json_file(indir,filename+'_analysis_2',good_fits)

    if return_data:
        return oa_single.x,oa_single.y,oa_single.fit_result

    return good_lws

