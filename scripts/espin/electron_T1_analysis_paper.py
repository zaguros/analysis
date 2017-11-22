import numpy as np
import os, sys
from analysis.lib.m2.ssro import sequence
from analysis.lib.tools import plot
from analysis.lib.tools import toolbox

from matplotlib import pyplot as plt
from analysis.lib.fitting import fit, common
from matplotlib import pyplot as plt

import matplotlib as matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 11}

matplotlib.rc('font', **font)


def get_T1_data_uncorrected(folder):
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    #a.get_electron_ROC()
    #a.get_cr_results('ssro')
    #a.plot_cr_vs_sweep()
    x = a.sweep_pts
    #y = a.p0
    y=a.ssro_results
    reps=a.reps
    print reps
    print y
    #y_err = a.u_p0
    #y_err=a.u_normalized_ssro
    #ax = a.plot_result_vs_sweepparam(ret='ax')
    return x,y,reps


def electron_T1_mul_3_uncorrected(older_than='20161111_091500',newer_than='20161110_224400', Amplitude=1, offset=1, T1=1e9, do_print = False,contains='T1',return_x='3'):

    Folder_list = toolbox.latest_data(contains=contains,older_than=older_than,newer_than=newer_than,return_all=True)
    if return_x=='3' or return_x=='all':    
        x_tot=np.zeros(3)
        y_tot=np.zeros(3)
        reps_tot=np.zeros(3)
        y_var_tot=np.zeros(3)
    # elif return_x == 'all':    
    #     x_tot=np.zeros(3)
    #     y_tot=np.zeros(3)
    #     reps_tot=np.zeros(3)
    #     y_var_tot=np.zeros(3)
    else:
        x_tot=np.zeros(2)
        y_tot=np.zeros(2)
        reps_tot=np.zeros(2)
        y_var_tot=np.zeros(2)


    for i in range(len(Folder_list)):
        print Folder_list[len(Folder_list)-i-1]
        Folder = Folder_list[len(Folder_list)-i-1]
        x,y,reps = get_T1_data_uncorrected(Folder)
        y_tot+=y
        reps_tot+=reps

    x_tot=x

    a = sequence.SequenceAnalysis(Folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')

    print reps_tot
    a.ssro_results=y_tot
    a.reps=reps_tot
    print a.reps


    a.normalized_ssro = a.ssro_results/((a.reps)/len(a.sweep_pts))
    a.u_normalized_ssro = \
    (a.normalized_ssro*(1.-a.normalized_ssro)/((a.reps)/len(a.sweep_pts)))**0.5  #this is quite ugly, maybe replace?

    a.get_electron_ROC()
    y_tot=a.p0
    y_var_tot=a.u_p0

    if return_x =='3':
        return x_tot[2],y_tot[2],y_var_tot[2]
    elif return_x == '5min':
        return x_tot[1],y_tot[0],y_var_tot[0]
    elif return_x == '2':
        return x_tot[1],y_tot[1],y_var_tot[1]
    elif return_x == 'all' or return_x == '2_all' :
        return x_tot,y_tot,y_var_tot


def get_T1_data(folder):
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()
    #a.get_cr_results('ssro')
    #a.plot_cr_vs_sweep()
    x = a.sweep_pts
    y = a.p0
    #y=a.normalized_ssro
    #print y
    y_err = a.u_p0
    #y_err=a.u_normalized_ssro
    #ax = a.plot_result_vs_sweepparam(ret='ax')
    return x,y,y_err






def electron_T1_mul(older_than='20161110_180000',newer_than='20161110_141200',mode='init0', Amplitude=1, offset=1, T1=1e9, do_print = True):

    Folder_list = toolbox.latest_data(contains='T1',older_than=older_than,newer_than=newer_than,return_all=True)
    Directory = Folder_list[0]
    Date = Folder_list[1]
    x_tot=[]
    y_tot=[]
    y_var_tot=[]

    for i in range(len(Folder_list)):
        #print Folder_list[len(Folder_list)-i-1]
        Folder = Folder_list[len(Folder_list)-i-1]
        x,y,y_var = get_T1_data(Folder)
        #print y
        x_tot.extend(list(x))
        y_tot.extend(list(y))
        y_var_tot.extend(list(y_var))
        


    ################################################data from 20s to 60s #####################################################
    x1,y1,y_var1 = electron_T1_mul_1(older_than='20161110_211700',newer_than='20161110_172900')


    x_tot.extend(list(x1))
    y_tot.extend(list(y1))
    y_var_tot.extend(list(y_var1))

    ############################################################## data at 5 mins
    x2,y2,y_var2 = electron_T1_mul_3_uncorrected(older_than='20161112_033700',newer_than='20161110_224600',return_x='2')
   
    x_tot.append(x2)
    y_tot.append(y2)
    y_var_tot.append(y_var2)
    ############################################################## data at 10 mins

    x3,y3,y_var3 = electron_T1_mul_3_uncorrected(older_than='20161113_200000',newer_than='20161112_033700',return_x='2')
   
    x_tot.append(x3)
    y_tot.append(y3)
    y_var_tot.append(y_var3)
    #################################################################

    x_tot=np.array(x_tot)/1e3
    y_tot=np.array(y_tot)
    y_var_tot=np.array(y_var_tot)
    #ax_tot=np.array(ax_tot)


    #print x_tot
    #print y_tot
    #print y_var_tot

    p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(offset, Amplitude, T1)
    # fit_result = fit.fit1d(x_tot,y_tot, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True)
    # #plt.plot(x_tot,y_tot)
    # plot.plot_fit1d(fit_result, np.linspace(0,x_tot[-1],201), ax= None, plot_data=False)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Fidelity')
    ax.errorbar(x_tot,y_tot,fmt='o',yerr=y_var_tot)
    ax.set_ylim([0.3,1.1])
    ax.set_xlim([0,650])
    


    fit_result = fit.fit1d(x_tot,y_tot, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True)
    # plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)
        
    

    plt.savefig(os.path.join(Folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(Folder, 'analyzed_result.png'),
    format='png')


    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.grid()
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Fidelity')
    ax1.errorbar(x_tot,y_tot,fmt='o',yerr=y_var_tot)
    ax1.set_ylim([0.3,1.1])
    ax1.set_xscale('log')
    
    

    plt.savefig(os.path.join(Folder, 'analyzed_result_log.pdf'),
    format='pdf')
    plt.savefig(os.path.join(Folder, 'analyzed_result_log.png'),
    format='png')

    print 'data saved in  ' + str(Folder)

    return x_tot, y_tot, y_var_tot



def electron_T1_mul_1(older_than='20161111_091500',newer_than='20161110_224400',mode='init0', Amplitude=1, offset=1, T1=1e9, do_print = False):

    Folder_list = toolbox.latest_data(contains='T1',older_than=older_than,newer_than=newer_than,return_all=True)
    Directory = Folder_list[0]
    Date = Folder_list[1]
    x_tot=np.zeros(3)
    y_tot=np.zeros(3)
    y_var_tot=np.zeros(3)

    for i in range(len(Folder_list)):
        #print Folder_list[len(Folder_list)-i-1]
        Folder = Folder_list[len(Folder_list)-i-1]
        x,y,y_var = get_T1_data(Folder)
        y_tot+=y
        y_var_tot+=y_var**(2.0)
        #ax_tot.extend(list(ax))





    x_tot=x
    y_tot=y_tot/len(Folder_list)
    y_var_tot=y_var_tot**(1/2.0)/len(Folder_list)

    # print 'x_tot is  ' +str(x_tot)
    # print 'y_tot is  ' +str(y_tot)
    # print 'error_tot is  ' +str(y_var_tot)

    # # p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(offset, Amplitude, T1)
    # # fit_result = fit.fit1d(x_tot,y_tot, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True)
    # # #plt.plot(x_tot,y_tot)
    # # plot.plot_fit1d(fit_result, np.linspace(0,x_tot[-1],201), ax= None, plot_data=False)

    # fig, ax = plt.subplots(1, 1)
    
    # ax.errorbar(x_tot, y_tot, yerr=y_var_tot, fmt='o')
    # #ax.title(' combined measurement ')
    # ax.set_ylim([0.5,1.2])
    # ax.set_xlim([-10000,310000])


    # plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    # format='pdf')
    # plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    # format='png')


    return x_tot,y_tot,y_var_tot



def electron_T1_mul_2(older_than='20161111_091500',newer_than='20161110_224400', Amplitude=1, offset=1, T1=1e9, do_print = False):

    Folder_list = toolbox.latest_data(contains='T1',older_than=older_than,newer_than=newer_than,return_all=True)
    Directory = Folder_list[0]
    Date = Folder_list[1]
    x_tot=np.zeros(2)
    y_tot=np.zeros(2)
    y_var_tot=np.zeros(2)

    for i in range(len(Folder_list)):
        print Folder_list[len(Folder_list)-i-1]
        Folder = Folder_list[len(Folder_list)-i-1]
        x,y,y_var = get_T1_data(Folder)
        y_tot+=y
        y_var_tot+=y_var**(2.0)
        #print x
        #print y
        #print y_var
        

    print 'len(Folder_list) is' + str(len(Folder_list))
    x_tot=x
    y_tot=y_tot/len(Folder_list)
    y_var_tot=y_var_tot**(1/2.0)/len(Folder_list)

    #print 'x_tot is  ' +str(x_tot)
    #print 'y_tot is  ' +str(y_tot)
    #print 'error_tot is  ' +str(y_var_tot)




    # # p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(offset, Amplitude, T1)
    # # fit_result = fit.fit1d(x_tot,y_tot, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True)
    # # #plt.plot(x_tot,y_tot)
    # # plot.plot_fit1d(fit_result, np.linspace(0,x_tot[-1],201), ax= None, plot_data=False)

    # fig, ax = plt.subplots(1, 1)
    
    # ax.errorbar(x_tot, y_tot, yerr=y_var_tot, fmt='o')
    # #ax.title(' combined measurement ')
    # ax.set_ylim([0.5,1.2])
    # ax.set_xlim([-10000,310000])


    # plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    # format='pdf')
    # plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    # format='png')


    return x_tot[1],y_tot[1],y_var_tot[1]
    #return x_tot,y_tot,y_var_tot


def electron_T1_mul_3(older_than='20161111_091500',newer_than='20161110_224400', Amplitude=1, offset=1, T1=1e9, do_print = False,contains='T1'):

    Folder_list = toolbox.latest_data(contains=contains,older_than=older_than,newer_than=newer_than,return_all=True)
    x_tot=np.zeros(3)
    y_tot=np.zeros(3)
    y_var_tot=np.zeros(3)

    for i in range(len(Folder_list)):
        print Folder_list[len(Folder_list)-i-1]
        Folder = Folder_list[len(Folder_list)-i-1]
        x,y,y_var = get_T1_data(Folder)
        y_tot+=y
        y_var_tot+=y_var**(2.0)
        #print x
        #print y
        #print y_var
        

    print 'len(Folder_list) is' + str(len(Folder_list))
    x_tot=x
    y_tot=y_tot/len(Folder_list)
    y_var_tot=y_var_tot**(1/2.0)/len(Folder_list)

    #print 'x_tot is  ' +str(x_tot)
    #print 'y_tot is  ' +str(y_tot)
    #print 'error_tot is  ' +str(y_var_tot)

    return x_tot[2],y_tot[2],y_var_tot[2]
    


def electron_T1_ms_1(older_than='20161117_160000',newer_than='20161117_131000',mode='init0', Amplitude=1, offset=1, T1=1e9, do_print = True):

    Folder_list = toolbox.latest_data(contains='T1',older_than=older_than,newer_than=newer_than,return_all=True)
    Folder_list1 = toolbox.latest_data(contains='T1',older_than='20161120_172000',newer_than='20161120_165340',return_all=True)
    Folder_list.extend(Folder_list1)
    x_tot=[]
    y_tot=[]
    y_var_tot=[]

    for i in range(len(Folder_list)):
        #print Folder_list[len(Folder_list)-i-1]
        Folder = Folder_list[len(Folder_list)-i-1]
        x,y,y_var = get_T1_data(Folder)
        #print y
        x_tot.extend(list(x))
        y_tot.extend(list(y))
        y_var_tot.extend(list(y_var))
 ################################################data from 20s to 60s #####################################################
    # x1,y1,y_var1 = electron_T1_mul_600_1000(older_than='20161110_211700',newer_than='20161110_172900')


    # x_tot.extend(list(x1))
    # y_tot.extend(list(y1))
    # y_var_tot.extend(list(y_var1))

    ############################################################## 
    x2,y2,y_var2 = electron_T1_mul_3(older_than='20161120_214000',newer_than='20161120_173800',contains='20s')
   
    x_tot.append(x2)
    y_tot.append(y2)
    y_var_tot.append(y_var2)
    ############################################################## 

    x3,y3,y_var3 = electron_T1_mul_3(older_than='20161120_133000',newer_than='20161120_050000',contains='40s')
   
    x_tot.append(x3)
    y_tot.append(y3)
    y_var_tot.append(y_var3)
    #################################################################
    x4,y4,y_var4 = electron_T1_mul_3(older_than='20161120_194000',newer_than='20161120_050000',contains='60s')
   
    x_tot.append(x4)
    y_tot.append(y4)
    y_var_tot.append(y_var4)

    ###############################################################
    x5,y5,y_var5 = electron_T1_mul_3_uncorrected(older_than='20161117_020000',newer_than='20161116_154000',contains='T1',return_x='5min')
   
    x_tot.append(x5)
    y_tot.append(y5)
    y_var_tot.append(y_var5)


    #################################################################
    x6,y6,y_var6 = electron_T1_mul_3_uncorrected(older_than='20161119_205000',newer_than='20161118_110000',contains='10_min',return_x='3')
   
    x_tot.append(x6)
    y_tot.append(y6)
    y_var_tot.append(y_var6)

    x_tot=np.array(x_tot)/1e3
    y_tot=np.array(y_tot)
    y_var_tot=np.array(y_var_tot)
    #ax_tot=np.array(ax_tot)


    #print x_tot
    #print y_tot
    #print y_var_tot

    p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(offset, Amplitude, T1)
    # fit_result = fit.fit1d(x_tot,y_tot, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True)
    # #plt.plot(x_tot,y_tot)
    # plot.plot_fit1d(fit_result, np.linspace(0,x_tot[-1],201), ax= None, plot_data=False)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Fidelity')
    ax.errorbar(x_tot,y_tot,fmt='o',yerr=y_var_tot)
    ax.set_ylim([0.3,1.1])
    ax.set_xlim([0,650])
    


    fit_result = fit.fit1d(x_tot,y_tot, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True)
    # plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)
        
    

    plt.savefig(os.path.join(Folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(Folder, 'analyzed_result.png'),
    format='png')


    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Fidelity')
    ax1.grid()
    ax1.errorbar(x_tot,y_tot,fmt='o',yerr=y_var_tot)
    ax1.set_ylim([0.3,1.1])
    ax1.set_xscale('log')
    
    

    plt.savefig(os.path.join(Folder, 'analyzed_result_log.pdf'),
    format='pdf')
    plt.savefig(os.path.join(Folder, 'analyzed_result_log.png'),
    format='png')

    print 'data saved in  ' + str(Folder)

    return x_tot, y_tot, y_var_tot






def electron_T1_msp1(older_than='20161121_190000',newer_than='20161121_171000',mode='init0', Amplitude=1, offset=1, T1=1e9, do_print = True):

    Folder_list = toolbox.latest_data(contains='T1',older_than=older_than,newer_than=newer_than,return_all=True)
    x_tot=[]
    y_tot=[]
    y_var_tot=[]

    for i in range(len(Folder_list)):
        print Folder_list[len(Folder_list)-i-1]
        Folder = Folder_list[len(Folder_list)-i-1]
        x,y,y_var = get_T1_data(Folder)
        #print y
        x_tot.extend(list(x))
        y_tot.extend(list(y))
        y_var_tot.extend(list(y_var))
 ################################################data from 20s to 60s #####################################################
    # x1,y1,y_var1 = electron_T1_mul_600_1000(older_than='20161110_211700',newer_than='20161110_172900')


    # x_tot.extend(list(x1))
    # y_tot.extend(list(y1))
    # y_var_tot.extend(list(y_var1))

    ############################################################## data at 5 mins
    x2,y2,y_var2 = electron_T1_mul_3(older_than='20161124_080000',newer_than='20161121_190000',contains='20s')
   
    x_tot.append(x2)
    y_tot.append(y2)
    y_var_tot.append(y_var2)
    ############################################################## data at 10 mins

    x3,y3,y_var3 = electron_T1_mul_3(older_than='20161124_080000',newer_than='20161121_190000',contains='40s')
   
    x_tot.append(x3)
    y_tot.append(y3)
    y_var_tot.append(y_var3)
    #################################################################
    x4,y4,y_var4 = electron_T1_mul_3(older_than='20161124_080000',newer_than='20161121_190000',contains='60s')
   
    x_tot.append(x4)
    y_tot.append(y4)
    y_var_tot.append(y_var4)

    ###############################################################
    x5,y5,y_var5 = electron_T1_mul_3_uncorrected(older_than='20161124_080000',newer_than='20161121_190000',contains='5_min')
   
    x_tot.append(x5)
    y_tot.append(y5)
    y_var_tot.append(y_var5)


    #################################################################
    x6,y6,y_var6 = electron_T1_mul_3_uncorrected(older_than='20161124_080000',newer_than='20161121_190000',contains='10_min')
   
    x_tot.append(x6)
    y_tot.append(y6)
    y_var_tot.append(y_var6)

    x_tot=np.array(x_tot)/1e3
    y_tot=np.array(y_tot)
    y_var_tot=np.array(y_var_tot)
    #ax_tot=np.array(ax_tot)


    # print x_tot
    # print y_tot
    # print y_var_tot

    p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(offset, Amplitude, T1)
    # fit_result = fit.fit1d(x_tot,y_tot, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True)
    # #plt.plot(x_tot,y_tot)
    # plot.plot_fit1d(fit_result, np.linspace(0,x_tot[-1],201), ax= None, plot_data=False)


    fig = plt.figure(1,figsize=(9,6))
    ax = fig.add_subplot(111)
    ax.grid()
    ax.set_xlabel('Total evolution time (s)')
    ax.set_ylabel('State Fidelity')
    ax.errorbar(x_tot,y_tot,fmt='o',yerr=y_var_tot)
    ax.set_ylim([0.3,1.1])
    ax.set_xlim([0,650])
    


    fit_result = fit.fit1d(x_tot,y_tot, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True)
    # plot.plot_fit1d(fit_result, np.linspace(0,x[-1],201), ax=ax, plot_data=False)
        
    

    plt.savefig(os.path.join(Folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(Folder, 'analyzed_result.png'),
    format='png')


    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Fidelity')
    ax1.grid()
    ax1.errorbar(x_tot,y_tot,fmt='o',yerr=y_var_tot)
    ax1.set_ylim([0.3,1.1])
    ax1.set_xscale('log')
    
    

    plt.savefig(os.path.join(Folder, 'analyzed_result_log.pdf'),
    format='pdf')
    plt.savefig(os.path.join(Folder, 'analyzed_result_log.png'),
    format='png')

    print 'data saved in  ' + str(Folder)

    return x_tot, y_tot, y_var_tot





def plot_summation(Amplitude=1, offset=0.333, T1=1000, do_print = True):

    x_tot1,y_tot1,y_var_tot1= electron_T1_mul(older_than='20161110_172800',newer_than='20161110_141500')
    x_tot2,y_tot2,y_var_tot2= electron_T1_ms_1()
    x_tot3,y_tot3,y_var_tot3= electron_T1_msp1()

    # print len(x_tot1)
    # print len(x_tot2)
    # print len(x_tot3)

    x_tot= (x_tot1+x_tot2+x_tot3)/3.0
    y_tot= (y_tot1+y_tot2+y_tot3)/3.0
    y_var_tot=((y_var_tot1**2+y_var_tot2**2+y_var_tot3**2)**0.5)
    print 'y_tot is' + str(y_tot) 
    print 'y_var_tot is ' + str(y_var_tot)

    fig = plt.figure(111,figsize=(4,2))
    ax1 = fig.add_subplot(111)
    #ax1.errorbar(x_tot,y_tot,fmt='o',yerr=y_var_tot)

    ax1.errorbar(x_tot1,y_tot1,fmt='o',color='r',label='ms=0',yerr=y_var_tot1,lw=1.7,markersize=3 )
    ax1.errorbar(x_tot2,y_tot2,fmt='s',color='b',label='ms=-1',yerr=y_var_tot2,lw=1.7,markersize=3)
    ax1.errorbar(x_tot3,y_tot3,fmt='^',color='g',label='ms=+1',yerr=y_var_tot3,lw=1.7,markersize=3)

    ax1.set_xlabel('Total evolution time (s)')
    ax1.set_ylabel('State fidelity')
    #ax1.grid()
    ax1.set_ylim([0.3,1.08])
    ax1.set_yticks([0.4,0.6,0.8,1])
    ax1.set_xscale('log')
    plt.legend(numpoints=1,loc=3,fontsize=10,frameon=False)



    p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(offset, Amplitude, T1)
    fit_result = fit.fit1d(x_tot,y_tot, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=[0])

    #y=A * exp(-x/tau) + a
    plot.plot_fit1d(fit_result, np.linspace(0,50*x_tot[-1],501), ax=ax1, add_txt = False, color='b', plot_data=False)
    plt.savefig('C:\Users\TUD277931\Dropbox\TaminiauLab\Projects\Coherence in multi-qubit systems\Paper\Figures\Fig 1\T1_tim_proposal_2.pdf',
    format='pdf',bbox_inches='tight',pad_inches=0.2,transparent=True)
    
    # plt.savefig(os.path.join(Folder, 'sum_log.png'),
    # format='png')

    #print 'data saved in  ' + str(Folder)

    return x_tot1,y_tot1,y_var_tot1

def plot_total_spin_population(Amplitude=0.69315, offset=0.3333, T1=1000,do_print=True):
    x1,y1,e1 = electron_T1_mul_3_uncorrected(older_than='20161127_092000',newer_than='20161126_074000', Amplitude=1, offset=1, T1=1e7,return_x='all')
    x2,y2,e2 = electron_T1_mul_3_uncorrected(older_than='20170407_204000',newer_than='20170406_174400', Amplitude=1, offset=1, T1=1e7,return_x='all')
    x3,y3,e3 = electron_T1_mul_3_uncorrected(older_than='20161113_200000',newer_than='20161112_033700', Amplitude=1, offset=1, T1=1e7,return_x='2_all')
 

    x2=np.delete(x2,1)
    y2=np.delete(y2,1)
    x1=np.delete(x1,1)
    y1=np.delete(y1,1)
    e1=np.delete(e1,1)
    e2=np.delete(e2,1)

    e_tot=((e1**2+e2**2+e3**2)**0.5)/3

    fig = plt.figure(211,figsize=(5,3))
    ax = fig.add_subplot(111)
    #ax.errorbar(x1*1.0e-3, y1, yerr=e1, fmt='o',color='b',label='ms=0 to ms=+1')
    #ax.errorbar(x2*1.0e-3, y2, yerr=e2, fmt='^',color='g',label='ms=0 to ms=-1')
    ax.errorbar(x3*1.0e-3, y3, yerr=e3, fmt='o',color='r',label='ms=0 to ms=0')

    ax.errorbar(x1*1.0e-3, (y1+y2+y3), yerr=e_tot, fmt='o',color='k',label='Total spin population ')
    

    y_tot=(y1+y2+y3)
    x=x1*1.0e-3
    print y1
    print y3
    print y_tot
    p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(offset, Amplitude, T1)
    fit_result = fit.fit1d(x,y_tot, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True,fixed=[0,1])

    #y=A * exp(-x/tau) + a
    plot.plot_fit1d(fit_result, np.linspace(0,50*x[-1],501), ax=ax, add_txt = False, color='b', plot_data=False)

    
    ax.set_xscale('log')
    ax.set_xlabel('Total evolution time (s)')
    ax.set_ylabel('Population')
    ax.set_xlim(1.0e-3,1.0e4)
    ax.set_ylim(0.7,1.06)
    plt.legend(fontsize=11,loc=3)

    plt.savefig('C:\Users\TUD277931\Dropbox\TaminiauLab\Projects\Coherence in multi-qubit systems\Paper\Figures\Fig 1\T1_total_spin_population.pdf',
    format='pdf',bbox_inches='tight',pad_inches=0.2,transparent=True)

    #plt.show()











def electron_T1_mul_5min(older_than='20161111_091500',newer_than='20161110_224400', Amplitude=1, offset=1, T1=1e9, do_print = False):

    Folder_list = toolbox.latest_data(contains='T1',older_than=older_than,newer_than=newer_than,return_all=True)
    Directory = Folder_list[0]
    Date = Folder_list[1]
    x_tot=np.zeros(2)
    y_tot=np.zeros(2)
    y_var_tot=np.zeros(2)

    for i in range(len(Folder_list)):
        print Folder_list[len(Folder_list)-i-1]
        Folder = Folder_list[len(Folder_list)-i-1]
        x,y,y_var = get_T1_data(Folder)
        y_tot+=y
        y_var_tot+=y_var**(2.0)
        # print x
        # print y
        # print y_var
        

    print 'len(Folder_list) is' + str(len(Folder_list))
    x_tot=x
    y_tot=y_tot/len(Folder_list)
    y_var_tot=y_var_tot**(1/2.0)/len(Folder_list)

    print 'x_tot is  ' +str(x_tot)
    print 'y_tot is  ' +str(y_tot)
    print 'error_tot is  ' +str(y_var_tot)




    # # p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(offset, Amplitude, T1)
    # # fit_result = fit.fit1d(x_tot,y_tot, None, p0=p0, fitfunc=fitfunc, do_print=do_print, ret=True)
    # # #plt.plot(x_tot,y_tot)
    # # plot.plot_fit1d(fit_result, np.linspace(0,x_tot[-1],201), ax= None, plot_data=False)

    # fig, ax = plt.subplots(1, 1)
    
    # ax.errorbar(x_tot, y_tot, yerr=y_var_tot, fmt='o')
    # #ax.title(' combined measurement ')
    # ax.set_ylim([0.5,1.2])
    # ax.set_xlim([-10000,310000])


    # plt.savefig(os.path.join(folder, 'analyzed_result.pdf'),
    # format='pdf')
    # plt.savefig(os.path.join(folder, 'analyzed_result.png'),
    # format='png')


    return x_tot[1],y_tot[0],y_var_tot[0]


# def electron_T1_anal_mul(older_than='20161110_144000',newer_than='20161110_141500',mode='init0', Amplitude=1, Offset=1, T1=1e7):
#     ''' Function to analyze T1 measurements. Loads the results and fits them to a simple exponential.
#     Inputs:
#     older_than = older than timestamp
#     newer_than = newer than timestamp
#     mode = init_0 init_1 contrast
#     '''

#     Folder_list = toolbox.latest_data(contains='T1',older_than=older_than,newer_than=newer_than,return_all=True)
#     Directory = Folder_list[0]
#     Date = Folder_list[1]
#     init0=['0','1']
#     #init0 = [k for k in Folder_list[1] if 'init0' in k]
#     # init1 = [k for k in Folder_list[2] if 'init_1' in k]
#     print 'folder_list is' + str(Folder_list)
#     print 'init0 is' + str(init0)

#     def get_T1_data(folder): #Assures that the data is in the right order
#         a = sequence.SequenceAnalysis(folder)
#         a.get_sweep_pts()
#         a.get_readout_results('ssro')
#         a.get_electron_ROC()
#         x = a.sweep_pts
#         y = a.p0
#         y_var = (a.u_p0)**2
#         minloc = -np.where(x == min(x))[0][0]
        
#         x = np.roll(x,minloc)
#         y = np.roll(y,minloc)
#         y_var = np.roll(y_var,minloc)

#         return x,y,y_var

#     if mode == 'init0' or mode == 'contrast':
#         x0_tot = np.zeros((6,))
#         y0_tot = np.zeros((6,))
#         y0_var_tot = np.zeros((6,))

#         for i in range(len(init0)):
#             Folder = Directory + "\\" + Date + "\\" + init0[i]
#             x,y,y_var = get_T1_data(Folder)
#             y0_tot += y
#             y0_var_tot += y_var

#         y0_tot /= len(init0)
#         y0_var_tot /= len(init0)


#     if mode == 'init1' or mode == 'contrast':
#         x1_tot = np.zeros((6,))
#         y1_tot = np.zeros((6,))
#         y1_var_tot = np.zeros((6),)

#         for i in range(len(init1)):
#             Folder = Directory + "\\" + Date + "\\" + init1[i]
#             x,y,y_var = get_T1_data(Folder)
#             y1_tot += y
#             y1_var_tot += y_var
        
#         y1_tot /= len(init0)
#         y1_var_tot /= (len(init0)) 

#     if mode == 'init1':
#         y_tot = y1_tot
#         y_var_tot = y1_var_tot
#     elif mode == 'init0':
#         y_tot = y1_tot
#         y_var_tot = y1_var_tot
#     elif mode == 'contrast':
#         y_diff = y0_tot - y1_tot
#         var_diff = (y0_var_tot + y1_var_tot) / 2
#     else:
#         raise Exception('Mode not specified')


#     p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(Amplitude, Offset, T1)


#     a = sequence.SequenceAnalysis(Directory + "\\" + Date + "\\" + init0[0])
#     a.get_sweep_pts()
#     a.sweep_pts = x / 1e6
#     a.sweep_name = 'Times (sec)'
#     a.get_readout_results('ssro')
#     a.get_electron_ROC()

#     a.p0 = y_diff
#     a.u_p0 = var_diff**0.5
#     ax = a.plot_result_vs_sweepparam(ret='ax')

#     fit_result = fit.fit1d(x / 1e6,y_diff, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True)
#     plot.plot_fit1d(fit_result, np.linspace(0,x[-1]/1e6,201), ax=ax, plot_data=False)