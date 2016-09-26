### Functions for analysing oscilloscope data
##SvD
###

import numpy as np 
from matplotlib import pyplot as plt
import os
import json
c = 2.997e8


#functions for analysis of linewidths from several files, 
#taking into account the up down modulation

def get_up_down_lws(fitted_lws,first_lw_direction):
    even_lws=fitted_lws[::2]
    odd_lws = fitted_lws[1::2]
    if len(even_lws)>len(odd_lws):
        even_lws=even_lws[:-1]
    updown_mean_lws=np.add(even_lws,odd_lws)/2        
    if first_lw_direction=='u':
        lws_up = even_lws
        lws_down  = odd_lws
    elif first_lw_direction=='d':
        lws_up = odd_lws
        lws_down  = even_lws
    else:
        print 'specify a valid modulation direction for the first linewidth'
        return
    return updown_mean_lws,lws_up,lws_down

def add_to_analysis_dict(analysis,lws_up,lws_down,updown_mean_lws):
    analysis['lws_up']=lws_up
    analysis['lws_down']=lws_down
    analysis['updown_mean_lws']=updown_mean_lws
    analysis['avg_up']=np.average(lws_up)
    analysis['u_avg_up']=np.std(lws_up)
    analysis['avg_down']=np.average(lws_down)
    analysis['u_avg_down']=np.std(lws_down)
    analysis['avg_updown_mean']=np.average(updown_mean_lws)
    analysis['u_avg_updown_mean']=np.std(updown_mean_lws)
    return analysis

def print_results(analysis):
    print 'avg_up', analysis['avg_up'],analysis['u_avg_up']
    print 'avg_down', analysis['avg_down'],analysis['u_avg_down']
    print 'avg_updown_mean', analysis['avg_updown_mean'],analysis['u_avg_updown_mean']

#functions for saving in json/txt
def default(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError('Not serializable')

def save_to_json_file(indir,filename,analysis_dict):
    f = open(os.path.join(indir,filename+'.json'), 'w')
    json.dump(analysis_dict,f,default=default)
    f.close()

def save_to_txt_file(indir,filename,nparray):
    f = open(os.path.join(indir,filename+'.txt'), 'w')
    np.savetxt(f,nparray)
    f.close()












########functions to further process saved linewidth data.
def get_avg_values_from_json(V,Lnr,day_dir,tag=''):
    V_str = '%1.3f'%V
    string_name = 'V'+str(V)
    #string_name = '%sL%d_V%sp%s'%(tag,Lnr,V_str[0],V_str[2:])
#     print string_name
    datapath = os.path.join(day_dir,string_name)
    analysis=0
    
#     print datapath
    for filename in os.listdir(datapath):
        if filename.endswith("analysis.json"):
#             print filename 
            f = open(os.path.join(datapath,filename), 'r')
            analysis = json.load(f)
            f.close()
            
    if analysis==0:
        print 'analysis file not found!'
        return datapath, 0
        
            
    avg_values = np.array([[analysis['avg_up'], analysis['u_avg_up'],\
        analysis['avg_down'],analysis['u_avg_down'],\
        analysis['avg_updown_mean'],analysis['u_avg_updown_mean']]])
    return datapath,avg_values


def get_avg_array_for_Vs(Vs,Lnr,day_dir,tag=''):
    for i,V in enumerate(Vs):
        datapath,avg_values=get_avg_values_from_json(V,Lnr,day_dir,tag=tag)
        if i==0:
            datapath,avg_array=get_avg_values_from_json(V,Lnr,day_dir,tag=tag)
        else:
            avg_array=np.append(avg_array,avg_values,axis=0)
    return datapath,avg_array

def get_LT_avg_values_from_json(V,L_dir,nr_files=2):
    string_name = 'V'+str(V)
    #print string_name
    datapath = os.path.join(L_dir,string_name)
    analysis=0
    
    #print datapath
    lws=np.array([])
    for filename in os.listdir(datapath):
        if filename.endswith("analysis.json"):
            #print filename 
            f = open(os.path.join(datapath,filename), 'r')
            analysis = json.load(f)
            f.close()
            lws = np.append(lws,analysis['lws'])
    return datapath,np.average(lws),np.std(lws)

def get_LT_avg_array_from_json(Vs,L_dir):
    for i,V in enumerate(Vs):
        datapath,avg_values,std_values=get_LT_avg_values_from_json(V,L_dir)
        if i==0:
            datapath,avg_array,std_array=get_LT_avg_values_from_json(V,L_dir)
        else:
            avg_array=np.append(avg_array,avg_values)
            std_array=np.append(std_array,std_values)
    return datapath,avg_array,std_array





##########functions to plot linewidths/Finesse vs cavity length

def plot_updown_lws(avgs,Ns,ax=None):
    if ax==None:
        fig,ax = plt.subplots()
    ax.errorbar(Ns,avgs[:,0],yerr=avgs[:,1],fmt='ro',label='ascending slope')
    ax.errorbar(Ns,avgs[:,2],yerr=avgs[:,3],fmt='go',label='descending slope')
    ax.errorbar(Ns,avgs[:,4],yerr=avgs[:,5],fmt='bo',label='up and down mean')
    ax.set_xlabel('cavity length (lambda/2)')
    ax.set_ylabel('linewidth (GHz)')
    ax.legend()
    return(ax)

def plot_mean_lws_from_avgs(avgs,Ns,ax=None,style='o',label='mean'):
    plot_mean_lws(avgs[:,4],avgs[:,5],Ns,ax=ax,style=style,label=label)
    return(ax)

def plot_mean_lws(avg_lws,u_avg_lws,Ns,ax=None,style='o',label='mean'):
    if ax==None:
        fig,ax = plt.subplots()
    ax.errorbar(Ns,avg_lws,yerr=u_avg_lws,fmt=style,label=label)
    ax.set_xlabel('cavity length (lambda/2)')
    ax.set_ylabel('linewidth (GHz)')
    ax.legend()
    return(ax)

def plot_finesse_vs_length_from_avgs(avgs,Ns,lambda_c=636.6e-9,ax=None,label='finesse',style='o'):
    ax = plot_finesse_vs_length(avgs[:,4],avgs[:,5],Ns,lambda_c=lambda_c,ax=ax,label=label,style=style)
    return(ax)
    
def plot_finesse_vs_length(avg_lws,u_avg_lws,Ns,lambda_c=636.6e-9,ax=None,label='finesse',style='o'):
    Finesses = c/(np.multiply(Ns,lambda_c))/np.multiply(avg_lws,1.e9)
    u_Finesses = c/(np.multiply(Ns,lambda_c))/np.multiply((avg_lws**2),1.e9)*u_avg_lws
    if ax==None:
        fig,ax = plt.subplots()
#     print Finesses,u_Finesses
    ax.errorbar(Ns,Finesses,yerr=u_Finesses,fmt=style,label=label)
    ax.set_xlabel('cavity length (lambda/2)')
    ax.set_ylabel('Finesse')
    ax.legend()
    return(ax)

def fit_lws(linewidths,u_linewidths,Ns,g_Finesse,datapath,ax=None,lambda_c=636.6e-9):
    #cavity_lengths =  Ns#636.6e-9*Ns/2 # in um
    g_a = 0#0.1*1.e9
    g_b = (c)/(2*g_Finesse)*1.e-9*2/lambda_c #*1.e9
    print g_b
    print Ns
    print linewidths
    p0, fitfunc, fitfunc_str = common.fit_inverse(g_a, g_b)
    fit_result = fit.fit1d(Ns,linewidths, None, 
                           p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=[0])
    if ax == None:
        fig,ax = plt.subplots(figsize=(8,4))
        ax.errorbar(Ns,linewidths,yerr=u_linewidths,fmt = 'bo')
    x = np.linspace(min(Ns),max(Ns),10*len(Ns))
    plot.plot_fit1d(fit_result,x, ax=ax, label='Fit',show_guess=True, 
                    plot_data=False,color='red', print_info= False)
    b_fit=fit_result['params_dict']['b']#cav length in um, linewidth in GHz
    u_b_fit = fit_result['error_dict']['b']
    print b_fit
    Finesse = c/(2*b_fit)*1.e-9*2/lambda_c
    u_Finesse = c/(2*b_fit**2)*1.e-9*2/lambda_c*u_b_fit
    ax.set_title(datapath)
    ax.legend(title='Finesse = %d $\pm$ %d'%(Finesse,u_Finesse))
    plt.savefig(os.path.join(datapath,'plot_lw_vs_L.png'))
    plt.show()
    plt.close()
    return Finesse



























