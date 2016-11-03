### Functions for analysing oscilloscope data
##SvD
###
import scipy
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
    analysis['u_avg_up']=scipy.stats.sem(lws_up)
    analysis['avg_down']=np.average(lws_down)
    analysis['u_avg_down']=scipy.stats.sem(lws_down)
    analysis['avg_updown_mean']=np.average(updown_mean_lws)
    analysis['u_avg_updown_mean']=scipy.stats.sem(updown_mean_lws)
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



########functions to make a nice fit and plot of single linewidth data.
###also

def fit_single(x_single,conversion=1.0):
    print g_dx*conversion
    print g_x01*conversion
    p0, fitfunc, fitfunc_str = common.fit_3lorentz_symmetric(g_a1, g_A1*conversion, g_x01*conversion, g_gamma1*conversion, g_dx*conversion, g_A2*conversion)
    fit_result = fit.fit1d(x_single,y_single, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)

    dx = fit_result['params_dict']['dx']
    gamma = fit_result['params_dict']['gamma1']
    u_gamma = fit_result['error_dict']['gamma1']

    x0 = fit_result['params_dict']['x01']
    A1 = fit_result['params_dict']['A1']
    A2 = fit_result['params_dict']['A2']
    a1 = fit_result['params_dict']['a1']
    #print dx
    #print dx,gamma,x0,A1,A2,a1
    return dx,gamma,u_gamma,x0,A1,A2,a1,fit_result


def plot_single(x_single,conversion=1.0,x_name='datapoints',x_label=''):
    dx,gamma,u_gamma,x0,A1,A2,a1,fit_result = fit_single(x_single,conversion=conversion)
    gamma_GHz = gamma/dx*6.
    u_gamma_GHz = u_gamma/dx*6.
    dx_GHz = 6.
    x0_GHz = 0.
    print 'cavity linewidth is',gamma_GHz,'+-',u_gamma_GHz
    x_single_GHz = (x_single-x0)/dx*6. #in GHz-> hardcodes a 6 GHz EOM frequency
    fig,ax = plt.subplots()
    ax.plot(x_single_GHz,y_single, 'o',color='orange',linewidth=4,label='data',zorder=1)
    xs_fit = np.linspace(x_single[0],x_single[-1],10*len(x_single))
    xs_fit_GHz = np.linspace(x_single_GHz[0],x_single_GHz[-1],10*len(x_single_GHz))
    ax.plot(xs_fit_GHz, fit_result['fitfunc'](xs_fit), color = 'darkblue', lw=1, label='fit',zorder=2 ) 

    halfmax = (A1/(math.pi*gamma)+a1)
    
    head_width = 0.01*rel_y_conversion
    head_length=12.5*conversion*rel_x_conversion/dx*6.
    arrow_color = 'darkblue'
    ax.arrow(-dx_GHz,halfmax/2,dx_GHz-head_length,0, width = head_width/5.,head_width=head_width, head_length=head_length,  fc=arrow_color,ec=arrow_color)
    ax.arrow(0,halfmax/2,-dx_GHz+head_length,0, width = head_width/5.,head_width=head_width, head_length=head_length, fc=arrow_color, ec=arrow_color)
    ax.text(-3*dx_GHz/4,halfmax/2+0.01*rel_y_conversion,'6 GHz',color=arrow_color)
    
    head_length2=7.*conversion*rel_x_conversion/dx*6.
    arrow_length2 =20*conversion*rel_x_conversion/dx*6.
    arrow_color2 = 'darkblue'
    #print x0-gamma/2.
    ax.arrow(-gamma_GHz/2.-arrow_length2,halfmax,arrow_length2-head_length2,0, width = head_width/5.,head_width=head_width, head_length=head_length2,  fc=arrow_color2,ec=arrow_color2)
    ax.arrow(gamma_GHz/2.+arrow_length2,halfmax,-arrow_length2+head_length2,0,width = head_width/5., head_width=head_width, head_length=head_length2, fc=arrow_color2, ec=arrow_color2)
    ax.text(x0+0.03*rel_x_conversion/dx*6.,halfmax+0.01*rel_y_conversion,'FWHM',color=arrow_color)
        
    #ax.legend()
    ax.set_xlabel('detuning (GHz)',fontsize=14)
    ax.set_ylabel('transmitted signal (V)',fontsize=14)
    ax.set_xlim([x_single_GHz[0],x_single_GHz[-1]])
    

    print x_single_GHz[0], x_single_GHz[-1]
    #print x_in_GHz 

    ax2= ax.twiny()
    ax2.set_xlim([x_single[0],x_single[-1]])
    ax2.set_xlabel(x_label,fontsize=14)


    plt.savefig(os.path.join(indir,filename+'single_example_lw_plot_vs_%s.eps'%(x_name)))
    plt.savefig(os.path.join(indir,filename+'single_example_lw_plot_vs_%s.png'%(x_name)))
    plt.show()
    return ax



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


def get_values_from_json_calc_sem(V,day_dir):
    V_str = '%1.3f'%V
    string_name = 'V'+str(V)
    datapath = os.path.join(day_dir,string_name)
    analysis=0
    
#     print datapath
    for filename in os.listdir(datapath):
        if filename.endswith("analysis.json"):
            use_filename=filename
#             print filename 
            f = open(os.path.join(datapath,filename), 'r')
            analysis = json.load(f)
            f.close
    if analysis==0:
        print 'analysis file not found!'
        return datapath, 0

    updownmeans = analysis['updown_mean_lws']
    u_avg_updownmeans  = scipy.stats.sem(updownmeans)
    u_avg_down = scipy.stats.sem(analysis['lws_down'])
    u_avg_up = scipy.stats.sem(analysis['lws_up'])
    analysis['u_avg_updown_mean']=u_avg_updownmeans
    analysis['u_avg_down'] = u_avg_down
    analysis['u_avg_up'] = u_avg_down

    f2 = open(os.path.join(datapath,use_filename), 'w')
    json.dump(analysis,f2,default=default)
    f2.close()

    return analysis

def get_values_from_json_calc_sem_for_L(Vs,day_dir):
    for V in Vs:
        get_values_from_json_calc_sem(V,day_dir)

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
    return datapath,np.average(lws),scipy.stats.sem(lws)

def get_LT_avg_array_from_json(Vs,L_dir):
    for i,V in enumerate(Vs):
        datapath,avg_values,u_avg_values=get_LT_avg_values_from_json(V,L_dir)
        if i==0:
            datapath,avg_array,u_avg_array=get_LT_avg_values_from_json(V,L_dir)
        else:
            avg_array=np.append(avg_array,avg_values)
            u_avg_array=np.append(u_avg_array,u_avg_values)
    return datapath,avg_array,u_avg_array





##########functions to plot linewidths/Finesse vs cavity length

def get_averaged_finesses(avg_lws,u_avg_lws,Ns,lambda_c=636.6e-9):
    finesses, u_finesses = finesses_from_lws(avg_lws,u_avg_lws,Ns,lambda_c = lambda_c)#,lamba_c = 636.6e-9)
    #print np.histogram(Ns)
    bins = np.histogram(Ns)
    range_Ns = np.arange(min(Ns),max(Ns))
    
    averaged_finesse=np.zeros(len(range_Ns))
    u_averaged_finesse = np.zeros(len(range_Ns))
    for i,N in enumerate(range_Ns):
        avg_finesse,u_avg_finesse =get_finesse_at_N(avg_lws,u_avg_lws,Ns,N)
        averaged_finesse[i]=avg_finesse
        u_averaged_finesse[i] = u_avg_finesse
    return averaged_finesse,u_averaged_finesse,range_Ns

def get_finesse_at_N(avg_lws,u_avg_lws,Ns,N,lambda_c = 636.6e-9):
    finesses,u_finesses = finesses_from_lws(avg_lws,u_avg_lws,Ns,lambda_c=lambda_c)
    finesses_at_N = finesses[np.where(Ns==N)]
    finesse = np.average(finesses_at_N)
    #print finesses_at_N
    #print 'finesse = ',finesse

    if len(finesses_at_N)==1:
        u_finesse = u_finesses[np.where(Ns==N)][0]
    else:
        u_finesse = scipy.stats.sem(finesses_at_N)
        #print 'sem= ',scipy.stats.sem(finesses_at_N)

    return finesse,u_finesse

def get_finesse_at_Ns(avg_lws,u_avg_lws,all_Ns,Ns):
    finesse_Ns = np.zeros(len(Ns))
    u_finesse_Ns = np.zeros(len(Ns))
    for i,N in enumerate(Ns):
        f_i,u_f_i = get_finesse_at_N(avg_lws,u_avg_lws,all_Ns,N)
        finesse_Ns[i],u_finesse_Ns[i] = f_i,u_f_i

    #print np.where(~np.isnan(finesse_Ns))
    finesse_Ns = finesse_Ns[np.where(~np.isnan(finesse_Ns))]

    #print finesse_Ns
    avg_finesse_Ns = np.average(finesse_Ns)
    u_avg_finesse_Ns = scipy.stats.sem(finesse_Ns)

    return avg_finesse_Ns,u_avg_finesse_Ns

def finesses_from_lws(avg_lws,u_avg_lws,Ns,lambda_c=636.6e-9):
    finesses = c/(np.multiply(Ns,lambda_c))/np.multiply(avg_lws,1.e9)
    u_finesses = c/(np.multiply(Ns,lambda_c))/np.multiply((avg_lws**2),1.e9)*u_avg_lws
    return finesses,u_finesses

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
    # Finesses = c/(np.multiply(Ns,lambda_c))/np.multiply(avg_lws,1.e9)
    # u_Finesses = c/(np.multiply(Ns,lambda_c))/np.multiply((avg_lws**2),1.e9)*u_avg_lws
    Finesses,u_Finesses=finesses_from_lws(avg_lws,u_avg_lws,Ns)
    if ax==None:
        fig,ax = plt.subplots()
    ax.errorbar(Ns,Finesses,yerr=u_Finesses,fmt=style,label=label)
    ax.set_xlabel('cavity length (lambda/2)')
    ax.set_ylabel('Finesse')
    ax.legend()
    return(ax)

def plot_avgd_finesse_vs_length_from_avgs(avgs,Ns,lambda_c=636.6e-9,ax=None,label='finesse',style='o',xerror=False):
    ax = plot_avgd_finesse_vs_length(avgs[:,4],avgs[:,5],Ns,lambda_c=lambda_c,ax=ax,label=label,style=style,xerror=xerror)
    return(ax)

def plot_avgd_finesse_vs_length(avg_lws,u_avg_lws,Ns,lambda_c=636.6e-9,ax=None,label='finesse',style='o',xerror=False):
    finesses,u_finesses,Ns = get_averaged_finesses(avg_lws,u_avg_lws,Ns,lambda_c=636.6e-9)
    if xerror==True:
        xerr = 1.
    else:
        xerr = None
    if ax==None:
        fig,ax = plt.subplots()
    ax.errorbar(Ns,finesses,yerr=u_finesses,xerr=xerr,fmt=style,label=label)
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




if __name__=="__main__":

    # ###################values for single sweep at RT for cahraacterisation paper uncomment to use
    # #indir = 'K:/ns\qt\Diamond\Projects\Cavities\Cavity characterisation paper\data\data_for_cav_char_paper/20160822\RT_OND_POS1_L3\V3'
    # #NICE_LWS00076, second lw at datapoint 1627
    # #the sweep has 500 ms in 125000 data points, so we can convert the x-axis:
    # indir = 'K:/ns\qt\Diamond\Projects\Cavities\Cavity characterisation paper\data\data_for_cav_char_paper/20160822\RT_OND_POS1_L3\V3'
    # filename = 'NICE_LWS00076'
    # f = open(os.path.join(indir,filename+'single_example_lw_data.txt'), 'r')
    # x_single,y_single = np.loadtxt(f)
    # f.close()

    # conversion = 500/125000.

    # x_single=x_single-x_single[0] #start with datapoint 0
    # t_single = x_single*conversion #ms

    # g_a1 = 0.014
    # g_A1 = 4.5
    # g_x01= 160
    # g_gamma1=10
    # g_dx=89
    # g_A2 = 2.0
    # fixed =[]

    # rel_x_conversion=1.
    # rel_y_conversion=1.




    ###################values for singel sweep at LT for charactaritasiotn paper. uncomment to use
    indir = 'K:/ns\qt\Diamond\Projects\Cavities\Cavity characterisation paper\data\data_for_cav_char_paper/20160905/LT_OND_POS2_L6/V2'
    filename = 'LT_LWS024'  #first lw at 246
    #the sweep has 500 ms in 125000 data points, so we can convert the x-axis:

    f = open(os.path.join(indir,filename+'single_example_lw_data.txt'), 'r')
    x_single,y_single = np.loadtxt(f)
    f.close()

    conversion = 200/125000.

    cut_nr_datapoints=45
    x_single=x_single-x_single[0] #start with datapoint 0
    x_single=x_single[cut_nr_datapoints:-cut_nr_datapoints]
    y_single=y_single[cut_nr_datapoints:-cut_nr_datapoints]

    t_single = x_single*conversion #ms

    g_a1 = 0.01
    g_A1 = 0.6
    g_x01= 100
    g_gamma1=6
    g_dx=26
    g_A2 = 0.2
    fixed =[]

    rel_x_conversion=200/500.
    rel_y_conversion=0.08/0.3



    plot_single(t_single,conversion=conversion,x_name='time',x_label='time in scan (ms)')
    plot_single(x_single,x_name='datapoints',x_label='datapoints')













