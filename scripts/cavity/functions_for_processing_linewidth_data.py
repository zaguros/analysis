### Functions for analysing oscilloscope data
##SvD
###
import scipy
import numpy as np 
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
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

def fit_single(x,y,conversion=1.0):
    #print g_dx*conversion
    #print 'g_x01',g_x01*conversion
    #print 'x',x 
    p0, fitfunc, fitfunc_str = common.fit_3lorentz_symmetric(g_a1, g_A1*conversion, g_x01*conversion, g_gamma1*conversion, g_dx*conversion, g_A2*conversion)
    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)

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


def plot_single(x,y,conversion=1.0,x_name='datapoints',x_label='',rel_x_conversion=1.,rel_y_conversion=1.,ret_data=False):
    dx,gamma,u_gamma,x0,A1,A2,a1,fit_result = fit_single(x,y,conversion=conversion)
    gamma_GHz = gamma/dx*6.
    u_gamma_GHz = u_gamma/dx*6.
    dx_GHz = 6.
    x0_GHz = 0.
    #print 'conversion',conversion
    #print 'rel_x_conversion',rel_x_conversion
    print 'cavity linewidth is',gamma_GHz,'+-',u_gamma_GHz
    x_GHz = (x-x0)/dx*6. #in GHz-> hardcodes a 6 GHz EOM frequency
    fig,ax = plt.subplots()
    ax.plot(x_GHz,y, 'o',color='orange',linewidth=4,label='data',zorder=1)
    xs_fit = np.linspace(x[0],x[-1],10*len(x))
    xs_fit_GHz = np.linspace(x_GHz[0],x_GHz[-1],10*len(x_GHz))
    ys_fit = fit_result['fitfunc'](xs_fit)
    ax.plot(xs_fit_GHz, ys_fit , color = 'darkblue', lw=1, label='fit',zorder=2 ) 

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
    ax.set_xlim([x_GHz[0],x_GHz[-1]])
    

    print x_GHz[0], x_GHz[-1]
    #print x_in_GHz 

    ax2= ax.twiny()
    ax2.set_xlim([x[0],x[-1]])
    ax2.set_xlabel(x_label,fontsize=14)


    plt.savefig(os.path.join(indir,filename+'single_example_lw_plot_vs_%s.eps'%(x_name)))
    plt.savefig(os.path.join(indir,filename+'single_example_lw_plot_vs_%s.png'%(x_name)))
    plt.show()

    if ret_data:
        return x,x_GHz, y,xs_fit_GHz,ys_fit
    return ax


def plot_both(x01,x1,y1,xfit1,yfit1,x02,x2,y2,xfit2,yfit2,dy2=0.):
    #scale y2
    normalise = (max(y1)-min(y1))/(max(y2)-min(y2))
    y2 = y2 * normalise
    yfit2 = yfit2 * normalise

    plt.figure(figsize=[4,4])

    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(top=0.75)

    par1 = host.twiny()
    par2 = host.twiny()

    offset = 60
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["top"] = new_fixed_axis(loc="top",
                                        axes=par2,
                                        offset=(0, offset))

    par2.axis["top"].toggle(all=True)

    host.set_xlim([x1[0],x1[-1]])

    par1.set_xlim([x01[0]-x01[0],x01[-1]-x01[0]])
    par2.set_xlim([x02[0]-x02[0],x02[-1]-x02[0]])

    host.set_xlabel('detuning (GHz)',fontsize=14)
    host.set_ylabel('transmitted signal (V)',fontsize=14)
    par1.set_xlabel('scan time 1 (ms)',fontsize=14)
    par2.set_xlabel('scan time 2 (ms)',fontsize=14)

    host.plot(x1,y1, 'o',color='orange',linewidth=4,label='1',zorder=1)
    host.plot(xfit1,yfit1, '-',color='orange',linewidth=4,zorder=1)
    host.plot(x2,y2+dy2, 'o',color='blue',linewidth=4,label='2',zorder=1)
    host.plot(xfit2,yfit2+dy2, '-',color='blue',linewidth=4,zorder=1)
    host.legend()


    # fig,ax = plt.subplots()


    # ax.plot(x1,y1, 'o',color='orange',linewidth=4,label='data 1',zorder=1)
    # ax.plot(xfit1,yfit1, '-',color='orange',linewidth=4,label='fit 1',zorder=1)
    # ax.plot(x2,y2+dy2, 'o',color='blue',linewidth=4,label='data 2',zorder=1)
    # ax.plot(xfit2,yfit2+dy2, '-',color='blue',linewidth=4,label='fit 2',zorder=1)
   
    # #ax.legend()
    # ax.set_xlabel('detuning (GHz)',fontsize=14)
    # ax.set_ylabel('transmitted signal (V)',fontsize=14)
    # ax.set_xlim([x1[0],x1[-1]])
    
    # #print x_in_GHz 

    # ax2= ax.twiny()
    # ax2.set_xlim([x01[0],x01[-1]])
    # ax2.set_xlabel('scan time 1 (ms)',fontsize=14)

    # ax3= ax.twiny()
    # ax3.set_xlim([x02[0],x02[-1]])
    # ax3.set_xlabel('scan time 2 (ms)',fontsize=14)
    print indir 
    plt.savefig(os.path.join(indir,filename+'both_example_lw_plot.eps'))
    plt.savefig(os.path.join(indir,filename+'both_example_lw_plot.png'))

    return host



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


def plot_finesse_vs_airlikecharacter((RT_OFFD_POS0_avgs,RT_OFFD_POS0_Ns),(LT_OFFD_POS3_avg_array,LT_OFFD_POS3_u_avg_array,LT_OFFD_POS3_Ns),
            (RT_OND_POS0_avgs,RT_OND_POS0_Ns),(RT_OND_POS1_avgs,RT_OND_POS1_Ns),(RT_OND_POS2_avgs,RT_OND_POS2_Ns),
            (LT_OND_POS2_avg_array,LT_OND_POS2_u_avg_array,LT_OND_POS2_Ns),(LT_OND_POS4_avg_array,LT_OND_POS4_u_avg_array,LT_OND_POS4_Ns)):
    Ns= np.arange(47,56)
    #these things require global variables of 
    avg_finesse_RT_OFFD = np.average(finesses_from_lws(RT_OFFD_POS0_avgs[:,4],RT_OFFD_POS0_avgs[:,5],RT_OFFD_POS0_Ns)[0])
    avg_finesse_LT_OFFD = finesses_from_lws(LT_OFFD_POS3_avg_array,LT_OFFD_POS3_u_avg_array,LT_OFFD_POS3_Ns)[0][0]
    u_avg_finesse_RT_OFFD = scipy.stats.sem(finesses_from_lws(RT_OFFD_POS0_avgs[:,4],RT_OFFD_POS0_avgs[:,5],RT_OFFD_POS0_Ns)[0])
    u_avg_finesse_LT_OFFD = finesses_from_lws(LT_OFFD_POS3_avg_array,LT_OFFD_POS3_u_avg_array,LT_OFFD_POS3_Ns)[1][0]

    F_RT_POS0_N, u_F_RT_POS0_N = get_finesse_at_Ns(RT_OND_POS0_avgs[:,4],RT_OND_POS0_avgs[:,5],RT_OND_POS0_Ns,Ns)
    F_RT_POS1_N, u_F_RT_POS1_N = get_finesse_at_Ns(RT_OND_POS1_avgs[:,4],RT_OND_POS1_avgs[:,5],RT_OND_POS1_Ns,Ns)
    F_RT_POS2_N, u_F_RT_POS2_N = get_finesse_at_Ns(RT_OND_POS2_avgs[:,4],RT_OND_POS2_avgs[:,5],RT_OND_POS2_Ns,Ns)
    F_LT_POS2_N, u_F_LT_POS2_N = get_finesse_at_Ns(LT_OND_POS2_avg_array,LT_OND_POS2_u_avg_array,LT_OND_POS2_Ns,Ns)
    F_LT_POS4_N, u_F_LT_POS4_N = get_finesse_at_Ns(LT_OND_POS4_avg_array,LT_OND_POS4_u_avg_array,LT_OND_POS4_Ns,Ns)

    ####BELOW IS HARDCODED THE 'airlike character' of the modes as obtained for 2d plots in spectrometer analysis!
    RT_OND_POS0_ans=np.array([0.33,0.32,0.35])*2
    RT_OND_POS1_ans=np.array([0.29,0.26,0.25,0.25])*2
    RT_OND_POS2_ans=np.array([0.42,0.43,0.39,0.38,0.35,0.35,0.35])*2
    LT_OND_POS2_ans=np.array([0.14,0.19,0.1])*2
    LT_OND_POS4_ans=np.array([0.36,0.38,0.4,0.38])*2


    RT_finesses = np.array([F_RT_POS0_N,F_RT_POS1_N,F_RT_POS2_N])
    u_RT_finesses = np.array([u_F_RT_POS0_N,u_F_RT_POS1_N,u_F_RT_POS2_N])
    RT_airnesses= [0.30,0.25,0.38]
    RT_airnesses =[np.average(an) for an in [RT_OND_POS0_ans,RT_OND_POS1_ans,RT_OND_POS2_ans]]
    u_RT_airnesses =[scipy.stats.sem(an) for an in [RT_OND_POS0_ans,RT_OND_POS1_ans,RT_OND_POS2_ans]]

    LT_finesses = np.array([F_LT_POS2_N,F_LT_POS4_N])
    u_LT_finesses = np.array([u_F_LT_POS2_N,u_F_LT_POS4_N])
    LT_airnesses= [0.1,0.38]
    LT_airnesses =[np.average(an) for an in [LT_OND_POS2_ans,LT_OND_POS4_ans]]
    u_LT_airnesses = [scipy.stats.sem((an)) for an in [LT_OND_POS2_ans,LT_OND_POS4_ans]]

    fig,ax = plt.subplots()

    ax.errorbar([1],avg_finesse_RT_OFFD,yerr=u_avg_finesse_RT_OFFD,fmt='ro',label='RT OFFD')
    ax.errorbar([1],avg_finesse_LT_OFFD,yerr=u_avg_finesse_LT_OFFD,fmt='bo',label='LT OFFD')
    ax.errorbar(RT_airnesses,RT_finesses,xerr=u_RT_airnesses,yerr=u_RT_finesses,fmt='ro',label='RT')
    ax.errorbar(LT_airnesses,LT_finesses,xerr=u_LT_airnesses,yerr=u_LT_finesses,fmt='bo',label='LT')
    ax.legend(loc='upper left')
    ax.set_ylabel('average finesse N %d - %d'%(Ns[0],Ns[-1]))
    ax.set_xlabel('relative air-like character')
    ax.set_xlim([0,1.05])
    ax.set_ylim([0,30000])

    return ax

if __name__=="__main__":

    ###################values for single sweep at RT for cahraacterisation paper uncomment to use
    #indir = 'K:/ns\qt\Diamond\Projects\Cavities\Cavity characterisation paper\data\data_for_cav_char_paper/20160822\RT_OND_POS1_L3\V3'
    #NICE_LWS00076, second lw at datapoint 1627
    #the sweep has 500 ms in 125000 data points, so we can convert the x-axis:
    indir = 'K:/ns\qt\Diamond\Projects\Cavities\Cavity characterisation paper\data\data_for_cav_char_paper/20160822\RT_OND_POS1_L3\V3'
    filename = 'NICE_LWS00076'
    f = open(os.path.join(indir,filename+'single_example_lw_data.txt'), 'r')
    x_single,y_single = np.loadtxt(f)
    f.close()

    conversion = 500/125000.

    x_single=x_single-x_single[0] #start with datapoint 0
    t_single = x_single*conversion #ms

 

    rel_x_conversion=1.
    rel_y_conversion=1.




    ###################values for singel sweep at LT for charactaritasiotn paper. uncomment to use
    indir = 'K:/ns\qt\Diamond\Projects\Cavities\Cavity characterisation paper\data\data_for_cav_char_paper/20160905/LT_OND_POS2_L6/V2'
    filename = 'LT_LWS024'  #first lw at 246
    #the sweep has 500 ms in 125000 data points, so we can convert the x-axis:

    f_LT = open(os.path.join(indir,filename+'single_example_lw_data.txt'), 'r')
    x_single_LT,y_single_LT = np.loadtxt(f_LT)
    f_LT.close()

    conversion_LT = 200/125000.

    cut_nr_datapoints=45
    x_single_LT=x_single_LT-x_single_LT[0] #start with datapoint 0
    x_single_LT=x_single_LT[cut_nr_datapoints:-cut_nr_datapoints]
    y_single_LT=y_single_LT[cut_nr_datapoints:-cut_nr_datapoints]

    t_single_LT = x_single_LT*conversion_LT #ms

    rel_x_conversion_LT=200/500.
    rel_y_conversion_LT=0.08/0.3
    


    g_a1 = 0.014
    g_A1 = 4.5
    g_x01= 160
    g_gamma1=10
    g_dx=89
    g_A2 = 2.0
    fixed =[]


    x0_RT,x_RT,y_RT,xfit_RT,yfit_RT=plot_single(t_single,y_single,conversion=conversion,x_name='time',
        x_label='time in scan (ms)',rel_x_conversion=rel_x_conversion,rel_y_conversion=rel_y_conversion,ret_data=True)
    #plot_single(x_single,y_single,x_name='datapoints',x_label='datapoints')
    
    g_a1 = 0.01
    g_A1 = 0.6
    g_x01= 100
    g_gamma1=6
    g_dx=26
    g_A2 = 0.2
    fixed =[]
    #print conversion_LT
    #print rel_x_conversion_LT

    x0_LT,x_LT,y_LT,xfit_LT,yfit_LT=plot_single(t_single_LT,y_single_LT,conversion=conversion_LT,x_name='time',
        x_label='time in scan (ms)',rel_x_conversion=rel_x_conversion_LT,rel_y_conversion=rel_y_conversion_LT,ret_data=True)
    #plot_single(x_single_LT,y_single_LT,x_name='datapoints',x_label='datapoints')


    plot_both(x0_RT,x_RT,y_RT,xfit_RT,yfit_RT,x0_LT,x_LT,y_LT,xfit_LT,yfit_LT,dy2=0.3)










