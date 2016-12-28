### Functions to plot a single scan from the oscilloscope. 
### Used for figure 2(b) the cavity characterisation paper
### SvD  October 2016
###


from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

########functions to make a  fit and plot of a single scan of linewidth data (three lorentzians).
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


def plot_single(datadir,filename,x,y,conversion=1.0,x_name='datapoints',x_label='',rel_x_conversion=1.,rel_y_conversion=1.,ret_data=False,save_plot=False):
    dx,gamma,u_gamma,x0,A1,A2,a1,fit_result = fit_single(x,y,conversion=conversion)

    dx_GHz = 6.#in GHz-> hardcodes a 6 GHz EOM frequency
    gamma_GHz = gamma/dx*dx_GHz
    u_gamma_GHz = u_gamma/dx*dx_GHz
    x0_GHz = 0.

    #print 'conversion',conversion
    #print 'rel_x_conversion',rel_x_conversion
    print 'cavity linewidth is',gamma_GHz,'+-',u_gamma_GHz
    x_GHz = (x-x0)/dx*dx_GHz
    fig,ax = plt.subplots()
    ax.plot(x_GHz,y, 'o',color='orange',linewidth=4,label='data',zorder=1)
    xs_fit = np.linspace(x[0],x[-1],10*len(x))
    xs_fit_GHz = np.linspace(x_GHz[0],x_GHz[-1],10*len(x_GHz))
    ys_fit = fit_result['fitfunc'](xs_fit)
    ax.plot(xs_fit_GHz, ys_fit , color = 'darkblue', lw=1, label='fit',zorder=2 ) 

    halfmax = (A1/(math.pi*gamma)+a1)
    
    head_width = 0.01*rel_y_conversion
    head_length=12.5*conversion*rel_x_conversion/dx*dx_GHz
    arrow_color = 'darkblue'
    ax.arrow(-dx_GHz,halfmax/2,dx_GHz-head_length,0, width = head_width/5.,head_width=head_width, head_length=head_length,  fc=arrow_color,ec=arrow_color)
    ax.arrow(0,halfmax/2,-dx_GHz+head_length,0, width = head_width/5.,head_width=head_width, head_length=head_length, fc=arrow_color, ec=arrow_color)
    ax.text(-3*dx_GHz/4,halfmax/2+0.01*rel_y_conversion,'6 GHz',color=arrow_color)
    
    head_length2=7.*conversion*rel_x_conversion/dx*dx_GHz
    arrow_length2 =20*conversion*rel_x_conversion/dx*dx_GHz
    arrow_color2 = 'darkblue'
    #print x0-gamma/2.
    ax.arrow(-gamma_GHz/2.-arrow_length2,halfmax,arrow_length2-head_length2,0, width = head_width/5.,head_width=head_width, head_length=head_length2,  fc=arrow_color2,ec=arrow_color2)
    ax.arrow(gamma_GHz/2.+arrow_length2,halfmax,-arrow_length2+head_length2,0,width = head_width/5., head_width=head_width, head_length=head_length2, fc=arrow_color2, ec=arrow_color2)
    ax.text(x0+0.03*rel_x_conversion/dx*dx_GHz,halfmax+0.01*rel_y_conversion,'FWHM',color=arrow_color)
        
    #ax.legend()
    ax.set_xlabel('detuning (GHz)',fontsize=14)
    ax.set_ylabel('transmitted signal (V)',fontsize=14)
    ax.set_xlim([x_GHz[0],x_GHz[-1]])
    

    #print x_GHz[0], x_GHz[-1]
    #print x_in_GHz 

    ax2= ax.twiny()
    ax2.set_xlim([x[0],x[-1]])
    ax2.set_xlabel(x_label,fontsize=14)

    if save_plot:
        plt.savefig(os.path.join(indir,filename+'single_example_lw_plot_vs_%s.eps'%(x_name)))
        plt.savefig(os.path.join(indir,filename+'single_example_lw_plot_vs_%s.png'%(x_name)))
    plt.show()

    if ret_data:
        return x,x_GHz, y,xs_fit_GHz,ys_fit
    return ax

#########Function specifically to plot two single scans above each other. 
#########used for the cavity characterisation paper to plot both RT and LT linewidth.
def plot_both(datadir,filename,x01,x1,y1,xfit1,yfit1,x02,x2,y2,xfit2,yfit2,dy2=0.,save_plot=False):
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

    if save_plot:
        print indir 
        plt.savefig(os.path.join(indir,filename+'both_example_lw_plot.eps'))
        plt.savefig(os.path.join(indir,filename+'both_example_lw_plot.png'))

    return host



if __name__=="__main__":
    ###################values for single sweep at RT for cahraacterisation paper 
    #indir = 'K:/ns\qt\Diamond\Projects\Cavities\Cavity characterisation paper\data\data_for_cav_char_paper/20160822\RT_OND_POS1_L3\V3'
    #NICE_LWS00076, second lw at datapoint 1627
    #the sweep has 500 ms in 125000 data points, so we can convert the x-axis:
    datadir_RT = 'K:/ns\qt\Diamond\Projects\Cavities\Cavity characterisation paper\data\data_for_cav_char_paper/20160822\RT_OND_POS1_L3\V3'
    filename_RT = 'NICE_LWS00076'
    f = open(os.path.join(datadir_RT,filename_RT+'single_example_lw_data.txt'), 'r')
    x_single,y_single = np.loadtxt(f)
    f.close()

    conversion = 500/125000.
    x_single = x_single-x_single[0] #start couting at datapoint 0
    t_single = x_single*conversion #ms
    rel_x_conversion=1.
    rel_y_conversion=1.


    ###################values for singel sweep at LT for charactaritasiotn paper. 
    datadir_LT = 'K:/ns\qt\Diamond\Projects\Cavities\Cavity characterisation paper\data\data_for_cav_char_paper/20160905/LT_OND_POS2_L6/V2'
    filename_LT = 'LT_LWS024'  #first lw at 246
    #the sweep has 500 ms in 125000 data points, so we can convert the x-axis:

    f_LT = open(os.path.join(datadir_LT,filename_LT+'single_example_lw_data.txt'), 'r')
    x_single_LT,y_single_LT = np.loadtxt(f_LT)
    f_LT.close()

    conversion_LT = 200/125000.

    cut_nr_datapoints=45
    x_single_LT=x_single_LT-x_single_LT[0] #start counting at datapoint 0
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


    x0_RT,x_RT,y_RT,xfit_RT,yfit_RT=plot_single(datadir_RT,filename_RT,t_single,y_single,conversion=conversion,x_name='time',
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

    x0_LT,x_LT,y_LT,xfit_LT,yfit_LT=plot_single(datadir_LT,filename_LT,t_single_LT,y_single_LT,conversion=conversion_LT,x_name='time',
        x_label='time in scan (ms)',rel_x_conversion=rel_x_conversion_LT,rel_y_conversion=rel_y_conversion_LT,ret_data=True)
    #plot_single(x_single_LT,y_single_LT,x_name='datapoints',x_label='datapoints')


    plot_both(datadir_RT,filename_RT,x0_RT,x_RT,y_RT,xfit_RT,yfit_RT,x0_LT,x_LT,y_LT,xfit_LT,yfit_LT,dy2=0.3)
