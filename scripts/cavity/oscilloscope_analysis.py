### Functions for analysing oscilloscope data

import numpy as np 
from matplotlib import pyplot as plt




def get_avg_values_from_json(V,Lnr,day_dir,tag=''):
    V_str = '%1.3f'%V
    string_name = 'V'+str(V)
    #string_name = '%sL%d_V%sp%s'%(tag,Lnr,V_str[0],V_str[2:])
    print string_name
    datapath = os.path.join(day_dir,string_name)
    analysis=0
    
    print datapath
    for filename in os.listdir(datapath):
        if filename.endswith("analysis.json"):
            print filename 
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
    print string_name
    datapath = os.path.join(L_dir,string_name)
    analysis=0
    
    print datapath
    lws=np.array([])
    for filename in os.listdir(datapath):
        if filename.endswith("analysis.json"):
            print filename 
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
  










def plot_updown_lws(avgs,Ns,ax=None):
    if ax==None:
        fig,ax = plt.subplots()
    ax.errorbar(Ns,avgs[:,0],yerr=avgs[:,1],fmt='ro',label='upgoing slope')
    ax.errorbar(Ns,avgs[:,2],yerr=avgs[:,3],fmt='go',label='downgoing slope')
    ax.errorbar(Ns,avgs[:,4],yerr=avgs[:,5],fmt='bo',label='up and down mean')
    ax.set_xlabel('cavity length (lambda/2)')
    ax.set_ylabel('linewidth (GHz)')
    ax.legend()
    return(ax)

def plot_mean_lws(avgs,Ns,ax=None,style='bo',label='mean'):
    if ax==None:
        fig,ax = plt.subplots()
    ax.errorbar(Ns,avgs[:,4],yerr=avgs[:,5],fmt=style,label=label)
    ax.set_xlabel('cavity length (lambda/2)')
    ax.set_ylabel('linewidth (GHz)')
    ax.legend()
    return(ax)



def plot_finesse_vs_length_from_avgs(avgs,Ns,lambda_c=636.6e-9):
    ax = plot_finesse_vs_length(avgs[:,4],avgs[:,5],Ns,lambda_c=lambda_c)

def plot_finesse_vs_length(avg_lws,u_avg_lws,Ns,lambda_c=636.6e-9):
    Finesses = c/(np.multiply(Ns,lambda_c))/np.multiply(avg_lws,1.e9)
    u_Finesses = c/(np.multiply(Ns,lambda_c))/np.multiply((avg_lws**2),1.e9)*u_avg_lws
    fig,ax = plt.subplots()
    print Finesses,u_Finesses
    ax.errorbar(Ns,Finesses,yerr=u_Finesses,fmt='bo',label='mean Finesse')
    ax.set_xlabel('cavity length (lambda/2)')
    ax.set_ylabel('Finess')
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













































indir = 'D:\measuring\data/20160905/LT_OND_POS2_L6/V1'
filename = 'LT_LWS025'
EOM_freq= 6.
reload(oscillo_ana)

oa_all = oscillo_ana.oscilloscope_analysis(indir=indir,filename=filename)
oa_all.get_data(use_timetrace=False,x_min=None,x_max=8000,data_col=3)
oa_all.plot_data(plot_mod=False,figsize=(16,4),fmt='o')

x_min=64000
x_max=x_min+5000
oa_single = oscillo_ana.oscilloscope_analysis(indir=indir,filename=filename)
# oa_single.get_data(use_timetrace=False,x_min=x_min,x_max=x_max,data_col=3)#1400)
# oa_single.plot_data(plot_mod=False,figsize=(16,4),fmt='o')
#

lws = np.array([])#lws_old## 
As = np.array([])#lws_old## 
gammas = np.array([])#lws_old## 
chisqs= np.array([])
x01s = np.array([])
u_lws =np.array([])
nr_lws=400
x0 =275
success =np.zeros(nr_lws)
windowsize=120# print len(x)

for i in np.arange(nr_lws):
#     if i>10:
#         break
    xi = int(x0+(len(x)*i/nr_lws))

    if i%2==1:
        xi=xi
    if (xi-windowsize)>0:
        x_min = xi-windowsize
    else: 
        x_min = 0
    if (xi+windowsize)<len(x)-1:
        x_max = xi+windowsize
    else:
        x_max = -1
    

    oa_single.get_data(use_timetrace=False,x_min = x_min, x_max=x_max,data_col=3)    
#     oa_single.plot_data()

    lw,u_lw,A1,u_A1,gamma1,u_gamma1,chisq,x01 = \
        oa_single.plot_and_fit_with_EOM(EOM_freq=EOM_freq,
           g_gamma1=8, g_dx=30,g_A2=0.25,g_A3=0.25,show_fit=False,print_fit=False)
        
    if chisq>0.2:
        print 20*'*'
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
    if (x01+windowsize)<len(x)-1:
        x_max = x01+windowsize
    else:
        x_max = None
    oa_single.get_data(use_timetrace=False,x_min = x_min, x_max=x_max,data_col=3)
#     oa_single.plot_data()    
    lw,u_lw,A1,u_A1,gamma1,u_gamma1,chisq,x02 = \
        oa_single.plot_and_fit_with_EOM(EOM_freq=EOM_freq,
           g_gamma1=8, g_dx=30,g_A2=0.2,g_A3=0.2,show_fit=True,print_fit=False,plot_name=str(int(round(x01,0))))
        
    u_lws = np.append(u_lws,u_lw)
    lws = np.append(lws,lw)
    As = np.append(As,A1)
    gammas = np.append(gammas,gamma1)
    chisqs = np.append(chisqs,chisq)
    x01s = np.append(x01s,x01)
    if chisq>0.00055:
        print 20*'*'
        print chisq,lw
        print 'fit not good enough! '
        print 20*'*'
        continue
    
    print 20*'*'
    print chisq,lw
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

print 'the good linewidth fits from %s :'%filename
print 'x0s = ',good_x01s
print 'lws = ',good_lws
print 'u_lws = ',good_u_lws
print 'chisq = ',good_chisqs

save_to_json_file(indir,filename+'_analysis',good_fits)