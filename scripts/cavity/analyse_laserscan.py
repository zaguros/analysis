"""
script to analyse laserscan data
"""
import numpy as np
import math

from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from matplotlib import pyplot as plt
from analysis.lib.fitting import fit,common


from analysis.scripts.cavity import cavity_general_analysis as cga; reload(cga)

new_y = np.zeros(44)

class laserscan_analysis(cga.cavity_analysis):

    def get_and_bin_data(self,binsize=100,**kw):
        avg_per_s=kw.pop('avg_per_s',False)
        del_points =kw.pop('del_points',0)
        self.get_laserscan_data()
        self.bin_data_per_ms(binsize=binsize)
        if avg_per_s:
            self.avg_bins_per_s(**kw)
        self.calc_sweeptime()
        return self.x_data,self.y_data,self.frq_data,self.y_data_per_ms_binned

    def poly_fit_frq(self,**kw):
        plot_fit = kw.pop('plot_fit',True)

        g_x0 = self.x_data[0]
        g_a0 = self.frq_data[0]
        g_a1 = (self.frq_data[-1]-self.frq_data[0])/(self.x_data[-1]-self.x_data[0])
        g_a2 = -0.2
        g_a3=2
        fixed=[]

        p0, fitfunc, fitfunc_str = common.fit_poly_shifted(g_x0,g_a0, g_a1, g_a2,g_a3)
        fit_result = fit.fit1d(self.x_data,self.frq_data, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True, fixed=fixed)
        if plot_fit:
            plot.plot_fit1d(fit_result, np.linspace(self.x_data[-1],self.x_data[0],10*len(self.x_data)),
                label='Fit',show_guess=False, plot_data=True,color='red', data_linestyle = '-', print_info= True)

        self.frq_data =fit_result['fitfunc'](self.x_data)
        return self.frq_data

    def plot_bins(self,**kw):
        fig,ax = plt.subplots(figsize=(10,4))
        colormap = plt.cm.copper
        plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, self.nr_bins/2)])
        for i in np.arange(self.nr_bins):
            ax.plot(self.frq_data,self.y_data_per_ms_binned[:,i], label = '%d'%(i))
        ax.set_title(self.folder+'/binsize = %d'%(self.binsize))
        ax.legend()
        fig.savefig(self.folder+'/data_per_ms_binned_binsize_%d.png'%(self.binsize))
        
        plt.show()
        plt.close()

    def find_nearest(self, array,value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]

    def shift_centre(self,**kw):
        """
        function that produces new x-data, based on the detuning from the centre frequency.
        Ouput:
        self.x_data
        """
        y=kw.pop('y',self.y_data)
        x=kw.pop('x',self.frq_data)
        FWHM,u_FWHM,x0 = self.fit_gaussian(**kw)
        self.dx_data = x-x0
        return self.dx_data

    def bin_shifted_data(self,**kw):
        """
        function that bins the dx-data in the frequency domain 
        (note: these are not the bins in time-after-pulsetube domain).
        Always use after calling 'shift_centre', to determine the detuning from centre frequency.
        Output:
        self.binsGHz  - the x-bins
        self.y_data_per_ms_binned_shifted  - the y's in the x-bins
        """
        binsizeGHz = kw.pop('binsizeGHz',1.)#in GHz
        rangeGHz = kw.pop('rangeGHz',50)
        self.binsGHz = np.linspace(-rangeGHz/2.,rangeGHz/2,rangeGHz/binsizeGHz+1)
        self.y_data_per_ms_binned_shifted  = np.zeros((len(self.binsGHz),self.nr_bins))

        for i,b in enumerate(self.binsGHz):
            bin_indices = np.where( (self.dx_data>b-binsizeGHz/2.)&(self.dx_data<b+binsizeGHz/2.) )
            for j in np.arange(self.nr_bins):
                self.y_data_per_ms_binned_shifted[i,j] = np.average(self.y_data_per_ms_binned[:,j][bin_indices])

        return self.binsGHz, self.y_data_per_ms_binned_shifted

    def crop_data_around_centre(self,i_window=10,**kw):
        """
        function that crops that data i_window datapoints around the centre. Not used at the moment.
        """
        ix_centre = int(np.where(self.dx_data==0.)[0])
        i_min = max(0,ix_centre-i_window)
        i_max = ix_centre+i_window+1 
        if i_max > len(self.dx_data):
            i_max = -1
        #use_ix = np.arange(ix_centre-i_window,ix_centre+i_window)
        #use_ix=np.where((self.dx_data>-dfreq_window/2.) & (self.dx_data<dfreq_window/2.))
        self.dx_data_centred=self.dx_data[i_min:i_max]
        self.y_data_centred = self.y_data[i_min:i_max]
        print self.dx_data_centred
        return self.dx_data_centred,self.y_data_centred

    def crop_binned_data_around_centre(self,i_window=10,**kw):
        """
        function that crops the binned data i_window datapoints around the fitted centre frequency. Not used at the moment.
        """
        ix_centre = int(np.where(self.dx_data==0.)[0])
        i_min = max(0,ix_centre-i_window)
        i_max = ix_centre+i_window+1 
        if i_max > len(self.dx_data):
            i_max = -1
        #use_ix=np.where((self.dx_data>-dfreq_window/2.) & (self.dx_data<dfreq_window/2.))
        self.dx_data_centred=self.dx_data[i_min:i_max]
        self.y_data_per_ms_binned_centred = np.zeros((len(self.dx_data_centred),self.nr_bins))
        for i in np.arange(self.nr_bins):
            self.y_data_per_ms_binned_centred[:,i] = self.y_data_per_ms_binned[i_min:i_max,i]
        return self.y_data_per_ms_binned_centred

    def fit_gaussian(self,ax=None,ret_ax =False,**kw):
        """
        function that fits a gaussian, and returns the FWHM (not sigma).
        x and y can be given as input, ad well as all the parameter guesses.
        """
        y=kw.pop('y',self.y_data)
        x=kw.pop('x',self.frq_data)
        g_sigma = kw.pop('g_sigma',10)
        plot_title = kw.pop('plot_title','')
        plot_fit=kw.pop('plot_fit',False)
        save_plot=kw.pop('save_plot',True)
        fixed=kw.pop('fixed',[])
        label=kw.pop('label','')
        g_a = np.nanmean(y)
        g_A = max(y)-g_a
        g_x0= x[np.argmax(y)]

        p0, fitfunc, fitfunc_str = common.fit_gauss(g_a, g_A, g_x0, g_sigma)
        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
        if 'sigma' not in fit_result['params_dict']:
            print 'WARNING: COULD NOT FIT sigma'
            return 0,10,0

        sigma = np.abs(fit_result['params_dict']['sigma'])
        u_sigma = np.abs(fit_result['error_dict']['sigma'])
        x0 = fit_result['params_dict']['x0']
        FWHM = 2*math.sqrt(2*math.log(2))*sigma
        u_FWHM = 2*math.sqrt(2*math.log(2))*u_sigma

        if plot_fit:
            if ax==None:
                fig,ax = plt.subplots(figsize=(10,4))
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],len(x)*10),
                ax=ax,label=label,show_guess=False, plot_data=False,print_info=False)
            ax.set_title('FWHM = %.1f +- %.1f GHz %s'%(FWHM,u_FWHM,plot_title))
            if save_plot:
                fig = ax.get_figure()
                fig.savefig(self.folder+'/gaussian_fit_%s.png'%(plot_title))

            if ret_ax:
                return FWHM,u_FWHM,ax

            plt.show(ax.get_figure())
            plt.close(ax.get_figure())
        
        return FWHM,u_FWHM,x0

    def fit_lorentzian(self,**kw):
        y=kw.pop('y',self.y_data)
        x=kw.pop('x',self.frq_data)
        g_gamma = kw.pop('g_gamma',15)
        plot_title = kw.pop('plot_title','')
        plot_fit=kw.pop('plot_fit',False)
        fixed=kw.pop('fixed',[])

        g_a = np.average(y)
        g_A = max(y)-g_a
        g_x0= self.frq_data[np.argmax(y)]


        p0, fitfunc, fitfunc_str = common.fit_lorentz(g_a, g_A, g_x0, g_gamma)
        fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
        if 'gamma' not in fit_result['params_dict']:
            fig,ax = plt.subplots(figsize=(10,4))      
            ax.plot(x,y)
            print 'WARNING: COULD NOT FIT gamma'
            return 0,10

        gamma = np.abs(fit_result['params_dict']['gamma'])
        u_gamma = np.abs(fit_result['error_dict']['gamma'])
        x0 = fit_result['params_dict']['x0']
        FWHM = gamma
        u_FWHM = u_gamma

        if plot_fit:
            fig,ax = plt.subplots(figsize=(10,4))
            plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],len(x)*10),
                ax=ax,label='Fit',show_guess=True, plot_data=True)
            ax.set_title('FWHM = %.1f +- %.1f GHz %s'%(FWHM,u_FWHM,plot_title))
            fig.savefig(self.folder+'/lorentzian_fit_%s.png'%(plot_title))
            plt.show(fig)
            plt.close(fig)

        return FWHM,u_FWHM

    def fit_bins(self,**kw):
        """
        function that fits a gaussian to all the bins in which the data is ordered.
        """
        plot_sigmas = kw.pop('plot_sigmas',True)
        binned_data=kw.pop('binned_data',self.y_data_per_ms_binned)
        g_sigmas = kw.pop('g_sigmas',np.ones(self.nr_bins)*10)
        #print g_sigmas

        FWHMs=np.array([])
        u_FWHMs=np.array([])

        for i in np.arange(self.nr_bins):
            FWHM,u_FWHM,x0 = self.fit_gaussian(x=self.frq_data, y=binned_data[:,i],plot_title ='bin_%d_binsize_%d.png'%(i,self.binsize), g_sigma=g_sigmas[i],**kw )
            
            FWHMs=np.append(FWHMs,FWHM)
            u_FWHMs=np.append(u_FWHMs,u_FWHM)

        if plot_sigmas:
            fig2,ax2 = plt.subplots(figsize=(10,4))
            ax2.errorbar(np.arange(self.nr_bins),FWHMs,yerr=u_FWHMs)
            ax2.set_xlabel('bin nr')
            ax2.set_ylabel('FWHM')
            #ax2.set_ylim([0,80])
            ax2.set_title(self.folder+'/FWHMs_vs_bins_binsize_%d.png'%self.binsize)
            fig2.savefig(self.folder+'/FWHMs_vs_bins_binsize_%d.png'%self.binsize)


        return FWHMs,u_FWHMs,x0

    def get_bin(self,bin_nr):
        return self.y_data_per_ms_binned[:,bin_nr]

    def get_linecut_per_ms(self,linecut_point=15,**kw):
        do_plot=kw.pop('do_plot',True)
        self.y_per_ms_linecut=self.y_data_per_ms[linecut_point,:]
        if do_plot:
            fig,ax = plt.subplots(figsize=(10,4))
            ax.plot(self.y_per_ms_linecut, label = '%d'%(linecut_point))
            ax.set_title(self.folder+'/PDsignal_per_mslinecut%d.png'%linecut_point)
            ax.set_ylabel('PD signal at single frequency')
            ax.set_xlabel('ms after sync')
            fig.savefig(self.folder+'/PDsignal_per_mslinecut%d.png'%linecut_point)
        return self.y_per_ms_linecut

    def plot_PD_signal_vs_V(self,ax=None,ret_ax=False,**kw):
        save_fig = kw.pop('save_fig',True)
        if ax==None:
            fig,ax = plt.subplots()
        ax.set_title(self.folder+'\ntotalsweeptime=%.1fs'%(self.total_sweep_time))
        ax.set_xlabel('laser tuning voltage (V)')
        ax.set_ylabel('PD signal (V)')
        ax.plot(self.x_data,self.y_data)
        if save_fig:
            fig = ax.get_figure()
            fig.savefig(self.folder+'/PDsignal_vs_laser_tuning_voltage.png')
        if ret_ax:
            return ax
        plt.show()
        plt.close()

    def plot_PD_signal_vs_frq(self,ax=None,ret_ax=False,**kw):
        save_fig = kw.pop('save_fig',True)
        if ax==None:
            fig,ax = plt.subplots()     
        ax.set_title(self.folder+'\ntotalsweeptime=%.4fs'%(self.total_sweep_time))
        ax.set_xlabel('frequency (GHz) + 470.4 THz')
        ax.set_ylabel('PD signal (V)')
        ax.plot(self.frq_data,self.y_data)
        if save_fig:
            fig = ax.get_figure()
            fig.savefig(self.folder+'/PDsignal_vs_frq.png')
        if ret_ax:
            return ax
        plt.show()
        plt.close()

    def plot_PD_signal_vs_dfrq(self,ax=None,ret_ax=False,**kw):
        save_fig = kw.pop('save_fig',True)
        if ax==None:
            fig,ax = plt.subplots()     
        ax.set_title(self.folder+'\ntotalsweeptime=%.4fs'%(self.total_sweep_time))
        ax.set_xlabel('delta frequency (GHz)')
        ax.set_ylabel('PD signal (V)')
        ax.plot(self.dx_data_centred,self.y_data_centred)
        if save_fig:
            fig = ax.get_figure()
            fig.savefig(self.folder+'/PDsignal_vs_dfrq.png')
        if ret_ax:
            return ax
        plt.show()
        plt.close()

    def plot_frq_vs_V(self,ax=None,ret_ax=False,**kw):
        save_fig = kw.pop('save_fig',True)
        if ax==None:
            fig,ax = plt.subplots()        
        ax.set_title(self.folder+'\ntotalsweeptime=%.1fs'%(self.total_sweep_time))
        ax.set_ylabel('frequency (GHz) + 470.4 THz')
        ax.set_xlabel('laser tuning voltage (V)')
        ax.plot(self.x_data,self.frq_data)
        if save_fig:
            fig = ax.get_figure()
            fig.savefig(self.folder+'/laser_tuning_voltage_vs_frq.png')
        if ret_ax:
            return ax
        plt.show()
        plt.close()


# def laserscan_analysis(folder,**kw):
#     fit_gaussian= kw.pop('fit_gaussian',False)
#     fit_lorentz_avg  = kw.pop('fit_lorentz_avg',False)
#     plot_frq_vs_V = kw.pop('plot_frq_vs_V',False)
#     plot_PD_signal_vs_V = kw.pop('plot_PD_signal_vs_V',False)
#     poly_fit_frq = kw.pop('poly_fit_frq',False)
#     ms_min= kw.pop('ms_min',None)
#     ms_max= kw.pop('ms_max',None)
#     do_plot=kw.pop('do_plot',False)
#     linecut_point = kw.pop('linecut_point',0)

#     a = cga.cavity_analysis(folder)
#     x,y,frq,y_per_ms = a.get_laserscan_data()
#     y_per_ms_binned = a.bin_data_per_ms(**kw)



    # print 'total sweep time = ', total_sweep_time 





    # if fit_lorentz_avg:

    #     global new_y
    #     new_y = np.add(new_y,y)

    #     avg_y = new_y/50

    #     fig,ax = plt.subplots()
    #     ax.set_xlabel('Voltage (V)')
    #     ax.set_ylabel('PD signal (V)')
    #     plt.plot(x,avg_y)
    #     fig.savefig(a.folder+'/PDsignal_vs_laser_tuning_voltage.png')
    #     plt.show()

    #     fig,ax = plt.subplots()

    #     ax.set_xlabel('Frequency (GHz) + 470.4 THz')
    #     ax.set_ylabel('PD signal (V)')
    #     plt.plot(frq,avg_y)
    #     fig.savefig(a.folder+'/PDsignal_vs_frq.png')
    #     plt.show()

    #     g_a = 0.7
    #     g_A = 3
    #     g_x0 = 4
    #     g_sigma = 4
    #     fixed = []



    #     p0, fitfunc, fitfunc_str = common.fit_lorentz(g_a, g_A, g_x0, g_sigma)
    #     fit_result = fit.fit1d(x,avg_y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, fixed=fixed)
    #     fig,ax = plt.subplots(figsize=(8,4))
    #     ax.set_xlabel('Voltage (V)')
    #     ax.set_ylabel('PD signal (V)')
    #     plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],10*len(x)),ax=ax, label='Fit',show_guess=True, plot_data=True,color='red', data_linestyle = '-', print_info= False)

    #     A = fit_result['params_dict']['A']
    #     a = fit_result['params_dict']['a']
    #     x0 = fit_result['params_dict']['x0']
    #     gamma = fit_result['params_dict']['gamma']

    #     linewidth = (frq[0]-frq[-1])/(x[-1]-x[0])*gamma

    #     ax.set_title('Linewidth = %.4fGHz'%linewidth)


###########################################


    # fig,ax = plt.subplots()
    # ax.set_title(folder+'\ntotalsweeptime=%.4fs'%(total_sweep_time))
    # ax.set_xlabel('frequency (GHz) + 470.4 THz')
    # ax.set_ylabel('PD signal (V)')
    # ax.set_ylim([0.014,0.035])
    # ax.plot(frq,y_per_ms_avg)
    # fig.savefig(a.folder+'/PDsignal_vs_frq.png')
    # plt.show()









        




