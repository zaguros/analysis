# analyse 2D data of spectrometer 
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import seaborn as sns
import glob
import os
import matplotlib.image as mpimg
import scipy 
import time

from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
reload(common)
import analysis.scripts.cavity.spectrometer_analysis as sa
reload(sa)


# parameters to vary per measurement Note: you might have to change the vmin and vmax of the colorbar inside the script! 

n_diamond = 2.419 #refractive index diamond
c = 2.99792458e8 #speed of light


class spectrometer_2D_analysis(sa.spectrometer_analysis):
    def __init__(self,folder, V_min, V_max,laser_wavelength=None):
        #initialise data-independent parameters
        self.laser_wavelength = laser_wavelength
        self.folder = folder
        self.V_min = V_min
        self.V_max = V_max
        self.plot_name = os.path.split(os.path.split(os.path.split(self.folder)[0])[0])[-1] +\
            '/'+os.path.split(os.path.split(self.folder)[0])[-1] + '/'+os.path.split(self.folder)[-1]

        #get data
        self.get_data()

        #initialise some useful data-dependent parameters
        self.V_range=np.abs(self.V_max-self.V_min)
        self.V_extent_correction = self.V_range/float((len(self.filenumbers)-1))/2.
        self.frq_range = abs(self.frequencies[-1]-self.frequencies[0])
        self.frq_extent_correction = self.frq_range/(float(len(self.frequencies)-1))/2.





    def get_data(self,min_frq=0,max_frq=1000):
        self.frequencies,self.filenumbers,self.intensities = self.load_data_from_folder()
        #the ability to select a certain frequency range        
        i_min = np.argmin(np.abs(self.frequencies-min_frq)) 
        i_max = np.argmin(np.abs(self.frequencies-max_frq))
        self.frequencies=self.frequencies[i_max:i_min]
        self.intensities=self.intensities[i_max:i_min,:]

        return self.frequencies,self.filenumbers,self.intensities 


    def subtract_offset(self):
        """
        function that subtracts the average offset (average over voltage of spectrometer data) from the data
        """
        offset = np.average(self.intensities,axis=1)
        offsets=np.array([offset,]*len(self.intensities[0])).transpose()
        self.intensities = self.intensities-offsets
        return self.intensities

    def plot_data_quickly(self,**kw):
        """
        Function that can be used for quick plotting of the data.
        HOWEVER, it is plotting the frequencies on the y-axis evenly spaced. this is WRONG.
        """
        title = kw.pop('title','/quick_2D_plot_imshow')
        cmap = kw.pop('cmap','YlGnBu')
        vmax = kw.pop('vmax',None)
        vmin = kw.pop('vmin',None)
        aspect = kw.pop('aspect','auto')
        fig,ax = plt.subplots()

        extent = [self.V_min-self.V_extent_correction,self.V_max+self.V_extent_correction,\
            self.frequencies[-1]-self.frq_extent_correction,self.frequencies[0]+self.frq_extent_correction]
        extent = [self.V_min,self.V_max,self.wavelengths[-1],self.wavelengths[0]]
        im = ax.imshow(self.intensities, extent= extent, vmax =vmax, vmin=vmin,cmap = cmap, aspect = aspect,interpolation='None')
        ax = self.set_axes_basics(ax)
        ax.set_title(self.plot_name+title)
        plt.colorbar(im)
        print 'Note that the y axis of this figure is distorted - data is NOT in reality equally spaced in frequency'

        try: 
            print 'saving figure as:'
            print os.path.join(self.folder, '2D_plot_imshow.png')
            fig.savefig(os.path.join(self.folder, '2D_plot_imshow.png'))
        except:
            print('could not save figure')

        plt.show()
        plt.close()

        return fig,ax

    def plot_data(self, **kw):
        title = kw.pop('title','/2D_plot')
        cmap = kw.pop('cmap','YlGnBu')
        vmax = kw.pop('vmax',None)
        vmin = kw.pop('vmin',None)
        aspect = kw.pop('aspect','auto')

        extent = [self.V_min-self.V_extent_correction,self.V_max+self.V_extent_correction,\
            self.frequencies[-1]-self.frq_extent_correction,self.frequencies[0]+self.frq_extent_correction]

        #add an extra point to the Vs array, to make the pcolor plot work
        Vs_xl = np.linspace(self.V_min-self.V_extent_correction,self.V_max+self.V_extent_correction,len(self.filenumbers)+1)
        #also adding an extra point to frequencies, 
        #but due to the unequal spacing in frequencies, I do it this way, which is not completely correct.
        #the added frequency should also be unequally spaced... 
        frqs_xl = np.append(self.frequencies+self.frq_extent_correction,self.frequencies[-1]-self.frq_extent_correction)
        x,y = np.meshgrid(Vs_xl,frqs_xl) 
        fig,ax = plt.subplots()
        im = ax.pcolormesh(x,y,self.intensities,vmax =vmax, vmin=vmin,cmap = cmap)
        ax = self.set_axes_basics(ax)
        ax.set_title(self.plot_name+title)
        plt.colorbar(im)

        ax.set_xlim([self.V_min-self.V_extent_correction,self.V_max+self.V_extent_correction])
        ax.set_ylim([self.frequencies[-1]-self.frq_extent_correction,self.frequencies[0]+self.frq_extent_correction])
        print ax.get_xlim()
        print ax.get_ylim()

        try: 
            print 'saving figure as:'
            print os.path.join(self.folder, '2D_plot.png')
            fig.savefig(os.path.join(self.folder, '2D_plot.png'))
        except:
            print('could not save figure')

        plt.show()
        plt.close()

        return fig,ax

    def set_axes_basics(self, ax):
        ax.set_xlabel("Voltage (V)", fontsize = 14)
        ax.set_ylabel("Frequency (THz)", fontsize = 14)

        ax.grid(False)
        ax.set_axis_bgcolor('white')

        return ax

    def set_axes_ticks(self, ax):

        ax.tick_params(which = 'both', direction = 'out')
        xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],int((self.V_max-self.V_min)/2+1))
        xticklabels = np.linspace(self.V_min,self.V_max,int((self.V_max-self.V_min)/2+1))
        xticklabels = np.round(xticklabels,1)

        yticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[-1],7)
        ytickslabels = np.linspace(self.frequencies[-1],self.frequencies[0],7)#in THz
        ytickslabels =np.round(ytickslabels,0).astype(int)

        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, rotation=0)

        ax.set_yticks(yticks)
        ax.set_yticklabels(ytickslabels)
        return ax

    def set_axes(self, ax):
        ax=self.set_axes_basics(ax)
        ax=self.set_axes_ticks(ax)
        return ax

    def plot_from_2D_data(self, **kw):
        self.plot_data(**kw)    

    def shifted_plot_from_2D_data(self, folder):     
        fig,ax = plt.subplots(figsize=(6,4))
        ax=self.set_axes(ax)

        for i,intensity in enumerate(np.transpose(self.intensities)):
            ax.plot(self.frequencies,intensity-i*15000)
        
        plt.show()

        try: 
            print os.path.join(self.folder, 'shifted_plots.png')
            fig.savefig(os.path.join(self.folder, 'shifted_plots.png'))
        except:
            print('could not save figure')

        plt.close()

    def remove_higher_order_modes(self,x0s,u_x0s,hom_max =10,**kw):
        """
        function that removes the peaks that are higher order modes. 
        If the separation dx0 beteween two peaks is wihtin the range (HOM_min <) dx0 < HOM_max,
        only the oeak with the lowest x0 is kept.
        """
        # nr_peaks = 
        report_fails = kw.pop('report_fails',False)
        success = np.zeros(len(x0s))

        indices_fund = np.where(np.abs(np.diff(x0s))>hom_max)
        x0s_fund = x0s[indices_fund]
        x0s = np.append(x0s_fund,x0s[-1])
        u_x0s_fund = u_x0s[indices_fund]
        u_x0s =  np.append(u_x0s_fund,u_x0s[-1])
        success[indices_fund]=1
        success[-1]=1
        if report_fails:
            print '%d higher order modes removed'%(len(success) - len(x0s))
        return x0s,u_x0s,success


    def peaks_from_1D_data(self, intensity,**kw):
        '''
        function that finds the peaks in 1D data with 
        x = wavelengths,
        y = intensity,
        output parameters:
        x0s - array of xs corresponding to peak locations
        u_x0s - uncertainty in xs from lorentzian fit
        '''
        plot_fit = kw.pop('plot_fit',False)
        plot_peak_locations = kw.pop('plot_peak_locations',False)
        save_fig = kw.pop('save_fig', False)
        fit_peaks = kw.pop('fit_peaks', True)
        remove_hom = kw.pop('remove_hom',False)
        hom_max = kw.pop('hom_max',10)

        indices,peak_frequencies,peak_intensity = self.approximate_peak_location(intensity,**kw)
        if fit_peaks:
            x0s,u_x0s,success = self.fit_peak(intensity,indices,peak_frequencies,peak_intensity,plot_fit = plot_fit,**kw)
        else:
            x0s=peak_frequencies
            u_x0s = np.zeros(len(peak_frequencies))
            success = np.ones(len(peak_frequencies))
        peak_intensity_x0s = peak_intensity[np.where(success >0)]

        if remove_hom:
            x0s,u_x0s,hom_success = self.remove_higher_order_modes(x0s,u_x0s,hom_max=hom_max,**kw)
            peak_intensity_x0s = peak_intensity_x0s[np.where(hom_success >0)]



        # print len(peak_intensity),len(peak_intensity_x0s)
        if plot_peak_locations:
            fig, ax = plt.subplots()
            ax.plot(self.frequencies,intensity)
            # ax.plot(peak_wavelengths,peak_intensity,'o', mfc=None, mec='r', mew=2, ms=8)
            ax.plot(x0s,peak_intensity_x0s,'+', mfc=None, mec='r', mew=2, ms=8)
            ax.set_title(self.plot_name)

            if save_fig:
                plt.savefig(os.path.join(self.folder,'plot_1D_peaks.png'))
            plt.show()
            plt.close()

        return x0s,u_x0s    

    def peaks_from_2D_data(self,**kw):
        '''
        function that finds the peaks in 2D data in folder.
        Output: the peak locations, if return_peak_locations is true.
        '''
        return_peak_locations = kw.pop('return_peak_locations',False)

        x=np.array([])
        u_x = np.array([])
        y=np.array([])

        print 'getting peak locations'
        for i,intensity in enumerate(np.transpose(self.intensities)):
            # print i,'/',len(np.transpose(self.intensities))
            # if i>0:
            #     break
            x0s,u_x0s = self.peaks_from_1D_data(intensity,**kw)
            x = np.append(x,x0s)
            u_x = np.append(u_x,u_x0s)
            for j in np.arange(len(x0s)):
                y = np.append(y,i)

        self.peak_x = x
        self.u_peak_x = u_x
        self.peak_y = y /(len(self.filenumbers)-1)*self.V_range

        if return_peak_locations:
            return self.peak_x, self.peak_y


    def plot_peaks(self,**kw):
        """
        function that plots the peaks found in the data in a scatter plot
        """
        save_fig = kw.pop('save_fig',False)
        ax = kw.pop('ax',None)

        # make a scatter plot of the peaks
        if ax == None:
            fig,ax = plt.subplots(figsize =(6,4))    

        ax.scatter(self.peak_y,self.peak_x)
        print ax.get_xlim(),ax.get_ylim()
        # ax.errorbar(self.peak_y, self.peak_x, ls='none',marker=None,yerr= self.u_peak_x)
        ax.set_title(self.plot_name+'/peaks.png')
        ax.set_xlim((self.V_min-self.V_extent_correction,self.V_max+self.V_extent_correction))
        ax.set_ylim((self.frequencies[-1]-self.frq_extent_correction, self.frequencies[0]+self.frq_extent_correction))
        print ax.get_xlim(),ax.get_ylim()
        # if self.laser_wavelength!=None:
        #     ax.plot([ax.get_xlim()[0],ax.get_xlim()[-1]],[c/self.laser_wavelength*1.e-12,c/self.laser_wavelength*1.e-12]) #laser wavelength in THz

        ax=self.set_axes(ax)

        if save_fig:
            try: 
                print 'saving figure to: ',os.path.join(self.folder, 'peaks.png')
                fig = ax.get_figure()
                fig.savefig(os.path.join(self.folder, 'peaks.png'))
            except:
                print('could not save figure')
        # plt.show()
        # plt.close()

        return ax


    def find_best_overlap_peaks_and_modes(self,diamond_thicknesses, air_lengths, conversion_factor=307.e-9, **kw):
        # x,y,fig,ax = self.peaks_from_2D_data(return_peak_locations=True,**kw)
        fig,ax = plt.subplots()
        ax = self.plot_peaks(ax=ax,**kw)

        ms_errors = np.zeros((len(diamond_thicknesses),len(air_lengths)))
        u_ms_errors = np.zeros((len(diamond_thicknesses),len(air_lengths)))

        t0=time.time() 
        for i,diamond_thickness in enumerate(diamond_thicknesses):
            for j,air_length in enumerate(air_lengths):
                cavity_length = diamond_thickness + air_length
                modes = self.diamond_air_modes(cavity_length=cavity_length,diamond_thickness=diamond_thickness,
                    conversion_factor=conversion_factor,nr_points=max(self.peak_y)+1)

                ms_error, u_ms_error = self.calculate_overlap_quality(self.peak_x,self.peak_y,modes,**kw)
                # print 15*'*'
                # print 'mean squared error fit',ms_error, '+-', u_ms_error
                # print 15*'*'
                ms_errors[i,j] = ms_error
                u_ms_errors[i,j] = u_ms_error
            if (i%5==0):
                t1=time.time()
                time_elapsed = t1-t0
                average_time = (t1-t0)/(i+1)
                print 'diamond thickness',i+1,'out of ', len(diamond_thicknesses), 'done.',\
                    'estimated time remaining: ', int((len(diamond_thicknesses)-i-1)*average_time), 's'
                

        ix_min_mean_square_overlap = np.unravel_index(ms_errors.argmin(), ms_errors.shape)
        print 'lowest mean square error (',round(ms_errors[ix_min_mean_square_overlap],3), '+-',round(u_ms_errors[ix_min_mean_square_overlap],3),') is found for:'
        print 'diamond thickness: ',diamond_thicknesses[ix_min_mean_square_overlap[0]]
        print 'air length: ',air_lengths[ix_min_mean_square_overlap[1]]
        print 'total cavity length:', diamond_thicknesses[ix_min_mean_square_overlap[0]]+air_lengths[ix_min_mean_square_overlap[1]]

        self.ms_errors = ms_errors
        self.u_ms_errors = u_ms_errors

        return ms_errors, u_ms_errors

    def overlap_peaks_and_modes(self, diamond_thickness=4.e-6,air_length = 5.e-6,
            conversion_factor = 307.e-9,nr_points=61, ax=None, **kw):
        '''
        function that plots the fitted peak locations in 2D data in folder, 
        and overlaps it with the analytically derived diamond and air modes.
        Input parameters:
        folder - directory containing 2D data
        diamond_thickness - diamond thickness used to obtain analytic result for resonance frequency
        cavity_length - cavity length used to obtain analytic result for resonance frequency
        conversion_factor - the piezo conversion factor. at RT:307 nm/V,. at LT: 123 nm/V in latest mstm
        nr_points - the number of points used for plotting analytic results of resonances
        mode_type - the type of the modes plotted = possible are 'diamond_air_modes' or 'air_modes'
        **keywords for peak-finding:
        plot_fit - whether to plot the fit of each resonance found in the data

        '''
        # x,y,fig,ax = self.peaks_from_2D_data(return_peak_locations=True,**kw)
        # ax = kw.pop('ax',None)
        if ax==None:
            fig,ax = plt.subplots()
            ax = self.plot_peaks(ax=ax,**kw)

        mode_type = kw.pop('type','diamond_air_modes') #can be 'diamond_air_modes', or 'air_modes'
        ret_ax = kw.pop('ret_ax',False)


        if mode_type == 'diamond_air_modes':
            print 'overlapping diamond air modes'
            modes,ax = self.plot_diamond_air_modes(air_length=air_length,diamond_thickness=diamond_thickness,
                ax=ax,conversion_factor=conversion_factor,nr_points=nr_points, return_modes=True)
        elif mode_type == 'air_modes':
            modes,ax = self.plot_air_modes(air_length=air_length,conversion_factor=conversion_factor,nr_points = nr_points,ax=ax)
            pass 
        # ax = self.plot_diamond_modes(diamond_thickness=diamond_thickness,ax = ax)
        # ax = self.plot_air_modes(cavity_length=cavity_length,diamond_thickness=diamond_thickness,
        #     ax=ax, conversion_factor=conversion_factor,nr_points=nr_points)

        ms_error, u_ms_error = self.calculate_overlap_quality(self.peak_x,self.peak_y,modes,**kw)
        print 15*'*'
        print 'mean squared error', round(ms_error,3), '+-', round(u_ms_error,3)
        print 15*'*'

        title ='d={}um_L={}um'.format(str(diamond_thickness*1e6),str(air_length*1.e6))

        ax.text(ax.get_xlim()[0] + (ax.get_xlim()[-1]-ax.get_xlim()[0])/4,ax.get_ylim()[0],title, size=14, backgroundcolor = 'w')

        #add an axis at the top with the cavity length 
        ax2 = ax.twiny()
        xticks = np.linspace(ax2.get_xlim()[0],ax2.get_xlim()[-1],int((self.V_max-self.V_min)/2+1))
        xticklabels2 =np.linspace(air_length*1.e6,air_length*1.e6+(conversion_factor*(self.V_max-self.V_min)*1.e6),int((self.V_max-self.V_min)/2+1))
        xticklabels2 = np.round(xticklabels2,2)

        ax2.set_xticks(xticks)
        ax2.set_xticklabels(xticklabels2,rotation=0)
        ax2.set_xlabel('air length (um)',fontsize = 14)

        try: 
            print os.path.join(self.folder, 'overlap_peaks_and_modes{}.png'.format(title))
            fig = ax.get_figure()
            fig.savefig(os.path.join(self.folder, 'overlap_peaks_and_modes{}.png'.format(title)))
        except:
            print('could not save figure')

        if ret_ax:
            return fig,ax

        plt.close(fig)

    def find_nearest(self, array,value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]

    def calculate_overlap_quality(self, x,y,modes, **kw):
        min_frequency = kw.pop('min_frequency', 400)
        max_frequency = kw.pop('max_frequency', 550)
        min_voltage = kw.pop('min_voltage',0.)

        nr_scans_to_disregard = int((min_voltage - self.V_min)/(self.V_max - self.V_min)*(max(y)+1))

        squared_errors = []
        tot_nr_errors = 0

        for i in np.arange(max(y)+1):
            if i>nr_scans_to_disregard:
                x_i = x[np.where((y>i-0.2)&(y<i+0.2))]
                nu_i = np.transpose(modes)[i]
                
                for x_ii in x_i: #important to compare to x_ii: the data.
                    if ((x_ii>min_frequency) and (x_ii < max_frequency)):
                        nearest_nu_ii = self.find_nearest(nu_i, x_ii)
                        tot_nr_errors+=1
                        squared_errors.append((nearest_nu_ii-x_ii)**2 )

        squared_errors = np.array(squared_errors)
        total_squared_errors = np.sum(squared_errors)
        mean_squared_error = total_squared_errors/tot_nr_errors
        u_mean_squared_error = np.sqrt(np.sum((squared_errors-mean_squared_error)**2))/tot_nr_errors

        return mean_squared_error, u_mean_squared_error

    def pure_diamond_modes(self, diamond_thickness=4.e-6):
        max_nr_modes = 100
        nu_diamond = np.zeros(max_nr_modes)
        for N in np.arange(max_nr_modes):
            nu_diamond[N] = (N * c / (2 * n_diamond *diamond_thickness))/1.e12 # in THz

        return nu_diamond

    def plot_diamond_modes(self, diamond_thickness=4.e-6,ax = None):
        return_fig = False
        if ax == None:
            return_fig = True
            fig,ax = plt.subplots()

        nu_diamond = self.pure_diamond_modes(diamond_thickness)

        for N,nu in enumerate(nu_diamond):
            ax.plot(ax.get_xlim(),[nu,nu], lw=2)
            # ax.text(ax.get_xlim()[-1],nu, 'N={}'.format(N))

        if return_fig:
            return fig,ax

        return ax

    def pure_air_modes(self, cavity_length=1.e-6,conversion_factor = -150.e-9,nr_points=31):
        delta_V = self.V_max - self.V_min
        delta_L = delta_V*(conversion_factor) # in m
        Ls = np.linspace(cavity_length,cavity_length+delta_L,nr_points)

        max_nrmgi_odes = 180
        nu_air = np.zeros((max_nr_modes,nr_points))
        for N in np.arange(max_nr_modes):
            for i,L in enumerate(Ls):
                nu_air[N,i] = (N * c / (2 * L))/1.e12 # in THz

        return nu_air

    def plot_air_modes(self, cavity_length=1.e-6,diamond_thickness=0.e-6,ax = None,conversion_factor = -150.e-9,nr_points=31):
        return_fig = False
        if ax == None:
            return_fig = True
            fig,ax = plt.subplots()

        nu_air = self.pure_air_modes(cavity_length=cavity_length-diamond_thickness,conversion_factor=conversion_factor,nr_points=nr_points)
        xs = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],nr_points)

        for N,nu in enumerate(nu_air):
            ax.plot(xs,nu, lw=2)
            # ax.text(ax.get_xlim()[0],nu[0], 'N={}'.format(N))

        if return_fig:
            return fig,ax

        return ax

    def diamond_air_mode_freq(self, N=1,cavity_length=1.e-6, diamond_thickness=4.e-6):
        Ltot = cavity_length+(n_diamond-1)*diamond_thickness
        Lred = cavity_length-(n_diamond+1)*diamond_thickness
        nu = c / (2*math.pi*Ltot) * \
            (math.pi*N - (-1)**N * math.asin( (n_diamond-1)/(n_diamond+1) * \
            math.sin( (N*math.pi*Lred/Ltot))))
        return nu

    def diamond_air_modes(self, cavity_length = 1.e-6, diamond_thickness = 4.e-6, conversion_factor = 100e-9,nr_points=31):
        delta_V = self.V_max - self.V_min
        delta_L = delta_V*(conversion_factor) # in m

        Ls = np.linspace(cavity_length,cavity_length+delta_L,nr_points)

        max_nr_modes = 180
        nu_diamond_air = np.zeros((max_nr_modes,nr_points))

        for N in np.arange(max_nr_modes):
            for i,L in enumerate(Ls):
                nu_diamond_air[N,i] = self.diamond_air_mode_freq(N=N,cavity_length=L,diamond_thickness=diamond_thickness)/1.e12 # in THz

        return nu_diamond_air

    def plot_diamond_air_modes(self, air_length=1.e-6,diamond_thickness=4.e-6,ax = None,conversion_factor = -150.e-9,nr_points=31, return_modes=False):
        return_fig = False
        if ax == None:
            return_fig = True
            fig,ax = plt.subplots()

        cavity_length= air_length+diamond_thickness
        nu_diamond_air = self.diamond_air_modes(cavity_length=cavity_length,diamond_thickness=diamond_thickness,conversion_factor=conversion_factor,nr_points=nr_points)
        xs = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],nr_points)


        for N,nu in enumerate(nu_diamond_air):
            ax.plot(xs,nu, lw=2)
            if (nu[0]<ax.get_ylim()[-1]) and (nu[0]>ax.get_ylim()[0]):
                ax.text(ax.get_xlim()[0],nu[0], 'N={}'.format(N))

        if return_fig:
            return fig,ax

        if return_modes:
            return nu_diamond_air,ax
        return ax

