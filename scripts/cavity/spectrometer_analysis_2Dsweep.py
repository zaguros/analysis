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
import json
import h5py

from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
reload(common)
from analysis.scripts.cavity.spectrometer_analysis import spectrometer_analysis


# parameters to vary per measurement Note: you might have to change the vmin and vmax of the colorbar inside the script! 

n_diamond = 2.419 #refractive index diamond
c = 2.99792458e8 #speed of light


class spectrometer_2D_analysis(spectrometer_analysis):
    def __init__(self,folder,load_analysis_from_file=False):
        self.folder = folder
        
        self.plot_name = os.path.split(os.path.split(os.path.split(self.folder)[0])[0])[-1] +\
            '/'+os.path.split(os.path.split(self.folder)[0])[-1] + '/'+os.path.split(self.folder)[-1]

        self.ana_pars = {}

        if load_analysis_from_file:
            print 'loading analysis from analysis.hdf5 file'
            self.load_analysis()
        else:
            print 'setting default analysis params'
            self.set_default_params()

    def set_default_params(self):
        spectrometer_analysis.set_default_params(self)

        # initialise some default settings
        self.ana_pars['vmin'] = None
        self.ana_pars['vmax'] = 1000

        self.ana_pars['V_min']=0.
        self.ana_pars['V_max']=10
        self.ana_pars['min_frq']=432
        self.ana_pars['max_frq']=496
        self.ana_pars['laser_wavelength']= None  

        self.ana_pars['remove_hom']=True  #remove the higher order modes 
        self.ana_pars['hom_max']=2.       #maximum distance between higher order modes for one to be removed
        self.ana_pars['fit_peaks']=True 

        self.ana_pars['conversion_factor']=123e-9#80e-9 -> quadratic combination at LT# 243e-9 at RT   # #123e-9#307e-9#RT #123e-9#LT
        self.ana_pars['quadratic_conversion_factor']=0.e-9#3.23e-9 at LT # 24.5e-9/2 at RT
        self.ana_pars['cubic_conversion_factor']=0.e-9 #-2.e-9/3 at RT
        self.ana_pars['mode_type']='diamond_air_modes'

        self.ana_pars['min_fit_voltage']=0 
        self.ana_pars['min_fit_frequency']=400
        self.ana_pars['max_fit_frequency']=550
        self.ana_pars['max_error_distance'] = 100. #this was not a great success

        self.ana_pars['diamond_thickness'] = None
        self.ana_pars['air_length'] = None

        self.sweep_air_lengths = None
        self.sweep_diamond_thicknesses = None

        ################################getting the data###########################################
    def get_data(self):
        """
        function that gets the data and calculates some derived parameters
        """
        self.frequencies,self.filenumbers,self.intensities = self.load_data_from_folder()

        #the ability to select a certain frequency range
        i_min = np.argmin(np.abs(self.frequencies-self.ana_pars['min_frq'])) 
        i_max = np.argmin(np.abs(self.frequencies-self.ana_pars['max_frq']))
        self.frequencies=self.frequencies[i_max:i_min]
        self.intensities=self.intensities[i_max:i_min,:]

        #initialise some useful parameters
        self.nr_files = len(self.filenumbers)

        self.Vs = np.linspace(self.ana_pars['V_min'],self.ana_pars['V_max'],self.nr_files)
        self.V_range=np.ptp(self.Vs)#np.abs(self.ana_pars['V_max']-self.ana_pars['V_min'])

        self.V_extent_correction = self.V_range/float((self.nr_files-1))/2.
        self.frq_range = abs(self.frequencies[-1]-self.frequencies[0])
        self.frq_extent_correction = self.frq_range/(float(len(self.frequencies)-1))/2.

        return self.frequencies,self.filenumbers,self.intensities 

    def subtract_offset(self):
        """
        function that subtracts the average offset (average over voltage of spectrometer data) from the data
        """
        offset = np.average(self.intensities,axis=1)
        offsets=np.array([offset,]*len(self.intensities[0])).transpose()
        self.intensities = self.intensities-offsets
        return self.intensities




        ############################################getting the peaks fromt eh data###########################

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
        fit_peaks = kw.pop('fit_peaks', self.ana_pars['fit_peaks'])
        remove_hom = kw.pop('remove_hom',self.ana_pars['remove_hom'])
        hom_max = kw.pop('hom_max',self.ana_pars['hom_max'])

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
            #plt.close()

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

        self.peak_frq = x
        self.u_peak_frq = u_x
        self.peak_V = y /(self.nr_files-1)*self.V_range+self.Vs[0]
        self.peak_filenr = y

        if return_peak_locations:
            return self.peak_frq, self.peak_V

        ############################### finding best overlap #######################################
    def find_best_overlap_peaks_and_modes(self,ax=None,**kw):
        mode_type = kw.pop('mode_type',self.ana_pars['mode_type']) #can be 'diamond_air_modes', or 'air_modes'
        diamond_thicknesses = kw.pop('sweep_diamond_thicknesses',self.sweep_diamond_thicknesses)
        air_lengths = kw.pop('sweep_air_lengths', self.sweep_air_lengths)


        ms_errors = np.zeros((len(diamond_thicknesses),len(air_lengths)))
        u_ms_errors = np.zeros((len(diamond_thicknesses),len(air_lengths)))

        t0=time.time() 
        for i,diamond_thickness in enumerate(diamond_thicknesses):
            for j,air_length in enumerate(air_lengths):
                cavity_length = diamond_thickness + air_length
                if mode_type == 'diamond_air_modes':
                    modes = self.diamond_air_modes(air_length = air_length,diamond_thickness=diamond_thickness,
                        nr_points=self.nr_files) 
                        #important to take same nr_points for the modes as filenumbers
                elif mode_type == 'air_modes':
                    modes = self.pure_air_modes(cavity_length=cavity_length,
                        nr_points=self.nr_files)
                else: 
                    print 'enter valid mode_type!'
                    return 0,0
                ms_error, u_ms_error = self.calculate_overlap_quality(modes,**kw)
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

        self.ana_pars['min_ms_error'] = ms_errors[ix_min_mean_square_overlap]
        self.ana_pars['u_min_ms_error'] = u_ms_errors[ix_min_mean_square_overlap]
        self.ana_pars['diamond_thickness'] = diamond_thicknesses[ix_min_mean_square_overlap[0]]
        self.ana_pars['air_length']=air_lengths[ix_min_mean_square_overlap[1]]
        self.ana_pars['total_cavity_length']=self.ana_pars['diamond_thickness']+self.ana_pars['air_length']
        print 'lowest mean square error (',round(self.ana_pars['min_ms_error'],3), '+-',round(self.ana_pars['u_min_ms_error'],3),') is found for:'
        print 'diamond thickness: ',self.ana_pars['diamond_thickness']
        print 'air length: ',self.ana_pars['air_length']
        # print 'total cavity length:', diamond_thicknesses[ix_min_mean_square_overlap[0]]+air_lengths[ix_min_mean_square_overlap[1]]

        self.ms_errors = ms_errors
        self.u_ms_errors = u_ms_errors

        return ms_errors, u_ms_errors


    def find_nearest(self, array,value):
        idx = (np.abs(array-value)).argmin()
        return idx, array[idx]

    def calculate_overlap_quality(self, modes, **kw):
        min_frequency = kw.pop('min_frequency', self.ana_pars['min_fit_frequency'])
        max_frequency = kw.pop('max_frequency', self.ana_pars['max_fit_frequency'])
        min_voltage = kw.pop('min_voltage', self.ana_pars['min_fit_voltage'])
        max_error_distance= kw.pop('max_error_distance', self.ana_pars['max_error_distance'])

        nr_scans_to_disregard = int((min_voltage - self.Vs[0])/(self.V_range)*(self.nr_files))
        nr_scans_to_disregard = int((min_voltage - self.Vs[0])/float((self.V_range))*(self.nr_files))
        squared_errors = []
        tot_nr_errors = 0
        for i in self.filenumbers:
            if i>nr_scans_to_disregard:
                #select the peaks in this file. use self.peak_filenr since it lists the filenumbers instead of voltages as y
                x_i = self.peak_frq[np.where(self.peak_filenr==i)]#self.peak_frq[np.where((self.peak_V>i-0.2)&(self.peak_V<i+0.2))]
                nu_i = np.transpose(modes)[i]#so it is important that modes has the same number of points as filenumbers
                for x_ii in x_i: #important to compare to x_ii: the data.
                    if ((x_ii>min_frequency) and (x_ii < max_frequency)):
                        _tmp, nearest_nu_ii = self.find_nearest(nu_i, x_ii)
                        if np.abs(nearest_nu_ii - x_ii) < max_error_distance:
                            tot_nr_errors+=1
                            squared_errors.append((nearest_nu_ii-x_ii)**2 )
                        #else:
                            #print 'ignoring since max_err _dist %.1f > %.1f'%(np.abs(nearest_nu_ii - x_ii) , max_error_distance)
        squared_errors = np.array(squared_errors)
        total_squared_errors = np.sum(squared_errors)
        mean_squared_error = total_squared_errors/tot_nr_errors
        u_mean_squared_error = np.sqrt(np.sum((squared_errors-mean_squared_error)**2))/tot_nr_errors
        return mean_squared_error, u_mean_squared_error


    ##################################plotting functions##################################

    def set_axes_basics(self, ax):
        ax.set_xlabel("Voltage (V)", fontsize = 14)
        ax.set_ylabel("Frequency (THz)", fontsize = 14)

        ax.grid(False)
        ax.set_axis_bgcolor('white')

        return ax

    def set_axes_ticks(self, ax):

        ax.tick_params(which = 'both', direction = 'out')
        xticks = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],int(self.V_range/2+1))
        xticklabels = np.linspace(self.Vs[0],self.Vs[-1],int(self.V_range/2+1))
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


    def plot_data_quickly(self,ret_ax =False,ax=None,**kw):
        """
        Function that can be used for quick plotting of the data.
        HOWEVER, it is plotting the frequencies on the y-axis evenly spaced. this is WRONG.
        """
        title = kw.pop('title','/quick_2D_plot_imshow')
        cmap = kw.pop('cmap','YlGnBu')
        vmax = kw.pop('vmax',ana_pars['vmax'])
        vmin = kw.pop('vmin',ana_pars['vmin'])
        aspect = kw.pop('aspect','auto')
        if ax==None:
            fig,ax = plt.subplots()


        extent = [self.Vs[0]-self.V_extent_correction,self.Vs[-1]+self.V_extent_correction,\
            self.frequencies[-1]-self.frq_extent_correction,self.frequencies[0]+self.frq_extent_correction]
        #extent = [self.Vs[0],self.Vs[-1],self.wavelengths[-1],self.wavelengths[0]]
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

        if ret_ax:
            return ax

        plt.show()
        plt.close()


    def plot_data(self, ret_ax=False,ax=None,**kw):
        title = kw.pop('title','/2D_plot')
        cmap = kw.pop('cmap','YlGnBu')
        vmax = kw.pop('vmax',self.ana_pars['vmax'])
        vmin = kw.pop('vmin',self.ana_pars['vmin'])
        aspect = kw.pop('aspect','auto')
        save_fig = kw.pop('save_fig',False)
        if ax==None:
            fig,ax = plt.subplots()


        extent = [self.Vs[0]-self.V_extent_correction,self.Vs[-1]+self.V_extent_correction,\
            self.frequencies[-1]-self.frq_extent_correction,self.frequencies[0]+self.frq_extent_correction]

        #add an extra point to the Vs array, to make the pcolor plot work
        Vs_xl = np.linspace(self.Vs[0]-self.V_extent_correction,self.Vs[-1]+self.V_extent_correction,len(self.filenumbers)+1)
        #also adding an extra point to frequencies, 
        #but due to the unequal spacing in frequencies, I do it this way, which is not completely correct.
        #the added frequency should also be unequally spaced... 
        frqs_xl = np.append(self.frequencies+self.frq_extent_correction,self.frequencies[-1]-self.frq_extent_correction)
        x,y = np.meshgrid(Vs_xl,frqs_xl) 
        im = ax.pcolormesh(x,y,self.intensities,vmax =vmax, vmin=vmin,cmap = cmap)
        ax = self.set_axes_basics(ax)
        ax.set_title(self.plot_name+title)
        plt.colorbar(im)

        ax.set_xlim([self.Vs[0]-self.V_extent_correction,self.Vs[-1]+self.V_extent_correction])
        ax.set_ylim([self.frequencies[-1]-self.frq_extent_correction,self.frequencies[0]+self.frq_extent_correction])

        if save_fig:
            try: 
                print 'saving figure as:'
                print os.path.join(self.folder, '2D_plot.png')
                fig.savefig(os.path.join(self.folder, '2D_plot.png'))
            except:
                print('could not save figure')

        if ret_ax:
            return ax

        plt.show()
        plt.close()


    def plot_peaks(self,ret_ax=False,**kw):
        """
        function that plots the peaks found in the data in a scatter plot
        """
        save_fig = kw.pop('save_fig',False)
        ax = kw.pop('ax',None)

        # make a scatter plot of the peaks
        if ax == None:
            fig,ax = plt.subplots(figsize =(6,4))    

        ax.scatter(self.peak_V,self.peak_frq)
        # ax.errorbar(self.peak_V, self.peak_frq, ls='none',marker=None,yerr= self.u_peak_frq)
        ax.set_title(self.plot_name+'/peaks.png')
        ax.set_xlim((self.Vs[0]-self.V_extent_correction,self.Vs[-1]+self.V_extent_correction))
        ax.set_ylim((self.frequencies[-1]-self.frq_extent_correction, self.frequencies[0]+self.frq_extent_correction))
        if self.ana_pars['laser_wavelength']!=None:
            ax.plot([ax.get_xlim()[0],ax.get_xlim()[-1]],[c/self.ana_pars['laser_wavelength']*1.e-12,c/self.ana_pars['laser_wavelength']*1.e-12]) #laser wavelength in THz
        ax=self.set_axes(ax)
        ax.grid(False)
        if save_fig:
            try: 
                print 'saving figure to: ',os.path.join(self.folder, 'peaks.png')
                fig = ax.get_figure()
                fig.savefig(os.path.join(self.folder, 'peaks.png'))
            except:
                print('could not save figure')
        
        if ret_ax:
            return ax

        plt.show()
        plt.close(fig)

    def plot_modes(self,diamond_thickness,air_length,ax=None,ret_ax =False,**kw):
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
        '''
        name = kw.pop('name','')
        mode_type = kw.pop('mode_type',self.ana_pars['mode_type']) #can be 'diamond_air_modes', or 'air_modes'
        save_fig = kw.pop('save_fig',True)

        if ax==None:
            fig,ax = plt.subplots()
            ax = self.set_axes()
        print mode_type
        if mode_type == 'diamond_air_modes':
            print 'overlapping diamond air modes'
            modes,ax = self.plot_diamond_air_modes(air_length=air_length,diamond_thickness=diamond_thickness,
                ax=ax,nr_points=self.nr_files, return_modes=True)
        elif mode_type == 'air_modes':
            modes,ax = self.plot_air_modes(air_length=air_length,nr_points = self.nr_files,ax=ax,return_modes=True)
            

        #ms_error, u_ms_error = self.calculate_overlap_quality(self.peak_frq,self.peak_V,modes,**kw)
        #print 15*'*'
        #print 'mean squared error', round(ms_error,3), '+-', round(u_ms_error,3)
        #print 15*'*'


        title ='d={}um_L={}um'.format(str(diamond_thickness*1e6),str(air_length*1.e6))

        ax.text(ax.get_xlim()[0] + (ax.get_xlim()[-1]-ax.get_xlim()[0])/4,ax.get_ylim()[0],title, size=14, backgroundcolor = 'w')

        #add an axis at the top with the cavity length 

        ax2 = ax.twiny()
        dLs = self.Vs_to_dLs_conversion()
        nr_xticks = int((self.V_range)/2+1)

        dL_at_ticks = dLs[::(len(dLs)/nr_xticks)]
        xticks = np.linspace(ax2.get_xlim()[0],ax2.get_xlim()[-1],len(dL_at_ticks))
        xticklabels2 = (air_length+dL_at_ticks)*1.e6
        xticklabels2 = np.round(xticklabels2,2)

        ax2.set_xticks(xticks)
        ax2.set_xticklabels(xticklabels2,rotation=0)
        ax2.set_xlabel('air length (um)',fontsize = 14)
        ax2.grid(False)

        if save_fig:
            try: 
                print os.path.join(self.folder, '%s_%s%s.png'%(mode_type,title,name))
                fig = ax.get_figure()
                fig.savefig(os.path.join(self.folder, '%s_%s%s.png'%(mode_type,title,name)))
            except:
                print('could not save figure')

        if ret_ax:
            return ax
        fig = ax.get_figure()
        plt.show(fig)
        plt.close(fig)


    def get_slope(self,  diamond_thickness=4.e-6,air_length = 5.e-6, nr_points=61, N=39):

         modes= self.diamond_air_modes(air_length=air_length,diamond_thickness=diamond_thickness,nr_points=nr_points)
         modeN=modes[N]

         y_diff = np.diff(modeN)
         index, frequency = self.find_nearest(modeN, c/self.laser_wavelength*1e-12)
         slope_value=1./2*(y_diff[index]+y_diff[index+1])

         final_slope = slope_value/(self.V_range/nr_points)


         return final_slope


    def plot_peaks_and_modes(self, diamond_thickness=4.e-6,air_length = 5.e-6,
            nr_points=61, ax=None,ret_ax=False, **kw):
        '''
        function that plots the fitted peak locations in 2D data in folder, 
        and overlaps it with the analytically derived diamond and air modes.
        Input parameters:
        diamond_thickness - diamond thickness used to obtain analytic result for resonance frequency
        cavity_length - cavity length used to obtain analytic result for resonance frequency
        conversion_factor - the piezo conversion factor. at RT:307 nm/V,. at LT: 123 nm/V in latest mstm
        nr_points - the number of points used for plotting analytic results of resonances
        mode_type - the type of the modes plotted = possible are 'diamond_air_modes' or 'air_modes'. default: diamond_air_modes
        '''
        diamond_thickness = kw.pop('diamond_thickness', self.ana_pars['diamond_thickness'])
        air_length = kw.pop('air_length',self.ana_pars['air_length'])

        fig, ax = plt.subplots()
        ax = self.plot_peaks(ax=ax,ret_ax = True,**kw)
        self.plot_modes(diamond_thickness,air_length,nr_points=self.nr_files,ax=ax,name='_overlap_peaks',**kw)

        if ret_ax:
            return ax
        plt.show(fig)
        plt.close(fig)


    def plot_2dplot_and_modes(self,  ax=None,ret_ax=False, **kw):
        '''
        function that plots the fitted peak locations in 2D data in folder, 
        and overlaps it with the analytically derived diamond and air modes.
        Input parameters:
        diamond_thickness - diamond thickness used to obtain analytic result for resonance frequency
        cavity_length - cavity length used to obtain analytic result for resonance frequency
        conversion_factor - the piezo conversion factor. at RT:307 nm/V,. at LT: 123 nm/V in latest mstm
        nr_points - the number of points used for plotting analytic results of resonances
        mode_type - the type of the modes plotted = possible are 'diamond_air_modes' or 'air_modes'. default: diamond_air_modes
        '''
        diamond_thickness = kw.pop('diamond_thickness', self.ana_pars['diamond_thickness'])
        air_length = kw.pop('air_length',self.ana_pars['air_length'])
        
        fig, ax = plt.subplots()
        ax = self.plot_data(ax=ax,ret_ax = True,**kw)
        self.plot_modes(diamond_thickness,air_length,nr_points=self.nr_files,ax=ax,name='_overlap_2dplot',**kw)

        if ret_ax:
            return ax
        plt.show(fig)
        plt.close(fig)

    def plot_2dplot_peaks_and_modes(self, ax=None,ret_ax=False, **kw):
        '''
        function that plots the fitted peak locations in 2D data in folder, 
        and overlaps it with the analytically derived diamond and air modes.
        Input parameters:
        diamond_thickness - diamond thickness used to obtain analytic result for resonance frequency
        cavity_length - cavity length used to obtain analytic result for resonance frequency
        conversion_factor - the piezo conversion factor. at RT:307 nm/V,. at LT: 123 nm/V in latest mstm
        nr_points - the number of points used for plotting analytic results of resonances
        mode_type - the type of the modes plotted = possible are 'diamond_air_modes' or 'air_modes'. default: diamond_air_modes
        '''
        diamond_thickness = kw.pop('diamond_thickness', self.ana_pars['diamond_thickness'])
        air_length = kw.pop('air_length',self.ana_pars['air_length'])
        
        fig, ax = plt.subplots()
        ax = self.plot_data(ax=ax,ret_ax = True,**kw)
        ax = self.plot_peaks(ax=ax,ret_ax = True,**kw)
        ax = self.plot_modes(diamond_thickness,air_length,nr_points=self.nr_files,ax=ax,name='_overlap_2dplot_peaks',**kw)

        if ret_ax:
            return ax

        plt.show(fig)
        plt.close(fig)


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

    def plot_V_intersects(self,**kw):
        ret_ax = kw.pop('ret_ax',False)
        plot_intersects = kw.pop('plot_intersects',True)
        fig, ax = plt.subplots()
        ax = self.plot_peaks(ax=ax,ret_ax = True,**kw)
        ax.plot(self.Vs_at_frq,self.frq*np.ones(len(self.Vs_at_frq)),'or')

        if ret_ax:
            return ax

        plt.savefig(os.path.join(self.folder,'V_intersects.png' ))

        if plot_intersects:
            plt.show(fig)
        plt.close(fig)


        ###################################getting and plotting all the modes ##############################

    def find_nearest_peak_frq_at_V(self,frq,V):
        dV = self.V_range/float(len(self.Vs)-1)
        V_idxs = np.where((self.peak_V<(V+dV/2.))&(self.peak_V>(V-dV/2.))) # comparing self.peak_V == V gives empty array
        _tmp,frq_nearest = self.find_nearest(self.peak_frq[V_idxs], frq)
        return frq_nearest

    def get_V_intersects(self,**kw):
        dfrq=kw.pop('dfrq',10.)
        if self.ana_pars['laser_wavelength'] ==None:
            self.frq = 470.#THz
        else:
            self.frq=c/self.ana_pars['laser_wavelength']*1.e-12

        idxs_near_frq = np.where( np.abs(self.peak_frq - self.frq) < dfrq)[0]
        Vs_near_frq = self.peak_V[idxs_near_frq]

        #remove others if there were multiple in the vicinity. 
        #could improve by checking which one is nearer
        del_is = np.array([])
        for i in np.where(np.diff(Vs_near_frq)<0.21)[0]: 
            del_is = np.append(del_is,i)
        idxs_near_frq = np.delete(idxs_near_frq,del_is)

        self.Vs_at_frq=np.array([])
        dV = self.V_range/float(len(self.Vs)-1)

        for i in idxs_near_frq:
            frq_nearest = self.peak_frq[i]
            V_nearest = self.peak_V[i]
            if frq_nearest > self.frq:
                V_first = V_nearest
                frq_first = frq_nearest
                V_last = V_nearest+dV
                print V_last
                if V_last>10.:
                    continue 
                frq_last = self.find_nearest_peak_frq_at_V(self.frq,V_last)
            elif frq_nearest < self.frq:
                V_first = V_nearest-dV
                print V_first
                if V_first<0.:
                    continue
                frq_first = self.find_nearest_peak_frq_at_V(self.frq,V_first)
                V_last = V_nearest
                frq_last = frq_nearest
           
            rel_distance_from_first = (self.frq-frq_first)/(frq_last-frq_first) #dividing two negative numbers
            V_at_frq = V_first + dV*rel_distance_from_first
            self.Vs_at_frq = np.append(self.Vs_at_frq,V_at_frq)

        print self.Vs_at_frq


    def get_conversion_factors(self,**kw):
        plot_fit = kw.pop('plot_fit',True)

        self.get_V_intersects(**kw)
        self.plot_V_intersects(**kw)

        dpeaks = np.diff(self.Vs_at_frq)
        centrepeaks = self.Vs_at_frq[:-1]+dpeaks/2
        local_conversionfactors = (c/(self.frq*1.e12))/2/dpeaks*1.e9 # in nm/V

        g_a0= 240
        g_a1 = 24
        g_a2 = -2
        fixed=[]

        p0, fitfunc, fitfunc_str = common.fit_poly(g_a0,g_a1,g_a2)#,0.1,1)
        fit_result = fit.fit1d(centrepeaks,local_conversionfactors, None, p0=p0, 
            fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed,show_guess=True)


        fig,ax= plt.subplots()
        ax.plot(centrepeaks,local_conversionfactors,'o')
        ax.set_xlabel('V')
        ax.set_ylabel('local conversion factor (dL/dV)')
        plot.plot_fit1d(fit_result, np.linspace(-2,10,101),ax=ax)
        plt.savefig(os.path.join(self.folder,'conversion_factor.png' ))

        if plot_fit:
            plt.show(fig)
        plt.close(fig)

        self.ana_pars['conversion_factor'] = fit_result['params_dict']['a0']*1.e-9
        self.ana_pars['quadratic_conversion_factor'] = fit_result['params_dict']['a1']*1.e-9/2.
        self.ana_pars['cubic_conversion_factor'] = fit_result['params_dict']['a2']*1.e-9/3.


    def Vs_to_dLs_conversion(self,**kw):
        conversion_factor = kw.pop('conversion_factor',self.ana_pars['conversion_factor'])
        quadratic_conversion_factor = kw.pop('quadratic_conversion_factor',self.ana_pars['quadratic_conversion_factor'])
        cubic_conversion_factor = kw.pop('cubic_conversion_factor',self.ana_pars['cubic_conversion_factor'])
        dLs = conversion_factor*self.Vs + quadratic_conversion_factor*(self.Vs**2)+ cubic_conversion_factor*(self.Vs**3)
        return dLs 


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


    def pure_air_modes(self, cavity_length=1.e-6,nr_points=31):
        # delta_V = self.V_max - self.V_min
        # delta_L = delta_V*(conversion_factor) # in m
        dLs = self.Vs_to_dLs_conversion()
        Ls = cavity_length + dLs 

        max_nr_modes = 180
        nu_air = np.zeros((max_nr_modes,nr_points))
        for N in np.arange(max_nr_modes):
            for i,L in enumerate(Ls):
                nu_air[N,i] = (N * c / (2 * L))/1.e12 # in THz

        return nu_air

    def plot_air_modes(self, air_length=1.e-6,ax = None,nr_points=31,return_modes=False):
        return_fig = False
        if ax == None:
            return_fig = True
            fig,ax = plt.subplots()

        nu_air = self.pure_air_modes(cavity_length=air_length,nr_points=nr_points)
        xs = np.linspace(ax.get_xlim()[0],ax.get_xlim()[-1],nr_points)

        for N,nu in enumerate(nu_air):
            ax.plot(xs,nu, lw=2)
            # ax.text(ax.get_xlim()[0],nu[0], 'N={}'.format(N))

        if return_fig:
            return fig,ax
        if return_modes:
            return nu_air,ax
        return ax

    def diamond_air_mode_freq(self, N=1,cavity_length=1.e-6, diamond_thickness=4.e-6):
        
        Ltot = cavity_length+(n_diamond-1)*diamond_thickness
        Lred = cavity_length-(n_diamond+1)*diamond_thickness
        nu = c / (2*math.pi*Ltot) * \
            (math.pi*N - (-1)**N * math.asin( (n_diamond-1)/(n_diamond+1) * \
            math.sin( (N*math.pi*Lred/Ltot))))
        return nu

    def diamond_air_modes(self, cavity_length = 1.e-6, diamond_thickness = 4.e-6, nr_points=31):
        # delta_V = self.V_max - self.V_min
        # delta_L = delta_V*(conversion_factor) # in m

        # Ls = np.linspace(cavity_length,cavity_length+delta_L,nr_points)
        dLs = self.Vs_to_dLs_conversion()
        Ls = cavity_length + dLs 

        max_nr_modes = 180
        nu_diamond_air = np.zeros((max_nr_modes,nr_points))
        cavity_length= air_length+diamond_thickness

        for N in np.arange(max_nr_modes):
            for i,L in enumerate(Ls):
                nu_diamond_air[N,i] = self.diamond_air_mode_freq(N=N,cavity_length=L,diamond_thickness=diamond_thickness)/1.e12 # in THz

        return nu_diamond_air

    def plot_diamond_air_modes(self, air_length=1.e-6,diamond_thickness=4.e-6,ax = None,nr_points=31, return_modes=False):
        return_fig = False
        if ax == None:
            return_fig = True
            fig,ax = plt.subplots()

        nu_diamond_air = self.diamond_air_modes(air_length=air_length,diamond_thickness=diamond_thickness,nr_points=nr_points)
        xs = np.linspace(self.Vs[0],self.Vs[-1],nr_points)

        for N,nu in enumerate(nu_diamond_air):
            ax.plot(xs,nu, lw=2)
            if (nu[0]<ax.get_ylim()[-1]) and (nu[0]>ax.get_ylim()[0]):
                ax.text(ax.get_xlim()[0],nu[0], 'N={}'.format(N))

        if return_fig:
            return fig,ax

        if return_modes:
            return nu_diamond_air,ax
        return ax

    def save_analysis(self, fname='analysis.hdf5'):
        print 'saving analysis'
        if not os.path.exists(os.path.join(self.folder, fname)):
            print 'creating analysis file'
            mode = 'w'    
        else:
            print 'using existng analysis file'
            mode = 'r+'            
            
        f = h5py.File(os.path.join(self.folder, fname), mode)
        g = f.require_group('2d_plot_result')

        ########################saving data arrays
        # if the group already exists, we have to delete it first
        if f.__contains__('/2d_plot_result/peak_frq'):
            del f['/2d_plot_result/peak_frq']
        if f.__contains__('/2d_plot_result/peak_V'):
            del f['/2d_plot_result/peak_V']
        if f.__contains__('/2d_plot_result/peak_filenr'):
            del f['/2d_plot_result/peak_filenr']

        try:
            f['/2d_plot_result/peak_frq'] = self.peak_frq
            f['/2d_plot_result/peak_V'] = self.peak_V
            f['/2d_plot_result/peak_filenr'] = self.peak_filenr
        except:
            print 'could not save peak locations'

        if f.__contains__('/2d_plot_result/sweep_air_lengths'):
            del f['/2d_plot_result/sweep_air_lengths']        
        if f.__contains__('/2d_plot_result/sweep_diamond_thicknesses'):
            del f['/2d_plot_result/sweep_diamond_thicknesses']

        try: 
            f['/2d_plot_result/sweep_air_lengths'] = self.sweep_air_lengths
            f['/2d_plot_result/sweep_diamond_thicknesses'] = self.sweep_diamond_thicknesses
        except:
            print 'could not save sweep parameters'

        if f.__contains__('/2d_plot_result/ms_errors'):
            del f['/2d_plot_result/ms_errors']        
        if f.__contains__('/2d_plot_result/u_ms_errors'):
            del f['/2d_plot_result/u_ms_errors']

        try:
            f['/2d_plot_result/ms_errors'] = self.ms_errors
            f['/2d_plot_result/u_ms_errors'] = self.u_ms_errors
        except:
            print 'could not save ms errors'

        ########################saving parameters in dictionary

        for k in self.ana_pars:
            if self.ana_pars[k]==None: #hdf5 file cannot contain attribute of type None
                g.attrs[k] = 'None'
            else:
                g.attrs[k] = self.ana_pars[k]
                

        f.close()

    def load_analysis(self,fname='analysis.hdf5'):
        f = h5py.File(os.path.join(self.folder, fname), 'r+')

        try:
            self.peak_frq = f['/2d_plot_result/peak_frq'][:]
            self.peak_V = f['/2d_plot_result/peak_V'][:]
            self.peak_filenr = f['/2d_plot_result/peak_filenr'][:]
        except:
            print 'could not load peak locations'

        try:
            self.sweep_air_lengths  = f['/2d_plot_result/sweep_air_lengths'][:]
            self.sweep_diamond_thicknesses = f['/2d_plot_result/sweep_diamond_thicknesses'][:]  
        except:
            print 'could not load sweep parameters'

        try:
            self.ms_errors = f['/2d_plot_result/ms_errors'][:]
            self.u_ms_errors = f['/2d_plot_result/u_ms_errors'][:]
        except:
            print 'could not load ms errors'

        for k in f['/2d_plot_result'].attrs:
            if f['/2d_plot_result'].attrs[k] == 'None': #hdf5 file cannot contain attribute of type None
                self.ana_pars[k] = None
            else:
                self.ana_pars[k] = f['/2d_plot_result'].attrs[k]

        f.close()