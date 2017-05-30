### General functions for analysing spectrometer data
import numpy as np 
import scipy
import pandas as pd
import os
from stat import S_ISREG, ST_MTIME, ST_MODE
import matplotlib.pyplot as plt

from analysis.scripts.cavity import peakdetect; reload(peakdetect)
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
c=3.e8#speed of light

class spectrometer_analysis(object):
    def __init__(self,folder):
        self.folder = folder
        self.ana_pars = {}
        self.set_default_params()

    def set_default_params(self):
        #set some defaults for ana_pars
        self.ana_pars['minimum_peak_height']=0 #the minimum height above average
        self.ana_pars['minimum_peak_distance']=60 #units: number of datapoints
        self.ana_pars['min_gamma']=0.0
        self.ana_pars['max_gamma']=2.4
        self.ana_pars['g_gamma']=0.2
        self.ana_pars['keep_same_height']=True #keep peaks closer together than min_peak_dist, if the same height(!)
        self.ana_pars['max_height_diff'] = 0.28#0.2#120 #if keep_same_height, the max smaller the smaller one is allowed to be.
        self.ana_pars['rescale_to_avg']=True

    def get_files_from_folder(self):
        """
        function that searches for csv files in a folder, and returns (creation time, filepath) for each
        Input parameters:
        folder - the folder in which to search for files
        Output parameters:
        files - filepath for each csv file in folder
        """
        csv_files = (os.path.join(self.folder, fn) for fn in os.listdir(self.folder) if fn.endswith('.csv'))

        self.files = np.array([])
        for f in csv_files:
            self.files = np.append(self.files,f)
        return self.files

    def sort_files_alphabetically(self):
        self.files = sorted(self.files)

    def load_data(self,filepath):
        """
        Loading csv -files into pandas dataframe and convert to array
        Input: 
        filepath  - the filepath of csv file
        Output: 
        frequencies  -  numpy array with frequencies
        intensity - nuympy array with the averaged intensity per wavelength
        """
        data = pd.read_csv(filepath, usecols=[2,5]) #creating a dataframe in pandas and importing the data
        #print data
        #the data from the spectrometer contains several intensity values per wavelength. 
        #This corresponds to rows in the camera.
        #Here we group by wavelength, and use the mean of the intensity per row.
        intensity=data.groupby('Wavelength')
        intensity=intensity.agg([np.mean])
        #first datapoint in spectrometer is wrong - remove it.
        self.wavelengths = np.array(intensity.index.tolist()[1:])

        self.frequencies = c/(self.wavelengths*1.e-9)/1.e12 #frequency in THz
        intensity = intensity.as_matrix()
        #print len(intensity)
        intensity = intensity[1:]
        return self.frequencies,intensity

    def load_data_from_folder(self):
        """
        Functions that loads all the data in csv file from a folder
        Input parameters:
        folder - the folder from which to load data
        Output parameters:
        frequencies - a numpy array with the frequencies
        filenumbers - a numpy array with the filenumbers
        intensities - a 2d numpy array with the intensities per (frequencies, filenumbers)
        """
        self.get_files_from_folder()
        nr_of_files = 0
        ii=0

        for path in sorted(self.files):
            if 'raw' in path:
                continue
            if ii == 0:
                self.frequencies,self.intensities = self.load_data(path)
                
            else:
                frequency,intensity = self.load_data(path)
                if len(intensity)!= len(self.intensities[:,0]):
                    print 'Warning: data file {} has unequal data length {:d}, compared to previous data files with length {:d}'.format(os.path.split(path)[1], len(intensity), len(intensities[:,0]))
                    continue
                self.intensities = np.concatenate((self.intensities,intensity),axis =1)
            ii+=1

            nr_of_files+=1

        self.filenumbers = np.arange(nr_of_files)

        return self.frequencies,self.filenumbers,self.intensities


    ##########Functions for peak finding and fitting##################
    def approximate_peak_location(self,intensity,**kw):
        """
        Function that detects a rough location of a peak, using scipy's argrelextrema.
        Input parameters:
        frequencies - frequency data
        intensity - intensity data 
        minimum_peak_distance - in number of datapoints, default - 30
        minimum_peak_height - the minimum height above average; default=0.2*(max intensity)
        kpsh - "keep same height" - keeps peaks with the same height, even if they are <mpd away

        Output parameters:
        peak_wavelengths - the wavelengths at which a peak is found
        """
        minimum_peak_height = kw.pop('minimum_peak_height',self.ana_pars['minimum_peak_height'])
        minimum_peak_distance = kw.pop('minimum_peak_distance',self.ana_pars['minimum_peak_distance'])
        kpsh = kw.pop('kpsh',self.ana_pars['keep_same_height'])
        max_height_diff = kw.pop('max_height_diff',self.ana_pars['max_height_diff'])#the maximum difference between peak height, in order that they are kept if kpsh
        rescale_to_avg = kw.pop('rescale_to_avg',self.ana_pars['rescale_to_avg'])

        if not rescale_to_avg:
            minimum_peak_height=np.average(intensity)+minimum_peak_height

        peak_frequencies = np.array([])
        peak_intensity = np.array([])
        indices = peakdetect.detect_peaks(intensity,mph=minimum_peak_height,mpd=minimum_peak_distance,
            kpsh=kpsh,mhd=max_height_diff, rescale_to_avg=rescale_to_avg)
        # print indices_maxima
        for ii in indices:
            peak_frequencies = np.append(peak_frequencies,self.frequencies[ii])
            peak_intensity = np.append(peak_intensity,intensity[ii])

        return indices,peak_frequencies, peak_intensity

    def fit_peak(self,intensity,indices,peak_frequencies,peak_intensity, 
            plot_fit = False, **kw):
        """
        This function fits every presumed peak location with a lorentzian. 
        If the fit fails it rejects it as a peak.
        Input parameters:
        intensity - 1 d array of intensity data
        indices - indices of the peaks
        peak_frequencies - frequencies of the peaks
        peak_intensity - intensity of the peaks
        g_gamma - the guess parameter of Lorentzians FWHM. Default: 0.2
        g_offset - the guess parameter of the offset of the Lorentzian. Default:0
        plot_fit - whether to plot each fit. default = False
        Output parameters:
        x0s - 1d array of the fitted peak locations
        u_x0s -  1d array of the uncertainty in the fitted peak locations
        """
        report_fails = kw.pop('report_fails',False)

        x0s =np.array([])
        u_x0s =np.array([])
        g_gamma = kw.pop('g_gamma',self.ana_pars['g_gamma'])
        max_gamma = kw.pop('max_gamma',self.ana_pars['max_gamma'])
        min_gamma = kw.pop('min_gamma',self.ana_pars['min_gamma'])
        
        frequency_range = np.abs(self.frequencies[-1]-self.frequencies[0])
        indices_around_peak = int((len(self.frequencies)/frequency_range)*g_gamma*6)
        success = np.zeros(len(indices))
        nr_fails = 0

        for i,ii,g_x0,g_A in zip(np.arange(len(indices)),indices, peak_frequencies, peak_intensity*g_gamma):
            if ii - indices_around_peak <0:
                i_min = 0
            else:
                i_min = ii - indices_around_peak
            if ii + indices_around_peak > len(self.frequencies)-1:
                i_max = -1
            else:
                i_max = ii + indices_around_peak


            intensity_around_peak = intensity[i_min:i_max]
            frequencies_around_peak = self.frequencies[i_min:i_max]
            g_offset = np.average(intensity_around_peak)
            fixed = []

            p0, fitfunc, fitfunc_str = common.fit_lorentz(g_offset, g_A, g_x0, g_gamma)
            fit_result = fit.fit1d(frequencies_around_peak,intensity_around_peak, None, p0=p0, 
                fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
            
            #if the fit failed, it outputs the variable 'success', that is an integer flag.
            #If it is 1,2,3,4 the fit succeeded, otherwise it failed. 
            #Break the loop if the fit failed
            if fit_result['success']==False:
                nr_fails += 1
                continue
                
            if fit_result == 5:  
                #print 'fit failed'
                nr_fails += 1
                continue

            res_rms = fit_result['residuals_rms']/np.average(frequencies_around_peak)
            gamma = fit_result['params_dict']['gamma']
            u_gamma = fit_result['error_dict']['gamma']
            x0 = fit_result['params_dict']['x0']
            u_x0 = fit_result['error_dict']['x0']
            A = fit_result['params_dict']['A']
            # print x0,A,gamma
     
            if plot_fit:
                # print indices_around_peak
                # print wavelengths_around_peak
                print 'plot fit'
                fig,ax = plt.subplots()
                # ax.plot(wavelengths,intensity)
                plot.plot_fit1d(fit_result, np.linspace(frequencies_around_peak[0],frequencies_around_peak[-1],len(frequencies_around_peak)), 
                    ax =ax, label='Fit',show_guess=True, plot_data=True)
                plt.show()
            # if A < 20:
            #     print 'peak intensity is too low; disregarding'
            #     continue

            if u_x0 > np.abs(frequencies_around_peak[-1]-frequencies_around_peak[0]):
                if plot_fit:
                    print 'uncertainty in peak position too large; disregarding: ', 'x0', x0, '+-',u_x0 
                nr_fails+=1
                continue 

            if u_gamma > np.abs(gamma):
                if plot_fit:
                    print 'uncertainty in gamma too large; disregarding: ', 'gamma', gamma, '+-',u_gamma 
                nr_fails+=1
                continue                 
            
            if ((np.abs(gamma)>max_gamma) or (np.abs(gamma)<min_gamma)): 
                if plot_fit:
                    print 'ignoring this peak, since gamma is not within specs:',min_gamma, '<',  gamma, '>', max_gamma
                nr_fails+=1
                continue

            if A*gamma < 0:
                if report_fails or plot_fit:
                    print 'ignoring since negative '
                nr_fails+=1
                continue

            success[i] = 1 #mark this peak as succesfully fitted
            x0s = np.append(x0s,x0)
            u_x0s = np.append(u_x0s,u_x0)
        if report_fails:
            print 'number of failed fits:', nr_fails
        return x0s, u_x0s, success

# import os
# import sys
# import numpy as np
# sys.path.append("H:\My Documents\measuring/")

# %matplotlib inline

# import analysis.scripts.cavity.oscilloscope_analysis as oa
# import analysis.scripts.cavity.fit_oscilloscope_data as od


# data_dir = "K:/ns\qt\Diamond\Projects\Cavities\data/20160426/EOM_lw_LT/"
# file_name = "NNNNNNNLOWT007.csv"
# filename = os.path.join(data_dir,file_name)
# EOM_freq = 6 
# reload(oa)
# reload(od)
# data = oa.load_data(filename)
# od.get_linewidth(data,EOM_freq)