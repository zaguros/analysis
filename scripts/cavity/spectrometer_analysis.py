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

def get_files_from_folder(folder):
    """
    function that searches for csv files in a folder, and returns (creation time, filepath) for each
    Input parameters:
    folder - the folder in which to search for files
    Output parameters:
    files - (creation time, filepath) for each csv file in folder
    """
    files = (os.path.join(folder, fn) for fn in os.listdir(folder) if fn.endswith('.csv'))
    files = ((os.stat(path), path) for path in files)
    files = ((stat[ST_MTIME], path)
               for stat, path in files if S_ISREG(stat[ST_MODE]))
 
    return files

def load_data(filepath):
    """
    Loading csv -files into pandas dataframe and convert to array
    Input: 
    filepath  - the filepath of csv file
    Output: 
    wavelengths  -  numpy array with wavelengths
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
    wavelengths = np.array(intensity.index.tolist()[1:])

    frequencies = c/(wavelengths*1.e-9)/1.e12 #frequency in THz
    
    intensity = intensity.as_matrix()
    #print len(intensity)
    intensity = intensity[1:]
    return frequencies,intensity

def load_data_from_folder(folder):
    """
    Functions that loads all the data in csv file from a folder
    Input parameters:
    folder - the folder from which to load data
    Output parameters:
    wavelengths - a numpy array with the wavelengths
    filenumbers - a numpy array with the filenumbers
    intensities - a 2d numpy array with the intensities per (wavelengths, filenumbers)
    """
    files = get_files_from_folder(folder)
    nr_of_files = 0
    ii=0
    for cdate,path in sorted(files):
        if 'raw' in path:
            continue
        if ii == 0:
            wavelengths,intensities = load_data(path)
        else:
        #try:  
            wavelength,intensity = load_data(path)
            if len(intensity)!= len(intensities[:,0]):
                print 'Warning: data file {} has unequal data length {:d}, compared to previous data files with length {:d}'.format(os.path.split(path)[1], len(intensity), len(intensities[:,0]))
                continue
            intensities = np.concatenate((intensities,intensity),axis =1)
        ii+=1
        #
            #wavelengths,intensities = load_data(path)

            # wavelengths = data[:,0]
            # intensities = np.array([data[:,1]])
        nr_of_files+=1

    filenumbers = np.arange(nr_of_files)

    return wavelengths,filenumbers,intensities


##########Functions for peak finding and fitting
def approximate_peak_location(wavelengths,intensity,**kw):
    """
    Function that detects a rough location of a peak, using scipy's argrelextrema.
    Input parameters:
    wavelengths - wavelength data
    intensity - intensity data 
    minimum_peak_distance - in number of datapoints, default - 30
    minimum_peak_height - default=0.2*(max intensity)

    Output parameters:
    peak_wavelengths - the wavelengths at which a peak is found
    """
    minimum_peak_height = kw.pop('minimum_peak_height',0.2*max(intensity))
    minimum_peak_distance = kw.pop('minimum_peak_distance',30)
    peak_wavelengths = np.array([])
    peak_intensity = np.array([])
    indices = peakdetect.detect_peaks(intensity,mph=minimum_peak_height,mpd=minimum_peak_distance)
    # print indices_maxima
    for ii in indices:
        peak_wavelengths = np.append(peak_wavelengths,wavelengths[ii])
        peak_intensity = np.append(peak_intensity,intensity[ii])

    return indices,peak_wavelengths, peak_intensity

def fit_peak(wavelengths,intensity,indices,peak_wavelengths,peak_intensity, 
        plot_fit = False, **kw):
    """
    This function fits every presumed peak location with a lorentzian. 
    If the fit fails it rejects it as a peak.
    Input parameters:
    wavelengths - wavelength data
    intensity - 1 d array of intensity data
    indices - indices of the peaks
    peak_wavelengths - wavelengths of the peaks
    peak_intensity - intensity of the peaks
    g_gamma - the guess parameter of Lorentzians FWHM. Default: 0.2
    g_offset - the guess parameter of the offset of the Lorentzian. Default:0
    plot_fit - whether to plot each fit. default = False
    Output parameters:
    x0s - 1d array of the fitted peak locations
    u_x0s -  1d array of the uncertainty in the fitted peak locations
    """
    x0s =np.array([])
    u_x0s =np.array([])
    g_gamma = kw.pop('g_gamma',0.5)
    g_offset= kw.pop('g_offset',0)
    max_gamma = kw.pop('max_gamma',None)
     
    wavelength_range = np.abs(wavelengths[-1]-wavelengths[0])
    indices_around_peak = int((len(wavelengths)/wavelength_range)*g_gamma*4)
    success = np.zeros(len(indices))
    nr_fails = 0

    for i,ii,g_x0,g_A in zip(np.arange(len(indices)),indices, peak_wavelengths, peak_intensity):
        if ii - indices_around_peak <0:
            i_min = 0
        else:
            i_min = ii - indices_around_peak
        if ii + indices_around_peak > len(wavelengths)-1:
            i_max = -1
        else:
            i_max = ii + indices_around_peak


        intensity_around_peak = intensity[i_min:i_max]
        wavelengths_around_peak = wavelengths[i_min:i_max]

        fixed = []

        p0, fitfunc, fitfunc_str = common.fit_lorentz(g_offset, g_A, g_x0, g_gamma)
        fit_result = fit.fit1d(wavelengths_around_peak,intensity_around_peak, None, p0=p0, 
            fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
        
        #if the fit failed, it outputs the variable 'success', that is an integer flag.
        #If it is 1,2,3,4 the fit succeeded, otherwise it failed. 
        #Break the loop if the fit failed
        if fit_result == 5:  
            #print 'fit failed'
            nr_fails += 1
            continue

        res_rms = fit_result['residuals_rms']/np.average(wavelengths_around_peak)
        gamma = fit_result['params_dict']['gamma']
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
            plot.plot_fit1d(fit_result, np.linspace(wavelengths_around_peak[0],wavelengths_around_peak[-1],len(wavelengths_around_peak)), 
                ax =ax, label='Fit',show_guess=True, plot_data=True)
            plt.show()
        # if A < 20:
        #     print 'peak intensity is too low; disregarding'
        #     continue

        if u_x0 > np.abs(wavelengths_around_peak[-1]-wavelengths_around_peak[0]):
            # print 'uncertainty in peak position too large; disregarding: ', 'x0', x0, '+-',u_x0 
            nr_fails+=1
            continue 
        
        if max_gamma!=None:
            if gamma>max_gamma: 
                # print 'ignoring this peak, since gamma is larger than max gamma:', gamma, '>', max_gamma
                nr_fails+=1
                continue

        success[i] = 1 #mark this peak as succesfully fitted
        x0s = np.append(x0s,x0)
        u_x0s = np.append(u_x0s,u_x0)

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