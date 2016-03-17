####################################################################################################################################################################
# A.M.J. Zwerver
#a.m.j.zwerver@student.tudelft.nl
# October 2015

import os
import h5py
import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
import pylab as plt
from analysis.lib import fitting
from analysis.lib.fitting import esr, common
from analysis.lib.fitting import fit,esr, common
from matplotlib import rc, cm
from scipy import interpolate
from scipy.optimize import curve_fit
import operator

from analysis.lib.m2 import m2

class cavity_analysis(m2.M2Analysis):

    '''
    The next three functions are the basis function to analyse the data. They gather the data that you want to collect, assign x and y 
    (Volt/Frequency and photodiode signal)
    and return the x and y values of the indicated scans.
    '''

    ####write the below function such that it works with the m2.m2Analysis framweork
    def load_data (self,folder, timestamp, scan_type, return_folder = False):
        '''This function selects the folder with the defined timestamp in the folder of the indicated day. 
        Thereafter it picks the hdf5 file from this folder (the file that contains the data) and selects the x and y data 
        (which is different for laser scans and piezo scans).'''

        all_folders = [f for f in os.listdir(folder) if timestamp in f and scan_type in f]
        test = all_folders[0]
        if len(all_folders) == 0:
            x = []
            y = []
            f = None
        
        else:
            curr_fold = os.path.join(folder, all_folders[0])
            all_files =[f for f in os.listdir(curr_fold) if os.path.isfile(os.path.join(curr_fold,f)) ]
            file_name = [f for f in all_files if '.hdf5' in f]

            if (scan_type =='piezo'):
                grp_name = '/piezo_scan'
                x_name = 'piezo_voltage'
            elif (scan_type == 'lr_scan'):
                grp_name = '/lr_scan'
                x_name = 'frequency_GHz'
            elif (scan_type =='fine_laser'):
                grp_name = '/fine_laser_scan'
                x_name = 'laser_tuning_voltage'
            
            f = h5py.File(os.path.join(curr_fold, file_name[0]),'r')
            ls_grp = f[grp_name]
            x = ls_grp[x_name].value
            y = ls_grp['PD_signal'].value
        if return_folder:
            return curr_fold, x, y, test
        else:
            return x, y, test

    def get_x_pts(self):
        self.nr_x_pts = self.f.attrs['nr_steps']
        self.x_pts = np.linspace(self.f.attrs['start_voltage'],self.f.attrs['end_voltage'],self.f.attrs['nr_steps'])

    def get_sweep_pts(self):
        self.sweep_name = self.f.attrs['sweep_name']
        self.sweep_pts = self.f.attrs['sweep_pts']
        self.nr_syncs_per_pt = self.f.attrs['nr_syncs_per_pt']
        print self.sweep_name, self.sweep_pts

    def get_sweep_data (self):
        """
        function that loads the data from a sweep measurement
        """
        self.nr_scans = self.f.attrs['nr_scans']
        self.nr_remainder = self.f.attrs['nr_remainder']
        self.nr_repetitions = self.f.attrs['nr_repetitions']

        self.sweep_data=np.zeros([len(self.sweep_pts),self.nr_repetitions,self.nr_x_pts])

        for i in np.arange(self.nr_syncs_per_pt):
            if (i == self.nr_syncs_per_pt-1): #the last one
                if self.nr_remainder>0:
                    self.nr_scans = self.nr_remainder
            for j in np.arange(len(self.sweep_pts)):
                first_rep = str(int(i*self.nr_scans+1))
                last_rep = str(int(i*self.nr_scans+self.nr_scans))
                data_index_name = 'sweep_pt_'+str(j)+'_reps_'+first_rep+'-'+last_rep  
                grp = self.f['raw_data_'+data_index_name]
                for k in np.arange(self.nr_scans):
                    self.data = grp['scannr_'+str(k+1)].value
                    #fill the sweep_data array for j = sweep pt, i*self.nr_scans+k = repetition nr
                    self.sweep_data[j,int(i*self.nr_scans+k)] = self.data


        self.avg_sweep_data = np.average(self.sweep_data,axis=1)

        return self.sweep_data

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

def load_piezo_scan (folder, timestamp, return_folder = False):
    "Loads the piezo data and orders it by voltage"
    data_folder, V, y, name = load_data (folder = folder, timestamp = timestamp, scan_type = 'piezo',return_folder=True)
    V = np.ravel(V)
    y = np.ravel(y)

    # ind = np.where(abs(V)>5000)
    # V = np.delete(V, ind)
    # y = np.delete(y, ind)
    if return_folder:
        return data_folder, V, y, name
    else:
        return V, y, name

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

def load_fine_lr_scan (folder, timestamp, return_folder = False):
    "Loads the fine laser data and orders it by voltage"
    data_folder, V, y, name = load_data (folder = folder, timestamp = timestamp, scan_type = 'fine_laser', return_folder = True)
    V = np.ravel(V)
    y = np.ravel(y)

    # ind = np.where(abs(V)>5000)
    # V = np.delete(V, ind)
    # y = np.delete(y, ind)
    if return_folder:
        return data_folder, V, y, name
    else:
        return V, y, name

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

def load_lr_scan (folder, timestamp):
    "loads the laser data, orders it by frequency and deletes all points above a certain absolute value, for the laser scans make a big jump sometimes"
    f, y, name = load_data (folder = folder, timestamp = timestamp, scan_type = 'lr_scan')
    f = np.ravel(f)
    y = np.ravel(y)

    ind = np.where(abs(f)>1000)
    f = np.delete(f, ind)
    y = np.delete(y, ind)
    return f, y, name

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

def show_fine_laser_plots(date, timestamp, folder, measurement_type, f_min, f_max, set_range_f):
    'This function plots all fine laser scans within the given timestamp range'
    T = []
    LW = []

    for i in range(len(timestamp)):
        time = timestamp[i]
        f,y, name = load_fine_lr_scan(folder = folder, timestamp = time)

        if len(f) == 0:
            continue

        T.append(1)
        f = f*35./3.
        
        ''' Assigns f_min and f_max '''
        if set_range_f == True:
            ind_min = np.where(f < f_min)
            ind_max = np.where(f > f_max)
            f = np.delete(f, ind_max[0])
            f = np.delete(f, ind_min[0])
            y = np.delete(y, ind_max[0])
            y = np.delete(y, ind_min[0])
        else:
            f_max = max(f)
            f_min = min(f)
        # print name
        fig,ax = plt.subplots(figsize=(6,4.7))
        ax.plot(f,y)
        plt.xlabel('Relative frequency [GHz]')
        plt.ylabel('Transmission')
        ax.set_xlim(f_min,f_max)
        # plt.title(date + ' - ' + time + '\n' + measurement_type + '\n Fine laser plot without fit')
        plt.title(name)
        plt.show()
        
    if len(T) == 0:
        print "No fine laser data with these timestamps"

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

def show_piezo_plots(date, timestamp, folder, measurement_type, V_min, V_max, set_range_V):
    'This function plots all piezo scans within the given timestamp range'
    T = []
    for i in range(len(timestamp)):
        time = timestamp[i]

        V,y, name = load_piezo_scan(folder = folder, timestamp = time)

        if len(V) == 0:
            continue

        T.append(1)
        ''' Assigns V_min and V_max '''
        if set_range_V == True:
            ind_min = np.where(V < V_min)
            ind_max = np.where(V > V_max)
            V = np.delete(V, ind_max[0])
            V = np.delete(V, ind_min)
            y = np.delete(y, ind_max[0])
            y = np.delete(y, ind_min)
        else:
            V_max = max(V)
            V_min = min(V)

        fig,ax = plt.subplots(figsize=(6,4.7))
        ax.plot(V,y)
        plt.xlabel('Voltage [V]')
        plt.ylabel('Transmission')
        ax.set_xlim(V_min,V_max)
        plt.title(date + ' - ' + time + '\n' + measurement_type + '\n Piezo plot without fit')
    if len(T) ==0:
        print "No piezo data with these timestamps"

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

def show_laser_plots(date, timestamp, folder, measurement_type, f_min, f_max, set_range_f):
    'This function plots all laser scans within the given timestamp range'
    
    T = []
    for i in range(len(timestamp)):
        time = timestamp[i]
        f,y, name = load_lr_scan (folder=folder, timestamp = time)

        if len(f) == 0:
            continue

        T.append(1)
        ''' Assigns f_min and f_max '''
        if set_range_f == True:
            ind_min = np.where(f < f_min)
            ind_max = np.where(f > f_max)
            f = np.delete(f, ind_min)
            f = np.delete(f, ind_max)
            y = np.delete(y, ind_min)
            y = np.delete(y, ind_max)

        else:
            f_min = min(f)
            f_max = max(f)


        fig,ax = plt.subplots(figsize=(6,4.7))
        ax.plot(f,y)
        plt.xlabel('frequency [GHz]')
        plt.ylabel('Transmission')
        ax.set_xlim(f_min,f_max)
        plt.title(date + ' - ' + time + '\n' + measurement_type + '\n Laser plot without fit')
    if len(T) ==0:
        print "No laser data with these timestamps"

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

def get_FSR(date, timestamp, folder, newpath, V_min, V_max, set_range_V, labda = 637, save_data = False):
    ''' Plots all laser scans within the given timestamp range and assigns the FSR_L (in V) and the conversion factor, given FSR_L = lambda/2'''
    
    print '### YOU\'RE ANALYSING THE FSR ### \n'
    T = []
    ALPHA = []
    Delta_V = []
    
    for i in range(len(timestamp)):

        time = timestamp[i]

        '''Determine where the peaks are and select the resonance peaks'''
        V,y, name = load_piezo_scan(folder = folder, timestamp = time)
        V=-1*V
        if len(V) == 0:
            continue



        if len(V) == 0:
            continue

        T.append(1)
        ''' Assigns V_min and V_max '''
        if set_range_V == True:
            ind_min = np.where(V < V_min)
            ind_max = np.where(V > V_max)
            V = np.delete(V, ind_max[0])
            V = np.delete(V, ind_min)
            y = np.delete(y, ind_max[0])
            y = np.delete(y, ind_min)
        else:
            V_max = max(V)
            V_min = min(V)

        T.append(1)
        peak = np.where(y > 0.05)
        resonance1 = peak[0][0]
        resonance2 = peak[0][-1]
        V1 = V[resonance1]
        V2 = V[resonance2]

        DV = V2-V1
        Delta_V.append(DV)

        # fsr = DV * 2500/12
        # Length = labda**2/(2*fsr)
        # FSR.append(fsr)
        # L.append(Length)

        # print 'FSR and LENGTH',fsr, Length
        ''' (m/labda = 2L), therefore, for two consecutive peaks (m and m+1), (Delta_L = labda/2). And since Delta_L = alpha*Delta_V: alpha = labda/(2*Delta_V)'''
        alpha = labda/(2.*(V2-V1))
        ALPHA.append(alpha)
        print 'DELTA VOLTAGE:', DV
        print 'alpha is ', alpha

        ''' Plot the Free Spectral Range '''        
        fig, ax = plt.subplots(figsize=(6,4.7))
        ax.plot(V,y,'b-')
        ax.set_xlim(V_min,V_max)
        plt.plot(V1,0,'ro')
        plt.plot(V2,0,'ro')
        plt.xlabel('Volt [V]',fontsize=14)
        plt.ylabel('transmission [arb. units]',fontsize=14)
        plt.rc('text', usetex = True)
        plt.title(date + ' - ' + time + '\n Free Spectral Range is ' + str(round(DV,2)) + 'V \n' + r'$\alpha$' + ' is ' + str(round(alpha,2)) + 'nm/V')
        # plt.title('Free Spectral Range',fontsize=14)
        # plt.axis('off')
        plt.show()
        
        if save_data == True:
            fig.savefig(newpath +'/' + "FSR " + time + ".png")

    average_Delta_V = sum(Delta_V)/len(Delta_V)
    variance_Delta_V = sum([(x-average_Delta_V)**2 for x in Delta_V])/len(Delta_V)
    st_dev_Delta_V = np.sqrt(variance_Delta_V)
    
    average_ALPHA = sum(ALPHA)/len(ALPHA)
    variance_alpha = sum([(x-average_ALPHA)**2 for x in ALPHA])/len(ALPHA)
    st_dev_alpha = np.sqrt(variance_alpha)
    
    W = 182.5*average_Delta_V
    st_W = 182.5*st_dev_Delta_V


    print 'Average alpha is ', round(average_ALPHA,2), 'with a standard deviation of ', round(st_dev_alpha,2)
    print 'Average Delta_V is ', round(W,2), 'with a standard deviation of ', round(st_W,4)




    if len(T) ==0:
        print "No piezo data with these timestamps"

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

def get_piezo_linewidth(date, timestamp, folder, measurement_type, V_min, V_max, newpath, set_range_V = False, save_data = False):
    '''
    This function loads the data for piezo scans, determines the peak, the mean of the peak and the standard deviation of this mean.
    The linewidth is defined as three times the standard deviation of the mean. The scan is plotted and an indication of the determined
    linewidth is given by two green dots.
    '''

    print '### YOU\'RE ANALYSING PIEZO SCANS ### \n'
    LW = []
    Quality = []
    t = []
    drift = []
    drift_error = []
    frequency = []
    standard_dev = []
    T = []

    for i in range(len(timestamp)):
        time = timestamp[i]
        V,y, name = load_piezo_scan(folder = folder, timestamp = time)
        # V,y = zip(*sorted(zip(V, y)))

        if len(V) == 0:
            continue

        T.append(1)
        '''
        If a certain range for V is selected, delete all the datapoints outside this range
        '''
        if set_range_V == True:
            ind_min = np.where(V < V_min)
            ind_max = np.where(V > V_max)
            # print len(V)
            V = np.delete(V, ind_max[0])
            V = np.delete(V, ind_min)
            y = np.delete(y, ind_max[0])
            y = np.delete(y, ind_min)
        else:
            V_max = max(V)
            V_min = min(V)

        ''' get the maximum value in the data and the index of this value (this should be the top of the linewidth peak) '''
        max_index, max_value = max(enumerate(y), key = operator.itemgetter(1))

        ''' Determine mean and standard deviation of the peaks '''
        peak = np.where(y > 0.1)
        if len(peak[0]) > 0:
            print time + ' Peak detected'
            X = [V[peak[0][i]] for i in range(len(peak[0]))]
            mean = sum(X)/len(X)
            variance = sum([(x-mean)**2 for x in X])/len(X)
            st_dev = np.sqrt(variance)
            st_dev_3 = 3*st_dev
            X_left = mean - (st_dev_3/2)
            X_right = mean + (st_dev_3/2)
            LW.append(st_dev_3)
            Ra = round(st_dev_3*111,3)
            Ro = round(st_dev_3*3,2)
            ''' Make a plot of the data '''
            fig,ax = plt.subplots(figsize=(6,4.7))
            ax.plot(V,y,'b', label = 'LT, tube low')
            # ax.legend()
            ax.legend(loc=2)
            # plt.plot(X_left,0, 'go')
            # plt.plot(X_right,0, 'go')
            plt.xlabel('Voltage [V]',fontsize=14)
            plt.ylabel('Transmission [arb. units]',fontsize=14)
            ax.set_xlim(V_min,V_max)
            # plt.title(date + ' - ' + time + '\n' + measurement_type + '\n Linewidth is %s V' %(round(st_dev_3,5)))
            plt.title('\n Range peak is %s $\pm$ %s nm' %(round(st_dev_3*111,2), round(st_dev_3*3,2)))

            if save_data == True:
                fig.savefig(newpath +'/' + "piezolinewidth " + time + ".png")

        else:
            print time + " No significant peak detected"
            continue
        
            
    ''' Determine linewidth by means of sigma '''

    if len(LW) > 0:
        average_sigma = sum(LW)/float(len(LW))
        variance_average_sigma = sum([(x-average_sigma)**2 for x in LW])/len(LW)
        st_dev_average_sigma = np.sqrt(variance_average_sigma)
        print 'AVERAGE SIGMA IS ', round(average_sigma,5), ' [V] WITH A ST DEV OF ', round(st_dev_average_sigma,5), ' [V]'
    elif len(T) == 0:
        print "No piezo data with these timestamps"
    else:
        print "\n !!! The average linewidth could not be determined, for there was no significant peak detected !!!"

        

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

def get_laser_linewidth(date, timestamp, folder, measurement_type, threshold, f_min, f_max, newpath, set_range_f = False, save_data = False):

    print '### YOU\'RE ANALYSING LASER SCANS ### \n'

    LW = []
    Quality = []
    t = []
    drift = []
    drift_error = []
    delete = []
    T = []

    for i in range(len(timestamp)):
        time = timestamp[i]
        f,y, name = load_lr_scan (folder=folder, timestamp = time)

        if len(f) == 0:
            continue

        T.append(1)

        print "f: ", len(f), f[0], f[-1]

        if set_range_f == True:
            ind_min = np.where(f < f_min)
            ind_max = np.where(f > f_max)
            f = np.delete(f, ind_min)
            f = np.delete(f, ind_max)
            y = np.delete(y, ind_min)
            y = np.delete(y, ind_max)
        else:
            f_min = min(f)
            f_max = max(f)


        print "f: ", len(f), f[0], f[-1]

        max_index, max_value = max(enumerate(y), key = operator.itemgetter(1))
        if max_value < threshold:
            continue

        '''
        Determine mean and standard deviation of the peaks
        '''
        peak = np.where(y >= threshold)
        X = [f[peak[0][i]] for i in range(len(peak[0]))]
        mean = sum(X)/len(X)
        variance = sum([(x-mean)**2 for x in X])/len(X)
        st_dev = np.sqrt(variance)
        st_dev_3 = 3*st_dev
        X_left = mean - (st_dev_3/2)
        X_right = mean + (st_dev_3/2)
        LW.append(st_dev_3)

        '''
        Make a plot of the data
        '''
        fig,ax = plt.subplots(figsize=(6,4.7))
        ax.plot(f,y,'r', label = 'piezo scan linewidth')
        ax.legend()
        plt.plot(X_left,0, 'go')
        plt.plot(X_right,0, 'go')
        plt.xlabel('frequency [GHz]')
        plt.ylabel('Transmission')
        ax.set_xlim(f_min,f_max)
        plt.title(date + ' - ' + time + '\n' + measurement_type + '\n Linewidth is %s GHz' %(round(st_dev_3,2)))

        if save_data == True:
            print newpath
            fig.savefig(newpath +'/' + "laserlinewidth " + time + ".png")

    ''' 
    Determine linewidth by means of sigma
    '''
    if len(LW) > 0:
        average_sigma = sum(LW)/float(len(LW))
        variance_average_sigma = sum([(x-average_sigma)**2 for x in LW])/len(LW)
        st_dev_average_sigma = np.sqrt(variance_average_sigma)
        print 'AVERAGE SIGMA IS ', round(average_sigma,5), 'WITH A ST DEV OF ', round(st_dev_average_sigma,5)
    else:
        print "\n !!! The average linewidth could not be determined, for there was no significant peak detected !!!"

    if len(T) ==0:
        print "No laser data with these timestamps"

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

def fit_laser_linewidth(date, timestamp, folder, measurement_type, threshold, f_min, f_max, newpath, set_range_f = False, save_data = False):

    print '### YOU\'RE FITTING LASER SCANS ### \n'
    
    LW = []
    Quality = []
    t = []
    drift = []
    drift_error = []
    delete = []
    T = []

    for i in range(len(timestamp)):
        time = timestamp[i]
        f,y, name = load_lr_scan (folder=folder, timestamp = time)

        if len(f) == 0:
            continue

        T.append(1)

        if set_range_f == True:
            ind_min = np.where(f < f_min)
            ind_max = np.where(f > f_max)
            f = np.delete(f, ind_min)
            f = np.delete(f, ind_max)
            y = np.delete(y, ind_min)
            y = np.delete(y, ind_max)
        else:
            f_min = min(f)
            f_max = max(f)

        max_index, max_value = max(enumerate(y), key = operator.itemgetter(1))
        if max_value < threshold:
            continue

        # f,y = zip(*sorted(zip(f, y)))

        offset = 0.01
        amplitude = max_value
        x0 = f[max_index]
        exponent = 1
        decay_constant = 3
        sigma = 2.413625
        gamma = 20

        """ Fit a Gaussian """
        fixed = [0]
        p0, fitfunc, fitfunc_str = common.fit_gauss(offset, amplitude, x0, sigma)
        fit_result = fit.fit1d(f,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)

        ''' 
            In order to calculate the linewidth a help-function is created which resembles the fitting function (by means of
            the fitting parameters). Of this function all indices where the function exceeds the half maximum are stored in
            Half_Max. Of this list the very first and very last values are considered to by the indices at the half maximum and
            the linewidth (in nm) is then calculated by subtracting the voltages of those points and multiplying that by 25 nm/V, 
            a spec that is provided by JPE. Subsequently Q can be calculated by the well-known formula Q=omega/(delta omega), for 
            which omega is the wavelength of the laser and delta omega the linewidth
        '''
        
        A = fit_result['params_dict']['A']
        a = offset
        x0 = fit_result['params_dict']['x0']
        sigma = fit_result['params_dict']['sigma']

        function = a + A*np.exp(-(f-x0)**2/(2*sigma**2))

        if A > 5*max_value:
            print 'fit failed'
            continue

        Half_Max = np.where(function > A/2)
        Half_Max_l = Half_Max[0][0]
        Half_Max_r = Half_Max[0][-1]

        Linewidth = abs(f[Half_Max_r]-f[Half_Max_l])
        LW.append(Linewidth)


        ''' Plot figure '''
        fig,ax = plt.subplots(figsize=(6,4.7))
        ax.plot(f,y, 'bo', label = 'data')
        plot.plot_fit1d(fit_result, np.linspace(f[0],f[-1],len(f)), ax=ax,label='Fit',show_guess=True, plot_data=False)
        ax.plot(f,function, 'g', label='testplot')
        ax.legend()
        plt.plot(f[Half_Max_r], A/2+offset, 'go')
        plt.plot(f[Half_Max_l], A/2+offset, 'go')
        plt.xlabel('frequency [GHz]')
        plt.ylabel('Transmission')
        plt.xlim(f_min,f_max)
        plt.title(date + ' - ' + time + '\n' + measurement_type + '\n Linewidth is %s GHz' %round(Linewidth,2))

        if save_data == True:
            fig.savefig(newpath +'/' + "laserfit " + time + ".png")

    ''' 
    Determine linewidth by means of sigma
    '''
    if len(LW) > 0:
        average_sigma = sum(LW)/float(len(LW))
        variance_average_sigma = sum([(x-average_sigma)**2 for x in LW])/len(LW)
        st_dev_average_sigma = np.sqrt(variance_average_sigma)
        print 'AVERAGE SIGMA IS ', round(average_sigma,5), ' [GHz] WITH A ST DEV OF ', round(st_dev_average_sigma,5), ' [GHz]'
    elif len(T) == 0:
        print "No fine laser data with these timestamps"
    else:
        print "\n !!! The average linewidth could not be determined, for there was no significant peak detected !!!"

        


def fit_fine_laser_plots(date, timestamp, folder,  f_min, f_max, set_range_f = False, save_data = False, threshold=0.1):
    'This function plots all fine laser scans within the given timestamp range'
    T = []
    LW = []

    for i in range(len(timestamp)):
        time = timestamp[i]
        data_folder, f,y, name = load_fine_lr_scan(folder = folder, timestamp = time, return_folder = True)

        if len(f) == 0:
            continue

        T.append(1)
        f = f*30./3.
        
        ''' Assigns f_min and f_max '''
        if set_range_f == True:
            ind_min = np.where(f < f_min)
            ind_max = np.where(f > f_max)
            f = np.delete(f, ind_max[0])
            f = np.delete(f, ind_min[0])
            y = np.delete(y, ind_max[0])
            y = np.delete(y, ind_min[0])
        else:
            f_max = max(f)
            f_min = min(f)

        max_index, max_value = max(enumerate(y), key = operator.itemgetter(1))
        if max_value < threshold:
            continue

        offset = 0.01
        amplitude = max_value
        x0 = f[max_index]
        sigma = 5
        gamma = 5

        
        fixed = []
        p0, fitfunc, fitfunc_str = common.fit_lorentz(offset, amplitude, x0, gamma)
        fit_result = fit.fit1d(f,y, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)


        #print fit_result['params_dict']
        A = fit_result['params_dict']['A']
        a = fit_result['params_dict']['a']
        x0 = fit_result['params_dict']['x0']
        #sigma = fit_result['params_dict']['sigma']
        gamma = fit_result['params_dict']['gamma']
        u_gamma = fit_result['error_dict']['gamma']

        # f2 = np.linspace(min(f),max(f),1001)
        function = a + 2*A/np.pi*gamma/(4*(f-x0)**2+gamma**2)
        if A > 5  *max_value:
            print 'fit failed'
            continue

        Linewidth = np.abs(gamma)
        u_Linewidth = u_gamma
        LW.append(Linewidth)

        ''' Plot figure '''
        fig,ax = plt.subplots(figsize=(6,4.7))
        ax.plot(f,y, 'bo', label = 'data')
        plot.plot_fit1d(fit_result, np.linspace(f[0],f[-1],len(f)), ax=ax,color='r',show_guess=True, plot_data=False, label = 'fit', add_txt=False)
        #ax.plot(f,function, 'g', label='testplot')
        ax.legend()

        ax.set_xlabel('Frequency [GHz]',fontsize=14)
        ax.set_ylabel('Transmission (arb. units)',fontsize=14)
        ax.tick_params(axis = 'both', which = 'major',labelsize = 14)
        ax.set_xlim(f_min,f_max)
        ax.set_title(data_folder+'\n Linewidth in frequency is %s $\pm$ %s GHz' %(round(Linewidth,2), round(u_Linewidth,2)))
        if save_data == True:
            try:
                fig.savefig(data_folder+'/fit.png')
                fig.savefig(data_folder+'/fit.pdf')
            except:
                print('figure not saved')

        plt.show()
        plt.close()
    if len(T) ==0:
        print "No fine laser data with these timestamps"

    average_LW = sum(LW)/len(LW)
    variance_LW = sum([(x-average_LW)**2 for x in LW])/len(LW)
    st_dev_LW = np.sqrt(variance_LW)
    print 'Average LW_piezo is ', round(average_LW,1), ' GHz with a standard deviation of ', round(st_dev_LW,1), ' GHz'







def fit_fine_laser_plots2(date, timestamp, folder, measurement_type, f_min, f_max, set_range_f = False, save_data = False, threshold = 0.1):
    'This function plots all fine laser scans within the given timestamp range'
    T = []
    LW1 = []
    LW2 = []

    for i in range(len(timestamp)):
        print i
        time = timestamp[i]
        data_folder, f,y, name = load_fine_lr_scan(folder = folder, timestamp = time, return_folder = True)

        if len(f) == 0:
            continue

        T.append(1)
        f = f*30./3.
        
        ''' Assigns f_min and f_max '''
        if set_range_f == True:
            ind_min = np.where(f < f_min)
            ind_max = np.where(f > f_max)
            f = np.delete(f, ind_max[0])
            f = np.delete(f, ind_min[0])
            y = np.delete(y, ind_max[0])
            y = np.delete(y, ind_min[0])
        else:
            f_max = max(f)
            f_min = min(f)
        print threshold
        max_index, max_value = max(enumerate(y), key = operator.itemgetter(1))
        if max_value < threshold:
            print "removing plot, max value < threshold"
            continue

        offset = 0.00
        amplitude1 = max_value
        x01 = f[max_index]
        sigma1 = 5

        amplitude2 = max_value 
        x02 = f[max_index]+10
        sigma2 = 2

        """ Fit a Gaussian """
        fixed = [0]
        p0, fitfunc, fitfunc_str = common.fit_2gauss(offset, amplitude1, x01, sigma1, amplitude2, x02, sigma2)
        fit_result = fit.fit1d(f,y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)


        A1 = fit_result['params_dict']['A1']
        a = offset
        x01 = fit_result['params_dict']['x01']
        sigma1 = fit_result['params_dict']['sigma1']
        A2 = fit_result['params_dict']['A2']
        x02 = fit_result['params_dict']['x02']
        sigma2 = fit_result['params_dict']['sigma2']

        # f2 = np.linspace(min(f),max(f),1001)
        function1 = a + A1*np.exp(-(f-x01)**2/(2*sigma1**2))
        function2 = a + A2*np.exp(-(f-x02)**2/(2*sigma2**2))

        if A1 > 5  *max_value:
            continue

        Half_Max1 = np.where(function1 > A1/2)
        Half_Max_l1 = Half_Max1[0][0]
        Half_Max_r1 = Half_Max1[0][-1]
        # print Half_Max,
        # print Half_Max_l, Half_Max_r
        Linewidth1 = abs(f[Half_Max_r1]-f[Half_Max_l1])
        LW1.append(Linewidth1)

        Half_Max2 = np.where(function2 > A2/2)
        Half_Max_l2 = Half_Max2[0][0]
        Half_Max_r2 = Half_Max2[0][-1]
        # print Half_Max,
        # print Half_Max_l, Half_Max_r
        Linewidth2 = abs(f[Half_Max_r2]-f[Half_Max_l2])
        LW2.append(Linewidth2)
        
        ''' Plot figure '''
        fig,ax = plt.subplots(figsize=(6,4.7))
        ax.plot(f,y, 'bo', label = 'data')
        plot.plot_fit1d(fit_result, np.linspace(f[0],f[-1],len(f)), ax=ax,label='Fit',show_guess=True, plot_data=False)
        ax.plot(f,function1, 'g', label='testplot')
        ax.plot(f,function2,'b')
        #ax.legend()
        plt.plot(f[Half_Max_r1], A1/2+offset, 'go')
        plt.plot(f[Half_Max_l1], A1/2+offset, 'go')
        plt.plot(f[Half_Max_r2], A2/2+offset, 'ro')
        plt.plot(f[Half_Max_l2], A2/2+offset, 'ro')
        plt.xlabel('frequency [GHz]')
        plt.ylabel('Transmission')
        plt.xlim(f_min,f_max)
        plt.title(date + ' - ' + time + '\n' + measurement_type + '\n' + 'Linewidth is ' + str(round(Linewidth1,2)) + 'GHz') #'GHz and \n Linewidth2 is ' + str(round(Linewidth2,2)) + 'GHz')

        # if save_data == True:
        #   fig.savefig(newpath +'/' + "laserfit " + time + ".png")

    if len(T) ==0:
        print "No fine laser data with these timestamps"




#! Not tested yet
def bla_fine_laser_plots(date, timestamp, folder, measurement_type, f_min, f_max, set_range_f, save_data):
    'This function plots all fine laser scans within the given timestamp range'
    T = []
    LW = []

    for i in range(len(timestamp)):
        time = timestamp[i]
        f,y, name = load_fine_lr_scan(folder = folder, timestamp = time)

        if len(f) == 0:
            continue

        T.append(1)
        f = f*35./3.
        
        ''' Assigns f_min and f_max '''
        if set_range_f == True:
            ind_min = np.where(f < f_min)
            ind_max = np.where(f > f_max)
            f = np.delete(f, ind_max[0])
            f = np.delete(f, ind_min[0])
            y = np.delete(y, ind_max[0])
            y = np.delete(y, ind_min[0])
        else:
            f_max = max(f)
            f_min = min(f)

        threshold = 0.02
        peak = np.where(y >= threshold)

        # print peak[0]
        # print range(peak[0])

        X2 = [f[index] for index in range(peak[0])]
        X3 = np.delete(f,X2)
        Y3 = np.delete(y,X2)


        fig,ax = plt.subplots(figsize=(6,4.7))
        ax.plot(X3,Y3)
        plt.xlabel('Relative frequency [GHz]')
        plt.ylabel('Transmission')
        ax.set_xlim(f_min,f_max)
        # plt.title(date + ' - ' + time + '\n' + measurement_type + '\n Fine laser plot without fit')
        plt.title(name)
        plt.show()
        
    if len(T) == 0:
        print "No fine laser data with these timestamps"


def fit_piezo_plots(date, timestamp, folder, V_min, V_max, set_range_V = False, save_data = False, threshold = 0.1,show_plots=True):
    'This function plots all length scans within the given timestamp range'
    print '### YOU\'RE ANALYSING PIEZO SCANS ### \n'
    N_points = 140 #the number of points around the resonant centre cavity length that are averaged for all scans 
    conversion_factor = 182 #nm / V . This is the calculated conversion factor from piezo voltage to distance moved
    LW = []
    Quality = []
    t = []
    peak_positions = []
    times = []
    frequency = []
    standard_dev = []
    T = np.array([])
    Ts = np.zeros(2*N_points)

    for i in range(len(timestamp)):
        time = timestamp[i]
        data_folder, V,y, name = load_piezo_scan(folder = folder, timestamp = time, return_folder = True)
        # V,y = zip(*sorted(zip(V, y)))

        if len(V) == 0:
            continue

        T = np.append(T,1)
        '''
        If a certain range for V is selected, delete all the datapoints outside this range
        '''
        if set_range_V == True:
            ind_min = np.where(V < V_min)
            ind_max = np.where(V > V_max)
            # print len(V)
            V = np.delete(V, ind_max[0])
            V = np.delete(V, ind_min)
            y = np.delete(y, ind_max[0])
            y = np.delete(y, ind_min)
        else:
            V_max = max(V)
            V_min = min(V)

        max_index, max_value = max(enumerate(y), key = operator.itemgetter(1))
        if max_value < threshold:
            print "removing plot, max value < threshold"
            continue

        offset = 0.01
        amplitude = 2*max_value
        x0 = V[max_index]
        gamma = 0.002



        """ Fit a Lorentzian """
        fixed = []
        p0, fitfunc, fitfunc_str = common.fit_lorentz(offset, amplitude, x0, gamma)
        fit_result = fit.fit1d(V,y, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)


        A = fit_result['params_dict']['A']
        a = offset
        x0 = fit_result['params_dict']['x0']
        gamma = fit_result['params_dict']['gamma']
        u_gamma = fit_result['error_dict']['gamma']

        # f2 = np.linspace(min(f),max(f),1001)
        function = a + 2*A/np.pi*gamma/(4*(V-x0)**2+gamma**2)
        if A > 5  *max_value:
            continue

        Linewidth = gamma#
        u_Linewidth = u_gamma
        LW.append(Linewidth)
        peak_positions = np.append(peak_positions, x0)
        times = np.append(times,time)
        
        dLs, V_zoom, y_zoom, zoom_success = zoom_around_peak(V, y, x0, N_points, conversion_factor = conversion_factor)

        if zoom_success: #if the ranging was not succesful, y_plot does not have the right dimension to fit on Ts, and should not be taken into account
            final_dLs = dLs
            Ts = Ts + y_zoom
        
        if (show_plots and zoom_success):
            plot_current_and_average_transmission(dLs,Ts,T,y_zoom)

        if show_plots:
            ''' Plot figure '''
            fig,ax = plt.subplots(figsize=(6,4))
            plot.plot_fit1d(fit_result, np.linspace(V_zoom[0],V_zoom[-1],len(V_zoom)*10 ), data_linestyle = '-',data_color = 'navy', data_lw = 2,ax=ax,color='r',show_guess=True, plot_data=True, label = 'fit', add_txt=False)
            ax.legend()
            ax.set_xlabel('Voltage [V]',fontsize=14)
            ax.set_ylabel('Transmission (arb. units)',fontsize=14)
            ax.tick_params(axis = 'both', which = 'major',labelsize = 14)
            ax.set_xlim(V_zoom[0],V_zoom[-1])
            ax.set_title(data_folder + ' \n Cavity linewidth in length is %s $\pm$ %s nm' %(round(Linewidth*conversion_factor,2), round(u_Linewidth*conversion_factor,3)))

            if save_data == True:
                fig.savefig(data_folder +'/'+ "piezofit " + ".png")
            plt.show()
            plt.close()

    if len(T) ==0:
        print "No piezo data with these timestamps"

    average_LW = sum(LW)/len(LW)
    variance_LW = sum([(x-average_LW)**2 for x in LW])/len(LW)
    st_dev_LW = np.sqrt(variance_LW)
    print  10*'*' + 'Average linewidth is ', round(average_LW*conversion_factor,3), '+-', round(st_dev_LW*conversion_factor,3), ' nm '+ 10*'*'


    #fit the averaged data -> if the linewidth is worse, the averaging is not done well.

    offset = 0.05
    amplitude = 2*max_value
    x0 = 0
    gamma = 0.3

    fixed = []
    p0, fitfunc, fitfunc_str = common.fit_lorentz(offset, amplitude, x0, gamma)
    fit_result = fit.fit1d(final_dLs,Ts/len(T), None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)

    A = fit_result['params_dict']['A']
    a = offset
    x0 = fit_result['params_dict']['x0']
    gamma = fit_result['params_dict']['gamma']
    u_gamma = fit_result['error_dict']['gamma']

    fig,ax = plt.subplots(figsize=(6,4))
    plot.plot_fit1d(fit_result, np.linspace(final_dLs[0],final_dLs[-1],10*len(final_dLs)), ax=ax,color='navy',show_guess=True, plot_data=True, data_color = 'darkorange',data_linestyle ='-',data_lw  =3,lw = 1.5,label = 'fit', add_txt=False)
    ax.legend()
    ax.set_xlabel('detuning in length (nm)',fontsize=14)
    ax.set_ylabel('Transmission (a.u.)',fontsize=14)
    ax.tick_params(axis = 'both', which = 'major',labelsize = 14)
    ax.set_xlim(final_dLs[0],final_dLs[-1])
    ax.set_title(data_folder + ' \n Cavity linewidth in length is %s $\pm$ %s nm' %(round(gamma,3), round(u_gamma,3)))
    fig.savefig(data_folder +'/'+ "average_transmission " + ".png")
    fig.savefig(data_folder +'/'+ "average_transmission " + ".eps")
    fig.savefig(data_folder +'/'+ "average_transmission " + ".pdf")
    plt.tight_layout()

    plt.show()
    plt.close()


    determine_drift(times,peak_positions, data_folder)

def fit_birefringence(date, timestamp, folder, V_min=0, V_max=1, set_range_V = False, save_data = False, threshold = 0.1,show_plots=True, sweep_polarization=False):
    """
    Function to fit linewidth data with a splitted resonance due to birefringence.

    """
    print '### YOU\'RE ANALYSING PIEZO SCANS ### \n'
    LW1 = []
    LW2 = []
    As1 = []
    As2 = []
    t = []
    T = []
    peak_positions1 = []
    peak_positions2 = []
    times = []
    frequency = []
    standard_dev = []
    polarizations = []
    N_points=150
    conversion_factor=182
    Ts = np.zeros(2*N_points)

    for i in range(len(timestamp)):
        time = timestamp[i]
        data_folder, V,y, name = load_piezo_scan(folder = folder, timestamp = time, return_folder = True)

        if sweep_polarization:
            polarization = name[-3:]

        if len(V) == 0:
            print 'your voltage array was empty!!'
            continue

        '''
        If a certain range for V is selected, delete all the datapoints outside this range
        '''
        if set_range_V == True:
            ind_min = np.where(V < V_min)
            ind_max = np.where(V > V_max)
            # print len(V)
            V = np.delete(V, ind_max[0])
            V = np.delete(V, ind_min)
            y = np.delete(y, ind_max[0])
            y = np.delete(y, ind_min)
        else:
            V_max = max(V)
            V_min = min(V)

        max_index, max_value = max(enumerate(y), key = operator.itemgetter(1))
        if max_value < threshold:
            print "removing plot, max value < threshold"
            continue

        offset = 0.003
        amplitude = 2*max_value
        x0 = V[max_index]
        gamma = 0.001
        amplitude2 = 0.5*amplitude
        gamma2 = gamma
        splitting = 0.004

        #Fit a Lorentzian with splitting fixed. 
        fixed = [5]

        fit_result, left_peak_high, discard_both =fit_two_lorentzian_peaks(V,y,offset,amplitude,x0,gamma,amplitude2,gamma2,splitting,fixed)

        if discard_both:
            print 'discarding both Lorentzian fits since they do not give a positive transmission' 
            continue


        if left_peak_high:
            A1 = fit_result['params_dict']['A1']
            x01 = fit_result['params_dict']['x01']
            gamma1 = fit_result['params_dict']['gamma1']
            u_gamma1 = fit_result['error_dict']['gamma1']
            A2 = fit_result['params_dict']['A2']
            if len(fixed) == 0:
                x02 = fit_result['params_dict']['x01'] + fit_result['params_dict']['dx2']
            else:
                x02 = fit_result['params_dict']['x01'] + splitting
            gamma2 = fit_result['params_dict']['gamma2']
            u_gamma2 = fit_result['error_dict']['gamma2']

        else:
            # interchange 1 & 2 compared to above: the left peak is now corresponding to A1, x01, gamma1 etcetera
            A2 = fit_result['params_dict']['A1']
            x02 = fit_result['params_dict']['x01']
            gamma2 = fit_result['params_dict']['gamma1']
            u_gamma2 = fit_result['error_dict']['gamma1']
            A1 = fit_result['params_dict']['A2']
            if len(fixed) == 0:
                x01 = fit_result['params_dict']['x01'] + fit_result['params_dict']['dx2']
            else:
                x01 = fit_result['params_dict']['x01'] - splitting
            gamma1 = fit_result['params_dict']['gamma2']
            u_gamma1 = fit_result['error_dict']['gamma2']

        a1 = fit_result['params_dict']['a1']

        LW1.append(np.abs(gamma1))
        LW2.append(np.abs(gamma2))
        # multiply A1 and A2 by 1000. We get a.u. anyway for now.
        As1.append(np.abs(A1*1000))
        As2.append(np.abs(A2*1000))
        peak_positions1 = np.append(peak_positions1, x01)
        peak_positions2 = np.append(peak_positions2, x02)
        times = np.append(times,time)

        if sweep_polarization:
            polarizations = np.append(polarizations,int(polarization))

        dLs, V_zoom, y_zoom, zoom_success = zoom_around_peak(V, y, x02, N_points, conversion_factor = conversion_factor)

        if zoom_success: #if the ranging was not succesful, y_plot does not have the right dimension to fit on Ts, and should not be taken into account
            final_dLs = dLs
            Ts = Ts + y_zoom
        
        if show_plots and zoom_success:
            plot_current_and_average_transmission(dLs,Ts,T,y_zoom)

        if show_plots:
            ''' Plot figure '''
            fig,ax = plt.subplots(figsize=(6,4.7))
            #ax.plot(V,y, 'b', label = 'data')
            plot.plot_fit1d(fit_result, np.linspace(V_zoom[0],V_zoom[-1],10*len(V_zoom)), ax=ax,color='r',show_guess=True, plot_data=True, label = 'fit', add_txt=False)
            ax.legend()
            ax.set_xlabel('Voltage [V]',fontsize=14)
            ax.set_ylabel('Transmission (arb. units)',fontsize=14)
            ax.tick_params(axis = 'both', which = 'major',labelsize = 14)
            ax.set_xlim(V_zoom[0],V_zoom[-1])
            ax.set_title(data_folder + ' \n Cavity linewidth in length is %s $\pm$ %s nm (%s $\pm$ %s V)' %(round(gamma1*182,2), round(u_gamma1*182,3),round(gamma1,4), round(u_gamma1,4)))

            if save_data == True:
                fig.savefig(data_folder +'/'+ "piezofit " + ".png")
            plt.show()
            plt.close()


    peak_differences = np.abs(peak_positions2 - peak_positions1)
    peak_differences = (peak_differences[np.where(peak_differences>3.e-3)])
    average_peak_diff, variance_peak_diff, st_dev_peak_diff = calculate_stats(peak_differences)

    average_LW1, variance_LW1, st_dev_LW1 = calculate_stats(LW1)
    average_LW2, variance_LW2, st_dev_LW2 = calculate_stats(LW2)

    average_A1, variance_A1, st_dev_A1 = calculate_stats(As1)
    average_A2, variance_A2, st_dev_A2 = calculate_stats(As2)

    print  10*'*' + 'Average peak difference is ', round(average_peak_diff,4), '+-', round(st_dev_peak_diff,4), ' V '+ 10*'*'
    print  10*'*' + 'Average linewidth of peak 1 is ', round(average_LW1*182,2), '+-', round(st_dev_LW1*182,2), ' nm  (' , round(average_LW1,2), '+-', round(st_dev_LW1,2), ' V )'+ 10*'*'
    print  10*'*' + 'Average linewidth of peak 2 is ', round(average_LW2*182,2), '+-', round(st_dev_LW2*182,2), ' nm  (', round(average_LW2,4), '+-', round(st_dev_LW2,4), ' V )'+10*'*'



    #fit the averaged data -> if the linewidth is worse, the averaging is not done well.

    offset = 0.003
    amplitude = 2*max_value
    x0 = 0
    gamma = average_LW1*conversion_factor
    splitting = 0.002*conversion_factor
    amplitude2 = amplitude
    gamma2 = average_LW2*conversion_factor
    #Fit a Lorentzian with splitting fixed. 
    fixed = []

    #p0, fitfunc, fitfunc_str = common.fit_lorentz(offset, amplitude, x0, gamma)
    fit_result,  left_peak_high, fit_failed = fit_two_lorentzian_peaks(final_dLs,Ts,offset,amplitude,x0,gamma,amplitude2,gamma2,splitting,fixed)

    a = offset
    x01 = fit_result['params_dict']['x01']
    gamma1 = fit_result['params_dict']['gamma1']
    u_gamma1 = fit_result['error_dict']['gamma1']
    gamma2 = fit_result['params_dict']['gamma2']
    u_gamma2 = fit_result['error_dict']['gamma2']

    fig,ax = plt.subplots(figsize=(6,4.5))
    plt.tight_layout()

    plot.plot_fit1d(fit_result, np.linspace(final_dLs[0],final_dLs[-1],10*len(finaL_dLs)), ax=ax,color='navy',show_guess=True, plot_data=True, data_color = 'darkorange',data_linestyle ='-',data_lw  =3,lw = 1.5,label = 'fit', add_txt=False)
    ax.legend()
    ax.set_xlabel('detuning in length (nm)',fontsize=14)
    ax.set_ylabel('Transmission (a.u.)',fontsize=14)
    ax.tick_params(axis = 'both', which = 'major',labelsize = 14)
    ax.set_xlim(dLs[0],dLs[-1])
    ax.set_title(data_folder + ' \n Cavity linewidth of peak 1 in length is %s $\pm$ %s nm. \n Cavity linewidth of peak 2 in length is %s $\pm$ %s nm.' %(round(gamma1,3), round(u_gamma1,3), round(gamma2,3), round(u_gamma2,3)))
    fig.savefig(data_folder +'/'+ "average_transmission " + ".png")
    fig.savefig(data_folder +'/'+ "average_transmission " + ".eps")
    fig.savefig(data_folder +'/'+ "average_transmission " + ".pdf")
    plt.show()
    plt.close()

    determine_drift(times,peak_positions1, data_folder)
    determine_drift(times,peak_positions2, data_folder)

    # print polarizations
    # print As1
    # print As2
    if sweep_polarization:
        plot_peak_intensities(As1,As2,polarizations,data_folder,save_data=save_data)


#######################################################
#some useful functions
#######################################################


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def fit_two_lorentzian_peaks(V,y,offset,amplitude,x0,gamma,amplitude2,gamma2,splitting,fixed):
    """
    Function that fits two lorentzian peaks to noisy data. 
    Input:
    V - the x array of data 
    y - the y array of data 
    offset - guess for the fit offset
    amplitude - guess for the amplitude of the highest peak
    x0 - guess for the location of the highest peak
    gamma - guess for the FWHM of the highest peak
    amplitude2 - guess for the amplitude of the lowest peak
    gamma2 - guess for the FWHM of the lowest peak
    splitting - guess for the splitting between the two peaks
    fixed - array of the parameters to be fixed in fit function. (keep order; e.g.to fix spliting: fixed = [5])
    """
    #First try fitting with the left peak the largest
    dx2_1 = splitting
    p0_1, fitfunc_1, fitfunc_str_1 = common.fit_2lorentz_splitting(offset, amplitude, x0, gamma, amplitude2, dx2_1, gamma2)
    fit_result_1 = fit.fit1d(V,y, None, p0=p0_1, fitfunc=fitfunc_1, do_print=False, ret=True,fixed=fixed)
    #Then try fitting with the right peak the largest
    dx2_2 = - splitting
    p0_2, fitfunc_2, fitfunc_str_2 = common.fit_2lorentz_splitting(offset, amplitude, x0, gamma, amplitude2, dx2_2, gamma2)
    fit_result_2 = fit.fit1d(V,y, None, p0=p0_2, fitfunc=fitfunc_2, do_print=False, ret=True,fixed=fixed)

    discard_1 = False
    discard_2 = False
    if fit_result_1['params_dict']['A2']*fit_result_1['params_dict']['gamma2'] < 0:
        # print 'discard first fit' 
        discard_1 = True
    if fit_result_2['params_dict']['A2']*fit_result_2['params_dict']['gamma2'] < 0:
        # print 'discard second fit' 
        discard_2 = True

    if discard_1 and discard_2:
        print 'discarding both Lorentzian fits since they do not give a positive transmission' 
        fit_failed = True
    else:
        fit_failed = False

    if (fit_result_1['chisq'] < fit_result_2['chisq']) or discard_2:
        fit_result = fit_result_1
        left_peak_high = True

    else:
        fit_result = fit_result_2
        left_peak_high = False

    return fit_result, left_peak_high, fit_failed

def zoom_around_peak(V, y, x0, N_points, conversion_factor = 182):
    """function that calculates the new arrays of data, zoomed in on a range around the resonance peak.
    input:
    V - the old voltage range
    y - the old signal data
    x0 - the fitted value of the peak (in V)
    N_points - the number of points to take into account around the peak
    conversion_factor - the piezo conversion factor from voltage to length

    Returns:
    dLs - the new range expressed in detuning in length from the peak position
    V_zoom - the new voltage range 
    y_zoom - the photodiode signal in the new zoomed in range
    success - False if the zoom in failed 
    """
    V_near, idx_near  = find_nearest(V,x0)
    idx_min = idx_near-N_points
    idx_max = idx_near+N_points

    if (idx_min < 0) or (idx_max > len(V)):        
        print("The linewidth is too far to the side of the plot to add it to the averaged plot")
        success = False
        V_zoom = V
        y_zoom = y
        dLs = V
    else:
        success = True
        V_zoom = V[idx_min:idx_max]
        y_zoom = y[idx_min:idx_max]

        dV = np.abs(V_zoom[0] - V_zoom[-1])
        dL = dV * conversion_factor
        dLs = np.linspace(-dL/2,dL/2,N_points*2)
    return dLs, V_zoom, y_zoom, success

def plot_current_and_average_transmission(dLs,Ts,T,y_zoom):
    fig,ax = plt.subplots(figsize=(6,4))
    ax.plot(dLs,Ts/len(T), '-',color = 'dimgray', linewidth = 2, label = 'average')
    ax.plot(dLs,y_zoom,'-',color = 'navy', linewidth = 2)
    ax.legend()
    ax.set_xlabel('detuning in length (nm)',fontsize=14)
    ax.set_xlim(dLs[0],dLs[-1])
    ax.set_ylabel('Transmission (a.u.)',fontsize=14)
    plt.show()
    plt.close()  

def determine_drift(times, peak_positions, data_folder):
    times_dec = []

    ftr = [3600,60,1]
    for time in times:
        time_dec = sum([a*b for a,b in zip(ftr, [ int(time[0:2]),int(time[2:4]),int(time[4:6]) ]) ])
        times_dec = np.append(times_dec, time_dec)

    times_dec = times_dec - times_dec[0]

    g_a = 2.13
    g_b = 0.0007
    fixed=[]
    p0, fitfunc, fitfunc_str = common.fit_line(g_a, g_b)
    fit_result = fit.fit1d(times_dec,peak_positions, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)

    b = fit_result['params_dict']['b']
    u_b = fit_result['error_dict']['b']

    fig,ax = plt.subplots(figsize=(6,4.7))

    ax = plot.plot_fit1d(fit_result, np.linspace(times_dec[0],times_dec[-1],10*len(times_dec)), ax=ax,color='r',show_guess=True, plot_data=True, label = 'fit', add_txt=False, ret='ax')
    ax.set_xlabel('time (s)')
    ax.set_ylabel('peak position (V)')
    ax.set_title(data_folder + ' \n Drift is {} $\pm$ {} mV/s '.format(round(b*1.e3,3), round(u_b*1.e3,3)))

    fig.savefig(data_folder +'/'+ "drift "  + ".png")

    plt.show()
    plt.close()

def calculate_stats(LW):
    average_LW = sum(LW)/len(LW)
    variance_LW = sum([(x-average_LW)**2 for x in LW])/len(LW)
    st_dev_LW = np.sqrt(variance_LW) 

    return average_LW, variance_LW, st_dev_LW

def plot_peak_intensities(As1,As2,polarizations,data_folder,save_data=True):
    p_old = 0
    p_indices = []
    avg_As1 = []
    avg_As2 = []
    u_avg_As1 = []
    u_avg_As2 = []
    j_old = 0
    index_old = 0

    for i,p in enumerate(polarizations):
        if p != p_old:
            p_indices = np.append(p_indices,int(i))
        p_old = p

    polarizations_values = np.array([polarizations[i] for i in p_indices])
    polarizations_min = np.min(polarizations_values)
    polarizations_max = np.max(polarizations_values)
    p_indices = np.append(p_indices,len(polarizations))

    for j,index in enumerate(p_indices):
        if j == 0:
            continue

        As1_j = As1[index_old:int(index)]
        As2_j = As2[index_old:int(index)]

        Average_As1_j, variance_As1_j, st_dev_As1_j = calculate_stats(As1_j)
        avg_As1 = np.append(avg_As1, Average_As1_j)
        u_avg_As1 = np.append(u_avg_As1, st_dev_As1_j)

        Average_As2_j, variance_As2_j, st_dev_As2_j = calculate_stats(As2_j)
        avg_As2 = np.append(avg_As2, Average_As2_j)
        u_avg_As2 = np.append(u_avg_As2, st_dev_As2_j)

        index_old = int(index)

    # print avg_As1
    # print avg_As2
    # print u_avg_As1
    # print u_avg_As2

    rel_As1 = np.divide(avg_As1,(avg_As1+avg_As2))
    rel_As2 = np.divide(avg_As2,(avg_As1+avg_As2))
    u_rel_As1 = np.sqrt( u_avg_As1**2 * ( 1/(avg_As1+avg_As2) - ( avg_As1/(avg_As1+avg_As2)**2) )**2 + u_avg_As2**2 * ( avg_As1/(avg_As1+avg_As2)**2 ) **2  )
    u_rel_As2 = np.sqrt( u_avg_As2**2 * ( avg_As2/(avg_As1+avg_As2)**2) **2 + u_avg_As2**2 * ( 1/(avg_As1+avg_As2) - ( avg_As2/(avg_As1+avg_As2)**2) **2 ) )
    # print rel_As1
    # print rel_As2
    # print u_rel_As1
    # print u_rel_As2


    # fitting the average peak intensities with cosines
    g_f =1/90.
    g_a = 1
    g_A = 1
    g_phi = 0.
    fixed =[]
    p0, fitfunc, fitfunc_str = common.fit_cos(g_f, g_a, g_A, g_phi)
    g_f = 1/90.
    g_a = 0.5
    g_A = 0.2
    g_phi = 0.
    fit_result1 = fit.fit1d(polarizations_values,avg_As1, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)
    p0, fitfunc, fitfunc_str = common.fit_cos(g_f, g_a, g_A, g_phi)
    fit_result2 = fit.fit1d(polarizations_values,avg_As2, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True,fixed=fixed)



    # Plot fits
    fig,ax3 = plt.subplots(figsize=(6,4.7))
    # with plot_fit1d, I cannot plot errorbars unless I add the value manually to the result dictionary
    fit_result1['yerr'] = u_avg_As1
    fit_result2['yerr'] = u_avg_As2
    plot.plot_fit1d(fit_result1, np.linspace(polarizations_min,polarizations_max,10*len(polarizations_values)), ax=ax3,color='c',show_guess=True, plot_data=True, label = 'fit peak 1', add_txt=False)
    plot.plot_fit1d(fit_result2, np.linspace(polarizations_min,polarizations_max,10*len(polarizations)), ax=ax3,color='r',show_guess=True, plot_data=True, label = 'fit peak 2', add_txt=False)
    
    #fitted period of the oscilation - average of both
    if len(fixed) != 0:
        period = 1/g_f
        u_period = 0.
    else:
        period = 2./(fit_result1['params_dict']['f']+fit_result2['params_dict']['f'])
        u_period = np.sqrt( (2/(fit_result1['params_dict']['f']+fit_result2['params_dict']['f'])**4)  * ( fit_result1['error_dict']['f']**2 + fit_result2['error_dict']['f']**2 ) )
    dphi = (fit_result1['params_dict']['phi']-fit_result2['params_dict']['phi'])
    u_dphi = np.sqrt(fit_result1['error_dict']['phi']**2 + fit_result2['error_dict']['phi']**2 )
    
    ax3.legend()
    ax3.set_xlabel('polarization [deg]',fontsize=14)
    ax3.set_ylabel('Peak intensities (a. u.)',fontsize=14)
    ax3.tick_params(axis = 'both', which = 'major',labelsize = 14)
    ax3.set_xlim(polarizations_min,polarizations_max)
    ax3.set_title(data_folder + ' \n period = {} +- {} deg, dphi = {} +- {} deg '.format(round(period,1),round(u_period,1),round(dphi,0),round(u_dphi,0)))

    if save_data:
        fig.savefig(data_folder +'/'+ "birefringence_fit " + ".png")
    plt.show()
    plt.close()




    g_f =1/90.
    g_a = 0.6
    g_A = 0.4
    g_phi = 0.
    fixed =[]
    p0, fitfunc, fitfunc_str = common.fit_cos(g_f, g_a, g_A, g_phi)
    fit_result3 = fit.fit1d(polarizations_values,rel_As1, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
    g_f = 1/90.
    g_a = 0.4
    g_A = 0.4
    g_phi = 0.
    p0, fitfunc, fitfunc_str = common.fit_cos(g_f, g_a, g_A, g_phi)
    fit_result4 = fit.fit1d(polarizations_values,rel_As2, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)

    # Plot fits
    fig,ax = plt.subplots(figsize=(6,4.7))
    # with plot_fit1d, I cannot plot errorbarsm unless I add the value manually to the result dictionary
    fit_result3['yerr'] = u_rel_As1
    fit_result4['yerr'] = u_rel_As2
    plot.plot_fit1d(fit_result3, np.linspace(polarizations_min,polarizations_max,10*len(polarizations_values)), ax=ax,color='c',show_guess=True, plot_data=True, label = 'peak 1', add_txt=False)
    plot.plot_fit1d(fit_result4, np.linspace(polarizations_min,polarizations_max,10*len(polarizations)), ax=ax,color='r',show_guess=True, plot_data=True, label = 'peak 2', add_txt=False)
    
    #average period of the oscilation
    if len(fixed) != 0:
        period = 1/g_f
        u_period = 0.
    else:
        period = 2./(fit_result3['params_dict']['f']+fit_result4['params_dict']['f'])
        u_period = np.sqrt( (2/(fit_result3['params_dict']['f']+fit_result4['params_dict']['f'])**4)  * ( fit_result3['error_dict']['f']**2 + fit_result4['error_dict']['f']**2 ) )
    dphi = (fit_result3['params_dict']['phi']-fit_result4['params_dict']['phi'])
    u_dphi = np.sqrt(fit_result3['error_dict']['phi']**2 + fit_result4['error_dict']['phi']**2 )
    
    ax.legend()
    ax.set_xlabel('polarization [deg]',fontsize=14)
    ax.set_ylabel('Relative peak intensities',fontsize=14)
    ax.tick_params(axis = 'both', which = 'major',labelsize = 14)
    ax.set_xlim(polarizations_min,polarizations_max)
    ax.set_title(data_folder + ' \n period = {} +- {} deg, dphi = {} +- {} deg '.format(round(period,1),round(u_period,1),round(dphi,0),round(u_dphi,0)))

    if save_data == True:
        fig.savefig(data_folder +'/'+ "relative_birefringence_fit " + ".png")
    plt.show()
    plt.close()