"""
This script uses a 2d fit (in fitting.fit) to fit 2D gaussians to a 2d counts scan
From this it estimates where in the scan unique NV's are located 
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import msvcrt

from analysis.lib.fitting import fit
from analysis.lib.fitting import common

import analysis.scripts.Fabrication.Display_scan2d as ds; reload(ds)


def plot_fit2d(fitres,meshx,meshy):
    """
    General function that plots data from a 2d fit
    """
    extent = [meshx[0,0],meshx[-1,-1],meshy[0,0],meshy[-1,-1]] #we need this for plotting
    fig,ax1 = plt.subplots()
    ax1.imshow(fitres['fitfunc'](meshx,meshy),extent = extent,origin='lower')
    plt.show()


def get_fit(data,meshx,meshy,g_offset,g_A,g_x0,g_y0,g_sigmax,g_sigmay,g_theta,**kw):
    """
    fit a single 2d scan with a general 2d gaussian
    """
    do_plot = kw.pop('do_plot',False)
    
    fitres = fit.fit2d((meshx,meshy),data,common.fit_2d_gaussian, g_offset,g_A, g_x0, g_y0, g_sigmax, g_sigmay,g_theta, ret = True )
    if type(fitres)!=dict:  #the fit failed; no NV
        return 0,0

    sigmax  = fitres['params_dict']['sigmax']
    sigmay  = fitres['params_dict']['sigmay']
    A = fitres['params_dict']['A']

    #rejection conditions
    ellipticity_condition = ((sigmax < sigmay*0.8) or (sigmax > sigmay*1.2)) #NV cannot be too elliptic
    amplitude_condition = A<(0.2*g_offset) #amplitude cannot be too small
    size_condition = ((sigmax>0.5) or (sigmay>0.5)) #NV cannot have a diameter larger than 0.5um.
    conditions = (ellipticity_condition or amplitude_condition) or size_condition
    
    if conditions:
        return 0,0      
    # print fitres['params_dict']
    
    if do_plot:
        plot_fit2d(fitres,meshx,meshy)
    
    return 1, fitres


def get_circular_fit(data,meshx,meshy,g_offset,g_A,g_x0,g_y0,g_sigma,**kw):
    """
    fit a single 2d scan with a circular 2d gaussian
    """
    do_plot = kw.pop('do_plot',False)

    fitres = fit.fit2d((meshx,meshy),data,common.fit_2d_gaussian_circular, g_offset,g_A, g_x0, g_y0, g_sigma, ret = True )
    if type(fitres)!=dict: #the fit failed; no NV
        return 0,0

    sigma  = fitres['params_dict']['sigma']
    A = fitres['params_dict']['A']

    #rejection conditions
    amplitude_condition = A<(0.2*g_offset) #amplitude cannot be too small
    size_condition = sigma>0.5 #NV cannot have a diameter larger than 0.5um.
    conditions = amplitude_condition or size_condition
    
    if conditions: 
        return 0,0

    if do_plot:
        plot_fit2d(fitres,meshx,meshy)
    
    return 1, fitres



def crop_array_around_value(array,value,patchsize):
    """
    Generally useful function that crops an array to 2*patchsize around value
    """
    i_closest = np.argmin(abs(array-value))
    if (i_closest-patchsize)>0:
        i_min = i_closest-patchsize
    else:
        i_min = 0
    if (i_closest-patchsize)<len(array):
        i_max = i_closest+patchsize
    else:
        i_max = -1
    
    cropped_array = array[i_min:i_max]
    return i_min,i_max,cropped_array


def get_data(folder):
    """
    Get the 2d scan data from a folder
    """
    dscan = ds.DisplayScan(folder = folder)
    dscan.get_data()

    data = dscan.countrates
    x = dscan.xvalues
    y = dscan.yvalues

    return x,y,data


def find_NV_locations(x,y,data,folder,g_sigma = 0.2,**kw):
    plot_NV_zoom=kw.pop('plot_NV_zoom',False)

    # guess parameters for the fit
    g_offset=np.average(data)
    g_A = data.max()
    g_sigmax = g_sigma
    g_sigmay = g_sigma
    g_theta=0

    #arrays of guesses for x0 and y0 - a fit attempt will be done for each
    xstepsize = g_sigmax*2 #the stepsize is 2*sigma, to ensure each NV is found
    ystepsize = g_sigmay*2
    g_x0s = np.arange(x[0]+g_sigmax,x[-1],xstepsize)
    g_y0s = np.arange(y[0]+g_sigmay,y[-1],ystepsize)

    #the fit is performed on a subset of the data around the intial guesses for x0 and y0
    patchsizex = int(len(x)/abs(max(x)-min(x))*g_sigmax*3) #number of indices around g_x0
    patchsizey = int(len(y)/abs(max(y)-min(y))*g_sigmay*3) #number of indices around g_x0
    # print patchsizex
    # print patchsizey

    #initialise some arrays or dictionaries to store results
    results=np.zeros((len(g_x0s),len(g_y0s)))
    x0s=np.array([0])
    y0s=np.array([0])
    NVs = {}
    NV_nr=0
    f=0

    total_nr_fit_attempts = len(g_x0s)*len(g_y0s)
    
    for i,g_x0 in enumerate(g_x0s):       
        ##################################################
        if (msvcrt.kbhit() and (msvcrt.getch() == 'q')):
            aborted=True
            break
        ##############################################
        i_min,i_max,x_around_x0 = crop_array_around_value(x,g_x0,patchsizex)
        
        for j,g_y0 in enumerate(g_y0s):
            ##################################################
            if (msvcrt.kbhit() and (msvcrt.getch() == 'q')):
                aborted=True
                break
            ##############################################
            j_min,j_max,y_around_y0 = crop_array_around_value(y,g_y0,patchsizey)
            #define meshgrid for x,y cropped around x0,y0
            meshx,meshy = np.meshgrid(x_around_x0,y_around_y0)
            data_x0y0 = data[j_min:j_max,i_min:i_max]
            
            #some progress indication for the impatient
            if f%100==0:
                print f,'out of ',total_nr_fit_attempts,'attempts to fit completed'
                print 'found', NV_nr, 'NVs so far'
            f+=1    
            
            #the actual fitting
            results[i,j],fitres = get_fit(data_x0y0,meshx,meshy,g_offset,g_A,g_x0,g_y0,g_sigmax,g_sigmay,g_theta)
            #results[i,j],fitres = get_circular_fit(data,meshx,meshy,g_offset,g_A,g_x0,g_y0,g_sigma) #general fit works better

            if results[i,j] == 1: #an NV was found
                #first check if it is a unique NV centre, or a previously found one.
                x0 = fitres['params_dict']['x0']
                min_dx0 = x0s[np.argmin(abs(x0s-x0))]-x0 #minimal distance to all previously found x0s

                y0 = fitres['params_dict']['y0']
                min_dy0 = y0s[np.argmin(abs(y0s-y0))]-y0 #minimal distance to all previously found y0s                  

                if ((abs(min_dx0)>g_sigmax) or (abs(min_dy0)>g_sigmay)): 
                    #if in either direction the NV is more than sigma away from a previously found one, it's unique
                    NV_nr += 1
                    'found NV! ', 'number', NV_nr
                    NVs['NV_'+str(NV_nr)] = fitres

                    if plot_NV_zoom:                      
                        extent = [meshx[0,0],meshx[-1,-1],meshy[0,0],meshy[-1,-1]]
                        fig,ax = plt.subplots()
                        ax.imshow(data_x0y0,extent = extent,origin='lower',interpolation='None')
                        plotname = 'NV_%d_x_%1.f_y_%1.f_y_'%(NV_nr,x0,y0)
                        title = folder+'\n'+plotname
                        ax.set_title(title)
                        try:
                            plt.savefig(os.path.join(folder,plotname+'.png'))
                        except:
                            print 'could not save fig'

                    # print 'x range',x_around_x0[0],g_x0,x_around_x0[-1] #for debugging
                    # print 'y range',y_around_y0[0],g_y0,y_around_y0[-1] #for debugging

                    #update x0s and y0s with our new NV
                    x0s=np.append(x0s,x0)
                    y0s=np.append(y0s,y0)

    return NVs

def plot_NV_locations(x,y,data,NVs):

    extent = (min(x),max(x), min(y),max(y))

    fig,ax = plt.subplots()
    ax.imshow(data, extent=extent,origin='lower',interpolation='None')

    for NV in NVs:
        fitres=NVs[NV]
        x0=fitres['params_dict']['x0']
        y0=fitres['params_dict']['y0']
        sigmax=fitres['params_dict']['sigmax']
        sigmay=fitres['params_dict']['sigmay']
        print 10*'*'
        print NV
        print 'x0 = %.1f um, y0 = %.1f um'%(x0,y0)
        print 'sigma_x = %.2f um, sigma_y = %.2f um'%(sigmax,sigmay)
        print 10*'*'
        ax.plot(x0,y0,'+w',mew=2,markersize=50)
        

    plt.show()  

    return ax

def zoom_around_NV(x,y,data,folder,NVs):
    for i in NVs:
        nr= i[3:]
        x0 = NVs[i]['params_dict']['x0']
        y0 = NVs[i]['params_dict']['y0']

        i_min,i_max,x_around_x0 = crop_array_around_value(x,x0,11)

        j_min,j_max,y_around_y0 = crop_array_around_value(y,y0,11)

        meshx,meshy = np.meshgrid(x_around_x0,y_around_y0)
        data_x0y0 = data[j_min:j_max,i_min:i_max]


        extent = [meshx[0,0],meshx[-1,-1],meshy[0,0],meshy[-1,-1]]
        fig,ax = plt.subplots()
        plotname = '_zoom_NV_%s_x_%.1f_y_%.1f'%(nr,x0,y0)
        ax.set_title(folder+'\n'+plotname)
        c=ax.imshow(data_x0y0,extent = extent,origin='lower',interpolation='None')
        plt.colorbar(c)
        plt.savefig(os.path.join(folder,plotname+'.png'))

