"""
analysis of DD files to extract T2 for a certain number of pulses
NK 2017
"""
import numpy as np
import os
from analysis.lib.tools import plot
from analysis.lib.tools import toolbox as tb
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt


###############################
"""     Helper functions    """
###############################

def get_tstamps_dsets(contains,older_than = '',newer_than = '',n_datasets = 5,subsets = 5,**kw):
    """
    expects multiple data sets (n_datasets)
    extracts the initial timestamp for all of these data sets.
    Each data set can have multiple subset folders (these are specified by subsets).
    Counting start from new and goes towards old.
    returns a list of folder_lists (divided into datasets and subsets)
    """


    i_dataset = 0
    f_list = []
    while tb.latest_data(contains=contains,older_than = older_than,newer_than = newer_than, folder=  None,
                raise_exc = False) != False and (i_dataset< n_datasets):

        subset_f_list = []
        for i_subset in np.array(range(subsets))+1:
            latest_tstamp, f_sub = tb.latest_data(contains=contains+str(i_subset),older_than = older_than,
                                            newer_than = newer_than, folder=   None,
                raise_exc = False,return_timestamp = True)

            subset_f_list.append(f_sub)
        older_than = latest_tstamp
        i_dataset +=1
        f_list.append(subset_f_list)

    return f_list

def compile_xy_values_of_datasets(f_list,ssro_tstamp = '112128',**kw):
    """
    besides returning all msmt values this function also returns 
    one data object in case one wants to review msmt params.
    """


    if ssro_tstamp == None: 
        ssro_calib_folder = tb.latest_data('SSROCalib')
    else:
        ssro_dstmp, ssro_tstmp = tb.verify_timestamp(ssro_tstamp)
        ssro_calib_folder = tb.latest_data(ssro_tstmp)

    x = []
    y = []
    y_u = []
    for f in f_list:
        a = mbi.MBIAnalysis(f)
        a.get_sweep_pts()
        a.get_readout_results(name = 'adwindata')
        a.get_electron_ROC(ssro_calib_folder)

        x.extend(a.sweep_pts.reshape(-1))
        y.extend(a.p0.reshape(-1))
        y_u.extend(a.u_p0.reshape(-1))

    return x,y,y_u,a

def sort_for_best_signal(x_list,y_list,y_u_list,**kw):
    """
    compares rows in the y matrix and returns only the best values
    all other values are returned as a sorted 1D list to plot as well
    uses some numpy and matrix transpose magic.
    """
    s = np.shape(x_list)[1] ### need length of the dataset in terms of sweep points
    best_y = np.amax(y_list.T,axis = 1)
    best_x = x_list.T[range(s),np.argmax(y_list.T,axis = 1)]
    best_y_u = y_u_list.T[range(s),np.argmax(y_list.T,axis = 1)]
    if kw.pop('print_best', False):
        print 'Best decoupling times:'
        print best_x

    return best_x,best_y,best_y_u



###############################
"""      Full analysis      """
###############################


def analyse_dataset(contains,**kw):

    f_list = get_tstamps_dsets(contains,**kw)
    x_list,y_list,y_u_list = [],[],[]

    print f_list

    for f_sub_list in f_list:
        x,y,y_u,a = compile_xy_values_of_datasets(f_sub_list,**kw)

        x_list.append(x); y_list.append(y); y_u_list.append(y_u)


    if kw.pop('plot_raw',True):
        ax = a.default_ax(figsize = (7,4))

        for x,y, y_u in zip(x_list,y_list,y_u_list):
            ax.errorbar(x,y,y_u,zorder = 0, fmt = 'o', color = '#d3d3d3')

    x_list,y_list,y_u_list = np.array(x_list),np.array(y_list),np.array(y_u_list) ## conver tto array

    best_x,best_y,best_y_u = sort_for_best_signal(x_list,y_list,y_u_list,**kw)
    ax.set_xlabel(a.sweep_name)
    ax.set_ylabel(r'$F(|0\rangle)$')
    ax.set_xlim([0,np.amax(best_x)*1.05])
   

    # To avoid dips in the analysis
    for jj in range (1,2):
        index_list=[]
        for ii in range(1,len(best_x)-10):
            #if y[ii+2] > 0.65:
            for kk in range (ii+1, len(best_x)-1):
                if y[kk] > 1.05*y[ii] :
                    index_list.append(ii)
        best_x=np.delete(best_x,[index_list])
        best_y=np.delete(best_y,[index_list])
        best_y_u =np.delete(best_y_u,[index_list])


    ax.errorbar(best_x,best_y,best_y_u,zorder = 5, fmt = 'bo')

    ########## time to fit the best results
    p0, fitfunc, fitfunc_str = common.fit_general_exponential(0.5,0.5, 0., 100.,2) ## one could start to make cool guesses based on the data but i refrain
    fit_result = fit.fit1d(best_x,best_y,None,p0=p0,fitfunc=fitfunc,do_print=True,ret=True,fixed = [0,2])

    plot.plot_fit1d(fit_result,np.linspace(0,np.amax(best_x)*1.05,201),ax=ax,plot_data=False)
    plt.savefig(os.path.join(a.folder, 'analyzed_result.pdf'),
    format='pdf')
    plt.savefig(os.path.join(a.folder, 'analyzed_result.png'),
    format='png')

    print 'plots are saved in'
    print a.folder

    plt.show()
