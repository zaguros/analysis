import numpy as np
import os,sys
# import qt

sys.path.append("/measuring/")
from analysis.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.fitting import fit, common
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt
import matplotlib.cm as cm
reload(common)
reload(plot)
""""
Takes 4 datapoints instead of 2 to combine them for amplitude
"""

def get_max_amplitude(px,py,mx,my,px_u,py_u,mx_u,my_u):

    ### init
    a_arr,a_u_arr = np.array([]),np.array([])

    ##calc

    b_arr = np.sqrt((px-mx)**2+(py-my)**2)/2
    b_u_arr =  np.sqrt((px-mx)**2*(px_u**2+mx_u**2) +(py-my)**2*(py_u**2+my_u**2))/(np.sqrt((px-mx)**2+(py-my)**2)*2)

    return b_arr,b_u_arr


def get_raw_data(carbon,**kw):
    """
    extracts the data for one carbon from a positive and negative file.
    returns arrays for the contrast along x, y and the respective uncertainties.
    """
    older_than = kw.pop('older_than',None)
    ssro_tstamp = kw.pop('ssro_tstamp',None)
    plotwidth = kw.pop('width',None)
    Ndiff = kw.pop('Ndiff',None)

    if ssro_tstamp == None:
        ssro_calib_folder = toolbox.latest_data(contains = 'SSRO')
        print ssro_calib_folder
    if ssro_tstamp != None:
        ssro_calib_folder = toolbox.data_from_time(ssro_tstamp, folder = None)

    
    #if Ndiff == None:
    #    search_string_pos = 'Sweep_carbon_Gate_positive_C'+str(carbon)+'_width_'+str(plotwidth)+'ns'
    #    search_string_neg = 'Sweep_carbon_Gate_negative_C'+str(carbon)+'_width_'+str(plotwidth)+'ns'

    #if Ndiff != None:
    if carbon == 1:
        search_string_pos = 'Sweep_carbon_Gate_positive_C'+str(carbon)+'_det_'+str(plotwidth)+'nas_Ndiff_'+str(Ndiff)
        search_string_neg = 'Sweep_carbon_Gate_negative_C'+str(carbon)+'_det_'+str(plotwidth)+'nas_Ndiff_'+str(Ndiff)
    if carbon == 2 or carbon == 5:
        search_string_pos = 'Sweep_carbon_Gate_positive_C'+str(carbon)+'_det_'+str(plotwidth)+'nas_Ndiff_'+str(Ndiff)
        search_string_neg = 'Sweep_carbon_Gate_negative_C'+str(carbon)+'_det_'+str(plotwidth)+'nas_Ndiff_'+str(Ndiff)
    
    #if plotwidth == None and Ndiff == None:
    #    search_string_pos = 'Sweep_carbon_Gate_positive_C'+str(carbon)
    #    search_string_neg = 'Sweep_carbon_Gate_negative_C'+str(carbon)

    #search_string_pos = 'Sweep_carbon_Gate_positive_C'+str(carbon)+'_width_'+str(plotwidth)+'nas_Ntot_'+str(Ndiff)
    #search_string_neg = 'Sweep_carbon_Gate_negative_C'+str(carbon)+'_width_'+str(plotwidth)+'nas_Ntot_'+str(Ndiff)
    


    
    folder_pos = toolbox.latest_data(contains = search_string_pos, older_than=older_than)
    folder_neg = toolbox.latest_data(contains = search_string_neg, older_than=older_than)

    #print(folder_pos)
    #print(folder_neg)

    a = mbi.MBIAnalysis(folder_pos)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)

    b = mbi.MBIAnalysis(folder_neg)
    b.get_sweep_pts()
    b.get_readout_results(name='adwindata')
    b.get_electron_ROC(ssro_calib_folder = ssro_calib_folder)

    a.p0 = 2*a.p0-1; a.u_p0 = 2*a.u_p0
    b.p0 = 2*b.p0-1; b.u_p0 = 2*b.u_p0

    ###combine positive & negative:

    a.p0 = (a.p0 - b.p0)/2
    a.u_p0 = ((a.u_p0**2 + b.u_p0**2)**0.5)/2


    mx_arr,mx_u_arr = np.array([]),np.array([])
    my_arr,my_u_arr = np.array([]),np.array([])

    px_arr,px_u_arr = np.array([]),np.array([])
    py_arr,py_u_arr = np.array([]),np.array([])

    gates,xlabels = np.array([]),np.array([]) ### used parameters

    ### sort into X and y lists.
    for (pt,val,val_u) in zip(a.sweep_pts.reshape(-1),a.p0.reshape(-1),a.u_p0.reshape(-1)):
        if 'mX' in pt:
            mx_arr = np.append(mx_arr,val)
            mx_u_arr = np.append(mx_u_arr,val_u)
            gates = np.append(gates,'N1 = '+str(pt[3:5])+',\ntau1 = '+str(pt[8:17])+', N2 = '+str(pt[18:20])+', \ntau2 = '+str(pt[23:32])+', phase = '+str(pt[32:]))
            xlabels= np.append(xlabels,str(pt[33:]))        
        elif 'mY' in pt:
            my_arr = np.append(my_arr,val)
            my_u_arr = np.append(my_u_arr,val_u)
        elif 'pY' in pt:
            py_arr = np.append(py_arr,val)
            py_u_arr = np.append(py_u_arr,val_u)
        elif 'pX' in pt:
            px_arr = np.append(px_arr,val)
            px_u_arr = np.append(px_u_arr,val_u)
        
    return gates,mx_arr,my_arr,mx_u_arr,my_u_arr,px_arr,py_arr,px_u_arr,py_u_arr,folder_pos,xlabels

def get_gate_fidelity(carbon,**kw):
    """
    gets data, plots it and prints the gate parameters for maximum bloch vector length.
    """

    older_than = kw.pop('older_than',None)
    ssro_tstamp = kw.pop('ssro_tstamp',None)
    #plotwidth = kw.pop('width',None)
    #Ndiff = kw.pop('Ndiff',None)
    pwlist=kw.pop('pwlist',[0])
    Ndifflist=kw.pop('Ndifflist',[0])
    for plotwidth in pwlist:

        fig = plt.figure(figsize=(12,6))
        fig.suptitle('Sweep of N1 - N2 (Ntot is constant) for width ='+str(plotwidth)+'ns.', fontsize=20)
        legendlst=[]
        bestblist=[]
        bestbindlist=[]
        gatelist=[]

        for ind,Ndiff in enumerate(Ndifflist):
                
            
            ax_ind = plt.subplot()
           
            gates,mx,my,mx_u,my_u,px,py,px_u,py_u,folder_pos,xlabels = get_raw_data(carbon,older_than = older_than,ssro_tstamp = ssro_tstamp,Ndiff=Ndiff,width=plotwidth)


            #print mx,my,px,py
            
            b,b_u = get_max_amplitude(px,py,mx,my,px_u,py_u,mx_u,my_u)

            bestblist.append(np.amax(b))
            bestbindlist.append(np.argmax(b))
            gatelist.append(gates)
         

            legendlst.append('N2 - N1 =' + str(Ndiff))
        

            errobarplot_= ax_ind.errorbar(np.arange(len(gates)),b,yerr=b_u)
            ax_ind.set_xticks(np.arange(len(gates)))
            ax_ind.set_xticks(np.array(range(len(gates))))
            ax_ind.set_xticklabels(xlabels, rotation=90)
            ax_ind.set_ylim([0,1.1])
            plt.xlabel('Phase between decoupling sequence of gates')
            plt.ylabel('Amplitude of measurement')
            
        
        
        plt.legend(legendlst, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
        plt.savefig(os.path.join(folder_pos, 'Sweep_of_gates_width_'+str(plotwidth)+'.png'), format='png',bbox_inches='tight', pad_inches=0.1)
        print 'best gate configuration at: ', gatelist[np.argmax(bestblist)][bestbindlist[np.argmax(bestblist)]]
        print 'Oscilation Amplitude: ', np.amax(bestblist)
        plt.show()
        plt.close('all')
        print('FIG saved in'+str(folder_pos))




