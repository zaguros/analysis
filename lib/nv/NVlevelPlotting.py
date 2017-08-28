import sys
import numpy as np
import scipy

sys.path.append(r'//Users/humphreys/Repositories/')

from analysis.lib.nv import nvlevels; reload(nvlevels)
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm


def plot_ES_energies(Ex,Ey,Strain_start,Strain_end,pts,B_field=[0.,0.,0.], ax=None):
    Strain = Ex-Ey
    ExZero = (Ex+Ey)/2 
    EyZero = (Ex+Ey)/2 

    result_list = np.zeros([6,pts])
    Strainarray = np.linspace(Strain_start,Strain_end,pts)
    for i in range(pts):
        result_list[:,i] = (nvlevels.get_ES_SpinComp_ExEy(ExZero+Strainarray[i]/2,EyZero-Strainarray[i]/2,B_field=B_field))[0]
    
    if ax == None:
        fig = plt.figure()
        ax = plt.subplot()

    legend_list = ["E-'","E+'",'Ey','Ex','A1','A2']
    for res,legend in zip(result_list,legend_list):
        plt.plot(Strainarray,res,label = legend)

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylabel('Frequency (GHz)')
    plt.xlabel('Strain splitting (GHz)')
    plt.title('Excited state energies')
    plt.show()
    plt.close("all")
   
def plot_GS_energies(B_start,B_end,pts,E_field=[0.,0.,0.], ax=None):

    result_list = np.zeros([3,pts])
    Barray = np.linspace(B_start,B_end,pts)
    for i in range(pts):
        result_list[:,i] = (nvlevels.get_GS_SpinComp(E_field=E_field,B_field=[0,0,Barray[i]]))[0] +1.92

    if ax == None:
        fig = plt.figure()
        ax = plt.subplot()

    legend_list = ["A_0", "E_-", "E_+"]
    for res,legend in zip(result_list,legend_list):
        plt.plot(Barray,res,label = legend)

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylabel('Frequency (GHz)')
    plt.xlabel('B field (Gauss)')
    plt.ylim([-0.2,plt.ylim()[1]])
    plt.title('Ground state energies')
    plt.show()
    plt.close("all")

def plot_eigenstates(Ex,Ey,Strain_start,Strain_end,pts,B_field=[0.,0.,0.]):
    
    EZero = (Ex+Ey)/2 

    Strainarray = np.linspace(Strain_start,Strain_end,pts)
        
    eigen_components = np.zeros([6,6,pts])

    for ii in range(pts):
        ES, v = nvlevels.get_ES_SpinComp_ExEy(EZero+Strainarray[ii]/2,EZero-Strainarray[ii]/2,B_field=B_field)
        eigen_components[:,:,ii] = np.abs(v.T)**2
           
    legend_list = ["E-'","E+'",'Ey','Ex','A1','A2']

    shape = np.shape(eigen_components)
    fig = plt.figure(figsize = (8,12))
    for x in range(shape[0]):
        ax = plt.subplot(shape[0],1,x+1)

        if x == 0:
            for y,label in zip(range(shape[1]),legend_list):
                plt.plot(Strainarray,eigen_components[x,y,:], label=label)
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) 
        else:
            for y in range(shape[1]):
                plt.plot(Strainarray,eigen_components[x,y,:])
        plt.ylabel('Component')
        if x == shape[0]-1:
            plt.xlabel('Strain splitting (GHz)')
    
    # fig = matplotlib.pyplot.gcf()
    # fig.set_size_inches(10.5, 15, forward=True)      
    plt.show()
    plt.close('all')

def plot_transitions(Ex,Ey,Strain_start,Strain_end,pts,B_field=[0.,0.,0.],m0=True,m1=True,p1=True,return_dict = False, ax=None):

    ExZero = (Ex+Ey)/2 
    EyZero = (Ex+Ey)/2 
    
    Strainarray = np.linspace(Strain_start,Strain_end,pts)
    
    if return_dict:
        transition_results = {}
        trans_keys = []
        if m0:
            trans_keys.append('ms0')
            transition_results['ms0'] = np.empty([2,pts])
        if m1:
            trans_keys.append('msm1')
            transition_results['msm1'] = np.empty([3,pts])
        if p1:
            trans_keys.append('msp1')
            transition_results['msp1'] = np.empty([3,pts])
    else:
        transition_results = np.empty([2*bool(m0)+3*bool(m1)+3*bool(p1),pts])   
    
    for i in range(pts):
            slice_list = nvlevels.get_transitions_ExEy(ExZero+Strainarray[i]/2,EyZero-Strainarray[i]/2,B_field=B_field,return_dict = return_dict, show_FB_E_transitions=False,show_FB_A_transitions=False, show_E_prime_flip_transitions=False, show_ms0_transitions = m0, show_m1_transitions = m1, show_p1_transitions = p1)
            if not return_dict:
                transition_results[:,i] = slice_list
            else:
                for key in trans_keys:
                    transition_results[key][:,i] = slice_list[key]

    plot_style = {'msp1' : '--', 'msm1': '-', 'ms0' : '-'}

    if ax == None:
        fig = plt.figure()
        ax = plt.subplot()

    if type(result_list) == np.ndarray:
        for res in result_list:
            plt.plot(Strainarray,res)
    elif type(result_list) == dict:
        for key in result_list:
            for i,res in enumerate(result_list[key]):
                if i == 0: 
                    baseline, = plt.plot(Strainarray,res,label = key, ls=plot_style[key])
                else:
                    plt.plot(Strainarray,res, c=baseline.get_color(), ls=plot_style[key])
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylabel('Frequency (GHz)')
    plt.xlabel('Strain splitting (GHz)')
    plt.show()
    plt.close("all")
    
def plot_transitions_w_strengths(Ex,Ey,Strain_start,Strain_end,pts,B_field=[0.,0.,0.],m0=True,m1=True,p1=True,ax=None,log_scale =False):
    
    Strain = Ex-Ey
    ExZero = (Ex+Ey)/2 
    EyZero = (Ex+Ey)/2 
    
    Strainarray = np.linspace(Strain_start,Strain_end,pts)

    trans_keys = []
    if m0:
        trans_keys.append('ms0')
    if m1:
        trans_keys.append('msm1')
    if p1:
        trans_keys.append('msp1')
        
    transitions = {}
    for key in trans_keys:
        transitions[key] = {}
        transitions[key]['strength'] = np.empty([6,pts])
        transitions[key]['freq'] = np.empty([6,pts])
        
    for ii in range(pts):
        slice_list = nvlevels.get_optical_transition_strengths_ExEy(ExZero+Strainarray[ii]/2,EyZero-Strainarray[ii]/2,B_field=B_field,
                    show_ms0_transitions=m0,show_m1_transitions=m1,show_p1_transitions=p1)
        
        for key in trans_keys:
            transitions[key]['strength'][:,ii] = slice_list[key]['strength']
            transitions[key]['freq'][:,ii] = slice_list[key]['freq']
            
    color_map_key = {'msp1' : 'Blues', 'msm1': 'Greens', 'ms0' : 'Reds'}

    if ax == None:
        fig = plt.figure(figsize = (8,6))

        ax = plt.subplot()

    for key in transitions:
        
        for ii in range(np.shape(transitions[key]['freq'])[0]):

            points = np.array([Strainarray, transitions[key]['freq'][ii,:]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            if not log_scale:
                lc = LineCollection(segments, cmap=plt.get_cmap(color_map_key[key]),norm=plt.Normalize(vmin=0,vmax=1))
                lc.set_array(transitions[key]['strength'][ii,:])
            else:
                lc = LineCollection(segments, cmap=plt.get_cmap(color_map_key[key]),norm=plt.Normalize(vmin=-3,vmax=0))
                transitions[key]['strength'][ii,(transitions[key]['strength'][ii,:] == 0)] = 1e-5
                lc.set_array(np.log10(transitions[key]['strength'][ii,:]))
            lc.set_linewidth(1)
            ax.add_collection(lc)

    plt.ylabel('Frequency (GHz)')
    plt.xlabel('Strain splitting (GHz)')

    ax.autoscale_view(True,True,True)
    
    plt.show()
    plt.close("all")