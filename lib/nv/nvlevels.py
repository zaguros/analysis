#Bas Hensen 2013

import numpy as np
import os
import scipy.constants as spc
#import sys

def get_levels(**kw):
    """
    Returns an array with the ES energies as a function of strain Ex, 
    also returned
    """
    Ex=kw.pop('strainvals', np.linspace(0,20,50))
    return Ex,np.array([np.sort(get_ES(E_field=[i,0,0], **kw)[0]) for i in Ex])

def get_ES_ExEy(Ex,Ey,fast=False,transitions=True):
    """
    Returns the six transition energies in GHz of the ES of the NV centre, 
    when given the Energies of the Ex and Ey transitions in GHz
    """
    
    strain=abs(Ex-Ey)/2.0
    offset=np.min([Ey,Ex])+strain
    if fast:
        return np.sort(get_ES_fast(offset,strain,transitions=transitions))
    return np.sort(get_ES(E_field=[strain,0,0],Ee0=offset-1.94,transitions=transitions)[0])

def get_ES_fast(f0,D,transitions=True):
    D=D/0.749
    ms1_off=0
    if transitions:
        ms1_off=2.87
    return np.array([f0 - ms1_off - 3.8821 - 0.01856*D - 0.06452*D**2 + 0.00284*D**3 - 4.925e-5*D**4,
                     f0 - ms1_off - 3.9021 + 0.0157*D  - 0.05397*D**2 + 0.00202*D**3 - 3.09e-5*D**4,
                     f0 - 0.749*D,
                     f0 + 0.749*D,
                     f0 - ms1_off + 5.1379 + 0.0159*D  + 0.06463*D**2 - 0.00287*D**3 + 4.959e-5*D**4,
                     f0 - ms1_off + 8.2579 - 0.00975*D + 0.05384*D**2 - 0.00202*D**3 + 3.083e-5*D**4])

def get_ES_ExEy_plottable(Ex,Ey,height):
    """
    Returns an array plottable with qt.Plot2D of the six transition energies 
    in GHz of the ES of the NV centre, when given the Energies of the Ex and 
    Ey transitions in GHz
    """
    x=get_ES_ExEy(Ex,Ey)
    y=np.zeros(3*len(x))
    for ii in range(len(x)):
        x=np.append(x,x[ii]-0.0001)
        x=np.append(x,x[ii]+0.0001)
        y[3*ii+1]=height
    return [np.sort(x),y]
    
def get_ES(E_field=[0.,0.,0.],B_field=[0.,0.,303.],Ee0=-1.94, **kw):
    """Returns the eigenvalues and eigenvectors of the NV excited state 
    pertubation matrix.
    inputs:
    - E-field xyz vector in GHz
    - B-field xyz vector in Gauss
    - Energy offset for the eigenvalues
    - boolean transitions - whether to return the transition energies 
    (ms0 energies increased by the zero-field splitting)s
    """
    # [1]: Doherty, M. W. et al. Physics Reports 528, 1-45 (2013)
    # [2]: Maze, J. R. et al. New J. Phys. 13, 025025 (2011)
    # see also:
    # Doherty, M. W., Manson, N. B., Delaney, P. and Hollenberg, L. C. L. New J. Phys. 13, 025019 (2011).
    # K:\ns\qt\Diamond\Reports and Theses\MSc\Bas Hensen\Hensen_msc_mail 2011-10-07.pdf

    ### THT: this actually does not handle correctly the magnetic field yet in the case of transitions due to the GS Zeeman splittings

    Ex = E_field[0]
    Ey = E_field[1]
    Ez = E_field[2]
    
    mu_B=spc.e*spc.hbar/(2*spc.m_e)/spc.h/1e9  #GHz/Tesla
    Bx = mu_B*B_field[0]*1e-4 #GHz
    By = mu_B*B_field[1]*1e-4 #GHz
    Bz = mu_B*B_field[2]*1e-4 #GHz
    
    #Bfield
    lambdaA2=.1                  #observed, [1]   
    g_es_par = 2.                #RT value, likely to be different at LT! [1]
    g_es_ort = 2.                #RT value, likely to be different at LT! [1]              

    lambda_par=5.3               #observed, [1] 
    #lambda_ort_2=1.5*lambda_par #unknown, calculated by [2]
    D1A1=2.88/3                  #observed, [1]
    D2A1=1.42/3                  #observed, [1]
    D2E1=1.55/2                  #observed, [1]
    D2E2=.2/np.sqrt(2)           #observed, [1] AKA lambda_es_ort

    w2=np.sqrt(2)
    
    Vss = np.matrix([[D2A1, 0, D2E2*w2, 0, 0, 0],
                    [0, D2A1,  0, D2E2*w2, 0, 0],
                    [D2E2*w2, 0, -2*D2A1, 0, 0, 0],
                    [0, D2E2*w2, 0, -2*D2A1, 0, 0],
                    [0, 0, 0, 0, D2A1 - 2*D2E1, 0],
                    [0, 0, 0, 0, 0, D2A1 + 2*D2E1]])
            
    Vso = np.diag([-lambda_par, -lambda_par, 0, 0, lambda_par, lambda_par])
    
    Ve = np.matrix([[Ez, 0, 0, 0, Ey, Ex],
                   [0, Ez, 0, 0, -Ex, Ey],
                   [0, 0, Ez + Ex, -Ey, 0, 0],
                   [0, 0, -Ey, Ez - Ex, 0, 0],
                   [Ey, -Ex, 0, 0, Ez, 0],
                   [Ex, Ey, 0, 0, 0, Ez]])
    Vb = np.matrix([[0,  1j*(g_es_par*Bz + lambdaA2*Bz), 1j*(g_es_ort*By)/w2,  1j*(g_es_ort*Bx)/w2, 0, 0],
                    [-1j*(g_es_par*Bz + lambdaA2*Bz), 0, 1j*(g_es_ort*Bx)/w2, -1j*(g_es_ort*By)/w2, 0, 0],
                    [-1j*(g_es_ort*By)/w2, -1j*(g_es_ort*Bx)/w2, 0,                 -1j*lambdaA2*Bz, 1j*(g_es_ort*By)/w2, -1j*(g_es_ort*Bx)/w2],
                    [-1j*(g_es_ort*Bx)/w2,  1j*(g_es_ort*By)/w2, -1j*lambdaA2*Bz,    0,             -1j*(g_es_ort*Bx)/w2, -1j*(g_es_ort*By)/w2],
                    [0, 0, -1j*(g_es_ort*By)/w2, 1j*(g_es_ort*Bx)/w2,  0, 1j*(g_es_par*Bz - lambdaA2*Bz)],
                    [0, 0, 1j*(g_es_ort*Bx)/w2,  1j*(g_es_ort*By)/w2, -1j*(g_es_par*Bz - lambdaA2*Bz), 0]])
      
   
   
    if kw.pop('transitions', False):
        print 'transitions kw deprecated, use function get_transitions instead'
        VGSoffset =  np.diag([0, 0, 3*D1A1, 3*D1A1, 0, 0])
    else:
        VGSoffset = 0.

    V = Vss + Vso + Ve + Vb + VGSoffset
    
    w,v=np.linalg.eig(V)
    
    return np.real(w+Ee0),v
 

def get_GS(E_field=[0.,0.,0.],B_field=[0.,0.,0.], **kw):

    Ex = E_field[0]
    Ey = E_field[1]
    Ez = E_field[2]
    
    mu_B=spc.e*spc.hbar/(2*spc.m_e)/spc.h/1e9  #GHz/Tesla
    Bx = mu_B*B_field[0]*1e-4 #GHz
    By = mu_B*B_field[1]*1e-4 #GHz
    Bz = mu_B*B_field[2]*1e-4 #GHz

    D1A1=2.88/3                  #observed, [1]
    g_gs_par = 2.                #approx
    g_gs_ort = 2.                #R
    
    Vss = np.diag([-2*D1A1, D1A1, D1A1])
    Ve = np.diag([Ez, Ez, Ez])
    Vb = np.matrix([[0, 1j*(g_gs_ort*By), -1j*(g_gs_ort*Bx)],
                    [-1j*(g_gs_ort*By), 0, -1j*(g_gs_par*Bz)],
                    [1j*(g_gs_ort*Bx), 1j*(g_gs_par*Bz), 0]])

    V = Vss + Ve + Vb
    
    w,v=np.linalg.eig(V)
    return np.real(w),v

def get_optical_transitions(show_FB_E_transitions=True, show_FB_A_transitions=False, show_E_prime_flip_transitions=False, **kw):

    E_GS=np.sort(get_GS(**kw)[0])
    E_ES=np.sort(get_ES(**kw)[0])
    allowed_transitions=np.array([E_ES[2]-E_GS[0],
                                 E_ES[3]-E_GS[0],
                                 E_ES[0]-E_GS[1],#E_ES[0]-E_GS[2],
                                 E_ES[1]-E_GS[2],#E_ES[1]-E_GS[1],
                                 E_ES[4]-E_GS[1],E_ES[4]-E_GS[2],
                                 E_ES[5]-E_GS[1],E_ES[5]-E_GS[2]])
    E_prime_flip_transitions = np.array([E_ES[0]-E_GS[2],
                                        E_ES[1]-E_GS[1]])
    FB_E_transitions=np.array([E_ES[2]-E_GS[1],E_ES[2]-E_GS[2],
                               E_ES[3]-E_GS[1],E_ES[3]-E_GS[2]])
    FB_A_transitions=np.array([E_ES[0]-E_GS[0],
                               E_ES[1]-E_GS[0],
                               E_ES[4]-E_GS[0],
                               E_ES[5]-E_GS[0]])
    transitions = allowed_transitions
    if show_FB_E_transitions: 
        transitions = np.append(transitions, FB_E_transitions)
    if show_FB_A_transitions: 
        transitions = np.append(transitions, FB_A_transitions)
    if show_E_prime_flip_transitions: 
        transitions = np.append(transitions, E_prime_flip_transitions)
    return transitions


def fit_laserscan_from_file(file, plot=False,plot_save=True, **kw):
    d=np.load(file)['data']
    x=d[:,1]
    y=d[:,2]
    fit_res=fit_laserscan(x,y,**kw)
    if plot:
        from matplotlib import pyplot as plt
        plt.figure()
        plt.plot(x,y)
        x2,y2=get_ES_ExEy_plottable(fit_res[0],fit_res[1],max(y))
        plt.plot(x2,y2)
        if plot_save:
            plt.savefig(os.path.join(os.path.split(file)[0],'fitted_scan.png'))
    return fit_res
    
def fit_laserscan(x,y,points=100,strain_range=(0,10),Ex_range=None,Ey_range=None,linewidths=[0.1,0.1,0.4,0.2,0.1,0.1],fast=False):
    """
    Fit a laserscan with the theoretical NV spectrum, with frequency offset 
    and strain as a free variable. 
    Input: x : frequency axis of the laser scan
           y : counts of the scan
           points: precision used in the fitting, eg the stepsize in both 
                   frequency offset and strain.
           strain_range: defines the search space for the strain fit, in GHz, in the form (min,max).
           Ey_range: defines the search space for the Ey frequency, in GHz, in the form (min,max).
                     if y_range is None, the frequency axis range (x is taken)
           linewidth: expected linewith used in the fitting procedure in GHz.
    Returns: the fitted peak Ex and Ey positions as [Ex,Ey] in GHz
    """
    
    l=len(x)
    fmax=max(x)
    fmin=min(x)
    ls=points
    
    if Ey_range==None:
        fEys=np.unique(x[y>(2*np.median(y))])
    else:
        fEys=np.linspace(Ey_range[0],Ey_range[1],ls)
        
    if Ex_range==None:
        fExs=np.copy(fEys)
    else:
        fExs=np.linspace(Ex_range[0],Ex_range[1],ls)
    lw=linewidths
    li=len(fEys)
    lj=len(fExs)
    A=np.zeros((li,lj,l),dtype=np.bool)
    for i,fEy in enumerate(fEys):
        for j,fEx in enumerate(fExs):
            if fEx<fEy or (fEx-fEy)<strain_range[0] or (fEx-fEy)>strain_range[1]:
                continue
            flev=get_ES_ExEy(fEx,fEy,fast)
            for k in range(6):
                A[i,j]=A[i,j]+ np.logical_and(x>(flev[k]-lw[k]/2),x<(flev[k]+lw[k]/2))
    result=np.sum(A*y,axis=2)
    im,jm=np.where(result==np.max(result))
    if len(im)>1:
        print 'Multiple laserscan fits found, returned largest strain result only:' 
        d=0.
        imm=0
        for ii,imi in enumerate(im):
            if (fExs[jm[ii]]-fEys[im[ii]])>d:
                d=(fExs[jm[ii]]-fEys[im[ii]])
                imm=ii
            print ii,':', 'Ey:', fEys[im[ii]], 'Ex:' , fExs[jm[ii]]
        return np.array([fExs[jm[imm]],fEys[im[imm]],])
    return np.ravel([fExs[jm],fEys[im],])

def get_ExEy_from_two_levels(f1,i1,f2,i2, precision=0.03, fast=True):
    
    for str_split in np.linspace(0,20,20/precision):
        levels=get_ES_ExEy(0,str_split,fast)
        offset=(f1-levels[i1])
        levels=levels+offset
        #print levels
        if abs(f2-levels[i2])<precision:
            return levels[2]+str_split, levels[2]
