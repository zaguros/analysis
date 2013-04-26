import numpy as np
import os
#import sys

def get_levels():
    """
    Returns an array with the ES energies as a function of strain Ex, 
    also returned
    """
    Ex=np.linspace(0,20,50)
    return Ex,np.array([np.sort(get_ES(E_field=[i,0,0])[0]) for i in Ex])

def get_ES_ExEy(Ex,Ey,fast=False):
    """
    Returns the six transition energies in GHz of the ES of the NV centre, 
    when given the Energies of the Ex and Ey transitions in GHz
    """
    
    strain=abs(Ex-Ey)/2.0
    offset=np.min([Ey,Ex])+strain
    if fast:
        return np.sort(get_ES_fast(offset,strain))
    return np.sort(get_ES(E_field=[strain,0,0],Ee0=offset-1.94)[0])

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
    return np.sort(x),y
    
def get_ES(E_field=[0.,0.,0.],B_field=[0.,0.,0.],Ee0=-1.94,transitions=True):
    """Returns the eigenvalues and eigenvectors of the NV excited state 
    pertubation matrix.
    inputs:
    - E-field xyz vector in GHz
    - B-field xyz vector in GHz
    - Energy offset for the eigenvalues
    - boolean transitions - whether to return the transition energies 
    (ms0 energies increased by the zero-field splitting)
    """
    
    Ex = E_field[0]
    Ey = E_field[1]
    Ez = E_field[2]
    Bx = B_field[0]
    By = B_field[1]
    Bz = B_field[2]
    lambdaA2=0.1              #observed
    lambda_par=5.3           #observed
    lambda_ort=1.5*lambda_par      #unknown, calculated by Maze, p9
    D1A1=2.87/3           #observed
    D2A1=1.42/3            #observed
    D2E1=1.55/2             #observed
    D2E2=.2/np.sqrt(2)        #observed

    w2=np.sqrt(2)
    
    Vss = np.matrix([[D2A1, 0, D2E2*w2, 0, 0, 0],
                    [0, D2A1, 0, D2E2 *w2, 0, 0],
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
    Vb = np.matrix([[0, 1j*(Bz + lambdaA2*Bz), 1j*(By)/w2, 1j*(Bx)/w2, 0, 0],
                    [-1j*(Bz + lambdaA2*Bz), 0, 1j*(Bx)/w2, -1j*(By)/w2, 0, 0],
                    [-1j*(By)/w2, -1j*(Bx)/w2, 0, -1j*lambdaA2*Bz, 1j*(By)/w2, -1j*(Bx)/w2],
                    [-1j*(Bx)/w2, 1j*(By)/w2, -1j*lambdaA2*Bz,    0, -1j*(Bx)/w2, -1j*(By)/w2],
                    [0, 0, -1j*(By)/w2, 1j*(Bx)/w2, 0, 0],
                    [0, 0, 1j*(Bx)/w2, 1j*(By)/w2, 0, 0]])
      
   
    VGSoffset =  np.diag([0, 0, 3*D1A1, 3*D1A1, 0, 0]) if transitions else 0
  
   
    V = Vss + Vso + Ve + Vb + VGSoffset
    
    w,v=np.linalg.eig(V)
    
    return np.real(w+Ee0),v
    
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
    
def fit_laserscan(x,y,points=100,strain_min=0,strain_max=10,linewidth=0.1,fast=False):
    """
    Fit a laserscan with the theoretical NV spectrum, with frequency offset 
    and strain as a free variable. the frequqncy offset search space is taken as 
    the frequqcy axis range.
    Input: x : frequency axis of the laser scan
           y : counts of the scan
           points: precision used in the fitting, eg the stepsize in both 
                   frequency offset and strain.
           strain_min/max: defines the search space for the strain fit, in GHz.
           linewidth: expected linewith used in the fitting procedure in GHz.
    Returns: the fitted peak Ex and Ey positions as [Ex,Ey] in GHz
    """
    
    l=len(x)
    fmax=max(x)
    fmin=min(x)  
    ls=points
    A=np.zeros((ls,ls,l),dtype=np.bool)

    fEys=np.linspace(fmin,fmax,ls)
    fstrs=np.linspace(strain_min,strain_max,ls)
    lw=linewidth
    for i,fEy in enumerate(fEys):
        for j,fstr in enumerate(fstrs):
            flev=get_ES_ExEy(fEy,fEy+fstr,fast)
            A[i,j]=np.logical_and(x>(flev[0]-lw/2),x<(flev[0]+lw/2))+ \
                        np.logical_and(x>(flev[1]-lw/2),x<(flev[1]+lw/2))+ \
                        np.logical_and(x>(flev[2]-lw/2),x<(flev[2]+lw/2))+ \
                        np.logical_and(x>(flev[3]-lw/2),x<(flev[3]+lw/2))+ \
                        np.logical_and(x>(flev[4]-lw/2),x<(flev[4]+lw/2))+ \
                        np.logical_and(x>(flev[5]-lw/2),x<(flev[5]+lw/2))

    result=np.sum(A*y,axis=2)
    im,jm=np.where(result==np.max(result))
    if len(im)>1:
        print 'Multiple laserscan fits found, returned largest strain result only:' 
        d=0.
        imm=0
        for ii,imi in enumerate(im):
            if fstrs[jm[ii]]>d:
                d=fstrs[jm[ii]]
                imm=ii
            print ii,':', 'Ey:', fEys[imi], 'Ex:' , fEys[imi]+fstrs[jm[ii]]
        return [fEys[im[imm]],fEys[im[imm]]+fstrs[jm[imm]]]
    return [fEys[im],fEys[im]+fstrs[jm]]
