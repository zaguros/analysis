"""
Evaluates the density matrix of a combined electron nuclear spin state.
Based upon 12 measurements for different initial measurements of the electron
"""


import numpy as np
import copy as cp
import os,h5py
from analysis.lib.tools import toolbox as tb
from analysis.lib.m2.ssro import mbi; reload(mbi)
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from analysis.lib.m2.ssro import ssro
from analysis.lib.m2 import m2
from analysis.lib.lde import sscorr; reload(sscorr)




####################
#                  # 
# Helper functions #
#                  #
####################

def generate_pauli_matrices():
    return [np.matrix([[1,0],[0,1]],dtype=complex),np.matrix([[0,1],[1,0]],dtype=complex),np.matrix([[0,-1j],[1j,0]],dtype=complex),np.matrix([[1,0],[0,-1]],dtype=complex)]

def get_RO_results(folder,ssro_calib_folder):
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder)
    y_a= ((a.p0.reshape(-1)[:])-0.5)*2
    y_err_a = 2*a.u_p0.reshape(-1)[:] 

    return y_a,y_err_a

def get_correlations(**kw):

    ### pull data
    ssro_calib_folder = kw.pop('ssro_calib_folder',None)
    ssro_calib_timestamp = kw.pop('ssro_calib_timestamp',None)
    search_string = kw.pop('search_string','el_13C_dm_')
    base_folder = kw.pop('base_folder',None)

    if ssro_calib_folder == None:
        if ssro_calib_timestamp == None: 
            ssro_calib_folder = tb.latest_data('SSROCalibration')
        else:
            ssro_dstmp, ssro_tstmp = tb.verify_timestamp(ssro_calib_timestamp)
            ssro_calib_folder = tb.latest_data(contains = ssro_tstmp,older_than = str(int(ssro_dstmp)+1)+'_'+ssro_tstmp,folder = base_folder)

    ### for basis assignment, see onenote 2016-08-24 or alternatively mathematica file E_13C_Bell_state.nb
    tomo_pulses_p = ['none','x','my'] ### used when looking for folders
    tomo_pulses_m = ['X','mx','y'] ### used when looking for folders

    f_list_pos_p,f_list_neg_p,f_list_pos_m,f_list_neg_m = [],[],[],[]
    exp_values,exp_vals_u = np.array([]),np.array([]) ### this will be a list with all combined expectation values such as XX or YZ
    sweep_pts = [] ### this will be a list of the bases associated expectation values

    ### carbon 1-qubit correlations
    c_x,c_y,c_z = 0.,0.,0.
    c_x_u,c_y_u,c_z_u = 0.,0.,0.

    ### this dictionary is used to translate the electron RO pulse into a measurement basis
    tomo_pulse_translation_dict = {'none': 'Z','X':'Z','x':'Y','mx':'Y','y':'X','my':'X'}

    for p,m in zip(tomo_pulses_p,tomo_pulses_m):
        f_list_pos_p.append(tb.latest_data(search_string + p+'_positive' ,folder = base_folder,return_timestamp = False,**kw))
        f_list_neg_p.append(tb.latest_data(search_string + p+'_negative' ,folder = base_folder,return_timestamp = False,**kw))
        f_list_pos_m.append(tb.latest_data(search_string + m+'_positive' ,folder = base_folder,return_timestamp = False,**kw))
        f_list_neg_m.append(tb.latest_data(search_string + m+'_negative' ,folder = base_folder,return_timestamp = False,**kw))

        #### now also calculate the measured contrast
        y_a,y_err_a = get_RO_results(f_list_pos_p[-1],ssro_calib_folder)
        y_b,y_err_b= get_RO_results(f_list_neg_p[-1],ssro_calib_folder)


        y_c,y_err_c = get_RO_results(f_list_pos_m[-1],ssro_calib_folder)
        y_d,y_err_d= get_RO_results(f_list_neg_m[-1],ssro_calib_folder)

        y_pos_electron   = (y_a - y_b)/2.
        y_pos_electron_u = 1./2*(y_err_a**2 + y_err_b**2)**0.5
        y_neg_electron   = (y_c - y_d)/2.
        y_neg_electron_u =  1./2*(y_err_c**2 + y_err_d**2)**0.5

        #### now do the final combination of results regarding the electron RO basis to get the combined expectation value
        y   = (y_pos_electron - y_neg_electron)/2.
        y_u = 1./2*(y_pos_electron_u**2 + y_neg_electron_u**2)**0.5


        ### Combine data
        exp_values = np.append(exp_values,y)
        exp_vals_u = np.append(exp_vals_u,y_u)

        #### we also need to calculate the single qubit expectation values.
        #### easy for the electron spin. One simply adds the electron outcomes and divides by 2.
        #### for the nuclear spin one has to add up several RO bases and look at the correlations there 
        #### (because experiments are not grouped by RO basis):
        exp_val_sum = (y_pos_electron + y_neg_electron)/2.
        exp_val_sum_u = (y_pos_electron_u**2 + y_neg_electron_u**2)/4.


        ### update the running average of the nuclear spin expectation values:
        c_x,c_y,c_z = c_x + exp_val_sum[0]/3.,c_y + exp_val_sum[1]/3., c_z + exp_val_sum[2]/3.
        c_x_u,c_y_u,c_z_u = c_x_u + exp_val_sum_u[0]/9.,c_y_u + exp_val_sum_u[1]/9., c_z_u + exp_val_sum_u[2]/9.

        ### add the electron values
        y_electron = np.sum(exp_val_sum)/3.
        y_electron_u = np.sum(exp_val_sum_u/9.)**0.5
        exp_values = np.append(exp_values,y_electron)
        exp_vals_u = np.append(exp_vals_u,y_electron_u)




        ### write up the combined analysis bases:
        a = mbi.MBIAnalysis(f_list_pos_p[-1])
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        e_base = np.array([tomo_pulse_translation_dict[p]]*4)
        c_base = a.adgrp.attrs['Tomography Bases'].flatten()
        comb_base = np.core.defchararray.add(e_base, np.append(c_base,'I')) ### add identity to 13C bases as this is evaluated on the fly
        sweep_pts = np.append(sweep_pts,comb_base)


    ### the nuclear spin single qubit correlations have been esimated in a running average (transposed correlations)
    ### and are now added to the results (square root still needs to be taken for the uncertainties)
    sweep_pts = np.append(sweep_pts,['IX','IY','IZ'])
    exp_values = np.append(exp_values,np.array([c_x,c_y,c_z]))
    print 'Nuclear spin correlations',c_x,c_y,c_z
    print sweep_pts
    print np.round(exp_values,2)

    exp_vals_u = np.append(exp_vals_u,np.array([np.sqrt(c_x_u),np.sqrt(c_y_u),np.sqrt(c_z_u)]))



    return f_list_pos_p[0],sweep_pts,exp_values,exp_vals_u

def carbon_ROC(exp,exp_u,folder):
    """
    setup dependent carbon read-correction coefficients are implemented here
    """
    if 'Pippin' in folder:
        ROC_coeff =  0.978264
        ROC_coeff_u = 0.00194222
    else:
        ROC_coeff = 0.972934 
        ROC_coeff_u = 0.0028265


    ### calc new uncertainty
    u = np.sqrt((ROC_coeff*exp_u)**2+(ROC_coeff_u*exp)**2)/ROC_coeff**2
    return exp/ROC_coeff,u

####################
#                  # 
#   Plotting etc.  #
#                  #
####################

def plot_dm(dm,dm_u_re = None,dm_u_im = None,plot_im = False):
    """
    routine for bar plotting
    """
    color = '#3594F2'
    alpha = 0.67


    xticks = [r'$|$Z,Z$\rangle$',r'$|$Z,-Z$\rangle$',r'$|$-Z,Z$\rangle$',r'$|$-Z,-Z$\rangle$']
    yticks = [r'$\langle$Z,Z$|$',r'$\langle$Z,-Z$|$',r'$\langle$-Z,Z$|$',r'$\langle$-Z,-Z$|$']
    fontsize = 10
    hf = plt.figure(figsize=plt.figaspect(0.3)/1.5)
    ha = plt.subplot(121, projection='3d')
    ha.grid(True)
    # plt.gca().patch.set_facecolor('white')
    ha.pbaspect = [1.0, 1.0, 0.255]
    xpos, ypos = np.array(range(4)),np.array(range(4))



    dx = 0.35 * np.ones(16)
    dy = dx.copy()


    a=np.arange(0,4,1)
    xpos,ypos = np.meshgrid(a,a)
    zpos = np.zeros((4,4))
    xpos = xpos.flatten()/2.
    ypos = ypos.flatten()/2.
    zpos = zpos.flatten()
    dz = np.reshape(np.asarray(np.abs(dm.real)), 16)


    #### now plot the error bars if given as input
    if dm_u_re != None:
        dm_err_re = np.reshape(np.asarray(dm_u_re), 16)
        for i in np.arange(0,len(xpos)):
            ha.plot([dx[i]/2+xpos[i],dx[i]/2+xpos[i]],[dy[i]/2+ypos[i],dy[i]/2+ypos[i]],[dz[i]-dm_err_re[i],dz[i]+dm_err_re[i]],marker="_",color = 'black')

    ha.bar3d(xpos, ypos, zpos, dx, dy,dz, color=color,alpha = alpha)

    # ha.set_title('abs(Real part)')
    ha.set_xticklabels(xticks,va = 'baseline',size=  fontsize,rotation = 40)
        #### fine adjustment of the x label positions... thanks stackexchange
    import types,matplotlib
    SHIFTX = 0.008 # Data coordinates
    SHIFTY = 0.003 # Data coordinates
    for label in ha.xaxis.get_majorticklabels():
        label.customShiftValueX = SHIFTX
        label.customShiftValueY = SHIFTY
        label.set_x = types.MethodType( lambda self, x: matplotlib.text.Text.set_x(self, x-self.customShiftValueX ), 
                                        label, matplotlib.text.Text )
        label.set_y = types.MethodType( lambda self, x: matplotlib.text.Text.set_y(self, x-self.customShiftValueY ), 
                                        label, matplotlib.text.Text )
    ha.set_yticklabels(yticks,size=  fontsize,
                   verticalalignment='baseline',
                   horizontalalignment='left',rotation =-15)
    SHIFT = 0.004 # Data coordinates
    for label in ha.yaxis.get_majorticklabels():
        label.customShiftValue = SHIFT
        label.set_x = types.MethodType( lambda self, x: matplotlib.text.Text.set_x(self, x-self.customShiftValue ), 
                                        label, matplotlib.text.Text )
    ha.set_zticklabels([0.0,0.1,0.2,0.3,0.4,0.5],size=  fontsize,
                   va='center',
                   ha ='left')
    ha.set_xticks([0.125,0.625,1.125,1.625])
    ha.set_yticks([0.125,0.625,1.125,1.625])
    ha.set_zlim([-0.0,0.5])

    if plot_im:
        dz = np.reshape(np.asarray(np.abs(dm.imag)), 16)

        hb = hf.add_subplot(122, projection='3d')
        hb.bar3d(xpos, ypos, zpos, dx, dy,dz, color=color,alpha = alpha)
        # hb.grid(False)
        #### now plot the error bars if given as input
        if dm_u_im != None:
            dm_err_im = np.reshape(np.asarray(dm_u_im), 16)
            for i in np.arange(0,len(xpos)):
                hb.plot([dx[i]/2+xpos[i],dx[i]/2+xpos[i]],[dy[i]/2+ypos[i],dy[i]/2+ypos[i]],[dz[i]-dm_err_im[i],dz[i]+dm_err_im[i]],marker="_",color = 'black')

        hb.set_title('abs(Imaginary part)')
        hb.set_xticklabels(xticks,va='center',size = fontsize)
        hb.set_yticklabels(yticks,size=  fontsize,
                       verticalalignment='baseline',
                       horizontalalignment='left')
        hb.set_zticklabels([0.0,0.1,0.2,0.3,0.4,0.5],size=  fontsize,
                       va='center',
                       ha ='left')
        hb.set_xticks([0.125,0.625,1.125,1.625])
        hb.set_yticks([0.125,0.625,1.125,1.625])
        hb.set_zlim([-0.0,0.5])

def electron_carbon_density_matrix(**kw):
    """
    Calculates a density matrix for our favourite nuclear spin-electron bell state
    Addendum: This function implements error propagation in the gaussian way and by assuming independence of measurement outcomes.
    """

    folder,sweep_pts,exp_values,exp_vals_u = get_correlations(**kw)
    print 'this is the folder', folder

    paulis = generate_pauli_matrices()
    ### initialize the dm via the correlations of a fully mixed state
    dm = np.kron(paulis[0],paulis[0])/4.
    dm_u_re,dm_u_im = np.zeros((4,4),dtype = float),np.zeros((4,4),dtype = float)
    ### basis definition
    t_dict = {'I' : 0, 'X':1, 'Y':2, 'Z':3}

    ### put dm together
    for t,exp,exp_u in zip(sweep_pts,exp_values,exp_vals_u):
            if t == 'II':
                continue
                
            sigma_kron = np.kron(paulis[t_dict[t[0]]],paulis[t_dict[t[1]]])


            ### carbon ROC necessary? YES
            if t[1] =='I':
                dm +=exp*sigma_kron/4.

                ### error calc is separate for real and imaginary part of the dm
                dm_u_re += ((exp_u/4.)**2)*np.abs(sigma_kron.real)
                dm_u_im += ((exp_u/4.)**2)*np.abs(sigma_kron.imag)

            else:
                exp,exp_u = carbon_ROC(exp,exp_u,folder)

                dm +=exp*sigma_kron/4.

                ### error calc is separate for real and imaginary part of the dm
                dm_u_re += ((exp_u/4.)**2)*np.abs(sigma_kron.real)
                dm_u_im += ((exp_u/4.)**2)*np.abs(sigma_kron.imag)

            # print t,(exp_u/4.)**2

    dm_u_re = np.sqrt(dm_u_re)
    dm_u_im = np.sqrt(dm_u_im)

    if kw.get('verbose',True):
        print 'Density matrix for the electron-carbon Bell state'
        print np.round(dm,decimals=3)
        print 'Expectation values from density matrix'
        print '      XX        YY          ZZ'
        print np.round(np.trace(np.dot(dm,np.kron(paulis[1],paulis[1]))),decimals=3),np.round(np.trace(np.dot(dm,np.kron(paulis[2],paulis[2]))),decimals=3),np.round(np.trace(np.dot(dm,np.kron(paulis[3],paulis[3]))),decimals=3)
        # print 'Error matrix'
        # print np.round(dm_p_u,decimals=3)
        print 'Eigenvalues'
        print np.linalg.eigh(dm)[0]

    if kw.pop('plot_errors',True):
        plot_dm(dm,dm_u_re,dm_u_im)
    else:
        plot_dm(dm)
